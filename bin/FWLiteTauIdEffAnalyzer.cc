
/** \class FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.1 2011/07/01 10:41:48 veelken Exp $
 *
 */

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

typedef std::vector<std::string> vstring;

struct regionEntryType
{
  regionEntryType(const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName)
    : selector_(0),
      histManager_(0)
  {
    edm::ParameterSet cfgSelector;
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators);
    cfgSelector.addParameter<std::string>("region", region);
    selector_ = new TauIdEffEventSelector(cfgSelector);
    edm::ParameterSet cfgHistManager;
    cfgHistManager.addParameter<std::string>("process", process);
    cfgHistManager.addParameter<std::string>("region", region);
    cfgHistManager.addParameter<std::string>("tauIdDiscriminator", tauIdName);
    std::string label;
    if      ( region.find("p") != std::string::npos ) label = "passed";
    else if ( region.find("f") != std::string::npos ) label = "failed";
    else                                              label = "all";
    cfgHistManager.addParameter<std::string>("label", label);
  }
  ~regionEntryType()
  {
    delete selector_;
    delete histManager_;
  }
  TauIdEffEventSelector* selector_;
  TauIdEffHistManager* histManager_;
};

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("TauIdEffEventSelector") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  const edm::ParameterSet& cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::InputTag srcMuTauPairs = cfg.getParameter<edm::InputTag>("srcMuTauPairs");
  edm::InputTag srcTrigger = cfg.getParameter<edm::InputTag>("srcTrigger");
  vstring hltPaths = cfg.getParameter<vstring>("hltPaths");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfg.getParameter<vInputTag>("weights");

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfg.getParameter<std::string>("process");
  vstring regions = cfg.getParameter<vstring>("regions");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfg.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");
      regionEntryType* regionEntry = new regionEntryType(process, *region, tauIdDiscriminators, tauIdName);
      regionEntries.push_back(regionEntry);
    }
  }

//--- book "dummy" histogram counting number of processed events
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 1, -0.5, +0.5);
  
  int numEvents_processed = 0;  
  int numMuTauPairs_selected = 0;
  
  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("TauIdEffEventSelector") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

//--- check that event has passed triggers
      edm::Handle<pat::TriggerEvent> hltEvent;
      evt.getByLabel(srcTrigger, hltEvent);
  
      bool isTriggered = false;
      for ( vstring::const_iterator hltPathName = hltPaths.begin();
	    hltPathName != hltPaths.end() && !isTriggered; ++hltPathName ) {
	if ( (*hltPathName) == "*" ) { // check for wildcard character "*" that accepts all events
	  isTriggered = true;
	} else {
	  const pat::TriggerPath* hltPath = hltEvent->path(*hltPathName);
	  if ( hltPath && hltPath->wasAccept() ) isTriggered = true;
	}
      }

      if ( isTriggered ) {
//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
	double evtWeight = 1.0;
	for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	      srcWeight != srcWeights.end(); ++srcWeight ) {
	  edm::Handle<double> weight;
	  evt.getByLabel(*srcWeight, weight);
	  evtWeight *= (*weight);
	}

//--- iterator over collection of muon + tau-jet pairs
	edm::Handle<PATMuTauPairCollection> muTauPairs;
	evt.getByLabel(srcMuTauPairs, muTauPairs);

	for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	      muTauPair != muTauPairs->end(); ++muTauPair ) {
	  bool isMuTauPair_selected;
	  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
		regionEntry != regionEntries.end(); ++regionEntry ) {
	    pat::strbitset evtSelFlags;
	    if ( (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	      (*regionEntry)->histManager_->fillHistograms(*muTauPair, evtWeight);
	      isMuTauPair_selected = true;
	    }
	  }
	  if ( isMuTauPair_selected ) ++numMuTauPairs_selected;
	}
      }
      
//--- fill "dummy" histogram counting number of processed events
      histogramEventCounter->Fill(0);

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;
    }

//--- close input file
    delete inputFile;
  }

  std::cout << " numEvents_processed = " << numEvents_processed << std::endl;
  std::cout << " numMuTauPairs_selected(ABCD) = " << numMuTauPairs_selected << std::endl;

  return 0;
}
