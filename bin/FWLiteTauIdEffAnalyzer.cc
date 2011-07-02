
/** \class FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.2 2011/07/02 10:42:44 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

typedef std::vector<std::string> vstring;

struct regionEntryType
{
  regionEntryType(fwlite::TFileService& fs,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName, const std::string& sysShift)
    : process_(process),
      region_(region),
      tauIdName_(tauIdName),
      selector_(0),
      histManager_(0),
      numMuTauPairs_selected_(0)
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
    if ( sysShift != "CENTRAL_VALUE" ) label.append("_").append(sysShift);
    cfgHistManager.addParameter<std::string>("label", label);
    histManager_ = new TauIdEffHistManager(cfgHistManager);
    histManager_->bookHistograms(fs);
  }
  ~regionEntryType()
  {
    delete selector_;
    delete histManager_;
  }
  std::string process_;
  std::string region_;
  std::string tauIdName_;
  TauIdEffEventSelector* selector_;
  TauIdEffHistManager* histManager_;
  int numMuTauPairs_selected_;
};

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauIdEffAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("TauIdEffEventSelector") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauIdEffAnalyzer = cfg.getParameter<edm::ParameterSet>("tauIdEffAnalyzer");

  edm::InputTag srcMuTauPairs = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcMuTauPairs");
  edm::InputTag srcTrigger = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcTrigger");
  vstring hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("hltPaths");
  edm::InputTag srcGoodMuons = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGoodMuons");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffAnalyzer.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffAnalyzer.exists("sysShift") ?
    cfgTauIdEffAnalyzer.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauIdEffAnalyzer.getParameter<std::string>("process");
  vstring regions = cfgTauIdEffAnalyzer.getParameter<vstring>("regions");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauIdEffAnalyzer.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");
      regionEntryType* regionEntry = new regionEntryType(fs, process, *region, tauIdDiscriminators, tauIdName, sysShift);
      regionEntries.push_back(regionEntry);
    }
  }

//--- book "dummy" histogram counting number of processed events
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 1, -0.5, +0.5);
  
  int numEvents_processed = 0;  
  
  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("TauIdEffEventSelector") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << " opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

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

      typedef std::vector<pat::Muon> PATMuonCollection;
      edm::Handle<PATMuonCollection> goodMuons;
      evt.getByLabel(srcGoodMuons, goodMuons);
      size_t numGoodMuons = goodMuons->size();

//--- require event to pass trigger requirements
//    and to contain only one "good quality" muon
      if ( isTriggered && numGoodMuons <= 1 ) {

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
	  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
		regionEntry != regionEntries.end(); ++regionEntry ) {
	    //std::cout << "checking region = " << (*regionEntry)->region_ << std::endl;
	    pat::strbitset evtSelFlags;
	    if ( (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	      //std::cout << "--> selection passed !!" << std::endl;
	      (*regionEntry)->histManager_->fillHistograms(*muTauPair, evtWeight);
	      ++(*regionEntry)->numMuTauPairs_selected_;
	    }
	  }
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

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed << std::endl;
  std::string lastTauIdName = "";
  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    if ( (*regionEntry)->tauIdName_ != lastTauIdName ) 
      std::cout << " numMuTauPairs_selected, " << (*regionEntry)->tauIdName_ << std::endl;
    std::cout << "  region " << (*regionEntry)->region_ << ": " << (*regionEntry)->numMuTauPairs_selected_ << std::endl;
    lastTauIdName = (*regionEntry)->tauIdName_;
  }

  clock.Show("FWLiteTauIdEffAnalyzer");

  return 0;
}
