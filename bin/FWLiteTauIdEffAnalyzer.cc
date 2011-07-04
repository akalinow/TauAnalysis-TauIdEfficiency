
/** \class FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.5 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.5 2011/07/03 10:15:55 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"

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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/Handle.h"

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
      histManagerTauEtaLt15_(0),
      histManagerTauEta15to19_(0),
      histManagerTauEta19to23_(0),
      histManagerTauPtLt25_(0),
      histManagerTauPt25to30_(0),
      histManagerTauPt30to40_(0),
      histManagerTauPtGt40_(0),
      histManagerSumEtLt75_(0),
      histManagerSumEt75to150_(0),
      histManagerSumEt150to225_(0),
      histManagerSumEtGt225_(0),
      histManagerNumVerticesLeq3_(0),
      histManagerNumVertices4to6_(0),
      histManagerNumVertices7to9_(0),
      histManagerNumVerticesGt10_(0),
      numMuTauPairs_selected_(0),
      numMuTauPairsWeighted_selected_(0.)
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

    histManager_                = addHistManager(fs, "",                cfgHistManager);

    histManagerTauEtaLt15_      = addHistManager(fs, "tauEtaLt15",      cfgHistManager);
    histManagerTauEta15to19_    = addHistManager(fs, "tauEta15to19",    cfgHistManager);
    histManagerTauEta19to23_    = addHistManager(fs, "tauEta19to23",    cfgHistManager);

    histManagerTauPtLt25_       = addHistManager(fs, "tauPtLt25",       cfgHistManager);
    histManagerTauPt25to30_     = addHistManager(fs, "tauPt25to30",     cfgHistManager);
    histManagerTauPt30to40_     = addHistManager(fs, "tauPt30to40",     cfgHistManager);
    histManagerTauPtGt40_       = addHistManager(fs, "tauPtGt40",       cfgHistManager);

    histManagerSumEtLt75_       = addHistManager(fs, "sumEtLt75",       cfgHistManager);
    histManagerSumEt75to150_    = addHistManager(fs, "sumEt75to150",    cfgHistManager);
    histManagerSumEt150to225_   = addHistManager(fs, "sumEt150to225",   cfgHistManager);
    histManagerSumEtGt225_      = addHistManager(fs, "sumEtGt225",      cfgHistManager);

    histManagerNumVerticesLeq3_ = addHistManager(fs, "numVerticesLeq3", cfgHistManager);
    histManagerNumVertices4to6_ = addHistManager(fs, "numVertices4to6", cfgHistManager);
    histManagerNumVertices7to9_ = addHistManager(fs, "numVertices7to9", cfgHistManager);
    histManagerNumVerticesGt10_ = addHistManager(fs, "numVerticesGt10", cfgHistManager);
  }
  ~regionEntryType()
  {
    delete selector_;
    delete histManager_;
  }
  TauIdEffHistManager* addHistManager(fwlite::TFileService& fs, const std::string& dirName, const edm::ParameterSet& cfg)
  {
    TauIdEffHistManager* retVal = new TauIdEffHistManager(cfg);

    if ( dirName != "" ) {
      TFileDirectory dir = fs.mkdir(dirName.data());
      retVal->bookHistograms(dir);
    } else {
      retVal->bookHistograms(fs);
    }

    return retVal;
  }
  void analyze(const PATMuTauPair& muTauPair, size_t numVertices, double evtWeight)
  {
    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, evtSelFlags) ) {
//--- fill histograms for "inclusive" tau id. efficiency measurement
      histManager_->fillHistograms(muTauPair, numVertices, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of tau-jet pseudo-rapidity
      double tauAbsEta = TMath::Abs(muTauPair.leg1()->eta());
      if      ( tauAbsEta < 1.5 ) histManagerTauEtaLt15_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( tauAbsEta < 1.9 ) histManagerTauEta15to19_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( tauAbsEta < 2.3 ) histManagerTauEta19to23_->fillHistograms(muTauPair, numVertices, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of tau-jet transverse momentum
      double tauPt = muTauPair.leg1()->pt();
      if      ( tauPt > 20.0 && tauPt < 25.0 ) histManagerTauPtLt25_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( tauPt > 25.0 && tauPt < 30.0 ) histManagerTauPt25to30_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( tauPt > 30.0 && tauPt < 40.0 ) histManagerTauPt30to40_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( tauPt > 40.0 && tauPt < 80.0 ) histManagerTauPtGt40_->fillHistograms(muTauPair, numVertices, evtWeight);
      
//--- fill histograms for tau id. efficiency measurement as function of sumEt
      double sumEt = muTauPair.met()->sumEt();
      if      ( sumEt <  75.0 ) histManagerSumEtLt75_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( sumEt < 150.0 ) histManagerSumEt75to150_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( sumEt < 225.0 ) histManagerSumEt150to225_->fillHistograms(muTauPair, numVertices, evtWeight);
      else                      histManagerSumEtGt225_->fillHistograms(muTauPair, numVertices, evtWeight);
      
//--- fill histograms for tau id. efficiency measurement as function of reconstructed vertex multiplicity
      if      ( numVertices <=  3 ) histManagerNumVerticesLeq3_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( numVertices <=  6 ) histManagerNumVertices4to6_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( numVertices <=  9 ) histManagerNumVertices7to9_->fillHistograms(muTauPair, numVertices, evtWeight);
      else if ( numVertices <= 20 ) histManagerNumVerticesGt10_->fillHistograms(muTauPair, numVertices, evtWeight);

      ++numMuTauPairs_selected_;
      numMuTauPairsWeighted_selected_ += evtWeight;
    }
  }
  std::string process_;
  std::string region_;
  std::string tauIdName_;

  TauIdEffEventSelector* selector_;

  TauIdEffHistManager* histManager_;

  TauIdEffHistManager* histManagerTauEtaLt15_;
  TauIdEffHistManager* histManagerTauEta15to19_;
  TauIdEffHistManager* histManagerTauEta19to23_;

  TauIdEffHistManager* histManagerTauPtLt25_;
  TauIdEffHistManager* histManagerTauPt25to30_;
  TauIdEffHistManager* histManagerTauPt30to40_;
  TauIdEffHistManager* histManagerTauPtGt40_;

  TauIdEffHistManager* histManagerSumEtLt75_;
  TauIdEffHistManager* histManagerSumEt75to150_;
  TauIdEffHistManager* histManagerSumEt150to225_;
  TauIdEffHistManager* histManagerSumEtGt225_;

  TauIdEffHistManager* histManagerNumVerticesLeq3_;
  TauIdEffHistManager* histManagerNumVertices4to6_;
  TauIdEffHistManager* histManagerNumVertices7to9_;
  TauIdEffHistManager* histManagerNumVerticesGt10_;

  int numMuTauPairs_selected_;
  double numMuTauPairsWeighted_selected_;
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
  edm::InputTag srcVertices = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcVertices");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffAnalyzer.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffAnalyzer.exists("sysShift") ?
    cfgTauIdEffAnalyzer.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";
  edm::InputTag srcEventCounter = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcEventCounter");

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
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 3, -0.5, +2.5);
  histogramEventCounter->GetXaxis()->SetBinLabel(1, "all Events (DBS)");      // CV: bin numbers start at 1 (not 0) !!
  histogramEventCounter->GetXaxis()->SetBinLabel(2, "processed by Skimming");
  histogramEventCounter->GetXaxis()->SetBinLabel(3, "analyzed in PAT-tuple");
  
  if ( cfgTauIdEffAnalyzer.exists("allEvents_DBS") ) {
    histogramEventCounter->SetBinContent(1, cfgTauIdEffAnalyzer.getParameter<int>("allEvents_DBS"));
  } else {
    histogramEventCounter->SetBinContent(1, -1.);
  }
  
  int numEvents_processed = 0; 
  double numEventsWeighted_processed = 0;
  
  edm::RunNumber_t lastLumiBlock_run = -1;
  edm::LuminosityBlockNumber_t lastLumiBlock_ls = -1;

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

//--- check if new luminosity section has started;
//    if so, retrieve number of events contained in this luminosity section before skimming
      if ( !(evt.id().run() == lastLumiBlock_run && evt.luminosityBlock() == lastLumiBlock_ls) ) {
	const fwlite::LuminosityBlock& ls = evt.getLuminosityBlock();
	edm::Handle<edm::MergeableCounter> numEvents_skimmed;
	ls.getByLabel(srcEventCounter, numEvents_skimmed);
	if ( numEvents_skimmed.isValid() ) histogramEventCounter->Fill(1, numEvents_skimmed->value);
	lastLumiBlock_run = evt.id().run();
	lastLumiBlock_ls = evt.luminosityBlock();
      }

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

//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
      double evtWeight = 1.0;
      for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	    srcWeight != srcWeights.end(); ++srcWeight ) {
	edm::Handle<double> weight;
	evt.getByLabel(*srcWeight, weight);
	evtWeight *= (*weight);
      }

//--- require event to pass trigger requirements
//    and to contain only one "good quality" muon
      typedef std::vector<pat::Muon> PATMuonCollection;
      edm::Handle<PATMuonCollection> goodMuons;
      evt.getByLabel(srcGoodMuons, goodMuons);
      size_t numGoodMuons = goodMuons->size();

      if ( isTriggered && numGoodMuons <= 1 ) {
	
	edm::Handle<reco::VertexCollection> vertices;
	evt.getByLabel(srcVertices, vertices);
	size_t numVertices = vertices->size();

//--- iterate over collection of muon + tau-jet pairs
	edm::Handle<PATMuTauPairCollection> muTauPairs;
	evt.getByLabel(srcMuTauPairs, muTauPairs);

	for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	      muTauPair != muTauPairs->end(); ++muTauPair ) {
	  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
		regionEntry != regionEntries.end(); ++regionEntry ) {
	    (*regionEntry)->analyze(*muTauPair, numVertices, evtWeight);
	  }
	}
      }
      
//--- fill "dummy" histogram counting number of processed events
      histogramEventCounter->Fill(2);

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;
    }

//--- close input file
    delete inputFile;
  }

//--- scale histograms to account for events lost, 
//    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs
  if ( histogramEventCounter->GetBinContent(1) > histogramEventCounter->GetBinContent(2) && 
       histogramEventCounter->GetBinContent(2) > 0.                                      ) {
    double factor = histogramEventCounter->GetBinContent(1)/histogramEventCounter->GetBinContent(2);
    std::cout << "--> scaling histograms by factor = " << factor 
	      << " to account for events lost," 
	      << " due to aborted skimming/crab or PAT-tuple production/lxbatch jobs." << std::endl;
    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {
      (*regionEntry)->histManager_->scaleHistograms(factor);
    }
  }

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed << std::endl;
  std::string lastTauIdName = "";
  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    if ( (*regionEntry)->tauIdName_ != lastTauIdName ) 
      std::cout << " numMuTauPairs_selected, " << (*regionEntry)->tauIdName_ << std::endl;
    std::cout << "  region " << (*regionEntry)->region_ << ":" 
	      << " " << (*regionEntry)->numMuTauPairs_selected_ 
	      << " (weighted = " << (*regionEntry)->numMuTauPairsWeighted_selected_ << ")" << std::endl;
    lastTauIdName = (*regionEntry)->tauIdName_;
  }

  clock.Show("FWLiteTauIdEffAnalyzer");

  return 0;
}
