
/** \executable FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.21 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.21 2011/08/15 17:11:04 veelken Exp $
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"
#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

#include <vector>
#include <string>
#include <fstream>

typedef std::vector<std::string> vstring;

struct histManagerEntryType
{
  histManagerEntryType(const edm::ParameterSet& cfg, bool fillGenMatchHistograms,
		       const std::string& binVariable = "", double min = -0.5, double max = +0.5)
    : binVariable_(binVariable),
      min_(min),
      max_(max),
      fillGenMatchHistograms_(fillGenMatchHistograms)
  {
    histManager_ = new TauIdEffHistManager(cfg);

    if ( fillGenMatchHistograms_ ) {
      const std::string& label = cfg.getParameter<std::string>("label");

      edm::ParameterSet cfgJetToTauFake(cfg);
      std::string labelJetToTauFake = std::string(label).append("_").append("JetToTauFake");
      cfgJetToTauFake.addParameter<std::string>("label", labelJetToTauFake);
      histManagerJetToTauFake_ = new TauIdEffHistManager(cfgJetToTauFake);

      edm::ParameterSet cfgMuToTauFake(cfg);
      std::string labelMuToTauFake = std::string(label).append("_").append("MuToTauFake");
      cfgMuToTauFake.addParameter<std::string>("label", labelMuToTauFake);
      histManagerMuToTauFake_ = new TauIdEffHistManager(cfgMuToTauFake);

      edm::ParameterSet cfgGenTau(cfg);
      std::string labelGenTau = std::string(label).append("_").append("GenTau");
      cfgGenTau.addParameter<std::string>("label", labelGenTau);
      histManagerGenTau_ = new TauIdEffHistManager(cfgGenTau);
    }
  }
  ~histManagerEntryType() {}
  void bookHistograms(TFileDirectory& dir)
  {
    histManager_->bookHistograms(dir);

    if ( fillGenMatchHistograms_ ) {
      histManagerJetToTauFake_->bookHistograms(dir);
      histManagerMuToTauFake_->bookHistograms(dir);
      histManagerGenTau_->bookHistograms(dir);
    }
  }
  void fillHistograms(double x, const PATMuTauPair& muTauPair, size_t numVertices, int genMatchType, double weight)
  {
    if ( x > min_ && x <= max_ ) {
      histManager_->fillHistograms(muTauPair, numVertices, weight);

      if ( fillGenMatchHistograms_ ) {
	if      ( genMatchType == kJetToTauFakeMatched ) histManagerJetToTauFake_->fillHistograms(muTauPair, numVertices, weight);
	else if ( genMatchType == kMuToTauFakeMatched  ) histManagerMuToTauFake_->fillHistograms(muTauPair, numVertices, weight);
	else if ( genMatchType == kGenTauHadMatched    ||
		  genMatchType == kGenTauOtherMatched  ) histManagerGenTau_->fillHistograms(muTauPair, numVertices, weight);
      }
    }
  }

  std::string binVariable_;

  double min_;
  double max_;

  TauIdEffHistManager* histManager_;

  bool fillGenMatchHistograms_;

  TauIdEffHistManager* histManagerJetToTauFake_;
  TauIdEffHistManager* histManagerMuToTauFake_;
  TauIdEffHistManager* histManagerGenTau_;
};


struct regionEntryType
{
  regionEntryType(fwlite::TFileService& fs,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName, const std::string& sysShift,
		  const edm::ParameterSet& cfgBinning, const std::string& svFitMassHypothesis, 
		  const std::string& tauChargeMode, bool disableTauCandPreselCuts, bool fillGenMatchHistograms,
		  const std::string& selEventsFileName)
    : process_(process),
      region_(region),
      tauIdDiscriminators_(tauIdDiscriminators),
      tauIdName_(tauIdName),
      sysShift_(sysShift),
      selector_(0),
      histogramsUnbinned_(0),
      numMuTauPairs_selected_(0),
      numMuTauPairsWeighted_selected_(0.),
      selEventsFile_(0)
  {
    edm::ParameterSet cfgSelector;
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators_);
    cfgSelector.addParameter<std::string>("region", region_);
    cfgSelector.addParameter<std::string>("tauChargeMode", tauChargeMode);
    cfgSelector.addParameter<bool>("disableTauCandPreselCuts", disableTauCandPreselCuts);

    selector_ = new TauIdEffEventSelector(cfgSelector);

    edm::ParameterSet cfgHistManager;
    cfgHistManager.addParameter<std::string>("process", process_);
    cfgHistManager.addParameter<std::string>("region", region_);
    cfgHistManager.addParameter<std::string>("tauIdDiscriminator", tauIdName_);
    if      ( region.find("p") != std::string::npos ) label_ = "passed";
    else if ( region.find("f") != std::string::npos ) label_ = "failed";
    else                                              label_ = "all";
    if ( sysShift_ != "CENTRAL_VALUE" ) label_.append("_").append(sysShift_);
    cfgHistManager.addParameter<std::string>("label", label_);
    cfgHistManager.addParameter<std::string>("svFitMassHypothesis", svFitMassHypothesis);

    histogramsUnbinned_ = new histManagerEntryType(cfgHistManager, fillGenMatchHistograms);
    histogramsUnbinned_->bookHistograms(fs);

    typedef std::vector<edm::ParameterSet> vParameterSet;
    vstring binVariableNames = cfgBinning.getParameterNamesForType<vParameterSet>();
    for ( vstring::const_iterator binVariableName = binVariableNames.begin();
	  binVariableName != binVariableNames.end(); ++binVariableName ) {
      vParameterSet cfgBinVariableBins = cfgBinning.getParameter<vParameterSet>(*binVariableName);
      for ( vParameterSet::const_iterator cfgBinVariableBin = cfgBinVariableBins.begin();
	    cfgBinVariableBin != cfgBinVariableBins.end(); ++cfgBinVariableBin ) {
	double min = cfgBinVariableBin->getParameter<double>("min");
	double max = cfgBinVariableBin->getParameter<double>("max");
	histManagerEntryType* histManagerEntry = new histManagerEntryType(cfgHistManager, fillGenMatchHistograms, 
									  *binVariableName, min, max);
	std::string dir_string = cfgBinVariableBin->getParameter<std::string>("subdir");
	TFileDirectory dir = fs.mkdir(dir_string);
	histManagerEntry->bookHistograms(dir);
	histogramEntriesBinned_.push_back(histManagerEntry);
      }
    }

    if ( selEventsFileName != "" ) {
      size_t idx = selEventsFileName.rfind(".");
      if ( idx != std::string::npos ) {
	std::string selEventsFileName_region = std::string(selEventsFileName, 0, idx);
	selEventsFileName_region.append("_").append(region_);
	selEventsFileName_region.append(std::string(selEventsFileName, idx));
	//std::cout << "selEventsFileName_region = " << selEventsFileName_region << std::endl;
	selEventsFile_ = new std::ofstream(selEventsFileName_region.data(), std::ios::out);
      } else throw cms::Exception("regionEntryType")
	  << "Invalid selEventsFileName = " << selEventsFileName << " !!\n";
    }
  }
  ~regionEntryType()
  {
    delete selector_;

    delete histogramsUnbinned_;

    for ( std::vector<histManagerEntryType*>::iterator it = histogramEntriesBinned_.begin();
	  it != histogramEntriesBinned_.end(); ++it ) {
      delete (*it);
    }
    
    delete selEventsFile_;
  }
  void analyze(const fwlite::Event& evt, const PATMuTauPair& muTauPair, size_t numVertices, int genMatchType, double evtWeight)
  {
    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, evtSelFlags) ) {
//--- fill histograms for "inclusive" tau id. efficiency measurement
      histogramsUnbinned_->fillHistograms(0., muTauPair, numVertices, genMatchType, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of 
//   o tau-jet transverse momentum
//   o tau-jet pseudo-rapidity
//   o reconstructed vertex multiplicity
//   o sumEt
//   o ...
      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = histogramEntriesBinned_.begin();
	    histManagerEntry != histogramEntriesBinned_.end(); ++histManagerEntry ) {
	double x = 0.;
	if      ( (*histManagerEntry)->binVariable_ == "tauPt"       ) x = muTauPair.leg2()->pt();
	else if ( (*histManagerEntry)->binVariable_ == "tauAbsEta"   ) x = TMath::Abs(muTauPair.leg2()->eta());
	else if ( (*histManagerEntry)->binVariable_ == "numVertices" ) x = numVertices;
	else if ( (*histManagerEntry)->binVariable_ == "sumEt"       ) x = muTauPair.met()->sumEt();
	else throw cms::Exception("regionEntryType::analyze")
	  << "Invalid binVariable = " << (*histManagerEntry)->binVariable_ << " !!\n";
	(*histManagerEntry)->fillHistograms(x, muTauPair, numVertices, genMatchType, evtWeight);
      }
 
      if ( selEventsFile_ ) 
	(*selEventsFile_) << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() << std::endl;

      ++numMuTauPairs_selected_;
      numMuTauPairsWeighted_selected_ += evtWeight;
    }
  }

  std::string process_;
  std::string region_;
  vstring tauIdDiscriminators_;
  std::string tauIdName_;
  std::string sysShift_;
  std::string label_;

  TauIdEffEventSelector* selector_;

  histManagerEntryType* histogramsUnbinned_;
  std::vector<histManagerEntryType*> histogramEntriesBinned_;

  int numMuTauPairs_selected_;
  double numMuTauPairsWeighted_selected_;

  std::ofstream* selEventsFile_;
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
    throw cms::Exception("FWLiteTauIdEffAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauIdEffAnalyzer = cfg.getParameter<edm::ParameterSet>("tauIdEffAnalyzer");

  edm::InputTag srcMuTauPairs = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcMuTauPairs");
  std::string svFitMassHypothesis = cfgTauIdEffAnalyzer.getParameter<std::string>("svFitMassHypothesis");
  std::string tauChargeMode = cfgTauIdEffAnalyzer.getParameter<std::string>("tauChargeMode");
  bool disableTauCandPreselCuts = cfgTauIdEffAnalyzer.getParameter<bool>("disableTauCandPreselCuts");
  edm::InputTag srcTrigger = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcTrigger");
  vstring hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("hltPaths");
  edm::InputTag srcGoodMuons = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGoodMuons");
  edm::InputTag srcVertices = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcVertices");
  edm::InputTag srcGenParticles = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGenParticles");
  bool fillGenMatchHistograms = cfgTauIdEffAnalyzer.getParameter<bool>("fillGenMatchHistograms");
  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch = cfgTauIdEffAnalyzer.getParameter<vint>("skipPdgIdsGenParticleMatch");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffAnalyzer.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffAnalyzer.exists("sysShift") ?
    cfgTauIdEffAnalyzer.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";
  edm::InputTag srcEventCounter = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcEventCounter");

  std::string selEventsFileName = ( cfgTauIdEffAnalyzer.exists("selEventsFileName") ) ? 
    cfgTauIdEffAnalyzer.getParameter<std::string>("selEventsFileName") : "";

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauIdEffAnalyzer.getParameter<std::string>("process");
  std::cout << " process = " << process << std::endl;
  std::string processType = cfgTauIdEffAnalyzer.getParameter<std::string>("type");
  std::cout << " type = " << processType << std::endl;
  bool isData = (processType == "Data");
  vstring regions = cfgTauIdEffAnalyzer.getParameter<vstring>("regions");
  edm::ParameterSet cfgBinning = cfgTauIdEffAnalyzer.getParameter<edm::ParameterSet>("binning");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauIdEffAnalyzer.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");
      regionEntryType* regionEntry = 
	new regionEntryType(fs, process, *region, tauIdDiscriminators, tauIdName, 
			    sysShift, cfgBinning, svFitMassHypothesis, 
			    tauChargeMode, disableTauCandPreselCuts, fillGenMatchHistograms, selEventsFileName);
      regionEntries.push_back(regionEntry);
    }
  }

  edm::ParameterSet cfgSelectorABCD;
  cfgSelectorABCD.addParameter<vstring>("tauIdDiscriminators", vstring());
  cfgSelectorABCD.addParameter<std::string>("region", "ABCD");
  cfgSelectorABCD.addParameter<std::string>("tauChargeMode", tauChargeMode);
  cfgSelectorABCD.addParameter<bool>("disableTauCandPreselCuts", disableTauCandPreselCuts);
  TauIdEffEventSelector* selectorABCD = new TauIdEffEventSelector(cfgSelectorABCD);

//--- book "dummy" histogram counting number of processed events
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 3, -0.5, +2.5);
  histogramEventCounter->GetXaxis()->SetBinLabel(1, "all Events (DBS)");      // CV: bin numbers start at 1 (not 0) !!
  histogramEventCounter->GetXaxis()->SetBinLabel(2, "processed by Skimming");
  histogramEventCounter->GetXaxis()->SetBinLabel(3, "analyzed in PAT-tuple");
  
  int allEvents_DBS = cfgTauIdEffAnalyzer.getParameter<int>("allEvents_DBS");
  if ( allEvents_DBS > 0 ) {
    histogramEventCounter->SetBinContent(1, allEvents_DBS);
  } else {
    histogramEventCounter->SetBinContent(1, -1.);
  }
  
  double xSection = cfgTauIdEffAnalyzer.getParameter<double>("xSection");
  double intLumiData = cfgTauIdEffAnalyzer.getParameter<double>("intLumiData");

  int    numEvents_processed                     = 0; 
  double numEventsWeighted_processed             = 0.;
  int    numEvents_passedTrigger                 = 0;
  double numEventsWeighted_passedTrigger         = 0.;
  int    numEvents_passedDiMuonVeto              = 0;
  double numEventsWeighted_passedDiMuonVeto      = 0.;
  int    numEvents_passedDiMuTauPairVeto         = 0;
  double numEventsWeighted_passedDiMuTauPairVeto = 0.;

  edm::RunNumber_t lastLumiBlock_run = -1;
  edm::LuminosityBlockNumber_t lastLumiBlock_ls = -1;

  double intLumiData_analyzed = 0.;
  edm::InputTag srcLumiProducer = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcLumiProducer");

  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauIdEffAnalyzer") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << "opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
      double evtWeight = 1.0;
      for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	    srcWeight != srcWeights.end(); ++srcWeight ) {
	edm::Handle<double> weight;
	evt.getByLabel(*srcWeight, weight);
	evtWeight *= (*weight);
      }

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

      //std::cout << "processing run = " << evt.id().run() << ":" 
      //	  << " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << std::endl;

//--- check if new luminosity section has started;
//    if so, retrieve number of events contained in this luminosity section before skimming
      if ( !(evt.id().run() == lastLumiBlock_run && evt.luminosityBlock() == lastLumiBlock_ls) ) {
	const fwlite::LuminosityBlock& ls = evt.getLuminosityBlock();
	edm::Handle<edm::MergeableCounter> numEvents_skimmed;
	ls.getByLabel(srcEventCounter, numEvents_skimmed);
	if ( numEvents_skimmed.isValid() ) histogramEventCounter->Fill(1, numEvents_skimmed->value);
	lastLumiBlock_run = evt.id().run();
	lastLumiBlock_ls = evt.luminosityBlock();

	if ( isData ) {
	  edm::Handle<LumiSummary> lumiSummary;
	  edm::InputTag srcLumiProducer("lumiProducer");
	  ls.getByLabel(srcLumiProducer, lumiSummary);
	  intLumiData_analyzed += lumiSummary->intgRecLumi();
	}
      }

//--- fill "dummy" histogram counting number of processed events
      histogramEventCounter->Fill(2);

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

      if ( !isTriggered ) continue;
      ++numEvents_passedTrigger;
      numEventsWeighted_passedTrigger += evtWeight;

//--- require event to contain only one "good quality" muon
      typedef std::vector<pat::Muon> PATMuonCollection;
      edm::Handle<PATMuonCollection> goodMuons;
      evt.getByLabel(srcGoodMuons, goodMuons);
      size_t numGoodMuons = goodMuons->size();
	
      if ( !(numGoodMuons <= 1) ) continue;
      ++numEvents_passedDiMuonVeto;
      numEventsWeighted_passedDiMuonVeto += evtWeight;

//--- require event to contain exactly one muon + tau-jet pair
//    passing the selection criteria for region "ABCD"
      edm::Handle<PATMuTauPairCollection> muTauPairs;
      evt.getByLabel(srcMuTauPairs, muTauPairs);

      unsigned numMuTauPairsABCD = 0;
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {
	pat::strbitset evtSelFlags;
	if ( selectorABCD->operator()(*muTauPair, evtSelFlags) ) ++numMuTauPairsABCD;
      }
      
      if ( !(numMuTauPairsABCD <= 1) ) continue;
      ++numEvents_passedDiMuTauPairVeto;
      numEventsWeighted_passedDiMuTauPairVeto += evtWeight;

//--- determine number of vertices reconstructed in the event
//   (needed to parametrize dependency of tau id. efficiency on number of pile-up interactions)
      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      size_t numVertices = vertices->size();

//--- iterate over collection of muon + tau-jet pairs:
//    check which region muon + tau-jet pair is selected in,
//    fill histograms for that region
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {

//--- determine type of particle matching reconstructed tau-jet candidate
//    on generator level (used in case of Ztautau or Zmumu Monte Carlo samples only,
//    in order to distinguish between jet --> tau fakes, muon --> tau fakes and genuine taus)
	int genMatchType = kUnmatched;
	if ( fillGenMatchHistograms ) {
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  evt.getByLabel(srcGenParticles, genParticles);	  
	  genMatchType = getGenMatchType(*muTauPair, *genParticles);
	}

	for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	      regionEntry != regionEntries.end(); ++regionEntry ) {	  
	  (*regionEntry)->analyze(evt, *muTauPair, numVertices, genMatchType, evtWeight);

	  //pat::strbitset evtSelFlags;
	  //if ( (*regionEntry)->region_ == "D1p" && (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	  //  std::cout << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() 
	  //	        << " (weight = " << evtWeight << ")" << std::endl;
	  //}
	}
      }
    }

//--- close input file
    delete inputFile;
  }

//--- scale histograms taken from Monte Carlo simulation
//    according to cross-section times luminosity
  if ( !isData ) {
    double mcScaleFactor = (intLumiData*xSection)/(double)allEvents_DBS;
    std::cout << " intLumiData = " << intLumiData << std::endl;
    std::cout << " xSection = " << xSection << std::endl;
    std::cout << " allEvents_DBS = " << allEvents_DBS << std::endl;
    std::cout << "--> scaling histograms by factor = " << mcScaleFactor
	      << " according to cross-section times luminosity." << std::endl;

//--- apply correction to scale-factor in order to account for events lost, 
//    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs
    double lostStatCorrFactor = 1.;
    if ( histogramEventCounter->GetBinContent(1) > histogramEventCounter->GetBinContent(2) && 
	 histogramEventCounter->GetBinContent(2) > 0.                                      ) {
      lostStatCorrFactor = histogramEventCounter->GetBinContent(1)/histogramEventCounter->GetBinContent(2);
      std::cout << "--> scaling histograms by additional factor = " << lostStatCorrFactor
		<< " to account for events lost," << std::endl; 
      std::cout << "    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs." << std::endl;
    }

    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {  
      (*regionEntry)->histogramsUnbinned_->histManager_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = (*regionEntry)->histogramEntriesBinned_.begin();
	    histManagerEntry != (*regionEntry)->histogramEntriesBinned_.end(); ++histManagerEntry ) {
	(*histManagerEntry)->histManager_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      }
    }
  }

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;
  std::cout << " numEvents_passedTrigger: " << numEvents_passedTrigger 
	    << " (weighted = " << numEventsWeighted_passedTrigger << ")" << std::endl;
  std::cout << " numEvents_passedDiMuonVeto: " << numEvents_passedDiMuonVeto 
	    << " (weighted = " << numEventsWeighted_passedDiMuonVeto << ")" << std::endl;
  std::cout << " numEvents_passedDiMuTauPairVeto: " << numEvents_passedDiMuTauPairVeto
	    << " (weighted = " << numEventsWeighted_passedDiMuTauPairVeto << ")" << std::endl;
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

//--- close ASCII files containing run + event numbers of events selected in different regions
  for ( std::vector<regionEntryType*>::iterator it = regionEntries.begin();
	it != regionEntries.end(); ++it ) {
    delete (*it);
  }
  
  if ( isData ) {
    std::cout << " intLumiData (recorded, IsoMu17 prescale corr.) = " << intLumiData << " pb" << std::endl;
    // CV: luminosity is recorded in some 'weird' units,
    //     needs to be multiplied by factor 0.10 in order to be in units of pb^-1
    std::cout << " intLumiData_analyzed (recorded) = " << intLumiData_analyzed*1.e-6*0.10 << " pb" << std::endl;
  }

  clock.Show("FWLiteTauIdEffAnalyzer");

  return 0;
}
