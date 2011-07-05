
/** \class FWLiteTauIdEffPreselNumbers
 *
 * Determine efficiency of "leading" track finding, leading track Pt and loose (PF)isolation requirements 
 * applied in preselection of tau-jet candidates considered for tau id. efficiency measurement.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: FWLiteTauIdEffPreselNumbers.cc,v 1.1 2011/07/05 12:30:16 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"

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
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffCutFlowTable.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

typedef std::vector<std::string> vstring;

enum { kUnmatched, kTauHadMatched, kFakeTauMatched };

struct cutFlowEntryType
{
  cutFlowEntryType(const edm::ParameterSet& cfg)
  {
    std::string label = cfg.getParameter<std::string>("label");

    edm::ParameterSet cfgTauHadMatched = cfg;
    std::string labelTauHadMatched = std::string(label).append("TauHadMatched");
    cfgTauHadMatched.addParameter<std::string>("label", labelTauHadMatched);
    cutFlowTauHadMatched_ = new TauIdEffCutFlowTable(cfgTauHadMatched);
    edm::ParameterSet cfgFakeTauMatched = cfg;
    std::string labelFakeTauMatched = std::string(label).append("FakeTauMatched");
    cfgFakeTauMatched.addParameter<std::string>("label", labelFakeTauMatched);
    cutFlowFakeTauMatched_ = new TauIdEffCutFlowTable(cfgFakeTauMatched);

    edm::ParameterSet cfgNoMatchingApplied = cfg;
    std::string labelNoMatchingApplied = std::string(label).append("NoMatchingApplied");
    cfgNoMatchingApplied.addParameter<std::string>("label", labelNoMatchingApplied);
    cutFlowNoMatchingApplied_ = new TauIdEffCutFlowTable(cfgNoMatchingApplied);
  }
  ~cutFlowEntryType()
  {
    delete cutFlowTauHadMatched_;
    delete cutFlowFakeTauMatched_;

    delete cutFlowNoMatchingApplied_;
  }
  void bookCutFlowTables(TFileDirectory& dir, int numBins, float* binning)
  {
    cutFlowTauHadMatched_->bookCutFlowTable(dir, numBins, binning);
    cutFlowFakeTauMatched_->bookCutFlowTable(dir, numBins, binning);
    cutFlowNoMatchingApplied_->bookCutFlowTable(dir, numBins, binning);
  }
  void fillCutFlowTables(double x, const std::vector<bool>& selectionFlags, int genMatchType, double evtWeight)
  {
    if      ( genMatchType == kTauHadMatched  ) cutFlowTauHadMatched_->fillCutFlowTable(x, selectionFlags, evtWeight);
    else if ( genMatchType == kFakeTauMatched ) cutFlowFakeTauMatched_->fillCutFlowTable(x, selectionFlags, evtWeight);

    cutFlowNoMatchingApplied_->fillCutFlowTable(x, selectionFlags, evtWeight);
  }

  TauIdEffCutFlowTable* cutFlowTauHadMatched_;
  TauIdEffCutFlowTable* cutFlowFakeTauMatched_;

  TauIdEffCutFlowTable* cutFlowNoMatchingApplied_;
};

struct regionEntryType
{
  regionEntryType(fwlite::TFileService& fs,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName, const std::string& sysShift)
    : process_(process),
      region_(region),
      tauIdDiscriminators_(tauIdDiscriminators),
      tauIdName_(tauIdName),
      sysShift_(sysShift),
      selector_(0),
      cutFlow_(0),
      cutFlowTauAbsEtaBinning_(0),
      cutFlowTauPtBinning_(0),
      cutFlowSumEtBinning_(0),
      cutFlowVtxMultiplicityBinning_(0)
  {
    edm::ParameterSet cfgSelector;
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators_);
    cfgSelector.addParameter<std::string>("region", region_);

    selector_ = new TauIdEffEventSelector(cfgSelector);

//--- disable preselection cuts applied on tau-jet candidates
    selector_->tauLeadTrackPtMin_      =  -1.e+3; 
    selector_->tauAbsIsoMax_           =  +1.e+3;
    selector_->muTauPairAbsDzMax_      =  +1.e+3;
    selector_->muTauPairChargeProdMin_ =  -1.e+3;
    selector_->muTauPairChargeProdMax_ =  +1.e+3; 
    
    tauIdFlags_.resize(6 + tauIdDiscriminators_.size());
    
    edm::ParameterSet cfgCutFlowTable;
    cfgCutFlowTable.addParameter<std::string>("process", process_);
    cfgCutFlowTable.addParameter<std::string>("region", region_);
    cfgCutFlowTable.addParameter<std::string>("tauIdDiscriminator", tauIdName_);
    label_ = "presel";
    if ( sysShift != "CENTRAL_VALUE" ) label_.append("_").append(sysShift_);
    cfgCutFlowTable.addParameter<std::string>("label", label_);

    selectionNames_.resize(6 + tauIdDiscriminators_.size());
    selectionNames_[0] = std::string(region);
    selectionNames_[1] = "leadTrackFinding";
    selectionNames_[2] = "leadTrackPtCut";
    selectionNames_[3] = "loosePFIso";
    selectionNames_[4] = "eVeto";
    selectionNames_[5] = "muVeto";
    int idx = 6;
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator, ++idx ) {
      selectionNames_[idx] = tauIdDiscriminator->data();
    }
    cfgCutFlowTable.addParameter<vstring>("selectionNames", selectionNames_);

    float unbinned[] = { -0.5, +0.5 };
    edm::ParameterSet cfgUnbinned = cfgCutFlowTable;
    cfgUnbinned.addParameter<std::string>("binVariable", "");
    cutFlow_ = addCutFlowEntry(fs, "", 1, unbinned, cfgUnbinned);

    float tauAbsEtaBins[] = { 0.0, 1.5, 1.9, 2.3 }; 
    edm::ParameterSet cfgTauAbsEtaBinning = cfgCutFlowTable;
    cfgTauAbsEtaBinning.addParameter<std::string>("binVariable", "tauAbsEta");
    cutFlowTauAbsEtaBinning_ = addCutFlowEntry(fs, "tauAbsEtaBinning", 3, tauAbsEtaBins, cfgTauAbsEtaBinning);

    float tauPtBins[] = { 20.0, 25.0, 30.0, 40.0, 100.0 };
    edm::ParameterSet cfgTauPtBinning = cfgCutFlowTable;
    cfgTauPtBinning.addParameter<std::string>("binVariable", "tauPt");
    cutFlowTauPtBinning_ = addCutFlowEntry(fs, "tauPtBinning", 4, tauPtBins, cfgTauPtBinning);

    float sumEtBins[] = { 0.0, 250.0, 350.0, 450.0, 1000.0 };
    edm::ParameterSet cfgSumEtBinning = cfgCutFlowTable;
    cfgSumEtBinning.addParameter<std::string>("binVariable", "sumEt");
    cutFlowSumEtBinning_ = addCutFlowEntry(fs, "sumEtBinning", 4, sumEtBins, cfgSumEtBinning);

    float vtxMultiplicityBins[] = { -0.5, 4.5, 6.5, 8.5, 20.5 };
    edm::ParameterSet cfgVtxMultiplicityBinning = cfgCutFlowTable;
    cfgVtxMultiplicityBinning.addParameter<std::string>("binVariable", "sumEt");
    cutFlowVtxMultiplicityBinning_ = addCutFlowEntry(fs, "vtxMultiplicityBinning", 4, vtxMultiplicityBins, cfgVtxMultiplicityBinning);
  }
  ~regionEntryType()
  {
    delete selector_;

    delete cutFlow_;

    delete cutFlowTauAbsEtaBinning_;
    delete cutFlowTauPtBinning_;
    delete cutFlowSumEtBinning_;
    delete cutFlowVtxMultiplicityBinning_;
  }
  cutFlowEntryType* addCutFlowEntry(
    fwlite::TFileService& fs, const std::string& dirName, int numBins, float* binning, const edm::ParameterSet& cfg)
  {
    cutFlowEntryType* retVal = new cutFlowEntryType(cfg);

    if ( dirName != "" ) {
      TFileDirectory dir = fs.mkdir(dirName.data());
      retVal->bookCutFlowTables(dir, numBins, binning);
    } else {
      retVal->bookCutFlowTables(fs, numBins, binning);
    }

    return retVal;
  }
  void analyze(const PATMuTauPair& muTauPair, int genMatchType, size_t numVertices, double evtWeight)
  {
    //std::cout << "<cutFlowEntryType::analyze>:" << std::endl;

    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, evtSelFlags) ) {

//--- set flags indicating whether tau-jet candidate passes 
//    "leading" track finding, leading track Pt and loose (PF)isolation requirements 
//    plus tau id. discriminators
      tauIdFlags_[0] = true;
      tauIdFlags_[1] = (muTauPair.leg2()->userFloat("hasLeadTrack")       > 0.5);
      tauIdFlags_[2] = (muTauPair.leg2()->userFloat("leadTrackPt")        > 5.0);      
      tauIdFlags_[3] = (muTauPair.leg2()->userFloat("preselLoosePFIsoPt") < 2.5);
      tauIdFlags_[4] = (muTauPair.leg2()->userFloat("PFElectronMVA")      < 0.6);
      tauIdFlags_[5] = (muTauPair.leg2()->userFloat("dRnearestMuon")      > 0.5);
      int idx = 6;
      for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	    tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator, ++idx ) {
	//std::cout << " tauIdDiscriminator = " << (*tauIdDiscriminator) << ":" 
	//	    << " " << muTauPair.leg2()->tauID(tauIdDiscriminator->data()) << std::endl;
	tauIdFlags_[idx] = (muTauPair.leg2()->tauID(tauIdDiscriminator->data()) > 0.5);
      }
      //std::cout << "tauIdFlags = " << format_vbool(tauIdFlags_) << std::endl;

//--- fill histograms for "inclusive" tau id. efficiency measurement
      cutFlow_->fillCutFlowTables(0., tauIdFlags_, genMatchType, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of tau-jet pseudo-rapidity
      double tauAbsEta = TMath::Abs(muTauPair.leg1()->eta());
      cutFlowTauAbsEtaBinning_->fillCutFlowTables(tauAbsEta, tauIdFlags_, genMatchType, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of tau-jet transverse momentum
      double tauPt = muTauPair.leg1()->pt();
      cutFlowTauPtBinning_->fillCutFlowTables(tauPt, tauIdFlags_, genMatchType, evtWeight);
      
//--- fill histograms for tau id. efficiency measurement as function of sumEt
      double sumEt = muTauPair.met()->sumEt();
      cutFlowSumEtBinning_->fillCutFlowTables(sumEt, tauIdFlags_, genMatchType, evtWeight);
      
//--- fill histograms for tau id. efficiency measurement as function of reconstructed vertex multiplicity
      cutFlowVtxMultiplicityBinning_->fillCutFlowTables(numVertices, tauIdFlags_, genMatchType, evtWeight);
    }
  }

  std::string process_;
  std::string region_;
  vstring tauIdDiscriminators_;
  std::string tauIdName_;
  std::string sysShift_;
  std::string label_;
  vstring selectionNames_;
  
  TauIdEffEventSelector* selector_;

  std::vector<bool> tauIdFlags_;
  
  cutFlowEntryType* cutFlow_;

  cutFlowEntryType* cutFlowTauAbsEtaBinning_;
  cutFlowEntryType* cutFlowTauPtBinning_;
  cutFlowEntryType* cutFlowSumEtBinning_;
  cutFlowEntryType* cutFlowVtxMultiplicityBinning_;
};

int getGenMatchType(const PATMuTauPair& muTauPair, const reco::GenParticleCollection& genParticles)
{
//--- check if reconstructed tau-jet candidate matches "true" hadronic tau decay on generator level,
//    is a "fake" tau (i.e. matches a quark/gluon/e/mu/photon on generator level)
//    or fails to be matched to any generator level object
//
//    NOTE: code to perform matching taken from TauAnalysis/Core/plugins/TauHistManager.cc
//
  int absMatchingGenParticlePdgId           = TMath::Abs(getMatchingGenParticlePdgId(muTauPair.leg2()->p4(), genParticles, 0, true));
  int absMatchingFinalStateGenParticlePdgId = TMath::Abs(getMatchingGenParticlePdgId(muTauPair.leg2()->p4(), genParticles, 0, false));

  std::string genTauDecayMode = "";
  if ( muTauPair.leg2()->genJet() != 0 ) {
    genTauDecayMode = JetMCTagUtils::genTauDecayMode(*muTauPair.leg2()->genJet());
  } else if ( absMatchingGenParticlePdgId == 15 ) { // special handling of tau --> electron/muon decays
    if      ( absMatchingFinalStateGenParticlePdgId == 11 ) genTauDecayMode = "electron";
    else if ( absMatchingFinalStateGenParticlePdgId == 13 ) genTauDecayMode = "muon";
  }

  if ( absMatchingGenParticlePdgId == 15 &&
       (genTauDecayMode == "oneProng0Pi0"    ||
	genTauDecayMode == "oneProng1Pi0"    ||
	genTauDecayMode == "oneProng2Pi0"    ||
	genTauDecayMode == "oneProngOther"   ||
	genTauDecayMode == "threeProng0Pi0"  ||
	genTauDecayMode == "threeProngOther" ||
	genTauDecayMode == "threeProngOther" ||
	genTauDecayMode == "rare"            ) ) return kTauHadMatched;
  else if ( (absMatchingGenParticlePdgId >= 1 && absMatchingGenParticlePdgId <= 6) ||
	     absMatchingGenParticlePdgId == 11 || absMatchingGenParticlePdgId == 13 ||
	     absMatchingGenParticlePdgId == 21 ||
	     absMatchingGenParticlePdgId == 22 ) return kFakeTauMatched;
  else                                           return kUnmatched;
}

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauIdEffPreselNumbers>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauIdEffPreselNumbers");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteTauIdEffPreselNumbers") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauIdEffPreselNumbers = cfg.getParameter<edm::ParameterSet>("tauIdEffPreselNumbers");

  edm::InputTag srcMuTauPairs = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcMuTauPairs");
  edm::InputTag srcGenParticles = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcGenParticles");
  edm::InputTag srcTrigger = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcTrigger");
  vstring hltPaths = cfgTauIdEffPreselNumbers.getParameter<vstring>("hltPaths");
  edm::InputTag srcGoodMuons = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcGoodMuons");
  edm::InputTag srcVertices = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcVertices");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffPreselNumbers.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffPreselNumbers.exists("sysShift") ?
    cfgTauIdEffPreselNumbers.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauIdEffPreselNumbers.getParameter<std::string>("process");
  vstring regions = cfgTauIdEffPreselNumbers.getParameter<vstring>("regions");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauIdEffPreselNumbers.getParameter<vParameterSet>("tauIds");
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

  int numEvents_processed = 0; 
  double numEventsWeighted_processed = 0;
  
  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauIdEffPreselNumbers") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << " opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

      //std::cout << "processing event " << evt.id().event() << ", run " << evt.id().run() << std::endl;

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
      //std::cout << "isTriggered = " << isTriggered << std::endl;

//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
      double evtWeight = 1.0;
      for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	    srcWeight != srcWeights.end(); ++srcWeight ) {
	edm::Handle<double> weight;
	evt.getByLabel(*srcWeight, weight);
	evtWeight *= (*weight);
      }
      //std::cout << "evtWeight = " << evtWeight << std::endl;

//--- require event to pass trigger requirements
//    and to contain only one "good quality" muon
      typedef std::vector<pat::Muon> PATMuonCollection;
      edm::Handle<PATMuonCollection> goodMuons;
      evt.getByLabel(srcGoodMuons, goodMuons);
      size_t numGoodMuons = goodMuons->size();
      //std::cout << "numGoodMuons = " << numGoodMuons << std::endl;

      if ( isTriggered && numGoodMuons <= 1 ) {

	edm::Handle<reco::VertexCollection> vertices;
	evt.getByLabel(srcVertices, vertices);
	size_t numVertices = vertices->size();

//--- iterate over collection of muon + tau-jet pairs
	edm::Handle<PATMuTauPairCollection> muTauPairs;
	evt.getByLabel(srcMuTauPairs, muTauPairs);

	edm::Handle<reco::GenParticleCollection> genParticles;
	evt.getByLabel(srcGenParticles, genParticles);

	int muTauPairIdx = 0;
	for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	      muTauPair != muTauPairs->end(); ++muTauPair, ++muTauPairIdx ) {
	  int genMatchType = getGenMatchType(*muTauPair, *genParticles);
	  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
		regionEntry != regionEntries.end(); ++regionEntry ) {   
/* 
	    if ( (*regionEntry)->region_ == "C1" ) {
	      pat::strbitset evtSelFlags;
	      if ( genMatchType == kTauHadMatched || (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	    	std::cout << "muTauPair #" << muTauPairIdx << std::endl;
	    	std::cout << " leg1: Pt = " << muTauPair->leg1()->pt() << "," 
	    		  << " eta = " << muTauPair->leg1()->eta() << ", phi = " << muTauPair->leg1()->phi() << std::endl;
	    	std::cout << " leg2: Pt = " << muTauPair->leg2()->pt() << "," 
	    		  << " eta = " << muTauPair->leg2()->eta() << ", phi = " << muTauPair->leg2()->phi();
	    	if      ( genMatchType == kTauHadMatched  ) std::cout << " ('true' hadronic tau decay)";
	    	else if ( genMatchType == kFakeTauMatched ) std::cout << " (tau fake)"; 
	    	else                                        std::cout << " (no gen. match)"; 
	    	std::cout << std::endl;
	      }
	    }
 */
	    (*regionEntry)->analyze(*muTauPair, genMatchType, numVertices, evtWeight);
	  }
	}
      }

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;
    }

//--- close input file
    delete inputFile;
  }

  std::cout << "<FWLiteTauIdEffPreselNumbers>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;
  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    std::cout << " region " << (*regionEntry)->region_ << ", " << (*regionEntry)->tauIdName_ << std::endl;
    TauIdEffCutFlowTable* cutFlowTableTauHadMatched = (*regionEntry)->cutFlow_->cutFlowTauHadMatched_;
    double effPreselection = cutFlowTableTauHadMatched->getCutFlowNumber(0, 5)/
                             cutFlowTableTauHadMatched->getCutFlowNumber(0, 0);
    std::cout << "  eff(preselection) = " << effPreselection << std::endl;
    double effTauId = cutFlowTableTauHadMatched->getCutFlowNumber(0, 5 + (*regionEntry)->tauIdDiscriminators_.size())/
                      cutFlowTableTauHadMatched->getCutFlowNumber(0, 5);
    std::cout << "  eff(tauId) = " << effTauId << std::endl;
    TauIdEffCutFlowTable* cutFlowTableFakeTauMatched = (*regionEntry)->cutFlow_->cutFlowFakeTauMatched_;
    double purity = cutFlowTableTauHadMatched->getCutFlowNumber(0, 5)/
                   (cutFlowTableTauHadMatched->getCutFlowNumber(0, 5) + cutFlowTableFakeTauMatched->getCutFlowNumber(0, 5));
    std::cout << "  purity = " << purity << std::endl;
  }

  clock.Show("FWLiteTauIdEffPreselNumbers");

  return 0;
}
