
/** \executable FWLiteMuonIsolationAnalyzer
 *
 * Study probability for muons in QCD background events to pass tight muon isolation,
 * provided muon passes loose muon isolation applied on trigger level
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: FWLiteMuonIsolationAnalyzer.cc,v 1.1 2011/10/19 14:37:00 veelken Exp $
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/TauIdEfficiency/interface/MuonIsolationHistManager.h"
#include "TauAnalysis/RecoTools/interface/PATObjectLUTvalueExtractorFromKNN.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

#include <vector>
#include <string>
#include <fstream>

typedef std::vector<double> vdouble;
typedef std::vector<std::string> vstring;

template <typename T>
double getUserFloat(const T& lepton, const std::string& userFloatName)
{ 
  if ( lepton.hasUserFloat(userFloatName) ) {
    return lepton.userFloat(userFloatName);
  } else {
    std::cerr << "Error in <getUserFloat>: no userFloat with name = " << userFloatName << " stored in lepton !!" << std::endl;
    std::cerr << " available userFloats = " << format_vstring(lepton.userFloatNames()) << std::endl;
    return -1.;
  }
}

struct histManagerEntryType
{
  histManagerEntryType(const edm::ParameterSet& cfg, 
		       const std::string& triggerPath, double muonIsoThreshold_loose, double muonIsoThreshold_tight)
    : triggerPath_(triggerPath),
      muonIsoThreshold_loose_(muonIsoThreshold_loose),
      muonIsoThreshold_tight_(muonIsoThreshold_tight),
      histManager_all_weighted_(0),
      muonIsoProbExtractor_(0),
      ntuple_(0)
  {
    histManager_all_ = new MuonIsolationHistManager(cfg);
    histManager_passed_ = new MuonIsolationHistManager(cfg);
    histManager_failed_ = new MuonIsolationHistManager(cfg);

    if ( cfg.exists("muonIsoProbExtractor") ) {
      edm::ParameterSet cfgMuonIsoProbExtractor = cfg.getParameter<edm::ParameterSet>("muonIsoProbExtractor");
      muonIsoProbExtractor_ = new PATMuonLUTvalueExtractorFromKNN(cfgMuonIsoProbExtractor);
      histManager_all_weighted_ = new MuonIsolationHistManager(cfg);
    }
  }
  ~histManagerEntryType() 
  {
    delete muonIsoProbExtractor_;
  }
  void bookHistograms(TFileDirectory& dir)
  {
    std::string subdirectory = Form("%s_loose%02.0f_tight%02.0f", triggerPath_.data(), muonIsoThreshold_loose_*10., muonIsoThreshold_tight_*10.);
    TFileDirectory subdir = ( subdirectory != "" ) ? dir.mkdir(subdirectory.data()) : dir;
    TFileDirectory subdir_all = subdir.mkdir("all");
    histManager_all_->bookHistograms(subdir_all);
    TFileDirectory subdir_passed = subdir.mkdir("passed");
    histManager_passed_->bookHistograms(subdir_passed);
    TFileDirectory subdir_failed = subdir.mkdir("failed");
    histManager_failed_->bookHistograms(subdir_failed);

    if ( histManager_all_weighted_ ) {
      TFileDirectory subdir_all_weighted = subdir.mkdir("all_weighted");
      histManager_all_weighted_->bookHistograms(subdir_all_weighted);
    }
    
    ntuple_ = subdir_all.make<TTree>("ntuple", "Muon Isolation Ntuple");
    ntuple_->Branch("muonPt",  &muonPt_,  "muonPt/F");
    ntuple_->Branch("muonEta", &muonEta_, "muonEta/F");
    ntuple_->Branch("muonIso", &muonIso_, "muonIso/F");
    ntuple_->Branch("sumEt",   &sumEt_,   "sumEt/F");
    ntuple_->Branch("weight",  &weight_,  "weight/F");
  }
  void fillHistograms(const PATMuTauPair& muTauPair, size_t numVertices, double weight)
  {
    const pat::Muon& muon = *muTauPair.leg1();
    if ( getUserFloat(muon, triggerPath_) > 0.5 ) {
      double muonPt = muon.pt();
      // compute deltaBeta corrected isolation Pt sum "by hand":
      //   muonIsoPtSum = pfChargedParticles(noPileUp) + pfNeutralHadrons + pfGammas - deltaBetaCorr, deltaBetaCorr = 0.5*pfChargedParticlesPileUp
      // ( User1Iso = pfAllChargedHadrons(noPileUp), User2Iso = pfAllChargedHadronsPileUp
      //   as defined in TauAnalysis/TauIdEfficiency/test/commissioning/produceMuonIsolationPATtuple_cfg.py )
      double muonIsoPtSum1 = muon.userIsolation(pat::User1Iso) 
	                    + TMath::Max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) 
                                            - 0.5*muon.userIsolation(pat::User2Iso));
      //double muonIsoPtSum2 = getUserFloat(muon, "pfLooseIsoPt04");
      //std::cout << "muonIsoPtSum1 = " << muonIsoPtSum1 << ", muonIsoPtSum2 = " << muonIsoPtSum2 << std::endl;
      if ( muonIsoPtSum1 < (muonIsoThreshold_loose_*muonPt) ) {
	histManager_all_->fillHistograms(muTauPair, numVertices, muonIsoPtSum1, weight);
	if ( muonIsoPtSum1 < (muonIsoThreshold_tight_*muonPt) ) {
	  histManager_passed_->fillHistograms(muTauPair, numVertices, muonIsoPtSum1, weight);
	} else {
	  histManager_failed_->fillHistograms(muTauPair, numVertices, muonIsoPtSum1, weight);
	}

	if ( histManager_all_weighted_ && muonIsoProbExtractor_ ) {
	  double muonIsoProbValue = (*muonIsoProbExtractor_)(muon);
	  histManager_all_weighted_->fillHistograms(muTauPair, numVertices, muonIsoPtSum1, weight*muonIsoProbValue);
	}

	muonPt_ = muonPt;
	muonEta_ = muon.eta();
	muonIso_ = muonIsoPtSum1;
	sumEt_ = muTauPair.met()->sumEt() - (muTauPair.leg1()->pt() + muTauPair.leg2()->pt());
	weight_ = weight;
	ntuple_->Fill();
      }
    }
  }

  std::string triggerPath_;

  double muonIsoThreshold_loose_;
  double muonIsoThreshold_tight_;

  MuonIsolationHistManager* histManager_all_;
  MuonIsolationHistManager* histManager_passed_;
  MuonIsolationHistManager* histManager_failed_;

  MuonIsolationHistManager* histManager_all_weighted_;
  PATMuonLUTvalueExtractorFromKNN* muonIsoProbExtractor_;

  TTree* ntuple_;
  Float_t muonPt_;
  Float_t muonEta_;
  Float_t muonIso_;
  Float_t sumEt_;
  Float_t weight_;
};

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteMuonIsolationAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteMuonIsolationAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteMuonIsolationAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgMuonIsolationAnalyzer = cfg.getParameter<edm::ParameterSet>("muonIsolationAnalyzer");

  edm::InputTag srcMuonsTightId = cfgMuonIsolationAnalyzer.getParameter<edm::InputTag>("srcMuonsTightId");
  edm::InputTag srcMuonsLooseId = cfgMuonIsolationAnalyzer.getParameter<edm::InputTag>("srcMuonsLooseId");
  edm::InputTag srcTauJetCandidates = cfgMuonIsolationAnalyzer.getParameter<edm::InputTag>("srcTauJetCandidates");
  edm::InputTag srcMuTauPairs = cfgMuonIsolationAnalyzer.getParameter<edm::InputTag>("srcMuTauPairs");
  edm::InputTag srcVertices = cfgMuonIsolationAnalyzer.getParameter<edm::InputTag>("srcVertices");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgMuonIsolationAnalyzer.getParameter<vInputTag>("weights");

  std::string directory = cfgMuonIsolationAnalyzer.getParameter<std::string>("directory");
  
  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- book histograms
  TFileDirectory dir = ( directory != "" ) ? fs.mkdir(directory) : fs;
  edm::ParameterSet cfgMuonIsolationHistManager;
  if ( cfgMuonIsolationAnalyzer.exists("muonIsoProbExtractor") ) {
    edm::ParameterSet cfgMuonIsoProbExtractor = cfgMuonIsolationAnalyzer.getParameter<edm::ParameterSet>("muonIsoProbExtractor");
    cfgMuonIsolationHistManager.addParameter<edm::ParameterSet>("muonIsoProbExtractor", cfgMuonIsoProbExtractor);
  }
  std::vector<histManagerEntryType*> histManagerEntries;
  vstring triggerPaths = cfgMuonIsolationAnalyzer.getParameter<vstring>("triggerPaths");
  vdouble muonIsoThresholds_loose = cfgMuonIsolationAnalyzer.getParameter<vdouble>("muonIsoThresholdsLoose");
  vdouble muonIsoThresholds_tight = cfgMuonIsolationAnalyzer.getParameter<vdouble>("muonIsoThresholdsTight");
  for ( vstring::const_iterator triggerPath = triggerPaths.begin();
	triggerPath != triggerPaths.end(); ++triggerPath ) {
    for ( vdouble::const_iterator muonIsoThreshold_loose = muonIsoThresholds_loose.begin();
	  muonIsoThreshold_loose != muonIsoThresholds_loose.end(); ++muonIsoThreshold_loose ) {
      for ( vdouble::const_iterator muonIsoThreshold_tight = muonIsoThresholds_tight.begin();
	    muonIsoThreshold_tight != muonIsoThresholds_tight.end(); ++muonIsoThreshold_tight ) {
	histManagerEntryType* histManagerEntry = 
	  new histManagerEntryType(cfgMuonIsolationHistManager, *triggerPath, *muonIsoThreshold_loose, *muonIsoThreshold_tight);
	histManagerEntry->bookHistograms(dir);
	histManagerEntries.push_back(histManagerEntry);
      }
    }
  }

  int    numEvents_processed                = 0; 
  double numEventsWeighted_processed        = 0.;
  int    numEvents_passedPresel             = 0;
  double numEventsWeighted_passedPresel     = 0.;
  int    numEvents_passedDiMuonVeto         = 0;
  double numEventsWeighted_passedDiMuonVeto = 0.;
  int    numEvents_passedMuTauPair          = 0;
  double numEventsWeighted_passedMuTauPair  = 0.;

  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteMuonIsolationAnalyzer") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << "opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

      // CV: due to problem with EDFilter configuration during PAT-tuple production,
      //     not all objects are available for each event --> test presence and skip event processing in case objects are not available
      //    (use mu + tau-jet pair as "test" object)
      edm::Handle<PATMuTauPairCollection> testObject;
      evt.getByLabel(srcMuTauPairs, testObject);
      if ( !testObject.isValid() ) continue;
      ++numEvents_passedPresel;
      //numEventsWeighted_passedPresel += evtWeight;

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
      
      edm::Handle<pat::MuonCollection> muonsLooseIdSel;
      evt.getByLabel(srcMuonsLooseId, muonsLooseIdSel);
      if ( muonsLooseIdSel->size() >= 2 ) continue;
      ++numEvents_passedDiMuonVeto;
      numEventsWeighted_passedDiMuonVeto += evtWeight;

      edm::Handle<pat::MuonCollection> muonsTightIdSel;
      evt.getByLabel(srcMuonsTightId, muonsTightIdSel);
      
      edm::Handle<pat::TauCollection> tauJetCandidates;
      evt.getByLabel(srcTauJetCandidates, tauJetCandidates);
      
      edm::Handle<PATMuTauPairCollection> muTauPairs;
      evt.getByLabel(srcMuTauPairs, muTauPairs);

      const PATMuTauPair* bestMuTauPair = 0;
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {
	if ( muTauPair->leg2()->pfJetRef()->pt() > 20. && TMath::Abs(muTauPair->leg2()->pfJetRef()->eta()) < 2.3 &&
	     TMath::Abs(muTauPair->leg1()->vertex().z() - muTauPair->leg2()->vertex().z()) < 0.2 ) {
	  if ( !bestMuTauPair ) bestMuTauPair = &(*muTauPair); // CV: simply take first object passing selection criteria for now...
	}
      }
      if ( !bestMuTauPair ) continue;
      ++numEvents_passedMuTauPair;
      numEventsWeighted_passedMuTauPair += evtWeight;

//--- determine number of vertices reconstructed in the event
//   (needed to parametrize dependency of tau id. efficiency on number of pile-up interactions)
      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      size_t numVertices = vertices->size();

      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = histManagerEntries.begin();
	    histManagerEntry != histManagerEntries.end(); ++histManagerEntry ) {
	(*histManagerEntry)->fillHistograms(*bestMuTauPair, numVertices, evtWeight);
      }
    }

//--- close input file
    delete inputFile;
  }

  for ( std::vector<histManagerEntryType*>::iterator it = histManagerEntries.begin();
	it != histManagerEntries.end(); ++it ) {  
    delete (*it);
  }

  std::cout << "<FWLiteMuonIsolationAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;
  std::cout << " numEvents_passedPresel: " << numEvents_passedPresel
	    << " (weighted = " << numEventsWeighted_passedPresel << ")" << std::endl;
  std::cout << " numEvents_passedDiMuonVeto: " << numEvents_passedDiMuonVeto 
	    << " (weighted = " << numEventsWeighted_passedDiMuonVeto << ")" << std::endl;
  std::cout << " numEvents_passedMuTauPair: " << numEvents_passedMuTauPair 
	    << " (weighted = " << numEventsWeighted_passedMuTauPair << ")" << std::endl;
  
  clock.Show("FWLiteMuonIsolationAnalyzer");

  return 0;
}
