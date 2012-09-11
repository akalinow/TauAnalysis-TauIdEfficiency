
/** \executable FWLiteTauPtResAnalyzer
 *
 * Study reconstructed tauPt vs. Pt of generated tau decay products
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.6 $
 *
 * $Id: FWLiteTauPtResAnalyzer.cc,v 1.6 2011/09/30 12:26:40 veelken Exp $
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
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/View.h"

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauPtResHistManager.h"
#include "TauAnalysis/TauIdEfficiency/interface/tauPtResAuxFunctions.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

#include <vector>
#include <string>
#include <fstream>

typedef std::vector<std::string> vstring;

TauPtResHistManager* addHistManager(const edm::ParameterSet& cfg, TFileDirectory& dir, const std::string& subdirectory)
{
  TFileDirectory subdir = ( subdirectory != "" ) ? dir.mkdir(subdirectory) : dir;
  TauPtResHistManager* histManager = new TauPtResHistManager(cfg);
  histManager->bookHistograms(subdir);
  return histManager;
}

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauPtResAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauPtResAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteTauPtResAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
 
  edm::ParameterSet cfgTauPtResAnalyzer = cfg.getParameter<edm::ParameterSet>("tauPtResAnalyzer");

  edm::InputTag srcTauJetCandidates = cfgTauPtResAnalyzer.getParameter<edm::InputTag>("srcTauJetCandidates");
  edm::InputTag srcGenParticles = cfgTauPtResAnalyzer.getParameter<edm::InputTag>("srcGenParticles");
  edm::InputTag srcTracks = cfgTauPtResAnalyzer.getParameter<edm::InputTag>("srcTracks");
  edm::InputTag srcVertices = cfgTauPtResAnalyzer.getParameter<edm::InputTag>("srcVertices");
  edm::InputTag srcCaloTowers = cfgTauPtResAnalyzer.getParameter<edm::InputTag>("srcCaloTowers");

  edm::ParameterSet cfgQualityCuts = cfgTauPtResAnalyzer.getParameter<edm::ParameterSet>("qualityCuts");
  reco::tau::RecoTauQualityCuts qualityCuts(cfgQualityCuts);

  std::string directory = cfgTauPtResAnalyzer.getParameter<std::string>("directory");
  
  std::string selEventsFileName = cfgTauPtResAnalyzer.getParameter<std::string>("selEventsFileName");

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- book histograms
  edm::ParameterSet cfgTauPtResHistManager;
  TFileDirectory dir = ( directory != "" ) ? fs.mkdir(directory) : fs;
  TauPtResHistManager* histManager                      = addHistManager(cfgTauPtResHistManager, dir, "");
  TauPtResHistManager* histManagerVtxMultiplicityLe4    = addHistManager(cfgTauPtResHistManager, dir, "vtxMultiplicityLe4");
  TauPtResHistManager* histManagerVtxMultiplicity5to12  = addHistManager(cfgTauPtResHistManager, dir, "vtxMultiplicity5to12");
  TauPtResHistManager* histManagerVtxMultiplicity13to20 = addHistManager(cfgTauPtResHistManager, dir, "vtxMultiplicity13to20");
  TauPtResHistManager* histManagerVtxMultiplicityGe21   = addHistManager(cfgTauPtResHistManager, dir, "vtxMultiplicityGe21");
  
  std::ofstream* selEventsFile = new std::ofstream(selEventsFileName.data(), std::ios::out);

  int    numEvents_processed         = 0; 
  double numEventsWeighted_processed = 0.;

  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauFakeRateAnalyzer") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << "opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

      std::cout << "processing run = " << evt.id().run() << ":" 
		<< " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << std::endl;

      double evtWeight = 1.0; // vertex multiplicity reweighting not yet implemented...

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      if ( !(vertices->size() >= 1) ) continue;
      size_t vtxMultiplicity = vertices->size();
      const reco::Vertex& theEventVertex = vertices->at(0);
      
      qualityCuts.setPV(reco::VertexRef(vertices, 0));

      edm::Handle<pat::TauCollection> tauJetCandidates;
      evt.getByLabel(srcTauJetCandidates, tauJetCandidates);

      edm::Handle<reco::TrackCollection> tracks;
      evt.getByLabel(srcTracks, tracks);

      edm::Handle<CaloTowerCollection> caloTowers;
      evt.getByLabel(srcCaloTowers, caloTowers);

      edm::Handle<reco::GenParticleCollection> genParticles;
      evt.getByLabel(srcGenParticles, genParticles);

//--- iterate over collection of tau-jet candidates:
//    check if tau-jet candidate passed/fails tau id. criteria,
//    fill histograms for selected tau-jet candidates
      for ( pat::TauCollection::const_iterator tauJetCand = tauJetCandidates->begin();
	    tauJetCand != tauJetCandidates->end(); ++tauJetCand ) {
	if ( tauJetCand->genJet() &&
	     TMath::Abs(tauJetCand->genJet()->eta()) < 2.3 &&
	     tauJetCand->tauID("decayModeFinding") > 0.5 &&
	     tauJetCand->tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 &&
	     tauJetCand->tauID("againstElectronLoose") > 0.5 &&
	     tauJetCand->tauID("againstMuonMedium") > 0.5 ) {	  
	  std::string genTauDecayMode = getGenTauDecayMode(*tauJetCand, *genParticles);
	  if (  (genTauDecayMode == "oneProng0Pi0" ||
		 genTauDecayMode == "oneProng1Pi0" ||
		 genTauDecayMode == "oneProng2Pi0" ||
		 genTauDecayMode == "threeProng0Pi0") &&
		tauJetCand->genJet()  ) {
	    if ( tauJetCand->pt() < (0.5*tauJetCand->genJet()->pt()) ) {
	      std::cout << "run = " << evt.id().run() << "," 
			<< " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << ":" << std::endl;
	      printPatTau(*tauJetCand);
	      std::cout << "tauPt = " << tauJetCand->pt() << std::endl;
	      std::cout << std::endl;
	      printRecoPFJet(*tauJetCand->pfJetRef(), theEventVertex);
	      std::cout << std::endl;

	      std::cout << "<printTracks>:" << std::endl;
	      size_t numTracks = tracks->size();
	      for ( size_t iTrack = 0; iTrack < numTracks; ++iTrack ) {
		reco::TrackRef track(tracks, iTrack);
		double dR = deltaR(track->eta(), track->phi(), tauJetCand->eta(), tauJetCand->phi());
		if ( dR < 0.5 && track->pt() > 2. ) printTrack(track, theEventVertex);
	      }
	      std::cout << std::endl;
	      
	      printCaloTowers(*caloTowers, tauJetCand->p4(), 0.5);
	      std::cout << std::endl;

	      (*selEventsFile) << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() << std::endl;
	    }
	  }

	  histManager->fillHistograms(*tauJetCand, *genParticles, theEventVertex, qualityCuts, evtWeight);

	  if      ( vtxMultiplicity <=  4 ) 
	    histManagerVtxMultiplicityLe4->fillHistograms(*tauJetCand, *genParticles, theEventVertex, qualityCuts, evtWeight);
	  else if ( vtxMultiplicity <= 12 ) 
	    histManagerVtxMultiplicity5to12->fillHistograms(*tauJetCand, *genParticles, theEventVertex, qualityCuts, evtWeight);
	  else if ( vtxMultiplicity <= 20 ) 
	    histManagerVtxMultiplicity13to20->fillHistograms(*tauJetCand, *genParticles, theEventVertex, qualityCuts, evtWeight);
	  else 
	    histManagerVtxMultiplicityGe21->fillHistograms(*tauJetCand, *genParticles, theEventVertex, qualityCuts, evtWeight);	  
	}
      }
    }

//--- close input file
    delete inputFile;
  }

  delete histManager;
  delete histManagerVtxMultiplicityLe4;
  delete histManagerVtxMultiplicity5to12;
  delete histManagerVtxMultiplicity13to20;
  delete histManagerVtxMultiplicityGe21;

//--- close ASCII file containing 
//     run:lumi-section:event 
//    numbers of events with PtRes < 0.5
  delete selEventsFile;

  std::cout << "<FWLiteTauPtResAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;

  clock.Show("FWLiteTauPtResAnalyzer");

  return 0;
}
