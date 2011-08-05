
/** \executable FWLiteTauPtResAnalyzer
 *
 * Study reconstructed tauPt vs. Pt of generated tau decay products
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: FWLiteTauPtResAnalyzer.cc,v 1.4 2011/08/02 15:16:09 veelken Exp $
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauPtResHistManager.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

#include <vector>
#include <string>

typedef std::vector<std::string> vstring;

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
    throw cms::Exception("FWLiteTauFakeRateAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauFakeRateAnalyzer = cfg.getParameter<edm::ParameterSet>("tauPtResAnalyzer");

  edm::InputTag srcTauJetCandidates = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcTauJetCandidates");
  edm::InputTag srcGenParticles = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcGenParticles");

  std::string directory = cfgTauFakeRateAnalyzer.getParameter<std::string>("directory");
  
  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- book histograms
  edm::ParameterSet cfgTauPtResHistManager;
  TauPtResHistManager histManager(cfgTauPtResHistManager);
  TFileDirectory dir = ( directory != "" ) ? fs.mkdir(directory) : fs;
  histManager.bookHistograms(dir);

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

      //std::cout << "processing run = " << evt.id().run() << ":" 
      //	  << " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << std::endl;

      double evtWeight = 1.0; // vertex multiplicity reweighting not yet implemented...

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

      edm::Handle<pat::TauCollection> tauJetCandidates;
      evt.getByLabel(srcTauJetCandidates, tauJetCandidates);

      edm::Handle<reco::GenParticleCollection> genParticles;
      evt.getByLabel(srcGenParticles, genParticles);

//--- iterate over collection of tau-jet candidates:
//    check if tau-jet candidate passed/fails tau id. criteria,
//    fill histograms for selected tau-jet candidates
      for ( pat::TauCollection::const_iterator tauJetCand = tauJetCandidates->begin();
	    tauJetCand != tauJetCandidates->end(); ++tauJetCand ) {
	if ( tauJetCand->genJet() &&
	     tauJetCand->genJet()->pt() > 20. && tauJetCand->genJet()->pt() < 40. &&
	     TMath::Abs(tauJetCand->genJet()->eta()) < 2.3 &&
	     tauJetCand->tauID("decayModeFinding")     > 0.5 &&
	     tauJetCand->tauID("byLooseIsolation")     > 0.5 &&
	     tauJetCand->tauID("againstElectronLoose") > 0.5 &&
	     tauJetCand->tauID("againstMuonTight")     > 0.5 ) {
	  histManager.fillHistograms(*tauJetCand, *genParticles, evtWeight);
	}
      }
    }

//--- close input file
    delete inputFile;
  }

  std::cout << "<FWLiteTauPtResAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;

  clock.Show("FWLiteTauPtResAnalyzer");

  return 0;
}
