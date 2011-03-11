#include "TauAnalysis/TauIdEfficiency/plugins/PATTauCollectionDiffAnalyzer.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

PATTauCollectionDiffAnalyzer::PATTauCollectionDiffAnalyzer(const edm::ParameterSet& cfg)
  : dqmError_(0)
{
  //std::cout << "<PATTauCollectionDiffAnalyzer::PATTauCollectionDiffAnalyzer>:" << std::endl;

  patTauSource1_ = cfg.getParameter<edm::InputTag>("patTauSource1");
  patTauSource2_ = cfg.getParameter<edm::InputTag>("patTauSource2");

  patTauSelector_ = ( cfg.exists("patTauSelection") ) ? 
    new StringCutObjectSelector<pat::Tau>(cfg.getParameter<std::string>("patTauSelection")) : 0;

  dRmatch_ = cfg.getParameter<double>("dRmatch");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

PATTauCollectionDiffAnalyzer::~PATTauCollectionDiffAnalyzer()
{
  delete patTauSelector_;
}

void PATTauCollectionDiffAnalyzer::beginJob()
{
  //std::cout << "<PATTauCollectionDiffAnalyzer::beginJob>:" << std::endl;

  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("PATTauCollectionDiffAnalyzer::beginJob") << " Failed to access dqmStore --> histograms will NOT be booked !!";
    dqmError_ = 1;
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory_);

  histogramsPatTaus1_.jetPt_     = dqmStore.book1D("jetPt1",     "jetPt1",     50,  0.,  100.);
  histogramsPatTaus1_.jetEta_    = dqmStore.book1D("jetEta1",    "jetEta1",    50, -2.5,  +2.5);
  histogramsPatTausSel1_.jetPt_  = dqmStore.book1D("selJetPt1",  "selJetPt1",  50,  0.,  100.);
  histogramsPatTausSel1_.jetEta_ = dqmStore.book1D("selJetEta1", "selJetEta1", 50, -2.5,  +2.5);
  
  histogramsPatTaus2_.jetPt_     = dqmStore.book1D("jetPt2",     "jetPt2",     50,  0.,  100.);
  histogramsPatTaus2_.jetEta_    = dqmStore.book1D("jetEta2",    "jetEta2",    50, -2.5,  +2.5);
  histogramsPatTausSel2_.jetPt_  = dqmStore.book1D("selJetPt2",  "selJetPt2",  50,  0.,  100.);
  histogramsPatTausSel2_.jetEta_ = dqmStore.book1D("selJetEta2", "selJetEta2", 50, -2.5,  +2.5);
}

reco::Candidate::LorentzVector getJetP4(const pat::Tau& patTau)
{
  if ( patTau.isPFTau() ) {
    return patTau.pfJetRef()->p4();
  } else if ( patTau.isCaloTau() ) {
    return patTau.caloTauTagInfoRef()->jetRef()->p4();
  } else {
    edm::LogError ("getJetP4") << " Undefined pat::Tau type !!";
    return reco::Candidate::LorentzVector(0,0,0,0);
  }
}

void fillHistograms(MonitorElement* hTauJetPt, MonitorElement* hTauJetEta, const pat::TauCollection& patTaus)
{
  for ( pat::TauCollection::const_iterator patTau = patTaus.begin();
	patTau != patTaus.end(); ++patTau ) {
    reco::Candidate::LorentzVector patTauJetP4 = getJetP4(*patTau);

    hTauJetPt->Fill(patTauJetP4.pt());
    hTauJetEta->Fill(patTauJetP4.eta());
  }
}

void matchTaus(const pat::TauCollection& patTaus_test, const std::string& patTauType_test,
	       const pat::TauCollection& patTaus_ref,  const std::string& patTauType_ref,
	       double dRmatch)
{
  for ( pat::TauCollection::const_iterator patTau_test = patTaus_test.begin();
	patTau_test != patTaus_test.end(); ++patTau_test ) {
    reco::Candidate::LorentzVector patTauJetP4_test = getJetP4(*patTau_test);

    bool isMatched = false;

    for ( pat::TauCollection::const_iterator patTau_ref = patTaus_ref.begin();
	  patTau_ref != patTaus_ref.end(); ++patTau_ref ) {
      reco::Candidate::LorentzVector patTauJetP4_ref = getJetP4(*patTau_ref);
      
      if ( reco::deltaR(patTauJetP4_test, patTauJetP4_ref) < dRmatch ) isMatched = true;
    }

    if ( !isMatched ) {
      std::cout << "<matchTaus>:" << std::endl;
      std::cout << " pat::Tau in collection " << patTauType_test << " not matched" 
		<< " by any entry in collection " << patTauType_ref << ":" << std::endl;
      std::cout << "tauPt  = " << patTau_test->pt() << std::endl;
      std::cout << "tauEta = " << patTau_test->eta() << std::endl;
      std::cout << "tauPhi = " << patTau_test->phi() << std::endl;
      std::cout << "jetPt  = " << patTauJetP4_test.pt() << std::endl;
      std::cout << "jetEta = " << patTauJetP4_test.eta() << std::endl;
      std::cout << "jetPhi = " << patTauJetP4_test.phi() << std::endl;
    }
  }
}

void compareTaus(const pat::TauCollection& patTaus1, const std::string& patTauType1,
		 const pat::TauCollection& patTaus2, const std::string& patTauType2,
		 double dRmatch)
{
  if ( patTaus1.size() != patTaus2.size() ) {
    std::cout << "<compareTaus>:" << std::endl;
    std::cout << " mismatch in number of pat::Tau collection entries !!" << std::endl;
    std::cout << "(" << patTauType1 << " = " << patTaus1.size() << "," 
	      << " " << patTauType2 << " = " << patTaus2.size() << ")" << std::endl;
  }

  matchTaus(patTaus1, patTauType1, patTaus2, patTauType2, dRmatch);
  matchTaus(patTaus2, patTauType2, patTaus1, patTauType1, dRmatch);
} 

pat::TauCollection selectTaus(StringCutObjectSelector<pat::Tau>* patTauSelector, const pat::TauCollection& patTaus)
{
  pat::TauCollection patTaus_selected;

  for ( pat::TauCollection::const_iterator patTau = patTaus.begin();
	patTau != patTaus.end(); ++patTau ) {
    bool isSelected = ( patTauSelector ) ? (*patTauSelector)(*patTau) : true;
    if ( isSelected ) patTaus_selected.push_back(*patTau);
  }

  return patTaus_selected;
}

void PATTauCollectionDiffAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{  
  std::cout << "<PATTauCollectionDiffAnalyzer::analyze>:" << std::endl; 

//--- check that configuration parameters contain no errors
  if ( dqmError_ ) {
    edm::LogError("PATTauCollectionDiffAnalyzer::analyze") << " Failed to access dqmStore --> histograms will NOT be filled !!";
    return;
  }

  edm::Handle<pat::TauCollection> patTaus1;
  evt.getByLabel(patTauSource1_, patTaus1);

  edm::Handle<pat::TauCollection> patTaus2;
  evt.getByLabel(patTauSource2_, patTaus2);

  fillHistograms(histogramsPatTaus1_.jetPt_, histogramsPatTaus1_.jetEta_, *patTaus1);
  fillHistograms(histogramsPatTaus2_.jetPt_, histogramsPatTaus2_.jetEta_, *patTaus2);

  std::cout << " --> comparing pat::Taus **before** applying selection..." << std::endl;
  compareTaus(*patTaus1, patTauSource1_.label(), *patTaus2, patTauSource2_.label(), dRmatch_);

  pat::TauCollection patTausSel1 = selectTaus(patTauSelector_, *patTaus1);
  pat::TauCollection patTausSel2 = selectTaus(patTauSelector_, *patTaus2);

  fillHistograms(histogramsPatTausSel1_.jetPt_, histogramsPatTausSel1_.jetEta_, patTausSel1);
  fillHistograms(histogramsPatTausSel2_.jetPt_, histogramsPatTausSel2_.jetEta_, patTausSel2);

  std::cout << " --> comparing pat::Taus **after** applying selection..." << std::endl;
  compareTaus(patTausSel1, patTauSource1_.label(), patTausSel2, patTauSource2_.label(), dRmatch_);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauCollectionDiffAnalyzer);

