#include "TauAnalysis/TauIdEfficiency/plugins/PATPFTauSelectorForTauIdEff.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

typedef std::vector<pat::Tau> PATTauCollection;

PATPFTauSelectorForTauIdEff::PATPFTauSelectorForTauIdEff(const edm::ParameterSet& cfg)
  : trackQualityCuts_(0),
    muonSelection_(0),				
    pfIsolationExtractor_(0),
    verbosity_(0)
{
  filter_ = cfg.getParameter<bool>("filter");

  src_ = cfg.getParameter<edm::InputTag>("src");

  jetEnergyCorrection_ = cfg.getParameter<std::string>("jetEnergyCorrection");

  minJetPt_ = cfg.getParameter<double>("minJetPt");
  maxJetEta_ = cfg.getParameter<double>("maxJetEta");

  trackQualityCuts_ = new reco::tau::RecoTauQualityCuts(cfg.getParameter<edm::ParameterSet>("trackQualityCuts"));
  minLeadTrackPt_ = cfg.getParameter<double>("minLeadTrackPt");
  maxDzLeadTrack_ = cfg.getParameter<double>("maxDzLeadTrack");
  maxLeadTrackPFElectronMVA_ = cfg.getParameter<double>("maxLeadTrackPFElectronMVA");
  applyECALcrackVeto_ = cfg.getParameter<bool>("applyECALcrackVeto");

  minDeltaRtoNearestMuon_ = cfg.getParameter<double>("minDeltaRtoNearestMuon");
  if ( cfg.exists("muonSelection") ) {
    muonSelection_ = new StringCutObjectSelector<pat::Muon>(cfg.getParameter<std::string>("muonSelection"));
  }
  srcMuon_ = cfg.getParameter<edm::InputTag>("srcMuon");

  pfIsolationExtractor_ = new ParticlePFIsolationExtractor<pat::Tau>(cfg.getParameter<edm::ParameterSet>("pfIsolation"));
  maxPFIsoPt_ = cfg.getParameter<double>("maxPFIsoPt");
  srcPFIsoCandidates_ = cfg.getParameter<edm::InputTag>("srcPFIsoCandidates");
  srcBeamSpot_   = cfg.getParameter<edm::InputTag>("srcBeamSpot");
  srcVertex_     = cfg.getParameter<edm::InputTag>("srcVertex");
  if ( cfg.exists("srcRhoFastJet") ) {
    srcRhoFastJet_ = cfg.getParameter<edm::InputTag>("srcRhoFastJet");
  }

  produces<PATTauCollection>();
}

PATPFTauSelectorForTauIdEff::~PATPFTauSelectorForTauIdEff() 
{
  delete muonSelection_;
  delete trackQualityCuts_;
  delete pfIsolationExtractor_;
}

bool PATPFTauSelectorForTauIdEff::filter(edm::Event& evt, const edm::EventSetup& es)
{
  if ( verbosity_ ) {
    std::cout << "<PATPFTauSelectorForTauIdEff::filter>:" << std::endl;
    std::cout << " src = " << src_.label() << std::endl;
  }

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(srcVertex_, vertices);
  reco::VertexRef theVertex;
  if ( vertices->size() >= 1 ) theVertex = reco::VertexRef(vertices, 0);
  
  edm::Handle<reco::PFCandidateCollection> pfIsoCandidates;
  evt.getByLabel(srcPFIsoCandidates_, pfIsoCandidates);

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_, beamSpot);

  double rhoFastJet = -1.;
  if ( srcRhoFastJet_.label() != "" ) {
    edm::Handle<double> rhoFastJetHandle;
    evt.getByLabel(srcRhoFastJet_, rhoFastJetHandle);
    if ( rhoFastJetHandle.isValid() ) rhoFastJet = (*rhoFastJetHandle);
  }

  typedef std::vector<pat::Muon> PATMuonCollection;
  edm::Handle<PATMuonCollection> muons;
  evt.getByLabel(srcMuon_, muons);
  
  edm::Handle<PATTauCollection> pfTaus_input;
  evt.getByLabel(src_, pfTaus_input);

  const JetCorrector* jetEnergyCorrector = JetCorrector::getJetCorrector(jetEnergyCorrection_, es);
  
  std::auto_ptr<PATTauCollection> pfTaus_output(new PATTauCollection());

  for ( PATTauCollection::const_iterator pfTau_input = pfTaus_input->begin();
	pfTau_input != pfTaus_input->end(); ++pfTau_input ) {
    reco::PFJetRef pfJet = pfTau_input->pfJetRef();

    reco::Candidate::LorentzVector p4PFJetUncorrected = pfJet->p4();

    double pfJetJEC = jetEnergyCorrector->correction(*pfJet, edm::RefToBase<reco::Jet>(pfJet), evt, es);
    if ( verbosity_ ) {
      std::cout << " PFJet: uncorrected Pt = " << p4PFJetUncorrected.pt() << "," 
    	        << " eta = " << p4PFJetUncorrected.eta() << ", phi = " << p4PFJetUncorrected.phi()
    	        << " --> pfJetJEC = " << pfJetJEC << std::endl;
    }
    reco::Candidate::LorentzVector p4PFJetCorrected(pfJetJEC*p4PFJetUncorrected);
    
//--- check that (PF)tau-jet candidate passes Pt and eta selection
    if ( !(p4PFJetCorrected.pt() > minJetPt_ && TMath::Abs(p4PFJetCorrected.eta()) < maxJetEta_) ) continue;

//--- select PFChargedHadrons passing track quality cuts
//    applied in PFTau reconstruction
    trackQualityCuts_->setPV(theVertex);
    std::vector<reco::PFCandidatePtr> pfChargedJetConstituents = reco::tau::pfChargedCands(*pfJet);
    std::vector<reco::PFCandidatePtr> selPFChargedHadrons;
    for ( std::vector<reco::PFCandidatePtr>::const_iterator pfChargedJetConstituent = pfChargedJetConstituents.begin();
	  pfChargedJetConstituent != pfChargedJetConstituents.end(); ++pfChargedJetConstituent ) {
      if ( trackQualityCuts_->filter(**pfChargedJetConstituent) ) selPFChargedHadrons.push_back(*pfChargedJetConstituent);
    }

//--- find highest Pt "leading" PFChargedHadron
    const reco::PFCandidate* leadPFChargedHadron = 0;
    for ( std::vector<reco::PFCandidatePtr>::const_iterator selPFChargedHadron = selPFChargedHadrons.begin();
	  selPFChargedHadron != selPFChargedHadrons.end(); ++selPFChargedHadron ) {
      if ( leadPFChargedHadron == 0 || (*selPFChargedHadron)->pt() > leadPFChargedHadron->pt() ) 
	leadPFChargedHadron = selPFChargedHadron->get();
    }

    if ( verbosity_ ) {
      std::cout << " leadPFChargedHadron: ";
      if ( leadPFChargedHadron ) std::cout << "Pt = " << leadPFChargedHadron->pt() << "," 
					   << " eta = " << leadPFChargedHadron->eta() << ", phi = " << leadPFChargedHadron->phi();
      else std::cout << " None";
      std::cout << std::endl;
    }

//--- require at least one PFChargedHadron passing quality cuts
//   (corresponding to "leading track finding" applied in PFTau reconstruction,
//    but without dR < 0.1 matching between leading track and jet-axis applied)
    if ( !leadPFChargedHadron ) continue;

//--- require "leading" PFChargedHadron to pass Pt cut
//   (corresponding to "leading track Pt cut" applied in PFTau reconstruction,
//    but without dR < 0.1 matching between leading track and jet-axis applied)
    if ( !(leadPFChargedHadron->pt() > minLeadTrackPt_) ) continue;

//--- require that (PF)Tau-jet candidate passes loose isolation criteria
    reco::VertexCollection theVertexCollection;
    theVertexCollection.push_back(*theVertex);
    double loosePFIsoPt = (*pfIsolationExtractor_)(*pfTau_input, leadPFChargedHadron->momentum(), 
						   *pfIsoCandidates, &theVertexCollection, beamSpot.product(), rhoFastJet);
    if ( verbosity_ ) std::cout << " loosePFIsoPt = " << loosePFIsoPt << std::endl;
    if ( loosePFIsoPt > maxPFIsoPt_ ) continue;

//--- require "leading" PFChargedHadron to pass cut on (anti-)PFElectron MVA output
//   (corresponding to discriminatorAgainstElectrons(Loose) applied in PFTau reconstruction)
    if ( verbosity_ ) std::cout << " PFElectronMVA = " << leadPFChargedHadron->mva_e_pi() << std::endl;
    if ( !(leadPFChargedHadron->mva_e_pi() < maxLeadTrackPFElectronMVA_) ) continue;
    bool isECALcrack = (TMath::Abs(p4PFJetCorrected.eta()) > 1.442 && TMath::Abs(p4PFJetCorrected.eta()) < 1.560);
    if ( applyECALcrackVeto_ && isECALcrack ) continue;

//--- require "leading" PFChargedHadron not to overlap with reconstructed muon
//   (corresponding to discriminatorAgainstMuons applied in PFTau reconstruction;
//    whether the cut against muons is looser or tighter can be controlled 
//    via the 'muonSelection' configuration parameter)
    bool isMuon = false;
    for ( PATMuonCollection::const_iterator muon = muons->begin();
	  muon != muons->end(); ++muon ) {
      if ( muonSelection_ == 0 || (*muonSelection_)(*muon) ) {
	double dR = deltaR(leadPFChargedHadron->p4(), muon->p4());
	if ( dR < minDeltaRtoNearestMuon_ ) {
	  if ( verbosity_ ) {
	    std::cout << " muon: Pt = " << muon->pt() << "," 
		      << " eta = " << muon->eta() << ", phi = " << muon->phi() 
		      << " --> dR = " << dR << std::endl;
	    std::cout << "(isGlobalMuon = " << muon->isGlobalMuon() << ","
		      << " isTrackerMuon = " << muon->isTrackerMuon() << ","
		      << " isStandAloneMuon = " << muon->isStandAloneMuon() << ")" << std::endl;
	  }
	  isMuon = true;
	}
      }
    }
    if ( verbosity_ ) std::cout << " isMuon = " << isMuon << std::endl;
    if ( isMuon ) continue;

//--- all cuts passed 
//   --> create selected (PF)tau-jet candidate
    pat::Tau pfTau_output(*pfTau_input);

//--- set four-vector of selected (PF)tau-jet candidate
//    to jet-energy corrected four-vector of associated PFJet
    pfTau_output.setP4(p4PFJetCorrected);

//--- store (PF)isolation Pt sum in userFloat variables
    pfTau_output.addUserFloat("preselLoosePFIsoPt", loosePFIsoPt);

//--- store Pt, eta, phi and mass of "original" PFTau four-vector in userFloat variables
    pfTau_output.addUserFloat("origTauPt",   pfTau_input->pt());
    pfTau_output.addUserFloat("origTauEta",  pfTau_input->eta());
    pfTau_output.addUserFloat("origTauPhi",  pfTau_input->phi());
    pfTau_output.addUserFloat("origTauMass", pfTau_input->mass());    

    pfTaus_output->push_back(pfTau_output);
  }

//--- sort (PF)tau-jet candidates in order of decreasing Pt
  std::sort(pfTaus_output->begin(), pfTaus_output->end(), pfTauPtComparator_);

  size_t numPFTaus_output = pfTaus_output->size();
  if ( verbosity_ ) std::cout << " numPFTaus_output = " << numPFTaus_output << std::endl;

  evt.put(pfTaus_output);

  bool retVal = ( filter_ && numPFTaus_output == 0 ) ? false : true;
  if ( verbosity_ ) std::cout << "--> returning retVal = " << retVal << std::endl;
  
  return retVal;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATPFTauSelectorForTauIdEff);
