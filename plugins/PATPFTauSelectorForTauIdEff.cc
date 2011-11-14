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
    save_(0),
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

  if ( cfg.exists("save") ) {
    std::cout << "<PATPFTauSelectorForTauIdEff::PATPFTauSelectorForTauIdEff>:" << std::endl;
    std::cout << " src = " << src_.label() << std::endl;
    std::string save_string = cfg.getParameter<std::string>("save");
    std::cout << "--> saving pat::Taus passing: " << save_string << std::endl;
    save_ = new StringCutObjectSelector<pat::Tau>(cfg.getParameter<std::string>("save"));
  }

  produceAll_ = ( cfg.exists("produceAll") ) ? 
    cfg.exists("produceAll") : false;

  produces<PATTauCollection>();
}

PATPFTauSelectorForTauIdEff::~PATPFTauSelectorForTauIdEff() 
{
  delete muonSelection_;
  delete trackQualityCuts_;
  delete pfIsolationExtractor_;
  delete save_;
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

    double pfJetJEC = jetEnergyCorrector->correction(*pfJet, evt, es);
    if ( verbosity_ ) {
      std::cout << " PFJet: uncorrected Pt = " << p4PFJetUncorrected.pt() << "," 
    	        << " eta = " << p4PFJetUncorrected.eta() << ", phi = " << p4PFJetUncorrected.phi()
    	        << " --> pfJetJEC = " << pfJetJEC << std::endl;
    }
    reco::Candidate::LorentzVector p4PFJetCorrected(pfJetJEC*p4PFJetUncorrected);
    
    bool isSaved = ( save_ ) ?
      (*save_)(*pfTau_input) : false;
    
//--- check that (PF)tau-jet candidate passes Pt and eta selection
    if ( !(produceAll_ || isSaved) && !(p4PFJetCorrected.pt() > minJetPt_ && TMath::Abs(p4PFJetCorrected.eta()) < maxJetEta_) ) continue;

//--- select PFChargedHadrons passing track quality cuts
//    applied in PFTau reconstruction
    trackQualityCuts_->setPV(theVertex);
    std::vector<reco::PFCandidatePtr> pfChargedJetConstituents = reco::tau::pfChargedCands(*pfJet);
    unsigned numTracks = 0;
    unsigned numSelTracks = 0;
    std::vector<reco::PFCandidatePtr> selPFChargedHadrons;
    for ( std::vector<reco::PFCandidatePtr>::const_iterator pfChargedJetConstituent = pfChargedJetConstituents.begin();
	  pfChargedJetConstituent != pfChargedJetConstituents.end(); ++pfChargedJetConstituent ) {
      if ( TMath::Abs((*pfChargedJetConstituent)->charge()) > 0.5 ) ++numTracks;
      if ( trackQualityCuts_->filter(**pfChargedJetConstituent) ) {
	++numSelTracks;
	selPFChargedHadrons.push_back(*pfChargedJetConstituent);
      }
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
    if ( !(produceAll_ || isSaved) && !leadPFChargedHadron ) continue;

//--- require "leading" (highest Pt) PFChargedHadron to pass Pt cut
//   (corresponding to "leading track Pt cut" applied in PFTau reconstruction,
//    but without dR < 0.1 matching between leading track and jet-axis applied)
    if ( !(produceAll_ || isSaved) && !(leadPFChargedHadron->pt() > minLeadTrackPt_) ) continue;

//--- require that (PF)Tau-jet candidate passes loose isolation criteria
    reco::VertexCollection theVertexCollection;
    theVertexCollection.push_back(*theVertex);
    double loosePFIsoPt = -1.;
    if ( leadPFChargedHadron ) 
      loosePFIsoPt = (*pfIsolationExtractor_)(*pfTau_input, leadPFChargedHadron->momentum(), 
					      *pfIsoCandidates, &theVertexCollection, beamSpot.product(), rhoFastJet);
    if ( verbosity_ ) std::cout << " loosePFIsoPt = " << loosePFIsoPt << std::endl;
    if ( !(produceAll_  || isSaved) && loosePFIsoPt > maxPFIsoPt_ ) continue;

//--- require "leading" PFChargedHadron to pass cut on (anti-)PFElectron MVA output
//   (corresponding to discriminatorAgainstElectrons(Loose) applied in PFTau reconstruction)
    double PFElectronMVA = -1.;
    if ( leadPFChargedHadron ) {
      if ( verbosity_ ) std::cout << " PFElectronMVA = " << leadPFChargedHadron->mva_e_pi() << std::endl;
      PFElectronMVA = leadPFChargedHadron->mva_e_pi();
    }
    bool isElectron = !(PFElectronMVA < maxLeadTrackPFElectronMVA_);
    if ( !(produceAll_ || isSaved) && isElectron ) continue;
    bool isECALcrack = (TMath::Abs(p4PFJetCorrected.eta()) > 1.442 && TMath::Abs(p4PFJetCorrected.eta()) < 1.560);
    if ( !(produceAll_ || isSaved) && applyECALcrackVeto_ && isECALcrack ) continue;

//--- require "leading" PFChargedHadron not to overlap with reconstructed muon
//   (corresponding to discriminatorAgainstMuons applied in PFTau reconstruction;
//    whether the cut against muons is looser or tighter can be controlled 
//    via the 'muonSelection' configuration parameter)
    const pat::Muon* nearestMuon = 0;
    double dRnearestMuon = 1.e+3;
    if ( leadPFChargedHadron ) {
      for ( PATMuonCollection::const_iterator muon = muons->begin();
	    muon != muons->end(); ++muon ) {
	if ( muonSelection_ == 0 || (*muonSelection_)(*muon) ) {
	  double dR = deltaR(leadPFChargedHadron->p4(), muon->p4());
	  if ( dR < dRnearestMuon ) {
	    nearestMuon = &(*muon);
	    dRnearestMuon = dR;	    
	  }
	}
      }
    }
    bool isMuon = (dRnearestMuon < minDeltaRtoNearestMuon_);
    if ( verbosity_ ) {
      if ( nearestMuon && dRnearestMuon < minDeltaRtoNearestMuon_ ) {
	std::cout << " muon: Pt = " << nearestMuon->pt() << "," 
		  << " eta = " << nearestMuon->eta() << ", phi = " << nearestMuon->phi() 
		  << " --> dR = " << dRnearestMuon << std::endl;
	std::cout << "(isGlobalMuon = " << nearestMuon->isGlobalMuon() << ","
		  << " isTrackerMuon = " << nearestMuon->isTrackerMuon() << ","
		  << " isStandAloneMuon = " << nearestMuon->isStandAloneMuon() << ")" << std::endl;
      }
      std::cout << " isMuon = " << isMuon << std::endl;
    }
    if ( !(produceAll_ || isSaved) && isMuon ) continue;

//--- all cuts passed 
//   --> create selected (PF)tau-jet candidate
    pat::Tau pfTau_output(*pfTau_input);

//--- set four-vector of selected (PF)tau-jet candidate
//    to jet-energy corrected four-vector of associated PFJet
    pfTau_output.setP4(p4PFJetCorrected);

//--- store number of tracks/number of tracks passing quality criteria
    pfTau_output.addUserFloat("numTracks",    numTracks);
    pfTau_output.addUserFloat("numSelTracks", numSelTracks);

//--- store Pt and charge of "leading" track
    if ( leadPFChargedHadron ) {
      pfTau_output.addUserFloat("hasLeadTrack",    1.);
      pfTau_output.addUserFloat("leadTrackPt",     leadPFChargedHadron->pt());
      pfTau_output.addUserFloat("leadTrackCharge", leadPFChargedHadron->charge());

//--- store (PF)isolation Pt sum in userFloat variables
      pfTau_output.addUserFloat("preselLoosePFIsoPt", loosePFIsoPt);
    } else {
      pfTau_output.addUserFloat("hasLeadTrack",    0.);
    }

//--- store electron/muon veto flags
    pfTau_output.addUserFloat("PFElectronMVA", PFElectronMVA);
    pfTau_output.addUserFloat("dRnearestMuon", dRnearestMuon);

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
