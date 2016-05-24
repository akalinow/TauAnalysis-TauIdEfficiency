#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

class PFTauMerger : public edm::stream::EDProducer<> {
public:
  explicit PFTauMerger(const edm::ParameterSet & iConfig);
  virtual ~PFTauMerger() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

private:
  
  edm::InputTag taus_;
  StringCutObjectSelector<pat::Tau, false> tausCut_;

  bool mergeTracks_;
  edm::InputTag tracks_;
  StringCutObjectSelector<pat::PackedCandidate, false> tracksCut_;

  edm::InputTag photons_;

  edm::EDGetTokenT<std::vector<pat::Tau> > tauToken_;
  edm::EDGetTokenT<std::vector<pat::Photon> > photonToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > trackToken_;
};


PFTauMerger::PFTauMerger(const edm::ParameterSet & iConfig) :
    taus_(iConfig.getParameter<edm::InputTag>("taus")),
    tausCut_(iConfig.existsAs<std::string>("tausCut") ? iConfig.getParameter<std::string>("tausCut") : ""),    
    mergeTracks_(iConfig.existsAs<bool>("mergeTracks") ? iConfig.getParameter<bool>("mergeTracks") : false),
    tracks_(mergeTracks_ ? iConfig.getParameter<edm::InputTag>("tracks") : edm::InputTag()),
    tracksCut_(iConfig.existsAs<std::string>("tracksCut") ? iConfig.getParameter<std::string>("tracksCut") : ""),
    photons_(iConfig.getParameter<edm::InputTag>("photons"))
{
  tauToken_ = consumes<std::vector<pat::Tau> >(taus_);
  photonToken_ = consumes<std::vector<pat::Photon> >(photons_);
  trackToken_ = consumes<std::vector<pat::PackedCandidate> > (tracks_);
  produces<std::vector<pat::Tau> >();
}

void 
PFTauMerger::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    edm::Handle<std::vector<pat::Tau> > taus;
    edm::Handle<std::vector<pat::Photon> > photons;
    edm::Handle<std::vector<pat::PackedCandidate> > tracks;

    iEvent.getByToken(photonToken_,photons);
    iEvent.getByToken(tauToken_,taus);
    if(mergeTracks_) iEvent.getByToken(trackToken_,tracks);

    std::auto_ptr<std::vector<pat::Tau> >  out(new std::vector<pat::Tau>());
    out->reserve(taus->size() + (mergeTracks_?tracks->size():0));

    // copy reco::Taus, turning on the CaloCompatibility flag if enabled and possible
    for (std::vector<pat::Tau>::const_iterator it = taus->begin(), ed = taus->end(); it != ed; ++it) {
        if(!tausCut_(*it)) continue;
	//out->push_back(*it);
    }
    // merge reco::Track avoiding duplication of innerTracks
    if(mergeTracks_){
        for (size_t i = 0; i < tracks->size(); i++) {
	    pat::PackedCandidateRef track(tracks, i);
            if(!tracksCut_(*track)) continue;
            // check if it is a tau	    
            bool isTau = false;
	    const pat::Tau *aTau = 0;
            for(std::vector<pat::Tau>::const_iterator tau = taus->begin(); tau < taus->end(); tau++){
	      if(reco::deltaR(tau->leadChargedHadrCand()->p4(),*track)<0.007){
		isTau = true;
		aTau = &(*tau);
		break;
                }
            }	    
            // make a pat::Tau from a track
            double energy = sqrt(track->p() * track->p() + 0.13957018);
            math::XYZTLorentzVector p4(track->px(), track->py(), track->pz(), energy);
	    /*
	    for (std::vector<pat::Photon>::const_iterator itPhoton = photons->begin(); itPhoton!=photons->end(); ++itPhoton){
	      float isolation = itPhoton->chargedHadronIso() + itPhoton->neutralHadronIso() + itPhoton->photonIso();
	      if(itPhoton->p4().E()/p4.E()<1000.1 &&
		 ((reco::deltaR(p4,*itPhoton)<0.07 && itPhoton->pt()>2) ||
		 (reco::deltaR(p4,*itPhoton)>0.07 && reco::deltaR(p4,*itPhoton)<0.5 && itPhoton->pt()>4 && isolation<1.0))){
		p4+=itPhoton->p4();
		//p4.SetE(p4.E()+itPhoton->p4().E());
	      }
	    }
	    */
	    reco::BaseTau aBaseTau(track->charge(), p4, track->vertex());
	    pat::Tau aPatTau(aBaseTau);
	    std::vector<pat::Tau::IdPair> aID;	   
            if(isTau){	     	      
	      aPatTau.setalternatLorentzVect(aTau->p4());
	      aID.push_back(pat::Tau::IdPair("decayMode",aTau->decayMode()));
	      aID.push_back(pat::Tau::IdPair("againstMuonTight3",aTau->tauID("againstMuonTight3")));
	      aID.push_back(pat::Tau::IdPair("againstMuonLoose3",aTau->tauID("againstMuonLoose3")));
	      aID.push_back(pat::Tau::IdPair("decayModeFindingNewDMs",aTau->tauID("decayModeFindingNewDMs")));
	      aID.push_back(pat::Tau::IdPair("decayModeFinding",aTau->tauID("decayModeFinding")));
	      aID.push_back(pat::Tau::IdPair("byLooseCombinedIsolationDeltaBetaCorr3Hits",aTau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")));
	    }
	    else{
	      aID.push_back(pat::Tau::IdPair("decayMode",-1));
	      aID.push_back(pat::Tau::IdPair("againstMuonTight3",0));
	      aID.push_back(pat::Tau::IdPair("againstMuonLoose3",0));
	      aID.push_back(pat::Tau::IdPair("decayModeFindingNewDMs",0));
	      aID.push_back(pat::Tau::IdPair("decayModeFinding",0));
	      aID.push_back(pat::Tau::IdPair("byLooseCombinedIsolationDeltaBetaCorr3Hits",0));
	    }	    
	    aPatTau.setTauIDs(aID);
            out->push_back(aPatTau);
        }
    }
    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFTauMerger);
