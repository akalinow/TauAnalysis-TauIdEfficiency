#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

class PFMuonMerger : public edm::stream::EDProducer<> {
public:
  explicit PFMuonMerger(const edm::ParameterSet & iConfig);
  virtual ~PFMuonMerger() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

private:
  
  edm::InputTag muons_;
  StringCutObjectSelector<pat::Muon, false> muonsCut_;

  bool mergeTracks_;
  edm::InputTag tracks_;
  StringCutObjectSelector<pat::PackedCandidate, false> tracksCut_;

  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > trackToken_;
};


PFMuonMerger::PFMuonMerger(const edm::ParameterSet & iConfig) :
    muons_(iConfig.getParameter<edm::InputTag>("muons")),
    muonsCut_(iConfig.existsAs<std::string>("muonsCut") ? iConfig.getParameter<std::string>("muonsCut") : ""),    
    mergeTracks_(iConfig.existsAs<bool>("mergeTracks") ? iConfig.getParameter<bool>("mergeTracks") : false),
    tracks_(mergeTracks_ ? iConfig.getParameter<edm::InputTag>("tracks") : edm::InputTag()),
    tracksCut_(iConfig.existsAs<std::string>("tracksCut") ? iConfig.getParameter<std::string>("tracksCut") : "")
{
  muonToken_ = consumes<std::vector<pat::Muon> >(muons_);
  trackToken_ = consumes<std::vector<pat::PackedCandidate> > (tracks_);
  produces<std::vector<pat::Muon> >();
}

void 
PFMuonMerger::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    edm::Handle<std::vector<pat::Muon> > muons;
    edm::Handle<std::vector<pat::PackedCandidate> > tracks;

    iEvent.getByToken(muonToken_,muons);
    if(mergeTracks_) iEvent.getByToken(trackToken_,tracks);

    std::auto_ptr<std::vector<pat::Muon> >  out(new std::vector<pat::Muon>());
    out->reserve(muons->size() + (mergeTracks_?tracks->size():0));

    // copy reco::Muons, turning on the CaloCompatibility flag if enabled and possible
    for (std::vector<pat::Muon>::const_iterator it = muons->begin(), ed = muons->end(); it != ed; ++it) {
        if(!muonsCut_(*it)) continue;
	out->push_back(*it);
    }
    // merge reco::Track avoiding duplication of innerTracks
    if(mergeTracks_){
        for (size_t i = 0; i < tracks->size(); i++) {
	    pat::PackedCandidateRef track(tracks, i);
            if(!tracksCut_(*track)) continue;
            // check if it is a muon
            bool isMuon = false;
            for(std::vector<pat::Muon>::const_iterator muon = muons->begin(); muon < muons->end(); muon++){
	      if(reco::deltaR(*muon->bestTrack(),*track)<0.007){
                    isMuon = true;
                    break;
                }
            }
            if(isMuon) continue;           
            // make a reco::Muon
            double energy = sqrt(track->p() * track->p() + 0.011163691);
            math::XYZTLorentzVector p4(track->px(), track->py(), track->pz(), energy);
	    reco::Muon aRecoMuon(track->charge(), p4, track->vertex());
            out->push_back(pat::Muon(aRecoMuon));
        }
    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFMuonMerger);
