#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

template<typename T>
class TriggerMatch : public edm::EDProducer {
    public:
        explicit TriggerMatch(const edm::ParameterSet & iConfig);
        virtual ~TriggerMatch() ;

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::EDGetTokenT<edm::View<reco::Candidate>> tags_;            
        edm::EDGetTokenT<edm::View<T>> objects_; 
        StringCutObjectSelector<T,true> objCut_; // lazy parsing, to allow cutting on variables not in reco::Candidate class
        double objDR2Tag_;
};

template<typename T>
TriggerMatch<T>::TriggerMatch(const edm::ParameterSet & iConfig) :
    tags_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("tags"))),
    objects_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("objects"))),
    objCut_(iConfig.existsAs<std::string>("objectSelection") ? iConfig.getParameter<std::string>("objectSelection") : "", true),
    objDR2Tag_(std::pow(iConfig.getParameter<double>("maxTagObjDR"),2))
{
    produces<edm::ValueMap<float> >();
}


template<typename T>
TriggerMatch<T>::~TriggerMatch()
{
}

template<typename T>
void 
TriggerMatch<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > tags;
    Handle<View<T> > objects;
    iEvent.getByToken(tags_,  tags);
    iEvent.getByToken(objects_, objects);
    
    // fill
    std::vector<const T *> selObjs;
    typename View<T>::const_iterator object, endobjects = objects->end();
    for (object = objects->begin(); object != endobjects; ++object) {
      if ( (objCut_(*object)) ) selObjs.push_back(&*object);
    }

    std::vector<float> values; 
    values.reserve(tags->size());

    // prepare vector for output    
    View<reco::Candidate>::const_iterator tagIter, endtags = tags->end(); 
    typename std::vector<const T *>::const_iterator selbegin = selObjs.begin(), selend = selObjs.end();
    for (tagIter = tags->begin(); tagIter!= endtags; ++tagIter) {
        const reco::Candidate &tag   = *tagIter;
	values.push_back(false);
        for (typename std::vector<const T *>::const_iterator ito = selbegin; ito != selend; ++ito) { 
	  if (reco::deltaR2(  tag.eta(),   tag.phi(), (*ito)->eta(), (*ito)->phi()) < objDR2Tag_){
	      values.back() = true;
            }
        }
    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(tags, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}


typedef TriggerMatch<pat::TriggerObjectStandAlone> TriggerObjectStandAloneMatch;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectStandAloneMatch);
