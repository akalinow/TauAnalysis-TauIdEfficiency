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
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

template<typename T>
class MTPair : public edm::EDProducer {
    public:
        explicit MTPair(const edm::ParameterSet & iConfig);
        virtual ~MTPair() ;

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;            
        edm::EDGetTokenT<edm::View<T>> objects_; 
        StringCutObjectSelector<T,true> objCut_; // lazy parsing, to allow cutting on variables not in reco::Candidate class
        double objDR2Tag_, objDR2Probe_;
        bool useProbe_;
};

template<typename T>
MTPair<T>::MTPair(const edm::ParameterSet & iConfig) :
    pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
    objects_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("objects"))),
    objCut_(iConfig.existsAs<std::string>("objectSelection") ? iConfig.getParameter<std::string>("objectSelection") : "", true),
    objDR2Tag_(std::pow(iConfig.getParameter<double>("minTagObjDR"),2)),
    objDR2Probe_(std::pow(iConfig.getParameter<double>("minProbeObjDR"),2)),
    useProbe_(iConfig.getParameter<bool>("useProbe"))
{
    produces<edm::ValueMap<float> >();
}


template<typename T>
MTPair<T>::~MTPair(){ }

template<typename T>
void 
MTPair<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > pairs;
    Handle<View<T> > objects;
    iEvent.getByToken(pairs_,  pairs);
    iEvent.getByToken(objects_, objects);
    
    // fill
    std::vector<const T *> selObjs;
    typename View<T>::const_iterator object, endobjects = objects->end();
    for (object = objects->begin(); object != endobjects; ++object) {
      if ( (objCut_(*object)) ) selObjs.push_back(&*object);
    }

    std::vector<float> values; 
    values.reserve(pairs->size());

    // prepare vector for output    
    View<reco::Candidate>::const_iterator pair, endpairs = pairs->end(); 
    typename std::vector<const T *>::const_iterator selbegin = selObjs.begin(), selend = selObjs.end();
    for (pair = pairs->begin(); pair != endpairs; ++pair) {
        if (pair->numberOfDaughters() != 2) throw cms::Exception("LogicError", "Pairs must have *two* daughters");
        const reco::Candidate &tag   = *pair->daughter(0);
        const reco::Candidate &probe = *pair->daughter(1);
	
        for (typename std::vector<const T *>::const_iterator ito = selbegin; ito != selend; ++ito) { 
            if (reco::deltaR2(  tag.eta(),   tag.phi(), (*ito)->eta(), (*ito)->phi()) > objDR2Tag_  &&
                reco::deltaR2(probe.eta(), probe.phi(), (*ito)->eta(), (*ito)->phi()) > objDR2Probe_) {
	     
	      float mtTag = sqrt(2*(*ito)->p4().pt()*tag.pt()*(1 - cos((*ito)->p4().phi() - tag.p4().phi())));
	      float mtProbe = sqrt(2*(*ito)->p4().pt()*probe.pt()*(1 - cos((*ito)->p4().phi() - probe.p4().phi())));
	      if(useProbe_) values.push_back(mtProbe);
	      else  values.push_back(mtTag);
            }
        }
    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(pairs, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}


typedef MTPair<reco::Candidate> CandMTPair;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CandMTPair);
