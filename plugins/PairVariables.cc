#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/TauIdEfficiency/interface/GenHelper.h"

template<typename T>
class PairVariables : public edm::EDProducer {
    public:
        explicit PairVariables(const edm::ParameterSet & iConfig);
        virtual ~PairVariables() ;

        virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
        edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;
        edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticles_;
        std::string variableName_;
       
};

template<typename T>
PairVariables<T>::PairVariables(const edm::ParameterSet & iConfig) :
    pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
    genParticles_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
    variableName_(iConfig.getParameter<std::string>("variableName"))
{
    produces<edm::ValueMap<float> >();
}


template<typename T>
PairVariables<T>::~PairVariables(){ }

template<typename T>
void 
PairVariables<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;

    // read input
    Handle<View<reco::Candidate> > pairs;
    Handle<View<reco::GenParticle> > genParticles;
    iEvent.getByToken(pairs_,  pairs);
    
    // fill
    std::vector<float> values; 
    values.reserve(pairs->size());

    // prepare vector for output    
    View<reco::Candidate>::const_iterator pair, endpairs = pairs->end(); 
    for (pair = pairs->begin(); pair != endpairs; ++pair) {
        if (pair->numberOfDaughters() != 2) throw cms::Exception("LogicError", "Pairs must have *two* daughters");
        const reco::CandidateBaseRef tag   = pair->daughter(0)->masterClone();
	const reco::CandidateBaseRef probe = pair->daughter(1)->masterClone();

	//StringObjectFunction<reco::Candidate,true> function("alternatLorentzVect().pt()");
	//float pt = function(*probe);
	
	if(variableName_=="alternativeMass"){
	  const pat::Tau* aTauCandidate = (const pat::Tau*)probe.get();
	  float massAlternative = (tag->p4() + aTauCandidate->alternatLorentzVect()).M();
	  values.push_back(massAlternative);
	}
	if(variableName_=="ZDecayMode"){
	  genhelper::HZDecay decay = genhelper::Other;
	  iEvent.getByToken(genParticles_, genParticles);
	  for (unsigned int iGen = 0; iGen < genParticles->size(); iGen++){
	    const reco::GenParticle genP = genParticles->at(iGen);
	    bool isLast = genhelper::IsLastCopy(genP);
	    int pdgId = genP.pdgId();	    
	    if (isLast && abs(pdgId)==23){
		decay = genhelper::GetHZDecay (&genP);	  
		break;
	    }
	  }
	  values.push_back(decay);
	}
    }

    // convert into ValueMap and store
    std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    ValueMap<float>::Filler filler(*valMap);
    filler.insert(pairs, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap);
}


typedef PairVariables<reco::Candidate> CandPairVariables;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CandPairVariables);
