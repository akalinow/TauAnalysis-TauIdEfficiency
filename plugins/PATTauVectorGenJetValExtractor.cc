#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorGenJetValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "TauAnalysis/GenSimTools/interface/genParticleAuxFunctions.h"

#include <string>

namespace {
   // Forward declaration
   template<typename T> const reco::GenJet* getGenJet(edm::Ptr<T> input);

   // Getter for when input collection is pat::taus
   template<> const reco::GenJet* getGenJet<pat::Tau>(edm::Ptr<pat::Tau> input) 
   { 
      return input->genJet();
   }

   // Getter for when input collection is already a collection of
   // GenJets.
   template<> const reco::GenJet* getGenJet<reco::GenJet>(edm::Ptr<reco::GenJet> input) 
   {
      return input.get();
   }
}
  
template<typename T>
VectorGenJetValExtractor<T>::VectorGenJetValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  std::string value_string = cfg.getParameter<std::string>("value");
//--- variables for "true" hadronic tau decays
  if      ( value_string == "genMatch"     ) value_ = kGenMatch;
  else if ( value_string == "genPt"        ) value_ = kGenPt;
  else if ( value_string == "genEta"       ) value_ = kGenEta;
  else if ( value_string == "genPhi"       ) value_ = kGenPhi;
  else if ( value_string == "genMass"      ) value_ = kGenMass;
  else if ( value_string == "genDecayMode" ) value_ = kGenDecayMode;
//--- variables for quark/gluon jets faking signature of hadronic tau decays
  else if ( value_string == "genPdgId"     ) {
    value_ = kGenPdgId;
    srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
    skipPdgIdsGenParticleMatch_ = cfg.getParameter<vint>("skipPdgIdsGenParticleMatch");
  } else {
    edm::LogError ("VectorGenJetValExtractor") << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

template<typename T>
std::vector<double> VectorGenJetValExtractor<T>::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;
  typedef edm::View<T> inputCollectionType;
  edm::Handle<inputCollectionType> input;
  evt.getByLabel(src_, input);

  edm::Handle<reco::GenParticleCollection> genParticles;
  if ( srcGenParticles_.label() != "" ) evt.getByLabel(srcGenParticles_, genParticles);

  unsigned nInput = input->size();
  for ( unsigned i = 0; i < nInput; ++i ) {
    edm::Ptr<T> inputPtr = input->ptrAt(i);

    double vec_i = -1.;

    // get the associated GenJet
    const reco::GenJet* genJet = getGenJet<T>(inputPtr);

    bool isMatched = ( genJet != NULL ) ? true : false;

    if      ( value_ == kGenMatch ) vec_i = isMatched;
    else if ( isMatched ) {
      if      ( value_ == kGenPt        ) vec_i = genJet->pt();
      else if ( value_ == kGenEta       ) vec_i = genJet->eta();
      else if ( value_ == kGenPhi       ) vec_i = genJet->phi();
      else if ( value_ == kGenMass      ) vec_i = genJet->mass();
      else if ( value_ == kGenDecayMode ) {	
	std::string genDecayMode_string = JetMCTagUtils::genTauDecayMode(*genJet);
//--- decode generated tau decay mode
//    ( as defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc )
	if      ( genDecayMode_string == "electron"        ) vec_i = 0;
	else if ( genDecayMode_string == "muon"            ) vec_i = 1;
	else if ( genDecayMode_string == "oneProng0Pi0"    ) vec_i = 2;
	else if ( genDecayMode_string == "oneProng1Pi0"    ) vec_i = 3;
	else if ( genDecayMode_string == "oneProng2Pi0"    ) vec_i = 4;
	else if ( genDecayMode_string == "oneProngOther"   ) vec_i = 5;
	else if ( genDecayMode_string == "threeProng0Pi0"  ) vec_i = 6;
	else if ( genDecayMode_string == "threeProng1Pi0"  ) vec_i = 7;
	else if ( genDecayMode_string == "threeProngOther" ) vec_i = 8;
	else if ( genDecayMode_string == "rare"            ) vec_i = 9;
	else {
	  edm::LogError ("VectorGenJetValExtractor::operator()") 
	    << " Undefined genDecayMode = " << genDecayMode_string << " --> returning -1 !!";
	  vec_i = -1;
	}
      } 
    } else if ( value_ == kGenPdgId ) vec_i = getMatchingGenParticlePdgId(inputPtr->p4(), *genParticles, &skipPdgIdsGenParticleMatch_);

    vec.push_back(vec_i);
  }

  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

typedef VectorGenJetValExtractor<pat::Tau> PATTauVectorGenJetValExtractor;
typedef VectorGenJetValExtractor<reco::GenJet> GenJetVectorValExtractor;

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, GenJetVectorValExtractor, "GenJetVectorValExtractor");
DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorGenJetValExtractor, "PATTauVectorGenJetValExtractor");
