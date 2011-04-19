#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorJetIdValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

#include <string>

template<typename T>
PATTauVectorJetIdValExtractor<T>::PATTauVectorJetIdValExtractor(const edm::ParameterSet& cfg)
  : jetId_(0)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  srcJet_ = cfg.getParameter<edm::InputTag>("srcJet");

  std::string value_string = cfg.getParameter<std::string>("value");
  std::string quality;
  if      ( value_string == "minimal"   ) quality = "MINIMAL";
  else if ( value_string == "loose"     ) quality = "LOOSE";
  else if ( value_string == "loose_AOD" ) quality = "LOOSE_AOD";
  else if ( value_string == "medium"    ) quality = "MEDIUM";
  else if ( value_string == "tight"     ) quality = "TIGHT";
  else {
    edm::LogError ("PATTauVectorDecayModeValExtractor") 
      << " Invalid Configuration Parameter 'value' = " << value_string << " !!";
  }

  edm::ParameterSet cfgJetId;
  cfgJetId.addParameter<std::string>("version", "PURE09");
  cfgJetId.addParameter<std::string>("quality", quality);
  
  jetId_ = new T(cfgJetId);
}

template<typename T>
PATTauVectorJetIdValExtractor<T>::~PATTauVectorJetIdValExtractor()
{
  delete jetId_;
}

template<typename T>
std::vector<double> PATTauVectorJetIdValExtractor<T>::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;

  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> patTaus;
  evt.getByLabel(src_, patTaus);

  edm::Handle<pat::JetCollection> patJets;
  evt.getByLabel(srcJet_, patJets);

  for ( patTauCollectionType::const_iterator patTau = patTaus->begin(); 
	patTau != patTaus->end(); ++patTau ) {

    double vec_i = -1.;

    const pat::Jet* patJet = getJet_Tau(*patTau, *patJets);
    
    if ( patJet ) {
      pat::strbitset bits = jetId_->getBitTemplate();
      vec_i = (*jetId_)(*patJet, bits);
    }

    vec.push_back(vec_i);
  }

  return vec;
}

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

typedef PATTauVectorJetIdValExtractor<JetIDSelectionFunctor> PATTauVectorCaloJetIdValExtractor;
typedef PATTauVectorJetIdValExtractor<PFJetIDSelectionFunctor> PATTauVectorPFJetIdValExtractor;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorCaloJetIdValExtractor, "PATTauVectorCaloJetIdValExtractor");
DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorPFJetIdValExtractor, "PATTauVectorPFJetIdValExtractor");

