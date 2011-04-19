#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorJetCorrMomValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

#include <string>

PATTauVectorJetCorrMomValExtractor::PATTauVectorJetCorrMomValExtractor(const edm::ParameterSet& cfg)
  : stringObjFunction_(cfg.getParameter<std::string>("value"))
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  srcJet_ = cfg.getParameter<edm::InputTag>("srcJet");
}

PATTauVectorJetCorrMomValExtractor::~PATTauVectorJetCorrMomValExtractor()
{
//--- nothing to be done yet...
}

std::vector<double> PATTauVectorJetCorrMomValExtractor::operator()(const edm::Event& evt) const
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
    
    if ( patJet ) vec_i = stringObjFunction_(*patJet);

    vec.push_back(vec_i);
  }

  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorJetCorrMomValExtractor, "PATTauVectorJetCorrMomValExtractor");

