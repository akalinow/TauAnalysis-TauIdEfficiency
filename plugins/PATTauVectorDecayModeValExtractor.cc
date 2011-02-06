#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorDecayModeValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <string>

PATTauVectorDecayModeValExtractor::PATTauVectorDecayModeValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  srcPFTau_ = cfg.getParameter<edm::InputTag>("srcPFTau");
  srcPFTauDecayMode_ = cfg.getParameter<edm::InputTag>("srcPFTauDecayMode");

  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "energy"   ) value_ = kEnergy;
  else if ( value_string == "pt"       ) value_ = kPt;
  else if ( value_string == "eta"      ) value_ = kEta;
  else if ( value_string == "phi"      ) value_ = kPhi;
  else if ( value_string == "mass"     ) value_ = kMass;
  else {
    edm::LogError ("PATTauVectorJetIdValExtractor") 
      << " Invalid Configuration Parameter 'value' = " << value_string << " !!";
    value_ = -1;
  }
}

const reco::PFTauDecayMode* getDecayMode_pfTau(const pat::Tau& tau, 
					       const edm::Handle<reco::PFTauCollection>& recoPFTaus, 
					       const reco::PFTauDecayModeAssociation& recoPFTauDecayModes) 
{
  const reco::PFTauDecayMode* recoPFTauDecayMode = 0;

  double dRmin = 1.e+3;

  size_t numPFTaus = recoPFTaus->size();
  for ( size_t idx = 0; idx < numPFTaus; ++idx ) {
    reco::PFTauRef recoPFTauRef(recoPFTaus, idx);

    double dR = deltaR(recoPFTauRef->p4(), tau.p4());
    if ( dR < 0.5 && dR < dRmin ) {
      recoPFTauDecayMode = &recoPFTauDecayModes[recoPFTauRef];
      dRmin = dR;
    }
  }

  return recoPFTauDecayMode;
}

std::vector<double> PATTauVectorDecayModeValExtractor::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;

  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> patTaus;
  evt.getByLabel(src_, patTaus);

  edm::Handle<reco::PFTauCollection> recoPFTaus;
  evt.getByLabel(srcPFTau_, recoPFTaus);
  
  edm::Handle<reco::PFTauDecayModeAssociation> recoPFTauDecayModes;
  evt.getByLabel(srcPFTauDecayMode_, recoPFTauDecayModes);

  unsigned numPatTaus = patTaus->size();
  for ( unsigned iTau = 0; iTau < numPatTaus; ++iTau ) {
    edm::Ptr<pat::Tau> patTauPtr = patTaus->ptrAt(iTau);

    double vec_i = -1.;

    const reco::PFTauDecayMode* recoPFTauDecayMode = getDecayMode_pfTau(*patTauPtr, recoPFTaus, *recoPFTauDecayModes);
    
    if ( recoPFTauDecayMode ) {
      if      ( value_ == kEnergy ) vec_i = recoPFTauDecayMode->energy();
      else if ( value_ == kPt     ) vec_i = recoPFTauDecayMode->pt();
      else if ( value_ == kEta    ) vec_i = recoPFTauDecayMode->eta();
      else if ( value_ == kPhi    ) vec_i = recoPFTauDecayMode->phi();
      else if ( value_ == kMass   ) vec_i = recoPFTauDecayMode->mass();
    }

    vec.push_back(vec_i);
  }

  return vec;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorDecayModeValExtractor, "PATTauVectorDecayModeValExtractor");
