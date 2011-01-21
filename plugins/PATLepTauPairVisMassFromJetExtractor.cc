#include "TauAnalysis/TauIdEfficiency/plugins/PATLepTauPairVisMassFromJetExtractor.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

template<typename T>
PATLepTauPairVisMassFromJetExtractor<T>::PATLepTauPairVisMassFromJetExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  index_ = ( cfg.exists("index") ) ? cfg.getParameter<unsigned>("index") : 0;
}

template<typename T>
PATLepTauPairVisMassFromJetExtractor<T>::~PATLepTauPairVisMassFromJetExtractor()
{
//--- nothing to be done yet...
}

template<typename T>
double PATLepTauPairVisMassFromJetExtractor<T>::operator()(const edm::Event& evt) const
{
  typedef CompositePtrCandidateT1T2MEt<T, pat::Tau> diTauType;

  typedef edm::View<diTauType> diTauCollectionType;
  edm::Handle<diTauCollectionType> diTauPairs;
  evt.getByLabel(src_, diTauPairs);

  if ( diTauPairs->size() > index_ ) {
    edm::Ptr<diTauType> diTauPairPtr = diTauPairs->ptrAt(index_);

    return (diTauPairPtr->leg1()->p4() + diTauPairPtr->leg2()->pfTauTagInfoRef()->pfjetRef()->p4()).mass();
  } else {
    return -1.;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

typedef PATLepTauPairVisMassFromJetExtractor<pat::Electron> PATElecTauPairVisMassFromJetExtractor;
typedef PATLepTauPairVisMassFromJetExtractor<pat::Muon> PATMuTauPairVisMassFromJetExtractor;

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, PATElecTauPairVisMassFromJetExtractor, "PATElecTauPairVisMassFromJetExtractor");
DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, PATMuTauPairVisMassFromJetExtractor, "PATMuTauPairVisMassFromJetExtractor");
