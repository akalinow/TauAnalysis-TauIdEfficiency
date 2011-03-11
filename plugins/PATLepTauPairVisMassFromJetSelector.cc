#include "TauAnalysis/TauIdEfficiency/plugins/PATLepTauPairVisMassFromJetSelector.h"

#include "DataFormats/Common/interface/Handle.h"

template<typename T>
PATLepTauPairVisMassFromJetSelector<T>::PATLepTauPairVisMassFromJetSelector(const edm::ParameterSet& cfg)
{
  std::cout << "<PATLepTauPairVisMassFromJetExtractor>:" << std::endl;

  src_ = cfg.getParameter<edm::InputTag>("src");

  minVisMass_ = ( cfg.exists("minVisMass") ) ?
    cfg.getParameter<double>("minVisMass") : -1.;
  std::cout << " minVisMass = " << minVisMass_ << std::endl;
  maxVisMass_ = ( cfg.exists("maxVisMass") ) ?
    cfg.getParameter<double>("maxVisMass") : +1.e+6;
  std::cout << " maxVisMass = " << maxVisMass_ << std::endl;
}

template<typename T>
PATLepTauPairVisMassFromJetSelector<T>::~PATLepTauPairVisMassFromJetSelector()
{
//--- nothing to be done yet...
}

template <typename T>
void PATLepTauPairVisMassFromJetSelector<T>::select(const edm::Handle<collection>& diTauCollection,
						    edm::Event& evt, const edm::EventSetup& es) 
{
  selected_.clear();

  for ( typename collection::const_iterator diTauCandidate = diTauCollection->begin(); 
	diTauCandidate != diTauCollection->end(); ++diTauCandidate ) {
    const edm::Ptr<pat::Tau> patTauPtr = diTauCandidate->leg2();
    reco::Candidate::LorentzVector p4Tau;
    if      ( patTauPtr->isPFTau()   ) p4Tau = patTauPtr->pfJetRef()->p4();
    else if ( patTauPtr->isCaloTau() ) p4Tau = patTauPtr->caloTauTagInfoRef()->calojetRef()->p4();
    else assert(0);

    double visMass = (diTauCandidate->leg1()->p4() + p4Tau).mass();
    
    if ( visMass > minVisMass_ && 
	 visMass < maxVisMass_ ) {
      selected_.push_back(&(*diTauCandidate));
    }
  }
}

#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"

typedef ObjectSelector<PATLepTauPairVisMassFromJetSelector<pat::Electron> > PATElecTauPairVisMassFromJetSelector;
typedef ObjectSelector<PATLepTauPairVisMassFromJetSelector<pat::Muon> > PATMuTauPairVisMassFromJetSelector;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATElecTauPairVisMassFromJetSelector);
DEFINE_FWK_MODULE(PATMuTauPairVisMassFromJetSelector);
