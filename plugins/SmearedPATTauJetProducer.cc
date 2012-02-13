#include "PhysicsTools/PatUtils/interface/SmearedJetProducerT.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "RecoMET/METAlgorithms/interface/SignAlgoResolutions.h"
#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"

#include "TauAnalysis/TauIdEfficiency/interface/PATTauJetCorrExtractor.h"

namespace SmearedJetProducer_namespace
{
  template <>
  class JetResolutionExtractorT<pat::Tau>
  {
   public:

    JetResolutionExtractorT(const edm::ParameterSet& cfg) 
      : jetResolutions_(cfg)
    {}
    ~JetResolutionExtractorT() {}
    
    double operator()(const pat::Tau& jet) const
    {
      if ( !jet.isPFTau() )
	throw cms::Exception("SmearedJetProducer::produce")
	  << " Tau-jets of type other than PF not supported yet !!\n";
      
      metsig::SigInputObj pfJetResolution = jetResolutions_.evalPFJet(&(*jet.pfJetRef()));
      if ( pfJetResolution.get_energy() > 0. ) {
	return jet.energy()*(pfJetResolution.get_sigma_e()/pfJetResolution.get_energy());
      } else {
	return 0.;
      }
    }
    
    metsig::SignAlgoResolutions jetResolutions_;
  };
}

typedef SmearedJetProducerT<pat::Tau, PATTauJetCorrExtractor> SmearedPATTauJetProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedPATTauJetProducer);
