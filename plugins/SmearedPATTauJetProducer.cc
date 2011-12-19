#include "PhysicsTools/PatUtils/interface/SmearedJetProducerT.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "RecoMET/METAlgorithms/interface/SignAlgoResolutions.h"
#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"

#include "JetMETCorrections/Type1MET/interface/JetCorrExtractorT.h"

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

template <>
class JetCorrExtractorT<pat::Tau>
{
public:
  
  reco::Candidate::LorentzVector operator()(const pat::Tau& rawJet, const std::string& jetCorrLabel, 
					    const edm::Event* evt = 0, const edm::EventSetup* es = 0, 
					    double jetCorrEtaMax = 9.9, 
					    const reco::Candidate::LorentzVector* rawJetP4_specified = 0)
  {
    if ( !rawJet.isPFTau() )
      throw cms::Exception("SmearedJetProducer::produce")
	<< " Tau-jets of type other than PF not supported yet !!\n";
    
    static JetCorrExtractorT<reco::PFJet> jetCorrExtractor_jet;
    return jetCorrExtractor_jet(*rawJet.pfJetRef(), jetCorrLabel, evt, es, jetCorrEtaMax);
  }
};

typedef SmearedJetProducerT<pat::Tau> SmearedPATTauJetProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(SmearedPATTauJetProducer);
