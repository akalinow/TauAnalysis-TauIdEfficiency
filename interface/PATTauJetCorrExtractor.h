#ifndef TauAnalysis_TauIdEfficiency_PATTauJetCorrExtractor_h
#define TauAnalysis_TauIdEfficiency_PATTauJetCorrExtractor_h

#include "JetMETCorrections/Type1MET/interface/JetCorrExtractorT.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

class PATTauJetCorrExtractor
{
public:
  
  reco::Candidate::LorentzVector operator()(const pat::Tau& rawJet, const std::string& jetCorrLabel, 
                                            const edm::Event* evt = 0, const edm::EventSetup* es = 0, 
                                            double jetCorrEtaMax = 9.9, 
                                            const reco::Candidate::LorentzVector* rawJetP4_specified = 0)
  {
    if ( !rawJet.isPFTau() )
      throw cms::Exception("PATTauJetCorrExtractor::operator()")
        << " Tau-jets of type other than PF not supported yet !!\n";
    
    static JetCorrExtractorT<reco::PFJet> jetCorrExtractor_jet;
    return jetCorrExtractor_jet(*rawJet.pfJetRef(), jetCorrLabel, evt, es, jetCorrEtaMax);
  }
};

#endif
