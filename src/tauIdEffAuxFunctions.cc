#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

const pat::Jet* getJet_Tau(const pat::Tau& tau, const pat::JetCollection& patJets) 
{
  const pat::Jet* retVal = 0;

  reco::Candidate::LorentzVector tauJetP4;
  if      ( tau.isPFTau()   ) tauJetP4 = tau.pfJetRef()->p4();
  else if ( tau.isCaloTau() ) tauJetP4 = tau.caloTauTagInfoRef()->jetRef()->p4();

  double dRmin = 1.e+3;
  
  for ( pat::JetCollection::const_iterator patJet = patJets.begin();
	patJet != patJets.end(); ++patJet ) {
    double dR = deltaR(patJet->correctedJet("Uncorrected").p4(), tauJetP4);
    if ( dR < 0.5 && dR < dRmin ) {
      retVal = &(*patJet);
      dRmin = dR;
    }
  }

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

int getGenMatchType(const PATMuTauPair& muTauPair, const reco::GenParticleCollection& genParticles,
		    double* genTauCharge, double* recTauCharge)
{
//--- check if reconstructed tau-jet candidate matches "true" hadronic tau decay on generator level,
//    is a "fake" tau (i.e. matches a quark/gluon/e/mu/photon on generator level)
//    or fails to be matched to any generator level object
//
//    NOTE: code to perform matching taken from TauAnalysis/Core/plugins/TauHistManager.cc
//
  //std::cout << "<getGenMatchType>:" << std::endl;

  const reco::GenParticle* matchingGenParticle = findGenParticle(muTauPair.leg2()->p4(), genParticles);
  int matchingGenParticleAbsPdgId = ( matchingGenParticle ) ?
    TMath::Abs(matchingGenParticle->pdgId()) : 0;
  //std::cout << " matchingGenParticleAbsPdgId = " << matchingGenParticleAbsPdgId << std::endl;
  
  if ( genTauCharge ) {
    if ( matchingGenParticle ) (*genTauCharge) = matchingGenParticle->charge();
    else                       (*genTauCharge) = 0.;
  }

  if ( recTauCharge ) (*recTauCharge) = muTauPair.leg2()->charge();

  std::string genTauDecayMode = ( matchingGenParticle && matchingGenParticleAbsPdgId == 15 ) ?
    getGenTauDecayMode(matchingGenParticle) : "";
  //std::cout << " genTauDecayMode = " << genTauDecayMode << std::endl;

  if      ( (matchingGenParticleAbsPdgId >=  1 && matchingGenParticleAbsPdgId <=  6) ||
	     matchingGenParticleAbsPdgId == 11 || matchingGenParticleAbsPdgId == 13  ||
	     matchingGenParticleAbsPdgId == 22                                            ) return kJetToTauFakeMatched;
  else if (  matchingGenParticleAbsPdgId == 11 || 
	    (matchingGenParticleAbsPdgId == 15 &&  genTauDecayMode == "muon")             ) return kMuToTauFakeMatched;
  else if (  matchingGenParticleAbsPdgId == 15 && (genTauDecayMode == "oneProng0Pi0"    ||
						   genTauDecayMode == "oneProng1Pi0"    ||
						   genTauDecayMode == "oneProng2Pi0"    ||
						   genTauDecayMode == "oneProngOther"   ||
						   genTauDecayMode == "threeProng0Pi0"  ||
						   genTauDecayMode == "threeProngOther" ||
						   genTauDecayMode == "threeProngOther" ||
						   genTauDecayMode == "rare"            ) ) return kGenTauHadMatched;    
  else if (  matchingGenParticleAbsPdgId == 15                                            ) return kGenTauOtherMatched;
  else                                                                                      return kUnmatched;
}

