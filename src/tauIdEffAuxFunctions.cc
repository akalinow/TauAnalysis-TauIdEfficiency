#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

#include "DataFormats/Math/interface/deltaR.h"

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

