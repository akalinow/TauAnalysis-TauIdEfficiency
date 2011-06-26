
#ifndef TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h
#define TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h

#include "DataFormats/Candidate/interface/Candidate.h" 

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

const pat::Jet* getJet_Tau(const pat::Tau&, const pat::JetCollection&);

#endif
