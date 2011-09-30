
#ifndef TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h
#define TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h

#include "DataFormats/Candidate/interface/Candidate.h" 

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

const pat::Jet* getJet_Tau(const pat::Tau&, const pat::JetCollection&);

enum { kUnmatched, kJetToTauFakeMatched, kMuToTauFakeMatched, kGenTauHadMatched, kGenTauOtherMatched };

int getGenMatchType(const PATMuTauPair&, const reco::GenParticleCollection&, double* = 0, double* = 0);

#endif
