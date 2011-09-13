#ifndef TauAnalysis_TauIdEfficiency_tauPtResAuxFunctions_h
#define TauAnalysis_TauIdEfficiency_tauPtResAuxFunctions_h

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <TMath.h>

std::string getGenTauDecayMode(const pat::Tau&, const reco::GenParticleCollection&);

std::string getPFCandidateType(const reco::PFCandidate&);
const reco::TrackBaseRef getTrack(const reco::PFCandidate&);

double getTauPtManCorr(const pat::Tau&, const reco::Vertex&, const reco::tau::RecoTauQualityCuts&, unsigned);

void printPatTau(const pat::Tau&);
void printRecoPFJet(const reco::PFJet&, const reco::Vertex&);

#endif
