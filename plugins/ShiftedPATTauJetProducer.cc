#include "PhysicsTools/PatUtils/interface/ShiftedJetProducerT.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TauAnalysis/TauIdEfficiency/interface/PATTauJetCorrExtractor.h"

typedef ShiftedJetProducerT<pat::Tau, PATTauJetCorrExtractor> ShiftedPATTauJetProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ShiftedPATTauJetProducer);
