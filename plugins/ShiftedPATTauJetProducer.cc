#include "PhysicsTools/PatUtils/interface/ShiftedJetProducerT.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

typedef ShiftedJetProducerT<pat::Tau> ShiftedPATTauJetProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ShiftedPATTauJetProducer);
