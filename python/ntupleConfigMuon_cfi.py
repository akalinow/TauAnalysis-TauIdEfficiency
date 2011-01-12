import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to Muons
#--------------------------------------------------------------------------------

muons_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATMuonVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patMuonsLoosePFIsoEmbedded06"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for Muon
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),

        # charge of Muon
        charge = cms.string("charge()"),

        # momentum of the "inner" track reconstructed in the silicon Pixel plus Strip detectors,
        # of the "outer" tracks reconstructed in the muon system,
        # and of the "global" track reconstructed by a combined fit of the inner plus outer tracks
        innerTrackPt = cms.string("? innerTrack().isNonnull() ? innerTrack().pt() : 0."),
        outerTrackPt = cms.string("? outerTrack().isNonnull() ? outerTrack().pt() : 0."),
        globalTrackPt = cms.string("? globalTrack().isNonnull() ? globalTrack().pt() : 0."),

        # Pt sum of tracks/Et sum of ECAL/HCAL energy deposits within isolation cone
        trackIso = cms.string("trackIso"),
        ecalIso = cms.string("ecalIso"),
        hcalIso = cms.string("hcalIso"),

        # Pt sum of charged + neutral hadrons and photons within isolation cones of size dR = 0.4/0.6
        ptSumLooseIsolation04 = cms.string("userFloat('pfLooseIsoPt04')"),
        ptSumLooseIsolation06 = cms.string("userFloat('pfLooseIsoPt06')")
    )
)

muons_genInfo = muons_recInfo.clone(
    columns = cms.PSet(
        # generator level information
        genMatch = cms.string("? genParticleRef().isNonnull() ? 1. : 0."),
        
        genPt = cms.string("? genParticleRef().isNonnull() ? genParticleRef().pt() : 0."),
        genEta = cms.string("? genParticleRef().isNonnull() ? genParticleRef().eta() : 0."),
        genPhi = cms.string("? genParticleRef().isNonnull() ? genParticleRef().phi() : 0.")
    )
)    
