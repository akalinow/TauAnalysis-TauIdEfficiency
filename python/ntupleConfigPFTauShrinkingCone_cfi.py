import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to PFTaus
# reconstructed by shrinking signal cone algorithm
#--------------------------------------------------------------------------------

pfTausShrinkingCone_template01 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausShrinkingCone"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),

        # tau id. discriminators based on leading track/PFChargedHadron                                                 
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byLeadTrackPtCut = cms.string("tauID('leadingTrackPtCut')"),
        byTrackIsolation = cms.string("tauID('trackIsolation')"),
        byEcalIsolation = cms.string("tauID('ecalIsolation')"),
        byIsolation = cms.string("tauID('byIsolation')"),
        
        # tau id. discriminators based on leading pion (PFChargedHadron or PFGamma)                                    
        byLeadPionPtCut = cms.string("tauID('leadingPionPtCut')"),
        byTrackIsolationUsingLeadingPion = cms.string("tauID('trackIsolationUsingLeadingPion')"),
        byEcalIsolationUsingLeadingPion = cms.string("tauID('ecalIsolationUsingLeadingPion')"),
        byIsolationUsingLeadingPion = cms.string("tauID('byIsolationUsingLeadingPion')"),
        
        # additional TaNC (neural network) based discriminators
        byTaNC = cms.string("tauID('byTaNC')"),
        byTaNCfrOnePercent = cms.string("tauID('byTaNCfrOnePercent')"),
        byTaNCfrHalfPercent = cms.string("tauID('byTaNCfrHalfPercent')"),
        byTaNCfrQuarterPercent = cms.string("tauID('byTaNCfrQuarterPercent')"),
        byTaNCfrTenthPercent = cms.string("tauID('byTaNCfrTenthPercent')"),                                                     
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')")                             
    )
)

pfTausShrinkingCone_template02 = pfTausShrinkingCone_template01.clone(
    pluginType = cms.string("PATTauVectorGenJetValExtractor"),

    columns = cms.PSet(
        # generator level information
        genMatch = cms.string("genMatch"),
        
        genPt = cms.string("genPt"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi"),
        
        genDecayMode = cms.string("genDecayMode")
    )
)    
