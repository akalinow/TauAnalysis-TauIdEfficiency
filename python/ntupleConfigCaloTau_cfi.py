import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to CaloTaus
#--------------------------------------------------------------------------------

caloTaus_template01 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patCaloTaus"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),
        
        # tau id. discriminators based on leading track
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byLeadTrackPtCut = cms.string("tauID('leadingTrackPtCut')"),
        byIsolation = cms.string("tauID('byIsolation')"),
        
        # discriminators against electrons
        againstElectron = cms.string("tauID('againstElectron')")
    )
)

caloTaus_template02 = caloTaus_template01.clone(
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
