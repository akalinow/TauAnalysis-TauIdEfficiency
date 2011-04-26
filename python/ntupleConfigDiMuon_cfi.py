import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to Mu+ Mu- pairs
#--------------------------------------------------------------------------------

diMuons_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("CandidateVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("goldenZmumuCandidatesGe2IsoMuons"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables of Mu+ Mu- system
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),
        
        # invariant mass of Mu+ Mu- pair
        mass = cms.string("mass()")
    )
)
