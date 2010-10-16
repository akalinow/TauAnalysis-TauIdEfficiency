import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for x,y,z coordinates of reconstructed primary event vertices
#
# NOTE: x and y coordinates are computed with respect to beam-spot
#
#--------------------------------------------------------------------------------

vertex_template = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    pluginType = cms.string("VertexVectorValExtractor"),

    src = cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),

    columns = cms.PSet(
        x = cms.string("vertexX"),
        y = cms.string("vertexY"),
        z = cms.string("vertexZ")
    )
)    
