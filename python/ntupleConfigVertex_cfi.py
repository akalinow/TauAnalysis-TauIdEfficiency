import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables for x,y,z coordinates of reconstructed primary event vertex
#
# NOTE: x and y coordinates are computed with respect to beam-spot
#
#--------------------------------------------------------------------------------

vertex_template = cms.PSet(
    vector = cms.bool(False),
    
    pluginType = cms.string("VertexValExtractor"),

    src = cms.InputTag("offlinePrimaryVertices"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),

    columns = cms.PSet(
        vertexX = cms.string("vertexX"),
        vertexY = cms.string("vertexY"),
        vertexZ = cms.string("vertexZ")
    )
)    
