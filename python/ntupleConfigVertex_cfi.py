import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi import *

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

vertexMultiplicity_template = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    pluginType = cms.string("VertexValExtractor"),

    src = cms.InputTag("offlinePrimaryVertices"),

    columns = cms.PSet(
        numVerticesPtGt5 = cms.string("numVerticesPtGt5"),
        numVerticesPtGt10 = cms.string("numVerticesPtGt10"),
        numVerticesPtGt15 = cms.string("numVerticesPtGt15"),
        numVerticesPtGt20 = cms.string("numVerticesPtGt20"),
    )
)

vertexMultReweight_template = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    pluginType = cms.string("VertexMultiplicityReweightExtractor"),

    src = vertexMultiplicityReweight.src,
    inputFileName = vertexMultiplicityReweight.inputFileName,
    lutName = vertexMultiplicityReweight.lutName,
    type = vertexMultiplicityReweight.type,

    columns = cms.PSet(
        vtxMultReweight = cms.string("vtxMultReweight")
    )
)
