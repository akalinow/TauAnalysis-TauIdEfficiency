import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to Muons
#--------------------------------------------------------------------------------

diTaus_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATMuTauPairVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),
    
    # Variables to compute for this source
    columns = cms.PSet(
	Mt = cms.string("mt1MET()"),  
        pZeta = cms.string("pZeta()"),
        pZetaVis = cms.string("pZetaVis()"),
    )
)
