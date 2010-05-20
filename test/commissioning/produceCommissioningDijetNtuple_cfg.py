import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

from May6thPDSkim2_SD_JetMETTau_files import data_files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_5_5/RelValZTT/GEN-SIM-RECO/START3X_V25-v1/0009/FAB3991D-8039-DF11-8E2E-002618FDA277.root',
        #'/store/data/Commissioning10/MinimumBias/RECO/May6thPDSkim2_SD_JetMETTau-v1/0137/FEADF2F9-D25D-DF11-91B3-002618943935.root',
        #data_files
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# Load the selections to apply on dijet events
process.load("TauAnalysis.TauIdEfficiency.TauIdDijetEventSelection_cfi")

from TauAnalysis.TauIdEfficiency.tools.sequenceBuilder \
        import buildDijetTauSequence

from TauAnalysis.TauIdEfficiency.patConfiguration.tauTypeDefintions \
        import patTauProducerOptions

# Build the sequence for shrinking cone taus
shrinkingConeSequenceInfo = buildDijetTauSequence(
    process, niceName = "shrinkingCone", 
    cleanerOverlapCheckerSource = "pfJetsTagAndProbes",
    patTauProdOptions = patTauProducerOptions['shrinkingCone']
)
# Output is a dict containing the relevant sequence and 
# a list of the important output collections
process.buildShrinkingConeTaus = shrinkingConeSequenceInfo['sequence']

# Build the sequence for the fixed cone taus
fixedConeSequenceInfo = buildDijetTauSequence(
    process, niceName = "fixedCone", cleanerOverlapCheckerSource = "pfJetsTagAndProbes",
    patTauProdOptions = patTauProducerOptions['fixedCone']
)
process.buildFixedConeTaus = fixedConeSequenceInfo['sequence']

#  TODO factorize this stuff out into separate function
# Build the ntuple skeleton
process.ntupleProducer = cms.EDAnalyzer(
    "ObjValEDNtupleProducer",
    ntupleName = cms.string("tauCommissioning"),
    sources = cms.PSet()
)

# Load the definitions of the variables to put into the ntuples
from TauAnalysis.TauIdEfficiency.ntupleVariables import common
from TauAnalysis.TauIdEfficiency.ntupleVariables import pftau
from TauAnalysis.TauIdEfficiency.ntupleVariables import tancDiscriminators

# Add the shrinkingCone collections to the ntuple
for ntupleInputName in shrinkingConeSequenceInfo['outputCollections']:
    #build this source
    newSource = cms.PSet(
        columns = cms.PSet(
            tancDiscriminators.variables,
            common.variables,
            pftau.variables,
        ),
        vector = cms.bool(True),
        pluginType = cms.string('PATTauVectorValExtractor'),
        src = cms.InputTag(ntupleInputName)
    )
    # Add to ntuple producer
    setattr(process.ntupleProducer.sources, ntupleInputName, newSource)

# Add the fixedCone collections to the ntuple
for ntupleInputName in fixedConeSequenceInfo['outputCollections']:
    #build this source
    newSource = cms.PSet(
        columns = cms.PSet(
            common.variables,
            pftau.variables,
        ),
        vector = cms.bool(True),
        pluginType = cms.string('PATTauVectorValExtractor'),
        src = cms.InputTag(ntupleInputName)
    )
    # Add to ntuple producer
    setattr(process.ntupleProducer.sources, ntupleInputName, newSource)

# Setup our path to run
process.p = cms.Path(
    process.pfJetTagAndProbeSequence
    *process.buildShrinkingConeTaus
    *process.buildFixedConeTaus
    *process.ntupleProducer
)

# Get trigger report
process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
)

# Save ntuple
process.out = cms.OutputModule(
    "PoolOutputModule",                                                                                                                                                        
    outputCommands = cms.untracked.vstring("drop *", "keep *_*ntupleProducer*_*_*" ),
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("commissioning_ntuple.root")      
)

process.end = cms.EndPath(process.out)
