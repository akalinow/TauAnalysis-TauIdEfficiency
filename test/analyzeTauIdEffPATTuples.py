import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os

channel = 'ZtoMuTau_tauIdEff'
#jobId = getJobId(channel)
jobId = '2011Jul23'

inputFilePath = '/data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Jul06_mauro/V4/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/'
outputFilePath = '/data1/veelken/tmp/muonPtGt20/V6/'

sampleToAnalyze = 'Ztautau_powheg'

inputFiles = os.listdir(inputFilePath)
#print(inputFiles)

# check if inputFile is PAT-tuple and
# matches sampleToAnalyze, jobId
inputFiles_sample = []
for inputFile in inputFiles:        
    if inputFile.find("tauIdEffMeasPATTuple") != -1 and \
      inputFile.find("".join(['_', sampleToAnalyze, '_'])) != -1 and \
      inputFile.find("".join(['_', jobId, '_'])) != -1:
        inputFiles_sample.append(os.path.join(inputFilePath, inputFile))

# find name of associated "process"
process_matched = None
processes = recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'].keys()
for process in processes:
    for sample in recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][process]['samples']:
        if sample == sampleToAnalyze:
            process_matched = process

if not process_matched:
    print("No process associated to sample %s --> skipping !!" % sampleToAnalyze)
    continue

processType = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['type']

outputFileName = os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_%s_%s_%s.root' % (sampleToAnalyze, sysUncertainty, jobId))

weights_string = ""
if not recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][process_matched]['type'] == 'Data':
    weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight", "'"])
    #weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#selectedPatMuonsForTauIdEffTrkIPcumulative#muonHLTeff", "'"])

allEvents_DBS = -1
xSection = 0.0
if not recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][process_matched]['type'] == 'Data':
    allEvents_DBS = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['events_processed']
    xSection = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['x_sec']
intLumiData = recoSampleDefinitionsTauIdEfficiency_7TeV['TARGET_LUMI']

#
#--------------------------------------------------------------------------------
#

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFiles_sample),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.tauIdEffAnalyzer = cms.PSet(
    process = cms.string(process_matched),
    type = cms.string(processType),

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',
        'A1p',
        'A1f',
        'B',
        'B1',
        'B1p',
        'B1f',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D',
        'D1',
        'D1p',
        'D1f'
    ),
    
    tauIds = cms.VPSet(
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byLooseIsolation'
            ),
            name = cms.string("tauDiscrHPSloose")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byMediumIsolation'
            ),
            name = cms.string("tauDiscrHPSmedium")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byTightIsolation'
            ),
            name = cms.string("tauDiscrHPStight")
        ),     
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byLooseIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPSlooseDBcorr")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byMediumIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPSmediumDBcorr")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byTightIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPStightDBcorr")
        ),     
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byLooseCombinedIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPScombLooseDBcorr")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byMediumCombinedIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPScombMediumDBcorr")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byTightCombinedIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPScombTightDBcorr")
        )             
    ),

    sysShift = cms.string('CENTRAL_VALUE'),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5', 'HLT_IsoMu17_v6', 'HLT_IsoMu17_v8', 'HLT_IsoMu17_v9', 'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(weights_string),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(allEvents_DBS_value),
    
    xSection = cms.double(xSection_value),
    
    intLumiData = cms.double(intLumiData_value),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
