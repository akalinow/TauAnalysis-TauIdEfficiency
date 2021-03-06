import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os
import re

inputFilePath = \
  '/data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Jul23/V6/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/tauIdEffMeasPATTuple_Ztautau_powheg_2011Jul23V6_5_75fd.root

  /data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Aug18/V7noZrecoilCorr/'
  'user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Aug18'
sampleToAnalyze = 'Ztautau_powheg'
jobId = '2011Aug18'

inputFile_regex = \
  r"tauIdEffMeasPATTuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)

outputFilePath = '/data1/veelken/tmp/muonPtGt20/V7noZrecoilCorr/'

# check if inputFile is PAT-tuple and
# matches sampleToAnalyze, jobId
inputFileNames = []
files = os.listdir(inputFilePath)
for file in files:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(file):
        inputFileNames.append(os.path.join(inputFilePath, file))
#print "inputFileNames = %s" % inputFileNames 
##inputFileNames = [ "".join([inputFilePath, "tauIdEffMeasPATTuple_Ztautau_powheg_2011Jul23V6_23_171e.root"]) ]

# find name of associated "process"
process_matched = None
processes = recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'].keys()
for process in processes:
    for sample in recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][process]['samples']:
        if sample == sampleToAnalyze:
            process_matched = process

if not process_matched:
    raise ValueError("No process associated to sample %s --> skipping !!" % sampleToAnalyze)

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFileNames),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string(os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_%s_%s.root' % (sampleToAnalyze, jobId)))
)

process.tauIdEffAnalyzer = cms.PSet(
    process = cms.string(process_matched),
    type = cms.string(recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['type']),

    regions = cms.vstring(
        ##'ABCD',
        ##'A',
        ##'A1',
        ##'A1p',
        ##'A1f',
        ##'B',
        'B1',
        'B1p',
        'B1f',
        ##'C',
        'C1',
        'C1p',
        'C1f',
        ##'C2',
        ##'C2p',
        ##'C2f',
        ##'D',
        ##'D1',
        'D1p',
        ##'D1f'
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

    binning = cms.PSet(
        sumEt = cms.VPSet(
            cms.PSet(
                subdir = cms.string('sumEtGt450'),
                min = cms.double(450.000000),
                max = cms.double(1000.000000)
            ),
            cms.PSet(
                subdir = cms.string('sumEt350to450'),
                min = cms.double(350.000000),
                max = cms.double(450.000000)
            ),
            cms.PSet(
                subdir = cms.string('sumEtLt250'),
                min = cms.double(0.000000),
                max = cms.double(250.000000)
            ),
            cms.PSet(
                subdir = cms.string('sumEt250to350'),
                min = cms.double(250.000000),
                max = cms.double(350.000000)
            ),
        ),
        tauPt = cms.VPSet(
            cms.PSet(
                subdir = cms.string('tauPt25to30'),
                min = cms.double(25.000000),
                max = cms.double(30.000000)
            ),
            cms.PSet(
                subdir = cms.string('tauPtLt25'),
                min = cms.double(20.000000),
                max = cms.double(25.000000)
            ),
            cms.PSet(
                subdir = cms.string('tauPtGt40'),
                min = cms.double(40.000000),
                max = cms.double(100.000000)
            ),
            cms.PSet(
                subdir = cms.string('tauPt30to40'),
                min = cms.double(30.000000),
                max = cms.double(40.000000)
            ),
        ),
        tauAbsEta = cms.VPSet(
            cms.PSet(
                subdir = cms.string('tauAbsEta19to23'),
                min = cms.double(1.900000),
                max = cms.double(2.300000)
            ),
            cms.PSet(
                subdir = cms.string('tauAbsEta14to19'),
                min = cms.double(1.500000),
                max = cms.double(1.900000)
            ),
            cms.PSet(
                subdir = cms.string('tauAbsEtaLt14'),
                min = cms.double(0.000000),
                max = cms.double(1.400000)
            ),
        ),
        numVertices = cms.VPSet(
            cms.PSet(
                subdir = cms.string('numVertices5to6'),
                min = cms.double(4.500000),
                max = cms.double(6.500000)
            ),
            cms.PSet(
                subdir = cms.string('numVerticesLeq4'),
                min = cms.double(-0.500000),
                max = cms.double(4.500000)
            ),
            cms.PSet(
                subdir = cms.string('numVerticesGt8'),
                min = cms.double(8.500000),
                max = cms.double(20.500000)
            ),
            cms.PSet(
                subdir = cms.string('numVertices7to8'),
                min = cms.double(6.500000),
                max = cms.double(8.500000)
            ),
        )
    ),

    sysShift = cms.string('CENTRAL_VALUE'),

    selEventsFileName = cms.string(os.path.join(outputFilePath, "selEvents_tauIdEff_%s.txt" % sampleToAnalyze)),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5',
        'HLT_IsoMu17_v6',
        'HLT_IsoMu17_v8',
        'HLT_IsoMu17_v9',
        'HLT_IsoMu17_v10',
        'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),
    svFitMassHypothesis = cms.string('psKine_MEt_logM_fit'),
    tauChargeMode = cms.string("tauSignalChargedHadronSum"),
    disableTauCandPreselCuts = cms.bool(False),
    #disableTauCandPreselCuts = cms.bool(True),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag('ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight'),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['events_processed']),
    
    xSection = cms.double(recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['x_sec']),
    
    intLumiData = cms.double(recoSampleDefinitionsTauIdEfficiency_7TeV['TARGET_LUMI']),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
