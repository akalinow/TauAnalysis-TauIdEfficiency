#!/usr/bin/env python

import copy
import os
import time

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdCommissioning_7TeV_grid_cfi import recoSampleDefinitionsTauIdCommissioning_7TeV
import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauFakeRateAnalysis import *
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import buildConfigFile_hadd
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *
from TauAnalysis.Configuration.tools.jobtools import make_bsub_script
from TauAnalysis.Configuration.tools.harvestingLXBatch import make_harvest_scripts
from TauAnalysis.Configuration.tools.harvesting import castor_source

version = 'patV2_0'

inputFilePath  = '/castor/cern.ch/user/v/veelken/TauIdCommissioning/'
harvestingFilePath = '/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/harvesting/TauIdCommissioning/'
outputFilePath = '/data1/veelken/tmp/tauFakeRateAnalysis/'

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
]

eventSelectionsToAnalyze = [
    # modify in case you want to submit jobs for some of the event selections only...
]

skipFWLiteTauFakeRateAnalyzer = False
#skipFWLiteTauFakeRateAnalyzer = True

#runPeriod = '2011RunA'
runPeriod = '2011RunB'

hltPaths_qcdDiJet = [
    'HLT_Jet30_v1',
    'HLT_Jet30_v2',
    'HLT_Jet30_v3',
    'HLT_Jet30_v4',
    'HLT_Jet30_v5',
    'HLT_Jet30_v6',
    'HLT_Jet30_v9',
    'HLT_Jet60_v1',
    'HLT_Jet60_v2',
    'HLT_Jet60_v3',
    'HLT_Jet60_v4',
    'HLT_Jet60_v5',
    'HLT_Jet60_v6',
    'HLT_Jet60_v9',
    'HLT_Jet80_v1',
    'HLT_Jet80_v2',
    'HLT_Jet80_v3',
    'HLT_Jet80_v4',
    'HLT_Jet80_v5',
    'HLT_Jet80_v6',
    'HLT_Jet110_v1',
    'HLT_Jet110_v2',
    'HLT_Jet110_v3',
    'HLT_Jet110_v4',
    'HLT_Jet110_v5',
    'HLT_Jet110_v6',
    'HLT_Jet110_v9',
    'HLT_Jet150_v1',
    'HLT_Jet150_v2',
    'HLT_Jet150_v3',
    'HLT_Jet150_v4',
    'HLT_Jet150_v5',
    'HLT_Jet150_v6'
]

hltPaths_qcdMuEnriched = [
    'HLT_Mu15_v2',
    'HLT_Mu15_v3',
    'HLT_Mu15_v4',
    'HLT_Mu15_v5',
    'HLT_Mu15_v6',
    'HLT_Mu15_v8',
    'HLT_Mu15_v9',
    'HLT_Mu15_v12',
    'HLT_Mu15_v13'
]

intLumiData = None
hltPaths = None
srcWeights = None
if runPeriod == '2011RunA':
    samplesToAnalyze.extend([
        'data_Jet_Run2011A_May10ReReco_v1',
        'data_Jet_Run2011A_PromptReco_v4',
        'data_Jet_Aug05ReReco_v1',
        'data_Jet_Run2011A_PromptReco_v6', 
        'data_SingleMu_Run2011A_May10ReReco_v1',
        'data_SingleMu_Run2011A_PromptReco_v4',
        'data_SingleMu_Run2011A_Aug05ReReco_v1',
        'data_SingleMu_Run2011A_PromptReco_v6'
    ])
    intLumiData = 2.12e+3 # runs 160431-173692
    hltPaths = {
        'qcdDiJet'      : hltPaths_qcdDiJet,
        'qcdMuEnriched' : hltPaths_qcdMuEnriched,
        'WplusJets'     : [
            'HLT_IsoMu17_v5',
            'HLT_IsoMu17_v6',
            'HLT_IsoMu17_v8',
            'HLT_IsoMu17_v9',
            'HLT_IsoMu17_v10',
            'HLT_IsoMu17_v11',
            'HLT_IsoMu17_v13'
        ],
        'Zmumu' : [
            'HLT_IsoMu17_v5',
            'HLT_IsoMu17_v6',
            'HLT_IsoMu17_v8',
            'HLT_IsoMu17_v9',
            'HLT_IsoMu17_v10',
            'HLT_IsoMu17_v11',
            'HLT_IsoMu17_v13'
        ]
    }
    srcWeights = {
        'Data' : [],
        'smMC' : [ 'vertexMultiplicityReweight3dRunA' ]
    }
elif runPeriod == '2011RunB':
    samplesToAnalyze.extend([
        'data_Jet_Run2011B_PromptReco_v1a', 
        'data_SingleMu_Run2011B_PromptReco_v1a'
    ])
    intLumiData = 2.53e+3 # runs 175860-179411
    hltPaths = {
        'qcdDiJet'      : hltPaths_qcdDiJet,
        'qcdMuEnriched' : hltPaths_qcdMuEnriched,
        'WplusJets'     : [
            'HLT_IsoMu24_v1',
            'HLT_IsoMu24_v2',
            'HLT_IsoMu24_v4',
            'HLT_IsoMu24_v5',
            'HLT_IsoMu24_v6',
            'HLT_IsoMu24_v7',
            'HLT_IsoMu24_v8',
            'HLT_IsoMu24_v9',
            'HLT_IsoMu24_v12',
            'HLT_IsoMu24_v13'
        ],
        'Zmumu' : {
            'Data' : [
                'HLT_Mu17_Mu8_v7',
                'HLT_Mu17_Mu8_v10'
            ],
            'smMC' : [
                'HLT_DoubleMu7_v1'
            ]
        }
    }
    srcWeights = {
        'Data' : [],
        ##'smMC' : [ 'vertexMultiplicityReweight3dRunB' ]
        'smMC' : []
    }
else:
    raise ValueError("Invalid runPeriod = %s !!" % runPeriod)

eventSelections = {
    'QCDj30' : {
        'jobNameInRecoSampleDef' : 'qcdDiJet',
        'tauJetCandSelection'    : [ "userFloat('probeJet30') > 0.5" ], 
        'inputFileNames'         : "tauCommissioningQCDdiJetPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausDijetTagAndProbeHPS',
        'legendEntry'            : 'QCDj30',
        'markerStyleData'        : 20,
        'markerStyleSim'         : 24,
        'color'                  : color_black.value()
    },
    'QCDj60' : {
        'jobNameInRecoSampleDef' : 'qcdDiJet',
        'tauJetCandSelection'    : [ "userFloat('probeJet60') > 0.5" ], 
        'inputFileNames'         : "tauCommissioningQCDdiJetPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausDijetTagAndProbeHPS',
        'legendEntry'            : 'QCDj60',
        'markerStyleData'        : 21,
        'markerStyleSim'         : 25,
        'color'                  : color_red.value()
    },
    'QCDj80' : {
        'jobNameInRecoSampleDef' : 'qcdDiJet',
        'tauJetCandSelection'    : [ "userFloat('probeJet80') > 0.5" ], 
        'inputFileNames'         : "tauCommissioningQCDdiJetPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausDijetTagAndProbeHPS',
        'legendEntry'            : 'QCDj80',
        'markerStyleData'        : 22,
        'markerStyleSim'         : 26,
        'color'                  : color_green.value()
    },
    'QCDj110' : {
        'jobNameInRecoSampleDef' : 'qcdDiJet',
        'tauJetCandSelection'    : [ "userFloat('probeJet110') > 0.5" ], 
        'inputFileNames'         : "tauCommissioningQCDdiJetPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausDijetTagAndProbeHPS',
        'legendEntry'            : 'QCDj110',
        'markerStyleData'        : 23,
        'markerStyleSim'         : 32,
        'color'                  : color_lightBlue.value()
    },
    'QCDj150' : {
        'jobNameInRecoSampleDef' : 'qcdDiJet',
        'tauJetCandSelection'    : [ "userFloat('probeJet150') > 0.5" ], 
        'inputFileNames'         : "tauCommissioningQCDdiJetPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausDijetTagAndProbeHPS',
        'legendEntry'            : 'QCDj150',
        'markerStyleData'        : 33,
        'markerStyleSim'         : 27,
        'color'                  : color_violett.value()
    },
    'QCDmu' : {
        'jobNameInRecoSampleDef' : 'qcdMuEnriched',
        'tauJetCandSelection'    : [], 
        'inputFileNames'         : "tauCommissioningQCDmuEnrichedPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausLoosePFIsoEmbedded06HPS',
        'legendEntry'            : 'QCD#mu',
        'markerStyleData'        : 33,
        'markerStyleSim'         : 27,
        'color'                  : color_darkBlue.value()
    },
    'Wmunu' : {
        'jobNameInRecoSampleDef' : 'WplusJets',
        'tauJetCandSelection'    : [], 
        'inputFileNames'         : "tauCommissioningWplusJetsEnrichedPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausLoosePFIsoEmbedded06HPS',
        'legendEntry'            : 'W #rightarrow #mu #nu',
        'markerStyleData'        : 33,
        'markerStyleSim'         : 27,
        'color'                  : color_red.value()
    },
    'Zmumu' : {
        'jobNameInRecoSampleDef' : 'Zmumu',
        'tauJetCandSelection'    : [], 
        'inputFileNames'         : "tauCommissioningZmumuEnrichedPATtuple.root",
        'srcTauJetCandidates'    : 'patPFTausLoosePFIsoEmbedded06HPS',
        'legendEntry'            : 'Z #rightarrow #mu^{+} #mu^{-}',
        'markerStyleData'        : 33,
        'markerStyleSim'         : 27,
        'color'                  : color_violett.value()
    }
}

tauIds = {
    # HPS isolation with no deltaBeta corrections applied
    # (separate isolation requirements wrt. PFChargedHadrons and PFGammas)
    'tauDiscrHPSloose'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseIsolation'
        ],
        'legendEntry' : "HPS Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPSmedium' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumIsolation'
        ],
        'legendEntry' : "HPS Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPStight' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightIsolation'
        ],
        'legendEntry' : "HPS Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    },
            
    # HPS isolation with deltaBeta corrections applied
    # (separate isolation requirements wrt. PFChargedHadrons and PFGammas)
    'tauDiscrHPSlooseDBcorr'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPSmediumDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPStightDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    },
    
    # HPS combined isolation discriminators
    # (based on isolation sumPt of PFChargedHadrons + PFGammas)
    'tauDiscrHPScombLooseDBcorr'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPScombMediumDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPScombTightDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    }
}

labels = [
    'CMS Preliminary 2011',
    '#sqrt{s} = 7 TeV, L = %1.2f fb^{-1}' % intLumiData/1.e+3
]    

srcMET = 'patPFMet' 

configFilePath = os.path.join(os.getcwd(), "lxbatch")
logFilePath    = os.path.join(os.getcwd(), "lxbatch_log")

if not os.path.exists(configFilePath):
    os.makedirs(configFilePath)
if not os.path.exists(logFilePath):
    os.makedirs(logFilePath)

execDir = "%s/bin/%s/" % (os.environ['CMSSW_BASE'], os.environ['SCRAM_ARCH'])

executable_FWLiteTauFakeRateAnalyzer = execDir + 'FWLiteTauFakeRateAnalyzer'
executable_bsub = 'bsub'
executable_waitForLXBatchJobs = 'python %s/src/TauAnalysis/Configuration/python/tools/waitForLXBatchJobs.py' % os.environ['CMSSW_BASE']
executable_rfcp = 'rfcp'
executable_rfrm = '- rfrm' # CV: ignore error code returned by 'rfrm' in case file on castor does not exist
executable_hadd = 'hadd -f'
executable_makeTauFakeRatePlots = execDir + 'makeTauFakeRatePlots'
executable_shell = '/bin/csh'

bsubQueue = "1nw"

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdCommissioning_7TeV['SAMPLES_TO_RUN']
if len(eventSelectionsToAnalyze) == 0:
    eventSelectionsToAnalyze = eventSelections.keys()

print "samplesToAnalyze = %s" % samplesToAnalyze
print "eventSelectionsToAnalyze = %s" % eventSelectionsToAnalyze

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauFakeRateAnalyzer macro on lxbatch
#
fileNames_FWLiteTauFakeRateAnalyzer         = {}
bsubJobNames_FWLiteTauFakeRateAnalyzer      = {}
bjobListFileNames_FWLiteTauFakeRateAnalyzer = {}
for sampleToAnalyze in samplesToAnalyze:
    fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze]         = {}
    bsubJobNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze]      = {}
    bjobListFileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze] = {}
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:
            tauJetCandSelection = eventSelections[eventSelectionToAnalyze]['tauJetCandSelection']
            srcTauJetCandidates = eventSelections[eventSelectionToAnalyze]['srcTauJetCandidates']            
            hltPaths_sample = None
            if if isinstance(hltPaths[eventSelectionToAnalyze], dict):
                processType = recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['type']
                hltPaths_sample = hltPaths[eventSelectionToAnalyze][processType]
            else:
                hltPaths_sample = hltPaths[eventSelectionToAnalyze]
            retVal_FWLiteTauFakeRateAnalyzer = \
              buildConfigFile_FWLiteTauFakeRateAnalyzer(sampleToAnalyze, eventSelectionToAnalyze, version,
                                                        os.path.join(inputFilePath, jobNameInRecoSampleDef, version, sampleToAnalyze),
                                                        tauIds, tauJetCandSelection, srcTauJetCandidates, srcMET, hltPaths_sample, 
                                                        configFilePath, logFilePath, harvestingFilePath,
                                                        recoSampleDefinitionsTauIdCommissioning_7TeV)
        
            if retVal_FWLiteTauFakeRateAnalyzer is not None:
                fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze] = retVal_FWLiteTauFakeRateAnalyzer
            else:
                fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze] = {
                    'inputFileNames'  : [],
                    'configFileNames' : [],
                    'outputFileNames' : [],
                    'logFileNames'    : []
                }
            fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]['bsubScriptFileNames'] = []
  
            bsubJobNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze] = []
            bjobListFileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze] = None

            if retVal_FWLiteTauFakeRateAnalyzer is None:
                continue

            for i in range(len(retVal_FWLiteTauFakeRateAnalyzer['inputFileNames'])):

                # The None in the tuple indicates that batch job has no dependencies on other batch jobs
                input_files_and_jobs = \
                  [ (None, os.path.join(inputFilePath, jobNameInRecoSampleDef, version, sampleToAnalyze,
                                        retVal_FWLiteTauFakeRateAnalyzer['inputFileNames'][i])) ]

                def log_file_maker(job_hash):
                    return os.path.join(logFilePath, retVal_FWLiteTauFakeRateAnalyzer['logFileNames'][i])

                # Build script for batch job submission
                jobName, bsubScript = make_bsub_script(
                    os.path.join(harvestingFilePath, retVal_FWLiteTauFakeRateAnalyzer['outputFileNames'][i]),
                    input_files_and_jobs,
                    log_file_maker,
                    "%s %s" % (executable_FWLiteTauFakeRateAnalyzer,
                               os.path.join(configFilePath, retVal_FWLiteTauFakeRateAnalyzer['configFileNames'][i])))

                #print "configFilePath = %s" % configFilePath
                #print "retVal_FWLiteTauFakeRateAnalyzer['logFileNames'][i] = %s" % retVal_FWLiteTauFakeRateAnalyzer['logFileNames'][i]
            
                bsubScriptFileName = \
                    os.path.join(configFilePath, retVal_FWLiteTauFakeRateAnalyzer['logFileNames'][i].replace(".log", ".sh"))
                bsubScriptFile = open(bsubScriptFileName, "w")
                bsubScriptFile.write(bsubScript)
                bsubScriptFile.close()

                fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]['bsubScriptFileNames'].append(
                bsubScriptFileName)

                bsubJobName = "tauFRana%s%s_%i" % (sampleToAnalyze, eventSelectionToAnalyze, i)
                bsubJobNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze].append(bsubJobName)

            bjobListFileName = \
              os.path.join(configFilePath, "batchJobs_FWLiteTauFakeRateAnalyzer_%s_%s.lst" % (sampleToAnalyze, eventSelectionToAnalyze))
            bjobListFile = open(bjobListFileName, "w")
            for bsubJobName in bsubJobNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]:
                bjobListFile.write("%s\n" % bsubJobName)
            bjobListFile.close()
        
            bjobListFileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze] = bjobListFileName
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to collect histograms
# for single sample and single event selection into single .root file
#
bsubFileNames_harvesting    = {}
bsubJobNames_harvesting     = {}
bsubJobNames_harvesting_all = []
for sampleToAnalyze in samplesToAnalyze:
    bsubFileNames_harvesting[sampleToAnalyze] = {}
    bsubJobNames_harvesting[sampleToAnalyze]  = {}
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:

            plot_regex = r"[a-zA-Z0-9._]+"
            skim_regex = r"dont match anything"
        
            def local_copy_mapper(sample):
                return os.path.join(
                  outputFilePath,
                  'analyzeTauFakeRateHistograms_%s_%s_%s_harvested.root' % (eventSelectionToAnalyze, sampleToAnalyze, version))

            inputFileInfos = []
            for inputFileName in fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]['outputFileNames']:
                inputFileInfo = {
                    'path'        : os.path.join(harvestingFilePath, inputFileName),
                    'size'        : 1,           # dummy
                    'time'        : time.localtime(),
                    'file'        : inputFileName,
                    'permissions' : 'mrw-r--r--' # "ordinary" file access permissions
                }
                #print "inputFileInfo = %s" % inputFileInfo
                inputFileInfos.append(inputFileInfo)

            retVal_make_harvest_scripts = make_harvest_scripts(
                plot_regex,
                skim_regex,
                channel = eventSelectionToAnalyze,
                sampleToAnalyze = sampleToAnalyze,
                job_id = version,
                input_files_info = inputFileInfos,
                harvester_command = executable_hadd,
                abort_on_rfcp_error = False,
                castor_output_directory = harvestingFilePath,
                script_directory = configFilePath,
                merge_script_name = \
                  os.path.join(configFilePath, "_".join(['submit', sampleToAnalyze, eventSelectionToAnalyze, 'merge']) + '.sh'),
                local_copy_mapper = local_copy_mapper,
                chunk_size = 2.e+9, # 2 GB
                check_old_files = False,
                max_bsub_concurrent_file_access = 250,
                verbosity = 0
            )

            bsubFileNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze] = retVal_make_harvest_scripts

            bsubJobName = "harvest%s%s" % (sampleToAnalyze, eventSelectionToAnalyze)
            bsubJobNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze] = bsubJobName
  
            bsubJobNames_harvesting_all.append(bsubJobName)

bjobListFileName_harvesting = os.path.join(configFilePath, "batchJobs_harvesting_all.lst")
bjobListFile_harvesting = open(bjobListFileName_harvesting, "w")
for sampleToAnalyze in samplesToAnalyze:
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:
            for bsubJobName in bsubFileNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze]['bsub_job_names']:        
                bjobListFile_harvesting.write("%s\n" % bsubJobName)
bjobListFile_harvesting.close()
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to collect histograms
# for all samples and event selections into single .root file
#
haddInputFileNames = []
for sampleToAnalyze in samplesToAnalyze:
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:
            for final_harvest_file in bsubFileNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze]['final_harvest_files']:
                # CV:
                #    (1) file name of final harvesting output file is stored at index[1] in final_harvest_file-tuple
                #       (cf. TauAnalysis/Configuration/python/tools/harvestingLXBatch.py)
                #    (2) assume that .root files containing histograms for single sample and single event selection
                #        are copied to local disk via rfcp prior to running 'hadd'
                haddInputFileNames.append(os.path.join(outputFilePath, os.path.basename(final_harvest_file[1])))
haddShellFileName = os.path.join(configFilePath, 'harvestTauFakeRateHistograms_%s.csh' % version)
haddOutputFileName = os.path.join(outputFilePath, 'analyzeTauFakeRateHistograms_all_%s.root' % version)
retVal_hadd = \
  buildConfigFile_hadd(executable_hadd, haddShellFileName, haddInputFileNames, haddOutputFileName)
haddLogFileName = retVal_hadd['logFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# make jet --> tau fake-rate plots as function of tauPt, tauEta,... for:
#  o QCD muon enriched, W + jets and Zmumu event selections
#  o QCD multi-jet events triggered by different Jet trigger Pt thresholds
#
fileNames_makeTauFakeRatePlots = []

evtSelQCDmu_WplusJets_Zmumu = [
    'QCDmu',
    'Wmunu',
    'Zmumu'
]    

evtSelQCDj = [
    'QCDj30',
    'QCDj60',
    'QCDj80',
    'QCDj110',
    'QCDj150'
]    

evtSelJobs = {
    'evtSelQCDmu_WplusJets_Zmumu' : evtSelQCDmu_WplusJets_Zmumu,
    'evtSelQCDj'                  : evtSelQCDj
}

for evtSelName, evtSelJob in evtSelJobs.items():

    # make plots for HPS isolation with no deltaBeta corrections applied
    discriminators_HPS = [
        'tauDiscrHPSloose',
        'tauDiscrHPSmedium',
        'tauDiscrHPStight'
    ]
    outputFileName_HPS = 'makeTauFakeRatePlots_HPS_%s.eps' % evtSelName
    retVal_makeTauFakeRatePlots = \
      buildConfigFile_makeTauFakeRatePlots(haddOutputFileName, eventSelections, evtSelJob, tauIds, discriminators_HPS, labels,
                                           outputFilePath, outputFileName_HPS, recoSampleDefinitionsTauIdCommissioning_7TeV)
    fileNames_makeTauFakeRatePlots.append(retVal_makeTauFakeRatePlots)

    # make plots for HPS isolation with no applied deltaBeta corrections
    discriminators_HPSdbCorr = [
        'tauDiscrHPSlooseDBcorr',
        'tauDiscrHPSmediumDBcorr',
        'tauDiscrHPStightDBcorr'
    ]
    outputFileName_HPSdbCorr = 'makeTauFakeRatePlots_HPSdbCorr_%s.eps' % evtSelName
    retVal_makeTauFakeRatePlots = \
      buildConfigFile_makeTauFakeRatePlots(haddOutputFileName, eventSelections, evtSelJob, tauIds, discriminators_HPSdbCorr, labels,
                                           outputFilePath, outputFileName_HPSdbCorr, recoSampleDefinitionsTauIdCommissioning_7TeV)
    fileNames_makeTauFakeRatePlots.append(retVal_makeTauFakeRatePlots)
          
    # make plots for HPS combined isolation discriminators
    discriminators_HPScombined = [
        'tauDiscrHPScombLooseDBcorr',
        'tauDiscrHPScombMediumDBcorr',
        'tauDiscrHPScombTightDBcorr'
    ]
    outputFileName_HPScombined = 'makeTauFakeRatePlots_HPScombined_%s.eps' % evtSelName
    retVal_makeTauFakeRatePlots = \
      buildConfigFile_makeTauFakeRatePlots(haddOutputFileName, eventSelections, evtSelJob, tauIds, discriminators_HPScombined, labels,
                                           outputFilePath, outputFileName_HPScombined, recoSampleDefinitionsTauIdCommissioning_7TeV)
    fileNames_makeTauFakeRatePlots.append(retVal_makeTauFakeRatePlots)
#--------------------------------------------------------------------------------

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

# done building config files, now build Makefile...
makeFileName = "Makefile_TauFakeRateAnalysis"
makeFile = open(makeFileName, "w")
makeFile.write("\n")
outputFileNames_makeTauFakeRatePlots = []
for fileNameEntry in fileNames_makeTauFakeRatePlots:
    outputFileNames_makeTauFakeRatePlots.append(fileNameEntry['outputFileName'])
makeFile.write("all: %s %s\n" %
  (haddOutputFileName,
   make_MakeFile_vstring(outputFileNames_makeTauFakeRatePlots)))
makeFile.write("\techo 'Finished running TauFakeRateAnalysis.'\n")
makeFile.write("\n")
for sampleToAnalyze in samplesToAnalyze:
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:
            if not skipFWLiteTauFakeRateAnalyzer:
                fileNameEntry = fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]
                if fileNameEntry is None or len(fileNameEntry['inputFileNames']) == 0:
                    continue
                bsubJobEntry = bsubJobNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]
                for i in range(len(fileNameEntry['inputFileNames'])):
                    makeFile.write("%s: %s\n" %
                      (fileNameEntry['outputFileNames'][i],
                       executable_FWLiteTauFakeRateAnalyzer))
                    makeFile.write("\t%s -q %s -J %s < %s\n" %
                      (executable_bsub,
                       bsubQueue,
                       bsubJobEntry[i],
                       fileNameEntry['bsubScriptFileNames'][i]))
            else:
                fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]['outputFileNames'] = []
            makeFile.write("\n")
            makeFile.write("%s: %s\n" %
              (bsubJobNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze],
               make_MakeFile_vstring(fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]['outputFileNames'])))
            if not skipFWLiteTauFakeRateAnalyzer:
                makeFile.write("\t%s %s\n" %
                  (executable_waitForLXBatchJobs,
                   bjobListFileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]))
            makeFile.write("\t%s %s\n" %
              (executable_shell,
               bsubFileNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze]['harvest_script_name']))
            makeFile.write("\n")
makeFile.write("%s: %s\n" %
  (haddOutputFileName,
   make_MakeFile_vstring(bsubJobNames_harvesting_all)))
makeFile.write("\t%s %s\n" %
  (executable_waitForLXBatchJobs,
   bjobListFileName_harvesting))
for haddInputFileName in haddInputFileNames:
    makeFile.write("\t%s %s %s\n" %
      (executable_rfcp,
       os.path.join(harvestingFilePath, os.path.basename(haddInputFileName)),
       outputFilePath))
makeFile.write("\t%s %s &> %s\n" %
  (executable_shell,
   haddShellFileName,
   haddLogFileName))
makeFile.write("\n")
for fileNameEntry in fileNames_makeTauFakeRatePlots:
    makeFile.write("%s: %s %s\n" %
      (fileNameEntry['outputFileName'],
       executable_makeTauFakeRatePlots,
       haddOutputFileName))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_makeTauFakeRatePlots,
       fileNameEntry['configFileName'],
       fileNameEntry['logFileName']))
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
for sampleToAnalyze in samplesToAnalyze:
    for eventSelectionToAnalyze in eventSelectionsToAnalyze:
        jobNameInRecoSampleDef = eventSelections[eventSelectionToAnalyze]['jobNameInRecoSampleDef']
        if jobNameInRecoSampleDef in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sampleToAnalyze]['jobs']:
            if not skipFWLiteTauFakeRateAnalyzer:
                fileNameEntry = fileNames_FWLiteTauFakeRateAnalyzer[sampleToAnalyze][eventSelectionToAnalyze]
                if fileNameEntry is None:
                    continue
                for outputFileName in fileNameEntry['outputFileNames']:
                    makeFile.write("\t%s %s\n" %
                      (executable_rfrm,
                       os.path.join(harvestingFilePath, outputFileName)))
            for final_harvest_file in bsubFileNames_harvesting[sampleToAnalyze][eventSelectionToAnalyze]['final_harvest_files']:
                # CV: file name of final harvesting output file is stored at index[1] in final_harvest_file-tuple
                #    (cf. TauAnalysis/Configuration/python/tools/harvestingLXBatch.py)    
                makeFile.write("\t%s %s\n" %
                  (executable_rfrm,
                   os.path.join(harvestingFilePath, final_harvest_file[1])))
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(haddInputFileNames))            
makeFile.write("\trm -f %s\n" % haddShellFileName)
makeFile.write("\trm -f %s\n" % haddOutputFileName)
for fileNameEntry in fileNames_makeTauFakeRatePlots:
    makeFile.write("\trm -f %s\n" % fileNameEntry['outputFileName'])
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)
