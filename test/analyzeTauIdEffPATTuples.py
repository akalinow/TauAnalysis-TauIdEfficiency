#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os
import subprocess

channel = 'ZtoMuTau_tauIdEff'
#jobId = getJobId(channel)
jobId = '2011Jul06_mauroV4'

inputFilePath = '/data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Jul06_mauro/V4/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/'
outputFilePath = '/data1/veelken/tmp/muonPtGt20/V4c/'

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'Ztautau_pythia',
    ##'Ztautau_embedded_part1',
    ##'Ztautau_embedded_part2',
    'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

sysUncertainties = [
    "sysTauJetEn", # needed for diTauVisMass/diTauVisMassFromJet
    "sysJetEnUp"   # needed for diTauMt
]

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'].keys()

inputFiles = os.listdir(inputFilePath)
#print(inputFiles)

configFileNames = []
outputFileNames = []

sysUncertainties_expanded = [ "CENTRAL_VALUE" ]
for sysUncertainty in sysUncertainties:
    sysUncertainties_expanded.append(sysUncertainty + "Up")
    sysUncertainties_expanded.append(sysUncertainty + "Down")
    
for sampleToAnalyze in samplesToAnalyze:

    isMC = False

    # check if inputFile is PAT-tuple and
    # matches sampleToAnalyze, jobId
    inputFiles_sample = []
    for inputFile in inputFiles:        
        if inputFile.find("tauIdEffMeasPATTuple") != -1 and \
           inputFile.find("".join(['_', sampleToAnalyze, '_'])) != -1 and \
           inputFile.find("".join(['_', jobId, '_'])) != -1:
            inputFiles_sample.append(os.path.join(inputFilePath, inputFile))

    #print(sampleToAnalyze)
    #print(inputFiles_sample)

    if len(inputFiles_sample) == 0:
        print("Sample %s has not input files --> skipping !!" % sampleToAnalyze)
        continue

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

    print("building config file for sample %s..." % sampleToAnalyze)

    processType = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleToAnalyze]['type']

    inputFiles_string = "'"
    for i, inputFile_sample in enumerate(inputFiles_sample):
        if i > 0:
            inputFiles_string += "', '"
        inputFiles_string += inputFile_sample
    inputFiles_string += "'"

    tauIds = []

    for sysUncertainty in sysUncertainties_expanded:

        outputFileName = None
        if sysUncertainty != "CENTRAL_VALUE":            
            outputFileName = 'analyzeTauIdEffHistograms_%s_%s_%s.root' % (sampleToAnalyze, sysUncertainty, jobId)
        else:
            outputFileName = 'analyzeTauIdEffHistograms_%s_%s.root' % (sampleToAnalyze, jobId)
        outputFileName_full = os.path.join(outputFilePath, outputFileName)

        srcMuTauPairs = None
        if sysUncertainty != "CENTRAL_VALUE":  
            srcMuTauPairs = composeModuleName([ 'selectedMuPFTauHPSpairsDzForTauIdEff', sysUncertainty, "cumulative" ])            
        else:
            srcMuTauPairs = 'selectedMuPFTauHPSpairsDzForTauIdEffCumulative'

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

        config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(%s),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.tauIdEffAnalyzer = cms.PSet(
    process = cms.string('%s'),
    type = cms.string('%s'),

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

    sysShift = cms.string('%s'),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5', 'HLT_IsoMu17_v6', 'HLT_IsoMu17_v8', 'HLT_IsoMu17_v9', 'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('%s'),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(%s),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
""" % (inputFiles_string, outputFileName_full,
       process_matched, processType, sysUncertainty, srcMuTauPairs, weights_string, allEvents_DBS, xSection, intLumiData)

        configFileName = None
        if sysUncertainty != "CENTRAL_VALUE":  
            configFileName = "analyzeTauIdEffPATtuple_%s_%s_cfg.py" % (sampleToAnalyze, sysUncertainty)
        else:
            configFileName = "analyzeTauIdEffPATtuple_%s_cfg.py" % sampleToAnalyze
        configFile = open(configFileName, "w")
        configFile.write(config)
        configFile.close()
        configFileNames.append(configFileName)

        outputFileNames.append(outputFileName_full)

executable = '../../../../bin/slc5_amd64_gcc434/FWLiteTauIdEffAnalyzer'

shellFileName = "analyzeTauIdEffPATtuples.csh"
shellFile = open(shellFileName, "w")
shellFile.write("#!/bin/csh -f\n")
shellFile.write("\n")
for configFileName in configFileNames:
    # CV: add '&' at end of command-line to run all FWLiteTauIdEffAnalyzer jobs in parallel
    logFileName = configFileName.replace("_cfg.py", ".log")
    shellFile.write('%s %s >&! %s &\n' % (executable, configFileName, logFileName))
shellFile.close()

print("Finished building config files. Now execute 'source %s'." % shellFileName)

haddFileName = "mergeTauIdEffHistograms.csh"
haddFile = open(haddFileName, "w")
haddFile.write("#!/bin/csh -f\n")
haddFile.write("\n")
haddOutputFileName = os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_all_%s.root' % jobId)
haddCommandLine = "hadd %s" % haddOutputFileName
for outputFileName in outputFileNames:
    haddCommandLine += " %s" % outputFileName
haddFile.write("%s\n" % haddCommandLine)
haddFile.close()

print("Once all FWLiteTauIdEffAnalyzer jobs have finished, execute 'source %s' to merge output files." % haddFileName)
