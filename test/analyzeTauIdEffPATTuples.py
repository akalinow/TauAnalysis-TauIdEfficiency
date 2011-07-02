#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId

import os
import subprocess

channel = 'ZtoMuTau_tauIdEff'
#jobId = getJobId(channel)
jobId = '2011Jul03'

inputFilePath = '/data1/veelken/tmp/'
outputFilePath = '/data1/veelken/tmp/'

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
]

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'].keys()

inputFiles = os.listdir(inputFilePath)
#print(inputFiles)

configFileNames = []

for sampleToAnalyze in samplesToAnalyze:

    isMC = False

    inputFiles_process = []
    for inputFile in inputFiles:        
        isPATtuple = (inputFile.find("tauIdEffMeasPATtuple") != -1)        
        isProcess_matched = False
        processes = recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][sampleToAnalyze]['samples']
        for process in processes:
            processSearchString = "".join(['_', process, '_'])
            if inputFile.find(processSearchString) != -1:
                isProcess_matched = True
        jobIdSearchString = "".join(['_', jobId, '_'])        
        isJobId_matched = (inputFile.find(jobIdSearchString) != -1)
        
        if isPATtuple and isProcess_matched and isJobId_matched:
            inputFiles_process.append(os.path.join(inputFilePath, inputFile))

    #print(sampleToAnalyze)
    #print(inputFiles_process)

    if len(inputFiles_process) == 0:
        print("process %s has not input files --> skipping !!" % sampleToAnalyze)
        continue

    print("building config file for process %s..." % sampleToAnalyze)

    inputFiles_string = "'"
    for i, inputFile_process in enumerate(inputFiles_process):
        if i > 0:
            inputFiles_string += "', '"
        inputFiles_string += inputFile_process
    inputFiles_string += "'"

    weights_string = ""
    if not recoSampleDefinitionsTauIdEfficiency_7TeV['MERGE_SAMPLES'][sampleToAnalyze]['type'] == 'Data':
        weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight", "'"])
        #weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#selectedPatMuonsForTauIdEffTrkIPcumulative#muonHLTeff", "'"])

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

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',
        'B',
        'B1',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D'
    ),
    
    tauIds = cms.VPSet(
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byLooseIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPSloose")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byMediumIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPSmedium")
        ),
        cms.PSet(
            discriminators = cms.vstring(
                'decayModeFinding',
                'byTightIsolationDeltaBetaCorr'
            ),
            name = cms.string("tauDiscrHPStight")
        )        
    ),

    sysShift = cms.string("CENTRAL_VALUE"),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5', 'HLT_IsoMu17_v6', 'HLT_IsoMu17_v8', 'HLT_IsoMu17_v9', 'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),

    weights = cms.VInputTag(%s)
)
""" % (inputFiles_string,
       os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_%s_%s.root' % (sampleToAnalyze, jobId)),
       sampleToAnalyze,
       weights_string)

    configFileName = "analyzeTauIdEffPATtuple_%s_cfg.py" % sampleToAnalyze
    configFile = open(configFileName, "w")
    configFile.write(config)
    configFile.close()
    configFileNames.append(configFileName)

executable = '../../../../bin/slc5_amd64_gcc434/FWLiteTauIdEffAnalyzer'

shellFileName = "analyzeTauIdEffPATtuples.csh"
shellFile = open(shellFileName, "w")
shellFile.write("#!/bin/csh -f\n")
shellFile.write("\n")
for configFileName in configFileNames:
    shellFile.write('%s %s' % (executable, configFileName))
shellFile.close()

print("Finished building config files. Now execute 'source %s'." % shellFileName)
