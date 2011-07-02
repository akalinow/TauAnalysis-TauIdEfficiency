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
    samplesToAnalyze = recoSampleDefinitionsTauIdEff_7TeV['MERGE_SAMPLES'].keys()

inputFiles = os.listdir(inputFilePath)

configFileNames = []

for sampleToAnalyze in samplesToAnalyze:
    processSearchStrings = "".join('_', recoSampleDefinitionsTauIdEff_7TeV['MERGE_SAMPLES'][sampleToAnalyze]['samples'], '_')
    jobIdSearchString = "".join('_', jobId, '_')

    inputFiles_process = []
    for inputFile in inputFiles:
        isPATtuple = (inputFile.find("tauIdEffMeasPATtuple") != -1)        
        isProcess_matched = False
        for processSearchString in processSearchStrings:
            if inputFile.find(processSearchString) != -1:
                isProcess_matched = True
        isJobId_matched = (inputFile.find(jobIdSearchString) != -1)
        
        if isPATtuple and isProcess_matched and isJobId_matched:
            inputFiles_process.append(inputFile)

    if len(inputFiles_process):
        print("process %s has not input files --> skipping !!" % sampleToAnalyze)
        continue

config = """
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(%s),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('analyzeTauIdEffHistograms_%.root')
)

process.tauIdEffAnalyzer = cms.PSet(
    process = cms.string('%s'),

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',
        'B',
        'B1',
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

    srcTrigger cms.InputTag('patTriggerEvent'),
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsForTauIdEffCumulative')
)
""" % (inputFiles_process, sampleToAnalyze, sampleToAnalyze)

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
