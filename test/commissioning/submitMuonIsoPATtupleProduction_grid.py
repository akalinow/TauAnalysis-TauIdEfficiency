#!/usr/bin/env python

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.submitAnalysisToGrid import submitAnalysisToGrid
from TauAnalysis.Configuration.userRegistry import getJobId

import os
import subprocess

channel = 'muonIsoPATtuple'
configFile = 'produceMuonIsolationPATtuple_cfg.py'
jobId = getJobId(channel)
outputFilePath = '/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/PATtuples/MuonIsolation/2011Oct10/v3/'

samplesToAnalyze = [
    # modify in case you want to submit crab jobs for some of the samples only...   
    'PPmuXptGt20Mu15'
]

# Define what output file name a skimmed sample will have
def output_mapper(channel, sample, jobId):
    output_file = "muonIsolationPATtuple_%s_%s.root" % (sample, jobId)
    return output_file

# Function to prepare customized config files specific to TauIdEff. skim 
def customizeConfigFile(sampleName, cfgFileName_original, cfgFileName_modified = None):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    isMC = "False"
    if recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'Data' and \
       recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'embeddedData':
        isMC = "True"
    cfg_modified = cfg_modified.replace("#isMC#", isMC)
    HLTprocessName = 'HLT'
    if 'hlt' in recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName].keys():
        HLTprocessName = recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['hlt'].getProcessName()
    cfg_modified = cfg_modified.replace("#HLTprocessName#", "'%s'" % HLTprocessName)

    if cfgFileName_modified is None:
        cfgFileName_modified = cfgFileName_original.replace("_cfg.py", "_customized_%s_cfg.py" % sampleName)
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

    return cfgFileName_modified

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsZtoMuTau_7TeV['SAMPLES_TO_ANALYZE']

for sampleToAnalyze in samplesToAnalyze:

    # prepare customized config file as basis for further modifications by "TauAnalysis machinery"...
    configFile_customized = customizeConfigFile(sampleToAnalyze, configFile)

    # apply further modifications and submit job to grid
    submitAnalysisToGrid(configFile = configFile_customized, channel = channel, jobId = jobId,
                         samples = recoSampleDefinitionsZtoMuTau_7TeV,
                         samplesToAnalyze = [ sampleToAnalyze ],
                         disableFactorization = True, disableSysUncertainties = True, disableZrecoilCorrections = True,
                         outputFilePath = outputFilePath, outputFileMap = output_mapper)





