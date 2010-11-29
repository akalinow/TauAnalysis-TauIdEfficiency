#!/usr/bin/env python

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_Nov22_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.submitAnalysisToLocal import submitAnalysisToLocal
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId, getPickEventsPath, getHarvestingFilePath

import os

channel = 'AHtoMuTau'
configFile = 'produceTauIdEffMeasNtuple_cfg.py'
analysisFilePath = getAnalysisFilePath(channel)
#jobId = getJobId(channel)
jobId = 'Run26'

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
]

logFilePath = "/data1/veelken/logs/"

##outputFilePath = getHarvestingFilePath(channel)
outputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas"

# Function that maps a sample name to its pickevents file
def local_sample_mapper(channel, sample, jobId):
    return os.path.join(
        getPickEventsPath(channel),
        "skim_%s_%s_%s*.root" % (channel, sample, jobId)
    )

# Define what output ntuple file name a sample will have
def output_mapper(channel, sample, jobId):
    output_file = os.path.join(
        outputFilePath,
        'local',
        "tauIdEffMeasEDNtuple_%s_%s.root" % (sample, jobId)
    )
    # Make the directory if it doesn't exist yet
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    return output_file

submitAnalysisToLocal(configFile = configFile, channel = channel, jobId = jobId,
                      samples = recoSampleDefinitionsZtoMuTau_7TeV,
                      samplesToAnalyze = samplesToAnalyze,
                      disableSysUncertainties = True,
                      # Options for local running
                      cfgdir = 'local', 
                      maxEvents = 10000, maxJobsConcurrently = 8, submit = False, 
		      logFilePath = logFilePath,
                      inputFileMap = local_sample_mapper,
                      outputFileMap = output_mapper,
                      processName = 'local',
                      saveFinalEvents = False)
