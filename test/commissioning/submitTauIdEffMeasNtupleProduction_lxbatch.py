#!/usr/bin/env python

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.submitAnalysisToLXBatch import submitAnalysisToLXBatch
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId, getPickEventsPath, getHarvestingFilePath
import TauAnalysis.Configuration.tools.castor as castor

import os
import subprocess

channel = 'ZtoMuTau_tauIdEff'
configFile = 'produceTauIdEffMeasNtuple_cfg.py'
analysisFilePath = getAnalysisFilePath(channel)
#jobId = getJobId(channel)
jobId = '2011Apr05'

version = "V1"

pfCandidateCollection = "particleFlow" # pile-up removal disabled
#pfCandidateCollection = "pfNoPileUp"   # pile-up removal enabled

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
]

outputFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas"

# Get all the skim files from the castor directory
skimFilePath = getPickEventsPath(channel)
skim_files = [ file_info['path'] for file_info in castor.nslsl(skimFilePath) ]

# Function that maps a sample name to its skim file
def input_mapper(channel, sample, jobId):
    for input_file in skim_files:
        if input_file.find('skim_' + sample + '_chunk') != -1:
            yield input_file

# Define what output ntuple file name a sample will have
def output_mapper(channel, sample, jobId):
    output_file = "tauIdEffMeasEDNtuple_%s_%s%s.root" % (sample, jobId, version)
    return output_file

# Function to prepare customized config files specific to Tau(ED)Ntuple production
def customizeConfigFile(sampleName, cfgFileName_original, cfgFileName_modified = None):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    isMC = "False"
    if recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'Data':
        isMC = "True"
    cfg_modified = cfg_modified.replace("#isMC#", isMC)
    HLTprocessName = 'HLT'
    if 'hlt' in recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName].keys():
        HLTprocessName = recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['hlt'].getProcessName()
    cfg_modified = cfg_modified.replace("#HLTprocessName#", "'%s'" % HLTprocessName)
    cfg_modified = cfg_modified.replace("#pfCandidateCollection#", "'%s'" % pfCandidateCollection)
    applyZrecoilCorrection = "False"
    if 'applyZrecoilCorrection' in recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName].keys() and \
      recoSampleDefinitionsZtoMuTau_7TeV['RECO_SAMPLES'][sampleName]['applyZrecoilCorrection']:
        applyZrecoilCorrection = "True"                
    cfg_modified = cfg_modified.replace("#applyZrecoilCorrection#", "%s" % applyZrecoilCorrection)

    if cfgFileName_modified is None:
        cfgFileName_modified = cfgFileName_original.replace("_cfg.py", "_customized_%s_cfg.py" % sampleName)
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

    return cfgFileName_modified

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsZtoMuTau_7TeV['SAMPLES_TO_ANALYZE']

shFileNames_modified = []

for sampleToAnalyze in samplesToAnalyze:

    # prepare customized config file as basis for further modifications by "TauAnalysis machinery"...
    configFile_customized = customizeConfigFile(sampleToAnalyze, configFile)

    # apply further modifications and submit job to lxbatch
    shFileName = \
      submitAnalysisToLXBatch(configFile = configFile_customized, channel = channel,
                              samples = recoSampleDefinitionsZtoMuTau_7TeV,
                              samplesToAnalyze = [ sampleToAnalyze ],
                              disableFactorization = True, disableSysUncertainties = True, disableZrecoilCorrections = True,
                              # Options for local running
                              cfgdir = 'lxbatch', 
                              inputFileMap = input_mapper,
                              outputFileMap = output_mapper,
                              outputDirectory = outputFilePath,
                              processName = 'lxbatch',
                              saveFinalEvents = False)

    # rename shell script and move to "./lxbatch" subdirectory
    shFileName_modified = shFileName.replace(jobId, "%s_%s" % (sampleToAnalyze, jobId))
    subprocess.call("mv %s lxbatch/%s" % (shFileName, shFileName_modified), shell = True)
    shFileNames_modified.append(shFileName_modified)

    # move customized config file to "./lxbatch" subdirectory
    #
    # NOTE: "TauAnalysis machinery" does not work if customized config file
    #       is created in "./lxbatch" subdirectory from the start
    #
    subprocess.call("mv %s lxbatch" % configFile_customized, shell = True)

# create "master" shell script
shFileName_master = "submit_lxbatch_analysis_all_%s.sh" % jobId
shFile_master = open(shFileName_master, "w")
for shFileName_modified in shFileNames_modified:
    shFile_master.write("source lxbatch/%s\n" % shFileName_modified)
shFile_master.close()

print "\n"
print "Run ./%s to submit **all** jobs" % shFileName_master
os.chmod(shFileName_master, 0755)
