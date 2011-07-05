#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.submitAnalysisToLXBatch import submitAnalysisToLXBatch
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId, getBatchHarvestLocation, getHarvestingFilePath
import TauAnalysis.Configuration.tools.castor as castor

import os
import subprocess

channel = 'ZtoMuTau_tauIdEff'
configFile = 'produceTauIdEffMeasPATTuple_cfg.py'
analysisFilePath = getAnalysisFilePath(channel)
#jobId = getJobId(channel)
jobId = '2011Jul01_mauro'

version = "V2"

pfCandidateCollection = "particleFlow" # pile-up removal disabled
#pfCandidateCollection = "pfNoPileUp"   # pile-up removal enabled

samplesToAnalyze = [
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'Ztautau_pythia',
    'Zmumu_pythia',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

#outputFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/patTuples_v6"
outputFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/"

# Get all the skim files from the castor directory
skimFilePath = getBatchHarvestLocation(channel)
skim_files = [ file_info['path'] for file_info in castor.nslsl(skimFilePath) ]

if not os.path.isdir("lxbatch_pattuple"):
    print 'Creating directory to store the lxbatch jobs: lxbatch_pattuple'
    os.mkdir('lxbatch_pattuple')

if not os.path.isdir("lxbatch_pat_log"):
    print 'Creating directory to store the lxbatch logs: lxbatch_pat_log'
    os.mkdir('lxbatch_pat_log')


# Function that maps a sample name to its skim file
def input_mapper(channel, sample, jobId):
    for input_file in skim_files:
        #print " unmatched file: %s" % input_file
        if input_file.find('skim_' + sample + '_chunk') != -1:
            #print "--> matched file: %s" % input_file
            yield input_file

# Define what output ntuple file name a sample will have
def output_mapper(channel, sample, jobId):
    output_file = "tauIdEffMeasPATTuple_%s_%s%s.root" % (sample, jobId, version)
    return output_file

# Function to prepare customized config files specific to Tau(ED)Ntuple production
def customizeConfigFile(sampleName, cfgFileName_original, cfgFileName_modified = None):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    isMC = "False"
    if recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'Data':
        isMC = "True"
    cfg_modified = cfg_modified.replace("#isMC#", isMC)
    HLTprocessName = 'HLT'
    if 'hlt' in recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName].keys():
        HLTprocessName = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['hlt'].getProcessName()
    cfg_modified = cfg_modified.replace("#HLTprocessName#", "'%s'" % HLTprocessName)
    cfg_modified = cfg_modified.replace("#pfCandidateCollection#", "'%s'" % pfCandidateCollection)
    applyZrecoilCorrection = "False"
    if 'applyZrecoilCorrection' in recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName].keys() and \
      recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['applyZrecoilCorrection']:
        applyZrecoilCorrection = "True"                
    cfg_modified = cfg_modified.replace("#applyZrecoilCorrection#", "%s" % applyZrecoilCorrection)

    if cfgFileName_modified is None:
        cfgFileName_modified = cfgFileName_original.replace("_cfg.py", "_customized_%s_cfg.py" % sampleName)
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

    return cfgFileName_modified

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['SAMPLES_TO_ANALYZE']

shFileNames_modified = []

for sampleToAnalyze in samplesToAnalyze:

    # prepare customized config file as basis for further modifications by "TauAnalysis machinery"...
    configFile_customized = customizeConfigFile(sampleToAnalyze, configFile)

    # apply further modifications and submit job to lxbatch
    shFileName = \
      submitAnalysisToLXBatch(configFile = configFile_customized, channel = channel,
                              samples = recoSampleDefinitionsTauIdEfficiency_7TeV,
                              samplesToAnalyze = [ sampleToAnalyze ],
                              disableFactorization = True, disableSysUncertainties = True, disableZrecoilCorrections = True,
                              # Options for local running
                              cfgdir = 'lxbatch_pattuple', 
                              inputFileMap = input_mapper,
                              outputFileMap = output_mapper,
                              outputDirectory = outputFilePath,
                              processName = 'lxbatch',
                              saveFinalEvents = False,
                              jobExtention = "_PATTuple")

    # rename shell script and move to "./lxbatch" subdirectory
    shFileName_modified = shFileName.replace(jobId, "%s_%s" % (sampleToAnalyze, jobId))
    subprocess.call("mv %s lxbatch_pattuple/%s" % (shFileName, shFileName_modified), shell = True)
    shFileNames_modified.append(shFileName_modified)

    # move customized config file to "./lxbatch" subdirectory
    #
    # NOTE: "TauAnalysis machinery" does not work if customized config file
    #       is created in "./lxbatch" subdirectory from the start
    #
    subprocess.call("mv %s lxbatch_pattuple" % configFile_customized, shell = True)

# create "master" shell script
shFileName_master = "submit_lxbatch_PATTuple_analysis_all_%s.sh" % jobId
shFile_master = open(shFileName_master, "w")
for shFileName_modified in shFileNames_modified:
    shFile_master.write("source lxbatch_pattuple/%s\n" % shFileName_modified)
shFile_master.close()

os.system("./swapLogsToPatDir.sh")

print "\n"
print "Run ./%s to submit **all** jobs" % shFileName_master
os.chmod(shFileName_master, 0755)

