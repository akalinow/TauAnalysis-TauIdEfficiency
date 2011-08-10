#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.submitAnalysisToLXBatch import submitAnalysisToLXBatch
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId, getBatchHarvestLocation, getHarvestingFilePath
import TauAnalysis.Configuration.tools.castor as castor

import os
import subprocess

channel = 'ZtoMuTau_tauIdEff'
configFile = 'produceTauPtResPATTuple_cfg.py'
analysisFilePath = getAnalysisFilePath(channel)
jobId = '2011Jul23'

version = 'V1c'

samplesToAnalyze = [
    'Ztautau_powheg'
]

outputFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/PATtuples/TauPtRes/V1c/"

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
    output_file = "tauPtResPATtuple_%s_%s%s.root" % (sample, jobId, version)
    return output_file

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['SAMPLES_TO_ANALYZE']

shFileNames_modified = []

for sampleToAnalyze in samplesToAnalyze:

    # customize config file and submit job to lxbatch
    shFileName = \
      submitAnalysisToLXBatch(configFile = configFile, channel = channel,
                              samples = recoSampleDefinitionsTauIdEfficiency_7TeV,
                              samplesToAnalyze = [ sampleToAnalyze ],
                              disableFactorization = True, disableSysUncertainties = True, disableZrecoilCorrections = True,
                              # Options for local running
                              cfgdir = 'lxbatch_pattuple', 
                              inputFileMap = input_mapper,
                              outputFileMap = output_mapper,
                              outputDirectory = outputFilePath,
			      queue = '1nw',
                              processName = 'lxbatch',
                              saveFinalEvents = False,
                              jobExtention = "_PATTuple")

    # rename shell script and move to "./lxbatch" subdirectory
    shFileName_modified = shFileName.replace(jobId, "%s_%s" % (sampleToAnalyze, jobId))
    subprocess.call("mv %s lxbatch_pattuple/%s" % (shFileName, shFileName_modified), shell = True)
    shFileNames_modified.append(shFileName_modified)

    # copy config file to "./lxbatch" subdirectory
    subprocess.call("cp %s lxbatch_pattuple" % configFile, shell = True)

# create "master" shell script
shFileName_master = "submit_lxbatch_PATTuple_analysis_all_%s.sh" % jobId
shFile_master = open(shFileName_master, "w")
for shFileName_modified in shFileNames_modified:
    shFile_master.write("source lxbatch_pattuple/%s\n" % shFileName_modified)
shFile_master.close()

print "\n"
print "Run ./%s to submit **all** jobs" % shFileName_master
os.chmod(shFileName_master, 0755)

