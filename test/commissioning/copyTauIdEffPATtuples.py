#!/usr/bin/env python

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.TauIdEfficiency.tools.castor_mirror2 as castor_mirror

import os
import subprocess
import shlex

jobId = "2011Oct30" # Christian's Ntuples
version = "V10_1tauEnRecovery"
#jobId = "2011Jun06" # Mauro's Ntuples
#version = "V2"

# Get all the skim files from the castor directory
sourceFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/%s/%s" % (jobId, version) # Christian's PAT-tuples
#sourceFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/Jun06Skim/edntuples_v3/" # Mauro's Ntuples
source_files = [ file_info['path'] for file_info in castor.nslsl(sourceFilePath) ]
print "source_files:"
print source_files

targetFilePath = "/data1/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/%s/%s/" % (jobId, version)

if not os.path.exists(targetFilePath):
    os.mkdir(targetFilePath)

samplesToCopy = [
    # modify in case you want to copy some of the samples only...
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_Aug05ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v6',
    'data_MET_Run2011B_PromptReco_v1',
    'data_MET_Run2011B_PromptReco_v1a',
    'Ztautau_powheg',
    #'Ztautau_embedded_part1',
    #'Ztautau_embedded_part2',
    'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

files_to_copy = []

for source_file in source_files:

    if not (source_file.find(jobId) != -1 and source_file.find(version) != -1):
	continue

    isSampleToCopy = False
    if len(samplesToCopy) == 0:
        isSampleToCopy = True
    for sampleToCopy in samplesToCopy:
        if source_file.find(sampleToCopy) != -1:
            isSampleToCopy = True
    if not isSampleToCopy:
        continue;

    target_file = source_file.replace(sourceFilePath, targetFilePath)
    
    print("copying %s --> %s" % (source_file, target_file))
    files_to_copy.append(source_file)

castor_mirror.mirror_files(castor_mirror.needs_local_copy(files_to_copy, [ targetFilePath ]), [ targetFilePath ], 3)
