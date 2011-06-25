#!/usr/bin/env python

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.TauIdEfficiency.tools.castor_mirror2 as castor_mirror

import subprocess
import shlex

# Get all the skim files from the castor directory
sourceFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/" # Christian's Ntuples
#sourceFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/Jun06Skim/edntuples_v3/" # Mauro's Ntuples
source_files = [ file_info['path'] for file_info in castor.nslsl(sourceFilePath) ]
#print "source_files:"
#print source_files

#targetFilePath = "/tmp/veelken/"
targetFilePath = "/data2/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/2011Jun18/"

jobId = "2011Jun18" # Christian's Ntuples
version = "V1"
#jobId = "2011Jun06" # Mauro's Ntuples
#version = "V2"

samplesToCopy = [
    # modify in case you want to submit jobs for some of the samples only...
]

files_to_copy = []

for source_file in source_files:

    if source_file.find("%s%s" % (jobId, version)) == -1:
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

castor_mirror.mirror_files(castor_mirror.needs_local_copy(files_to_copy, [ targetFilePath ]), [ targetFilePath ], 50)
