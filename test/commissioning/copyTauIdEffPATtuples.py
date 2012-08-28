#!/usr/bin/env python

import TauAnalysis.Configuration.tools.castor as castor
import TauAnalysis.TauIdEfficiency.tools.castor_mirror2 as castor_mirror

import os
import subprocess
import shlex

#jobId = getJobId(channel)
jobId = '2012May12'

version = 'v1_8'

# Get all the skim files from the castor directory
sourceFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_5_2_x/PATtuples/TauIdEffMeas/%s/%s" % (jobId, version)
source_files = [ file_info['path'] for file_info in castor.nslsl(sourceFilePath) ]
#print "source_files:"
#print source_files

targetFilePath = "/data1/veelken/CMSSW_5_2_x/PATtuples/TauIdEffMeas/%s/%s/" % (jobId, version)

def createFilePath_recursively(filePath):
    filePath_items = filePath.split('/')
    currentFilePath = "/"
    for filePath_item in filePath_items:
        currentFilePath = os.path.join(currentFilePath, filePath_item)
        if len(currentFilePath) <= 1:
            continue
        if not os.path.exists(currentFilePath):
            #sys.stdout.write("creating directory %s\n" % currentFilePath)
            os.mkdir(currentFilePath)

createFilePath_recursively(targetFilePath)

samplesToCopy = [
    # modify in case you want to copy some of the samples only...
    'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to193621',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs193752to194076v2',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194108to194479',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194790to195016',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195099to195947',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195948to196509',
    'ZplusJets_madgraph2',
    'Ztautau_pythia',
    'WplusJets_madgraph',
    'WplusJets_madgraph_extension',
    'PPmuXptGt20Mu15',      
    'TTplusJets_madgraph2'
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
    
    #print("mirroring %s --> %s" % (source_file, target_file))
    files_to_copy.append(source_file)

castor_mirror.mirror_files(castor_mirror.needs_local_copy(files_to_copy, [ targetFilePath ]), [ targetFilePath ], 1)
