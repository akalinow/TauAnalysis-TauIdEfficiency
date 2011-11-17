#!/usr/bin/env python

import os
import subprocess
import time

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdCommissioning_7TeV_grid_cfi import *
from TauAnalysis.Configuration.submitToGrid2 import submitToGrid

print("<submitCommissioningPATTtupleProductionJobs_grid>:")

castorFilePath = "/castor/cern.ch/user/v/veelken/TauIdCommissioning/"
#crabFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/crab/"
#crabFilePath = "/data1/veelken/CMSSW_4_2_x/crab/TauIdEfficiency/"
crabFilePath = "/tmp/veelken/crab/"

pfCandidateCollection = "particleFlow" # pile-up removal disabled
#pfCandidateCollection = "pfNoPileUp"   # pile-up removal enabled

version = 'patV2_0'

crab_template = """
[CRAB]

jobtype = cmssw
scheduler = glite
use_server = 0

[CMSSW]

datasetpath = DATASETPATH
dbs_url = DBS_URL
lumi_mask = LUMI_MASK
runselection = RUNSELECTION
pset = PSET
output_file = OUTPUT_FILE
total_number_of_events = -1
events_per_job = 25000
total_number_of_lumis = -1
lumis_per_job = LUMIS_PER_JOB

[USER]

ui_working_dir = UI_WORKING_DIR
return_data = 0
copy_data = 1
storage_element = srm-cms.cern.ch
storage_path = /srm/managerv2?SFN=/castor/cern.ch
user_remote_dir = USER_REMOTE_DIR

[GRID]

SE_black_list = SE_BLACK_LIST
"""

def setOption(crabConfig, sample, placeholder, option, optionValue = None):
     if optionValue is None and sample.get(option) is not None:
          optionValue = sample[option]
     if optionValue is not None:
          return crabConfig.replace(placeholder, optionValue)
     else:
          return crabConfig.replace("%s = %s" % (option, placeholder), "")

shFileName = os.path.join(crabFilePath, "submitCommissioningPATTupleProductionJobs_grid.csh")
shFile = open(shFileName, "w")
shFile.write("#!/bin/csh -f\n")
shFile.write("\n")

for sampleName in SAMPLES_TO_RUN:
     for jobType in JOBS_TO_RUN:
          if jobType in RECO_SAMPLES[sampleName].get('jobs'):

               if RECO_SAMPLES[sampleName]['type'] in JOB_OPTIONS[jobType]['submitTypes']:               
                    print("--> submitting PAT-tuple production job for sample = %s, job = %s..." % (sampleName, jobType))
               else:
                    print("--> skipping submission of PAT-tuple production job for sample = %s, job = %s..." % (sampleName, jobType))
                    continue

               cfgFileName = CONFIG_FILES[jobType]
               cfgFile = open(cfgFileName, "r")
               cfg = cfgFile.read()
               cfgFile.close()

               cfg = cfg.replace("#__", "")
               isMC = "False"
               if RECO_SAMPLES[sampleName]['type'] == 'smMC':
                    isMC = "True"
               cfg = cfg.replace("#isMC#", isMC)
               cfg = cfg.replace("#HLTprocessName#", "'%s'" % RECO_SAMPLES[sampleName]['hlt'])
               cfg = cfg.replace("#pfCandidateCollection#", "'%s'" % pfCandidateCollection)
               applyEventSelection = "True"
               if not JOB_OPTIONS[jobType]['applyEventSelection']:
                    applyEventSelection = "False"
               cfg = cfg.replace("#applyEventSelection#", applyEventSelection)
               
               cfgFileName_modified = \
                 os.path.join(crabFilePath, cfgFileName.replace("_cfg.py", "_%s_%s_%s_cfg.py" % (sampleName, jobType, version)))
               cfgFile_modified = open(cfgFileName_modified, "w")
               cfgFile_modified.write(cfg)
               cfgFile_modified.close()

               crabFileName = crabFilePath + "crab_%s_%s_%s.cfg" % (sampleName, jobType, version)
               crabConfig = str(crab_template)
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "DATASETPATH", 'datasetpath')
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "DBS_URL", 'dbs_url')
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "LUMI_MASK", 'lumi_mask')
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "RUNSELECTION", 'runselection')
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "PSET", 'pset', cfgFileName_modified)
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "OUTPUT_FILE", 'output_file', ROOT_FILE_NAMES[jobType])
               print(" isMC = %s" % isMC)
               if isMC == "False":
                    crabConfig = crabConfig.replace("total_number_of_events = -1", "")
                    crabConfig = crabConfig.replace("events_per_job = 25000", "")
                    crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "LUMIS_PER_JOB", 'lumis_per_job') 
               else:
                    crabConfig = crabConfig.replace("total_number_of_lumis = -1", "")
                    crabConfig = crabConfig.replace("lumis_per_job = LUMIS_PER_JOB", "")
               ui_working_dir = os.path.join(crabFilePath, "crabdirProduceFakeRatePATtuple_%s_%s_%s" % (sampleName, jobType, version))
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "UI_WORKING_DIR", 'ui_working_dir', ui_working_dir)
               user_remote_dir = os.path.join(castorFilePath, jobType) + '/' + version + '/' + sampleName
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "USER_REMOTE_DIR", 'user_remote_dir', user_remote_dir.replace("/castor/cern.ch", ""))
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "SE_BLACK_LIST", 'SE_black_list')
               crabFile = open(crabFileName, "w")
               crabFile.write(crabConfig)
               crabFile.close()
               
               # create directory structure for output files
               def createDir(outputFilePath):
                    subprocess.call("rfmkdir %s"     % outputFilePath, shell = True)
                    subprocess.call("rfchmod 777 %s" % outputFilePath, shell = True)
               
               createDir(os.path.join(castorFilePath, jobType))
               createDir(os.path.join(castorFilePath, jobType) + '/' + version)
               createDir(os.path.join(castorFilePath, jobType) + '/' + version + '/' + sampleName)

               submitToGrid(configFile = None, jobInfo = None, crabOptions = None,
                            crabFileName_full = crabFileName, ui_working_dir = ui_working_dir)
