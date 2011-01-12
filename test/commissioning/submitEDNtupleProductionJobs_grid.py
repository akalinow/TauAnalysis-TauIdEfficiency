#!/usr/bin/env python

import os
import subprocess
import time

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdCommissioning_7TeV_grid_cfi import *

print("<submitEDNtupleProductionJobs_grid>:")

castorFilePath = "/castor/cern.ch/user/v/veelken/TauIdCommissioning/"
#crabFilePath = "/afs/cern.ch/user/v/veelken/scratch0/CMSSW_3_8_7/src/TauAnalysis/TauIdEfficiency/test/commissioning/crab/"
crabFilePath = "/data1/veelken/CMSSW_3_8_x/crab/TauIdEfficiency/"

pfCandidateCollection = "particleFlow" # pile-up removal disabled
#pfCandidateCollection = "pfNoPileUp"   # pile-up removal enabled

version = "v4_0q"

crab_template = """
[CRAB]

jobtype = cmssw
scheduler = glidein
use_server = 1

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

shFileName = os.path.join(crabFilePath, "submitEDNtupleProductionJobs_grid.csh")
shFile = open(shFileName, "w")
shFile.write("#!/bin/csh -f\n")
shFile.write("\n")

for sampleName in SAMPLES_TO_RUN:
     for jobType in JOBS_TO_RUN:
          if jobType in RECO_SAMPLES[sampleName].get('jobs'):

               print("--> submitting EDNtuple production job for sample = %s, job = %s..." % (sampleName, jobType))

               cfgFileName = CONFIG_FILES[jobType]
               cfgFile = open(cfgFileName, "r")
               cfg = cfgFile.read()
               cfgFile.close()

               cfg = cfg.replace("#__", "")
               isMC = "False"
               if RECO_SAMPLES[sampleName]['type'] == 'MC':
                    isMC = "True"
               cfg = cfg.replace("#isMC#", isMC)
               cfg = cfg.replace("#HLTprocessName#", "'%s'" % RECO_SAMPLES[sampleName]['hlt'])
               cfg = cfg.replace("#pfCandidateCollection#", "'%s'" % pfCandidateCollection)

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
               ui_working_dir = os.path.join(crabFilePath, "crabdirProduceEDNtuple_%s_%s_%s" % (sampleName, jobType, version))
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "UI_WORKING_DIR", 'ui_working_dir', ui_working_dir)
               user_remote_dir = os.path.join(castorFilePath, jobType) + '/' + version + '/' + sampleName
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "USER_REMOTE_DIR", 'user_remote_dir', user_remote_dir.replace("/castor/cern.ch", ""))
               crabConfig = setOption(crabConfig, RECO_SAMPLES[sampleName], "SE_BLACK_LIST", 'SE_black_list')
               crabFile = open(crabFileName, "w")
               crabFile.write(crabConfig)
               crabFile.close()

               shFile.write("#-------------------- %s, %s --------------------\n" % (sampleName, jobType))
               shFile.write("rfmkdir %s\n" % (os.path.join(castorFilePath, jobType)))
               shFile.write("rfchmod 777 %s\n" % (os.path.join(castorFilePath, jobType)))
               shFile.write("rfmkdir %s\n" % (os.path.join(castorFilePath, jobType) + '/' + version))
               shFile.write("rfchmod 777 %s\n" % (os.path.join(castorFilePath, jobType) + '/' + version))
               shFile.write("rfmkdir %s\n" % (os.path.join(castorFilePath, jobType) + '/' + version + '/' + sampleName))
               shFile.write("rfchmod 777 %s\n" % (os.path.join(castorFilePath, jobType) + '/' + version + '/' + sampleName))
               shFile.write("\n")
               shFile.write("mkdir -p %s\n" % ui_working_dir)
               shFile.write("crab -create -cfg %s\n" % crabFileName)
               shFile.write("crab -submit -c %s\n" % ui_working_dir)
               shFile.write("\n")

shFile.close()

subprocess.call("chmod +x %s" % shFileName, shell = True)

# need to wait until config files have finished writing...
time.sleep(1)

print("Finished building config files. Now execute 'source %s'." % shFileName)
#subprocess.call("source %s" % shFileName, shell = True)
