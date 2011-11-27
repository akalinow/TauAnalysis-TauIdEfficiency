#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.tools.jobtools import make_bsub_script
from TauAnalysis.Configuration.tools.harvestingLXBatch import make_harvest_scripts
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId, getBatchHarvestLocation, getHarvestingFilePath
import TauAnalysis.Configuration.tools.castor as castor

import os
import re
import subprocess

channel = 'ZtoMuTau_tauIdEff'

#jobId = getJobId(channel)
jobId = '2011Oct30'

version = 'V10_4tauEnRecovery'

lxbatch_queue = '1nw'

pfCandidateCollection = "particleFlow" # pile-up removal disabled
#pfCandidateCollection = "pfNoPileUp"   # pile-up removal enabled

samplesToAnalyze = [
    #'data_SingleMu_Run2011A_May10ReReco_v1',
    #'data_SingleMu_Run2011A_PromptReco_v4',
    #'data_SingleMu_Run2011A_Aug05ReReco_v1',
    #'data_SingleMu_Run2011A_PromptReco_v6',
    #'data_MET_Run2011B_PromptReco_v1s1',
    'Ztautau_powheg',
    #'Ztautau_embedded_Run2011A_May10ReReco',
    #'Ztautau_embedded_Run2011A_PromptReco_v4',
    #'Ztautau_embedded_Run2011A_Aug05ReReco_v1',
    #'Ztautau_embedded_Run2011A_PromptReco_v6',
    #'Ztautau_embedded_Run2011B_PromptReco_v1',
    #'Zmumu_powheg',
    #'PPmuXptGt20Mu15',
    #'WplusJets_madgraph',
    #'TTplusJets_madgraph'
]

samplesToAnalyze_noTauSel = [
    'Ztautau_powheg'
]    

numInputFilesPerJob = {
    'data_SingleMu_Run2011A_May10ReReco_v1'    : 10,
    'data_SingleMu_Run2011A_PromptReco_v4'     : 10,
    'data_SingleMu_Run2011A_Aug05ReReco_v1'    :  5,
    'data_SingleMu_Run2011A_PromptReco_v6'     :  5,
    'data_MET_Run2011B_PromptReco_v1s1'        :  3,
    'Ztautau_powheg'                           :  3,
    'Ztautau_embedded_Run2011A_May10ReReco'    :  5,
    'Ztautau_embedded_Run2011A_PromptReco_v4'  :  5,
    'Ztautau_embedded_Run2011A_Aug05ReReco_v1' :  5,
    'Ztautau_embedded_Run2011A_PromptReco_v6'  :  5,
    'Ztautau_embedded_Run2011B_PromptReco_v1'  :  5,
    'Zmumu_powheg'                             :  5,
    'PPmuXptGt20Mu15'                          :  5,
    'WplusJets_madgraph'                       : 10,
    'TTplusJets_madgraph'                      :  3
}    

# Get all the skim files from the castor directory
inputFilePath = getAnalysisFilePath(channel)
print("inputFilePath = %s" % inputFilePath)

skipExistingPATtuples = True

#outputFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/patTuples_v6"
#outputFilePath = "/castor/cern.ch/user/m/mverzett/tagprobe/"
outputFilePath = "/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/%s/%s" % (jobId, version)

configFile = 'produceTauIdEffMeasPATTuple_cfg.py'
configFile_noTauSel = 'produceTauIdEffMeasPATTuple_noTauSel_cfg.py'

executable_bsub = 'bsub'

if not os.path.isdir("lxbatch"):
    os.mkdir('lxbatch')

configFilePath = os.path.join(os.getcwd(), "lxbatch")

if not os.path.isdir("lxbatch_log"):
    os.mkdir('lxbatch_log')

logFilePath = os.path.join(os.getcwd(), "lxbatch_log")

# Function that maps a sample name to its skim file
def input_mapper(inputFileNames, sample, jobId):
    inputFile_regex = \
      r"[a-zA-Z0-9_./]*tauIdEffSample_%s_%s_RECO_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sample, jobId)
    inputFile_matcher = re.compile(inputFile_regex)
    for input_file in inputFileNames:        
        #print "sample = %s:  unmatched file: %s" % (sample, input_file)
        if inputFile_matcher.match(input_file):
            #print "sample = %s: --> matched file: %s" % (sample, input_file)
            yield "".join(["rfio:", input_file])

# Define what output ntuple file name a sample will have
def output_mapper(sample, suffix, version, jobNumber):
    output_file = "tauIdEffMeasPATTuple_%s%s_%s_%s.root" % (sample, suffix, version, jobNumber)
    return output_file

# Split list of inputFileNames into groups of size N
def chunks(inputFileNames, N):
    for idx in range(0, len(inputFileNames), N):
        yield inputFileNames[idx:idx + N]

# Function to prepare customized config files specific to Tau(ED)Ntuple production
def customizeConfigFile(sampleName, jobId, version, inputFileNames, outputFileName, cfgFileName_original, cfgFileName_modified = None):
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    isMC = "False"
    if recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'Data' and \
       recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['type'] != 'embeddedData':
        isMC = "True"
    cfg_modified = cfg_modified.replace("#isMC#", isMC)
    isEmbedded = "False"
    if recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['type'] == 'embeddedData':
        isEmbedded = "True"
    cfg_modified = cfg_modified.replace("#isEmbedded#", isEmbedded)
    HLTprocessName = 'HLT'
    if 'hlt' in recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName].keys():
        HLTprocessName = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'][sampleName]['hlt'].getProcessName()
    cfg_modified = cfg_modified.replace("#HLTprocessName#", "'%s'" % HLTprocessName)
    cfg_modified = cfg_modified.replace("#pfCandidateCollection#", "'%s'" % pfCandidateCollection)

    cfg_modified += "\n"
    cfg_modified += "process.source.fileNames = cms.untracked.vstring(%s)\n" % \
                      [ "".join([ "file:", inputFileName ]) for inputFileName in inputFileNames ]
    cfg_modified += "process.patTupleOutputModule.fileName = cms.untracked.string('%s')\n" % outputFileName
    
    if cfgFileName_modified is None:
        cfgFileName_modified = cfgFileName_original.replace("_cfg.py", "_customized_%s_cfg.py" % sampleName)
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

    return cfgFileName_modified

if len(samplesToAnalyze) == 0 and len(samplesToAnalyze_noTauSel) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['SAMPLES_TO_ANALYZE']
samplesToAnalyze_all = []
samplesToAnalyze_all.extend(samplesToAnalyze)
samplesToAnalyze_all.extend(samplesToAnalyze_noTauSel)

bsubFileNames       = {}
bsubScriptFileNames = {}
bsubJobNames        = {}

for sampleToAnalyze in samplesToAnalyze_all:

    print "checking sample %s" % sampleToAnalyze

    bsubFileNames[sampleToAnalyze]       = {}
    bsubScriptFileNames[sampleToAnalyze] = {}
    bsubJobNames[sampleToAnalyze]        = {}

    inputFileNames = [ file_info['path'] for file_info in castor.nslsl(inputFilePath) ]
    #print " inputFileNames = %s" % inputFileNames
    
    inputFileNames_matched = [ os.path.basename(input_file) for input_file in input_mapper(inputFileNames, sampleToAnalyze, jobId) ]
    #print "inputFileNames_matched = %s" % inputFileNames_matched
    print "--> found %i inputFiles" % len(inputFileNames_matched)

    for jobNumber, inputFileNames_chunk in enumerate(chunks(inputFileNames_matched, numInputFilesPerJob[sampleToAnalyze])):
        # Build script for batch job submission;
        # the None in the tuple indicates that batch job has no dependencies on other batch jobs
        input_files_and_jobs = \
          [ (None, os.path.join(inputFilePath, inputFileName)) for inputFileName in inputFileNames_chunk ]

        for jobType in [ "", "_noTauSel" ]:
            configFile_original = None
            if jobType == "":
                if sampleToAnalyze in samplesToAnalyze:
                    configFile_original = configFile
            elif jobType == "_noTauSel":
                if sampleToAnalyze in samplesToAnalyze_noTauSel:
                    configFile_original = configFile_noTauSel
            else:
                raise ValueError("Invalid jobType = %s !!" % jobType)

            #print "jobType = %s: configFile_original = %s" % (jobType, configFile_original)
            
            if configFile_original is None:
                continue

            outputFileName = output_mapper(sampleToAnalyze, jobType, version, jobNumber)
            #print "outputFileName = %s" % outputFileName
            
            configFileName = os.path.join(
                configFilePath,
                configFile_original.replace("_cfg.py", "_%s%s_%s_cfg.py" % (sampleToAnalyze, jobType, jobNumber)))
            #print " configFileName = %s" % configFileName
            customizeConfigFile(sampleToAnalyze, jobNumber, version,
                                inputFileNames_chunk, outputFileName, configFile_original, configFileName)

            logFileName = os.path.basename(configFileName.replace('_cfg.py', '.log'))
            #print " logFileName = %s" % logFileName
            def log_file_maker(job_hash):
                return os.path.join(logFilePath, logFileName)

            bsubId = "%s_%i" % (jobType, jobNumber)

            jobName, bsubScript = make_bsub_script(
                os.path.join(outputFilePath, outputFileName),
                input_files_and_jobs,
                log_file_maker,
                "cmsRun %s" % os.path.join(configFilePath, configFileName))
            bsubFileNames[sampleToAnalyze][bsubId] = [ outputFileName ]

            bsubScriptFileName = os.path.join(configFilePath, logFileName.replace(".log", ".sh"))
            bsubScriptFile = open(bsubScriptFileName, "w")
            bsubScriptFile.write(bsubScript)
            bsubScriptFile.close()

            bsubScriptFileNames[sampleToAnalyze][bsubId] = bsubScriptFileName
            
            bsubJobName = "tauIdEffPATtuple_%s%s" % (sampleToAnalyze, bsubId)
            bsubJobNames[sampleToAnalyze][bsubId] = bsubJobName

# create "master" shell script
shFileName_master = "submitTauIdPATTupleProduction_lxbatch_%s_%s.sh" % (jobId, version)
shFile_master = open(shFileName_master, "w")
existingOutputFiles = [ file_info['file'] for file_info in castor.nslsl(outputFilePath) ]
for sampleToAnalyze in samplesToAnalyze_all:
    for bsubId in bsubScriptFileNames[sampleToAnalyze].keys():
        outputFilesExist = True
        for outputFileName in bsubFileNames[sampleToAnalyze][bsubId]:
            if not outputFileName in existingOutputFiles:
                outputFilesExist = False
        if skipExistingPATtuples and outputFilesExist:
            print "Output files for sample = %s, bsubId = %s exist --> skipping !!" % (sampleToAnalyze, bsubId)
        else:
            shFile_master.write("%s -q %s -J %s < %s\n" %
              (executable_bsub,
               lxbatch_queue,
               bsubJobNames[sampleToAnalyze][bsubId],
               bsubScriptFileNames[sampleToAnalyze][bsubId]))
    shFile_master.write("\n")
shFile_master.close()

print "\n"
print "Run ./%s to submit **all** jobs" % shFileName_master
os.chmod(shFileName_master, 0755)

