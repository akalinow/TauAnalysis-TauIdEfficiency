#!/usr/bin/env python

import os
import re

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdCommissioning_7TeV_grid_cfi import recoSampleDefinitionsTauIdCommissioning_7TeV
from TauAnalysis.Configuration.tools.harvestingLXBatch import make_harvest_scripts
from TauAnalysis.Configuration.tools.harvesting import castor_source, clean_by_crab_id

castorFilePath = '/castor/cern.ch/user/v/veelken/TauIdCommissioning/'
version = 'patV2_1'

SAMPLES_TO_ANALYZE = [
    # modify in case you want to submit jobs for some of the samples only...
    'data_Jet_Run2011A_May10ReReco_v1',
    'data_Jet_Run2011A_PromptReco_v4',
    'data_Jet_Aug05ReReco_v1',
    'data_Jet_Run2011A_PromptReco_v6', 
    'data_Jet_Run2011B_PromptReco_v1a', 
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_Aug05ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v6',
    'data_SingleMu_Run2011B_PromptReco_v1a',
    'ZplusJets',
    'WplusJets',
    'qcdDiJetPtHat15to30s4',
    'qcdDiJetPtHat30to50s4',
    'qcdDiJetPtHat50to80s4',
    'qcdDiJetPtHat80to120s4',
    'qcdDiJetPtHat120to170s4',
    'qcdDiJetPtHat170to300s4',
    'qcdDiJetPtHat300to470s4',
    'PPmuXptGt20Mu15',
    'TTplusJets'
]

JOBS_TO_ANALYZE = [
    # modify in case you want to submit jobs for some of the event selections:
    #  o QCD multi-jet:     'qcdDiJet'
    #  o QCD muon enriched: 'qcdMuEnriched'
    #  o W + jets:          'WplusJets'
    #  o Z --> mu+ mu-:     'Zmumu'
    # only...
]   

if len(SAMPLES_TO_ANALYZE) == 0:
    SAMPLES_TO_ANALYZE = recoSampleDefinitionsTauIdCommissioning_7TeV['SAMPLES_TO_RUN']
if len(JOBS_TO_ANALYZE) == 0:
    JOBS_TO_ANALYZE = recoSampleDefinitionsTauIdCommissioning_7TeV['JOBS_TO_RUN']

logFilePath = os.path.join(os.getcwd(), "lxbatch")

harvest_scripts = []
    
for sample in SAMPLES_TO_ANALYZE:
    for evtSel in JOBS_TO_ANALYZE:

        # CV: skip combinations of sample & evtSel that were not run
        #    (to avoid exception from castor that inputFilePath does not exists)
        if not evtSel in recoSampleDefinitionsTauIdCommissioning_7TeV['RECO_SAMPLES'][sample]['jobs']:
            continue
                
        inputFilePath = os.path.join(castorFilePath, evtSel, version, sample) + '/' # CV: add trailing '/'
        outputFilePath = inputFilePath
        print "harvesting files in inputFilePath = %s," % inputFilePath \
             + " copying harvested files to outputFilePath = %s..." % outputFilePath
       
        plot_regex = r"dont match anything"
        skim_regex = r"%s" % recoSampleDefinitionsTauIdCommissioning_7TeV['ROOT_FILE_NAMES'][evtSel].replace(
            ".root", "_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root")
        
        def matches_either(files):
            # Check if the file matches either of the regexes we are interested in.
            # We do this to skip extra files in the directories before we pass them to
            # clean_by_crab_id
            skim_matcher = re.compile(skim_regex)
            for file in files:
                 #print " unmatched file: %s" % file['path']
                 if skim_matcher.match(file['file']):
                     #print "--> matched file: %s" % file['path']
                     yield file  

        def local_copy_mapper(sample):
            # Define where we want to copy the final output locally 
            return os.path.join(
                outputFilePath,
                recoSampleDefinitionsTauIdCommissioning_7TeV['ROOT_FILE_NAMES'][evtSel].replace(".root" + "_harvested.root"))

        retVal_make_harvest_scripts = make_harvest_scripts(
            plot_regex,
            skim_regex,
            sampleToAnalyze = sample,
            job_id = "_".join([evtSel, version]),
            input_source = matches_either(castor_source(inputFilePath)),
            castor_output_directory = outputFilePath,
            script_directory = logFilePath,
            merge_script_name = os.path.join("lxbatch", "_".join(['submit', sample, evtSel, version, 'merge']) + '.sh'),
            local_copy_mapper = local_copy_mapper,
            chunk_size = 2.e+9, # 2 GB
            run_harvesting = False,
            check_old_files = False,
            max_bsub_concurrent_file_access = 500,
            verbosity = 0
        )

        harvest_scripts.append(retVal_make_harvest_scripts['merge_script_name'])

# create "master" shell script
shFileName_master = os.path.join("lxbatch", "submit_lxbatch_harvestCommissioningPATTuples_all_%s.sh" % version)
shFile_master = open(shFileName_master, "w")
for harvest_script in harvest_scripts:
    shFile_master.write("source %s\n" % harvest_script)
shFile_master.close()

print "\n"
print "Run ./%s to submit **all** jobs" % shFileName_master
os.chmod(shFileName_master, 0755)
