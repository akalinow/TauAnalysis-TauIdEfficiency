import FWCore.ParameterSet.Config as cms
import copy
import itertools

import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as styles

import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as ZtoMuTau

# List of samples to run in the analysis
SAMPLES_TO_ANALYZE = [
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_Aug05ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v6',
    #'data_MET_Run2011B_PromptReco_v1',
    #'data_MET_Run2011B_PromptReco_v1a',
    'data_MET_Run2011B_PromptReco_v1c',
    #'DYtautauM10to20_powheg',
    #'Ztautau_pythia',
    'Ztautau_powheg',
    'Ztautau_embedded_Run2011A_May10ReReco',
    'Ztautau_embedded_Run2011A_PromptReco_v4',
    'Ztautau_embedded_Run2011A_Aug05ReReco_v1',
    'Ztautau_embedded_Run2011A_PromptReco_v6',
    'Ztautau_embedded_Run2011B_PromptReco_v1',
    #'qqZll',
    #'DYmumuM10to20_pythia',
    #'Zmumu_pythia',
    'Zmumu_powheg',
    #'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    #'WW',
    #'WZ',
    #'ZZ',
    'TTplusJets_madgraph'
]

# Conversions to pico barns
_millibarns = 1.0e+9
_microbarns = 1.0e+6
_nanobarns  = 1.0e+3
_picobarns =  1.0
_femtobarns = 1.0e-3

# Integrated luminosity to normalize
# (computed by lumiCalc)
TARGET_LUMI = (
       33.60 # HLT_IsoMu17_v5  (160431-163261)
  +   168.61 # HLT_IsoMu17_v6  (163270-163869)
  +     0.00 # HLT_IsoMu17_v7  (never used ?)
  +   139.03 # HLT_IsoMu17_v8  (165088-165633)
  +   542.74 # HLT_IsoMu17_v9  (165970-167043)
  +     4.29 # HLT_IsoMu17_v10 (166346, used only for one run)
  +   243.01 # HLT_IsoMu17_v11 (167078-167913)
  +     0.00 # HLT_IsoMu17_v12 (never used ?)
  +   391.54 # HLT_IsoMu17_v13 (170826-173198)
##  +     3.39 # HLT_IsoMu17_v14 (173236-177053, WARNING: HLT_IsoMu17 trigger got prescaled !!)  
)/_picobarns

RECO_SAMPLES = copy.deepcopy(ZtoMuTau.RECO_SAMPLES)
TauIdEfficiencySpecific_RECO_SAMPLES = {
    'data_SingleMu_Run2011A_May10ReReco_v1' : {
        'datasetpath' : '/SingleMu/Run2011A-May10ReReco-v1/AOD',
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt",
        'runselection' : "160329-163869",
        'number_of_jobs' : 750,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_Mu15_v2'    : '160431:MIN-163261:MAX',
            'HLT_Mu15_v3'    : '163270:MIN-163869:MAX',
            'HLT_IsoMu17_v5' : '160431:MIN-163261:MAX',
            'HLT_IsoMu17_v6' : '163270:MIN-163869:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_SingleMu_Run2011A_PromptReco_v4' : {
        'datasetpath' : "/SingleMu/Run2011A-PromptReco-v4/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "165071-167913",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_Mu15_v4'     : '165088:MIN-167043:MAX',
            'HLT_Mu15_v5'     : '167078:MIN-167913:MAX',
            'HLT_IsoMu17_v8'  : '165088:MIN-165633:MAX',
            'HLT_IsoMu17_v9'  : '165970:MIN-167043:MAX',
            'HLT_IsoMu17_v10' : '166346:MIN-166346:MAX',
            'HLT_IsoMu17_v11' : '167078:MIN-167913:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_SingleMu_Run2011A_Aug05ReReco_v1' : {
        'datasetpath' : "/SingleMu/Run2011A-05Aug2011-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt",
        'runselection' : "170053-172619",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu17_v13' : '170826:MIN-173198:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_SingleMu_Run2011A_PromptReco_v6' : {
        'datasetpath' : "/SingleMu/Run2011A-PromptReco-v6/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-176023_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "172620-175770",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu17_v13' : '170826:MIN-173198:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_MET_Run2011B_PromptReco_v1' : {
        'datasetpath' : "/MET/Run2011B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-179431_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "175832-179411",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_Mu15_L1ETM20_v3'    : '178420:MIN-179411:MAX',
            'HLT_IsoMu15_L1ETM20_v3' : '178420:MIN-179411:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_MET_Run2011B_PromptReco_v1a' : {
        'datasetpath' : "/MET/Run2011B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "179434-180252",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_Mu15_L1ETM20_v3'    : '178420:MIN-179889:MAX',
            'HLT_IsoMu15_L1ETM20_v3' : '178420:MIN-179889:MAX',
            'HLT_Mu15_L1ETM20_v4'    : '179959:MIN-180252:MAX',
            'HLT_IsoMu15_L1ETM20_v4' : '179959:MIN-180252:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
    'data_MET_Run2011B_PromptReco_v1c' : {
        'datasetpath' : "/MET/Run2011B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "175832-180252",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_42_V20::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_Mu15_L1ETM20_v3'    : '178420:MIN-179889:MAX',
            'HLT_IsoMu15_L1ETM20_v3' : '178420:MIN-179889:MAX',
            'HLT_Mu15_L1ETM20_v4'    : '179959:MIN-180252:MAX',
            'HLT_IsoMu15_L1ETM20_v4' : '179959:MIN-180252:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    },
}
RECO_SAMPLES.update(TauIdEfficiencySpecific_RECO_SAMPLES)

# Define samples that get merged together
#
# NOTE: the merge sample name will become part of the histogram name
#       when running FWLiteTauIdEffAnalyzer
#
MERGE_SAMPLES = {
    'Data_2011RunA' : {
        'samples' : [
            'data_SingleMu_Run2011A_May10ReReco_v1',
            'data_SingleMu_Run2011A_PromptReco_v4',
            'data_SingleMu_Run2011A_Aug05ReReco_v1',
            'data_SingleMu_Run2011A_PromptReco_v6'
        ],
        'type' : 'Data'
    },
    'Data_2011RunB' : {
        'samples' : [
            #'data_MET_Run2011B_PromptReco_v1',
            #'data_MET_Run2011B_PromptReco_v1a'
            'data_MET_Run2011B_PromptReco_v1c'
        ],
        'type' : 'Data'
    },
    'Ztautau' : {
        'samples' : [
            ##'Ztautau_pythia'
            'Ztautau_powheg'
        ],
        'type' : plotter.process_Ztautau.config_dqmHistPlotter.type.value()
    },
    'Ztautau_embedded' : {
        'samples' : [
            'Ztautau_embedded_part1',
            'Ztautau_embedded_part2'
        ],
        'type' : plotter.process_Ztautau.config_dqmHistPlotter.type.value()
    },
    'Zmumu' : {
        'samples' : [
            ##'Zmumu_pythia'
            'Zmumu_powheg'
        ],
        'type' : plotter.process_Zmumu.config_dqmHistPlotter.type.value()
    },
    'QCD' : {
        'samples' : [
            'PPmuXptGt20Mu15'
        ],
        'type' : plotter.process_PPmuXptGt20.config_dqmHistPlotter.type.value()
    },
    'WplusJets' : {
        'samples' : [
            'WplusJets_madgraph'
        ],
        'type' : plotter.process_WplusJets.config_dqmHistPlotter.type.value()
    },
    'TTplusJets' : {
        'samples' : [
            'TTplusJets_madgraph'
        ],
        'type' : plotter.process_TTplusJets.config_dqmHistPlotter.type.value()
    }
}

ALL_SAMPLES = {}
# Update to use the defaults if necessary
for sample in RECO_SAMPLES.keys():
    defaults = copy.copy(ZtoMuTau.SAMPLE_DEFAULTS)
    defaults.update(RECO_SAMPLES[sample])
    RECO_SAMPLES[sample] = defaults
    # Combine MERGE and RECO samples in ALL samples
    # for simple access
    ALL_SAMPLES.update(MERGE_SAMPLES)
    ALL_SAMPLES.update(RECO_SAMPLES)

recoSampleDefinitionsTauIdEfficiency_7TeV = {
    'SAMPLES_TO_ANALYZE' : SAMPLES_TO_ANALYZE,
    'SAMPLE_DEFAULTS'    : ZtoMuTau.SAMPLE_DEFAULTS,
    'TARGET_LUMI'        : TARGET_LUMI,
    'RECO_SAMPLES'       : RECO_SAMPLES,
    'MERGE_SAMPLES'      : MERGE_SAMPLES,
    'ALL_SAMPLES'        : ALL_SAMPLES
}
