import FWCore.ParameterSet.Config as cms
import copy
import itertools

import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as styles

import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as ZtoMuTau

# List of samples to run in the analysis
SAMPLES_TO_ANALYZE = [
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_May10ReReco_v1',
    #'DYtautauM10to20_powheg',
    'Ztautau_pythia',
    #'Ztautau_powheg',
    #'qqZll',
    #'DYmumuM10to20_pythia',
    'Zmumu_pythia',
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
TARGET_LUMI = 0.90*869.1/_picobarns # for runs 160404 - 167496 ("golden" quality) corrected for ~90% PFTau trigger efficiency

RECO_SAMPLES = copy.deepcopy(ZtoMuTau.RECO_SAMPLES)
TauIdEfficiencySpecific_RECO_SAMPLES = {
    'data_SingleMu_Run2011A_May10ReReco_v1' : {
        'datasetpath' : '/SingleMu/Run2011A-May10ReReco-v1/AOD',
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON.txt",
        'runselection' : "160329-163869",
        'number_of_jobs' : 750,
        'conditions' : 'GR_R_42_V14::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
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
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167784_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection' : "165071-167784",
        'number_of_jobs' : 2000,
        'conditions' : 'GR_R_42_V14::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu17_v8'  : '165088:MIN-165633:MAX',
            'HLT_IsoMu17_v9'  : '165970:MIN-167043:MAX',
            'HLT_IsoMu17_v10' : '166346:MIN-166346:MAX',
            'HLT_IsoMu17_v11' : '167078:MIN-167784:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT"),
        'SE_black_list' : 'T2_US_UCSD'
    }
}
RECO_SAMPLES.update(TauIdEfficiencySpecific_RECO_SAMPLES)

# Define samples that get merged together
#
# NOTE: the merge sample name will become part of the histogram name
#       when running FWLiteTauIdEffAnalyzer
#
MERGE_SAMPLES = {
    'Data' : {
        'samples' : [
            'data_SingleMu_Run2011A_May10ReReco_v1',
            'data_SingleMu_Run2011A_PromptReco_v4'
        ],
        'type' : 'Data'
    },
    'Ztautau' : {
        'samples' : [
            'Ztautau_pythia',
        ],
        'type' : plotter.process_Ztautau.config_dqmHistPlotter.type.value()
    },
    'Zmumu' : {
        'samples' : [
            'Zmumu_pythia'
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
    'SAMPLE_DEFAULTS' : ZtoMuTau.SAMPLE_DEFAULTS,
    'TARGET_LUMI' : TARGET_LUMI,
    'RECO_SAMPLES' : RECO_SAMPLES,
    'MERGE_SAMPLES' : MERGE_SAMPLES,
    'ALL_SAMPLES' : ALL_SAMPLES
}
