import FWCore.ParameterSet.Config as cms
import copy
import itertools

import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as styles

import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as ZtoMuTau

# List of samples to run in the analysis
SAMPLES_TO_ANALYZE = [
    'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to193621',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs193752to194076v2',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194108to194479',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194790to195016',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195099to195947',
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195948to196509',
    'Ztautau_pythia',
    'Zmumu_pythia',
    'DYmumuM10to20_pythia',
    'ZplusJets_madgraph2',
    'WplusJets_madgraph',
    'WplusJets_madgraph_extension',
    'PPmuXptGt20Mu15',
    'PPmuXptGt20Mu15v2',
    'WW',
    'WZ',
    'ZZ',
    'TTplusJets_madgraph2'
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
       29.9 # HLT_IsoMu15_eta2p1_L1ETM20_v3 (190645-190738)
    + 668.5 # HLT_IsoMu15_eta2p1_L1ETM20_v4 (191057-193621)
    + 890.6 # HLT_IsoMu15_eta2p1_L1ETM20_v5 (193998-194479)
)/_picobarns

RECO_SAMPLES = copy.deepcopy(ZtoMuTau.RECO_SAMPLES)
TauIdEfficiencySpecific_RECO_SAMPLES = {
    'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to193621' : {
        'datasetpath' : "/TauPlusX/Run2012A-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "190456-193621",
        'number_of_jobs' : 500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v3' : '190645:MIN-190738:MAX',
            'HLT_IsoMu15_eta2p1_L1ETM20_v4' : '191057:MIN-193621:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    },
    'data_TauPlusX_Run2012B_PromptReco_v1_runs193752to194076v2' : {
        'datasetpath' : "/TauPlusX/Run2012B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "193752-194076",
        'number_of_jobs' : 1500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v5' : '193752:MIN-194076:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    },
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194108to194479' : {
        'datasetpath' : "/TauPlusX/Run2012B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "194108-194479",
        'number_of_jobs' : 1500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v5' : '194108:MIN-194479:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    },
    'data_TauPlusX_Run2012B_PromptReco_v1_runs194790to195016' : {
        'datasetpath' : "/TauPlusX/Run2012B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-195016_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "194790-195016",
        'number_of_jobs' : 1500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v5' : '194790:MIN-195016:MAX'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    },
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195099to195947' : {
        'datasetpath' : "/TauPlusX/Run2012B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-195947_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "195099-195947",
        'number_of_jobs' : 1500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v5' : '195099:MIN-:MAX195947'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    },
    'data_TauPlusX_Run2012B_PromptReco_v1_runs195948to196509' : {
        'datasetpath' : "/TauPlusX/Run2012B-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-196509_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "195948-196509",
        'number_of_jobs' : 1500,
        'conditions' : 'GR_R_52_V7C::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v5' : '195948:MIN-:MAX196509'
        },
        'enableSysUncertainties' : False,
        'enableFakeRates' : True,
        'hlt' : cms.InputTag("TriggerResults", "", "HLT")
    }
}
RECO_SAMPLES.update(TauIdEfficiencySpecific_RECO_SAMPLES)

# Define samples that get merged together
#
# NOTE: the merge sample name will become part of the histogram name
#       when running FWLiteTauIdEffAnalyzer
#
# Define samples that get merged together
#
# NOTE: the merge sample name will become part of the histogram name
#       when running FWLiteTauIdEffAnalyzer
#
MERGE_SAMPLES = {
    'Data_2012RunA' : {
        'samples' : [
            'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to193621',
            'data_TauPlusX_Run2012B_PromptReco_v1_runs193752to194076v2',
            'data_TauPlusX_Run2012B_PromptReco_v1_runs194108to194479',
            'data_TauPlusX_Run2012B_PromptReco_v1_runs194790to195016',
            'data_TauPlusX_Run2012B_PromptReco_v1_runs195099to195947',
            'data_TauPlusX_Run2012B_PromptReco_v1_runs195948to196509'
        ],
        'type' : 'Data'
    },    
    'Ztautau' : {
        'samples' : [
            'Ztautau_pythia'
        ],
        'type' : plotter.process_Ztautau.config_dqmHistPlotter.type.value()
    },
    'Zmumu' : {
        'samples' : [
            'Zmumu_pythia',
            'DYmumuM10to20_pythia'
        ],
        'type' : plotter.process_Zmumu.config_dqmHistPlotter.type.value()
    },
    'ZplusJets' : {
        'samples' : [
            'ZplusJets_madgraph2'
        ],
        'type' : plotter.process_Zmumu.config_dqmHistPlotter.type.value()
    },
    'QCD' : {
        'samples' : [
            'PPmuXptGt20Mu15v2'
        ],
        'type' : plotter.process_PPmuXptGt20.config_dqmHistPlotter.type.value()
    },
    'WplusJets' : {
        'samples' : [
            'WplusJets_madgraph',
            'WplusJets_madgraph_extension'
        ],
        'type' : plotter.process_WplusJets.config_dqmHistPlotter.type.value()
    },
    'TTplusJets' : {
        'samples' : [
            'TTplusJets_madgraph2'
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
