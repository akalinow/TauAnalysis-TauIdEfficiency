import FWCore.ParameterSet.Config as cms
import copy
import itertools

import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as styles

import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as ZtoMuTau

# List of samples to run in the analysis
SAMPLES_TO_ANALYZE = [
    'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to191859',
    'ZplusJets_madgraph',
    #'WplusJets_madgraph',
    'Wenu_pythia',
    'Wmunu_pythia',
    'Wtaunu_pythia',  
    'PPmuXptGt20Mu15',      
    'WW',
    'WZ',
    'ZZ',
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
       21.30 # HLT_IsoMu15_eta2p1_L1ETM20_v3 (190645-190738)
    + 346.31 # HLT_IsoMu15_eta2p1_L1ETM20_v4 (191057-191810)
)/_picobarns

RECO_SAMPLES = copy.deepcopy(ZtoMuTau.RECO_SAMPLES)
TauIdEfficiencySpecific_RECO_SAMPLES = {
    'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to191859' : {
        'datasetpath' : "/TauPlusX/Run2012A-PromptReco-v1/AOD",
        'dbs_url' :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask' : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-191859_8TeV_PromptReco_Collisions12_JSON.txt",
        'runselection' : "190456-191859",
        'number_of_jobs' : 2500,
        'conditions' : 'GR_R_52_V7::All',
        'events_processed' : -1,
        'skim_eff' : 1.0,
        'type' : 'Data',
        'drawOption' : styles.drawOption_Data,
        'hlt_paths' : {
            'HLT_IsoMu15_eta2p1_L1ETM20_v3' : '190645:MIN-190738:MAX',
            'HLT_IsoMu15_eta2p1_L1ETM20_v4' : '191057:MIN-191810:MAX'
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
MERGE_SAMPLES = copy.deepcopy(ZtoMuTau.MERGE_SAMPLES)
TauIdEfficiencySpecific_MERGE_SAMPLES = {
    'Data_2012RunA' : {
        'samples' : [
            'data_TauPlusX_Run2012A_PromptReco_v1_runs190456to191859'
        ],
        'type' : 'Data'
    }
}
MERGE_SAMPLES.update(TauIdEfficiencySpecific_RECO_SAMPLES)

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
