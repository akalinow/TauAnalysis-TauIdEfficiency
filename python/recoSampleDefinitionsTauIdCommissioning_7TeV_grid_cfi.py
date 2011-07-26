import FWCore.ParameterSet.Config as cms

import TauAnalysis.Configuration.plotterProcessDefinitions_cfi as plotter
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as styles

# List of samples to run in the analysis
SAMPLES_TO_RUN = [
    'data_Jet_Run2011A_May10ReReco_v1',
    'data_Jet_Run2011A_PromptReco_v4',    
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'ZplusJets',
    'WplusJets',
    'qcdDiJetPtHat15to30',
    'qcdDiJetPtHat30to50',
    'qcdDiJetPtHat50to80',
    'qcdDiJetPtHat80to120',
    'qcdDiJetPtHat120to170',
    'qcdDiJetPtHat170to300',
    'PPmuXptGt20Mu15',
    'TTplusJets'
]

JOBS_TO_RUN = [
    'qcdDiJet',
    'qcdMuEnriched',
    'WplusJets',
    'Zmumu'
]

CONFIG_FILES = {
    'qcdDiJet'      : "produceCommissioningQCDdiJetPATTuple_cfg.py",
    'qcdMuEnriched' : "produceCommissioningQCDmuEnrichedPATTuple_cfg.py",
    'WplusJets'     : "produceCommissioningWplusJetsEnrichedPATTuple_cfg.py",
    'Zmumu'         : "produceCommissioningZmumuEnrichedPATTuple_cfg.py"
}

ROOT_FILE_NAMES = {
    'qcdDiJet'      : "tauCommissioningQCDdiJetPATtuple.root",
    'qcdMuEnriched' : "tauCommissioningQCDmuEnrichedPATtuple.root",
    'WplusJets'     : "tauCommissioningWplusJetsEnrichedPATtuple.root",
    'Zmumu'         : "tauCommissioningZmumuEnrichedPATtuple.root"
}

JOB_OPTIONS = {
    'qcdDiJet' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'qcdMuEnriched' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'WplusJets' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'Zmumu' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    }
}

#--------------------------------------------------------------------------------
# NOTE: Tau fake-rate PAT-tuple production requires RECO event content,
#       in order to rerun the Calo/TCTau reconstruction (based on TrackExtra objects)
#--------------------------------------------------------------------------------

_millibarns = 1.0e+9
_picobarns =  1.0
_femtobarns = 1.0e-3

TARGET_LUMI = 1.1/_femtobarns # for runs 160404 - 167913 ("golden" quality) and unprescaled triggers

RECO_SAMPLES = {
    # JetMET datasets
    'data_Jet_Run2011A_May10ReReco_v1' : {
        'datasetpath'      : '/Jet/Run2011A-May10ReReco-v1/AOD',
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'        : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'     : "160329-161312",
        'lumis_per_job'    : "25",
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'Data',
        'hlt'              : 'HLT'
    },
    'data_Jet_Run2011A_PromptReco_v4' : {
        'datasetpath'      : '/Jet/Run2011A-PromptReco-v4/AOD',
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'        : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'     : "165071-167913",
        'lumis_per_job'    : "25",
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'Data',
        'hlt'              : 'HLT'
    },
    # Muon datasets
    'data_SingleMu_Run2011A_May10ReReco_v1' : {
        'datasetpath'      : '/SingleMu/Run2011A-May10ReReco-v1/AOD',
        'dbs_url'          :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'        : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.txt",
        'runselection'     : "160329-163869",
        'lumis_per_job'    : "25",
        'jobs'             : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'Data',
        'hlt'              : 'HLT'
    },
    'data_SingleMu_Run2011A_PromptReco_v4' : {
        'datasetpath'      : "/SingleMu/Run2011A-PromptReco-v4/AOD",
        'dbs_url'          :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'        : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'     : "165071-167913",
        'lumis_per_job'    : "25",
        'jobs'             : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'Data',
        'hlt'              : 'HLT'
    },
    # Monte Carlo samples
    'ZplusJets' : {
        'datasetpath'      : "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 36277961,
        'x_sec'            : 2946*_picobarns, # NLO cross-section for Z --> l+ l-, M(l+ l-) > 50 GeV
                                              # taken from http://alcaraz.web.cern.ch/alcaraz/CROSS_SECTIONS.txt
        'jobs'             : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'WplusJets' : {
        'datasetpath'      : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'          :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 49527177,
        'x_sec'            : 31314*_picobarns,
        'jobs'             : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat15to30' : {
        'datasetpath'      : "/QCD_Pt-15to30_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 11000000,
        'x_sec'            : 8.16e8*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat30to50' : {
        'datasetpath'      : "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6583068,
        'x_sec'            : 5.31e7*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat50to80' : {
        'datasetpath'      : "/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6600000,
        'x_sec'            : 6.36e6*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat80to120' : {
        'datasetpath'      : "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6589956,
        'x_sec'            : 7.84e5*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat120to170' : {
        'datasetpath'      : "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6127528,
        'x_sec'            : 1.15e5*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat170to300' : {
        'datasetpath'      : "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6220160,
        'x_sec'            : 2.43e4*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'qcdDiJetPtHat300to470' : {
        'datasetpath'      : "/QCD_Pt-300to470_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 6432669,
        'x_sec'            : 1.17e3*_picobarns,
        'jobs'             : [ 'qcdDiJet' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'PPmuXptGt20Mu15' : {
        'datasetpath'      : "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 20416038,
        'x_sec'            : 0.2966*_millibarns*2.855e-4, # x-sec * gen filter efficiency
        'jobs'             : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    },
    'TTplusJets' : {
        'datasetpath'      : "/TTJets_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'          : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'events_processed' : 3701947,
        'x_sec'            : 157.5*_picobarns, # NLO cross-section from
                                               # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
        'jobs'             : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'             : 'MC',
        'hlt'              : 'HLT'
    }
}

# Define samples that get merged together
#
# NOTE:
#      (1) the merge sample name will become part of the histogram name
#          when running FWLiteTauFakeRateAnalyzer
#      (2) add Data samples triggered by Jet triggers and single Muon triggers
#         --> either one or the other will exist...
#      (3) add QCD multi-jet and PPmuX samples
#         --> either one or the other will exist...
#
MERGE_SAMPLES = {
    'Data' : {
        'samples' : [
            'data_Jet_Run2011A_May10ReReco_v1',
            'data_Jet_Run2011A_PromptReco_v4',    
            'data_SingleMu_Run2011A_May10ReReco_v1',
            'data_SingleMu_Run2011A_PromptReco_v4'
        ],
        'type' : 'Data',
        'legendEntry' : "Data",
        'drawOption' : styles.drawOption_Data
    },
    'ZplusJets' : {
        'samples' : [
            'ZplusJets'
        ],
        'type' : plotter.process_Ztautau.config_dqmHistPlotter.type.value(),
        'legendEntry' : plotter.process_Ztautau.config_dqmHistPlotter.legendEntry.value(),
        'drawOption' : styles.drawOption_ZplusJets
    },
    'WplusJets' : {
        'samples' : [
            'WplusJets_madgraph'
        ],
        'type' : plotter.process_WplusJets.config_dqmHistPlotter.type.value(),
        'legendEntry' : plotter.process_WplusJets.config_dqmHistPlotter.legendEntry.value(),
        'drawOption' : styles.drawOption_WplusJets
    },
    'QCD' : {
        'samples' : [
            'qcdDiJetPtHat15to30',
            'qcdDiJetPtHat30to50',
            'qcdDiJetPtHat50to80',
            'qcdDiJetPtHat80to120',
            'qcdDiJetPtHat120to170',
            'qcdDiJetPtHat170to300',
            'PPmuXptGt20Mu15'
        ],
        'type' : plotter.process_PPmuXptGt20.config_dqmHistPlotter.type.value(),
        'legendEntry' : plotter.process_PPmuXptGt20.config_dqmHistPlotter.legendEntry.value(),
        'drawOption' : styles.drawOption_QCD
    },
    'TTplusJets' : {
        'samples' : [
            'TTplusJets_madgraph'
        ],
        'type' : plotter.process_TTplusJets.config_dqmHistPlotter.type.value(),
        'legendEntry' : plotter.process_TTplusJets.config_dqmHistPlotter.legendEntry.value(),
        'drawOption' : styles.drawOption_TTplusJets
    }
}

recoSampleDefinitionsTauIdCommissioning_7TeV = {
    'SAMPLES_TO_RUN'  : SAMPLES_TO_RUN,    
    'JOBS_TO_RUN'     : JOBS_TO_RUN,
    'CONFIG_FILES'    : CONFIG_FILES, 
    'ROOT_FILE_NAMES' : ROOT_FILE_NAMES,
    'JOB_OPTIONS'     : JOB_OPTIONS,
    'TARGET_LUMI'     : TARGET_LUMI,
    'RECO_SAMPLES'    : RECO_SAMPLES,
    'MERGE_SAMPLES'   : MERGE_SAMPLES
}
