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
# NOTE: Tau(ED)Ntuple production requires RECO event content,
#       in order to rerun the Calo/TCTau reconstruction (based on TrackExtra objects)
#--------------------------------------------------------------------------------

RECO_SAMPLES = {
    # JetMET datasets
    'data_Jet_Run2011A_May10ReReco_v1' : {
        'datasetpath'   : '/Jet/Run2011A-May10ReReco-v1/AOD',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "160329-161312",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_Jet_Run2011A_PromptReco_v4' : {
        'datasetpath'   : '/Jet/Run2011A-PromptReco-v4/AOD',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "165071-167913",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    # Muon datasets
    'data_SingleMu_Run2011A_May10ReReco_v1' : {
        'datasetpath'   : '/SingleMu/Run2011A-May10ReReco-v1/AOD',
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.txt",
        'runselection'  : "160329-163869",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_SingleMu_Run2011A_PromptReco_v4' : {
        'datasetpath'   : "/SingleMu/Run2011A-PromptReco-v4/AOD",
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "165071-167913",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    # Monte Carlo samples
    'ZplusJets' : {
        'datasetpath'   : "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'WplusJets' : {
        'datasetpath'   : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat15to30' : {
        'datasetpath'   : "/QCD_Pt-15to30_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat30to50' : {
        'datasetpath'   : "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat50to80' : {
        'datasetpath'   : "/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat80to120' : {
        'datasetpath'   : "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat120to170' : {
        'datasetpath'   : "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'qcdDiJetPtHat170to300' : {
        'datasetpath'   : "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'PPmuXptGt20Mu15' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'TTplusJets' : {
        'datasetpath'   : "/TTJets_TuneZ2_7TeV-madgraph-tauola/Summer11-PU_S4_START42_V11-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet', 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    }
}
