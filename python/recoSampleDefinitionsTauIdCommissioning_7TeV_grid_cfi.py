# List of samples to run in the analysis
SAMPLES_TO_RUN = [
    'data_Jet_Run2011A_PromptReco_v1',
    'data_Mu_Run2011A_PromptReco_v1',
    'data_Mu_Run2011A_PromptReco_v2',
    'Ztautau',
    'ZplusJets',
    'qcdDiJetPtHat15to30',
    'qcdDiJetPtHat30to50',
    'qcdDiJetPtHat50to80',
    'qcdDiJetPtHat80to120',
    'qcdDiJetPtHat120to170',
    'qcdDiJetPtHat170to300',
    'PPmuXptGt20Mu15',
    'WplusJets'
]

JOBS_TO_RUN = [
    #'qcdDiJet',
    #'qcdMuEnriched',
    'WplusJets',
    'Zmumu'
]

CONFIG_FILES = {
    'qcdDiJet'               : "produceCommissioningQCDdiJetNtuple_cfg.py",
    'qcdMuEnriched'          : "produceCommissioningQCDmuEnrichedNtuple_cfg.py",
    'WplusJets'              : "produceCommissioningWplusJetsEnrichedNtuple_cfg.py",
    'Zmumu'                  : "produceCommissioningZmumuEnrichedNtuple_cfg.py"
}

ROOT_FILE_NAMES = {
    'qcdDiJet'               : "tauIdEffEDNtuple_qcdDiJet.root",
    'qcdMuEnriched'          : "tauIdEffEDNtuple_qcdMuEnriched.root",
    'WplusJets'              : "tauIdEffEDNtuple_wPlusJetsEnriched.root",
    'Zmumu'                  : "tauIdEffEDNtuple_zMuMuEnriched.root"
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
        'submitTypes'         : [ 'MC' ]
    }
}

#--------------------------------------------------------------------------------
# NOTE: Tau(ED)Ntuple production requires RECO event content,
#       in order to rerun the Calo/TCTau reconstruction (based on TrackExtra objects)
#--------------------------------------------------------------------------------

RECO_SAMPLES = {
    # JetMET secondary datasets
    'data_Jet_Run2011A_PromptReco_v1' : {
        'datasetpath'   : '/Jet/Run2011A-PromptReco-v1/AOD',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163369_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "160404-161216",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    # Muon secondary datasets
    'data_Mu_Run2011A_PromptReco_v1' : {
        'datasetpath'   : '/SingleMu/Run2011A-PromptReco-v1/AOD',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163369_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "160404-161216",
        'lumis_per_job' : "50",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_Mu_Run2011A_PromptReco_v2' : {
        'datasetpath'   : '/SingleMu/Run2011A-PromptReco-v2/AOD',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-163369_7TeV_PromptReco_Collisions11_JSON.txt",
        'runselection'  : "162718-163796",
        'lumis_per_job' : "50",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets', 'Zmumu' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    # Monte Carlo samples
    'Ztautau' : {
        'datasetpath'   : "/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Spring11-PU_S1_START311_V1G1-v2/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'ZplusJets' : {
        'datasetpath'   : "/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'Zmumu' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat15to30' : {
        'datasetpath'   : "/QCD_Pt_15to30_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat30to50' : {
        'datasetpath'   : "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat50to80' : {
        'datasetpath'   : "/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat80to120' : {
        'datasetpath'   : "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat120to170' : {
        'datasetpath'   : "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'qcdDiJetPtHat170to300' : {
        'datasetpath'   : "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'PPmuXptGt20Mu15' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    },
    'WplusJets' : {
        'datasetpath'   : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Spring11-PU_S1_START311_V1G1-v1/AODSIM",
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'WplusJets' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI311X'
    }
}

# add jobs with no event selection applied (for MC only)
for sample in RECO_SAMPLES.keys():
    if RECO_SAMPLES[sample]['type'] == 'MC':
        jobs          = RECO_SAMPLES[sample]['jobs']
        jobs_mcNoCuts = []
        for job in jobs:
            jobs_mcNoCuts.append(job)
            jobs_mcNoCuts.append(job + "_mcNoCuts")
            RECO_SAMPLES[sample]['jobs'] = jobs_mcNoCuts
