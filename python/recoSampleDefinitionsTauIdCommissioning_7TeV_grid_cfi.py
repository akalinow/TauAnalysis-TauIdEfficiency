# List of samples to run in the analysis
SAMPLES_TO_RUN = [
    'data_MinBias_Commissioning_Jun14ReReco',
    'data_JetMETTau_Run2010Ai_Sep17ReReco',
    'data_JetMET_Run2010Aii_Sep17ReReco',
    'data_Jet_Run2010B_Nov4ReReco',
    'data_Mu_Run2010A_Nov4ReReco',
    'data_Mu_Run2010B_Nov4ReReco',
    'Ztautau',
    'ZtautauPU156bx',
    'Ztautau_powhegZ2',
    'qcdDiJetPtHat5to15',
    'qcdDiJetPtHat5to15PU156bx',
    'qcdDiJetPtHat15to30',
    'qcdDiJetPtHat15to30PU156bx',
    'qcdDiJetPtHat30to50',
    'qcdDiJetPtHat30to50PU156bx',
    'qcdDiJetPtHat50to80',
    'qcdDiJetPtHat50to80PU156bx',
    'qcdDiJetPtHat80to120',
    'qcdDiJetPtHat80to120PU156bx',
    'qcdDiJetPtHat120to170',
    'qcdDiJetPtHat120to170PU156bx',
    'qcdDiJetPtHat170to300',
    'qcdDiJetPtHat170to300PU156bx',
    'PPmuXptGt20Mu10',
    'PPmuXptGt20Mu10PU156bx',
    'PPmuXptGt20Mu15',
    'PPmuXptGt20Mu15PU156bx',
    'WplusJets',
    'WplusJetsPU156bx'
]

JOBS_TO_RUN = [
    'qcdDiJet',
    #'qcdDiJet_mcNoCuts',
    'qcdMuEnriched',
    #'qcdMuEnriched_mcNoCuts',
    'WplusJets',
    #'WplusJets_mcNoCuts'
]

CONFIG_FILES = {
    'qcdDiJet'               : "produceCommissioningQCDdiJetNtuple_cfg.py",
    'qcdDiJet_mcNoCuts'      : "produceCommissioningQCDdiJetNtuple_cfg.py",
    'qcdMuEnriched'          : "produceCommissioningQCDmuEnrichedNtuple_cfg.py",
    'qcdMuEnriched_mcNoCuts' : "produceCommissioningQCDmuEnrichedNtuple_cfg.py",
    'WplusJets'              : "produceCommissioningWplusJetsEnrichedNtuple_cfg.py",
    'WplusJets_mcNoCuts'     : "produceCommissioningWplusJetsEnrichedNtuple_cfg.py"
}

ROOT_FILE_NAMES = {
    'qcdDiJet'               : "tauIdEffEDNtuple_qcdDiJet.root",
    'qcdDiJet_mcNoCuts'      : "tauIdEffEDNtuple_qcdDiJet.root",
    'qcdMuEnriched'          : "tauIdEffEDNtuple_qcdMuEnriched.root",
    'qcdMuEnriched_mcNoCuts' : "tauIdEffEDNtuple_qcdMuEnriched.root",
    'WplusJets'              : "tauIdEffEDNtuple_wPlusJetsEnriched.root",
    'WplusJets_mcNoCuts'     : "tauIdEffEDNtuple_wPlusJetsEnriched.root"
}

JOB_OPTIONS = {
    'qcdDiJet' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'qcdDiJet_mcNoCuts' :{
        'applyEventSelection' : False,
        'submitTypes'         : [ 'MC' ]
    },
    'qcdMuEnriched' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'qcdMuEnriched_mcNoCuts' : {
        'applyEventSelection' : False,
        'submitTypes'         : [ 'MC' ]
    },
    'WplusJets' : {
        'applyEventSelection' : True,
        'submitTypes'         : [ 'MC', 'Data' ]
    },
    'WplusJets_mcNoCuts' : {
        'applyEventSelection' : False,
        'submitTypes'         : [ 'MC' ]
    }
}

RECO_SAMPLES = {
    # JetMET secondary datasets
    'data_MinBias_Commissioning_Jun14ReReco' : {
        'datasetpath'   : '/MinimumBias/Commissioning10-Jun14thReReco_v1/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.txt",
        'runselection'  : "132440-135802",
        'lumis_per_job' : "50",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_JetMETTau_Run2010Ai_Sep17ReReco' : {
        'datasetpath'   : '/JetMETTau/Run2010A-Sep17ReReco_v2/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.txt",
        'runselection'  : "135821-141887",
        'lumis_per_job' : "5",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_JetMET_Run2010Aii_Sep17ReReco' : {
        'datasetpath'   : '/JetMET/Run2010A-Sep17ReReco_v2/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON.txt",
        'runselection'  : "141950-144114",
        'lumis_per_job' : "5",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    'data_Jet_Run2010B_Nov4ReReco' : {
        'datasetpath'   : '/Jet/Run2010B-Nov4ReReco_v1/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON.txt",
        'runselection'  : "146428-149442",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'Data',
        'hlt'           : 'HLT',
        'SE_black_list' : 'T2_BE_IIHE'
    },
    # Muon secondary datasets
    'data_Mu_Run2010A_Nov4ReReco' : {
        'datasetpath'   : '/Mu/Run2010A-Nov4ReReco_v1/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON.txt",
        'runselection'  : "136033-144114",
        'lumis_per_job' : "50",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets' ],
        'type'          : 'Data',
        'hlt'           : 'HLT',
        'SE_black_list' : 'T2_IT_Legnaro'
    },
    'data_Mu_Run2010B_Nov4ReReco' : {
        'datasetpath'   : '/Mu/Run2010B-Nov4ReReco_v1/RECO',
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'lumi_mask'     : "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions10/7TeV/Reprocessing/Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_JSON.txt",
        'runselection'  : "146428-149442",
        'lumis_per_job' : "25",
        'jobs'          : [ 'qcdMuEnriched', 'WplusJets' ],
        'type'          : 'Data',
        'hlt'           : 'HLT'
    },
    # Monte Carlo samples
    'Ztautau' : {
        'datasetpath'   : "/DYtoTauTau_M_20_TuneD6T_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'ZtautauPU156bx' : {
        'datasetpath'   : "/DYtoTauTau_M_20_TuneD6T_7TeV-pythia6-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'Ztautau_powhegZ2' : {
        'datasetpath'   : "/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v3/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ], 
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat5to15' : {
        'datasetpath'   : "/QCD_Pt_5to15_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat5to15PU156bx' : {
        'datasetpath'   : "/QCD_Pt_5to15_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat15to30' : {
        'datasetpath'   : "/QCD_Pt_15to30_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat15to30PU156bx' : {
        'datasetpath'   : "/QCD_Pt_15to30_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat30to50' : {
        'datasetpath'   : "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat30to50PU156bx' : {
        'datasetpath'   : "/QCD_Pt_30to50_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat50to80' : {
        'datasetpath'   : "/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat50to80PU156bx' : {
        'datasetpath'   : "/QCD_Pt_50to80_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat80to120' : {
        'datasetpath'   : "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat80to120PU156bx' : {
        'datasetpath'   : "/QCD_Pt_80to120_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat120to170' : {
        'datasetpath'   : "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat120to170PU156bx' : {
        'datasetpath'   : "/QCD_Pt_120to170_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'qcdDiJetPtHat170to300' : {
        'datasetpath'   : "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'qcdDiJetPtHat170to300PU156bx' : {
        'datasetpath'   : "/QCD_Pt_170to300_TuneZ2_7TeV_pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdDiJet' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'PPmuXptGt20Mu10' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-10_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X'
    },
    'PPmuXptGt20Mu10PU156bx' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-10_TuneZ2_7TeV-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'PPmuXptGt20Mu15' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38X',
        'SE_black_list' : 'T2_US_Florida,T2_RU_ITEP'
    },
    'PPmuXptGt20Mu15PU156bx' : {
        'datasetpath'   : "/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       : "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'qcdMuEnriched' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
    },
    'WplusJets' : {
        'datasetpath'   : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'WplusJets' ],
        'type'          : 'MC',
        'hlt'           : 'HLT'
    },
    'WplusJetsPU156bx' : {
        'datasetpath'   : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/GEN-SIM-RECO",
        'dbs_url'       :  "http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet",
        'jobs'          : [ 'WplusJets' ],
        'type'          : 'MC',
        'hlt'           : 'REDIGI38XPU'
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
