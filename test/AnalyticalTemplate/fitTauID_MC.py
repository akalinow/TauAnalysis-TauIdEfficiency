import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#filePath = "/dmj/hep/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/18_01_2017/"
#filePath = "/cms/cms/akalinow/CMS/TauID/Crab/Data/TauID_TnP_2017/v6_Mu2Tau_2017/"
#filePath = "/mnt/shared/cms/akalinow/CMS/TauID/Crab/Data/TauID_TnP_2017/Mu2Tau_2017_v9/"
filePath = "/mnt/shared/cms/akalinow/CMS/TauID/Crab/Data/TauID_TnP_2018/Mu2Tau_2018_v2/"

filePath += "tnpZ_MCwithWeights_filtered.root"

efficiencyPSetTemplate = cms.PSet(
    #UnbinnedVariables = cms.vstring("mass", "weight", "alternatLorentzVectPt", "alternatLorentzVectEta", "tag_pt", "tag_triggerMatch", "tag_dB", "pair_dz", "pair_deltaR", "pair_probeMultiplicity", "pair_BestZ", "pair_MET", "pair_MTtag", "pair_MTprobe", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    UnbinnedVariables = cms.vstring("mass", "weight"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        abseta = cms.vdouble(0.0, 0.4, 0.8, 1.2, 1.7, 2.3),
        ## flags and conditions required at the denominator,
    ),
    BinToPDFmap = cms.vstring("Zll_Model"), ## PDF to use, as defined below
)

againstMuonLoose3_Zmumu = efficiencyPSetTemplate.clone()
againstMuonLoose3_Zmumu.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0.5,1.0))
##againstMuonLoose3_Zmumu.BinnedVariables._Parameterizable__addParameter("pair_ZDecayMode",cms.vdouble(6.5,7.5))

againstMuonLoose3_Zmumu.BinToPDFmap = cms.vstring("Zll_Model_LooseEta0","*abseta_bin1*","Zll_Model_LooseEta1", "*abseta_bin2*","Zll_Model_LooseEta2", "*abseta_bin3*","Zll_Model_LooseEta3", "*abseta_bin4*","Zll_Model_LooseEta4")

againstMuonLoose3_Zll = efficiencyPSetTemplate.clone()
againstMuonLoose3_Zll.BinToPDFmap = cms.vstring("Zll_Model_LooseEta0","*abseta_bin1*","Zll_Model_LooseEta1", "*abseta_bin2*","Zll_Model_LooseEta2",  "*abseta_bin3*","Zll_Model_LooseEta3",  "*abseta_bin4*","Zll_Model_LooseEta4")

againstMuonTight3_Zmumu = againstMuonLoose3_Zmumu.clone()
againstMuonTight3_Zmumu.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")
againstMuonTight3_Zmumu.BinToPDFmap = cms.vstring("Zll_Model_TightEta0","*abseta_bin1*","Zll_Model_TightEta1", "*abseta_bin2*","Zll_Model_TightEta2", "*abseta_bin3*","Zll_Model_TightEta3", "*abseta_bin4*","Zll_Model_TightEta4")

againstMuonTight3_Zll = efficiencyPSetTemplate.clone()
againstMuonTight3_Zll.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")
againstMuonTight3_Zll.BinToPDFmap = cms.vstring("Zll_Model_TightEta0","*abseta_bin1*","Zll_Model_TightEta1", "*abseta_bin2*","Zll_Model_TightEta2", "*abseta_bin3*","Zll_Model_TightEta3", "*abseta_bin4*","Zll_Model_TightEta4")

Zll_Model_LooseEta0_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LooseEta1_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LooseEta2_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LooseEta3_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LooseEta4_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightEta0_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightEta1_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightEta2_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightEta3_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightEta4_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

###############################################################
###############################################################
process.TnP_Muon_ID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output
    InputFileNames = cms.vstring("file:"+filePath),
    OutputFileName = cms.string("TnP_MuonToTau_MisID_MC.root"),
    InputTreeName = cms.string("fitter_tree"), 
    InputDirectoryName = cms.string("tpTree"),  
    WeightVariable = cms.string("weight"),
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-muon Mass", "60", "120", "GeV/c^{2}"),
        abseta = cms.vstring("muon |#eta|", "0", "2.3", ""),
        alternatLorentzVectPt = cms.vstring("probe tau pT", "20", "1500", ""),
        alternatLorentzVectEta = cms.vstring("probe tau eta", "-2.3", "2.3", ""),
        tag_pt  = cms.vstring("tag pT", "0", "1500", ""),
        tag_triggerMatch = cms.vstring("Tag matched to HLT item", "0.5", "1.0", ""),
        tag_dB  = cms.vstring("dB", "0.0", "0.004", ""),
        pair_dz = cms.vstring("#Deltaz between two muons", "-0.01", "0.01", "cm"),
        pair_deltaR = cms.vstring("#DeltaR between two muons", "0.5", "5", ""),
        pair_probeMultiplicity = cms.vstring("Probe multiplicity", "1", "1", ""),
        pair_BestZ = cms.vstring("pair closest to Z mass", "0.5", "1.0", ""),
        pair_MET = cms.vstring("MET", "0", "2500", ""),
        pair_MTtag = cms.vstring("MTtag", "0", "40", ""),
        pair_MTprobe = cms.vstring("MTprobe", "0", "4000", ""),
        decayModeFindingNewDMs = cms.vstring("Decay mode finding NewDMs", "0.5", "1.0", ""),
        byTightIsolationMVArun2v1DBnewDMwLT2017v2 = cms.vstring("Tight MVARun2_2017v2 id.", "0.5", "1.0", ""),
        mcTrue = cms.vstring("Match to gen muons", "0.0", "1.0", ""),
        pair_ZDecayMode = cms.vstring("Z decay mode", "0", "10", ""),
        pair_alternativeMass = cms.vstring("mass from tau and muon", "-70", "1200", ""),
	weight = cms.vstring("weight","0","10","")
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        againstMuonTight3 = cms.vstring("againstMuonTight3", "dummy[pass=1,fail=0]"),
        againstMuonLoose3 = cms.vstring("againstMuonLoose3", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        againstMuonLoose3_Zmumu = againstMuonLoose3_Zmumu,
        againstMuonLoose3_Zll = againstMuonLoose3_Zll,
        againstMuonTight3_Zmumu = againstMuonTight3_Zmumu,
        againstMuonTight3_Zll = againstMuonTight3_Zll,
    ),
    PDFs = cms.PSet(
        Zll_Model_LooseEta0 = Zll_Model_LooseEta0_Template,
        Zll_Model_LooseEta1 = Zll_Model_LooseEta1_Template,
        Zll_Model_LooseEta2 = Zll_Model_LooseEta2_Template,
        Zll_Model_LooseEta3 = Zll_Model_LooseEta3_Template,
        Zll_Model_LooseEta4 = Zll_Model_LooseEta4_Template,

        Zll_Model_TightEta0 = Zll_Model_TightEta0_Template,
        Zll_Model_TightEta1 = Zll_Model_TightEta1_Template,
        Zll_Model_TightEta2 = Zll_Model_TightEta2_Template,
        Zll_Model_TightEta3 = Zll_Model_TightEta3_Template,
        Zll_Model_TightEta4 = Zll_Model_TightEta4_Template,
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



