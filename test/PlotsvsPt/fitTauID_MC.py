import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

filePath = "/cms/cms/akalinow/CMS/TauID/Crab/Data/TauID_TnP_2017/v6_Mu2Tau_2017/"
filePath += "tnpZ_MC.root"

efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass","weight", "alternatLorentzVectPt", "alternatLorentzVectEta", "tag_pt", "tag_triggerMatch", "tag_dB", "pair_dz", "pair_deltaR", "pair_probeMultiplicity", "pair_BestZ", "pair_MET", "pair_MTtag", "pair_MTprobe", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        pt = cms.vdouble(10, 20, 30, 40, 50, 100),
        ## flags and conditions required at the denominator,
    ),
    BinToPDFmap = cms.vstring("Zll_Model"), ## PDF to use, as defined below
)

againstMuonLoose3_Zmumu = efficiencyPSetTemplate.clone()
againstMuonLoose3_Zmumu.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0.5,1.0))
againstMuonLoose3_Zmumu.BinToPDFmap = cms.vstring("Zll_Model_LoosePt0","*pt_bin1*","Zll_Model_LoosePt1", "*pt_bin2*","Zll_Model_LoosePt2", "*pt_bin3*","Zll_Model_LoosePt3", "*pt_bin4*","Zll_Model_LoosePt4")

againstMuonLoose3_Zll = efficiencyPSetTemplate.clone()
againstMuonLoose3_Zll.BinToPDFmap = cms.vstring("Zll_Model_LoosePt0","*pt_bin1*","Zll_Model_LoosePt1", "*pt_bin2*","Zll_Model_LoosePt2",  "*pt_bin3*","Zll_Model_LoosePt3",  "*pt_bin4*","Zll_Model_LoosePt4")

againstMuonTight3_Zmumu = againstMuonLoose3_Zmumu.clone()
againstMuonTight3_Zmumu.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")
againstMuonTight3_Zmumu.BinToPDFmap = cms.vstring("Zll_Model_TightPt0","*pt_bin1*","Zll_Model_TightPt1", "*pt_bin2*","Zll_Model_TightPt2", "*pt_bin3*","Zll_Model_TightPt3", "*pt_bin4*","Zll_Model_TightPt4")

againstMuonTight3_Zll = efficiencyPSetTemplate.clone()
againstMuonTight3_Zll.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")
againstMuonTight3_Zll.BinToPDFmap = cms.vstring("Zll_Model_TightPt0","*pt_bin1*","Zll_Model_TightPt1", "*pt_bin2*","Zll_Model_TightPt2", "*pt_bin3*","Zll_Model_TightPt3", "*pt_bin4*","Zll_Model_TightPt4")

Zll_Model_LoosePt0_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LoosePt1_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LoosePt2_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LoosePt3_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_LoosePt4_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightPt0_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightPt1_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightPt2_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightPt3_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceFixedParams:backgroundFail",
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
            )

Zll_Model_TightPt4_Template = cms.vstring(
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceFixedParams:signalPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceFixedParams:signalFail",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceFixedParams:backgroundPass",
            "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceFixedParams:backgroundFail",
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
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-muon Mass", "60", "120", "GeV/c^{2}"),
        abseta = cms.vstring("muon |#eta|", "0", "0.8", ""),
        pt  = cms.vstring("probe pT", "0", "100", ""),
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
        decayModeFinding = cms.vstring("Decay mode finding", "0.5", "1.0", ""),
        #decayModeFindingNewDMs = cms.vstring("Decay mode finding NewDMs", "0.5", "1.0", ""),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vstring("Combined loose isolation", "0.5", "1.0", ""),
        mcTrue = cms.vstring("Match to gen muons", "0.0", "1.0", ""),
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

        #againstMuonTight3_Zmumu = againstMuonTight3_Zmumu,
        #againstMuonTight3_Zll = againstMuonTight3_Zll,
    ),
    PDFs = cms.PSet(
        Zll_Model_LoosePt0 = Zll_Model_LoosePt0_Template,
        Zll_Model_LoosePt1 = Zll_Model_LoosePt1_Template,
        Zll_Model_LoosePt2 = Zll_Model_LoosePt2_Template,
        Zll_Model_LoosePt3 = Zll_Model_LoosePt3_Template,
        Zll_Model_LoosePt4 = Zll_Model_LoosePt4_Template,

        Zll_Model_TightPt0 = Zll_Model_TightPt0_Template,
        Zll_Model_TightPt1 = Zll_Model_TightPt1_Template,
        Zll_Model_TightPt2 = Zll_Model_TightPt2_Template,
        Zll_Model_TightPt3 = Zll_Model_TightPt3_Template,
        Zll_Model_TightPt4 = Zll_Model_TightPt4_Template,
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(20),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



