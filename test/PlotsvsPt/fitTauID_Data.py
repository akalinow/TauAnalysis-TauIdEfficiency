import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/15_09_2016/"
filePath += "tnpZ_Data.root"

efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass","alternatLorentzVectPt", "alternatLorentzVectEta","abseta", "tag_pt", "tag_triggerMatch", "tag_dB", "pair_dz", "pair_deltaR", "pair_probeMultiplicity", "pair_MET", "pair_MTtag", "pair_MTprobe", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        pt = cms.vdouble(10, 20, 30, 40, 50, 100),
        #run = cms.vdouble(271036, 274422, 275125),
        ## flags and conditions required at the denominator,
    ),
    BinToPDFmap = cms.vstring("Zll_Model"), ## PDF to use, as defined below
)

againstMuonLoose3_Data = efficiencyPSetTemplate.clone()
againstMuonLoose3_Data.BinToPDFmap = cms.vstring("Data_Model_LoosePt0","*pt_bin1*","Data_Model_LoosePt1", "*pt_bin2*","Data_Model_LoosePt2", "*pt_bin3*","Data_Model_LoosePt3", "*pt_bin4*","Data_Model_LoosePt4")

againstMuonTight3_Data = efficiencyPSetTemplate.clone()
againstMuonTight3_Data.BinToPDFmap = cms.vstring("Data_Model_TightPt0","*pt_bin1*","Data_Model_TightPt1", "*pt_bin2*","Data_Model_TightPt2", "*pt_bin3*","Data_Model_TightPt3", "*pt_bin4*","Data_Model_TightPt4")
againstMuonTight3_Data.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

Data_Model_LoosePt0_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_LoosePt1_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_LoosePt2_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_LoosePt3_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_LoosePt4_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_TightPt0_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin0__Zmumu_Model_Pt0/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin0__Ztautau_Model_Pt0/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_TightPt1_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin1__Zmumu_Model_Pt1/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin1__Ztautau_Model_Pt1/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"           
)

Data_Model_TightPt2_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin2__Zmumu_Model_Pt2/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin2__Ztautau_Model_Pt2/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_TightPt3_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin3__Zmumu_Model_Pt3/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin3__Ztautau_Model_Pt3/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_TightPt4_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceHistPdf:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/mcTrue_bin0__pt_bin4__Zmumu_Model_Pt4/workspaceHistPdf:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceHistPdf:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/mcTrue_bin0__pt_bin4__Ztautau_Model_Pt4/workspaceHistPdf:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
    )
###############################################################
###############################################################
process.TnP_Muon_ID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output
    InputFileNames = cms.vstring("file:"+filePath),
    OutputFileName = cms.string("TnP_MuonToTau_MisID_Data.root"),
    InputTreeName = cms.string("fitter_tree"), 
    InputDirectoryName = cms.string("tpTree"),  
    ## Variables for binning
    Variables = cms.PSet(
        run   = cms.vstring("Run", "0", "999999", ""),
        mass   = cms.vstring("Tag-muon Mass", "60", "110", "GeV/c^{2}"),
        abseta = cms.vstring("muon |#eta|", "0.0", "0.8", ""),        
        pt  = cms.vstring("probe pT", "0", "100", ""),
        alternatLorentzVectPt = cms.vstring("probe tau pT", "20", "1500", ""),
        alternatLorentzVectEta = cms.vstring("probe tau eta", "2.4", "2.4", ""),
        tag_pt  = cms.vstring("tag pT", "0", "1500", ""),
        tag_triggerMatch = cms.vstring("Tag matched to HLT item", "0.5", "1.0", ""),
        tag_dB  = cms.vstring("dB", "0.0", "0.004", ""),
        pair_dz = cms.vstring("#Deltaz between two muons", "-0.05", "0.05", "cm"),
        pair_deltaR = cms.vstring("#DeltaR between two muons", "0.5", "5", ""),
        pair_probeMultiplicity = cms.vstring("Probe multiplicity", "1", "1", ""),
        pair_MET = cms.vstring("MET", "0", "25", ""),
        pair_MTtag = cms.vstring("MTtag", "0", "40", ""),
        pair_MTprobe = cms.vstring("MTprobe", "0", "4000", ""),
        decayModeFinding = cms.vstring("Decay mode finding", "0.5", "1.0", ""),
        #decayModeFindingNewDMs = cms.vstring("Decay mode finding NewDMs", "0.5", "1.0", ""),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vstring("Combined loose isolation", "0.5", "1.0", ""),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        againstMuonTight3 = cms.vstring("againstMuonTight3", "dummy[pass=1,fail=0]"),
        againstMuonLoose3 = cms.vstring("againstMuonLoose3", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        againstMuonLoose3 = againstMuonLoose3_Data,
        againstMuonTight3 = againstMuonTight3_Data,
    ),
    PDFs = cms.PSet(
        Data_Model_LoosePt0 = Data_Model_LoosePt0_Template,
        Data_Model_LoosePt1 = Data_Model_LoosePt1_Template,
        Data_Model_LoosePt2 = Data_Model_LoosePt2_Template,
        Data_Model_LoosePt3 = Data_Model_LoosePt2_Template,
        Data_Model_LoosePt4 = Data_Model_LoosePt2_Template,

        Data_Model_TightPt0 = Data_Model_TightPt0_Template,
        Data_Model_TightPt1 = Data_Model_TightPt1_Template,
        Data_Model_TightPt2 = Data_Model_TightPt2_Template,
        Data_Model_TightPt3 = Data_Model_TightPt3_Template,
        Data_Model_TightPt4 = Data_Model_TightPt4_Template,
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(20),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



