import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/18_01_2017/"
filePath += "tnpZ_Data.root"

efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass","alternatLorentzVectPt", "alternatLorentzVectEta", "tag_pt", "tag_triggerMatch", "tag_dB", "pair_dz", "pair_deltaR", "pair_probeMultiplicity", "pair_BestZ", "pair_MET", "pair_MTtag", "pair_MTprobe", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        abseta = cms.vdouble(0.0, 0.4, 0.8, 1.2, 1.7, 2.3),
        #run = cms.vdouble(271036, 274422, 275125),
        ## flags and conditions required at the denominator,
    ),
    BinToPDFmap = cms.vstring("Zll_Model"), ## PDF to use, as defined below
)

againstMuonLoose3_Data = efficiencyPSetTemplate.clone()
againstMuonLoose3_Data.BinToPDFmap = cms.vstring("Data_Model_LooseEta0","*abseta_bin1*","Data_Model_LooseEta1", "*abseta_bin2*","Data_Model_LooseEta2", "*abseta_bin3*","Data_Model_LooseEta3", "*abseta_bin4*","Data_Model_LooseEta4")

againstMuonTight3_Data = efficiencyPSetTemplate.clone()
againstMuonTight3_Data.BinToPDFmap = cms.vstring("Data_Model_TightEta0","*abseta_bin1*","Data_Model_TightEta1", "*abseta_bin2*","Data_Model_TightEta2", "*abseta_bin3*","Data_Model_TightEta3", "*abseta_bin4*","Data_Model_TightEta4")
againstMuonTight3_Data.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

Data_Model_LooseEta0_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0.1,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_LooseEta1_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.15,0.1,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"           
)

Data_Model_LooseEta2_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0.1,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_LooseEta3_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0.1,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_LooseEta4_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonLoose3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.2,0.1,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_TightEta0_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin0__mcTrue_bin0__Zmumu_Model_Eta0/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin0__mcTrue_bin0__Ztautau_Model_Eta0/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
            )

Data_Model_TightEta1_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin1__mcTrue_bin0__Zmumu_Model_Eta1/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin1__mcTrue_bin0__Ztautau_Model_Eta1/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"           
)

Data_Model_TightEta2_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin2__mcTrue_bin0__Zmumu_Model_Eta2/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin2__mcTrue_bin0__Ztautau_Model_Eta2/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_TightEta3_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin3__mcTrue_bin0__Zmumu_Model_Eta3/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin3__mcTrue_bin0__Ztautau_Model_Eta3/workspaceFixedParams:backgroundZtautauPass",
    "Exponential::backgroundWJetPass(mass, lp3[-0.235,-1,0])",
    "SUM::backgroundPass(vFracWJetBkgPass[0.1,0,1]*backgroundWJetPass, backgroundZtautauPass)",
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Data_Model_TightEta4_Template = cms.vstring(
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalPass",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Zmumu/abseta_bin4__mcTrue_bin0__Zmumu_Model_Eta4/workspaceFixedParams:signalFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundFail",
    "#import TnP_MuonToTau_MisID_MC_Templates.root:tpTree/againstMuonTight3_Ztautau/abseta_bin4__mcTrue_bin0__Ztautau_Model_Eta4/workspaceFixedParams:backgroundZtautauPass",
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
        mass   = cms.vstring("Tag-muon Mass", "70", "120", "GeV/c^{2}"),
        abseta = cms.vstring("muon |#eta|", "0", "2.4", ""),
        alternatLorentzVectPt = cms.vstring("probe tau pT", "20", "1500", ""),
        alternatLorentzVectEta = cms.vstring("probe tau eta", "-2.4", "2.4", ""),
        tag_pt  = cms.vstring("tag pT", "0", "1500", ""),
        tag_triggerMatch = cms.vstring("Tag matched to HLT item", "0.5", "1.0", ""),
        tag_dB  = cms.vstring("dB", "0.0", "0.004", ""),
        pair_dz = cms.vstring("#Deltaz between two muons", "-0.05", "0.05", "cm"),
        pair_deltaR = cms.vstring("#DeltaR between two muons", "0.5", "5", ""),
        pair_probeMultiplicity = cms.vstring("Probe multiplicity", "1", "1", ""),
        pair_BestZ = cms.vstring("pair closest to Z mass", "0.5", "1.0", ""),
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
        Data_Model_LooseEta0 = Data_Model_LooseEta0_Template,
        Data_Model_LooseEta1 = Data_Model_LooseEta1_Template,
        Data_Model_LooseEta2 = Data_Model_LooseEta2_Template,
        Data_Model_LooseEta3 = Data_Model_LooseEta2_Template,
        Data_Model_LooseEta4 = Data_Model_LooseEta2_Template,

        Data_Model_TightEta0 = Data_Model_TightEta0_Template,
        Data_Model_TightEta1 = Data_Model_TightEta1_Template,
        Data_Model_TightEta2 = Data_Model_TightEta2_Template,
        Data_Model_TightEta3 = Data_Model_TightEta3_Template,
        Data_Model_TightEta4 = Data_Model_TightEta4_Template,
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(50),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



