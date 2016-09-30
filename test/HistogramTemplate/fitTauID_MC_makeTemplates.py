import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/15_09_2016/"
filePath += "tnpZ_MC.root"

efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass","alternatLorentzVectPt", "alternatLorentzVectEta","tag_pt", "tag_triggerMatch", "tag_dB", "pair_dz", "pair_deltaR", "pair_probeMultiplicity", "pair_MET", "pair_MTtag", "pair_MTprobe", "decayModeFinding", "byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        abseta = cms.vdouble(0.0, 0.4, 0.8, 1.2, 1.7, 2.3),
        ## flags and conditions required at the denominator,
    ),
    BinToPDFmap = cms.vstring("Zmumu_Model"), ## PDF to use, as defined below
)

againstMuonLoose3_Zmumu = efficiencyPSetTemplate.clone()
againstMuonLoose3_Zmumu.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0.5,1.0))
againstMuonLoose3_Zmumu.BinToPDFmap = cms.vstring("Zmumu_Model_Eta0","*abseta_bin1*","Zmumu_Model_Eta1", "*abseta_bin2*","Zmumu_Model_Eta2", "*abseta_bin3*","Zmumu_Model_Eta3", "*abseta_bin4*","Zmumu_Model_Eta4")

againstMuonLoose3_Ztautau = efficiencyPSetTemplate.clone()
againstMuonLoose3_Ztautau.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0,0.4))
againstMuonLoose3_Ztautau.BinToPDFmap = cms.vstring("Ztautau_Model_Eta0","*abseta_bin1*","Ztautau_Model_Eta1", "*abseta_bin2*","Ztautau_Model_Eta2", "*abseta_bin3*","Ztautau_Model_Eta3", "*abseta_bin4*","Ztautau_Model_Eta4")

againstMuonTight3_Zmumu = againstMuonLoose3_Zmumu.clone()
againstMuonTight3_Zmumu.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

againstMuonTight3_Ztautau = againstMuonLoose3_Ztautau.clone()
againstMuonTight3_Ztautau.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

Zmumu_Model = cms.vstring(
    "Voigtian::signal1Fail(mass, mean1Fail[92, 85,95], width[2.495], sigma1Fail[1, 0.1,2])",
    "Voigtian::signal2Fail(mass, mean2Fail[92, 85,95], width[2.495], sigma2Fail[3, 2,10])",
    "SUM::signalFail(vFracFail[0.9, 0,1]*signal1Fail, signal2Fail)",
    
    "Voigtian::signal1Pass(mass, mean1Pass[92, 85,95], width[2.495], sigma1Pass[2, 1,10])",
    "CBShape::signal2Pass(mass, mean2Pass[92, 80,95], sigma2Pass[5, 2, 10], alpha[3, 1,5], n[2, 0,5])",
    "SUM::signalPass(vFracPass[0.9, 0,1]*signal1Pass, signal2Pass)",
    
    "Exponential::backgroundPass1(mass, lp1[-0.033])",
    "SUM::backgroundPass(0*backgroundPass1)",
    
    "Exponential::backgroundFail1(mass, lf[-0.039])",
    "SUM::backgroundFail(0*backgroundFail1)",
    
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)

Ztautau_Model = cms.vstring(
    "Voigtian::signal1Fail(mass, mean1Fail[90.95], width[2.495], sigma1Fail[1.32])",
    "SUM::signalFail(0*signal1Fail)",
    
    "Voigtian::signal1Pass(mass, mean1Pass[90.96], width[2.495], sigma1Pass[1.3])",
    "SUM::signalPass(0*signal1Pass)",

    "Gaussian::backgroundPass1(mass, mean1p[55, 50,60], sigma1p[12,10,20])",
    #"Exponential::backgroundPass2(mass, lp2[-0.168, -1,0.1])",
    "Chebychev::backgroundPass2(mass, cPass[0,-1,1])",
    "SUM::backgroundZtautauPass(vFracBkgPass[0.2, 0,1]*backgroundPass1, backgroundPass2)",
    "SUM::backgroundPass(backgroundZtautauPass)",
    
    "Gaussian::backgroundFail1(mass, mean1p[55, 50,60], sigma1p[12,10,20])",
    "Exponential::backgroundFail2(mass, lf[-0.039, -1, 0.1])",
    "SUM::backgroundFail(vFracBkgFail[0.2, 0,1]*backgroundFail1,backgroundFail2)",
    
    "efficiency[0.001,0,0.01]",
    "signalFractionInPassing[0.9]"
)
###############################################################
###############################################################
process.TnP_Muon_ID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output
    InputFileNames = cms.vstring("file:"+filePath),
    OutputFileName = cms.string("TnP_MuonToTau_MisID_MC_Templates.root"),
    InputTreeName = cms.string("fitter_tree"), 
    InputDirectoryName = cms.string("tpTree"),  
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-muon Mass", "60", "110", "GeV/c^{2}"),
        abseta = cms.vstring("muon |#eta|", "0", "2.4", ""),
        alternatLorentzVectPt = cms.vstring("probe tau pT", "20", "1500", ""),
        alternatLorentzVectEta = cms.vstring("probe tau eta", "-2.4", "2.4", ""),
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
        mcTrue = cms.vstring("Match to gen muons", "0.0", "1.0", ""),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        againstMuonTight3 = cms.vstring("againstMuonTight3", "dummy[pass=1,fail=0]"),
        againstMuonLoose3 = cms.vstring("againstMuonLoose3", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        againstMuonLoose3_Zmumu = againstMuonLoose3_Zmumu,
        againstMuonLoose3_Ztautau = againstMuonLoose3_Ztautau,

        againstMuonTight3_Zmumu = againstMuonTight3_Zmumu,
        againstMuonTight3_Ztautau = againstMuonTight3_Ztautau,
    ),
    PDFs = cms.PSet(
        Zmumu_Model_Eta0 = Zmumu_Model,
        Zmumu_Model_Eta1 = Zmumu_Model,
        Zmumu_Model_Eta2 = Zmumu_Model,
        Zmumu_Model_Eta3 = Zmumu_Model,
        Zmumu_Model_Eta4 = Zmumu_Model,

        Ztautau_Model_Eta0 = Ztautau_Model,
        Ztautau_Model_Eta1 = Ztautau_Model,
        Ztautau_Model_Eta2 = Ztautau_Model,
        Ztautau_Model_Eta3 = Ztautau_Model,
        Ztautau_Model_Eta4 = Ztautau_Model,
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(20),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(True),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



