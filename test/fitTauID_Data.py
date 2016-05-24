import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/SingleMuon_Run2015D_16Dec2015_v1_v11/SingleMuon/SingleMuon_Run2015D_16Dec2015_v1_v11/160520_134214/0000/"

filePath+="/tnpZ_Data.root"

efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        #pt     = cms.vdouble( 10, 100 ),
        abseta = cms.vdouble( 0.0, 1.2, 1.7, 2.3),
        ## flags and conditions required at the denominator,
        pair_probeMultiplicity = cms.vdouble(1.0,1.1),
        pair_dz = cms.vdouble(-0.05,0.05),             
        pair_deltaR = cms.vdouble(0.5,5.),
        decayModeFinding = cms.vdouble(0.5,1.1),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vdouble(0.5,1.1),
    ),
    BinToPDFmap = cms.vstring("cbPlusPolyEta0","*abseta_bin1*","cbPlusPolyEta1", "*abseta_bin2*","cbPlusPolyEta2")
)

againstMuonLoose3_pt_abseta = efficiencyPSetTemplate.clone()

againstMuonTight3_pt_abseta = efficiencyPSetTemplate.clone()
againstMuonTight3_pt_abseta.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")
againstMuonTight3_pt_abseta.BinToPDFmap = cms.vstring("cbPlusPolyTightEta0","*abseta_bin1*","cbPlusPolyTightEta1", "*abseta_bin2*","cbPlusPolyTightEta2")

process.TnP_Muon_ID = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    ## Input, output
    InputFileNames = cms.vstring("file:"+filePath),
    OutputFileName = cms.string("TnP_MuonToTau_MisID_Data.root"),
    InputTreeName = cms.string("fitter_tree"), 
    InputDirectoryName = cms.string("tpTree"),  
    ## Variables for binning
    Variables = cms.PSet(
        mass   = cms.vstring("Tag-muon Mass", "70", "120", "GeV/c^{2}"),
        pt     = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("muon |#eta|", "0", "2.4", ""),
        pair_dz = cms.vstring("#Deltaz between two muons", "-100", "100", "cm"),
        pair_deltaR = cms.vstring("#DeltaR between two muons", "0", "10", ""),
        pair_probeMultiplicity = cms.vstring("Probe multiplicity", "0", "5", ""),
        decayModeFinding = cms.vstring("Decay mode finding", "-0.1", "1.1", ""),
        decayModeFindingNewDMs = cms.vstring("Decay mode finding NewDMs", "-0.1", "1.1", ""),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vstring("Combined loose isolation", "-0.1", "1.1", ""),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        againstMuonTight3 = cms.vstring("againstMuonTight3", "dummy[pass=1,fail=0]"),
        againstMuonLoose3 = cms.vstring("againstMuonLoose3", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        againstMuonLoose3_pt_abseta = againstMuonLoose3_pt_abseta,
        againstMuonTight3_pt_abseta = againstMuonTight3_pt_abseta,
    ),
    ## PDF for signal and background (double voigtian + exponential background)
    PDFs = cms.PSet(
        cbPlusPolyEta0 = cms.vstring(
            "Voigtian::signal1Fail(mass, mean1Fail[91], width[2.495], sigma1Fail[1.137])",
            "Voigtian::signal2Fail(mass, mean2Fail[86.3], width[2.495], sigma2Fail[4.71])",
            "SUM::signalFail(vFracFail[0.926]*signal1Fail, signal2Fail)",
            
            "Voigtian::signal1Pass(mass, mean1Pass[90.93], width[2.495], sigma1Pass[1.09])",
            "CBShape::signal2Pass(mass, mean2Pass[83.2], sigma2Pass[7.0], alpha[3], n[1])",
            "SUM::signalPass(vFracPass[0.83]*signal1Pass, signal2Pass)",
            
            "Exponential::backgroundPass1(mass, lp1[-0.033])",
            "Exponential::backgroundPass2(mass, lp2[-0.168])",
            "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
            "SUM::backgroundPass12(vFracBkgPass[0.2]*backgroundPass1, backgroundPass2)",
            "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
            
            "Exponential::backgroundFail1(mass, lf[-0.039])",
            "SUM::backgroundFail(backgroundFail1)",
            
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
        ),
        cbPlusPolyEta1 = cms.vstring(
            "Voigtian::signal1Fail(mass, mean1Fail[90.95], width[2.495], sigma1Fail[1.32])",
            "Voigtian::signal2Fail(mass, mean2Fail[85.98], width[2.495], sigma2Fail[4.9])",
            "SUM::signalFail(vFracFail[0.927]*signal1Fail, signal2Fail)",
            
            "Voigtian::signal1Pass(mass, mean1Pass[90.96], width[2.495], sigma1Pass[1.3])",
            "CBShape::signal2Pass(mass, mean2Pass[85], sigma2Pass[4.4], alpha[4], n[0])",
            "SUM::signalPass(vFracPass[0.91]*signal1Pass, signal2Pass)",
            
            "Exponential::backgroundPass1(mass, lp1[-0.033])",
            "Exponential::backgroundPass2(mass, lp2[-0.168])",
            "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
            "SUM::backgroundPass12(vFracBkgPass12[0.2]*backgroundPass1, backgroundPass2)",
            "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
            
            "Exponential::backgroundFail1(mass, lf[-0.039])",
            "SUM::backgroundFail(backgroundFail1)",
            
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
        ),
        cbPlusPolyEta2 = cms.vstring(
            "Voigtian::signal1Fail(mass, mean1Fail[90.92], width[2.495], sigma1Fail[1.617])",
            "Voigtian::signal2Fail(mass, mean2Fail[86.02], width[2.495], sigma2Fail[5.18])",
            "SUM::signalFail(vFracFail[0.93]*signal1Fail, signal2Fail)",
            
            "Voigtian::signal1Pass(mass, mean1Pass[90.7], width[2.495], sigma1Pass[1.6])",
            "CBShape::signal2Pass(mass, mean2Pass[80.0], sigma2Pass[7], alpha[3], n[0.1])",
            "SUM::signalPass(vFracPass[0.78]*signal1Pass, signal2Pass)",
            
            "Exponential::backgroundPass1(mass, lp1[0.02])",
            "Exponential::backgroundPass2(mass, lp2[-0.140])",
            "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
            "SUM::backgroundPass12(vFracBkgPass[0.17]*backgroundPass1, backgroundPass2)",
            "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
            
            "Exponential::backgroundFail1(mass, lf[-0.037])",
            "SUM::backgroundFail(backgroundFail1)",
            
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
        ),
         cbPlusPolyTightEta0 = cms.vstring(
            "Voigtian::signal1Fail(mass, mean1Fail[90.98], width[2.495], sigma1Fail[1.27])",
            "Voigtian::signal2Fail(mass, mean2Fail[84.88], width[2.495], sigma2Fail[4.67])",
            "SUM::signalFail(vFracFail[0.944]*signal1Fail, signal2Fail)",
            
            "Voigtian::signal1Pass(mass, mean1Pass[90.79], width[2.495], sigma1Pass[1.36])",
            "CBShape::signal2Pass(mass, mean2Pass[82.3], sigma2Pass[5.8], alpha[5], n[1])",
            "SUM::signalPass(vFracPass[0.73]*signal1Pass, signal2Pass)",
            
            "Exponential::backgroundPass1(mass, lp1[-0.015])",
            "Exponential::backgroundPass2(mass, lp2[-0.176])",
            "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
            "SUM::backgroundPass12(vFracBkgPass[0.150]*backgroundPass1, backgroundPass2)",
             "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
            
            "Exponential::backgroundFail1(mass, lf[-0.039])",
            "SUM::backgroundFail(backgroundFail1)",
            
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
        ),
         cbPlusPolyTightEta1 = cms.vstring(
             "Voigtian::signal1Fail(mass, mean1Fail[90.949], width[2.495], sigma1Fail[1.427])",
             "Voigtian::signal2Fail(mass, mean2Fail[85.48], width[2.495], sigma2Fail[4.86])",
             "SUM::signalFail(vFracFail[0.936]*signal1Fail, signal2Fail)",
            
             "Voigtian::signal1Pass(mass, mean1Pass[91.5], width[2.495], sigma1Pass[1])",
             "CBShape::signal2Pass(mass, mean2Pass[88], sigma2Pass[6.1], alpha[5], n[0])",
             "SUM::signalPass(vFracPass[0.63]*signal1Pass, signal2Pass)",
             
             "Exponential::backgroundPass1(mass, lp1[0.1])",
             "Exponential::backgroundPass2(mass, lp2[-0.153])",
             "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
             "SUM::backgroundPass12(vFracBkgPass[0.06]*backgroundPass1, backgroundPass2)",
             "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
             
             "Exponential::backgroundFail1(mass, lf[-0.039, -1,0])",
             "SUM::backgroundFail(backgroundFail1)",
             
             "efficiency[0.001,0,0.01]",
             "signalFractionInPassing[0.9]"
        ),
        cbPlusPolyTightEta2 = cms.vstring(
             "Voigtian::signal1Fail(mass, mean1Fail[90.915], width[2.495], sigma1Fail[1.615])",
             "Voigtian::signal2Fail(mass, mean2Fail[86.01], width[2.495], sigma2Fail[5.17])",
             "SUM::signalFail(vFracFail[0.923]*signal1Fail, signal2Fail)",
            
             "Voigtian::signal1Pass(mass, mean1Pass[90.5], width[2.495], sigma1Pass[0.5])",
             "CBShape::signal2Pass(mass, mean2Pass[85], sigma2Pass[7], alpha[4], n[0])",
             "SUM::signalPass(vFracPass[0.51]*signal1Pass, signal2Pass)",
             
             "Exponential::backgroundPass1(mass, lp1[0.02])",
             "Exponential::backgroundPass2(mass, lp2[-0.141])",
             "Exponential::backgroundPass3(mass, lp3[-0.235,-1,0])",
             "SUM::backgroundPass12(vFracBkgPass[0.17]*backgroundPass1, backgroundPass2)",
             "SUM::backgroundPass(vFracBkgPass123[0.2,0,1]*backgroundPass12, backgroundPass3)",
             
             "Exponential::backgroundFail1(mass, lf[-0.037])",
             "SUM::backgroundFail(backgroundFail1)",
             
             "efficiency[0.001,0,0.01]",
             "signalFractionInPassing[0.9]"
        ),
    ),
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(15),
    saveDistributionsPlot = cms.bool(True),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(False),
)

process.p1 = cms.Path(process.TnP_Muon_ID)

