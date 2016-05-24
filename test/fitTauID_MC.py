import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v5_ext4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v5_ext4/160426_105344/0000/"
#filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/WJetsToLNu_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_v5/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/WJetsToLNu_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_v5/160426_124914/0000/"

#v11
filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v11_ext4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v11_ext4/160520_133905/0000/"

#filePath = "./"
filePath += "tnpZ_MC.root"


efficiencyPSetTemplate = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    EfficiencyCategoryAndState = cms.vstring("againstMuonLoose3", "pass"), ## Numerator definition
    BinnedVariables = cms.PSet(
        ## Binning in continuous variables
        #pt     = cms.vdouble( 10, 100 ),
        abseta = cms.vdouble(0.0, 1.2, 1.21),
        ## flags and conditions required at the denominator,
        #tag_pt = cms.vdouble(20,100),
        tag_triggerMatch = cms.vdouble(1.0,1.1),
        pair_probeMultiplicity = cms.vdouble(1.0,1.1),
        pair_dz = cms.vdouble(-0.05,0.05),             
        pair_deltaR = cms.vdouble(0.5,5.),
        decayModeFinding = cms.vdouble(0.5,1.1),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vdouble(0.5,1.1),
    ),
    BinToPDFmap = cms.vstring("cbPlusPolyEta0"), ## PDF to use, as defined below
)

againstMuonLoose3_pt_abseta0 = efficiencyPSetTemplate.clone()

againstMuonLoose3_pt_abseta0_mcTrue = efficiencyPSetTemplate.clone()
againstMuonLoose3_pt_abseta0_mcTrue.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0.5,1.1))
#againstMuonLoose3_pt_abseta0_mcTrue.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(-0.1,0.4))

againstMuonTight3_pt_abseta0 = efficiencyPSetTemplate.clone()
againstMuonTight3_pt_abseta0.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

againstMuonTight3_pt_abseta0_mcTrue = efficiencyPSetTemplate.clone()
againstMuonTight3_pt_abseta0_mcTrue.BinnedVariables._Parameterizable__addParameter("mcTrue",cms.vdouble(0.5,1.1))
againstMuonTight3_pt_abseta0_mcTrue.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")


againstMuonLoose3_pt_abseta12 = againstMuonLoose3_pt_abseta0.clone()
againstMuonLoose3_pt_abseta12.BinnedVariables.abseta = cms.vdouble(1.2, 1.7, 2.3)
againstMuonLoose3_pt_abseta12.BinToPDFmap = cms.vstring("cbPlusPolyEta1","*abseta_bin0*","cbPlusPolyEta1", "*abseta_bin1*","cbPlusPolyEta2")

againstMuonLoose3_pt_abseta12_mcTrue = againstMuonLoose3_pt_abseta0_mcTrue.clone()
againstMuonLoose3_pt_abseta12_mcTrue.BinnedVariables.abseta = cms.vdouble(1.2, 1.7, 2.3)
againstMuonLoose3_pt_abseta12_mcTrue.BinToPDFmap = cms.vstring("cbPlusPolyEta1","*abseta_bin0*","cbPlusPolyEta1", "*abseta_bin1*","cbPlusPolyEta2")
##againstMuonLoose3_pt_abseta12_mcTrue.BinnedVariables.abseta = cms.vdouble(1.7, 2.3)
##againstMuonLoose3_pt_abseta12_mcTrue.BinToPDFmap = cms.vstring("cbPlusPolyEta2")

againstMuonTight3_pt_abseta12 = againstMuonTight3_pt_abseta0.clone()
againstMuonTight3_pt_abseta12.BinnedVariables.abseta = cms.vdouble(1.2, 1.7, 2.3)

againstMuonTight3_pt_abseta12_mcTrue = againstMuonTight3_pt_abseta0_mcTrue.clone()
againstMuonTight3_pt_abseta12_mcTrue.BinnedVariables.abseta = cms.vdouble(1.2, 1.7, 2.3)

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
        mass   = cms.vstring("Tag-muon Mass", "75", "110", "GeV/c^{2}"),
        #pt     = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("muon |#eta|", "0", "2.4", ""),
        tag_triggerMatch = cms.vstring("Tag matched to HLT item", "0", "1.1", ""),
        #tag_pt = cms.vstring("Tag pt", "20", "100", ""),
        pair_dz = cms.vstring("#Deltaz between two muons", "-100", "100", "cm"),
        pair_deltaR = cms.vstring("#DeltaR between two muons", "0", "10", ""),
        pair_probeMultiplicity = cms.vstring("Probe multiplicity", "0", "5", ""),
        decayModeFinding = cms.vstring("Decay mode finding", "-0.1", "1.1", ""),
        #decayModeFindingNewDMs = cms.vstring("Decay mode finding NewDMs", "-0.1", "1.1", ""),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.vstring("Combined loose isolation", "-0.1", "1.1", ""),
        mcTrue = cms.vstring("Match to gen muons", "-0.1", "1.1", ""),
    ),
    ## Flags you want to use to define numerator and possibly denominator
    Categories = cms.PSet(
        againstMuonTight3 = cms.vstring("againstMuonTight3", "dummy[pass=1,fail=0]"),
        againstMuonLoose3 = cms.vstring("againstMuonLoose3", "dummy[pass=1,fail=0]"),
    ),
    ## What to fit
    Efficiencies = cms.PSet(
        againstMuonLoose3_pt_abseta0 = againstMuonLoose3_pt_abseta0,        
        againstMuonLoose3_pt_abseta0_mcTrue = againstMuonLoose3_pt_abseta0_mcTrue,
        #againstMuonTight3_pt_abseta0 = againstMuonTight3_pt_abseta0,
        #againstMuonTight3_pt_abseta0_mcTrue = againstMuonTight3_pt_abseta0_mcTrue,
        
        againstMuonLoose3_pt_abseta12 = againstMuonLoose3_pt_abseta12,        
        againstMuonLoose3_pt_abseta12_mcTrue = againstMuonLoose3_pt_abseta12_mcTrue,
        #againstMuonTight3_pt_abseta12 = againstMuonTight3_pt_abseta12,
        #againstMuonTight3_pt_abseta12_mcTrue = againstMuonTight3_pt_abseta12_mcTrue 
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
            "SUM::backgroundPass(vFracBkgPass[0.2]*backgroundPass1, backgroundPass2)",
            
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
            "SUM::backgroundPass(vFracBkgPass[0.2]*backgroundPass1, backgroundPass2)",
            
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
            "SUM::backgroundPass(vFracBkgPass[0.17]*backgroundPass1, backgroundPass2)",
            
            "Exponential::backgroundFail1(mass, lf[-0.037])",
            "SUM::backgroundFail(backgroundFail1)",
            
            "efficiency[0.001,0,0.01]",
            "signalFractionInPassing[0.9]"
        )
    ),    
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(15),
    saveDistributionsPlot = cms.bool(False),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(False),
)

process.p1 = cms.Path(process.TnP_Muon_ID)



