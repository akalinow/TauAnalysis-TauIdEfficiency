import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/SingleMuon_Run2015D_16Dec2015_v1_v3/SingleMuon/SingleMuon_Run2015D_16Dec2015_v1_v3/160418_184453/0000/"
filePath = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/SingleMuon_Run2015D_16Dec2015_v1_v5/SingleMuon/SingleMuon_Run2015D_16Dec2015_v1_v5/160426_125121/0000/"
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
    BinToPDFmap = cms.vstring("vpvPlusExpo"), ## PDF to use, as defined below
)

againstMuonLoose3_pt_abseta = efficiencyPSetTemplate.clone()

againstMuonTight3_pt_abseta = efficiencyPSetTemplate.clone()
againstMuonTight3_pt_abseta.EfficiencyCategoryAndState = cms.vstring("againstMuonTight3", "pass")

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
        vpvPlusExpo = cms.vstring(
            "Voigtian::signal1(mass, mean1[90,80,100], width[2.495], sigma1[2,1,3])",
            "Voigtian::signal2(mass, mean2[90,80,100], width,        sigma2[4,2,10])",
            "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
            "Exponential::backgroundPass1(mass, lp1[-0.1,-1,0.1])",
            "Exponential::backgroundPass2(mass, lp2[-0.1,-1,0.1])",
            "SUM::backgroundPass(vFracBkg[0.8,0,1]*backgroundPass1, backgroundPass2)",
            "Exponential::backgroundFail(mass, lf[-0.1,-1,0.1])",
            "efficiency[0.01,0,0.1]",
            "signalFractionInPassing[0.9]"
        ),      
    ),
    ## How to do the fit
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(True),
    NumCPU = cms.uint32(1), ## leave to 1 for now, RooFit gives funny results otherwise
    SaveWorkspace = cms.bool(False),
)

process.p1 = cms.Path(process.TnP_Muon_ID)

