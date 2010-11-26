import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

_PT_BINS = [10, 15, 20, 30, 40, 50, 60, 80, 120]
_ETA_BINS = [2.1/10*i for i in range(11)]

sources = {
    'data' : {
        'file' : "/data1/friis/ZmumuEff/tagAndProbe_data.root",
        'mc' : False,
    },
    'mc' : {
        'file' : "/data1/friis/ZmumuEff/tagAndProbe_mc.root",
        'mc' : True,
    },
    'mcpu' : {
        'file' : "/data1/friis/ZmumuEff/tagAndProbe_mcpu.root",
        'mc' : True,
    },
    'mcboth' : {
        'file' : [
            "/data1/friis/ZmumuEff/tagAndProbe_mcpu.root",
            "/data1/friis/ZmumuEff/tagAndProbe_mc.root",
        ],
        'mc' : True,
    },
}

source = sys.argv[2]
source_info = sources[source]

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(source_info['file']),
    InputDirectoryName = cms.string("tpTree"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("efficiency_%s.root" % source),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    fixVars = cms.vstring("mean"),

    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", "")
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        #passing = cms.vstring("isMuon", "dummy[pass=1,fail=0]"),
        AbsIso = cms.vstring("Passing isolation", "dummy[pass=1,fail=0]"),
        AbsIsoPFnoPileUp = cms.vstring("Passing isolation (PU removal)", "dummy[pass=1,fail=0]"),
        Glb = cms.vstring("GlobalMuon", "dummy[pass=1,fail=0]"),
        CustomTrigger = cms.vstring("Trigger", "dummy[pass=1,fail=0]"),
        VBTF = cms.vstring("VBTF", "dummy[pass=1,fail=0]"),
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values
    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
            "Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
            "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
            "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        gaussPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        cbPlusExpo = cms.vstring(
            "CBShape::signal(mass, mean[90,80,100], sigma[2.5, 0.5, 10], alpha[1, -10, 10], n[3, 0, 10])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting.
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        pt_absIso = cms.PSet(
            #specifies the efficiency of which category and state to measure
            EfficiencyCategoryAndState = cms.vstring("AbsIso","pass"),
            #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            UnbinnedVariables = cms.vstring("mass"),
            #specifies the binning of parameters
            BinnedVariables = cms.PSet(
                # Require that the muon was identified as a global muon and
                # passed the ID cuts
                Glb = cms.vstring("pass"),
                VBTF = cms.vstring("pass"),
                pt = cms.vdouble(_PT_BINS)
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
            BinToPDFmap = cms.vstring("cbPlusExpo")
        ),
        #eta_absIso = cms.PSet(
            ##specifies the efficiency of which category and state to measure
            #EfficiencyCategoryAndState = cms.vstring("AbsIso","pass"),
            ##specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon and
                ## passed the ID cuts
                #Glb = cms.vstring("pass"),
                #VBTF = cms.vstring("pass"),
                #abseta = cms.vdouble(_ETA_BINS)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
        #),
        #pt_vbtf = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("VBTF", "pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon and
                ## passed the ID cuts
                #Glb = cms.vstring("pass"),
                #pt = cms.vdouble(_PT_BINS)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
        #),
        #eta_vbtf = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("VBTF", "pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon and
                ## passed the ID cuts
                #Glb = cms.vstring("pass"),
                #abseta = cms.vdouble(_ETA_BINS)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
        #),
        #pt_trigger = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("CustomTrigger", "pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon and
                ## passed the ID cuts
                #Glb = cms.vstring("pass"),
                #VBTF = cms.vstring("pass"),
                #AbsIso = cms.vstring("pass"),
                #pt = cms.vdouble(_PT_BINS)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
        #),
        #eta_trigger = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("CustomTrigger", "pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon and
                ## passed the ID cuts
                #Glb = cms.vstring("pass"),
                #VBTF = cms.vstring("pass"),
                #AbsIso = cms.vstring("pass"),
                #abseta = cms.vdouble(_ETA_BINS)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusExpo")
        #),
        #pt_trigger = cms.PSet(
            ##specifies the efficiency of which category and state to measure
            #EfficiencyCategoryAndState = cms.vstring("CustomTrigger","pass"),
            ##specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon passed all offline cuts
                #Glb = cms.vstring("pass"),
                #VBTF = cms.vstring("pass"),
                #AbsIso = cms.vstring("pass"),
                ## Require that the muon fired the trigger?
                #pt = cms.vdouble(15, 100),
                ##pt = cms.vdouble(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
        #pt_id = cms.PSet(
            ##specifies the efficiency of which category and state to measure
            #EfficiencyCategoryAndState = cms.vstring("VBTF","pass"),
            ##specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            #UnbinnedVariables = cms.vstring("mass"),
            ##specifies the binning of parameters
            #BinnedVariables = cms.PSet(
                ## Require that the muon was identified as a global muon
                #Glb = cms.vstring("true"),
                ## Require that the muon fired the trigger?
                #pt = cms.vdouble(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
            #),
            ##first string is the default followed by binRegExp - PDFname pairs
            #BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
        #pt_mcTrue = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            #BinnedVariables = cms.PSet(
                #mcTrue = cms.vstring("true"),
                #pt = cms.vdouble(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)
            #),
            ##unspecified binToPDFmap means no fitting
            #BinToPDFmap = cms.vstring()
        #),
        #pt_eta = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            #BinnedVariables = cms.PSet(
                #pt = cms.vdouble(20, 30, 40, 50, 60, 70, 80, 90, 100, 110,120),
                #eta = cms.vdouble(-2.4,-1.2, 0.0, 1.2, 2.4)
            #),
            #BinToPDFmap = cms.vstring("gaussPlusLinear")
        #),
        #pt_eta_mcTrue = cms.PSet(
            #EfficiencyCategoryAndState = cms.vstring("passing","pass"),
            #UnbinnedVariables = cms.vstring("mass"),
            #BinnedVariables = cms.PSet(
                #mcTrue = cms.vstring("true"),
                #pt = cms.vdouble(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120),
                #eta = cms.vdouble(-2.4, -1.2, 0.0, 1.2, 2.4)
            #),
            #BinToPDFmap = cms.vstring()
        #)
    )
)

if source_info['mc']:
    params = process.TagProbeFitTreeAnalyzer.Efficiencies.parameterNames_()
    for param in params:
        object = getattr(process.TagProbeFitTreeAnalyzer.Efficiencies, param)
        if isinstance(object, cms.PSet):
            print "Updating param:", param
            object.BinnedVariables.mcTrue = cms.vstring("true")

process.fit = cms.Path(process.TagProbeFitTreeAnalyzer)
