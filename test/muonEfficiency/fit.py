import FWCore.ParameterSet.Config as cms
import sys
import copy
import itertools

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

_PT_BINS = [15, 20, 30, 40, 50, 60, 80, 120]
_ALL_PT = [15, 120]
_ETA_BINS = [2.1/10*i for i in range(11)]
# Special binning for trigger measurement, taken from Manuel Zeise's EWK talk
_ETA_BINS_TRG = [0, 0.9, 1.2, 2.1]
_VTX_BINS = [0.5, 1.5, 2.5, 4.5, 10.5]

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
    OutputFileName = cms.string("/data1/friis/ZmumuEff/efficiency_%s.root" % source),
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
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        tag_nVertices = cms.vstring("N_{vtx}", "0", "6", "")
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[pass=1,fail=0]"),
        #passing = cms.vstring("isMuon", "dummy[pass=1,fail=0]"),
        AbsIso = cms.vstring("Absolute isolation", "dummy[pass=1,fail=0]"),
        LooseAbsIso = cms.vstring("Loose absolute isolation", "dummy[pass=1,fail=0]"),
        RelIso = cms.vstring("Relative isolation", "dummy[pass=1,fail=0]"),
        LooseRelIso = cms.vstring("Loose relative isolation", "dummy[pass=1,fail=0]"),
        Glb = cms.vstring("GlobalMuon", "dummy[pass=1,fail=0]"),
        TM = cms.vstring('Tracker Muon', "dummy[pass=1,fail=0]"),
        STA = cms.vstring('Standalone Muon', "dummy[pass=1,fail=0]"),
        CustomTrigger = cms.vstring("Trigger", "dummy[pass=1,fail=0]"),
        VBTF = cms.vstring("VBTF", "dummy[pass=1,fail=0]"),
        innerTrack = cms.vstring("Inner track",  "dummy[pass=1,fail=0]"),
        outerTrack = cms.vstring("Outer track",  "dummy[pass=1,fail=0]"),
        PeriodA = cms.vstring("Period A",  "dummy[pass=1,fail=0]"),
        PeriodB = cms.vstring("Period B",  "dummy[pass=1,fail=0]"),
        PeriodC = cms.vstring("Period C",  "dummy[pass=1,fail=0]"),
        HLTMu9 = cms.vstring("HLT_Mu9",   "dummy[pass=1,fail=0]"),
        IsoMu9Mu11 = cms.vstring("HLT_IsoMu9 OR HLT_Mu11",   "dummy[pass=1,fail=0]"),
        IsoMu13Mu15 = cms.vstring("HLT_IsoMu13 OR HLT_Mu15",   "dummy[pass=1,fail=0]"),
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
    Efficiencies = cms.PSet()
)

process.TagProbeFitTreeAnalyzerSta = process.TagProbeFitTreeAnalyzer.clone(
    InputDirectoryName = cms.string("tpTreeSta"),
    OutputFileName = cms.string("/data1/friis/ZmumuEff/efficiencySta_%s.root" % source)
)

process.TagProbeFitTreeAnalyzerInner = process.TagProbeFitTreeAnalyzer.clone(
    InputDirectoryName = cms.string("tpTreeInner"),
    OutputFileName = cms.string("/data1/friis/ZmumuEff/efficiencyInner_%s.root" % source)
)

def append_pset(x_axis, eff_info, append_to):
    new_pset = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring(eff_info[1],"pass"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = cms.PSet(
        ),
        BinToPDFmap = cms.vstring('cbPlusExpo')
    )
    # Add all the selections
    for selection in eff_info[2]:
        setattr(new_pset.BinnedVariables, selection, cms.vstring("pass"))
    setattr(new_pset.BinnedVariables, x_axis[1], cms.vdouble(x_axis[2]))
    new_name = x_axis[0] + "_" + eff_info[0] + "_" + x_axis[1]
    print "Build efficiency:", new_name
    setattr(append_to, new_name, new_pset)

variables = [
    ('avg', 'pt', _ALL_PT),
    ('pt', 'pt', _PT_BINS),
    ('eta', 'abseta', _ETA_BINS),
    ('etatrig', 'abseta', _ETA_BINS_TRG),
]

iso_vars = [
    ('vtx', 'tag_nVertices', _VTX_BINS),
]

efficiencies = [
    ('iso', 'AbsIso', ['Glb', 'VBTF']),
    ('looseiso', 'LooseAbsIso', ['Glb', 'VBTF']),
    ('reliso', 'RelIso', ['Glb', 'VBTF']),
    ('loosereliso', 'LooseRelIso', ['Glb', 'VBTF']),
    ('id', 'VBTF', ['Glb']),
    ('hltMu9', ['Glb', 'VBTF', 'AbsIso']),
]

trigger_efficiencies = [
    ('trigA', 'HLTMu9', ['Glb', 'VBTF', 'AbsIso', 'PeriodA']),
    ('trigB', 'IsoMu9Mu11', ['Glb', 'VBTF', 'AbsIso', 'PeriodB']),
    ('trigC', 'IsoMu13Mu15', ['Glb', 'VBTF', 'AbsIso', 'PeriodC']),
]


# FIXME is TM
inner_efficiencies = [
    ('innerTrack', 'innerTrack', ['outerTrack']),
]

# FIXME
outer_efficiencies = [
    ('standAlone', 'outerTrack', ['innerTrack']),
    ('linking', 'Glb', ['outerTrack', 'innerTrack']),
]

for eff in efficiencies:
    for variable in variables + iso_vars:
        append_pset(variable, eff, process.TagProbeFitTreeAnalyzer.Efficiencies)

if not source_info['mc']:
    for eff in trigger_efficiencies:
        for variable in variables:
            append_pset(variable, eff,
                        process.TagProbeFitTreeAnalyzer.Efficiencies)

for eff in inner_efficiencies:
    for variable in variables:
        append_pset(variable, eff,
                    process.TagProbeFitTreeAnalyzerSta.Efficiencies)

for eff in outer_efficiencies:
    for variable in variables:
        append_pset(variable, eff,
                    process.TagProbeFitTreeAnalyzerInner.Efficiencies)

if source_info['mc']:
    params = process.TagProbeFitTreeAnalyzer.Efficiencies.parameterNames_()
    for param in params:
        object = getattr(process.TagProbeFitTreeAnalyzer.Efficiencies, param)
        if isinstance(object, cms.PSet):
            print "Updating param:", param
            object.BinnedVariables.mcTrue = cms.vstring("pass")

process.fit = cms.Path(
    process.TagProbeFitTreeAnalyzer*
    process.TagProbeFitTreeAnalyzerSta*
    process.TagProbeFitTreeAnalyzerInner
)
