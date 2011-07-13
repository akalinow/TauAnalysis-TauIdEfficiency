import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('/data1/veelken/tmp/muonPtGt20/V4b/compTauIdEffFinalNumbers_output.root')
)

process.makeTauIdEffFinalPlots = cms.PSet(

    tauIds = cms.VPSet(
        cms.PSet(
            name = cms.string('tauDiscrHPScombLooseDBcorr'),
            legendEntry = cms.string("HPS comb. Loose"),
            markerStyleData = cms.uint32(20),
            markerStyleSim = cms.uint32(24),
            color = cms.uint32(418)
        ),
        cms.PSet(
            name = cms.string('tauDiscrHPScombMediumDBcorr'),
            legendEntry = cms.string("HPS comb. Medium"),
            markerStyleData = cms.uint32(21),
            markerStyleSim = cms.uint32(25),
            color = cms.uint32(807)
        ),
        cms.PSet(
            name = cms.string('tauDiscrHPScombTightDBcorr'),
            legendEntry = cms.string("HPS comb. Tight"),
            markerStyleData = cms.uint32(22),
            markerStyleSim = cms.uint32(26),
            color = cms.uint32(618)
        )
    ),

    fitVariables = cms.vstring(
        'diTauVisMass',
        #'diTauVisMassFromJet' # CV: diTauVisMass always computed from PFJet momenta if using PAT-tuple workflow
    ),

    xAxisBinning = cms.vdouble(0.0, 4.5, 6.5, 8.5, 20.5),
    xAxisTitle = cms.string("Num. Vertices"),

    values = cms.VPSet(
        cms.PSet(
            directory = cms.string(''),
            xBinCenter = cms.double(2.0),
        ),
        cms.PSet(
            directory = cms.string(''),
            xBinCenter = cms.double(5.5),
        ),
        cms.PSet(
            directory = cms.string(''),
            xBinCenter = cms.double(7.5),
        ),
        cms.PSet(
            directory = cms.string(''),
            xBinCenter = cms.double(14.0),
        ),
    ),

    outputFileName = cms.string('/data1/veelken/tmp/muonPtGt20/V4b/makeTauIdEffFinalPlots_HPScombined.eps')
)
