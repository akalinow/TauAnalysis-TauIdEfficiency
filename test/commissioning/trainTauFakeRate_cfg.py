#!/usr/bin/env cmsRun
'''

Tau Fake Rate MVA trainer

Author: Evan K. Friis (UC Davis)

'''

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

# Register options
options.register(
    'xml', '',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "XML file with MVA configuration")

options.register(
    'numeratorfile', '',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Numerator file")

options.register(
    'denominatorfile', '',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Denominator file")

options.register(
    'output', '',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "Output sqlite file")

options.parseArguments()

process = cms.Process("TrainFakeRateKNN")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# Input files
process.source = cms.Source("EmptySource")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500), # every 100th only
    limit = cms.untracked.int32(-1)       # or limit to 10 printouts...
))
process.MessageLogger.statistics.append('cout')

#######################################################
# Database BS
#######################################################

process.PoolDBOutputService = cms.Service(
    "PoolDBOutputService",
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet( messageLevel = cms.untracked.int32(0) ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:%s' % options.outputFile.replace(
        '.root', '')),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('TauTagMVAComputerRcd'),
        tag = cms.string('FakeRate')
    ))
)

process.MVATrainerSave = cms.EDAnalyzer(
    "TauMVATrainerSave",
    toPut = cms.vstring('train'),
    toCopy = cms.vstring()
)

process.looper = cms.Looper(
    "TauMVATrainerLooper",
    trainers = cms.VPSet(cms.PSet(
        calibrationRecord = cms.string("train"),
        saveState = cms.untracked.bool(True),
        trainDescription = cms.untracked.string(options.xml),
        loadState = cms.untracked.bool(False),
        doMonitoring = cms.bool(True),
    ))
)

process.train = cms.EDAnalyzer(
    "TauFakeRateTrainer",
    numerator = cms.string(options.numeratorfile),
    denominator = cms.string(options.denominatorfile),
)

process.main = cms.Path(process.train)

process.outpath = cms.EndPath(process.MVATrainerSave)

process.sched = cms.Schedule(process.main, process.outpath)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

