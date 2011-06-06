#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# -- Conditions
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                                'file:/nfs/data4/verzetti/tagprobe/testNtuples/tauIdEffMeasTestPattupleV3.root'
)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load('TauAnalysis.TauIdEfficiency.plotMacroConfigTauIDs') #TauID PSet
process.load('TauAnalysis.TauIdEfficiency.plotMacroConfigVariables') #Variable PSet


process.makePlotsChannel = cms.EDAnalyzer('HistoPlotter',
                                          outputFile = cms.string('testPlotter.root'),
                                          channel = cms.string('data'),
                                          evtCounter = cms.InputTag("totalEventsProcessed"),
                                          sysUncertainties =  cms.vstring("CENTRAL_VALUE",
                                                                          "SysTauJetEnUp",
                                                                          "SysTauJetEnDown",
                                                                          "SysJetEnUp",
                                                                          "SysJetEnDown"
#                                                                          "SysZllRecoilCorrectionUp",
#                                                                          "SysZllRecoilCorrectionDown"
                                                                          ),#vstring
                                          regions = cms.vstring( "ABCD",
                                                                "A",
                                                                "B",
                                                                "B1",  # QCD enriched control region (SS, loose muon isolation, Mt && Pzeta cuts applied)
                                                                "C",
                                                                "C1",
                                                                "C1p",
                                                                "C1f",
                                                                "C2",
                                                                "C2p",
                                                                "C2f",
                                                                "D",   # generic background control region (SS, tight muon isolation)
                                                                "D1",
                                                                "D1p",
                                                                "D1f",
                                                                "D2",
                                                                "D2p",
                                                                "D2f"),#vstring 
                                          tauids = cms.VPSet(
                                              process.oldHPSLoose,
                                              process.oldHPSMedium,
                                              process.oldHPSTight,
                                              process.combinedTancOnePercent,
                                              process.combinedTancHalfPercent,
                                              process.combinedTancQuarterPercent,
                                              ),#vpset
                                          plotVariables = cms.VPSet(
                                              process.diTauVisMassFromJet
                                              ),#vspet
                                          muonCollections = cms.PSet(
                                              standardMuons = cms.InputTag('selectedPatMuonsForTauIdEffTrkCumulative'),
                                              standAloneMuons = cms.InputTag('patMuonsStandAlone'),
                                              ),#pset
                                          plotSys = cms.bool(True)#bool
                                          )

process.p = cms.Path(
    process.makePlotsChannel
    )
