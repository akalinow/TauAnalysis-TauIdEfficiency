import FWCore.ParameterSet.Config as cms

process = cms.Process("runPATTauIsolationAnalyzer")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START38_V12::All')

process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patDefaultSequence.remove(process.patJets)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.tauPFIsolationAnalyzerBeforeTauId = cms.EDAnalyzer("TauPFIsolationAnalyzer",
    pluginName = cms.string('tauPFIsolationAnalyzerBeforeTauId'),

    dqmDirectory_store = cms.string('TauPFIsolationQuantities/beforeTauId'),

    leptonSource = cms.InputTag("patTaus"),
    ptMin = cms.double(20.),
    ptMax = cms.double(100.),
    etaMin = cms.double(-2.1),
    etaMax = cms.double(+2.1),

    dRisoCone = cms.double(0.5),

    pfCandidateSource = cms.InputTag("particleFlow"),

    genLeptonMatch = cms.PSet(
        genParticleSource = cms.InputTag("genParticles"),
        dRmatch = cms.double(0.5),
        genTauJetSource = cms.InputTag("tauGenJets"),                                              
        tauDecayModes = cms.vstring(
##--- define names of tau lepton decay modes
##   ( using names defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc )            
            "oneProng0Pi0",
            "oneProng1Pi0",
            "oneProng2Pi0",
            "oneProngOther",
            "threeProng0Pi0",
            "threeProng1Pi0",
            "threeProngOther",
            "rare"
        )
    ),

    normalization = cms.string("leptons")                                                       
)

process.muonPFIsolationAnalyzer = cms.EDAnalyzer("MuonPFIsolationAnalyzer",
    pluginName = cms.string('muonPFIsolationHistManager'),

    dqmDirectory_store = cms.string('MuonPFIsolationQuantities'),

    leptonSource = cms.InputTag("patMuons"),
    ptMin = cms.double(15.),
    ptMax = cms.double(100.),
    etaMin = cms.double(-2.1),
    etaMax = cms.double(+2.1),

    dRisoCone = cms.double(0.5),

    pfCandidateSource = cms.InputTag("particleFlow"),                                            

    genLeptonMatch = cms.PSet(
        genParticleSource = cms.InputTag("genParticles"),
        dRmatch = cms.double(0.5)
    ),

    normalization = cms.string("leptons")   
)

process.electronPFIsolationAnalyzer = cms.EDAnalyzer("ElectronPFIsolationAnalyzer",
    pluginName = cms.string('electronPFIsolationHistManager'),

    dqmDirectory_store = cms.string('ElectronPFIsolationQuantities'),

    leptonSource = cms.InputTag("patElectrons"),
    ptMin = cms.double(15.),
    ptMax = cms.double(100.),
    etaMin = cms.double(-2.1),
    etaMax = cms.double(+2.1),

    dRisoCone = cms.double(0.5),

    pfCandidateSource = cms.InputTag("particleFlow"),                                                     

    genLeptonMatch = cms.PSet(
        genParticleSource = cms.InputTag("genParticles"),
        dRmatch = cms.double(0.5)
    ),

    normalization = cms.string("leptons")   
)

process.selectedPatTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("patTaus"),                                       
    cut = cms.string('tauID("leadingTrackFinding") > 0.5 & tauID("leadingTrackPtCut") > 0.5 & tauID("byTaNCfrQuarterPercent") > 0.5'),
    filter = cms.bool(False)
)

process.tauPFIsolationAnalyzerAfterTauId = process.tauPFIsolationAnalyzerBeforeTauId.clone(
    pluginName = cms.string('tauPFIsolationAnalyzerAfterTauId'),

    dqmDirectory_store = cms.string('TauPFIsolationQuantities/afterTauId'),

    leptonSource = cms.InputTag("selectedPatTaus")
)

process.dqmSimpleFileSaver = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('patTauIsolationAnalyzer.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data1/veelken/CMSSW_3_6_x/skims/Ztautau_1_1_sXK.root'
    )
)

process.p = cms.Path(
    process.genParticlesForJets + process.ak5GenJets
   + process.patDefaultSequence
   + process.tauPFIsolationAnalyzerBeforeTauId + process.muonPFIsolationAnalyzer + process.electronPFIsolationAnalyzer
   + process.selectedPatTaus + process.tauPFIsolationAnalyzerAfterTauId
   + process.dqmSimpleFileSaver
)

# print-out all python configuration parameter information
#print process.dumpPython()
