import FWCore.ParameterSet.Config as cms

process = cms.Process("prodMuonIsoPATtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ##'file:/data1/veelken/CMSSW_4_2_x/skims/skimByHLTpath_IsoMu12_run167830_AOD_11_1_lSt.root'
        'file:/data1/veelken/CMSSW_4_2_x/skims/ppMuXpt15_1_1_Sl5.root'
    ),
    #skipEvents = cms.untracked.uint32(0),
    #eventsToProcess = cms.untracked.VEventRange('1:17:16907789')            
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#--------------------------------------------------------------------------------
# read AOD/RECO input files from castor directory

import TauAnalysis.Configuration.tools.castor as castor
inputFilePath = '/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/GoldenZmumu/simDYtoMuMu/'
inputFileNames = []
if inputFilePath.find('/castor/') != -1:
    inputFileNames = [ 'rfio:%s' % file_info['path'] for file_info in castor.nslsl(inputFilePath) ]
else:
    inputFileNames = [ 'file:%s' % os.path.join(inputFilePath, file_name) for file_name in os.listdir(inputFilePath) ]

sample = 'simDYtoMuMu'

#print "inputFileNames = %s" % inputFileNames

import re
inputFile_regex = \
  r"[a-zA-Z0-9_/:.]*goldenZmumuEvents_%s_AOD_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % sample
inputFile_matcher = re.compile(inputFile_regex)

inputFileNames_matched = []
for inputFileName in inputFileNames:
    if inputFile_matcher.match(inputFileName):
	inputFileNames_matched.append(inputFileName)

#print "inputFileNames_matched = %s" % inputFileNames_matched

setattr(process.source, "fileNames", cms.untracked.vstring(inputFileNames_matched))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define configuration parameter default values

isMC = True # use for MC
##isMC = False # use for Data
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "HLT" # use for Summer'11 MC
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
if isMC:
    process.GlobalTag.globaltag = cms.string('START42_V13::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V20::All')
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# produce PAT objects
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# configure PAT trigger matching    
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process, hltProcess = HLTprocessName, outputModule = '')

# select primary event vertex
process.load("TauAnalysis/RecoTools/recoVertexSelection_cff")

# produce 'pfNoPileUp' and 'pfPileUp' collections of particle-flow candidates
process.load("CommonTools/ParticleFlow/pfNoPileUp_cff")
process.pfPileUp.Enable = cms.bool(True)
process.pfPileUp.checkClosestZVertex = cms.bool(True)
process.pfPileUp.Vertices = cms.InputTag('selectedPrimaryVertexPosition')

process.load("CommonTools/ParticleFlow/pfParticleSelection_cff")

# compute muon IsoDeposits and add muon isolation sums to pat::Muon objects
process.load("RecoMuon/MuonIsolation/muonPFIsolation_cff")
import PhysicsTools.PatAlgos.tools.helpers as patutils
patutils.massSearchReplaceAnyInputTag(process.muonPFIsolationDepositsSequence, cms.InputTag('muons1stStep'), cms.InputTag('muons'))
process.patMuons.isoDeposits = cms.PSet(
    # CV: strings for IsoDeposits defined in PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc
    pfChargedHadrons = cms.InputTag("muPFIsoDepositCharged"),
    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutral"),
    pfPhotons = cms.InputTag("muPFIsoDepositGamma"),
    user = cms.VInputTag(
        cms.InputTag("muPFIsoDepositChargedAll"),
        cms.InputTag("muPFIsoDepositPU")
    )
)

process.patMuons.userIsolation = cms.PSet(
    # CV: strings for Isolation values defined in PhysicsTools/PatAlgos/src/MultiIsolator.cc
    pfChargedHadron = cms.PSet(
        deltaR = cms.double(0.4),
        src = process.patMuons.isoDeposits.pfChargedHadrons,
        vetos = process.muPFIsoValueCharged04.deposits[0].vetos,
        skipDefaultVeto = process.muPFIsoValueCharged04.deposits[0].skipDefaultVeto
    ),
    pfNeutralHadron = cms.PSet(
        deltaR = cms.double(0.4),
        src = process.patMuons.isoDeposits.pfNeutralHadrons,
        vetos = process.muPFIsoValueNeutral04.deposits[0].vetos,
        skipDefaultVeto = process.muPFIsoValueNeutral04.deposits[0].skipDefaultVeto
    ),
    pfGamma = cms.PSet(
        deltaR = cms.double(0.4),
        src = process.patMuons.isoDeposits.pfPhotons,
        vetos = process.muPFIsoValueGamma04.deposits[0].vetos,
        skipDefaultVeto = process.muPFIsoValueGamma04.deposits[0].skipDefaultVeto
    ),
    user = cms.VPSet(
        cms.PSet(
            deltaR = cms.double(0.4),
            src = process.patMuons.isoDeposits.user[0],
            vetos = process.muPFIsoValueChargedAll04.deposits[0].vetos,
            skipDefaultVeto = process.muPFIsoValueChargedAll04.deposits[0].skipDefaultVeto
        ),
        cms.PSet(
            deltaR = cms.double(0.4),
            src = process.patMuons.isoDeposits.user[1],
            vetos = process.muPFIsoValuePU04.deposits[0].vetos,
            skipDefaultVeto = process.muPFIsoValuePU04.deposits[0].skipDefaultVeto
        )
    )
)

# "clean" tau-jet candidate collection (i.e. remove tau-jets overlapping with muons)
# and require Pt > 15 GeV && |eta| < 2.5
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.cleanPatTaus.preselection = cms.string('')
process.cleanPatTaus.checkOverlaps = cms.PSet(
    muons = cms.PSet(
       src                 = cms.InputTag("patMuons"),
       algorithm           = cms.string("byDeltaR"),
       preselection        = cms.string("isGlobalMuon"),
       deltaR              = cms.double(0.5),
       checkRecoComponents = cms.bool(False),
       pairCut             = cms.string(""),
       requireNoOverlaps   = cms.bool(True)
    )
)
process.cleanPatTaus.finalCut = cms.string('pfJetRef().pt() > 15.0 & abs(pfJetRef().eta()) < 2.5')

# produce pat::MEt objects for "raw" particle-flow MEt
process.patPFMETs = process.patMETs.clone(
    metSource = cms.InputTag('pfMet'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(True)
)

from PhysicsTools.PatAlgos.tools.coreTools import *
if not isMC:
    # remove MC matching from standard PAT sequences
    removeMCMatching(process, ["All"], outputInProcess = False)
    process.patPFMETs.addGenMET = cms.bool(False)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define skimming criteria

hltIsoMu12paths = [
    'HLT_IsoMu12_v1', 
    'HLT_IsoMu12_v2',
    'HLT_IsoMu12_v3',
    'HLT_IsoMu12_v4',
    'HLT_IsoMu12_v5',
    'HLT_IsoMu12_v6',
    'HLT_IsoMu12_v7',
    'HLT_IsoMu12_v8',
    'HLT_IsoMu12_v9'
]

hltIsoMu15paths = [
    'HLT_IsoMu15_v1',
    'HLT_IsoMu15_v2',
    'HLT_IsoMu15_v3',
    'HLT_IsoMu15_v4',
    'HLT_IsoMu15_v5',
    'HLT_IsoMu15_v6',
    'HLT_IsoMu15_v7',
    'HLT_IsoMu15_v8',                             
    'HLT_IsoMu15_v9',
    'HLT_IsoMu15_v10',
    'HLT_IsoMu15_v11',
    'HLT_IsoMu15_v12',
    'HLT_IsoMu15_v13',
    'HLT_IsoMu15_v14'
]

hltIsoMu17paths = [
    'HLT_IsoMu17_v1',
    'HLT_IsoMu17_v2',
    'HLT_IsoMu17_v3',
    'HLT_IsoMu17_v4',
    'HLT_IsoMu17_v5',
    'HLT_IsoMu17_v6',
    'HLT_IsoMu17_v7',
    'HLT_IsoMu17_v8',                             
    'HLT_IsoMu17_v9',
    'HLT_IsoMu17_v10',
    'HLT_IsoMu17_v11',
    'HLT_IsoMu17_v12',
    'HLT_IsoMu17_v13',
    'HLT_IsoMu17_v14'
]

hltIsoMuPaths = []
hltIsoMuPaths.extend(hltIsoMu12paths)
hltIsoMuPaths.extend(hltIsoMu15paths)
hltIsoMuPaths.extend(hltIsoMu17paths)

process.load('TauAnalysis.TauIdEfficiency.filterQCDmuEnriched_cfi')
process.hltMu = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('hltMu'),             
        pluginType = cms.string('TriggerResultEventSelector'),
        src = cms.InputTag('TriggerResults::HLT'),
        triggerPaths = cms.vstring(hltIsoMuPaths)
    )
)

# require muons to satisfy Pt > 20 GeV, |eta| < 2.1 and to pass VBTF muon id. criteria
##process.selectedPatMuonsTight = cms.EDFilter("PATMuonSelector",
##    src = cms.InputTag('patMuons'),
##    cut = cms.string(
##        "pt > 15. & abs(eta) < 2.1 & isTrackerMuon & isGlobalMuon & numberOfMatches >= 2 & globalTrack.isNonnull &" \
##       + "globalTrack.hitPattern.numberOfValidMuonHits >= 1 & globalTrack.hitPattern.numberOfValidPixelHits >= 1 &"
##       + "globalTrack.hitPattern.numberOfValidTrackerHits >= 10 & globalTrack.normalizedChi2 < 10 & (globalTrack.ptError/globalTrack.pt) < 0.1 &" \
##       + "abs(userFloat('dxyWrtPV')) < 0.045 & abs(userFloat('dzWrtPV')) < 0.2"
##    ),
##    filter = cms.bool(True)
##)
process.preselectedPatMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag('patMuons'),
    cut = cms.string("pt > 15. & abs(eta) < 2.1"),
    filter = cms.bool(True)
)
process.selectedPatMuonsTight = cms.EDFilter("PATMuonIdSelector",
    src = cms.InputTag('preselectedPatMuons'),
    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    beamSpotSource = cms.InputTag('offlineBeamSpot'),                                         
    filter = cms.bool(True)                                         
)

# produce second collection of loosely selected pat::Muons, 
# used for vetoing di-muon events on analysis level
process.selectedPatMuonsLoose = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag('patMuons'),
    cut = cms.string("isGlobalMuon | isStandAloneMuon"),
    filter = cms.bool(True)
)

# embedd into pat::Muon objects information whether or not reconstructed muons are matched to trigger primitives
process.selectedPatMuonsTightTriggerMatched = cms.EDProducer("PATMuonTriggerEmbedder",
    src = cms.InputTag("selectedPatMuonsTight"),
    matched = cms.InputTag("patTrigger"),
    triggerPaths = cms.vstring(
        'HLT_IsoMu12_v*',
        'HLT_IsoMu15_v*',
        'HLT_IsoMu17_v*'
    ),
    dRmax = cms.double(0.5)
)

process.muTauJetPairs = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsTightTriggerMatched'),
    srcLeg2 = cms.InputTag('cleanPatTaus'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag('patPFMETs'),
    srcPrimaryVertex = cms.InputTag("selectedPrimaryVertexPosition"),
    srcBeamSpot = cms.InputTag("offlineBeamSpot"),
    recoMode = cms.string(""),
    doSVreco = cms.bool(False),
    nSVfit = cms.PSet(),
    scaleFuncImprovedCollinearApprox = cms.string('1'),
    doPFMEtSign = cms.bool(False),
    verbosity = cms.untracked.int32(0)
)

if isMC:
    process.dataQualityFilters.remove(process.hltPhysicsDeclared)
    process.dataQualityFilters.remove(process.dcsstatus)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# update InputTags for HLT trigger result object
# in case running on reprocessed Monte Carlo samples

if HLTprocessName != "HLT":
    process.hltMu.selector.src = cms.InputTag('TriggerResults::' + HLTprocessName)
    process.patTrigger.processName = HLTprocessName
    process.patTriggerEvent.processName = HLTprocessName
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# produce MC-to-Data reweight factors


if isMC:
    process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
    process.load("TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi")
    process.load("TauAnalysis.RecoTools.vertexMultiplicityVsRhoPFNeutralReweight_cfi")
    process.mcToDataReweightFactors = cms.Sequence(
        process.selectedPrimaryVertexQuality
       + process.selectedPrimaryVertexPosition
       + process.vertexMultiplicityReweight 
       + process.produceVertexMultiplicityVsRhoPFNeutralReweights
    )
else:
    process.mcToDataReweightFactors = cms.Sequence()
#--------------------------------------------------------------------------------

#-------------------------------------------------------------------------------- 
# add event counter for Mauro's "self baby-sitting" technology
process.totalEventsProcessed = cms.EDProducer("EventCountProducer")
process.eventCounterPath = cms.Path(process.totalEventsProcessed)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save PAT-tuple
#
process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep EventAux_*_*_*',
            'keep LumiSummary_*_*_*',                       
            'keep edmMergeableCounter_*_*_*',                                      
            'keep *_offlinePrimaryVertices_*_*',
            'keep *_offlinePrimaryVerticesWithBS_*_*',
            'keep *_selectedPrimaryVertexHighestPtTrackSum_*_*',       
            'keep *_selectedPatMuonsTightTriggerMatched_*_*',                     
            'keep *_selectedPatMuonsLoose_*_*', 
	    'keep *_cleanPatTaus_*_*',
            'keep *_muTauJetPairs_*_*',
            'keep *_patPFMETs_*_*',
            'keep *' # CV: only for testing
        )
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
            'muonIsoSkimPath'
        )
    ),
    fileName = cms.untracked.string("/data1/veelken/tmp/muonIsoStudy/muonIsolationPATtuple.root")
)

from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.patTupleOutputModule.outputCommands.extend(patTriggerEventContent)

if isMC:
    process.patTupleOutputModule.outputCommands.extend(
        cms.untracked.vstring(
            'keep *_addPileupInfo_*_*',
            'keep *_vertexMultiplicityReweight_*_*',
            'keep *_vertexMultiplicityVsRhoPFNeutralReweight_*_*'
        )
    )
#--------------------------------------------------------------------------------

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hltMu 
   + process.selectPrimaryVertex
   + process.pfNoPileUpSequence
   + process.pfParticleSelectionSequence
   + process.muonPFIsolationSequence
   + process.PFTau
   + process.patDefaultSequence 
   + process.patPFMETs
   + process.preselectedPatMuons
   + process.selectedPatMuonsTight
   + process.selectedPatMuonsLoose    
   + process.selectedPatMuonsTightPFIso04
   + process.selectedPatMuonsTightTriggerMatched
   + process.muTauJetPairs
   + process.mcToDataReweightFactors
)

process.muonIsoSkimPath = cms.Path(
    process.hltMu
   + process.dataQualityFilters
   + process.preselectedPatMuons
   + process.selectedPatMuonsTight
   + process.selectedPatMuonsLoose
)

process.o = cms.EndPath(process.patTupleOutputModule)

process.schedule = cms.Schedule(
    process.eventCounterPath,
    process.p,
    process.muonIsoSkimPath,
    process.o
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

processDumpFile = open('produceMuonIsolationPATtuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
