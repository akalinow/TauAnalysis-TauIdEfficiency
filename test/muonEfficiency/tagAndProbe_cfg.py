import FWCore.ParameterSet.Config as cms
import glob
import sys

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Load our information
import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi \
        as samples
hlt_dict = samples.RECO_SAMPLES['data_Mu_Run2010B_Nov4ReReco']['hlt_paths']

sample_map = {
    'data' : {
        'files' : ['file:' + file for file in glob.glob(
            '/data1/friis/ZmumuSkim/dataB*.root')],
        'trigger' : 'HLT',
        'paths' : ['*'],
        'add_mc' : False,
    },
    'mc' : {
        'files' : ['file:' + file for file in glob.glob(
            '/data1/friis/ZmumuSkim/mc0*.root')],
        'trigger' : 'REDIGI38X',
        'paths' : ['HLT_Mu9'],
        'add_mc' : True,
    },
    'mcpu' : {
        'files' : ['file:' + file for file in glob.glob(
            '/data1/friis/ZmumuSkim/mcpu0*.root')],
        'trigger' : 'REDIGI38XPU',
        'paths' : ['HLT_Mu9'],
        'add_mc' : True,
    },
}

source = None
if len(sys.argv) > 2:
    source = sys.argv[2]
else:
    source = sys.argv[1]


print "Source is:", source
sample_info = sample_map[source]

print sample_info['files']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        sample_info['files']
    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V8::All')

# HLT selection
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.fastFilter = hltHighLevel.clone(
    HLTPaths = sample_info['paths'],
)

# The main path - first require the HLT
process.tagAndProbe = cms.Path(
    process.fastFilter
)

## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====
from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    muons     = cms.InputTag("muons"),
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 5 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 5"),
    tracksCut    = cms.string("pt > 5"),
)
process.tagAndProbe += process.mergedMuons


# Build and trigger match pat::Muons
## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
## with some customization
process.muonMatchHLTL2.maxDeltaR = 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
# Add to process
process.tagAndProbe += process.patMuonsWithTriggerSequence

# Add our extra triggers to the collection
process.triggerMatcher = process.muonTriggerMatchHLT.clone(
    pathNames = cms.vstring(hlt_dict.keys()),
)
process.patTriggerMatchers1Mu += process.triggerMatcher
process.patMuonsWithTrigger.matches += [cms.InputTag("triggerMatcher")]

# Add run-ranged based trigger information to pat::Muons
process.patMuonsWithRunRangeTrigger = cms.EDProducer(
    "PATMuonTriggerRangeInfoProducer",
    src = cms.InputTag("patMuonsWithTrigger"),
    userIntName = cms.string("triggerRange"),
)
# Parse trigger matching info
if not sample_info['add_mc']:
    # Only do the run ranges for data
    process.patMuonsWithRunRangeTrigger.config = cms.VPSet()
    for bit, range in hlt_dict.iteritems():
        to_add = cms.PSet(
            hltAcceptPath = cms.string(bit),
            runrange = cms.EventRange(range)
        )
        process.patMuonsWithRunRangeTrigger.config.append(to_add)
else:
    # Otherwise just match HLT_Mu9
    process.patMuonsWithRunRangeTrigger.hltAcceptPaths = cms.vstring('HLT_Mu9')

process.patMuonsWithRunRangeSequence = cms.Sequence(
    process.patMuonsWithTriggerSequence
    *process.patMuonsWithRunRangeTrigger
)

process.tagAndProbe += process.patMuonsWithRunRangeSequence
print process.patMuonsWithRunRangeSequence
print process.tagAndProbe

# Change the input to the pat muons
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import \
        changeRecoMuonInput
changeRecoMuonInput(process, "mergedMuons")


from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import \
        patMuonPFIsolationSelector

# Embed information about the PF isolation variables about the tau
process.patMuonsWithEmbeddedIso = cms.EDProducer(
    "PATMuonPFIsolationEmbedder",
    patMuonPFIsolationSelector,
    src = cms.InputTag("patMuonsWithRunRangeTrigger"),
    userFloatName = cms.string('absIso'),
)

process.patMuonsWithEmbeddedIso.chargedHadronIso.ptMin = 1.0
process.patMuonsWithEmbeddedIso.neutralHadronIso.ptMin = 2000.
process.patMuonsWithEmbeddedIso.photonIso.ptMin = 1.5
# This dont' matter here, but leave it for reference
process.patMuonsWithEmbeddedIso.sumPtMax = 1.0
process.patMuonsWithEmbeddedIso.sumPtMethod = "absolute"
process.tagAndProbe += process.patMuonsWithEmbeddedIso

# Do another one, using PF no pileup
# Actually we already use PF no pileup! This is redundant.
process.patMuonsWithEmbeddedIsoPFNoPileUp = \
        process.patMuonsWithEmbeddedIso.clone(
            src = cms.InputTag("patMuonsWithEmbeddedIso"),
            userFloatName = cms.string('absIsopfNoPileup'),

)
process.tagAndProbe += process.patMuonsWithEmbeddedIsoPFNoPileUp

import MuonAnalysis.TagAndProbe.common_variables_cff as common
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

# Muons that pass the pt and ID requirement
process.tagMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("patMuonsWithEmbeddedIsoPFNoPileUp"),
    cut = cms.string(
        "pt > 15 && "
        "userInt('triggerRange') != 0 && "
        + common.MuonIDFlags.VBTF.value() )
)
process.tagAndProbe += process.tagMuons

process.tagAndProbe += process.nverticesModule

process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

# Build probe muons (no real cut now)
process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithEmbeddedIsoPFNoPileUp"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)
process.tagAndProbe += process.probeMuons

# Make tag probe pairs
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('40 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.tagAndProbe += process.tpPairs


# Update trigger process information
process.fastFilter.TriggerResultsTag = cms.InputTag(
    "TriggerResults","", sample_info['trigger'])
process.patTrigger.processName = sample_info['trigger']

# Build tag probe tree
process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = common.AllVariables,
    flags = cms.PSet(
        common.TrackQualityFlags,
        common.MuonIDFlags,
        common.HighPtTriggerFlags,
        ## Isolation
        # FIX ME, add our Iso here?
        Isol    = cms.string(
            "(isolationR03.emEt + isolationR03.hadEt + isolationR03.sumPt)/pt < 0.15"),
        IsolTk3 = cms.string("isolationR03.sumPt < 3"),
        CustomTrigger = cms.string("userInt('triggerRange') > 0"),
        AbsIso = cms.string("userFloat('absIso') < 1.0"),
        AbsIsoPFnoPileUp = cms.string("userFloat('absIsopfNoPileup') < 1.0"),
    ),
    tagVariables = cms.PSet(
        nVertices = cms.InputTag("nverticesModule"),
    ),
    tagFlags = cms.PSet(),
    isMC           = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

if sample_info['add_mc']:
    # Build matching if this is MC
    process.tagMuonsMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        src = cms.InputTag("tagMuons"),
        matched = cms.InputTag("genParticles"),
        pdgId = cms.vint32(13),
        distMin = cms.double(0.3),
    )
    process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(
        src = "probeMuons")
    process.tagAndProbe += process.tagMuonsMCMatch
    process.tagAndProbe += process.probeMuonsMCMatch

    # Update the tree
    process.tpTree.isMC = True
    process.tpTree.tagMatches = cms.InputTag("tagMuonsMCMatch")
    process.tpTree.probeMatches = cms.InputTag("probeMuonsMCMatch")
    process.tpTree.motherPdgId = cms.vint32(22, 23)
    process.tpTree.makeMCUnbiasTree       = cms.bool(True)
    process.tpTree.checkMotherInUnbiasEff = cms.bool(True)
    process.tpTree.allProbes              = cms.InputTag("probeMuons")


process.tagAndProbe += process.tpPairs
process.tagAndProbe += process.tpTree

# FIXME add MC info to STA tree

##############################################################
# Making STA collection
# ############################################################

#Now make another collection for standalone muons
# Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)

# Match to trigger
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet, \
        massSearchReplaceAnyInputTag
# Note we use our run range sequence
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(
    process, process.patMuonsWithTriggerSequence, "Sta")
process.muonMatchHLTL2Sta.maxDeltaR = 0.5
process.muonMatchHLTL3Sta.maxDeltaR = 0.5

massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta,
                             "mergedMuons", "muonsSta")

## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("outerTrack.isNonnull"), # no real cut now
)
process.tpPairsSta = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsSta@-")

## Now I have to define the passing probes for tracking
process.tkTracks = cms.EDProducer(
    "ConcreteChargedCandidateProducer",
    src = cms.InputTag("generalTracks"),
    particleType = cms.string("mu+"),
)
process.staToTkMatch = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("probeMuonsSta"),
    matched = cms.InputTag("tkTracks"),
    algorithm = cms.string("byDirectComparison"),
    srcTrack     = cms.string("muon"),    srcState = cms.string("atVertex"),
    matchedTrack = cms.string("tracker"), matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(0.3),
    maxDeltaPtRel    = cms.double(2),   # |pt(sta) - pt(tk)|/pt(tk)
    maxDeltaLocalPos = cms.double(100),
    sortBy           = cms.string("deltaR"),
)
process.staPassingTk = cms.EDProducer("MatchedCandidateSelector",
    src   = cms.InputTag("probeMuonsSta"),
    match = cms.InputTag("staToTkMatch"),
)

process.tpTreeSta = process.tpTree.clone(
    tagProbePairs = "tpPairsSta",
    variables = cms.PSet(
        common.KinematicVariables,
        common.TriggerVariables,
        ## extra standalone muon quality variables
        outerHits      = cms.string("outerTrack.hitPattern.numberOfHits"),
        outerValidHits = cms.string("outerTrack.numberOfValidHits"),
        outerStationsAny   = cms.string("outerTrack.hitPattern.muonStationsWithAnyHits"),
        outerStationsValid = cms.string("outerTrack.hitPattern.muonStationsWithValidHits"),
        ## track matching variables
        tk_deltaR     = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaPtRel = cms.InputTag("staToTkMatch","deltaPtRel"),
        tk_deltaEta   = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaPhi   = cms.InputTag("staToTkMatch","deltaPhi"),
    ),
    flags = cms.PSet(
        common.HighPtTriggerFlags,
        hasTrack = cms.InputTag("staPassingTk"),
        L1DoubleMuOpen       = common.LowPtTriggerFlagsPhysics.L1DoubleMuOpen,
        L1DoubleMuOpen_Tight = common.LowPtTriggerFlagsPhysics.L1DoubleMuOpen_Tight,
        L2DoubleMu0          = common.LowPtTriggerFlagsPhysics.L2DoubleMu0,
    )
)

process.tpTreeSta.variables.l1pt = process.tpTreeSta.variables.l1pt.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1q  = process.tpTreeSta.variables.l1q.value( ).replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.variables.l1dr = process.tpTreeSta.variables.l1dr.value().replace("muonL1Info","muonL1InfoSta")
process.tpTreeSta.tagFlags = process.tpTreeSta.flags.clone(hasTrack = cms.string(""))

process.tnpSimpleSequenceSta = cms.Sequence(
    ( process.tkTracks * process.staToTkMatch * process.staPassingTk ) +
    process.tpPairsSta      +
    process.tpTreeSta
)

process.tagAndProbeSta = cms.Path(
    process.fastFilter                     +
    process.muonsSta                       +
    process.patMuonsWithTriggerSequenceSta
)

process.tagAndProbeSta += process.tagMuons
process.tagAndProbeSta += process.nverticesModule
process.tagAndProbeSta += process.probeMuonsSta

if sample_info['add_mc']:
    # Build matching if this is MC
    process.probeMuonsMCMatchSta = process.tagMuonsMCMatch.clone(
        src = "probeMuonsSta")
    process.tagAndProbeSta += process.tagMuonsMCMatch
    process.tagAndProbeSta += process.probeMuonsMCMatchSta

    # Update the tree
    process.tpTreeSta.isMC = True
    process.tpTreeSta.tagMatches = cms.InputTag("tagMuonsMCMatch")
    process.tpTreeSta.probeMatches = cms.InputTag("probeMuonsMCMatchSta")
    process.tpTreeSta.motherPdgId = cms.vint32(22, 23)
    process.tpTreeSta.makeMCUnbiasTree       = cms.bool(True)
    process.tpTreeSta.checkMotherInUnbiasEff = cms.bool(True)
    process.tpTreeSta.allProbes              = cms.InputTag("probeMuonsSta")

process.tagAndProbeSta += process.tnpSimpleSequenceSta

process.TFileService = cms.Service("TFileService", fileName = cms.string(
    "/data1/friis/ZmumuEff/tagAndProbe_%s.root" % source))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
