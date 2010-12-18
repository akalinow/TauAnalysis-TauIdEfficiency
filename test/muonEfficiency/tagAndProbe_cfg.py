import FWCore.ParameterSet.Config as cms
import glob
import sys
import copy
import itertools

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Load our information
import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi \
        as samples
hlt_dict = samples.RECO_SAMPLES['data_Mu_Run2010B_Nov4ReReco']['hlt_paths']

def chunk(inputs, chunks):
    output = []
    for input in inputs:
        output.append(input)
        if len(output) == input/chunks:
            yield output
            output = []

sample_map = {
    'data' : {
        'files' : ['file:' + file for file in itertools.chain(
            #glob.glob('/data1/friis/ZmumuSkim/dataB*.root'),
            #glob.glob('/data1/friis/ZmumuSkim/data[.0-9]*.root')
            glob.glob('/data1/friis/ZmumuSkim/data*.root')
        )],
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
            '/data1/friis/ZmumuSkim/mcpu*.root')],
        'trigger' : 'REDIGI38XPU',
        'paths' : ['HLT_Mu9'],
        'add_mc' : True,
    },
}


jobs_file = open('jobs.txt', 'w')
for index, file in enumerate(sample_map['data']['files']):
    jobs_file.write('data %i\n' % index)

for index, file in enumerate(sample_map['mcpu']['files']):
    jobs_file.write('mcpu %i\n' % index)


print "Run:"
print "cat jobs.txt | xargs -n 2 -P 5 cmsRun tagAndProbe_cfg.py"
print "to submit"

source = None
jobindex = None
if len(sys.argv) > 3:
    source = sys.argv[2]
    jobindex = sys.argv[3]
else:
    source = sys.argv[1]
    jobindex = sys.argv[2]
jobindex = int(jobindex)

print "Source is:", source, jobindex
sample_info = sample_map[source]

print "File is:", sample_info['files'][jobindex]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        sample_info['files'][jobindex]
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
        # Skip cross triggers
        if bit.find('PFTau') != -1:
            print "Skipping xtrigger %s period %s" % (bit, range)
            continue
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

_CROSS_TRIGGER_PERIOD = '148822:MIN-149442:MAX'
# Period A is when HLT_Mu9 is unprescaled
_PERIOD_A = '132440:MIN-147116:MAX'
# Period B is when IsoMu9 and Mu11 are unprescaled
_PERIOD_B = '147196:MIN-148058:MAX'
# Period C is when IsoMu13 and Mu15 are unprescaled
_PERIOD_C = '148059:MIN-149442:MAX'

# Produce a flag indicating whether or not we are in the cross trigger period
process.patMuonsAddXTriggerPeriodInfo = cms.EDProducer(
    "PATMuonTriggerRangeInfoProducer",
    src = cms.InputTag("patMuonsWithRunRangeTrigger"),
    userIntName = cms.string("xTriggerPeriod"),
    config = cms.VPSet(
        cms.PSet(
            hltAcceptPath = cms.string('*'),
            runrange = cms.EventRange(_CROSS_TRIGGER_PERIOD)
        ),
    )
)

process.tagAndProbe += process.patMuonsAddXTriggerPeriodInfo

process.patMuonsAddPeriodAInfo = cms.EDProducer(
    "PATMuonTriggerRangeInfoProducer",
    src = cms.InputTag("patMuonsAddXTriggerPeriodInfo"),
    userIntName = cms.string("periodA"),
    config = cms.VPSet(
        cms.PSet(
            hltAcceptPath = cms.string('*'),
            runrange = cms.EventRange(_PERIOD_A)
        ),
    )
)
process.tagAndProbe += process.patMuonsAddPeriodAInfo

process.patMuonsAddPeriodBInfo = cms.EDProducer(
    "PATMuonTriggerRangeInfoProducer",
    src = cms.InputTag("patMuonsAddPeriodAInfo"),
    userIntName = cms.string("periodB"),
    config = cms.VPSet(
        cms.PSet(
            hltAcceptPath = cms.string('*'),
            runrange = cms.EventRange(_PERIOD_B)
        ),
    )
)
process.tagAndProbe += process.patMuonsAddPeriodBInfo

process.patMuonsAddPeriodCInfo = cms.EDProducer(
    "PATMuonTriggerRangeInfoProducer",
    src = cms.InputTag("patMuonsAddPeriodBInfo"),
    userIntName = cms.string("periodC"),
    config = cms.VPSet(
        cms.PSet(
            hltAcceptPath = cms.string('*'),
            runrange = cms.EventRange(_PERIOD_C)
        ),
    )
)
process.tagAndProbe += process.patMuonsAddPeriodCInfo

# Change the input to the pat muons
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import \
        changeRecoMuonInput
changeRecoMuonInput(process, "mergedMuons")


from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import \
        patMuonPFIsolationSelector

# Embed information about the PF isolation variables about the tau
process.patMuonsWithEmbeddedIsoForRel = cms.EDProducer(
    "PATMuonPFIsolationEmbedder",
    copy.deepcopy(patMuonPFIsolationSelector),
    src = cms.InputTag("patMuonsAddPeriodCInfo"),
    userFloatName = cms.string('isoActivityForRel'),
)
process.tagAndProbe += process.patMuonsWithEmbeddedIsoForRel

# Raise thresholds
process.patMuonsWithEmbeddedIso = cms.EDProducer(
    "PATMuonPFIsolationEmbedder",
    patMuonPFIsolationSelector,
    src = cms.InputTag("patMuonsWithEmbeddedIsoForRel"),
    userFloatName = cms.string('isoActivityForAbs'),
)

process.patMuonsWithEmbeddedIso.chargedHadronIso.ptMin = 1.0
process.patMuonsWithEmbeddedIso.neutralHadronIso.ptMin = 2000.
process.patMuonsWithEmbeddedIso.photonIso.ptMin = 1.5
# This dont' matter here, but leave it for reference
process.patMuonsWithEmbeddedIso.sumPtMax = 1.0
process.patMuonsWithEmbeddedIso.sumPtMethod = "absolute"
process.tagAndProbe += process.patMuonsWithEmbeddedIso

import MuonAnalysis.TagAndProbe.common_variables_cff as common
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

# Muons that pass the pt and ID requirement
process.tagMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("patMuonsWithEmbeddedIso"),
    cut = cms.string(
        "pt > 15 && "
        "userInt('triggerRange') != 0 && "
        + common.MuonIDFlags.VBTF.value() + " && " +
        "%s < 0.1" % common.IsolationVariables.combRelIso.value()
    )
)
process.tagAndProbe += process.tagMuons

process.atLeastOneTag = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tagMuons"),
    minNumber = cms.uint32(1),
)
process.tagAndProbe += process.atLeastOneTag

# Load our vertex selection
process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
process.tagAndProbe += process.selectPrimaryVertex

process.load("TauAnalysis.RecoTools.vertexMultiplicityReweight_cfi")
process.tagAndProbe += process.selectedPrimaryVerticesTrackPtSumGt10

# Use only vertices w/ pt > 10 in our vertex count
process.nverticesModule.objects = cms.InputTag(
    "selectedPrimaryVerticesTrackPtSumGt10")
process.tagAndProbe += process.nverticesModule

# Build probe muons (no real cut now)
process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithEmbeddedIso"),
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

def trigger_match(*paths):
    return " | ".join(
        "!triggerObjectMatchesByPath('%s').empty()" % path for path in paths)

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
        innerTrack = cms.string('innerTrack.isNonnull'),
        outerTrack = cms.string('outerTrack.isNonnull'),
        AbsIso = cms.string("userFloat('isoActivityForAbs') < 1.0"),
        LooseAbsIso = cms.string("userFloat('isoActivityForAbs') < 2.5"),
        RelIso = cms.string("userFloat('isoActivityForRel') < (0.10*pt)"),
        LooseRelIso = cms.string("userFloat('isoActivityForRel') < (0.15*pt)"),
        PeriodA = cms.string("userInt('periodA') > 0"),
        PeriodB = cms.string("userInt('periodB') > 0"),
        PeriodC = cms.string("userInt('periodC') > 0"),
        PtThresh = cms.string("pt > 15"),
        EtaCut = cms.string("abs(eta) < 2.1"),
        HLTMu9 = cms.string(trigger_match('HLT_Mu9')),
        IsoMu9Mu11 = cms.string(trigger_match('HLT_IsoMu9', 'HLT_Mu11')),
        IsoMu13Mu15 = cms.string(
            trigger_match('HLT_IsoMu13_v3', 'HLT_IsoMu13_v4', 'HLT_Mu15_v1')),
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

print process.patMuonsWithTriggerSequenceSta

## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
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
        PtThresh = cms.string("pt > 15"),
        EtaCut = cms.string("abs(eta) < 2.1"),
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
process.tagAndProbeSta += process.atLeastOneTag
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

##############################################################
# Making Inner collection
# ############################################################

#Now make another collection for standalone muons
# Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsInner = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("inner"),
)

# Match to trigger
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet, \
        massSearchReplaceAnyInputTag
# Note we use our run range sequence
process.patMuonsWithTriggerSequenceInner = cloneProcessingSnippet(
    process, process.patMuonsWithTriggerSequence, "Inner")
process.muonMatchHLTL2Inner.maxDeltaR = 0.5
process.muonMatchHLTL3Inner.maxDeltaR = 0.5

massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceInner,
                             "mergedMuons", "muonsInner")

## Define probes and T&P pairs
process.probeMuonsInner = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerInner"),
    cut = cms.string("innerTrack.isNonnull"), # no real cut now
)
process.tpPairsInner = process.tpPairs.clone(decay = "tagMuons@+ probeMuonsInner@-")

## Now I have to define the passing probes for tracking
process.tkTracks = cms.EDProducer(
    "ConcreteChargedCandidateProducer",
    src = cms.InputTag("generalTracks"),
    particleType = cms.string("mu+"),
)

process.tpTreeInner = process.tpTree.clone(
    tagProbePairs = "tpPairsInner",
    variables = cms.PSet(
        common.KinematicVariables,
        common.TriggerVariables,
    ),
    flags = cms.PSet(
        common.HighPtTriggerFlags,
        common.MuonIDFlags,
        PtThresh = cms.string("pt > 15"),
        EtaCut = cms.string("abs(eta) < 2.1"),
    )
)

process.tnpSimpleSequenceInner = cms.Sequence(
    ( process.tkTracks * process.staToTkMatch * process.staPassingTk ) +
    process.tpPairsInner      +
    process.tpTreeInner
)

process.tagAndProbeInner = cms.Path(
    process.fastFilter                     +
    process.muonsInner                       +
    process.patMuonsWithTriggerSequenceInner
)

process.tagAndProbeInner += process.tagMuons
process.tagAndProbeInner += process.atLeastOneTag
process.tagAndProbeInner += process.nverticesModule
process.tagAndProbeInner += process.probeMuonsInner

if sample_info['add_mc']:
    # Build matching if this is MC
    process.probeMuonsMCMatchInner = process.tagMuonsMCMatch.clone(
        src = "probeMuonsInner")
    process.tagAndProbeInner += process.tagMuonsMCMatch
    process.tagAndProbeInner += process.probeMuonsMCMatchInner

    # Update the tree
    process.tpTreeInner.isMC = True
    process.tpTreeInner.tagMatches = cms.InputTag("tagMuonsMCMatch")
    process.tpTreeInner.probeMatches = cms.InputTag("probeMuonsMCMatchInner")
    process.tpTreeInner.motherPdgId = cms.vint32(22, 23)
    process.tpTreeInner.makeMCUnbiasTree       = cms.bool(True)
    process.tpTreeInner.checkMotherInUnbiasEff = cms.bool(True)
    process.tpTreeInner.allProbes              = cms.InputTag("probeMuonsInner")

process.tagAndProbeInner += process.tnpSimpleSequenceInner

process.TFileService = cms.Service("TFileService", fileName = cms.string(
    "/data2/friis/ZmumuEff/tagAndProbe_%s_%i.root" % (source, jobindex)))

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
