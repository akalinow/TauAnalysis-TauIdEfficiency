import FWCore.ParameterSet.Config as cms

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os
import re

inputFilePath = \
  '/data2/veelken/CMSSW_4_2_x/PATtuples/MuonIsolation/2011Oct10/v4/user/v/veelken/CMSSW_4_2_x/PATtuples/MuonIsolation/2011Oct10/v4/'

sampleToAnalyze = 'PPmuXptGt20Mu15'
jobId = '2011Oct10v3'

inputFile_regex = \
  r"muonIsolationPATtuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)

outputFilePath = '/data1/veelken/tmp/muonIsoStudy/v4/'

# check if inputFile is PAT-tuple and
# matches sampleToAnalyze, jobId
inputFileNames = []
files = os.listdir(inputFilePath)
for file in files:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(file):
        inputFileNames.append(os.path.join(inputFilePath, file))
#print "inputFileNames = %s" % inputFileNames 
##inputFileNames = [
##    inputFileNames[0] # for TESTING only !!
##]

# find name of associated "process"
process_matched = None
processes = recoSampleDefinitionsZtoMuTau_7TeV['MERGE_SAMPLES'].keys()
for process in processes:
    for sample in recoSampleDefinitionsZtoMuTau_7TeV['MERGE_SAMPLES'][process]['samples']:
        if sample == sampleToAnalyze:
            process_matched = process

if not process_matched:
    raise ValueError("No process associated to sample %s --> skipping !!" % sampleToAnalyze)

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFileNames),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string(os.path.join(outputFilePath, 'analyzeMuonIsoHistograms_%s_%s.root' % (sampleToAnalyze, jobId)))
)

process.muonIsolationAnalyzer = cms.PSet(
  
    directory = cms.string(''),

    triggerPaths = cms.vstring('HLT_IsoMu12', 'HLT_IsoMu15'),
    muonIsoThresholdsLoose = cms.vdouble(0.5, 1.0),
    muonIsoThresholdsTight = cms.vdouble(0.10),

    muonIsoProbExtractor = cms.PSet(
        inputFileName = cms.FileInPath('TauAnalysis/TauIdEfficiency/data_nocrab/train_kNNmuonIsolation_kNN.weights.xml'),
        parametrization = cms.VPSet(            
            cms.PSet(
                name = cms.string('logMuonPt'),
                expression = cms.string('log(pt)')
            ),
            cms.PSet(
                name = cms.string('absMuonEta'),
                expression = cms.string('abs(eta)')
            )
        ),
        selection = cms.string(
            '(userIsolation("pat::User1Iso")' + \
            ' + max(0., userIsolation("pat::PfNeutralHadronIso") + userIsolation("pat::PfGammaIso")' + \
            '          - 0.5*userIsolation("pat::User2Iso"))) > 0.20*pt'
        )
    ),

    srcMuonsTightId = cms.InputTag('selectedPatMuonsTightPFIso04'),
    srcMuonsLooseId = cms.InputTag('selectedPatMuonsLoose'),
    srcTauJetCandidates = cms.InputTag('cleanPatTausTriggerMatched'),
    srcMuTauPairs = cms.InputTag('muTauJetPairs'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    
    weights = cms.VInputTag('vertexMultiplicityReweight')
)
