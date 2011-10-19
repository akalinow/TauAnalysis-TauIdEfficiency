import FWCore.ParameterSet.Config as cms

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os
import re

inputFilePath = \
  '/data2/veelken/CMSSW_4_2_x/PATtuples/MuonIsolation/2011Oct10/user/v/veelken/CMSSW_4_2_x/PATtuples/MuonIsolation/2011Oct10/'

sampleToAnalyze = 'PPmuXptGt20Mu15'
jobId = '2011Oct10ii'

inputFile_regex = \
  r"muonIsolationPATtuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)

outputFilePath = '/data1/veelken/tmp/muonIsoStudy/'

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
##    
##    '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/muonIsolationPATtuple.root'
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

    srcMuonsTightId = cms.InputTag('selectedPatMuonsTightPFIso04'),
    srcMuonsLooseId = cms.InputTag('selectedPatMuonsLoose'),
    srcTauJetCandidates = cms.InputTag('cleanPatTausTriggerMatched'),
    srcMuTauPairs = cms.InputTag('muTauJetPairs'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    
    weights = cms.VInputTag('vertexMultiplicityReweight')
)
