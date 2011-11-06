import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
import TauAnalysis.Configuration.tools.castor as castor
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import \
  make_inputFileNames_vstring, make_tauIds_string, getStringRep_bool

import os
import re

#--------------------------------------------------------------------------------
#
# define auxiliary functions
#
def make_drawOptions_string(drawOptions, namesToPlot):
    retVal = ""
    for nameToPlot in namesToPlot:
        retVal += "        " + "cms.PSet(" + "\n"
        retVal += "        " + "    name = cms.string('%s'),\n" % nameToPlot
        if 'avTriggerPrescale' in drawOptions[nameToPlot]:
            retVal += "        " + "    avTriggerPrescale = cms.double(%f),\n" % drawOptions[nameToPlot]['avTriggerPrescale']
        retVal += "        " + "    legendEntry = cms.string('%s'),\n" % drawOptions[nameToPlot]['legendEntry']
        retVal += "        " + "    markerStyleData = cms.uint32(%u),\n" % drawOptions[nameToPlot]['markerStyleData']
        retVal += "        " + "    markerStyleSim = cms.uint32(%u),\n" % drawOptions[nameToPlot]['markerStyleSim']
        retVal += "        " + "    color = cms.uint32(%u)\n" % drawOptions[nameToPlot]['color']
        retVal += "        " + ")," + "\n"
    return retVal
#--------------------------------------------------------------------------------

def buildConfigFile_FWLiteTauFakeRateAnalyzer(sampleToAnalyze, evtSel, version, inputFilePath, tauIds, 
                                              tauJetCandSelection, srcTauJetCandidates, srcMET, hltPaths, 
                                              configFilePath, logFilePath, outputFilePath, recoSampleDefinitions):

    """Build cfg.py file to run FWLiteTauFakeRateAnalyzer macro to run on PAT-tuples,
       and fill histograms for passed/failed samples"""

    #print "inputFilePath = %s" % inputFilePath

    inputFileNames = None
    if inputFilePath.find('/castor/') != -1:
        inputFileNames = [ file_info['path'] for file_info in castor.nslsl(inputFilePath) ]
    else:
        inputFileNames = os.listdir(inputFilePath)
    #print "inputFileNames = %s" % inputFileNames

    # check if inputFile is PAT-tuple and
    # matches sampleToAnalyze, jobId
    inputFileNames_sample = []
    for inputFileName in inputFileNames:        
        if inputFileName.find("chunk") != -1 and \
           inputFileName.find("".join(['_', sampleToAnalyze, '_'])) != -1:
            # CV: assume that input file gets copied to local directory before FWLiteTauFakeRateAnalyzer macro gets started
            inputFileNames_sample.append(os.path.basename(inputFileName))

    #print(sampleToAnalyze)
    #print(inputFiles_sample)

    if len(inputFileNames_sample) == 0:
        print("Sample %s, evtSel = %s has no input files --> skipping !!" % (sampleToAnalyze, evtSel))
        return

    # find name of associated "process"
    process_matched = None
    processes = recoSampleDefinitions['MERGE_SAMPLES'].keys()
    for process in processes:
        for sample in recoSampleDefinitions['MERGE_SAMPLES'][process]['samples']:
            if sample == sampleToAnalyze:
                process_matched = process

    if not process_matched:
        print("No process associated to sample %s --> skipping !!" % sampleToAnalyze)
        return

    print("building config file(s) for sample %s, evtSel %s..." % (sampleToAnalyze, evtSel))

    processType = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['type']

    tauIds_string = make_tauIds_string(tauIds)

    hltPaths_string = make_inputFileNames_vstring(hltPaths[processType])

    configFileNames = []
    outputFileNames = []
    logFileNames    = []

    for inputFileName_sample in inputFileNames_sample:

        inputFileName_regex = r"[a-zA-Z0-9_./]*skim_(?P<sample>\w+?)_chunk_(?P<jobId>\d*)_(?P<hash>[a-zA-Z0-9]*).root"
        inputFileName_matcher = re.compile(inputFileName_regex)
        match = inputFileName_matcher.match(inputFileName_sample)
        if not match:
            raise ValueError("Failed to parse fileName = %s !!" % inputFileName_sample)
        jobId = match.group('jobId')

        outputFileName = 'analyzeTauFakeRateHistograms_%s_%s_%s_chunk_%s.root' % (evtSel, sampleToAnalyze, version, jobId)

        weights_string = ""
        if not recoSampleDefinitions['MERGE_SAMPLES'][process_matched]['type'] == 'Data':
            weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight", "'"])

        allEvents_DBS = -1
        xSection = 0.0
        if not recoSampleDefinitions['MERGE_SAMPLES'][process_matched]['type'] == 'Data':
            allEvents_DBS = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['events_processed']
            xSection = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['x_sec']
        intLumiData = recoSampleDefinitions['TARGET_LUMI']

        config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('%s'),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.tauFakeRateAnalyzer = cms.PSet(
    process = cms.string('%s'),
    type = cms.string('%s'),

    evtSel = cms.string('%s'),

    regions = cms.vstring(
        'P',
        'F',
        'A'
    ),
    
    tauIds = cms.VPSet(
%s
    ),
    
    srcTauJetCandidates = cms.InputTag('%s'),
    tauJetCandSelection = cms.vstring(
%s
    ),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(%s),
    
    srcMET = cms.InputTag('%s'),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(%s),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
""" % (inputFileName_sample, outputFileName,
       process_matched, processType, evtSel,
       tauIds_string, srcTauJetCandidates, tauJetCandSelection, hltPaths_string, srcMET, weights_string,
       allEvents_DBS, xSection, intLumiData)

        outputFileNames.append(outputFileName)

        configFileName = "analyzeTauFakeRatePATtuple_%s_%s_%s_cfg.py" % (evtSel, sampleToAnalyze, jobId)
        configFileName_full = os.path.join(configFilePath, configFileName)    
        configFile = open(configFileName_full, "w")
        configFile.write(config)
        configFile.close()
        configFileNames.append(configFileName)

        logFileName = configFileName.replace('_cfg.py', '.log')
        logFileName_full = os.path.join(logFilePath, logFileName)
        logFileNames.append(logFileName)

    retVal = {}
    retVal['inputFileNames']  = inputFileNames_sample
    retVal['configFileNames'] = configFileNames
    retVal['outputFileNames'] = outputFileNames
    retVal['logFileNames']    = logFileNames

    #print "retVal['inputFileNames']  = %s" % inputFileNames_sample
    #print "retVal['outputFileNames'] = %s" % outputFileNames

    return retVal

def buildConfigFile_makeTauFakeRatePlots(inputFileName,
                                         eventSelections, eventSelectionsToPlot, tauIds, tauIdsToPlot, labels,
                                         outputFilePath, outputFileName, recoSampleDefinitions):

    """Make plots of jet --> tau fake-rate as function of jetPt, jetEta,..."""

    processes_string = ""
    processesToPlot = []
    for processName, processConfig in recoSampleDefinitions['MERGE_SAMPLES'].items():
        processes_string += "        " + "cms.PSet(" + "\n"
        processes_string += "        " + "    name = cms.string('%s'),\n" % processName
        processes_string += "        " + "    type = cms.string('%s'),\n" % processConfig['type']
        processes_string += "        " + "    legendEntry = cms.string('%s'),\n" % processConfig['legendEntry']
        drawOptions = processConfig['drawOption']
        # represent all attributes of drawOption object in string format
        for drawOptionAttrName in dir(drawOptions):
            drawOptionAttr = getattr(drawOptions, drawOptionAttrName)
            if isinstance(drawOptionAttr, cms._ParameterTypeBase):
                processes_string += "        " + "    %s = cms.%s(%s),\n" % \
                  (drawOptionAttrName, drawOptionAttr.__class__.__name__, drawOptionAttr.pythonValue())
        processes_string += "        " + ")," + "\n"
        processesToPlot.append(processName)

    eventSelections_string = make_drawOptions_string(eventSelections, eventSelectionsToPlot)
    
    tauIds_string = make_drawOptions_string(tauIds, tauIdsToPlot)

    outputFileName_full = os.path.join(outputFilePath, outputFileName)        

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('%s')
)

process.makeTauFakeRatePlots = cms.PSet(

    processes = cms.VPSet(
%s
    ),
    processesToPlot = cms.vstring(
%s
    ),

    eventSelections = cms.VPSet(
%s
    ),
    eventSelectionsToPlot = cms.vstring(
%s
    ),

    tauIds = cms.VPSet(
%s
    ),
    tauIdsToPlot = cms.vstring(
%s
    ),

    labels = cms.vstring(
%s
    ),

    outputFileName = cms.string('%s')
)
""" % (inputFileName,
       processes_string, make_inputFileNames_vstring(processesToPlot),
       eventSelections_string, make_inputFileNames_vstring(eventSelectionsToPlot),
       tauIds_string, make_inputFileNames_vstring(tauIdsToPlot),
       make_inputFileNames_vstring(labels),
       outputFileName_full)

    configFileName = outputFileName;
    configFileName = configFileName.replace('.eps', '_cfg.py')
    configFileName = configFileName.replace('.pdf', '_cfg.py')
    configFileName = configFileName.replace('.png', '_cfg.py')
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)  

    retVal = {}
    retVal['configFileName'] = configFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal
