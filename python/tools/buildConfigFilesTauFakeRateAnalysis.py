import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
import TauAnalysis.Configuration.tools.castor as castor

import os
import re

#--------------------------------------------------------------------------------
#
# define auxiliary functions
#
def make_inputFileNames_vstring(list_of_strings):
    retVal = "'"
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += "', '"
        retVal += string_i
    retVal += "'"
    return retVal

def make_tauIds_string(tauIds):
    retVal = ""
    for tauIdName, tauId in tauIds.items():
        retVal += "        " + "cms.PSet(" + "\n"
        retVal += "        " + "    discriminators = cms.vstring(" + "\n"
        for discriminator in tauId['discriminators']:
            retVal += "        " + "        " + "'" + discriminator + "'," + "\n"
        retVal += "        " + "    )," + "\n"
        retVal += "        " + "    name = cms.string('" + tauIdName + "')" + "\n"
        retVal += "        " + ")," + "\n"
    return retVal
#--------------------------------------------------------------------------------

def buildConfigFile_FWLiteTauFakeRateAnalyzer(sampleToAnalyze, evtSel, version, inputFilePath, tauIds, srcTauJetCandidates, srcMET,
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
            if inputFilePath.find('/castor/') != -1:
                # CV: assume that input file gets copied to local directory before FWLiteTauFakeRateAnalyzer macro gets started
                #inputFileNames_sample.append(os.path.join("rfio:" + inputFilePath, inputFileName))
                inputFileNames_sample.append(inputFileName)
            else:
                inputFileNames_sample.append(os.path.join(inputFilePath, inputFileName))

    #print(sampleToAnalyze)
    #print(inputFiles_sample)

    if len(inputFileNames_sample) == 0:
        print("Sample %s has no input files --> skipping !!" % sampleToAnalyze)
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

    print("building config file(s) for sample %s..." % sampleToAnalyze)

    processType = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['type']

    tauIds_string = make_tauIds_string(tauIds)

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
        outputFileName_full = os.path.join(outputFilePath, outputFileName)

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
        'A'
    ),
    
    tauIds = cms.VPSet(
%s
    ),
    
    srcTauJetCandidates = cms.InputTag('%s'),
    
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
""" % (inputFileName_sample, outputFileName_full,
       process_matched, processType, evtSel,
       tauIds_string, srcTauJetCandidates, srcMET, weights_string, allEvents_DBS, xSection, intLumiData)

        outputFileNames.append(outputFileName_full)

        configFileName = "analyzeTauFakeRatePATtuple_%s_%s_%s_cfg.py" % (evtSel, sampleToAnalyze, jobId)
        configFileName_full = os.path.join(configFilePath, configFileName)    
        configFile = open(configFileName_full, "w")
        configFile.write(config)
        configFile.close()
        configFileNames.append(configFileName_full)

        logFileName = configFileName.replace('_cfg.py', '.log')
        logFileName_full = os.path.join(logFilePath, logFileName)
        logFileNames.append(logFileName_full)

    retVal = {}
    retVal['inputFileNames']  = inputFileNames_sample
    retVal['configFileNames'] = configFileNames
    retVal['outputFileNames'] = outputFileNames
    retVal['logFileNames']    = logFileNames

    #print "retVal['inputFileNames']  = %s" % inputFileNames_sample
    #print "retVal['outputFileNames'] = %s" % outputFileNames

    return retVal

def buildConfigFile_makeTauFakeRatePlots(inputFileName,
                                         eventSelections, eventSelectionsToPlot, tauIds, tauIdNamesToPlot,
                                         outputFilePath, outputFileName):

    """Make plots of jet --> tau fake-rate as function of jetPt, jetEta,..."""

    retVal = {}
    retVal['configFileName'] = "blah"
    retVal['outputFileName'] = "blah"
    retVal['logFileName']    = "blah"

    return retVal

    tauIds_string = ""
    for tauIdNameToPlot in tauIdNamesToPlot:
        tauIds_string += "        " + "cms.PSet(" + "\n"
        tauIds_string += "        " + "    name = cms.string('%s'),\n" % tauIdNameToPlot 
        tauIds_string += "        " + "    legendEntry = cms.string('%s'),\n" % tauIds[tauIdNameToPlot]['legendEntry']
        tauIds_string += "        " + "    markerStyleData = cms.uint32(%u),\n" % tauIds[tauIdNameToPlot]['markerStyleData']
        tauIds_string += "        " + "    markerStyleSim = cms.uint32(%u),\n" % tauIds[tauIdNameToPlot]['markerStyleSim']
        tauIds_string += "        " + "    color = cms.uint32(%u)\n" % tauIds[tauIdNameToPlot]['color']
        tauIds_string += "        " + ")," + "\n"

    def xAxisBinValue(binning_item):
        # CV: index [1] refers to binOptions
        #    (index [0] refers to binName)
        if isinstance(binning_item[1], dict):
            return binning_item[1]['min']
        else:
            return -1
    
    xAxisBinning_set = set()
    for binName, binOptions in sorted(binning.items(), key = xAxisBinValue):
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            xAxisBinning_set.add(binOptions['min'])
            xAxisBinning_set.add(binOptions['max'])
    xAxisBinning_string = ""
    for xAxisBin in sorted(xAxisBinning_set):
        xAxisBinning_string += "%f," % xAxisBin
    
    values_string = ""
    for binName, binOptions in sorted(binning.items(), key = xAxisBinValue):
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            values_string += "        " + "cms.PSet(" + "\n"
            values_string += "        " + "    directory = cms.string('%s'),\n" % binName
            values_string += "        " + "    xBinCenter = cms.double(%f),\n" % (0.5*(binOptions['min'] + binOptions['max']))
            values_string += "        " + ")," + "\n"

    outputFileName_full = os.path.join(outputFilePath, outputFileName)        

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('%s')
)

process.makeTauFakeRatePlots = cms.PSet(

    tauIds = cms.VPSet(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),

    xAxisBinning = cms.vdouble(%s),
    xAxisTitle = cms.string('%s'),

    values = cms.VPSet(
%s    
    ),

    outputFileName = cms.string('%s')
)
""" % (inputFileName, tauIds_string, fitVariables, xAxisBinning_string, binning['xAxisTitle'], values_string, outputFileName_full)

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
