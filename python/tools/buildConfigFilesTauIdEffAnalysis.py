import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os

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

def make_binning_string(binning):
    retVal = ""
    for binVariable in binning.keys():
        retVal += "        " + binVariable + " = cms.VPSet(" + "\n"
        for binName, binOptions in binning[binVariable].items():
            if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
                retVal += "            " + "cms.PSet(" + "\n"
                retVal += "                " + "subdir = cms.string('%s'),\n" % binName
                retVal += "                " + "min = cms.double(%f),\n" % binOptions['min']
                retVal += "                " + "max = cms.double(%f)\n" % binOptions['max']
                retVal += "            " + ")," + "\n"
        retVal += "        " + ")," + "\n"
    return retVal;

def getStringRep_bool(flag):
    retVal = None
    if flag:
        retVal = "True"
    else:
        retVal = "False"
    return retVal
#--------------------------------------------------------------------------------

def buildConfigFile_FWLiteTauIdEffAnalyzer(sampleToAnalyze, jobId, inputFilePath, tauIds, binning, sysUncertainties, outputFilePath,
                                           recoSampleDefinitions, passed_region, failed_region, tauChargeMode, disableTauCandPreselCuts):

    """Build cfg.py file to run FWLiteTauIdEffAnalyzer macro to run on PAT-tuples,
       apply event selections and fill histograms for A/B/C/D regions"""

    print "<buildConfigFile_FWLiteTauIdEffAnalyzer>:"
    print " processing sample %s" % sampleToAnalyze

    # CV: check that tauChargeMode parameter matches either of the modes 
    #     define in TauAnalysis/TauIdEfficiency/src/TauIdEffEventSelector.cc
    if not tauChargeMode == "tauLeadTrackCharge" or tauChargeMode == "tauSignalChargedHadronSum":
        raise ValueError("Invalid configuration parameter 'tauChargeMode' = %s !!" % tauChargeMode)

    disableTauCandPreselCuts_string = getStringRep_bool(disableTauCandPreselCuts)

    inputFileNames = os.listdir(inputFilePath)
    #print(inputFileNames)

    # check if inputFile is PAT-tuple and
    # matches sampleToAnalyze, jobId
    inputFile_regex = \
      r"tauIdEffMeasPATTuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)
    inputFile_matcher = re.compile(inputFile_regex)
    inputFileNames = []
    for inputFileName in inputFileNames:
        if inputFile_matcher.match(inputFileName):
            inputFileNames_sample.append(os.path.join(inputFilePath, inputFileName))

    #print(sampleToAnalyze)
    #print(inputFileNames_sample)
    print " found %i input files." % len(inputFileNames_sample)
    
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

    print(" building config file...")

    processType = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['type']

    sysUncertainties_expanded = [ "CENTRAL_VALUE" ]
    if processType != 'Data':
        for sysUncertainty in sysUncertainties:
            if sysUncertainty != "sysAddPUsmearing":
                sysUncertainties_expanded.append(sysUncertainty + "Up")
                sysUncertainties_expanded.append(sysUncertainty + "Down")
            else:
                sysUncertainties_expanded.append(sysUncertainty)
    print " sysUncertainties = %s" %  sysUncertainties_expanded     

    inputFileNames_string = make_inputFileNames_vstring(inputFileNames_sample)

    tauIds_string = make_tauIds_string(tauIds)

    binning_string = make_binning_string(binning)

    configFileNames = []
    outputFileNames = []
    logFileNames    = []

    for sysUncertainty in sysUncertainties_expanded:

        outputFileName = None
        if sysUncertainty != "CENTRAL_VALUE":            
            outputFileName = 'analyzeTauIdEffHistograms_%s_%s_%s.root' % (sampleToAnalyze, sysUncertainty, jobId)
        else:
            outputFileName = 'analyzeTauIdEffHistograms_%s_%s.root' % (sampleToAnalyze, jobId)
        outputFileName_full = os.path.join(outputFilePath, outputFileName)

        srcMuTauPairs = None
        if sysUncertainty != "CENTRAL_VALUE":  
            srcMuTauPairs = composeModuleName([ 'selectedMuPFTauHPSpairsDzForTauIdEff', sysUncertainty, "cumulative" ])            
        else:
            srcMuTauPairs = 'selectedMuPFTauHPSpairsDzForTauIdEffCumulative'

        weights_string = ""
        if not recoSampleDefinitions['MERGE_SAMPLES'][process_matched]['type'] == 'Data':
            weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight", "'"])
            #weights_string += "".join(["'", "ntupleProducer:tauIdEffNtuple#selectedPatMuonsForTauIdEffTrkIPcumulative#muonHLTeff", "'"])

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
    fileNames   = cms.vstring(%s),
    
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.tauIdEffAnalyzer = cms.PSet(
    process = cms.string('%s'),
    type = cms.string('%s'),

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',
        'A1p',
        'A1f',
        'B',
        'B1',
        'B1p',
        'B1f',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D',
        'D1',
        'D1p',
        'D1f'
    ),

    passed_region = cms.string('%s'),
    failed_region = cms.string('%s'),

    regionQCDtemplateFromData_passed = cms.string('B1p'),
    regionQCDtemplateFromData_failed = cms.string('B1f'),
    
    tauIds = cms.VPSet(
%s
    ),
    
    binning = cms.PSet(
%s
    ),
    
    sysShift = cms.string('%s'),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5', 'HLT_IsoMu17_v6', 'HLT_IsoMu17_v8', 'HLT_IsoMu17_v9', 'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('%s'),
    svFitMassHypothesis = cms.string('psKine_MEt_logM_fit'),
    tauChargeMode = cms.string('%s'),
    disableTauCandPreselCuts = cms.bool(%s),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(%s),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
""" % (inputFileNames_string, outputFileName_full,
       process_matched, processType, passed_region, failed_region, tauIds_string, binning_string, sysUncertainty,
       srcMuTauPairs, tauChargeMode, disableTauCandPreselCuts_string, weights_string, allEvents_DBS, xSection, intLumiData)

        outputFileNames.append(outputFileName_full)

        configFileName = None
        if sysUncertainty != "CENTRAL_VALUE":  
            configFileName = "analyzeTauIdEffPATtuple_%s_%s_%s_cfg.py" % (sampleToAnalyze, sysUncertainty, jobId)
        else:
            configFileName = "analyzeTauIdEffPATtuple_%s_%s_cfg.py" % (sampleToAnalyze, jobId)
        configFileName_full = os.path.join(outputFilePath, configFileName)    
        configFile = open(configFileName_full, "w")
        configFile.write(config)
        configFile.close()
        configFileNames.append(configFileName_full)

        logFileName = configFileName.replace('_cfg.py', '.log')
        logFileName_full = os.path.join(outputFilePath, logFileName)
        logFileNames.append(logFileName_full)

    retVal = {}
    retVal['configFileNames'] = configFileNames
    retVal['outputFileNames'] = outputFileNames
    retVal['logFileNames']    = logFileNames

    return retVal

def buildConfigFile_fitTauIdEff(fitMethod, jobId, directory, inputFileName, tauIds, fitVariables, outputFilePath,
                                regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, makeControlPlots):

    """Fit Ztautau signal plus background templates to Mt and visMass distributions
       observed in regions A/B/C/D, in order to determined Ztautau signal contribution
       in regions C1p (tau id. passed sample), C1f (tau id. failed sample)"""

    outputFileName = None
    if directory != '':
        outputFileName = '%s_%s_%s.root' % (fitMethod, jobId, directory)
    else:
        outputFileName = '%s_%s.root' % (fitMethod, jobId)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    makeControlPlots_string = getStringRep_bool(makeControlPlots)

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('%s')
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.%s = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string('%s'),

    runClosureTest = cms.bool(False),
    #runClosureTest = cms.bool(True),

    #takeQCDfromData = cms.bool(False),
    takeQCDfromData = cms.bool(True),

    # CV: fitting fake-rates of background processes
    #     in C2f/C2p regions causes bias of fit result (2011/06/28)
    fitTauIdEffC2 = cms.bool(False),
    #fitTauIdEffC2 = cms.bool(True),
    
    runSysUncertainties = cms.bool(False),
    #runSysUncertainties = cms.bool(True),

    numPseudoExperiments = cms.uint32(10000),

    regions = cms.vstring(
        'ABCD',
        'A',
        'A1',  # QCD enriched control region (OS, loose muon isolation, Mt && Pzeta cuts applied)
        'B',
        'B1',  # QCD enriched control region (SS, loose muon isolation, Mt && Pzeta cuts applied)
        'B1p',
        'B1f',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D',   # generic background control region (SS, tight muon isolation)
        #'D1',
        #'D1p',
        #'D1f',
        #'D2',
        #'D2p',
        #'D2f'
    ),
    regionQCDtemplateFromData_passed = cms.string('%s'),
    regionQCDtemplateFromData_failed = cms.string('%s'),
    
    tauIds = cms.vstring(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),

    sysUncertainties = cms.vstring(
        "sysTauJetEn", # needed for diTauVisMass/diTauVisMassFromJet
        "sysJetEnUp"   # needed for diTauMt
    ),

    makeControlPlots = cms.bool(%s)
)
""" % (inputFileName, outputFileName_full,
       fitMethod, directory, regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed,
       tauIds, fitVariables, makeControlPlots_string)

    configFileName = outputFileName.replace('.root', '_cfg.py')
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

def buildConfigFile_FWLiteTauIdEffPreselNumbers(inputFilePath, sampleZtautau, jobId, tauIds, binning, outputFilePath,
                                                tauChargeMode, disableTauCandPreselCuts):

    """Compute preselection efficiencies and purities in regions C1p, C1f"""

    # CV: check that tauChargeMode parameter matches either of the modes 
    #     define in TauAnalysis/TauIdEfficiency/src/TauIdEffEventSelector.cc
    if not tauChargeMode == "tauLeadTrackCharge" or tauChargeMode == "tauSignalChargedHadronSum":
        raise ValueError("Invalid configuration parameter 'tauChargeMode' = %s !!" % tauChargeMode)

    disableTauCandPreselCuts_string = getStringRep_bool(disableTauCandPreselCuts)

    inputFileNames_Ztautau = []
    inputFileNames = os.listdir(inputFilePath)
    for inputFileName in inputFileNames:
        if inputFileName.find(sampleZtautau) != -1 and inputFileName.find(jobId) != -1:
            inputFileNames_Ztautau.append(os.path.join(inputFilePath, inputFileName))

    inputFileNames_string = make_inputFileNames_vstring(inputFileNames_Ztautau)        

    tauIds_string = make_tauIds_string(tauIds)

    binning_string = make_binning_string(binning)

    outputFileName = 'compTauIdEffPreselNumbers_%s_%s.root' % (sampleZtautau, jobId)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
    
    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(%s),
        
    maxEvents   = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.tauIdEffPreselNumbers = cms.PSet(
    process = cms.string('Ztautau'),

    regions = cms.vstring(
        'C1',
        'C1p',
        'C1f', # for tau id. efficiency measurement
        'D1p'  # for measurement of tau charge misidentification rate
    ),
    
    tauIds = cms.VPSet(
%s
    ),

    binning = cms.PSet(
%s
    ),

    sysShift = cms.string('CENTRAL_VALUE'),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(
        'HLT_IsoMu17_v5', 'HLT_IsoMu17_v6', 'HLT_IsoMu17_v8', 'HLT_IsoMu17_v9', 'HLT_IsoMu17_v11'
    ),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),
    tauChargeMode = cms.string('%s'),
    disableTauCandPreselCuts = cms.bool(%s),
    
    srcGenParticles = cms.InputTag('genParticles'),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag('ntupleProducer:tauIdEffNtuple#addPileupInfo#vtxMultReweight')
)
""" % (inputFileNames_string, outputFileName_full,
       tauIds_string, binning_string, tauChargeMode, disableTauCandPreselCuts_string)

    configFileName = outputFileName.replace('.root', '_cfg.py')
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

def buildConfigFile_compTauIdEffFinalNumbers(inputFileName, directory, jobId, tauIds, fitVariables, outputFilePath,
                                             passed_region, failed_region):

    """Compute final tau id. efficiency values and uncertainties"""

    outputFileName = None
    if directory != '':
        outputFileName = 'compTauIdEffFinalNumbers_%s_%s.root' % (directory, jobId)
    else:
        outputFileName = 'compTauIdEffFinalNumbers_%s.root' % jobId
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring('%s')
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('%s')
)

process.compTauIdEffFinalNumbers = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string('%s'),
    
    tauIds = cms.vstring(
%s
    ),

    fitVariables = cms.vstring(
%s
    )

    passed_region = cms.string('%s'),
    failed_region = cms.string('%s')
)
""" % (inputFileName, outputFileName_full,
       directory, tauIds, fitVariables, passed_region, failed_region)

    configFileName = outputFileName.replace('.root', '_cfg.py')
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

def buildConfigFile_makeTauIdEffFinalPlots(inputFileName, tauIds, tauIdNamesToPlot, binning, fitVariables, outputFilePath, outputFileName,
                                           expEff_label, measEff_label):

    """Make plots of tau id. efficiency as function of tauPt, tauEta,..."""

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

process.makeTauIdEffFinalPlots = cms.PSet(

    tauIds = cms.VPSet(
%s
    ),

    expEff_label  = cms.string('%s'),
    measEff_label = cms.string('%s'),

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
""" % (inputFileName, tauIds_string, expEff_label, measEff_label,
       fitVariables, xAxisBinning_string, binning['xAxisTitle'], values_string, outputFileName_full)

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

def buildConfigFile_hadd(haddCommand, shellFileName_full, inputFileNames, outputFileName_full):

    """Build shell script to run 'hadd' command in order to add all histograms
       in files specified by inputFileNames argument and write the sum to file outputFileName"""

    shellFile = open(shellFileName_full, "w")
    shellFile.write("#!/bin/csh -f\n")
    shellFile.write("\n")
    haddCommandLine = "%s %s" % (haddCommand, outputFileName_full)
    for inputFileName in inputFileNames:
        haddCommandLine += " %s" % inputFileName
    shellFile.write("%s\n" % haddCommandLine)
    shellFile.close()

    logFileName_full = shellFileName_full.replace('.csh', '.log')

    retVal = {}
    retVal['shellFileName']  = shellFileName_full
    retVal['outputFileName'] = outputFileName_full
    retVal['logFileName']    = logFileName_full

    return retVal
