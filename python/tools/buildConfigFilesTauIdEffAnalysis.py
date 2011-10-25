import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName

import os
import re

#--------------------------------------------------------------------------------
#
# define auxiliary functions
#
def make_inputFileNames_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += ","
        retVal += "".join([ "'", string_i, "'" ])        
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
                                           recoSampleDefinitions, regions, passed_region, failed_region, hltPaths,
                                           tauChargeMode, disableTauCandPreselCuts):

    """Build cfg.py file to run FWLiteTauIdEffAnalyzer macro to run on PAT-tuples,
       apply event selections and fill histograms for A/B/C/D regions"""

    print "<buildConfigFile_FWLiteTauIdEffAnalyzer>:"
    print " processing sample %s" % sampleToAnalyze

    # CV: check that tauChargeMode parameter matches either of the modes 
    #     define in TauAnalysis/TauIdEfficiency/src/TauIdEffEventSelector.cc
    if not (tauChargeMode == "tauLeadTrackCharge" or tauChargeMode == "tauSignalChargedHadronSum"):
        raise ValueError("Invalid configuration parameter 'tauChargeMode' = %s !!" % tauChargeMode)

    disableTauCandPreselCuts_string = getStringRep_bool(disableTauCandPreselCuts)

    inputFileNames = os.listdir(inputFilePath)
    #print(inputFileNames)

    # check if inputFile is PAT-tuple and
    # matches sampleToAnalyze, jobId
    inputFile_regex = \
      r"tauIdEffMeasPATTuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)
    inputFileNames_sample = []
    for inputFileName in inputFileNames:
        inputFile_matcher = re.compile(inputFile_regex)
        if inputFile_matcher.match(inputFileName):
            inputFileNames_sample.append(os.path.join(inputFilePath, inputFileName))

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

    regions_string = make_inputFileNames_vstring(regions)

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

    hltPaths_string = make_inputFileNames_vstring(hltPaths)

    srcGenParticles = ''
    fillGenMatchHistograms = False
    if processType != 'Data' and \
      (process_matched.find('Ztautau') != -1 or process_matched.find('Zmumu') != -1 or process_matched.find('ZplusJets') != -1):
        srcGenParticles = 'genParticles'
        fillGenMatchHistograms = True

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

        srcMuTauPairs = 'selectedMuPFTauHPSpairsDzForTauIdEff'
        if recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['applyZrecoilCorrection']:
            srcMuTauPairs = composeModuleName([ srcMuTauPairs, "ZllRecoilCorrected" ])
        if sysUncertainty != "CENTRAL_VALUE":  
            srcMuTauPairs = composeModuleName([ srcMuTauPairs, sysUncertainty, "cumulative" ])            
        else:
            srcMuTauPairs = composeModuleName([ srcMuTauPairs, "cumulative" ])

        weights_string = ""
        if not processType == 'Data':
            weights_string += "".join([
                "'", "vertexMultiplicityReweight3d", "'", ","
                "'", "vertexMultiplicityVsRhoPFNeutralReweight", "'"
            ])

        allEvents_DBS = -1
        xSection = 0.0
        if not processType == 'Data':
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
%s
    ),
    
    tauIds = cms.VPSet(
%s
    ),
    
    binning = cms.PSet(
%s
    ),
    
    sysShift = cms.string('%s'),

    srcTrigger = cms.InputTag('patTriggerEvent'),
    hltPaths = cms.vstring(%s),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('%s'),
    svFitMassHypothesis = cms.string('psKine_MEt_logM_fit'),
    tauChargeMode = cms.string('%s'),
    disableTauCandPreselCuts = cms.bool(%s),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    srcGenParticles = cms.InputTag('%s'),
    fillGenMatchHistograms = cms.bool(%s),
    skipPdgIdsGenParticleMatch = cms.vint32(12, 14, 16),

    weights = cms.VInputTag(%s),

    muonIsoProbExtractor = cms.PSet(
        inputFileName = cms.FileInPath('TauAnalysis/TauIdEfficiency/data/train_kNNmuonIsolation_kNN.weights.xml'),
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

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('totalEventsProcessed'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
""" % (inputFileNames_string, outputFileName_full,
       process_matched, processType,
       regions_string, passed_region, failed_region, tauIds_string, binning_string, sysUncertainty, hltPaths_string,
       srcMuTauPairs, tauChargeMode, disableTauCandPreselCuts_string, srcGenParticles, getStringRep_bool(fillGenMatchHistograms),
       weights_string, allEvents_DBS, xSection, intLumiData)

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

def buildConfigFile_fitTauIdEff(fitMethod, jobId, directory, inputFileName, tauIds, fitVariables, sysUncertainties, outputFilePath,
                                regions, passed_region, failed_region, 
                                regionQCDtemplateFromData_all,
                                regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, makeControlPlots):

    """Fit Ztautau signal plus background templates to Mt and visMass distributions
       observed in regions A/B/C/D, in order to determined Ztautau signal contribution
       in regions C1p (tau id. passed sample), C1f (tau id. failed sample)"""

    outputFileName = '%s_%s_%s.root' % (fitMethod, jobId, directory)
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    regions_string = make_inputFileNames_vstring(regions)

    loadSysUncertainties_string = make_inputFileNames_vstring(sysUncertainties)
    varySysUncertainties = []
    if "sysTauJetEn" in sysUncertainties:
        varySysUncertainties.append("sysTauJetEn")
    varySysUncertainties_string = make_inputFileNames_vstring(varySysUncertainties)

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

    # CV: fitting fake-rates of background processes
    #     in C2f/C2p regions causes bias of fit result (2011/06/28)
    fitTauIdEffC2 = cms.bool(False),
    #fitTauIdEffC2 = cms.bool(True),

    regions = cms.vstring(
%s
    ),

    # regions (in Data) from which templates for QCD background are taken
    regionTakeQCDtemplateFromData_all    = cms.string('%s'),
    regionTakeQCDtemplateFromData_passed = cms.string('%s'),
    regionTakeQCDtemplateFromData_failed = cms.string('%s'),
    #takeQCDfromData = cms.bool(False),
    takeQCDfromData = cms.bool(True),

    # define "passed" and "failed" regions
    passed_region = cms.string('%s'),
    failed_region = cms.string('%s'),
    
    tauIds = cms.vstring(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),

    #allowTemplateMorphing = cms.bool(True), # WARNING: template morphing runs **very** slow !!
    allowTemplateMorphing = cms.bool(False),
    morphSysUncertainty = cms.string(
        "sysTauJetEn"
    ),
    sysVariedByNsigma = cms.double(3.0),
    morphQCDinABD = cms.bool(False),
    
    loadSysUncertainties = cms.vstring(
%s
    ),

    runPseudoExperiments = cms.bool(False),
    #runPseudoExperiments = cms.bool(True),
    numPseudoExperiments = cms.uint32(10000),
    varySysUncertainties = cms.vstring(
%s
    ),

    makeControlPlots = cms.bool(%s)
)
""" % (inputFileName, outputFileName_full,
       fitMethod, directory,
       regions_string,
       regionQCDtemplateFromData_all, regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed,
       passed_region, failed_region, 
       tauIds, fitVariables, loadSysUncertainties_string, varySysUncertainties_string, makeControlPlots_string)

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
                                                hltPaths, keyword_compTauIdEffPreselNumbers):

    """Compute preselection efficiencies and purities in regions C1p, C1f"""

    inputFileNames_Ztautau = []
    inputFileNames = os.listdir(inputFilePath)
    for inputFileName in inputFileNames:
        if inputFileName.find(sampleZtautau) != -1 and inputFileName.find(jobId) != -1:
            inputFileNames_Ztautau.append(os.path.join(inputFilePath, inputFileName))

    inputFileNames_string = make_inputFileNames_vstring(inputFileNames_Ztautau)        

    tauIds_string = make_tauIds_string(tauIds)

    binning_string = make_binning_string(binning)

    hltPaths_string = make_inputFileNames_vstring(hltPaths)

    outputFileName = '%s_%s_%s.root' % (keyword_compTauIdEffPreselNumbers, sampleZtautau, jobId)
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

process.%s = cms.PSet(
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
    hltPaths = cms.vstring(%s),
    
    srcGoodMuons = cms.InputTag('patGoodMuons'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),
    
    srcGenParticles = cms.InputTag('genParticles'),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(
        'vertexMultiplicityReweight3d',
        'vertexMultiplicityVsRhoPFNeutralReweight'
    )
)
""" % (inputFileNames_string, outputFileName_full, keyword_compTauIdEffPreselNumbers,
       tauIds_string, binning_string, hltPaths_string)

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
                                             keyword_compTauIdEffFinalNumbers, passed_region, failed_region):

    """Compute final tau id. efficiency values and uncertainties"""

    outputFileName = '%s_%s_%s.root' % (keyword_compTauIdEffFinalNumbers, directory, jobId)
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring('%s')
)
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('%s')
)

process.%s = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'fitTauIdEff' jobs to be run in parallel as there are bins)
    directory = cms.string('%s'),
    
    tauIds = cms.vstring(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),

    passed_region = cms.string('%s'),
    failed_region = cms.string('%s')
)
""" % (inputFileName, outputFileName_full, keyword_compTauIdEffFinalNumbers,
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
