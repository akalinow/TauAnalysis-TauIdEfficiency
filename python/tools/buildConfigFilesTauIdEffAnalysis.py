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

def buildConfigFile_FWLiteTauIdEffAnalyzer(sampleToAnalyze, jobId, inputFilePath, fwliteInput_firstRun, fwliteInput_lastRun, 
                                           tauIds, binning, sysUncertainties, outputFilePath,
                                           recoSampleDefinitions, regions, intLumiData, hltPaths, srcWeights,
                                           tauChargeMode, disableTauCandPreselCuts,
                                           muonPtMin, tauLeadTrackPtMin, tauAbsIsoMax, caloMEtPtMin, pfMEtPtMin,
                                           plot_hltPaths, fillControlPlots, requireUniqueMuTauPair = False):

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
      r"tauIdEffMeasPATTuple_%s_%s_(?P<hash>[a-zA-Z0-9]*).root" % (sampleToAnalyze, jobId)
    inputFileNames_sample = []
    fwliteInput_fileNames = ""
    for inputFileName in inputFileNames:
        inputFile_matcher = re.compile(inputFile_regex)
        if inputFile_matcher.match(inputFileName):
            inputFileNames_sample.append(os.path.join(inputFilePath, inputFileName))
            fwliteInput_fileNames += "process.fwliteInput.fileNames.append('%s')\n" % os.path.join(inputFilePath, inputFileName)

    print " found %i input files." % len(inputFileNames_sample)
    
    if len(inputFileNames_sample) == 0:
        print("Sample %s has no input files --> skipping !!" % sampleToAnalyze)
        return

    # find name of associated "process"
    process_matched = None
    processes = recoSampleDefinitions['MERGE_SAMPLES'].keys()
    for process in processes:
        process_samples = []
        if 'samples' in recoSampleDefinitions['MERGE_SAMPLES'][process].keys():
            process_samples = recoSampleDefinitions['MERGE_SAMPLES'][process]['samples']
        else:
            process_samples.append(process)
        for sample in process_samples:
            if sample == sampleToAnalyze:
                process_matched = process                
    if process_matched and process_matched.startswith('Data'):
        process_matched = 'Data'

    print("sample = %s: process_matched = %s" % (sampleToAnalyze, process_matched))
    
    if not process_matched:
        print("No process associated to sample %s --> skipping !!" % sampleToAnalyze)
        return

    print(" building config file...")

    processType = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['type']

    regions_string = make_inputFileNames_vstring(regions)

    sysUncertainties_expanded = [ "CENTRAL_VALUE" ]
    if processType != 'Data':
        for sysUncertainty in sysUncertainties:
            sysUncertainties_expanded.append(sysUncertainty + "Up")
            sysUncertainties_expanded.append(sysUncertainty + "Down")
    print " sysUncertainties = %s" %  sysUncertainties_expanded     

    tauIds_string = make_tauIds_string(tauIds)

    binning_string = make_binning_string(binning)

    hltPaths_string = make_inputFileNames_vstring(hltPaths[processType])
    plot_hltPaths_string = make_inputFileNames_vstring(plot_hltPaths)
    weights_string = make_inputFileNames_vstring(srcWeights[processType])

    srcGenParticles = ''
    fillGenMatchHistograms = False
    if processType != 'Data' and \
      (process_matched.find('Ztautau')   != -1 or
       process_matched.find('Zmumu')     != -1 or
       process_matched.find('ZplusJets') != -1 or
       process_matched.find('WplusJets') != -1):
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
        if sysUncertainty not in [ "CENTRAL_VALUE", "CaloMEtResponseUp", "CaloMEtResponseDown" ]:  
            srcMuTauPairs = composeModuleName([ srcMuTauPairs, sysUncertainty, "cumulative" ])            
        else:
            srcMuTauPairs = composeModuleName([ srcMuTauPairs, "cumulative" ])
        srcJets = None
        if not processType == 'Data':
            srcJets = 'patJetsSmearedAK5PF'
        else:
            srcJets = 'patJetsCalibratedAK5PF'

        allEvents_DBS = -1
        xSection = 0.0
        if not processType == 'Data':
            allEvents_DBS = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['events_processed']
            xSection = recoSampleDefinitions['RECO_SAMPLES'][sampleToAnalyze]['x_sec']

        config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),

    firstRun = cms.int32(%i),
    lastRun = cms.int32(%i),
    
    maxEvents = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)

%s
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('%s')
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

    srcHLTresults = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),
    
    srcGoodMuons = cms.InputTag('selectedPatMuonsForTauIdEffPFRelIsoCumulative'),
    
    srcMuTauPairs = cms.InputTag('%s'),
    requireUniqueMuTauPair = cms.bool(%s),
    srcCaloMEt = cms.InputTag('patCaloMetNoHF'),
    srcJets = cms.InputTag('%s'),
    #svFitMassHypothesis = cms.string('psKine_MEt_logM_fit'),
    svFitMassHypothesis = cms.string('psKine_MEt_int'),
    tauChargeMode = cms.string('%s'),
    disableTauCandPreselCuts = cms.bool(%s),
    eventSelCuts = cms.PSet(
        muonPtMin         = cms.double(%f),
        tauLeadTrackPtMin = cms.double(%f),
        tauAbsIsoMax      = cms.double(%f),
        caloMEtPtMin      = cms.double(%f),
        pfMEtPtMin        = cms.double(%f)
    ),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    srcGenParticles = cms.InputTag('%s'),
    fillGenMatchHistograms = cms.bool(%s),
    skipPdgIdsGenParticleMatch = cms.vint32(12, 14, 16),

    plot_hltPaths = cms.vstring(%s),

    weights = cms.VInputTag(%s),
    # CV: restrict event weights to avoid problem with too low Monte Carlo event statistics
    #     and large pile-up reweighting factors in Spring'12 MC production
    minWeight = cms.double(0.),
    maxWeight = cms.double(3.),

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
    #applyMuonIsoWeights = cms.bool(False),
    applyMuonIsoWeights = cms.bool(True),

    # CV: only book histograms needed to run tau id. efficiency fit;
    #     drop all control plots
    fillControlPlots = cms.bool(%s),

    # CV: 'srcEventCounter' is defined in TauAnalysis/Skimming/test/skimTauIdEffSample_cfg.py
    srcEventCounter = cms.InputTag('processedEventsSkimming'),
    allEvents_DBS = cms.int32(%i),
    
    xSection = cms.double(%f),
    
    intLumiData = cms.double(%f),

    srcLumiProducer = cms.InputTag('lumiProducer')
)
""" % (fwliteInput_firstRun, fwliteInput_lastRun, fwliteInput_fileNames, outputFileName_full,
       process_matched, processType,
       regions_string, tauIds_string, binning_string, sysUncertainty, hltPaths_string,
       srcMuTauPairs, getStringRep_bool(requireUniqueMuTauPair), srcJets, tauChargeMode, disableTauCandPreselCuts_string,
       muonPtMin, tauLeadTrackPtMin, tauAbsIsoMax, caloMEtPtMin, pfMEtPtMin,
       srcGenParticles, getStringRep_bool(fillGenMatchHistograms), plot_hltPaths_string,
       weights_string, getStringRep_bool(fillControlPlots), allEvents_DBS, xSection, intLumiData)

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


def buildConfigFile_smoothTauIdEffTemplates(jobId, directory, inputFileName, tauIds, fitVariables, sysUncertainties, outputFilePath,
                                            regions, fitIndividualProcesses, makeControlPlots, outputFilePath_plots):

    """Build config file for smoothing W+jets and QCD by analytic fits"""

    # CV: for normalization purposes, always add 'EventCounter';
    #     add 'diTauMt' in order to take QCD template for region 'D' from data
    fitVariables_extended = [ 'EventCounter', 'diTauMt' ]
    fitVariables_extended.extend(fitVariables)

    sysUncertainties_extended = [ 'CENTRAL_VALUE' ]
    for sysUncertainty in sysUncertainties:
        sysUncertainties_extended.append("%sUp" % sysUncertainty)
        sysUncertainties_extended.append("%sDown" % sysUncertainty)

    if directory != "":
        outputFileName = 'smoothTauIdEffTemplates_%s_%s.root' % (jobId, directory)
    else:
        outputFileName = 'smoothTauIdEffTemplates_%s.root' % jobId
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    histogramsToSmooth_string = ""
    for process_and_region in [ [ 'WplusJets', 'A'     ],
                                [ 'WplusJets', 'A1'    ],
                                [ 'WplusJets', 'A_mW'  ],
                                [ 'WplusJets', 'A1_mW' ],
                                [ 'WplusJets', 'B'     ],
                                [ 'WplusJets', 'C1f'   ],
                                [ 'WplusJets', 'C1p'   ],
                                [ 'WplusJets', 'D'     ],
                                [ 'QCD',       'A'     ],
                                [ 'QCD',       'A1'    ],
                                [ 'QCD',       'A_mW'  ],
                                [ 'QCD',       'A1_mW' ],
                                [ 'QCD',       'B'     ] ]:
        process = process_and_region[0]
        region  = process_and_region[1]
        if not region in regions:
            continue
        for tauId in tauIds:            
            for fitVariable in fitVariables_extended:
                for sysUncertainty in sysUncertainties_extended:
                    for genMatch_option in [ None, "JetToTauFake" ]:
                    
                        histogramName = "%s_%s_%s_%s" % (process, region, fitVariable, tauId)
                        if region.find("p") != -1:
                            histogramName += "_passed"
                        elif region.find("f") != -1:
                            histogramName += "_failed"
                        else:
                            histogramName += "_all"
                        if sysUncertainty != 'CENTRAL_VALUE':
                            histogramName += "_%s" % sysUncertainty
                        if genMatch_option is not None:
                            if process == 'WplusJets':
                                histogramName += "_%s" % genMatch_option
                            else:
                                continue

                        fitFunctionType = None                    
                        if fitVariable == "diTauVisMass":
                            fitFunctionType = "LG1"
                        elif fitVariable == "diTauMt":
                            if region in [ 'D' ]:
                                fitFunctionType = "CB1"
                            elif process == 'WplusJets' and region in [ 'A', 'A1', 'A_mW', 'A1_mW', 'B', 'C1f', 'C1p' ]:
                                fitFunctionType = "CB2"
                            elif process == 'QCD' and region in [ 'A', 'A1', 'A_mW', 'A1_mW', 'B' ]:
                                fitFunctionType = "CB3"
                            else:
                                raise ValueError("Undefined combination of region = %s and process = %s !!" % (region, process))
                        else:
                            continue
 
                        histogramsToSmooth_string += "process.smoothTauIdEffTemplates.histogramsToSmooth.append(\n"
                        histogramsToSmooth_string += "    cms.PSet(\n"
                        histogramsToSmooth_string += "        histogramName = cms.string('%s'),\n" % histogramName
                        histogramsToSmooth_string += "        fitFunctionType = cms.string('%s')\n" % fitFunctionType
                        histogramsToSmooth_string += "    )\n"
                        histogramsToSmooth_string += ")\n"

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

process.smoothTauIdEffTemplates = cms.PSet(
    directory = cms.string('%s'),
    
    histogramsToSmooth = cms.VPSet(),

    makeControlPlots = cms.bool(%s),
    controlPlotFilePath = cms.string('%s')
)

%s
""" % (inputFileName, outputFileName_full,
       directory,       
       getStringRep_bool(makeControlPlots), outputFilePath_plots,
       histogramsToSmooth_string)
    
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

def buildConfigFile_makeTauIdEffQCDtemplate(jobId, directory, inputFileName, tauIds, fitVariables, sysUncertainties, outputFilePath,
                                            regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, regionQCDtemplateFromData_D, 
                                            regionWplusJetsSideband_passed, regionWplusJetsSideband_failed, regionWplusJetsSideband_D,  
                                            histQCDtemplateFromData_passed, histQCDtemplateFromData_failed, histQCDtemplateFromData_D):

    """Build config file for correcting QCD template obtained from control region in Data
       for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds"""

    fitVariables_makeTauIdEffQCDtemplate = []
    fitVariables_makeTauIdEffQCDtemplate.extend(fitVariables)
    fitVariables_always = [
        'diTauVisMass',
        'diTauMt'
    ]
    for fitVariable_always in fitVariables_always:
        if fitVariables_makeTauIdEffQCDtemplate.count(fitVariable_always) == 0:
            fitVariables_makeTauIdEffQCDtemplate.append(fitVariable_always)

    if directory != "":
        outputFileName = 'makeTauIdEffQCDtemplate_%s_%s.root' % (jobId, directory)
    else:
        outputFileName = 'makeTauIdEffQCDtemplate_%s.root' % jobId
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    sysUncertainties_string = make_inputFileNames_vstring(sysUncertainties)

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

process.makeTauIdEffQCDtemplate = cms.PSet(

    # CV: set to '' if determining tau id. efficiency for whole Tag & Probe sample,
    #     set to name of one individual bin in case you want to measure the tau id. efficiency as function of tauPt, tauEta,...
    #    (needs as many 'makeTauIdEffQCDtemplate' jobs to be run in parallel as there are bins)
    directory = cms.string('%s'),

    # regions (in Data) from which templates for QCD background are taken
    regionTakeQCDtemplateFromData_passed = cms.string('%s'),
    regionTakeQCDtemplateFromData_failed = cms.string('%s'),
    regionTakeQCDtemplateFromData_D      = cms.string('%s'),

    # regions (in Data) from which W + jets background contribution to QCD control region is estimated
    regionWplusJetsSideband_passed       = cms.string('%s'),
    regionWplusJetsSideband_failed       = cms.string('%s'),
    regionWplusJetsSideband_D            = cms.string('%s'),

    # define "all", "passed" and "failed" regions
    regionStoreQCDtemplate_passed        = cms.string('%s'),
    regionStoreQCDtemplate_failed        = cms.string('%s'),
    regionStoreQCDtemplate_D             = cms.string('%s'),
    
    tauIds = cms.vstring(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),
    
    sysUncertainties = cms.vstring(
%s
    )
)
""" % (inputFileName, outputFileName_full,
       directory,
       regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, regionQCDtemplateFromData_D,
       regionWplusJetsSideband_passed, regionWplusJetsSideband_failed, regionWplusJetsSideband_D, 
       histQCDtemplateFromData_passed, histQCDtemplateFromData_failed, histQCDtemplateFromData_D,
       tauIds, fitVariables_makeTauIdEffQCDtemplate, sysUncertainties_string)
    
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

def buildConfigFile_fitTauIdEff(fitMethod, jobId, directory, inputFileName, tauId, fitVariable,
                                templateMorphingMode, sysUncertainties, outputFilePath,
                                regionsToFit, passed_region, failed_region, 
                                regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, regionQCDtemplateFromData_D,
                                fitIndividualProcesses, intLumiData, runClosureTest, makeControlPlots, outputFilePath_plots):

    """Fit Ztautau signal plus background templates to Mt and visMass distributions
       observed in regions A/B/C/D, in order to determined Ztautau signal contribution
       in regions C1p (tau id. passed sample), C1f (tau id. failed sample)"""

    if directory != "":  
        outputFileName = '%s_%s_%s_%s_%s.root' % (fitMethod, tauId, fitVariable, jobId, directory)
    else:
        outputFileName = '%s_%s_%s_%s.root' % (fitMethod, tauId, fitVariable, jobId)
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    regionsToFit_string = make_inputFileNames_vstring(regionsToFit)

    fitIndividualProcesses_string = getStringRep_bool(fitIndividualProcesses)

    sysUncertainties_string = make_inputFileNames_vstring(sysUncertainties)

    makeControlPlots_string = getStringRep_bool(makeControlPlots)

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

    runClosureTest = cms.bool(%s),
    
    regions = cms.vstring(
%s
    ),

    # regions (Data or MC) from which templates for QCD background are taken
    regionQCDtemplate_passed = cms.string('%s'),
    regionQCDtemplate_failed = cms.string('%s'),
    regionQCDtemplate_D      = cms.string('%s'),

    # define "passed" and "failed" regions
    region_passed = cms.string('%s'),
    region_failed = cms.string('%s'),
    
    tauId = cms.string('%s'),

    fitVariable = cms.string('%s'),

    fitIndividualProcesses = cms.bool(%s),

    templateMorphingMode = cms.string('%s'), 

    sysUncertainties = cms.vstring(
%s
    ),
    sysVariedByNsigma = cms.double(3.0),

    runPseudoExperiments = cms.bool(False),
    #runPseudoExperiments = cms.bool(True),
    numPseudoExperiments = cms.uint32(10000),

    intLumiData = cms.double(%f),

    makeControlPlots = cms.bool(%s),
    controlPlotFilePath = cms.string('%s')
)
""" % (inputFileName, outputFileName_full,
       fitMethod, directory, getStringRep_bool(runClosureTest),
       regionsToFit_string,
       regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, regionQCDtemplateFromData_D, 
       passed_region, failed_region, 
       tauId, fitVariable, fitIndividualProcesses_string, templateMorphingMode, sysUncertainties_string,
       intLumiData*1.e-3, makeControlPlots_string, outputFilePath_plots)
    
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
                                                hltPaths, srcWeights, keyword_compTauIdEffPreselNumbers,
                                                muonPtMin, tauLeadTrackPtMin, tauAbsIsoMax, caloMEtPtMin, pfMEtPtMin):

    """Compute preselection efficiencies and purities in regions C1p, C1f"""

    inputFileNames = os.listdir(inputFilePath)
    #print(inputFileNames)

    inputFile_regex = \
      r"tauIdEffMeasPATTuple_%s_%s_(?P<hash>[a-zA-Z0-9]*).root" % (sampleZtautau, jobId)
    inputFileNames_sample = []
    fwliteInput_fileNames = ""
    for inputFileName in inputFileNames:
        inputFile_matcher = re.compile(inputFile_regex)
        if inputFile_matcher.match(inputFileName):
            inputFileNames_sample.append(os.path.join(inputFilePath, inputFileName))
            fwliteInput_fileNames += "process.fwliteInput.fileNames.append('%s')\n" % os.path.join(inputFilePath, inputFileName)

    print " found %i input files." % len(inputFileNames_sample)

    tauIds_string = make_tauIds_string(tauIds)

    binning_string = make_binning_string(binning)

    hltPaths_string = make_inputFileNames_vstring(hltPaths['smMC'])
    weights_string = make_inputFileNames_vstring(srcWeights['smMC'])

    outputFileName = '%s_%s_%s.root' % (keyword_compTauIdEffPreselNumbers, sampleZtautau, jobId)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
    
    config = \
"""
import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
        
    maxEvents = cms.int32(-1),
    
    outputEvery = cms.uint32(1000)
)

%s
    
process.fwliteOutput = cms.PSet(
    fileName = cms.string('%s')
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

    srcHLTresults = cms.InputTag('TriggerResults::HLT'),
    hltPaths = cms.vstring(%s),
    
    srcGoodMuons = cms.InputTag('selectedPatMuonsForTauIdEffPFRelIsoCumulative'),
    
    srcMuTauPairs = cms.InputTag('selectedMuPFTauHPSpairsDzForTauIdEffCumulative'),
    srcCaloMEt = cms.InputTag('patCaloMetNoHF'),
    srcJets = cms.InputTag('patJetsSmearedAK5PF'),
    
    srcGenParticles = cms.InputTag('genParticles'),

    eventSelCuts = cms.PSet(
        muonPtMin         = cms.double(%f),
        tauLeadTrackPtMin = cms.double(%f),
        tauAbsIsoMax      = cms.double(%f),
        caloMEtPtMin      = cms.double(%f),
        pfMEtPtMin        = cms.double(%f)
    ),

    srcVertices = cms.InputTag('offlinePrimaryVertices'),

    weights = cms.VInputTag(%s)
)
""" % (fwliteInput_fileNames, outputFileName_full, keyword_compTauIdEffPreselNumbers,
       tauIds_string, binning_string, hltPaths_string, 
       muonPtMin, tauLeadTrackPtMin, tauAbsIsoMax, caloMEtPtMin, pfMEtPtMin,
       weights_string)

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
                                             keyword_compTauIdEffFinalNumbers, passed_region, failed_region,
                                             fitIndividualProcesses):

    """Compute final tau id. efficiency values and uncertainties"""

    outputFileName = '%s_%s_%s.root' % (keyword_compTauIdEffFinalNumbers, directory, jobId)
    outputFileName = outputFileName.replace('__', '_')
    outputFileName_full = os.path.join(outputFilePath, outputFileName)

    tauIds_string = make_tauIds_string(tauIds)

    fitIndividualProcesses_string = getStringRep_bool(fitIndividualProcesses)

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
    
    tauIds = cms.VPSet(
%s
    ),

    fitVariables = cms.vstring(
%s
    ),

    passed_region = cms.string('%s'),
    failed_region = cms.string('%s'),

    fitIndividualProcesses = cms.bool(%s)
)
""" % (inputFileName, outputFileName_full, keyword_compTauIdEffFinalNumbers,
       directory, tauIds_string, fitVariables, passed_region, failed_region, fitIndividualProcesses_string)

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

def buildConfigFile_makeTauIdEffFinalPlots(inputFileName, tauIds, tauIdNamesToPlot, binning, fitVariables,
                                           outputFilePath, outputFileName, expEff_label, measEff_label, intLumiData):

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
    fileNames = cms.vstring('%s')
)

process.makeTauIdEffFinalPlots = cms.PSet(

    tauIds = cms.VPSet(
%s
    ),

    expEff_label = cms.string('%s'),
    measEff_label = cms.string('%s'),

    intLumiData = cms.double(%f),

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
""" % (inputFileName, tauIds_string, expEff_label, measEff_label, intLumiData*1.e-3,
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
    # CV: delete output file in case it exists 
    shellFile.write("rm -f %s\n" % outputFileName_full)
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
