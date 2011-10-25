#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import *

import os

channel = 'ZtoMuTau_tauIdEff'
#jobId = getJobId(channel)
jobId = '2011Aug18'

version = 'V9'

inputFilePath = '/data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/%s/%s/' % (jobId, version) \
               + 'user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/%s/' % jobId
outputFilePath = '/data1/veelken/tmp/muonPtGt20_noPzetaCut/%s/' % version

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v4',
    'Ztautau_powheg',
    'Ztautau_embedded_part1',
    'Ztautau_embedded_part2',
    'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

# define sample name and jobId of Ztautau sample
# used to compute preselection efficiencies and purities in C1p and C1f/D1p regions
sampleZtautau = 'Ztautau_powheg'

hltPaths = [
    #'HLT_Mu15_v1',
    #'HLT_Mu15_v2',
    #'HLT_Mu15_v3',
    #'HLT_Mu15_v4',
    #'HLT_Mu15_v5',
    #'HLT_Mu15_v6'
    'HLT_IsoMu17_v5',
    'HLT_IsoMu17_v6',
    'HLT_IsoMu17_v8',
    'HLT_IsoMu17_v9',
    'HLT_IsoMu17_v10',
    'HLT_IsoMu17_v11'
]

fitVariables = [
    'diTauVisMass'
]    

mode = 'tauIdEfficiency'
#mode = 'tauChargeMisIdRate'

sysUncertainties = [
    "sysTauJetEn",     # needed for diTauVisMass/diTauVisMassFromJet
    "sysJetEn",        # needed for diTauMt
    "sysAddPUsmearing" # additional MET smearing, not correlated with nay particles reconstructed in the event
]

tauIds = {
    # HPS isolation with no deltaBeta corrections applied
    # (separate isolation requirements wrt. PFChargedHadrons and PFGammas)
    'tauDiscrHPSloose'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseIsolation'
        ],
        'legendEntry' : "HPS Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPSmedium' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumIsolation'
        ],
        'legendEntry' : "HPS Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPStight' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightIsolation'
        ],
        'legendEntry' : "HPS Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    },
            
    # HPS isolation with deltaBeta corrections applied
    # (separate isolation requirements wrt. PFChargedHadrons and PFGammas)
    'tauDiscrHPSlooseDBcorr'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPSmediumDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPStightDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS #delta#beta Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    },
    
    # HPS combined isolation discriminators
    # (based on isolation sumPt of PFChargedHadrons + PFGammas)
    'tauDiscrHPScombLooseDBcorr'  : {
        'discriminators' : [
            'decayModeFinding',
            'byLooseCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Loose",
        'markerStyleData' : 20,
        'markerStyleSim' : 24,
        'color' : 418
    },
    'tauDiscrHPScombMediumDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byMediumCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Medium",
        'markerStyleData' : 21,
        'markerStyleSim' : 25,
        'color' : 807
    },
    'tauDiscrHPScombTightDBcorr' : {
        'discriminators' : [
            'decayModeFinding',
            'byTightCombinedIsolationDeltaBetaCorr'
        ],
        'legendEntry' : "HPS comb. Tight",
        'markerStyleData' : 22,
        'markerStyleSim' : 26,
        'color' : 618
    }
}

binning = {
    # NOTE: these strings need to match what is hard-coded in regionEntryType::regionEntryType constructor,
    #       defined in TauAnalysis/TauIdEfficiency/bin/FWLiteTauIdEffAnalyzer.cc 
    'tauPt' : {
        'tauPtLt25' : {
            'min' :  20.0,
            'max' :  25.0
        },
        'tauPt25to30' : {
            'min' :  25.0,
            'max' :  30.0
        },
        'tauPt30to40' : {
            'min' :  30.0,
            'max' :  40.0
        },
        'tauPtGt40' : {
            'min' :  40.0,
            'max' : 100.0
        },
        'xAxisTitle' : "P_{T}^{#tau}"
    },
    'tauAbsEta' : {
        'tauAbsEtaLt14' : {
            'min' :  0.0,
            'max' :  1.4
        },
        'tauAbsEta14to19' : {
            'min' :  1.4,
            'max' :  1.9
        },
        'tauAbsEta19to23' : {
            'min' :  1.9,
            'max' :  2.3
        },
        'xAxisTitle' : "|#eta_{#tau}|"
    },
    'numVertices' : {
        'numVerticesLeq4' : {
            'min' :  -0.5,
            'max' :   4.5
        },
        'numVertices5to6' : {
            'min' :   4.5,
            'max' :   6.5
        },
        'numVertices7to8' : {
            'min' :   6.5,
            'max' :   8.5
        },
        'numVerticesGt8' : {
            'min' : 8.5,
            'max' :  20.5
        },
        'xAxisTitle' : "Num. Vertices"
    },
    'sumEt' : {
        'sumEtLt250' : {
            'min' :    0.0,
            'max' :  250.0
        },
        'sumEt250to350' : {
            'min' :  250.0,
            'max' :  350.0
        },
        'sumEt350to450' : {
            'min' :  350.0,
            'max' :  450.0
        },
        'sumEtGt450' : {
            'min' :  450.0,
            'max' : 1000.0
        },
        'xAxisTitle' : "#Sigma E_{T}"
    }
}

execDir = "%s/bin/%s/" % (os.environ['CMSSW_BASE'], os.environ['SCRAM_ARCH'])

executable_compTauIdEffPreselNumbers = None
suffix_noTauSel                      = None
keyword_compTauIdEffPreselNumbers    = None
passed_region                        = None
failed_region                        = None
regionQCDtemplateFromData_all        = None
regionQCDtemplateFromData_passed     = None
regionQCDtemplateFromData_failed     = None
fitMethod                            = None
tauChargeMode                        = None
disableTauCandPreselCuts             = None
executable_compTauIdEffFinalNumbers  = None
keyword_compTauIdEffFinalNumbers     = None
expEff_label                         = None                
measEff_label                        = None  
if mode == 'tauIdEfficiency':
    executable_compTauIdEffPreselNumbers = 'nice ' + execDir + 'FWLiteTauIdEffPreselNumbers'
    keyword_compTauIdEffPreselNumbers    = 'compTauIdEffPreselNumbers'
    suffix_noTauSel                      = '_noTauSel'
    regions                              = [
        'ABCD',
        'A',
        'A1',  # QCD enriched control region (OS, loose muon isolation, Mt && Pzeta cuts applied)
        'A1_mW',
        #'A1_mW_tW',
        'A1p',
        'A1f',
        'AWj',
        'AWj_mW',
        #'AWj_mW_tW',
        'B',
        'B1',  # QCD enriched control region (SS, loose muon isolation, Mt && Pzeta cuts applied)
        'B1_mW',
        #'B1_mW_tW',
        'B1p',
        'B1f',
        'BWj',
        'BWj_mW',
        #'BWj_mW_tW',
        'C',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'C2p',
        'C2f',
        'D',   # generic background control region (SS, tight muon isolation)
        'D1',
        'D1p',
        'D1f'
    ]
    all_region                           = 'C1'
    passed_region                        = 'C1p'
    failed_region                        = 'C1f'
    regionQCDtemplateFromData_all        = 'B1'
    regionQCDtemplateFromData_passed     = 'C1_qcd'
    regionQCDtemplateFromData_failed     = 'D1_qcd'
    fitMethod                            = 'fitTauIdEff_wConstraints'
    tauChargeMode                        = 'tauLeadTrackCharge'
    disableTauCandPreselCuts             = False
    executable_compTauIdEffFinalNumbers  = execDir + 'compTauIdEffFinalNumbers'
    keyword_compTauIdEffFinalNumbers     = 'compTauIdEffFinalNumbers'
    expEff_label                         = 'expEff'
    measEff_label                        = 'measEff'   
elif mode == 'tauChargeMisIdRate':
    executable_compTauIdEffPreselNumbers = 'nice ' + execDir + 'FWLiteTauChargeMisIdPreselNumbers'
    keyword_compTauIdEffPreselNumbers    = 'compTauChargeMisIdPreselNumbers'
    suffix_noTauSel                      = ''
    regions                              = [
        'B1',  # control region used to obtain QCD template from Data
        'B1p',
        'B1f',
        'B2',
        'C1',
        'C1p',
        'C1f',
        'C2',
        'D1',
        'D1p',
        'D1f',
        'D2'
    ]
    all_region                           = ''
    passed_region                        = 'C1p'
    failed_region                        = 'D1p'
    regionQCDtemplateFromData_all        = 'B1'
    regionQCDtemplateFromData_passed     = 'B1p'
    regionQCDtemplateFromData_failed     = 'B1f'
    fitMethod                            = 'fitTauIdEff'
    tauChargeMode                        = 'tauSignalChargedHadronSum'
    disableTauCandPreselCuts             = True
    executable_compTauIdEffFinalNumbers  = execDir + 'compTauChargeMisIdFinalNumbers'
    keyword_compTauIdEffFinalNumbers     = 'compTauChargeMisIdFinalNumbers'
    expEff_label                         = 'expRate'
    measEff_label                        = 'measRate'   
else:
    raise ValueError("Invalid mode = %s !!" % mode)

executable_FWLiteTauIdEffAnalyzer = 'nice ' + execDir + 'FWLiteTauIdEffAnalyzer'
executable_hadd = 'hadd -f'
executable_fitTauIdEff = execDir + fitMethod
executable_makeTauIdEffFinalPlots = execDir + 'makeTauIdEffFinalPlots'
executable_shell = '/bin/csh'

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'].keys()

outputFilePath = os.path.join(outputFilePath, mode)
if not os.path.exists(outputFilePath):
    os.mkdir(outputFilePath)

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauIdEffAnalyzer macro
#
configFileNames_FWLiteTauIdEffAnalyzer = []
outputFileNames_FWLiteTauIdEffAnalyzer = []
logFileNames_FWLiteTauIdEffAnalyzer    = []
for sampleToAnalyze in samplesToAnalyze:
    retVal_FWLiteTauIdEffAnalyzer = \
      buildConfigFile_FWLiteTauIdEffAnalyzer(sampleToAnalyze, "".join([ jobId, version ]), inputFilePath, tauIds,
                                             binning, sysUncertainties, outputFilePath,
                                             recoSampleDefinitionsTauIdEfficiency_7TeV,
                                             regions, passed_region, failed_region, hltPaths,
                                             tauChargeMode, disableTauCandPreselCuts)

    if retVal_FWLiteTauIdEffAnalyzer is None:
        continue
    
    configFileNames_FWLiteTauIdEffAnalyzer.extend(retVal_FWLiteTauIdEffAnalyzer['configFileNames'])
    outputFileNames_FWLiteTauIdEffAnalyzer.extend(retVal_FWLiteTauIdEffAnalyzer['outputFileNames'])
    logFileNames_FWLiteTauIdEffAnalyzer.extend(retVal_FWLiteTauIdEffAnalyzer['logFileNames'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "harvest" histograms
# produced by FWLiteTauIdEffAnalyzer macro
#
haddShellFileName_stage1 = os.path.join(outputFilePath, 'harvestTauIdEffHistograms_stage1_%s.csh' % "".join([ jobId, version ]))
haddInputFileNames_stage1 = outputFileNames_FWLiteTauIdEffAnalyzer
haddOutputFileName_stage1 = os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_all_%s.root' % "".join([ jobId, version ]))
retVal_hadd_stage1 = \
  buildConfigFile_hadd(executable_hadd, haddShellFileName_stage1, haddInputFileNames_stage1, haddOutputFileName_stage1)
haddLogFileName_stage1 = retVal_hadd_stage1['logFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running fitTauIdEff macro
#
configFileNames_fitTauIdEff = []
outputFileNames_fitTauIdEff = []
logFileNames_fitTauIdEff    = []
retVal_fitTauIdEff = \
  buildConfigFile_fitTauIdEff(fitMethod, "".join([ jobId, version ]), '', haddOutputFileName_stage1, tauIds.keys(),
                              fitVariables, sysUncertainties, outputFilePath,
                              regions, passed_region, failed_region,
                              regionQCDtemplateFromData_all,
                              regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, True)
configFileNames_fitTauIdEff.append(retVal_fitTauIdEff['configFileName'])
outputFileNames_fitTauIdEff.append(retVal_fitTauIdEff['outputFileName'])
logFileNames_fitTauIdEff.append(retVal_fitTauIdEff['logFileName'])
for binVariable in binning.keys():
    for binName, binOptions in binning[binVariable].items():
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            retVal_fitTauIdEff = \
              buildConfigFile_fitTauIdEff(fitMethod, "".join([ jobId, version ]), binName, haddOutputFileName_stage1, tauIds.keys(),
                                          fitVariables, sysUncertainties, outputFilePath,
                                          regions, passed_region, failed_region,
                                          regionQCDtemplateFromData_all,
                                          regionQCDtemplateFromData_passed, regionQCDtemplateFromData_failed, False)
            configFileNames_fitTauIdEff.append(retVal_fitTauIdEff['configFileName'])
            outputFileNames_fitTauIdEff.append(retVal_fitTauIdEff['outputFileName'])
            logFileNames_fitTauIdEff.append(retVal_fitTauIdEff['logFileName'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauIdEffPreselNumbers macro
#
retVal_FWLiteTauIdEffPreselNumbers = \
  buildConfigFile_FWLiteTauIdEffPreselNumbers(inputFilePath, sampleZtautau, "".join([ jobId, version, suffix_noTauSel ]), tauIds,
                                              binning, outputFilePath, hltPaths, keyword_compTauIdEffPreselNumbers)
configFileName_FWLiteTauIdEffPreselNumbers = retVal_FWLiteTauIdEffPreselNumbers['configFileName']
outputFileName_FWLiteTauIdEffPreselNumbers = retVal_FWLiteTauIdEffPreselNumbers['outputFileName']
logFileName_FWLiteTauIdEffPreselNumbers = retVal_FWLiteTauIdEffPreselNumbers['logFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "harvest" fit results
# with preselection efficiencies and purities in regions C1p/C1f
#
haddShellFileName_stage2 = os.path.join(outputFilePath, 'harvestTauIdEffHistograms_stage2_%s.csh' % "".join([ jobId, version ]))
haddInputFileNames_stage2 = []
haddInputFileNames_stage2.extend(outputFileNames_fitTauIdEff)
haddInputFileNames_stage2.append(outputFileName_FWLiteTauIdEffPreselNumbers)
haddOutputFileName_stage2 = os.path.join(outputFilePath, 'compTauIdEffFinalNumbers_input_%s.root' % "".join([ jobId, version ]))
retVal_hadd_stage2 = \
  buildConfigFile_hadd(executable_hadd, haddShellFileName_stage2, haddInputFileNames_stage2, haddOutputFileName_stage2)
haddLogFileName_stage2 = retVal_hadd_stage2['logFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running compTauIdEffFinalNumbers macro
#
configFileNames_compTauIdEffFinalNumbers = []
outputFileNames_compTauIdEffFinalNumbers = []
logFileNames_compTauIdEffFinalNumbers = []
retVal_compTauIdEffFinalNumbers = \
  buildConfigFile_compTauIdEffFinalNumbers(haddOutputFileName_stage2, '', "".join([ jobId, version ]), tauIds.keys(),
                                           fitVariables, outputFilePath,
                                           keyword_compTauIdEffFinalNumbers, passed_region, failed_region)
configFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['configFileName'])
outputFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['outputFileName'])
logFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['logFileName'])
for binVariable in binning.keys():
    for binName, binOptions in binning[binVariable].items():
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            retVal_compTauIdEffFinalNumbers = \
              buildConfigFile_compTauIdEffFinalNumbers(haddOutputFileName_stage2, binName, "".join([ jobId, version ]), tauIds.keys(),
                                                       fitVariables, outputFilePath,
                                                       keyword_compTauIdEffFinalNumbers, passed_region, failed_region)
            configFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['configFileName'])
            outputFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['outputFileName'])
            logFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['logFileName'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "harvest" final tau id. efficiency numbers
# for different tauPt, tauEta,... bins
#
haddShellFileName_stage3 = os.path.join(outputFilePath, 'harvestTauIdEffHistograms_stage3_%s.csh' % "".join([ jobId, version ]))
haddInputFileNames_stage3 = outputFileNames_compTauIdEffFinalNumbers
haddOutputFileName_stage3 = os.path.join(outputFilePath, 'compTauIdEffFinalNumbers_all_%s.root' % "".join([ jobId, version ]))
retVal_hadd_stage3 = \
  buildConfigFile_hadd(executable_hadd, haddShellFileName_stage3, haddInputFileNames_stage3, haddOutputFileName_stage3)
haddLogFileName_stage3 = retVal_hadd_stage3['logFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# make final tau id. efficiency plots as function of tauPt, tauEta,...
#
configFileNames_makeTauIdEffFinalPlots = []
outputFileNames_makeTauIdEffFinalPlots = []
logFileNames_makeTauIdEffFinalPlots    = []
for binVariable in binning.keys():

    # make plots for HPS isolation with no deltaBeta corrections applied
    discriminators_HPS = [
        'tauDiscrHPSloose',
        'tauDiscrHPSmedium',
        'tauDiscrHPStight'
    ]
    outputFileName_HPS = 'makeTauIdEffFinalPlots_HPS_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPS,
                                             binning[binVariable], fitVariables, outputFilePath, outputFileName_HPS,
                                             expEff_label, measEff_label)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])
    logFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['logFileName'])

    # make plots for HPS isolation with no applied deltaBeta corrections
    discriminators_HPSdbCorr = [
        'tauDiscrHPSlooseDBcorr',
        'tauDiscrHPSmediumDBcorr',
        'tauDiscrHPStightDBcorr'
    ]
    outputFileName_HPSdbCorr = 'makeTauIdEffFinalPlots_HPSdbCorr_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPSdbCorr,
                                             binning[binVariable], fitVariables, outputFilePath, outputFileName_HPSdbCorr,
                                             expEff_label, measEff_label)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])
    logFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['logFileName'])
          
    # make plots for HPS combined isolation discriminators
    discriminators_HPScombined = [
        'tauDiscrHPScombLooseDBcorr',
        'tauDiscrHPScombMediumDBcorr',
        'tauDiscrHPScombTightDBcorr'
    ]
    outputFileName_HPScombined = 'makeTauIdEffFinalPlots_HPScombined_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPScombined,
                                             binning[binVariable], fitVariables, outputFilePath, outputFileName_HPScombined,
                                             expEff_label, measEff_label)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])
    logFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['logFileName'])
#--------------------------------------------------------------------------------

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

# done building config files, now build Makefile...
makeFileName = "Makefile_TauIdEffMeasAnalysis_%s" % "".join([ jobId, version ])
makeFile = open(makeFileName, "w")
makeFile.write("\n")
makeFile.write("all: %s %s\n" %
  (haddOutputFileName_stage3,
   make_MakeFile_vstring(outputFileNames_makeTauIdEffFinalPlots)))
makeFile.write("\techo 'Finished running TauIdEffMeasAnalysis.'\n")
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_FWLiteTauIdEffAnalyzer):
    makeFile.write("%s: %s\n" %
      (outputFileName,
       executable_FWLiteTauIdEffAnalyzer))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_FWLiteTauIdEffAnalyzer,
       configFileNames_FWLiteTauIdEffAnalyzer[i],
       logFileNames_FWLiteTauIdEffAnalyzer[i]))
makeFile.write("\n")
makeFile.write("%s: %s\n" %
  (haddOutputFileName_stage1,
   make_MakeFile_vstring(outputFileNames_FWLiteTauIdEffAnalyzer)))
makeFile.write("\t%s %s &> %s\n" %
  (executable_shell,
   haddShellFileName_stage1,
   haddLogFileName_stage1))
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_fitTauIdEff):
    makeFile.write("%s: %s %s\n" %
      (outputFileName,
       executable_fitTauIdEff,
       haddOutputFileName_stage1))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_fitTauIdEff,
       configFileNames_fitTauIdEff[i],
       logFileNames_fitTauIdEff[i]))
makeFile.write("\n")
makeFile.write("%s: %s\n" %
  (outputFileName_FWLiteTauIdEffPreselNumbers,
   executable_compTauIdEffPreselNumbers))
makeFile.write("\t%s %s &> %s\n" %
  (executable_compTauIdEffPreselNumbers,
   configFileName_FWLiteTauIdEffPreselNumbers,
   logFileName_FWLiteTauIdEffPreselNumbers))
makeFile.write("\n")
makeFile.write("%s:  %s %s\n" %
  (haddOutputFileName_stage2,
   make_MakeFile_vstring(outputFileNames_fitTauIdEff),
   outputFileName_FWLiteTauIdEffPreselNumbers))
makeFile.write("\t%s %s &> %s\n" %
  (executable_shell,
   haddShellFileName_stage2,
   haddLogFileName_stage2))
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_compTauIdEffFinalNumbers):
    makeFile.write("%s: %s %s\n" %
      (outputFileName,
       executable_compTauIdEffFinalNumbers,
       haddOutputFileName_stage2))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_compTauIdEffFinalNumbers,
       configFileNames_compTauIdEffFinalNumbers[i],
       logFileNames_compTauIdEffFinalNumbers[i]))
makeFile.write("\n")
makeFile.write("%s: %s %s\n" %
  (haddOutputFileName_stage3,
   make_MakeFile_vstring(outputFileNames_compTauIdEffFinalNumbers),
   make_MakeFile_vstring(outputFileNames_compTauIdEffFinalNumbers)))
makeFile.write("\t%s %s &> %s\n" %
  (executable_shell,
   haddShellFileName_stage3,
   haddLogFileName_stage3))
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_makeTauIdEffFinalPlots):
    makeFile.write("%s: %s %s\n" %
      (outputFileName,
       executable_makeTauIdEffFinalPlots,
       haddOutputFileName_stage3))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_makeTauIdEffFinalPlots,
       configFileNames_makeTauIdEffFinalPlots[i],
       logFileNames_makeTauIdEffFinalPlots[i]))
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(outputFileNames_FWLiteTauIdEffAnalyzer))
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(haddInputFileNames_stage1))
makeFile.write("\trm -f %s\n" % haddShellFileName_stage1)
makeFile.write("\trm -f %s\n" % haddOutputFileName_stage1)
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(outputFileNames_fitTauIdEff))
makeFile.write("\trm -f %s\n" % outputFileName_FWLiteTauIdEffPreselNumbers)
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(haddInputFileNames_stage2))
makeFile.write("\trm -f %s\n" % haddShellFileName_stage2)
makeFile.write("\trm -f %s\n" % haddOutputFileName_stage2)
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(outputFileNames_compTauIdEffFinalNumbers))
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(haddInputFileNames_stage3))
makeFile.write("\trm -f %s\n" % haddShellFileName_stage3)
makeFile.write("\trm -f %s\n" % haddOutputFileName_stage3)
makeFile.write("\trm -f %s\n" % make_MakeFile_vstring(outputFileNames_makeTauIdEffFinalPlots))
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)
