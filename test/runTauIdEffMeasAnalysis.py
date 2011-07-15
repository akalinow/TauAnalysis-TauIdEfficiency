#!/usr/bin/env python

from TauAnalysis.TauIdEfficiency.recoSampleDefinitionsTauIdEfficiency_7TeV_grid_cfi import recoSampleDefinitionsTauIdEfficiency_7TeV
from TauAnalysis.Configuration.userRegistry import getJobId
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import *

import os

channel = 'ZtoMuTau_tauIdEff'
#jobId = getJobId(channel)
jobId = '2011Jul06_mauroV4'

inputFilePath = '/data2/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/2011Jul06_mauro/V4/user/v/veelken/CMSSW_4_2_x/PATtuples/TauIdEffMeas/'
outputFilePath = '/data1/veelken/tmp/muonPtGt20/V4d/'

samplesToAnalyze = [
    # modify in case you want to submit jobs for some of the samples only...
    'data_SingleMu_Run2011A_PromptReco_v4',
    'data_SingleMu_Run2011A_May10ReReco_v1',
    'Ztautau_pythia',
    'Ztautau_embedded_part1',
    'Ztautau_embedded_part2',
    'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

# define sample name and jobId of Ztautau sample
# used to compute preselection efficiencies and purities in C1p and C1f regions
sampleZtautau = 'Ztautau_pythia'
jobId_noTauSel = '2011Jul06_mauroV4_noTauSel'

fitVariables = [
    'diTauVisMass',
    #'diTauVisMassFromJet' # CV: diTauVisMass always computed from PFJet momenta if using PAT-tuple workflow
]    

sysUncertainties = [
    "sysTauJetEn", # needed for diTauVisMass/diTauVisMassFromJet
    "sysJetEnUp"   # needed for diTauMt
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
    'tauEta' : {
        'tauEtaLt15' : {
            'min' :  0.0,
            'max' :  1.5
        },
        'tauEta15to19' : {
            'min' :  1.5,
            'max' :  1.9
        },
        'tauEta19to23' : {
            'min' :  1.9,
            'max' :  2.3
        },
        'xAxisTitle' : "#eta_{#tau}"
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

fitMethod = 'fitTauIdEff_wConstraints'
#fitMethod = 'fitTauIdEff'

executable_FWLiteTauIdEffAnalyzer = 'FWLiteTauIdEffAnalyzer'
executable_hadd = 'hadd'
executable_fitTauIdEff = fitMethod
executable_FWLiteTauIdEffPreselNumbers = 'FWLiteTauIdEffPreselNumbers'
executable_compTauIdEffFinalNumbers = 'compTauIdEffFinalNumbers'
executable_makeTauIdEffFinalPlots = 'makeTauIdEffFinalPlots'

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = recoSampleDefinitionsTauIdEfficiency_7TeV['RECO_SAMPLES'].keys()

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauIdEffAnalyzer macro
#
configFileNames_FWLiteTauIdEffAnalyzer = []
outputFileNames_FWLiteTauIdEffAnalyzer = []
for sampleToAnalyze in samplesToAnalyze:
    retVal_FWLiteTauIdEffAnalyzer = \
      buildConfigFile_FWLiteTauIdEffAnalyzer(sampleToAnalyze, jobId, inputFilePath, tauIds, sysUncertainties, outputFilePath,
                                             recoSampleDefinitionsTauIdEfficiency_7TeV)
    configFileNames_FWLiteTauIdEffAnalyzer.extend(retVal_FWLiteTauIdEffAnalyzer['configFileNames'])
    outputFileNames_FWLiteTauIdEffAnalyzer.extend(retVal_FWLiteTauIdEffAnalyzer['outputFileNames'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "harvest" histograms
# produced by FWLiteTauIdEffAnalyzer macro
#
haddShellFileName_stage1 = os.path.join(outputFilePath, 'mergeTauIdEffHistograms_stage1_%s.csh' % jobId)
haddOutputFileName_stage1 = os.path.join(outputFilePath, 'analyzeTauIdEffHistograms_all_%s.root' % jobId)
retVal_hadd_stage1 = \
  buildConfigFile_hadd(haddShellFileName_stage1, outputFileNames_FWLiteTauIdEffAnalyzer, haddOutputFileName_stage1)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running fitTauIdEff macro
#
configFileNames_fitTauIdEff = []
outputFileNames_fitTauIdEff = []
retVal_fitTauIdEff = \
  buildConfigFile_fitTauIdEff(fitMethod, jobId, '', haddOutputFileName_stage1, tauIds.keys(), fitVariables, outputFilePath, True)
configFileNames_fitTauIdEff.append(retVal_fitTauIdEff['configFileName'])
outputFileNames_fitTauIdEff.append(retVal_fitTauIdEff['outputFileName'])
for binVariable in binning.keys():
    for binName, binOptions in binning[binVariable].items():
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            retVal_fitTauIdEff = \
              buildConfigFile_fitTauIdEff(fitMethod, jobId, binName, haddOutputFileName_stage1, tauIds.keys(), fitVariables, outputFilePath, False)
            configFileNames_fitTauIdEff.append(retVal_fitTauIdEff['configFileName'])
            outputFileNames_fitTauIdEff.append(retVal_fitTauIdEff['outputFileName'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauIdEffPreselNumbers macro
#
retVal_FWLiteTauIdEffPreselNumbers = \
  buildConfigFile_FWLiteTauIdEffPreselNumbers(inputFilePath, sampleZtautau, jobId_noTauSel, tauIds, outputFilePath)
configFileName_FWLiteTauIdEffPreselNumbers = retVal_FWLiteTauIdEffPreselNumbers['configFileName']
outputFileName_FWLiteTauIdEffPreselNumbers = retVal_FWLiteTauIdEffPreselNumbers['outputFileName']
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "merge" fit results
# with preselection efficiencies and purities in regions C1p/C1f
#
haddShellFileName_stage2 = os.path.join(outputFilePath, 'mergeTauIdEffHistograms_stage2_%s.csh' % jobId)
haddInputFileNames_stage2 = []
haddInputFileNames_stage2.extend(outputFileNames_fitTauIdEff)
haddInputFileNames_stage2.append(outputFileName_FWLiteTauIdEffPreselNumbers)
haddOutputFileName_stage2 = os.path.join(outputFilePath, 'compTauIdEffFinalNumbers_input_%s.root' % jobId)
retVal_hadd_stage2 = \
  buildConfigFile_hadd(haddShellFileName_stage2, haddInputFileNames_stage2, haddOutputFileName_stage2)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build config files for running compTauIdEffFinalNumbers macro
#
configFileNames_compTauIdEffFinalNumbers = []
outputFileNames_compTauIdEffFinalNumbers = []
retVal_compTauIdEffFinalNumbers = \
  buildConfigFile_compTauIdEffFinalNumbers(haddOutputFileName_stage2, '', jobId, tauIds.keys(), fitVariables, outputFilePath)
configFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['configFileName'])
outputFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['outputFileName'])
for binVariable in binning.keys():
    for binName, binOptions in binning[binVariable].items():
        if isinstance(binOptions, dict) and binOptions.get('min') is not None and binOptions.get('max') is not None:
            retVal_compTauIdEffFinalNumbers = \
              buildConfigFile_compTauIdEffFinalNumbers(haddOutputFileName_stage2, binName, jobId, tauIds.keys(), fitVariables, outputFilePath)
            configFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['configFileName'])
            outputFileNames_compTauIdEffFinalNumbers.append(retVal_compTauIdEffFinalNumbers['outputFileName'])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to "merge" final tau id. efficiency numbers
# for different tauPt, tauEta,... bins
#
haddShellFileName_stage3 = os.path.join(outputFilePath, 'mergeTauIdEffHistograms_stage3_%s.csh' % jobId)
haddInputFileNames_stage3 = outputFileNames_compTauIdEffFinalNumbers
haddOutputFileName_stage3 = os.path.join(outputFilePath, 'compTauIdEffFinalNumbers_all_%s.root' % jobId)
retVal_hadd_stage3 = \
  buildConfigFile_hadd(haddShellFileName_stage3, haddInputFileNames_stage3, haddOutputFileName_stage3)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# make final tau id. efficiency plots as function of tauPt, tauEta,...
#
configFileNames_makeTauIdEffFinalPlots = []
outputFileNames_makeTauIdEffFinalPlots = []
for binVariable in binning.keys():

    # make plots for HPS isolation with no deltaBeta corrections applied
    discriminators_HPS = [
        'tauDiscrHPSloose',
        'tauDiscrHPSmedium',
        'tauDiscrHPStight'
    ]
    outputFileName_HPS = 'makeTauIdEffFinalPlots_HPS_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPS, binning[binVariable], fitVariables, outputFilePath, outputFileName_HPS)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])

    # make plots for HPS isolation with no applied deltaBeta corrections
    discriminators_HPSdbCorr = [
        'tauDiscrHPSlooseDBcorr',
        'tauDiscrHPSmediumDBcorr',
        'tauDiscrHPStightDBcorr'
    ]
    outputFileName_HPSdbCorr = 'makeTauIdEffFinalPlots_HPSdbCorr_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPSdbCorr, binning[binVariable], fitVariables, outputFilePath, outputFileName_HPSdbCorr)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])
          
    # make plots for HPS combined isolation discriminators
    discriminators_HPScombined = [
        'tauDiscrHPScombLooseDBcorr',
        'tauDiscrHPScombMediumDBcorr',
        'tauDiscrHPScombTightDBcorr'
    ]
    outputFileName_HPScombined = 'makeTauIdEffFinalPlots_HPScombined_%s.eps' % binVariable
    retVal_makeTauIdEffFinalPlots = \
      buildConfigFile_makeTauIdEffFinalPlots(haddOutputFileName_stage3, tauIds, discriminators_HPScombined, binning[binVariable], fitVariables, outputFilePath, outputFileName_HPScombined)
    configFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['configFileName'])
    outputFileNames_makeTauIdEffFinalPlots.append(retVal_makeTauIdEffFinalPlots['outputFileName'])
#--------------------------------------------------------------------------------

# done building config files, now build Makefile...
makeFileName = "Makefile_TauIdEffMeasAnalysis_%s" % jobId
makeFile = open(makeFileName, "w")
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_FWLiteTauIdEffAnalyzer):
    makeFile.write("%s:\n" % outputFileName)
    makeFile.write("\t%s %s\n" % (executable_FWLiteTauIdEffAnalyzer, configFileNames_FWLiteTauIdEffAnalyzer[i]))
makeFile.write("\n")
makeFile.write("%s: %s\n" % (haddOutputFileName_stage1,
                             make_inputFileNames_vstring(outputFileNames_FWLiteTauIdEffAnalyzer)))
makeFile.write("\tsource %s\n" % haddShellFileName_stage1)
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_fitTauIdEff):
    makeFile.write("%s: %s\n" % (outputFileName,
                                 haddOutputFileName_stage1))
    makeFile.write("\t%s %s\n" % (executable_fitTauIdEff, configFileNames_fitTauIdEff[i]))
makeFile.write("\n")
makeFile.write("%s:\n" % outputFileName_FWLiteTauIdEffPreselNumbers)
makeFile.write("\t%s %s\n" % (executable_FWLiteTauIdEffPreselNumbers, configFileName_FWLiteTauIdEffPreselNumbers))
makeFile.write("\n")
makeFile.write("%s: %s %s\n" % (haddOutputFileName_stage2,
                                make_inputFileNames_vstring(outputFileNames_fitTauIdEff),
                                outputFileName_FWLiteTauIdEffPreselNumbers))
makeFile.write("\tsource %s\n" % haddShellFileName_stage2)
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_compTauIdEffFinalNumbers):
    makeFile.write("%s: %s\n" % (outputFileName,
                                 haddOutputFileName_stage2))
    makeFile.write("\t%s %s\n" % (executable_compTauIdEffFinalNumbers, configFileNames_compTauIdEffFinalNumbers[i]))
makeFile.write("\n")
makeFile.write("%s: %s %s\n" % (haddOutputFileName_stage3,
                                make_inputFileNames_vstring(outputFileNames_compTauIdEffFinalNumbers),
                                make_inputFileNames_vstring(outputFileNames_compTauIdEffFinalNumbers)))
makeFile.write("\tsource %s\n" % haddShellFileName_stage3)
makeFile.write("\n")
for i, outputFileName in enumerate(outputFileNames_makeTauIdEffFinalPlots):
    makeFile.write("%s: %s\n" % (outputFileName,
                                 haddOutputFileName_stage3))
    makeFile.write("\t%s %s\n" % (executable_makeTauIdEffFinalPlots, configFileNames_makeTauIdEffFinalPlots[i]))
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)
