#!/usr/bin/env python

import copy
import os
import re

from TauAnalysis.Configuration.tools.jobtools import make_bsub_script
from TauAnalysis.Configuration.tools.harvesting import castor_source
from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import buildConfigFile_hadd

jobId = "2011Jul23" # jobId defined by submitTauIdEffMeasPATTupleProduction_noTauSel_lxbatch.py
version = 'V1c'

inputFilePath  = '/data2/veelken/CMSSW_4_2_x/PATtuples/TauPtRes/V1c/user/v/veelken/CMSSW_4_2_x/PATtuples/TauPtRes/V1c/'
outputFilePath = '/data1/veelken/tmp/tauPtResStudies/V1cA/'

sampleZtautau = 'Ztautau_powheg'

tauPtResOptions = [
    'Default',
    'Seed05Add00Strip05',
    'noEleTrackQcutsDefault',
    'noEleTrackQcutsSeed00Add00Strip05',
    'noEleTrackQcutsSeed00Add00Strip10',
    'noEleTrackQcutsSeed05Add00Strip05',
    'noEleTrackQcutsSeed05Add00Strip10',
    'noEleTrackQcutsSeed05Add00Strip15',
    'noEleTrackQcutsSeed05Add00Strip20',
    'noEleTrackQcutsSeed05Add00Strip25',
    'noEleTrackQcutsSeed05Add05Strip05',
    'noEleTrackQcutsSeed05Add05Strip10'
]

execDir = "%s/bin/%s/" % (os.environ['CMSSW_BASE'], os.environ['SCRAM_ARCH'])

executable_FWLiteTauPtResAnalyzer = execDir + 'FWLiteTauPtResAnalyzer'
executable_hadd = 'hadd -f'
executable_shell = '/bin/csh'

# check if inputFile is PAT-tuple and
# matches sampleZtautau, version
inputFileNames = os.listdir(inputFilePath)
inputFile_regex = \
  r"tauPtResPATtuple_%s_%s_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root" % (sampleZtautau, "".join([jobId, version]))
inputFileNames_Ztautau = []
for inputFileName in inputFileNames:
    inputFile_matcher = re.compile(inputFile_regex)
    if inputFile_matcher.match(inputFileName):
        inputFileNames_Ztautau.append(os.path.join(inputFilePath, inputFileName))

print " found %i input files." % len(inputFileNames_Ztautau)

if len(inputFileNames_Ztautau) == 0:
    raise ValueError("Sample %s has no input files --> skipping !!" % sampleZtautau)

from TauAnalysis.TauIdEfficiency.tools.buildConfigFilesTauIdEffAnalysis import make_inputFileNames_vstring
inputFileNames_string = make_inputFileNames_vstring(inputFileNames_Ztautau)

#--------------------------------------------------------------------------------
#
# build config files for running FWLiteTauPtResAnalyzer macro 
#
fileNames_FWLiteTauPtResAnalyzer = {}
for tauPtResOption in tauPtResOptions:
    outputFileName = 'analyzeTauPtResHistograms_%s_%s_%s.root' % (tauPtResOption, sampleZtautau, version)
    outputFileName_full = os.path.join(outputFilePath, outputFileName)
    directory = tauPtResOption
    srcTauJetCandidates = 'patTaus%s' % tauPtResOption
    selEventsFileName = 'selEvents_tauPtResLt0_5_%s_%s_%s.txt' % (tauPtResOption, sampleZtautau, version)
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

process.tauPtResAnalyzer = cms.PSet(
    directory = cms.string('%s'),

    srcTauJetCandidates = cms.InputTag('%s'),
    srcGenParticles = cms.InputTag('genParticles'),

    srcVertices = cms.InputTag('offlinePrimaryVerticesWithBS'),

    selEventsFileName = cms.string('%s')
)
""" % (inputFileNames_string, outputFileName_full, directory, srcTauJetCandidates, os.path.join(outputFilePath, selEventsFileName))

    configFileName = "analyzeTauPtResPATtuple_%s_%s_%s_cfg.py" % (tauPtResOption, sampleZtautau, version)
    configFileName_full = os.path.join(outputFilePath, configFileName)    
    configFile = open(configFileName_full, "w")
    configFile.write(config)
    configFile.close()

    logFileName = configFileName.replace('_cfg.py', '.log')
    logFileName_full = os.path.join(outputFilePath, logFileName)
    
    fileNames_FWLiteTauPtResAnalyzer[tauPtResOption] = {}
    fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['outputFileName'] = outputFileName
    fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['configFileName'] = configFileName_full
    fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['logFileName']    = logFileName_full
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# build shell script for running 'hadd' in order to collect histograms
# for all "seed", "add" and strip Pt thresholds into single .root file
#
haddInputFileNames = []
for tauPtResOption in tauPtResOptions:
    haddInputFileNames.append(os.path.join(outputFilePath,
                                           os.path.basename(fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['outputFileName'])))
haddShellFileName = os.path.join(outputFilePath, 'harvestTauPtResHistograms_%s_%s.csh' % (sampleZtautau, version))
haddOutputFileName = os.path.join(outputFilePath, 'analyzeTauPtResHistograms_all_%s_%s.root' % (sampleZtautau, version))
retVal_hadd = \
  buildConfigFile_hadd(executable_hadd, haddShellFileName, haddInputFileNames, haddOutputFileName)
haddLogFileName = retVal_hadd['logFileName']
#--------------------------------------------------------------------------------

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal

# done building config files, now build Makefile...
makeFileName = "Makefile_TauPtResAnalysis"
makeFile = open(makeFileName, "w")
makeFile.write("\n")
outputFileNames_makeTauPtResPlots = []
for tauPtResOption in tauPtResOptions:
    outputFileNames_makeTauPtResPlots.append(
      os.path.join(outputFilePath, fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['outputFileName']))
makeFile.write("all: %s %s\n" %
  (make_MakeFile_vstring(outputFileNames_makeTauPtResPlots),
   haddOutputFileName))
makeFile.write("\techo 'Finished running TauFakeRateAnalysis.'\n")
makeFile.write("\n")
for tauPtResOption in tauPtResOptions:
    makeFile.write("%s: %s\n" %
      (os.path.join(outputFilePath, fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['outputFileName']),
       executable_FWLiteTauPtResAnalyzer))
    makeFile.write("\t%s %s &> %s\n" %
      (executable_FWLiteTauPtResAnalyzer,
       fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['configFileName'],
       fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['logFileName']))
makeFile.write("\n")
makeFile.write("%s: %s\n" %
  (haddOutputFileName,
   make_MakeFile_vstring(haddInputFileNames)))
makeFile.write("\t%s %s &> %s\n" %
  (executable_shell,
   haddShellFileName,
   haddLogFileName))
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
for tauPtResOption in tauPtResOptions:
    makeFile.write("\trm -f %s\n" % os.path.join(outputFilePath, fileNames_FWLiteTauPtResAnalyzer[tauPtResOption]['outputFileName']))
makeFile.write("\trm -f %s\n" % haddOutputFileName)    
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)
