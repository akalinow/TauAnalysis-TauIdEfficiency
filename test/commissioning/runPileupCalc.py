#!/usr/bin/env python

import os
import re
import shlex
import subprocess
import sys

#--------------------------------------------------------------------------------
# Make histogram of expected pile-up interactions
# for given JSON file and HLT path
# (needed for pile-up reweighting)
#
# Author: Christian Veelken (LLR)
# 
# Examples:
# 
# ./runPileupCalc.py JSONfile lumiJSONfile HLTpath
#
#--------------------------------------------------------------------------------

print("<runPileupCalc.py>:")

if not len(sys.argv) == 4:
    raise ValueError("Usage: runPileupCalc.py JSONfile lumiJSONfile HLTpath")

inputFileName_JSON = sys.argv[1]
print " inputFileName_JSON = %s" % inputFileName_JSON
inputFileName_lumiJSON = sys.argv[2]
print " inputFileName_lumiJSON = %s" % inputFileName_lumiJSON

HLT_path = sys.argv[3]
print " HLT_path = %s" % HLT_path

runRange_regex = \
  r"[a-zA-Z0-9_-]*_(?P<firstRun>[0-9]+)-(?P<lastRun>[0-9]+)[a-zA-Z0-9_-]*.txt"
runRange_matcher = re.compile(runRange_regex)
runRange_match = runRange_matcher.match(os.path.basename(inputFileName_JSON))
firstRun = int(runRange_match.group('firstRun'))
print " firstRun = %i " % firstRun
lastRun = int(runRange_match.group('lastRun'))
print " lastRun = %i " % lastRun

outputFileName_csv = "my" + os.path.basename(inputFileName_JSON).replace(".txt", "_%s.csv" % HLT_path)
outputFileName_JSONforHLTpath =  "my" + os.path.basename(inputFileName_JSON).replace(".txt", "_%s.txt" % HLT_path)
outputFileName_histogram3d = "expPUpoissonMean_runs%ito%i_%s.root" % (firstRun, lastRun, HLT_path)
outputFileName_histogram3d = outputFileName_histogram3d.replace("_HLT_", "_")
outputFileName_histogram1d = "expPUpoissonDist_runs%ito%i_%s.root" % (firstRun, lastRun, HLT_path)
outputFileName_histogram1d = outputFileName_histogram1d.replace("_HLT_", "_")

def runCommand(commandLine):
    sys.stdout.write("%s\n" % commandLine)
    args = shlex.split(commandLine)
    retVal = subprocess.Popen(args, stdout = subprocess.PIPE)
    retVal.wait()
    return retVal

commandLine_lumiCalc2 = \
  "lumiCalc2.py lumibyls -i %s --hltpath %s -o %s" % \
    (inputFileName_JSON, HLT_path, outputFileName_csv)
runCommand(commandLine_lumiCalc2)

commandLine_pileupReCalc_HLTpath = \
  "pileupReCalc_HLTpaths.py -i %s --inputLumiJSON %s -o %s" % \
    (outputFileName_csv, inputFileName_lumiJSON, outputFileName_JSONforHLTpath)
runCommand(commandLine_pileupReCalc_HLTpath)

commandLine_pileupCalc3d = \
  "pileupCalc.py -i %s --inputLumiJSON %s --calcMode true --minBiasXsec 69400 --maxPileupBin 50 --numPileupBins 500 %s" % \
    (inputFileName_JSON, outputFileName_JSONforHLTpath, outputFileName_histogram3d)
runCommand(commandLine_pileupCalc3d)

commandLine_pileupCalc1d = \
  "pileupCalc.py -i %s --inputLumiJSON %s --calcMode observed --minBiasXsec 69400 --maxPileupBin 60 --numPileupBins 60 %s" % \
    (inputFileName_JSON, outputFileName_JSONforHLTpath, outputFileName_histogram1d)
runCommand(commandLine_pileupCalc1d)
