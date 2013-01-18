#!/usr/bin/env python

import os
import re
import sys
from collections import defaultdict

print("<analyzeLumiCalcOutput>:")

if len(sys.argv) != 2:
    raise ValueError("Usage: analyzeLumiCalcOutput.py fileName")

inputFileName = sys.argv[1]

summaryRowMatcher_regexp = \
   r"\|\s*(?P<hltPath>[a-zA-Z0-9_]+)\s*\|\s*[0-9]+\s*" \
  + "\|\s*(?P<recValue>[0-9.]+)\((?P<recUnits>[/a-z]+)\)\s*" \
  + "\|\s*(?P<effValue>[0-9.]+)\((?P<effUnits>[/a-z]+)\)\s*|" 
summaryRowMatcher = re.compile(summaryRowMatcher_regexp)

units = {
    '/mb' : 1.0e-12,
    '/ub' : 1.0e-9,
    '/nb' : 1.0e-6,
    '/pb' : 1.0e-3,
    '/fb' : 1.0
}

recLumiSum = 0.
effLumiSum = 0.
hltPaths = defaultdict(float)
    
inputFile = open(inputFileName, "r")
isSummaryReached = False
for line in inputFile.readlines():
    if line.find("Total") != -1:
        isSummaryReached = True
    if not isSummaryReached:
        continue
    summaryRowMatch = summaryRowMatcher.match(line)
    if summaryRowMatch:
        recValue = summaryRowMatch.group('recValue')
        recUnits = summaryRowMatch.group('recUnits')
        if not (recValue and recUnits):
            continue
        if recUnits not in units:
            raise ValueError("Undefined units = %s !!" % recUnits)
        recLumi = float(recValue)*units[recUnits]
        #print("recLumi = %f" % recLumi)
        recLumiSum += recLumi
        
        effValue = summaryRowMatch.group('effValue')
        effUnits = summaryRowMatch.group('effUnits')
        if not (effValue and effUnits):
            continue
        if effUnits not in units:
            raise ValueError("Undefined units = %s !!" % effUnits)
        effLumi = float(effValue)*units[effUnits]
        #print("effLumi = %f" % effLumi)
        effLumiSum += effLumi

        hltPath = summaryRowMatch.group('hltPath')
        hltPaths[hltPath] += effLumi

print(" recLumi = %f /fb" % recLumiSum)
print(" effLumi = %f /fb" % effLumiSum)
print("--> average Prescale factor = %f" % (recLumiSum/effLumiSum))
print("hltPaths:")
for hltPathName, effLumiSum_hltPath in hltPaths.items():
    print(" %s: %f %s" % (hltPathName, effLumiSum_hltPath*units['/fb'], "/fb")) 
