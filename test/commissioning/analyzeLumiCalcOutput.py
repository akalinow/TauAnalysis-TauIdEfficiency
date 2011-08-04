#!/usr/bin/env python

import os
import re
import sys

print("<analyzeLumiCalcOutput>:")

if len(sys.argv) != 2:
    raise ValueError("Usage: findBadCastorFiles.py fileName")

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

recLumi_ref = None
effLumiSum = 0.
hltPaths = []
    
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

        if recLumi_ref is not None:
            if recLumi != recLumi_ref:
                raise ValueError("Recorded luminosity is ambiguous !!")
        recLumi_ref = recLumi
        
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
        hltPaths.append(hltPath)

print("hltPaths: %s" % hltPaths)
print(" recLumi = %f" % recLumi_ref)
print(" effLumi = %f" % effLumiSum)
print("--> average Prescale factor = %f" % (recLumi_ref/effLumiSum))
