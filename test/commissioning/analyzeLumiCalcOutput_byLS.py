#!/usr/bin/env python

import os
import re
import sys
from collections import defaultdict

print("<analyzeLumiCalcOutput_byLS>:")

if len(sys.argv) != 3:
    raise ValueError("Usage: analyzeLumiCalcOutput_byLS.py fileName_lumiCalc fileName_lumiContext")

inputFileName_lumiCalc = sys.argv[1]
inputFileName_lumiContext = sys.argv[2]

units = {
    '/mb' : 1.0e-12,
    '/ub' : 1.0e-9,
    '/nb' : 1.0e-6,
    '/pb' : 1.0e-3,
    '/fb' : 1.0
}

def readLumiContextOutput(inputFileName):
    
    retVal = {}
    
    header_regexp = "\|\s+Run\s+\|\s+LS\s+\|\s+dfrac\s+\|\s+\(bitname,count,presc\)\s+\|\s*"
    header_matcher = re.compile(header_regexp)
    commentLine_regexp = "[-]+\s*"
    commentLine_matcher = re.compile(commentLine_regexp)
    row_regexp = "\|\s*(?P<run>[0-9]*)\s*\|\s*(?P<ls>[0-9:]*)\s*\|\s*[0-9.]*\s*\|\s*(?P<l1Seed>[(a-zA-Z0-9_),\s*]*)\s*\|\s*"
    row_matcher = re.compile(row_regexp)

    inputFile = open(inputFileName, "r")
    isHeaderReached  = False
    current_run      = None
    current_ls_start = None
    current_ls_end   = None
    for line in inputFile.readlines():

        if header_matcher.match(line):
            isHeaderReached = True
 
        # CV: skip all lines before header
        if not isHeaderReached:
            continue

        # CV: skip all comment lines
        if commentLine_matcher.match(line):
            continue
        
        row_match = row_matcher.match(line)
        if row_match:
            run_string = row_match.group('run')
            if run_string != "":
                current_run = int(run_string)
            ls_string = row_match.group('ls')
            if ls_string != "":
                ls_items = ls_string.split(":")
                if len(ls_items) == 1:
                    current_ls_start = int(ls_items[0])
                    current_ls_end   = int(ls_items[0])
                elif len(ls_items) == 2:
                    current_ls_start = int(ls_items[0])
                    current_ls_end   = int(ls_items[1])
                else:
                    raise ValueError("Failed to parse line = '%s' !!" % line)
            l1Seed_string = row_match.group('l1Seed')
            startPos = l1Seed_string.find("(")
            while startPos != -1:
                endPos = l1Seed_string.find(")", startPos + 1)                
                if endPos == -1:
                    raise ValueError("Failed to parse line = '%s' !!" % line)
                l1Seed_items = l1Seed_string[startPos + 1:endPos].split(",")
                if len(l1Seed_items) != 3:
                    raise ValueError("Failed to parse line = '%s' !!" % line)
                l1Seed_name = l1Seed_items[0]
                l1Seed_counts = l1Seed_items[1]
                l1Seed_prescale = l1Seed_items[2]

                for current_ls in range(current_ls_start, current_ls_end + 1):
                    if current_run not in retVal.keys():
                        retVal[current_run] = {}
                    if current_ls not in retVal[current_run].keys():
                        retVal[current_run][current_ls] = {}
                    retVal[current_run][current_ls][l1Seed_name] = int(l1Seed_prescale)

                startPos = l1Seed_string.find("(", endPos + 1)

    inputFile.close()

    return retVal

def readLumiCalcOutput_byLS(inputFileName):

    retVal = {}

    header_regexp = "\|\s+Run:Fill\s+\|\s+LS\s+\|\s+HLTpath\s+\|\s+L1bit\s+\|\s+HLTpresc\s+\|\s+L1presc\s+\|\s+Recorded\((?P<unit>[/a-z]+)\)\s+\|\s+Effective\(/nb\)\s+\|\s*"
    header_matcher = re.compile(header_regexp)
    commentLine_regexp = "[-]+\s*"
    commentLine_matcher = re.compile(commentLine_regexp)
    row_regexp = "\|\s*(?P<run>[0-9:]+)\s*\|\s*(?P<ls>[0-9:]+)\s*\|\s*[a-zA-Z0-9_]*\s*\|\s*[a-zA-Z0-9_]*\s*\|\s*[0-9]+\s*\|\s*[0-9]+\s*\|\s*(?P<luminosity>[0-9.]+)\s*\|\s*[0-9.]+\s*\|\s*"
    row_matcher = re.compile(row_regexp)

    inputFile = open(inputFileName, "r")
    isHeaderReached  = False
    unit             = None
    unit_value       = None
    current_run      = None
    current_ls_start = None
    current_ls_end   = None
    for line in inputFile.readlines():

        header_match = header_matcher.match(line)
        if header_match:
            unit = header_match.group("unit")
            if not unit in units.keys():
                raise ValueError("Invalid unit = '%s' defined in line = %s !!" % (unit, line))
            unit_value = units[unit]
            isHeaderReached = True
 
        # CV: skip all lines before header
        if not isHeaderReached:
            continue

        # CV: skip all comment lines
        if commentLine_matcher.match(line):
            continue
        
        row_match = row_matcher.match(line)
        if row_match:
            run_string = row_match.group('run')
            run_items = run_string.split(":")
            current_run = int(run_items[0])
            ls_string = row_match.group('ls')
            current_ls = None
            if ls_string != "":
                ls_items = ls_string.split(":")
                if len(ls_items) == 1:
                    current_ls = int(ls_items[0])
                elif len(ls_items) == 2:
                    if ls_items[0] != ls_items[1]:
                        raise ValueError("Current version of analyzeLumiCalcOutput_byLS script requires unique lumi-sections !!")
                    current_ls = int(ls_items[0])
                else:
                    raise ValueError("Failed to parse line = '%s' !!" % line)

            if current_run not in retVal.keys():
                retVal[current_run] = {}
            if current_ls not in retVal[current_run].keys():
                retVal[current_run][current_ls] = {}
            luminosity_string = row_match.group('luminosity')
            retVal[current_run][current_ls] = float(luminosity_string)*unit_value

    inputFile.close()

    return retVal

luminosities = readLumiCalcOutput_byLS(inputFileName_lumiCalc)
prescales = readLumiContextOutput(inputFileName_lumiContext)

lumi_per_l1Seed = {}

for run in luminosities.keys():
    for ls in luminosities[run].keys():
        if luminosities[run][ls] <= 0.:
            continue
        
        if not (run in prescales.keys() and ls in prescales[run].keys()):
            raise ValueError("Failed to read L1 prescales for run = %i, ls = %i !!" % (run, ls))

        for l1Seed in prescales[run][ls]:
            if not l1Seed in lumi_per_l1Seed.keys():
                lumi_per_l1Seed[l1Seed] = 0.
            ##if prescales[run][ls][l1Seed] > 0:
            ##    lumi_per_l1Seed[l1Seed] = lumi_per_l1Seed[l1Seed] + luminosities[run][ls]/float(prescales[run][ls][l1Seed])
            if prescales[run][ls][l1Seed] == 1:
                lumi_per_l1Seed[l1Seed] = lumi_per_l1Seed[l1Seed] + luminosities[run][ls]/float(prescales[run][ls][l1Seed])
                
print("l1Seeds:")
for l1Seed in lumi_per_l1Seed.keys():
    print "%s: %f %s" % (l1Seed, lumi_per_l1Seed[l1Seed]*units['/fb'], "/fb")
