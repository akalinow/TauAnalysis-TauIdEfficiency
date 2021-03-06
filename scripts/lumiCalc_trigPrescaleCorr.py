#!/usr/bin/env python
from __future__ import with_statement
import sys, os
import json
import csv
from TauAnalysis.TauIdEfficiency.tools.prescales import RunRegistry

'''

lumiCalc.py

Compute the effective luminosity of a dataset, given the selected runs and
prescales for those runs.

Author: Matthias Edelhoff (Aachen)
Contributors: Evan K. Friis (UC Davis)

$Id: lumiCalc_trigPrescaleCorr.py,v 1.10 2010/06/28 15:19:21 edelhoff Exp $

Takes as input: 
    
    o CSV file mapping lumi sections to their corresponding integrated
    luminosities

    o A json file containing the mask that maps runs to good lumisection ranges
    corresponding to a dataset. (i.e. the file produced crab -report)

    o A text file containing the list of ROOT files associated with this
    dataset

'''

def saveUpdate(dictA, dictB):
    for key in dictB:
            if key in dictA and not dictB[key] == dictA[key]:
                raise StandardError, \
                        "there are overlaping and divergent data for %s"%(key)
    dictA.update( dictB )
    return dictA

def cleanNTupleFiles( path ):
    """
    clean crab retries: take only last retry located in given path and return list of .root files
    path can be:
    o rfio:/castor... for castor use
    o a list of strings
    """
    from TauAnalysis.TauIdEfficiency.tools.castor import nsls
    import re
    rawFiles = []
    if path.startswith("rfio:/castor"):
        rawFiles = [ i.replace(path[5:], "") for i in nsls(path[5:])]
    elif "__iter__" in dir(path):
        rawFiles = [i for i in path]
        
    expr = re.compile("(.*)_([0-9]*)_([0-9]*)_(?:...).root")
    fileMap = {}
    for fileName in rawFiles:
        skip = False
        try:
            groups = expr.search(fileName).groups()
        except AttributeError:
            skip = True
            print "WARNING: bad fromating for '%s'. skipping"%fileName
        if not skip:
            if not groups[0] in fileMap:
                fileMap[groups[0]] = {}
            if not int(groups[1]) in fileMap[groups[0]]:
                fileMap[groups[0]][int(groups[1])] = int(groups[2])
            else:
                if fileMap[groups[0]][int(groups[1])] < int(groups[2]):
                    fileMap[groups[0]][int(groups[1])] = int(groups[2])

    result =[]
    for name in fileMap:
        for number in fileMap[name]:
                result.append("%s_%s_%s.root"%(name,number,fileMap[name][number]))
    return result
            

def readMask(path):
    " Load a LumiSel (CRAB-style) JSON file "
    parse_result = None
    with open(path, "r") as maskFile:
        parse_result = json.load(maskFile)
    # convert all strings to ints
    result = {}
    for parse_key, parse_value in parse_result.iteritems():
        result[int(parse_key)] = parse_value
    return result

def readNTuples(path):
    """ Read a list of input Ntuple root files and their trigger path

    The file must have the trigger path on the first line, followed by the
    associated root files (one per line)
    
    Returns (HLTPath, fileList)
    
    """
    with open(path, "r") as fileList:
        result = [i.strip() for i in fileList.read().splitlines()]
    return (result[0], result[1:])

def readCSV(pathToFile, skipHeader = 0):
    ''' Yields dictionaries for each row of a CSV file

    The field names (keys) for each column are taken from the header line
    (after skip header) of the file.
    '''
    with open(pathToFile, 'r') as csvFile:
        # Skip header lines
        for skipper in range(skipHeader):
            csvFile.readline()
        reader = csv.DictReader(csvFile)
        # Strip any whitespace from the field names
        reader.fieldnames[:] = [field.strip() for field in reader.fieldnames]
        for row in reader:
            yield row

def readLumiCSV(lumiCSVPath, skipHead = 4):
    " Populate (run/lumisec) -> int. luminosity map "
    result = {}
    for row in readCSV(lumiCSVPath, skipHead):
        run_lumi_key = (int(row['Run']), int(row['LS']))
        # Double check this has't been inserted already
        if run_lumi_key in result:
            raise KeyError, "duplicated Run/Lumi key in file %s"%lumiCSVPath
        result[run_lumi_key] = {
            'hfLumi': float(row['HF Lumi']),
            'vtxLumi': float(row['Vtx Lumi'])
        }
    return result

def calcLumi(lumiTable, mask, run_registry, path, formula = "vtxLumi", maxRelDelta=1.):
    #does the lumisec list inclusive? think so...
    from math import fabs
    result = 0.0
    for run in mask:
        lumi_i = 0.0
        for lumisecRange in mask[run]:
            for lumisec in range(lumisecRange[0], lumisecRange[1]+1):
                key = (run,lumisec)
                if key in lumiTable:
                    lumisecLumi = eval(formula, lumiTable[key])
                    relDelta = (
                        fabs(lumiTable[key]["vtxLumi"]-lumiTable[key]["hfLumi"]) 
                        /lumisecLumi) if lumisecLumi != 0.0 else 0.0
                    if relDelta > maxRelDelta:
                        print "WARNING: large (%.2f) diffenrece in measured lumis detected in %s: vtx = %.2f hf = %.2f"%(relDelta,key,lumiTable[key]["vtxLumi"],lumiTable[key]["hfLumi"])
                    lumi_i+= lumisecLumi
        prescale = 1
        if path:
            prescale = run_registry.get_prescale(run, path)
        result += lumi_i/prescale
    return result

def main(argv=None):
    from optparse import OptionParser
    if argv == None:
        argv = sys.argv[1:]
    parser = OptionParser()
    parser.add_option("-d", "--datasetName", dest="datasetNames", action="append", default=[], 
                      help="appendable list of dataset names. Expects DATASETNAME_JSON.txt for used mask and DATASETNAME.list for a list if paths to nTuple Files")
    parser.add_option("-m", "--mcDatasetName", dest="mcDatasetNames", action="append", default=[], 
                      help="appendable list of dataset names. Name ust be key of mcLumiMap.json")
    parser.add_option("-M", "--mcLumiMap", dest="mcLumiMapPath", default="$CMSSW_BASE/src/TauAnalysis/TauIdEfficiency/test/commissioning/mcLumiMap.json", 
                      help="path to output file")


    parser.add_option("-l", "--lumi", dest="lumiCSVs", action="append", default=[], 
                      help="appendable list of luminosity by lumisection csv files to use")
    parser.add_option("-o", "--output", dest="outPath", default="lumiMap.json", 
                      help="path to output file")
    (opts, args) = parser.parse_args(argv)

    if opts.datasetNames == [] and opts.mcDatasetNames == []: raise StandardError, "it makes no sense not to specify at least on dataset name"
    if opts.lumiCSVs == []: opts.lumiCSVs = ["lumi_by_LS.csv"]
    
    if os.path.exists(opts.outPath): raise StandardError, "Output file '%s' exists. Please move it out of the way!"%opts.outPath
    
    opts.mcLumiMapPath = os.path.expandvars(opts.mcLumiMapPath)
    
    lumiMap = {}
    if len(opts.datasetNames) > 0:
        lumiTable = {}
        for lumiCSV in opts.lumiCSVs:
            newLumis = readLumiCSV(lumiCSV)
            lumiTable = saveUpdate(lumiTable, newLumis)

        # Load masks for the datasets
        masks = {}
        for datasetName in opts.datasetNames:
            masks[datasetName] = readMask("%s_JSON.txt"%datasetName)

        def loop_over_runs():
            for dataset, mask in masks.iteritems():
                for run in mask.keys():
                    yield run

        min_run_number = min(loop_over_runs())
        max_run_number = max(loop_over_runs())
        print min_run_number
        print max_run_number
        
        # Build our run registry to access prescale info from DB
        run_registry = RunRegistry(firstRun=min_run_number, lastRun=max_run_number, 
                                       groupName="Collisions%")

        for datasetName in opts.datasetNames:
            print "Building lumi map for %s dataset" % datasetName
            mask = masks[datasetName]
            triggerPath, nTuples = readNTuples("%s.list"%datasetName)
            intLumi = calcLumi(lumiTable, mask, run_registry, triggerPath)
            nTuples = cleanNTupleFiles(nTuples)
            lumiMap[datasetName] = {'files': nTuples, 'int_lumi': intLumi}

        
    if len(opts.mcDatasetNames) > 0:
        for mcDatasetName in opts.mcDatasetNames:
            lumiMap[mcDatasetName] = {}
            with open(opts.mcLumiMapPath, "r") as mcLumiMapFile:
                mcLumiMap = json.load(mcLumiMapFile)
                lumiMap[mcDatasetName].update(mcLumiMap[mcDatasetName])
                lumiMap[mcDatasetName]["files"] = cleanNTupleFiles(lumiMap[mcDatasetName]["directory"])                

    # Write output to json file
    with open(opts.outPath, 'w') as outputFile:
        json.dump(lumiMap, outputFile, sort_keys=True, indent=4)
        
if __name__ == '__main__':
    main()
