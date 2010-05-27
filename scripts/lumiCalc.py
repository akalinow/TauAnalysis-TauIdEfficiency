#!/usr/bin/env python
import sys, os

def saveUpdate(dictA, dictB):
    for key in dictB:
            if key in dictA and not dictB[key] == dictA[key]:
                raise StandardError, "there are overlaping and divergent data for %s"%( key)
    dictA.update( dictB )
    return dictA

def readMask(path):
    maskFile = open(path,"r")
    result = eval(maskFile.read())
    maskFile.close()
    return result

def readNTuples(path):
    maskFile = open(path,"r")
    result = [i.strip() for i in maskFile.read().splitlines()]
    maskFile.close()
    return result

def readLumiCSV(path, skipHead = 5):
    #todo: add some validity checks
    row= { "run": 0, "lumisec":1, "hfLumi": 4, "vtxLumi":5} 
    lumiCSVFile=open(path,"r")
    lumiCSV = lumiCSVFile.read().splitlines()
    lumiCSVFile.close()
    result = {}
    for line in lumiCSV[skipHead:]:
        run = int(line.split(",")[row["run"]])
        lumisec = int(line.split(",")[row["lumisec"]])
        hfLumi = float(line.split(",")[row["hfLumi"]])
        vtxLumi = float(line.split(",")[row["vtxLumi"]])
        result["%s:%s"%(run,lumisec)] = {"hfLumi": hfLumi, "vtxLumi":vtxLumi}
    return result

def readPrescaleCSV(path, skipHead = 0):
    #todo: add some validity checks
    #todo: avoid duplication from readLumiCSV
    row = {"run":0, "prescale":1}
    prescaleCSVFile=open(path,"r")
    prescaleCSV = prescaleCSVFile.read().splitlines()
    prescaleCSVFile.close()
    result = {}
    for line in prescaleCSV[skipHead:]:
        run = int(line.split(",")[row["run"]])
        prescale = float(line.split(",")[row["prescale"]])
        result[run]= {"prescale":prescale}
    return result

def checkPrescales(prescaleTable, mask):
    lastGiven = 1.0
    result = {}
    for run in [ int(i) for i in mask.keys()]:
        if not run in prescaleTable:
            given = raw_input("missing prescale for run %s (empty line = %s):"%(run, lastGiven))
            if given == "": given = lastGiven
            else: lastGiven = given
            result[run]= {"prescale":float(given)}
    return result

def calcLumi(lumiTable, prescaleTable, mask, formula = "vtxLumi", maxRelDelta=1.):
    #does the lumisec list inclusive? think so...
    from math import fabs
    result = 0.0
    for run in mask:
        lumi_i = 0.0
        for lumisecRange in mask[run]:
            for lumisec in range(lumisecRange[0], lumisecRange[1]+1):
                key = "%s:%s"%(run,lumisec)
                if key in lumiTable:
                    lumisecLumi = eval(formula, lumiTable[key])
                    relDelta = fabs(lumiTable[key]["vtxLumi"]-lumiTable[key]["hfLumi"])/lumisecLumi if lumisecLumi != 0.0 else 0.0
                    if relDelta > maxRelDelta:
                        print "WARNING: large (%.2f) diffenrece in measured lumis detected in %s: vtx = %.2f hf = %.2f"%(relDelta,key,lumiTable[key]["vtxLumi"],lumiTable[key]["hfLumi"])
                    lumi_i+= lumisecLumi
        result += lumi_i/prescaleTable[int(run)]["prescale"]
    return result

def main(argv=None):
    from optparse import OptionParser
    from pprint import pformat
    if argv == None:
        argv = sys.argv[1:]
    parser = OptionParser()
    parser.add_option("-d", "--datasetName", dest="datasetNames", action="append", default=[], 
                      help="appendable list of dataset names. Expects DATASETNAME_JSON.txt for used mask and DATASETNAME.list for a list if paths to nTuple Files")
    parser.add_option("-l", "--lumi", dest="lumiCSVs", action="append", default=[], 
                      help="appendable list of luminosity by lumisection csv files to use")
    parser.add_option("-p", "--prescale", dest="prescaleCSVs", action="append", default=[], 
                      help="appendable list of prescales by runnumber csv files to use")
    parser.add_option("-o", "--output", dest="outPath", default="lumiMap.json", 
                      help="path to output file")
    (opts, args) = parser.parse_args(argv)

    if opts.datasetNames == []: raise StandardError, "it makes no sense not to specify at least on dataset name"
    if opts.lumiCSVs == []: opts.lumiCSVs = ["lumi_by_LS.csv"]
    if opts.prescaleCSVs == []: opts.prescaleCSVs = ["prescales_by_RUN.csv"]
    
    if os.path.exists(opts.outPath): raise StandardError, "Output file '%s' exists. Please move it out of the way!"%opts.outPath

    lumiTable = {}
    for lumiCSV in opts.lumiCSVs:
        newLumis = readLumiCSV(lumiCSV)
        lumiTable = saveUpdate(lumiTable, newLumis)

    prescaleTable = {}
    for prescaleCSV in opts.prescaleCSVs:
        newPrescales = readPrescaleCSV(prescaleCSV)
        prescaleTable = saveUpdate(prescaleTable, newPrescales)

    prescaleTableAdditions = {}
    lumiMap = {}
    for datasetName in opts.datasetNames:
        mask = readMask("%s_JSON.txt"%datasetName)
        nTuples = readNTuples("%s.list"%datasetName)
        additions = checkPrescales(prescaleTable, mask)
        saveUpdate(prescaleTable, additions)
        prescaleTableAdditions.update(additions)
        intLumi = calcLumi(lumiTable, prescaleTable, mask)
        lumiMap[datasetName] = [nTuples, intLumi]

    
    outFile = open(opts.outPath,"w")
    outFile.write(pformat(lumiMap))
    outFile.close()
    
    if not prescaleTableAdditions == {}:
        prescaleAddionString = ""
        for run in prescaleTableAdditions:
            prescaleAddionString +="%s, %f\n"%(run, prescaleTableAdditions[run]["prescale"])
        prescaleAddionFile = open(opts.outPath.replace(".json","_prescaleAdditions.csv"),"w")
        prescaleAddionFile.write(prescaleAddionString)
        prescaleAddionFile.close()

        
if __name__ == '__main__':
    main()
