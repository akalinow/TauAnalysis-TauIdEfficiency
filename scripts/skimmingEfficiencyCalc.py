#!/usr/bin/env python
from __future__ import with_statement
try:
    from lxml import etree
except ImportError:
    #we deal with old pythons here
    import xml.etree.cElementTree as etree
"""
skimmingEfficiencyCalc.py

calculate the skimming efficiency from a given crabdir (you have to downloat the reports first)

Author: Matthias Edelhoff (Aachen)

Takes as input:

   o Input crab dir(s) to extract total number of skimmed events and total number of input events from

"""
def getFrameworkReports(crabDir):
    import os
    from glob import glob

    result = []
    for xmlPath in glob(os.path.join(crabDir,"res/crab_fjr_[0-9]*.xml")):
        with open(xmlPath,"r") as xmlFile:
            skip = False
            try:
                report = etree.parse(xmlFile)
            except SyntaxError, err:
                skip = True
                print "WARNING SyntaxError in '%s': %s. Skipping this file"%(xmlPath,err)
            #for some reason this exitCode does not coincide with crab -status
            #exitCodeNode = report.find("/ExitCode")
            #if not int(exitCodeNode.get("Value")) == 0:
            #    skip = True
            #    print "WARNING Exit code = %s in '%s'. skipping this file"%(int(exitCodeNode.get("Value")), xmlPath)
            if not skip:
                result.append(report)
    return result

def calcEfficiency(reports):
    skimmedEvents = 0
    readEvents = 0
    for report in reports:
        totalEventsNodes = report.findall("/File/TotalEvents")
        if len(totalEventsNodes) != 1:
            raise StandardError, "to exactly one totalEvents node in a jobReport"
        skimmedEvents += float(totalEventsNodes[0].text)

        for node in report.findall("/InputFile/EventsRead"):
            readEvents += int(node.text)
    return skimmedEvents / readEvents
            

def main(argv=None):
    from optparse import OptionParser
    import sys
    if argv == None:
        argv = sys.argv[1:]
    parser = OptionParser()
    parser.add_option("-c", "--crabDirs", dest="crabDirs", action="append", default=[],
                      help="crab directory which contains the framework reports (you have to do crab -get first!)")
    (opts, args) = parser.parse_args(argv)
    for crabDir in opts.crabDirs:
        reports = getFrameworkReports(crabDir)
        if not len(reports) > 0:
            raise StandardError, "no valid framework report found for '%s'"%crabDir
        eff = calcEfficiency(reports)
        print "Efficiency for %s: %s"%(crabDir, eff)

                     
                      
if __name__ == '__main__':
    main()
