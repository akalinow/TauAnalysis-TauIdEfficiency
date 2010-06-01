#!/usr/bin/env python

'''

prescales

Provides interface to retrieve HLT prescales for a given run and path from
Config databases.

Original author: Emmanuelle Perez

Modifications by Evan Friis

'''

import xmlrpclib
import xml.dom.minidom
import subprocess
import re
import json
from sys import stderr, exit

class RunRegistry(object):
    ''' Interface to run registry '''
    def __init__(self, firstRun, lastRun, groupName, 
                 runUrl="http://pccmsdqm04.cern.ch/runregistry/xmlrpc"):
        self.firstRun = firstRun
        self.lastRun = lastRun
        self.server = xmlrpclib.ServerProxy(runUrl)
        stderr.write(
            "Querying run registry for range [%d, %d], group name like %s ...\n"
            % (firstRun, lastRun, groupName))
        # Retrieve information about HLT from run registry
        xml_data =xml.dom.minidom.parseString(
            self.server.DataExporter.export(
            'RUN', 'GLOBAL', 'xml_datasets', 
            "{runNumber} >= %d AND {runNumber} <= %d AND {groupName} like '%s' AND {datasetName} = '/Global/Online/ALL'"  
            % (firstRun, lastRun, groupName)))
        xml_runs = xml_data.documentElement.getElementsByTagName("RUN_DATASET")
        # Associate information about HLT to each RUN
        self.runs = {}
        for xml_run in xml_runs:
            self.runs[
                xml_run.getElementsByTagName("RUN_NUMBER")[0].firstChild.nodeValue
            ] = xml_run.getElementsByTagName("RUN_HLTKEY")[0].firstChild.nodeValue

    def get_prescale(self, run, path):
        # get appropriate HLT config
        matcher = re.compile(r'\|\s*%s\s*\|\s*(?P<prescale>\d+)\s*\|' % path.strip())
        config = self.runs[run]
        cmd = [
            'edmConfigFromDB',
            '--orcoff',
            '--format',
            'summary.ascii',
            '--paths',
            path,
            '--configName',
            config ]
        # Run query
        query = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        query.wait()
        matching_prescale = list(
            match.group('prescale') for match in 
            (matcher.match(line) for line in query.stdout) 
            if match and match.group('prescale')
        )
        # No paths matched
        if not matching_prescale:
            return 0
        # Return last one (equivalent to tail -1 in perez' impl)
        else: return int(matching_prescale[-1])


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="usage: %prog [options] Trigger_Path")
    parser.add_option("--firstRun",  dest="firstRun",  help="first run", type="int", metavar="RUN", default="1")
    parser.add_option("--lastRun",   dest="lastRun",   help="last run",  type="int", metavar="RUN", default="9999999")
    parser.add_option("--groupName", dest="groupName", help="select runs of name like NAME", metavar="NAME", default="Collisions%")
    parser.add_option("--rrurl",     dest="rrurl",     help="run registry xmlrpc url", metavar="URL", default="http://pccmsdqm04.cern.ch/runregistry/xmlrpc")
    parser.add_option("--jsonOut",   dest="jsonOut",   help="dump prescales in JSON format on FILE", metavar="FILE")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_usage()
        exit(2)
    path = args[0]

    stderr.write("Building run registry\n")
    run_registry = RunRegistry(options.firstRun, options.lastRun, options.groupName,
                               options.rrurl)

    json_out = {}
    for run in run_registry.runs:
        prescale = run_registry.get_prescale(run, path)
        json_out[run] = prescale
        print run, path, prescale

    if options.jsonOut:
        with open(options.jsonOut, 'w') as jsonFile:
            jsonFile.write(json.dumps(json_out))

