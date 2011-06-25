import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasNtuple_template_cfg import *

processDumpFile = open('produceTauIdEffMeasNtuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
