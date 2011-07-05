import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific_template_cfg import *

processDumpFile = open('produceTauIdEffMeasPATTuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
