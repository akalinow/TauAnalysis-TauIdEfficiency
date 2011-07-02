import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific_template_cfg import *

##for processAttrName in dir(process):
##    processAttr = getattr(process, processAttrName)
##    if isinstance(processAttr, cms.EDFilter):
##        if processAttr.type_() == "PATTauSelector":
##            print "--> Disabling cut %s" % processAttrName
##            setattr(processAttr, "cut", cms.string('pt > 10. & abs(eta) < 2.5'))

for processAttrName in dir(process):
    processAttr = getattr(process, processAttrName)
    if isinstance(processAttr, cms.EDFilter):
        if processAttr.type_() == "PATPFTauSelectorForTauIdEff":
            print "--> Disabling cut %s" % processAttrName
            setattr(processAttr, "produceAll", cms.bool(True))

if hasattr(process.patTupleOutputModule, "SelectEvents"):
    delattr(process.patTupleOutputModule, "SelectEvents")

processDumpFile = open('produceTauIdEffMeasPATTuple_noTauPresel.dump' , 'w')
print >> processDumpFile, process.dumpPython()



