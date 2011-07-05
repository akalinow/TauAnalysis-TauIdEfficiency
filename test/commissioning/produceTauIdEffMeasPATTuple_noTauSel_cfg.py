import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific_template_cfg import *

for processAttrName in dir(process):
    processAttr = getattr(process, processAttrName)
    if isinstance(processAttr, cms.EDFilter):
        if processAttr.type_() == "PATPFTauSelectorForTauIdEff":
            print "--> Disabling cut %s" % processAttrName
            setattr(processAttr, "produceAll", cms.bool(True))

if hasattr(process.patTupleOutputModule, "SelectEvents"):
    delattr(process.patTupleOutputModule, "SelectEvents")

process.source.fileNames = cms.untracked.vstring(
    'file:/data1/veelken/CMSSW_4_2_x/skims/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_AOD.root'
)   
process.patTupleOutputModule.fileName = cms.untracked.string(
    '/data1/veelken/CMSSW_4_2_x/PATtuples/tauIdEffMeasPATtupleGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2.root'
)

processDumpFile = open('produceTauIdEffMeasPATTuple_noTauPresel.dump' , 'w')
print >> processDumpFile, process.dumpPython()



