import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasPATTupleNoTauSel")

#--------------------------------------------------------------------------------
# define configuration parameter default values

##isMC = True # use for MC
isMC = False # use for Data
isEmbedded = False # use for everything except for Ztautau samples produced via MCEmbedding technique
#isEmbedded = True # use for Ztautau samples produced via MCEmbedding technique
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "HLT" # use for Summer'11 MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__isEmbedded = #isEmbedded#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#
#--------------------------------------------------------------------------------

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific import produceTauIdEffMeasPATTuple
produceTauIdEffMeasPATTuple(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection)

for processAttrName in dir(process):
    processAttr = getattr(process, processAttrName)
    if isinstance(processAttr, cms.EDFilter):
        if processAttr.type_() == "PATPFTauSelectorForTauIdEff":
            print "--> Disabling cut %s" % processAttrName
            setattr(processAttr, "produceAll", cms.bool(True))

if hasattr(process.patTupleOutputModule, "SelectEvents"):
    delattr(process.patTupleOutputModule, "SelectEvents")

processDumpFile = open('produceTauIdEffMeasPATTuple_noTauSel.dump' , 'w')
print >> processDumpFile, process.dumpPython()



