import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasPATTuple")

#--------------------------------------------------------------------------------
# define configuration parameter default values

isMC = True # use for MC
##isMC = False # use for Data
isEmbedded = False # use for everything except for Ztautau samples produced via MCEmbedding technique
#isEmbedded = True # use for Ztautau samples produced via MCEmbedding technique
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "HLT" # use for Summer'11 MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
#applyZrecoilCorrection = False
applyZrecoilCorrection = True
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__isEmbedded = #isEmbedded#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#__applyZrecoilCorrection = #applyZrecoilCorrection#
#
#--------------------------------------------------------------------------------

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific import produceTauIdEffMeasPATTuple
produceTauIdEffMeasPATTuple(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection, applyZrecoilCorrection)

processDumpFile = open('produceTauIdEffMeasPATTuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
