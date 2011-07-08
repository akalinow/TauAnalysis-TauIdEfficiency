import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasNtuple")

#--------------------------------------------------------------------------------
# define configuration parameter default values

##isMC = True # use for MC
isMC = False # use for Data
##HLTprocessName = "HLT" # use for 2011 Data
HLTprocessName = "HLT" # use for Summer'11 MC
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
applyZrecoilCorrection = False
#applyZrecoilCorrection = True
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#__applyZrecoilCorrection = #applyZrecoilCorrection#
#
#--------------------------------------------------------------------------------

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasNtupleSpecific import produceTauIdEffMeasNtuple
produceTauIdEffMeasNtuple(process, isMC, HLTprocessName, pfCandidateCollection, applyZrecoilCorrection)

processDumpFile = open('produceTauIdEffMeasNtuple.dump' , 'w')
print >> processDumpFile, process.dumpPython()
