import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.ObjValEDNtupleProducer_cfi import ntupleProducer

def buildPatTauNtuple(process, ntupleName, collections, variables, appendTo=None):
    # check if we are adding to an existing ntuple
    output = appendTo
    if output is None:
        output = ntupleProducer.clone()


