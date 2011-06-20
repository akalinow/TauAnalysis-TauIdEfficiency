import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasNtuple_template_cfg import *

for processAttrName in dir(process):
    processAttr = getattr(process, processAttrName)
    if isinstance(processAttr, cms.EDFilter):
        if processAttr.type_() == "PATTauSelector":
            print "--> Disabling cut %s" % processAttrName
            setattr(processAttr, "cut", cms.string('pt > 10. & abs(eta) < 2.5'))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/m/mverzett/tagprobe/skims/TauIdEffMeas_Harvested_Jun06/skim_Ztautau_powheg_chunk_2_8025.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)


process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
                                              cms.PSet(
                                                  outputCommands = cms.untracked.vstring(
                                                      'drop *',
                                                      "keep *_*_*vtxMultReweight*_*",
                                                      'keep edmMergeableCounter_*_*_*',
                                                      'keep *TriggerEvent_*_*_*',
                                                      'keep patMuons_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
                                                      'keep *_patMuonsStandAlone_*_*',
                                                      'keep patTaus_selectedPatPFTausShrinkingConePFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSPFRelIsoCumulative_*_*',
                                                      'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsForTauIdEffCumulative_*_*',
                                                      'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
                                                      'keep *_ak5PFJets_*_*',
                                                      'keep *_ak5CaloJets_*_*',
                                                      'keep *Vertex*_*_*_*',
                                                      'keep *GenJet*_*_*_*',
#                                                      'keep *recoGenMET_*_*_*',
                                                      'keep *_genParticles_*_*',
                                                      'keep *_generator_*_*', 
#                                                      'keep *_*GenJets_*_*',
                                                      'keep *PFCandidate*_*_*_*',
                                                      'keep *_addPileupInfo_*_*'
                                                      )                               
                                                  ),
                                              process.tauIdEffSampleEventSelection,
                                              fileName = cms.untracked.string("tauIdEffMeasPATtuple_NoPreselection.root")
                                              )

processDumpFile = open('produceTauIdEffMeasNtuple_noTauPresel.dump' , 'w')
print >> processDumpFile, process.dumpPython()
