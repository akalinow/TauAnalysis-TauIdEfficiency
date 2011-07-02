import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTuple_template_cfg import *

#-------------------------------------------------------------------------------------------------------------------------
# Add few Ntuple variables to PAT-tuple
# (simple doubles)
#-------------------------------------------------------------------------------------------------------------------------
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTauIdEffMeas_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
    ntupleName = cms.string("tauIdEffNtuple"),
    sources = cms.PSet(
        # muon trigger efficiencies
        muonTriggerEff = process.tauIdEffMeas_muonTriggerEff.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        )                                    
    )
)

if isMC:
    # add reweighting factors to be applied to Monte Carlo simulated events
    # in order to match vertex multiplicity distribution in Data                                             
    setattr(process.ntupleProducer.sources, "vertexMultReweight", process.vertexMultReweight_template)
#-------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save PAT-tuple
#
process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep edmMergeableCounter_*_*_*',
            'keep *TriggerEvent_*_*_*',
            'keep patMuons_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
            'keep *_patGoodMuons_*_*',
            'keep patTaus_selectedPatPFTausShrinkingConePFRelIsoCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSPFRelIsoCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSpTaNCPFRelIsoCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsForTauIdEffCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsForTauIdEffCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative_*_*',
            'keep *Vertex*_*_*_*',
            'keep *_ak5PFJets_*_*',
            'keep *_ak5CaloJets_*_*',
            'keep *PFCandidate*_*_*_*'
        )               
    ),
    process.tauIdEffSampleEventSelection,
    fileName = cms.untracked.string("tauIdEffMeasPATtuple.root")      
)

if isMC:
    process.patTupleOutputModule.outputCommands.extend(
      cms.untracked.vstring(
            "keep *_*_*vtxMultReweight*_*",
            'keep *GenJet*_*_*_*',
            'keep *recoGenMET_*_*_*',
            'keep *_genParticles_*_*',
            'keep *_generator_*_*',
            'keep *_*GenJets_*_*',
            'keep *_addPileupInfo_*_*'
        )
    )
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patDefaultSequence
   + process.producePatTupleTauIdEffMeasSpecific
   #+ process.printEventContent
   #+ process.printGenParticleList
   + process.ntupleProducer
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.patTupleOutputModule)

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    process.muonPFTauFixedConeSkimPath,
    process.muonPFTauShrinkingConeSkimPath,
    process.muonPFTauHPSskimPath,
    process.muonPFTauHPSpTaNCskimPath,
    process.o
)
