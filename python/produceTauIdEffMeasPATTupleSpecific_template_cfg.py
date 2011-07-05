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
            'keep EventAux_*_*_*',
            'keep edmMergeableCounter_*_*_*',
            'keep patMuons_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
            'keep patMuons_patGoodMuons_*_*',
            #'keep patTaus_selectedPatPFTausFixedConeForTauIdEffCumulative_*_*',
            #'keep patTaus_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep patTaus_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep patTaus_selectedPatPFTausShrinkingConeForTauIdEffCumulative_*_*',
            #'keep patTaus_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep patTaus_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSForTauIdEffCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSpTaNCForTauIdEffCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSpTaNCForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep patTaus_selectedPatPFTausHPSpTaNCForTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauFixedConePairsDzForTauIdEffCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsDzForTauIdEffCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsDzForTauIdEffCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *CompositePtrCandidateT1T2MEt*_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *Vertex*_*_*_*',
            'keep *_*_*muonHLTeff*_*',
            'keep *_patPFMETs_*_*'                                                           
        )               
    ),
    process.tauIdEffSampleEventSelection,
    fileName = cms.untracked.string("tauIdEffMeasPATtuple.root")      
)

from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.patTupleOutputModule.outputCommands.extend(patTriggerEventContent)

if isMC:
    process.patTupleOutputModule.outputCommands.extend(
      cms.untracked.vstring(
            "keep *_*_*vtxMultReweight*_*",
            'keep *_addPileupInfo_*_*',
            'keep *_genParticles_*_*',
            'keep *_tauGenJets_*_*',
            'keep *_tauGenJetsSelectorAllHadrons_*_*',
            'keep *_genMetTrue_*_*'
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
