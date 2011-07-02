import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific_template_cfg import *

process.patTupleOutputModule.outputCommands = cms.untracked.vstring(
    'drop *',
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
    'keep *_*_*muonHLTeff*_*'
)

from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent
process.patTupleOutputModule.outputCommands.extend(patTriggerEventContent)

if isMC:
    process.patTupleOutputModule.outputCommands.extend(
        cms.untracked.vstring(
            'keep *_*_*vtxMultReweight*_*',
            'keep *_addPileupInfo_*_*'
        )
    )
