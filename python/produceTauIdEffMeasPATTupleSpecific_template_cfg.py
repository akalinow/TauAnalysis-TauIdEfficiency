import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTuple_template_cfg import *

#--------------------------------------------------------------------------------
#
# CV: add (new) SVfit algorithm;
#     for speed reasons, run SVfit in 'fit' mode only
#    (using combination of PS + MET likelihoods + logM regularization term
#     to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
#
muTauPairProducers = [
    ##"allMuPFTauHPSpairsForTauIdEff",
    ##"allMuPFTauHPSpTaNCpairsForTauIdEff"
]
for muTauPairProducer in muTauPairProducers:
    muTauPairProducerModule = getattr(process, muTauPairProducer)

    muTauPairProducerModule.doSVreco = cms.bool(True)
    
    muTauPairProducerModule.nSVfit = cms.PSet()
    muTauPairProducerModule.nSVfit.psKine_MEt_logM_fit = cms.PSet()
    muTauPairProducerModule.nSVfit.psKine_MEt_logM_fit.config = copy.deepcopy(process.nSVfitConfig_template)
    muTauPairProducerModule.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg1 = cms.PSet(
        src = muTauPairProducerModule.srcLeg1,
        likelihoodFunctions = cms.VPSet(process.nSVfitMuonLikelihoodPhaseSpace),
        builder = process.nSVfitTauToMuBuilder
    )
    muTauPairProducerModule.nSVfit.psKine_MEt_logM_fit.config.event.resonances.A.daughters.leg2 = cms.PSet(
        src = muTauPairProducerModule.srcLeg2,
        likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace),
        builder = process.nSVfitTauToHadBuilder
    )
    muTauPairProducerModule.nSVfit.psKine_MEt_logM_fit.algorithm = cms.PSet(
        pluginName = cms.string("nSVfitAlgorithmByLikelihoodMaximization"),
        pluginType = cms.string("NSVfitAlgorithmByLikelihoodMaximization"),
        minimizer  = cms.vstring("Minuit2", "Migrad"),
        maxObjFunctionCalls = cms.uint32(5000),
        verbosity = cms.int32(0)
    )
    
nSVfitProducers = [
    ##[ 'nSVfitProducerByLikelihoodMaximizationPFTauHPS',      'selectedPatPFTausHPSforTauIdEffCumulative'      ],
    ##[ 'nSVfitProducerByLikelihoodMaximizationPFTauHPSpTaNC', 'selectedPatPFTausHPSpTaNCforTauIdEffCumulative' ]
]
for nSVfitProducer in nSVfitProducers:
    nSVfitProducerModule = copy.deepcopy(process.nSVfitProducerByLikelihoodMaximization)
    setattr(process, nSVfitProducer[0], nSVfitProducerModule)

    nSVfitProducerModule.config.event.resonances.A.daughters.leg1 = cms.PSet(
        src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
        likelihoodFunctions = cms.VPSet(process.nSVfitMuonLikelihoodPhaseSpace),
        builder = process.nSVfitTauToMuBuilder
    )
    nSVfitProducerModule.config.event.resonances.A.daughters.leg2 = cms.PSet(
        src = cms.InputTag(nSVfitProducer[1]),
        likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace),
        builder = process.nSVfitTauToHadBuilder
    )
    nSVfitProducerModule.algorithm = cms.PSet(
        pluginName = cms.string("nSVfitAlgorithmByLikelihoodMaximization"),
        pluginType = cms.string("NSVfitAlgorithmByLikelihoodMaximization"),
        minimizer  = cms.vstring("Minuit2", "Migrad"),
        maxObjFunctionCalls = cms.uint32(5000),
        verbosity = cms.int32(0)
    )

    process.producePatTupleTauIdEffMeasSpecific += nSVfitProducerModule
#--------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------
#
# Add few Ntuple variables to PAT-tuple
# (simple doubles)
#
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
            'keep *_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
            'keep *_patGoodMuons_*_*',
            #'keep *_selectedPatPFTausFixedConeForTauIdEffCumulative_*_*',
            #'keep *_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep *_selectedPatPFTausShrinkingConeForTauIdEffCumulative_*_*',
            #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *_selectedPatPFTausHPSforTauIdEffCumulative_*_*',
            'keep *_selectedPatPFTausHPSforTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *_selectedPatPFTausHPSforTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *_selectedPatPFTausHPSpTaNCforTauIdEffCumulative_*_*',
            'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffCumulative_*_*',
            #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffCumulative_*_*',
            #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *_selectedMuPFTauHPSpairsDzForTauIdEffCumulative_*_*',
            'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffCumulative_*_*',
            'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
            'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
            'keep *_offlinePrimaryVertices_*_*',
            'keep *_offlinePrimaryVerticesWithBS_*_*',
            'keep *_selectedPrimaryVertexHighestPtTrackSum_*_*',                                         
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

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer") 

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patDefaultSequence
   + process.producePatTupleTauIdEffMeasSpecific
   + process.ntupleProducer
   ##+ process.printEventContent
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.patTupleOutputModule)

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    ##process.muonPFTauFixedConeSkimPath,
    ##process.muonPFTauShrinkingConeSkimPath,
    process.muonPFTauHPSskimPath,
    process.muonPFTauHPSpTaNCskimPath,
    process.o
)
