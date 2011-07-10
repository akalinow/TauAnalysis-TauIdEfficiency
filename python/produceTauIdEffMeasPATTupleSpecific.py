import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTuple_base import produceTauIdEffMeasPATTuple_base
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

def produceTauIdEffMeasPATTuple(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection, applyZrecoilCorrection):

    # produce PAT objects common between PAT-tuple and Ntuple production
    produceTauIdEffMeasPATTuple_base(process, isEmbedded, isMC, HLTprocessName, pfCandidateCollection, applyZrecoilCorrection)

    #--------------------------------------------------------------------------------
    #
    # CV: add (new) SVfit algorithm;
    #     for speed reasons, run SVfit in 'fit' mode only
    #    (using combination of PS + MET likelihoods + logM regularization term
    #     to reconstruct mass of tau lepton pair, as described in CMS AN-11-165)
    #
    muTauPairProducers = [
        #"allMuPFTauHPSpairsForTauIdEff",
        #"allMuPFTauHPSpairsForTauIdEffSysJetEnUp",
        #"allMuPFTauHPSpairsForTauIdEffSysJetEnDown",
        #"allMuPFTauHPSpairsForTauIdEffSysTauJetEnUp",
        #"allMuPFTauHPSpairsForTauIdEffSysTauJetEnDown",
        #"allMuPFTauHPSpTaNCpairsForTauIdEff",
        #"allMuPFTauHPSpTaNCpairsForTauIdEffSysJetEnUp",
        #"allMuPFTauHPSpTaNCpairsForTauIdEffSysJetEnDown",
        #"allMuPFTauHPSpTaNCpairsForTauIdEffSysTauJetEnUp",
        #"allMuPFTauHPSpTaNCpairsForTauIdEffSysTauJetEnDown"
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

        muTauPairProducerModule.doPFMEtSign = cms.bool(True)
    
        muTauPairProducerModule.pfMEtSign = cms.PSet(
            srcPFJets = cms.InputTag('ak5PFJets'),
            srcPFCandidates = cms.InputTag('particleFlow'),
            resolution = process.METSignificance_params,
            dRoverlapPFJet = cms.double(0.3),
            dRoverlapPFCandidate = cms.double(0.1)
        )
    #--------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------------
    #
    # add few Ntuple variables to PAT-tuple
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
    # save PAT-tuple
    #
    process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
        cms.PSet(
            outputCommands = cms.untracked.vstring(
                'drop *',
                'keep EventAux_*_*_*',
                'keep LumiSummary_*_*_*',                       
                'keep edmMergeableCounter_*_*_*',
                'keep *_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
                'keep *_patGoodMuons_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEffCumulative_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEffSysJetEnUpCumulative_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEffSysJetEnDownCumulative_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnUpCumulative_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEffSysTauJetEnDownCumulative_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEffCumulative_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysJetEnUpCumulative_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysJetEnDownCumulative_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnUpCumulative_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEffSysTauJetEnDownCumulative_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEffCumulative_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEffSysJetEnUpCumulative_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEffSysJetEnDownCumulative_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEffSysTauJetEnUpCumulative_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEffSysTauJetEnDownCumulative_*_*',
                'keep *_selectedPatPFTausHPSpTaNCforTauIdEffCumulative_*_*',
                'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysJetEnUpCumulative_*_*',
                'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysJetEnDownCumulative_*_*',                                            
                'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysTauJetEnUpCumulative_*_*',
                'keep *_selectedPatPFTausHPSpTaNCforTauIdEffSysTauJetEnDownCumulative_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffCumulative_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysJetEnUpCumulative_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysJetEnDownCumulative_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffCumulative_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysJetEnUpCumulative_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysJetEnDownCumulative_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEffCumulative_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauEnUpCumulative_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauEnDownCumulative_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
                'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffCumulative_*_*',
                'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysJetEnUpCumulative_*_*',
                'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysJetEnDownCumulative_*_*',
                'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnUpCumulative_*_*',
                'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEffSysTauJetEnDownCumulative_*_*',
                'keep *_offlinePrimaryVertices_*_*',
                'keep *_offlinePrimaryVerticesWithBS_*_*',
                'keep *_selectedPrimaryVertexHighestPtTrackSum_*_*',                                         
                'keep *_*_*muonHLTeff*_*',
                'keep *_patPFMETs_*_*',
                'keep *_smearedMETpFTauHPSsysJetEnUp_*_*',
                'keep *_smearedMETpFTauHPSsysJetEnDown_*_*',
                'keep *_smearedMETpFTauHPSsysTauJetEnUp_*_*',
                'keep *_smearedMETpFTauHPSsysTauJetEnDown_*_*',
                'keep *_smearedMETpFTauHPSpTaNCsysJetEnUp_*_*',
                'keep *_smearedMETpFTauHPSpTaNCsysJetEnDown_*_*',
                'keep *_smearedMETpFTauHPSpTaNCsysTauJetEnUp_*_*',
                'keep *_smearedMETpFTauHPSpTaNCsysTauJetEnDown_*_*',
                # CV: additional collections needed to run nSVfit algorithm                                                        
                #'keep recoTracks_generalTracks_*_*',                                             
                #'keep *_ak5PFJets_*_*',
                #'keep *_particleFlow_*_*'
            )               
        ),
        process.tauIdEffSampleEventSelection,
        fileName = cms.untracked.string("tauIdEffMeasPATtuple.root")      
    )

    process.patTupleOutputModule.outputCommands.extend(patTriggerEventContent)

    if isMC:
        process.patTupleOutputModule.outputCommands.extend(
          cms.untracked.vstring(
                'keep *_*_*vtxMultReweight*_*',
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

    #--------------------------------------------------------------------------------  
    #
    # CV: keep Z --> tau+ tau- --> muon + tau-jet events
    #     passing Pt and eta cuts on generator level
    #    (for studying preselection efficiencies)
    #
    if isMC:
        process.load('PhysicsTools.JetMCAlgos.TauGenJets_cfi')
        process.load('TauAnalysis.GenSimTools.gen_decaysFromZs_cfi')

        process.genMuonWithinAccFilter = cms.EDFilter("PATCandViewCountFilter",
            src = cms.InputTag('genMuonsFromZtautauDecaysWithinAcceptance'),
            minNumber = cms.uint32(1),
            maxNumber = cms.uint32(1000)
        )

        process.genHadTauWithinAccFilter = cms.EDFilter("PATCandViewCountFilter",
            src = cms.InputTag('genHadronsFromZtautauDecaysWithinAcceptance'),
            minNumber = cms.uint32(1),
            maxNumber = cms.uint32(1000)
        )

        process.genZtoMuTauWithinAccSkimPath = cms.Path(
            process.tauGenJets
           + process.produceGenDecayProductsFromZs
           + process.genMuonWithinAccFilter + process.genHadTauWithinAccFilter
        )
    
        extSkimPaths = process.patTupleOutputModule.SelectEvents.SelectEvents.value()
        extSkimPaths.append('genZtoMuTauWithinAccSkimPath')
        process.patTupleOutputModule.SelectEvents.SelectEvents = cms.vstring(extSkimPaths)
    #-------------------------------------------------------------------------------- 

    # define order in which different paths are run
    if isMC:
        process.schedule = cms.Schedule(
            process.p,
            ##process.muonPFTauFixedConeSkimPath,
            ##process.muonPFTauShrinkingConeSkimPath,
            process.muonPFTauHPSskimPath,
            process.muonPFTauHPSpTaNCskimPath,
            process.genZtoMuTauWithinAccSkimPath,
            process.o
        )
    else:
        process.schedule = cms.Schedule(
            process.p,
            ##process.muonPFTauFixedConeSkimPath,
            ##process.muonPFTauShrinkingConeSkimPath,
            process.muonPFTauHPSskimPath,
            process.muonPFTauHPSpTaNCskimPath,
            process.o
        )
