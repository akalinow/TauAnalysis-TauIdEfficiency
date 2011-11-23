import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTuple_base import produceTauIdEffMeasPATTuple_base
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

def produceTauIdEffMeasPATTuple(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection):

    # produce PAT objects common between PAT-tuple and Ntuple production
    produceTauIdEffMeasPATTuple_base(process, isMC, isEmbedded, HLTprocessName, pfCandidateCollection)

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

        pfJetCollection = 'selectedPatJetsAK5PFAntiOverlapWithMuonsVeto'
        if isMC:
            pfJetCollection = 'selectedPatJetsAK5PFsmearedAntiOverlapWithMuonsVeto'
        muTauPairProducerModule.pfMEtSign = cms.PSet(
            srcPFJets = cms.InputTag(pfJetCollection),
            srcPFCandidates = cms.InputTag('particleFlow'),
            resolution = process.METSignificance_params,
            dRoverlapPFJet = cms.double(0.3),
            dRoverlapPFCandidate = cms.double(0.1)
        )
    #--------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------------------------------------------------
    #
    # add Data to Monte-Carlo correction factors
    # (simple doubles)
    #
    if isMC:
        # add Data/MC correction factor for muon trigger efficiency
        process.muonTriggerEfficiencyCorrection = cms.EDProducer("PATMuonEfficiencyCorrectionProducer",
           inputFileName = cms.FileInPath("TauAnalysis/TauIdEfficiency/data/singleMuHLTeff.root"),
           lutName = cms.string('hEff'),
           parametrization = cms.PSet(
               src = cms.VInputTag(composeModuleName(['selectedPatMuonsForTauIdEffPFRelIso', "cumulative"])),
               x = cms.string("eta()"),
               y = cms.string("pt()"),
               z = cms.string("0.")
           ),
           noObjectSubstituteValue = cms.double(0.) # weight returned in case all 'src' collections do not contain any entries
        )
        process.producePatTupleTauIdEffMeasSpecific += process.muonTriggerEfficiencyCorrection
        
        # add reweighting factors (rho-Neutral correction)
        # to correct for Data/MC differences in modeling out-out-time pile-up
        process.load("TauAnalysis.RecoTools.vertexMultiplicityVsRhoPFNeutralReweight_cfi")
        process.producePatTupleTauIdEffMeasSpecific += process.produceVertexMultiplicityVsRhoPFNeutralReweights
    #-------------------------------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------
    # add event counter for Mauro's "self baby-sitting" technology
    process.processedEventsPATtupleProduction = cms.EDProducer("EventCountProducer")
    process.eventCounterPath = cms.Path(process.processedEventsPATtupleProduction)
    #--------------------------------------------------------------------------------

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
                'keep *_gtDigis_*_*',
                'keep *_hltL1GtObjectMap_*_*',
                'keep *_TriggerResults_*_*',                                            
                'keep *_selectedPatMuonsForTauIdEffPFRelIsoCumulative_*_*',
                'keep *_selectedPatMuonsForTauIdEffZmumuHypotheses_*_*',
                'keep *_selectedDiMuPairForTauIdEffZmumuHypotheses_*_*',
                'keep *_patJetsAK5PFnotOverlappingWithLeptonsForMEtUncertainty*_*_*',
                #'keep *_selectedPatPFTausFixedConeForTauIdEff*_*_*',
                #'keep *_selectedPatPFTausShrinkingConeForTauIdEff*_*_*',
                'keep *_selectedPatPFTausHPSforTauIdEff*_*_*',
                #'keep *_selectedPatPFTausHPSpTaNCforTauIdEff*_*_*',
                #'keep *_selectedMuPFTauFixedConePairsDzForTauIdEff*_*_*',
                #'keep *_selectedMuPFTauShrinkingConePairsDzForTauIdEff*_*_*',
                'keep *_selectedMuPFTauHPSpairsDzForTauIdEff*_*_*',
                #'keep *_selectedMuPFTauHPSpTaNCpairsDzForTauIdEff*_*_*',
                'keep *_offlinePrimaryVertices_*_*',
                'keep *_offlinePrimaryVerticesWithBS_*_*',
                'keep *_selectedPrimaryVertexHighestPtTrackSum_*_*',
                # keep CaloMEt objects                                                            
                'keep *_patMETs*_*_*',
                'keep *_patCaloMet*_*_*',                                  
                # keep PFMEt objects                                               
                'keep *_patPFMet*_*_*',
                'keep *_patType1CorrectedPFMet*_*_*',
                'keep *_patType1p2CorrectedPFMet*_*_*',
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
                'keep *_smearedPatJetsAK5PF*_*_*',
                'keep *_muonTriggerEfficiencyCorrection_*_*',
                'keep *_vertexMultiplicityReweight3dRunA_*_*',
                'keep *_vertexMultiplicityReweight3dRunB_*_*',
                'keep *_vertexMultiplicityVsRhoPFNeutralReweight_*_*',
                'keep *_addPileupInfo_*_*',
                'keep *_genParticles_*_*',
                'keep *_tauGenJets_*_*',
                'keep *_tauGenJetsSelectorAllHadrons_*_*',
                'keep *_genMetTrue_*_*',
                # CV: additional collections needed to run nSVfit algorithm  
                #'keep *_smearedAK5PFJets_*_*'
            )
        )
    #--------------------------------------------------------------------------------

    process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer") 

    process.p = cms.Path(
        process.prePatProductionSequence
       + process.producePatTupleTauIdEffMeasSpecific
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
            #process.muonPFTauFixedConeSkimPath,
            #process.muonPFTauShrinkingConeSkimPath,
            process.muonPFTauHPSskimPath,
            #process.muonPFTauHPSpTaNCskimPath,
            process.genZtoMuTauWithinAccSkimPath,
            process.o
        )
    else:
        process.schedule = cms.Schedule(
            process.p,
            #process.muonPFTauFixedConeSkimPath,
            #process.muonPFTauShrinkingConeSkimPath,
            process.muonPFTauHPSskimPath,
            #process.muonPFTauHPSpTaNCskimPath,
            process.o
        )
