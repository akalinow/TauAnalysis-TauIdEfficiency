import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
from TauAnalysis.CandidateTools.tools.objSelConfigurator import objSelConfigurator
from TauAnalysis.CandidateTools.tools.objProdConfigurator import objProdConfigurator
import PhysicsTools.PatAlgos.tools.helpers as configtools

def buildSequenceTauIdEffMeasSpecific(process,
                                      patMuonCollectionName = "selectedPatMuonsForTauIdEffPFRelIsoCumulative",
                                      tauIdAlgorithmName = None, patTauCollectionName = "patTaus", applyTauJEC = False,
                                      savePatTaus = None,
                                      patMEtCollectionName = "patType1CorrectedPFMet",
                                      isMC = False, isEmbedded = False,
                                      runSVfit = False):

    #print("<buildSequenceTauIdEffMeasSpecific>:")
    #print(" patTauCollectionName = %s" % patTauCollectionName)
    #print(" isMC = %s" % isMC)

    # check that tauIdAlgorithmName is defined, non-null and composed of two parts
    if tauIdAlgorithmName is None or len(tauIdAlgorithmName) != 2:
        raise ValueError("Undefined of invalid 'tauIdAlgorithmName' Parameter !!")

    patTauProducerName = "".join(patTauCollectionName)

    sequenceName = composeModuleName(["sequenceTauIdEffMeasSpecific", "".join(tauIdAlgorithmName)])
    sequence = cms.Sequence()
    setattr(process, sequenceName, sequence)
    
    #--------------------------------------------------------------------------------
    # define loose Tau-jet candidate selection
    #--------------------------------------------------------------------------------

    patTauSelectionModules = []

    process.load("TauAnalysis.RecoTools.patLeptonSelection_cff")
    selectedPatTausAntiOverlapWithMuonsVetoName = \
        "selectedPat%ss%sAntiOverlapWithMuonsVeto" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1])
    selectedPatTausAntiOverlapWithMuonsVeto = process.selectedPatTausForMuTauAntiOverlapWithMuonsVeto.clone(
        dRmin = cms.double(0.5),
        srcNotToBeFiltered = cms.VInputTag(composeModuleName([patMuonCollectionName, "cumulative"]))
    )
    setattr(process, selectedPatTausAntiOverlapWithMuonsVetoName, selectedPatTausAntiOverlapWithMuonsVeto)
    patTauSelectionModules.append(selectedPatTausAntiOverlapWithMuonsVeto)

    process.load("RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi")
    selectedPatPFTausForTauIdEffName = \
        composeModuleName(["selectedPat%ss%s" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1]), "ForTauIdEff"])
    selectedPatPFTausForTauIdEff = cms.EDFilter("PATPFTauSelectorForTauIdEff",
        minJetPt = cms.double(15.0),
        maxJetEta = cms.double(2.3),
        trackQualityCuts = process.PFTauQualityCuts.signalQualityCuts,
        minLeadTrackPt = cms.double(5.0),
        maxDzLeadTrack = cms.double(0.2),
        maxLeadTrackPFElectronMVA = cms.double(0.6),
        #applyECALcrackVeto = cms.bool(True),
        applyECALcrackVeto = cms.bool(False),                                        
        minDeltaRtoNearestMuon = cms.double(0.5),
        muonSelection = cms.string("isGlobalMuon() | isTrackerMuon() | isStandAloneMuon()"),
        srcMuon = cms.InputTag('patMuons'),
        pfIsolation = cms.PSet(
            chargedParticleIso = cms.PSet(
                ptMin = cms.double(1.0),        
                dRvetoCone = cms.double(0.15),
                dRisoCone = cms.double(0.6)
            ),
            neutralHadronIso = cms.PSet(
                ptMin = cms.double(1000.),        
                dRvetoCone = cms.double(0.15),        
                dRisoCone = cms.double(0.)
            ),
            photonIso = cms.PSet(
                ptMin = cms.double(1.5),        
                dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
                dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
                dRvetoCone = cms.double(0.15),
                dRisoCone = cms.double(0.6)
            ),
            pileUpCorr = cms.PSet(
                method = cms.string("deltaBeta"),
                chargedToNeutralFactor = cms.double(0.50)
            )                                        
        ),
        maxPFIsoPt = cms.double(5.0),
        srcPFIsoCandidates = cms.InputTag('particleFlow'),
        srcPFNoPileUpCandidates = cms.InputTag('pfNoPileUp'),                             
        srcBeamSpot = cms.InputTag('offlineBeamSpot'),
        srcVertex = cms.InputTag('offlinePrimaryVertices'),
        filter = cms.bool(False)                                                  
    )
    # CV: comment-out for now, in order to make sure there is no bias
    #     on the tau id. efficiency measurement
    #if savePatTaus is not None:
    #    setattr(selectedPatPFTausForTauIdEff, "save", cms.string(savePatTaus))
    setattr(process, selectedPatPFTausForTauIdEffName, selectedPatPFTausForTauIdEff)
    patTauSelectionModules.append(selectedPatPFTausForTauIdEff)

    # for MC   apply L1FastJet + L2 + L3 jet-energy corrections,
    # for Data apply L1FastJet + L2 + L3 + L2/L3 residual corrections
    #
    # CV: Ztautau samples produced via MCEmbedding technique are technically "Data',
    #     L2/L3 residual jet energy corrections **must not** be applied, however,
    #     since the tau-jet response is taken from the Monte Carlo simulation
    #
    if isMC or isEmbedded:
        setattr(selectedPatPFTausForTauIdEff, "jetEnergyCorrection", cms.string('ak5PFL1FastL2L3'))
    else:
        setattr(selectedPatPFTausForTauIdEff, "jetEnergyCorrection", cms.string('ak5PFL1FastL2L3Residual'))

    patTauSelConfigurator = objSelConfigurator(
        patTauSelectionModules,
        src = patTauCollectionName,
        pyModuleName = __name__,
        doSelIndividual = False
    )
    patTauSelectionSequenceName = "selectPat%ss%sForTauIdEff" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1])
    patTauSelectionSequence = patTauSelConfigurator.configure(process = process)
    setattr(process, patTauSelectionSequenceName, patTauSelectionSequence)
    sequence += patTauSelectionSequence

    selTauCollectionName = selectedPatPFTausForTauIdEffName
    
    patTauCollections = []
    for patTauSelectionModule in patTauSelectionModules:
        patTauCollections.append(composeModuleName([patTauSelectionModule.label(), "cumulative"]))
    patTauForTauIdEffCounter = cms.EDAnalyzer("PATCandViewCountAnalyzer",
        src = cms.VInputTag(patTauCollections),
        minNumEntries = cms.int32(1),
        maxNumEntries = cms.int32(1000)                                   
    )
    patTauForTauIdEffCounterName = "pat%s%sForTauIdEffCounter" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1])
    setattr(process, patTauForTauIdEffCounterName, patTauForTauIdEffCounter)
    sequence += patTauForTauIdEffCounter

    #--------------------------------------------------------------------------------
    # produce collections of Tau-jet candidates shifted Up/Down in energy
    #
    # NOTE: Tau-jet energy shifts/smearing needs to be applied
    #       **after** PATPFTauSelectorForTauIdEff module is run,
    #       as the PATPFTauSelectorForTauIdEff module sets Tau-jet momentum to PFJet momentum,
    #       regardless of any shifts/smearing applied to Tau-jet momentum
    #
    #--------------------------------------------------------------------------------

    tauSystematics = None
    if isMC:
        # CV: shift Tau-jet energy by 3 standard-deviations,
        #     so that template morphing remains an interpolation and no extrapolation is needed
        patTausJECshiftUpModuleName = composeModuleName([patTauSelectionModule.label(), "TauJetEnUp", "cumulative"])
        patTausJECshiftUpModule = cms.EDProducer("ShiftedPATTauJetProducer",
            src = cms.InputTag(composeModuleName([selTauCollectionName, "cumulative"])),
            jetCorrPayloadName = cms.string('AK5PF'),
            jetCorrUncertaintyTag = cms.string('Uncertainty'),
            addResidualJES = cms.bool(False),                                     
            shiftBy = cms.double(+3.)
        )
        setattr(process, patTausJECshiftUpModuleName, patTausJECshiftUpModule)
        sequence += patTausJECshiftUpModule

        patTausJECshiftDownModuleName = composeModuleName([patTauSelectionModule.label(), "TauJetEnDown", "cumulative"])
        patTausJECshiftDownModule = patTausJECshiftUpModule.clone(
            shiftBy = cms.double(-3.)
        )
        setattr(process, patTausJECshiftDownModuleName, patTausJECshiftDownModule)
        sequence += patTausJECshiftDownModule

        process.load("RecoMET.METProducers.METSigParams_cfi")
        patTausJERshiftUpModuleName = composeModuleName([patTauSelectionModule.label(), "TauJetResUp", "cumulative"])
        patTausJERshiftUpModule = cms.EDProducer("SmearedPATTauJetProducer",
            src = cms.InputTag(composeModuleName([selTauCollectionName, "cumulative"])),
            inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
            lutName = cms.string('pfJetResolutionMCtoDataCorrLUT'),
            jetResolutions = process.METSignificance_params,
            srcGenJets = cms.InputTag('ak5GenJetsNoNu'),
            dRmaxGenJetMatch = cms.string('TMath::Min(0.5, 0.1 + 0.3*TMath::Exp(-0.05*genJetPt - 10.))'), 
            smearBy = cms.double(0.),                                     
            shiftBy = cms.double(+3.)                          
        )
        setattr(process, patTausJERshiftUpModuleName, patTausJERshiftUpModule)
        sequence += patTausJERshiftUpModule

        patTausJERshiftDownModuleName = composeModuleName([patTauSelectionModule.label(), "TauJetResDown", "cumulative"])
        patTausJERshiftDownModule = patTausJERshiftUpModule.clone(
            shiftBy = cms.double(-3.)
        )
        setattr(process, patTausJERshiftDownModuleName, patTausJERshiftDownModule)
        sequence += patTausJERshiftDownModule

    #--------------------------------------------------------------------------------
    # define selection of Muon + Tau-jet candidate pairs
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.CandidateTools.muTauPairProduction_cff")
    allMuTauPairsModuleName = composeModuleName(["allMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    allMuTauPairsModule = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag(composeModuleName([patMuonCollectionName, "cumulative"])),
        srcLeg2 = cms.InputTag(composeModuleName([selTauCollectionName,  "cumulative"])),
        dRmin12 = cms.double(0.5),
        srcMET = cms.InputTag(patMEtCollectionName),
        doSVreco = cms.bool(runSVfit),
        doPFMEtSign = cms.bool(runSVfit),
        doMtautauMin = cms.bool(runSVfit)
    )
    if not runSVfit:
        if hasattr(allMuTauPairsModule, "nSVfit"):
            delattr(allMuTauPairsModule, "nSVfit")
        if hasattr(allMuTauPairsModule, "pfMEtSign"):
            delattr(allMuTauPairsModule, "pfMEtSign")
    else:
        # CV: for speed reasons, run NSVfit algorithm in "fit" mode only, not in "integration" mode
        if hasattr(allMuTauPairsModule, "nSVfit") and hasattr(allMuTauPairsModule.nSVfit, "psKine_MEt_logM_int"):
            delattr(allMuTauPairsModule.nSVfit, "psKine_MEt_logM_int")
    
    if isMC:
        setattr(allMuTauPairsModule, "srcGenParticles", cms.InputTag('genParticles'))
    else:
        setattr(allMuTauPairsModule, "srcGenParticles", cms.InputTag(''))
    setattr(process, allMuTauPairsModuleName, allMuTauPairsModule)
    sequence += allMuTauPairsModule
    muTauPairProdConfigurator_systematics = None
    if isMC:
        muTauPairProdConfigurator_systematics = {
            "TauJetEnUp" : {
                "srcLeg2" : cms.InputTag(patTausJECshiftUpModuleName),
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetEnUp"]))
            },
            "TauJetEnDown" : {
                "srcLeg2" : cms.InputTag(patTausJECshiftDownModuleName),
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetEnDown"]))
            },
            "TauJetResUp" : {
                "srcLeg2" : cms.InputTag(patTausJERshiftUpModuleName),
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetResUp"]))
            },
            "TauJetResDown" : {
                "srcLeg2" : cms.InputTag(patTausJERshiftDownModuleName),
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetResDown"]))
            },
            "JetEnUp" : {
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetEnUp"]))
            },
            "JetEnDown" : {
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "JetEnDown"]))
            },
            "UnclusteredEnUp" : {
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "UnclusteredEnUp"]))
            },
            "UnclusteredEnDown" : {
                "srcMET"  : cms.InputTag(composeModuleName([patMEtCollectionName, "UnclusteredEnDown"]))
            }
        }
    muTauPairProdConfigurator = objProdConfigurator(
        allMuTauPairsModule,
        pyModuleName = __name__
    )
    if isMC:
        setattr(muTauPairProdConfigurator, "systematics", muTauPairProdConfigurator_systematics)   

    prodMuTauPairSequenceName = composeModuleName(["prodMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    prodMuTauPairSequence = muTauPairProdConfigurator.configure(process = process)
    setattr(process, prodMuTauPairSequenceName, prodMuTauPairSequence)
    sequence += prodMuTauPairSequence

    muTauPairSystematicsForTauIdEff = {
        "TauJetEnUp"        : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "TauJetEnUp"])),
        "TauJetEnDown"      : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "TauJetEnDown"])),
        "TauJetResUp"       : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "TauJetResUp"])),
        "TauJetResDown"     : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "TauJetResDown"])),
        "JetEnUp"           : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "JetEnUp"])),
        "JetEnDown"         : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "JetEnDown"])),
        "UnclusteredEnUp"   : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "UnclusteredEnUp"])),
        "UnclusteredEnDown" : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "UnclusteredEnDown"]))
    }

    #print("muTauPairSystematicsForTauIdEff:", muTauPairSystematicsForTauIdEff)

    process.load("TauAnalysis.CandidateTools.muTauPairSelection_cfi")
    selectedMuTauPairsAntiOverlapVetoModuleName = \
      composeModuleName(["selectedMu", "".join(tauIdAlgorithmName), "PairsAntiOverlapVetoForTauIdEff"])
    selectedMuTauPairsAntiOverlapVetoModule = process.selectedMuTauPairsAntiOverlapVeto.clone(
        cut = cms.string('dR12 > 0.5')
    )
    setattr(process, selectedMuTauPairsAntiOverlapVetoModuleName, selectedMuTauPairsAntiOverlapVetoModule)
    selectedMuTauPairsDzModuleName = \
      composeModuleName(["selectedMu", "".join(tauIdAlgorithmName), "PairsDzForTauIdEff"])
    selectedMuTauPairsDzModule = selectedMuTauPairsAntiOverlapVetoModule.clone(
        #cut = cms.string('abs(leg1().vertex().z() - leg2().vertex().z()) < 0.2')
        cut = cms.string('abs(leg1().vertex().z() - leg2().vertex().z()) < 1.e+3')
    )
    setattr(process, selectedMuTauPairsDzModuleName, selectedMuTauPairsDzModule)    
    muTauPairSelConfigurator = objSelConfigurator(
        [ selectedMuTauPairsAntiOverlapVetoModule,
          selectedMuTauPairsDzModule ],
        src = allMuTauPairsModuleName,
        pyModuleName = __name__,
        doSelIndividual = False
    )
    if isMC:
        setattr(muTauPairSelConfigurator, "systematics", muTauPairSystematicsForTauIdEff)
    muTauPairSelectionSequenceName = composeModuleName(["selectMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    #print "muTauPairSelectionSequenceName = %s" % muTauPairSelectionSequenceName
    muTauPairSelectionSequence = muTauPairSelConfigurator.configure(process = process)
    setattr(process, muTauPairSelectionSequenceName, muTauPairSelectionSequence)
    sequence += muTauPairSelectionSequence

    muTauPairForTauIdEffCounter = cms.EDAnalyzer("PATCandViewCountAnalyzer",
        src = cms.VInputTag(
            allMuTauPairsModuleName,
            composeModuleName([selectedMuTauPairsAntiOverlapVetoModuleName, "cumulative"]),                         
            composeModuleName([selectedMuTauPairsDzModuleName, "cumulative"])
        ),      
        minNumEntries = cms.int32(1),
        maxNumEntries = cms.int32(1000)                                   
    )
    muTauPairForTauIdEffCounterName = composeModuleName(["mu", "".join(tauIdAlgorithmName), "PairForTauIdEffCounter"])
    setattr(process, muTauPairForTauIdEffCounterName, muTauPairForTauIdEffCounter)
    sequence += muTauPairForTauIdEffCounter

    retVal = {}
    retVal["sequence"] = sequence
    retVal["patTauCollections"] = [
        composeModuleName([selTauCollectionName,                                                          "cumulative"])
    ]
    retVal["muTauPairCollections"] = [
        composeModuleName([selectedMuTauPairsDzModuleName,                                                "cumulative"])
    ]
    if isMC:
        retVal["patTauCollections"].extend([
            composeModuleName([selTauCollectionName,                                 "TauJetEnUp",        "cumulative"]),
            composeModuleName([selTauCollectionName,                                 "TauJetEnDown",      "cumulative"]),
            composeModuleName([selTauCollectionName,                                 "TauJetResUp",       "cumulative"]),
            composeModuleName([selTauCollectionName,                                 "TauJetResDown",     "cumulative"])
        ])
        retVal["muTauPairCollections"].extend([
            composeModuleName([selectedMuTauPairsDzModuleName,                       "TauJetEnUp",        "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "TauJetEnDown",      "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "TauJetResUp",       "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "TauJetResDown",     "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "JetEnUp",           "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "JetEnDown",         "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "UnclusteredEnUp",   "cumulative"]),
            composeModuleName([selectedMuTauPairsDzModuleName,                       "UnclusteredEnDown", "cumulative"])        
        ])
    return retVal
    
