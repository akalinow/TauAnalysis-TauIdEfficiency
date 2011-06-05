import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.tools.configureZllRecoilCorrection import configureZllRecoilCorrection
from TauAnalysis.RecoTools.patLeptonSystematics_cff import *
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.CandidateTools.tools.composeModuleName import composeModuleName
from TauAnalysis.CandidateTools.tools.objSelConfigurator import objSelConfigurator
from TauAnalysis.CandidateTools.tools.objProdConfigurator import objProdConfigurator

def buildSequenceTauIdEffMeasSpecific(process,
                                      patMuonCollectionName = "selectedPatMuonsForTauIdEffTrkIP",
                                      tauIdAlgorithmName = None, patTauCollectionName = "patTaus",
                                      patMEtCollectionName = "patPFMETs",
                                      patTauPreselectionCriteria = None,
                                      patTauProducerPrototype = None,
                                      patTauCleanerPrototype = None,
                                      patTauTriggerMatcherProtoType = None,
                                      addGenInfo = False, applyZrecoilCorrection = False, runSVfit = False):

    #print("<buildSequenceTauIdEffMeasSpecific>:")
    #print(" patTauCollectionName = %s" % patTauCollectionName)

    # check that tauIdAlgorithmName is defined, non-null and composed of two parts
    if tauIdAlgorithmName is None or len(tauIdAlgorithmName) != 2:
        raise ValueError("Undefined of invalid 'tauIdAlgorithmName' Parameter !!")
    # check that patTauPreselectionCriteria are defined and non-null
    if patTauPreselectionCriteria is None:
        raise ValueError("Undefined 'patTauPreselectionCriteria' Parameter !!")

    patTauProducerName = "".join(patTauCollectionName)

    sequenceName = composeModuleName(["sequenceTauIdEffMeasSpecific", "".join(tauIdAlgorithmName)])
    sequence = cms.Sequence()
    setattr(process, sequenceName, sequence)

    #--------------------------------------------------------------------------------
    # produce collections of Tau-jet candidates shifted Up/Down in energy
    #--------------------------------------------------------------------------------

    patTausJECshiftUpModuleName = "%sJECshiftUp" % patTauCollectionName
    patTausJECshiftUpModule = patTausJECshiftUp.clone(
        src = cms.InputTag(patTauCollectionName)
    )
    setattr(process, patTausJECshiftUpModuleName, patTausJECshiftUpModule)
    sequence += patTausJECshiftUpModule

    patTausJECshiftDownModuleName = "%sJECshiftDown" % patTauCollectionName
    patTausJECshiftDownModule = patTausJECshiftDown.clone(
        src = cms.InputTag(patTauCollectionName)
    )
    setattr(process, patTausJECshiftDownModuleName, patTausJECshiftDownModule)
    sequence += patTausJECshiftDownModule

    tauSystematics = {   
        "sysTauJetEnUp"   : cms.InputTag(patTausJECshiftUpModuleName),
        "sysTauJetEnDown" : cms.InputTag(patTausJECshiftDownModuleName)
    }

    #--------------------------------------------------------------------------------
    # define loose Tau-jet candidate selection
    #--------------------------------------------------------------------------------

    patTauSelectionModules = []

    process.load("TauAnalysis.RecoTools.patLeptonSelection_cff")
    selectedPatTausAntiOverlapWithMuonsVetoName = \
        "selectedPat%ss%sAntiOverlapWithMuonsVeto" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1])
    selectedPatTausAntiOverlapWithMuonsVeto = process.selectedPatTausForMuTauAntiOverlapWithMuonsVeto.clone(
        dRmin = cms.double(0.7),
        srcNotToBeFiltered = cms.VInputTag(composeModuleName([patMuonCollectionName, "cumulative"]))
    )
    setattr(process, selectedPatTausAntiOverlapWithMuonsVetoName, selectedPatTausAntiOverlapWithMuonsVeto)
    patTauSelectionModules.append(selectedPatTausAntiOverlapWithMuonsVeto)

    allPatTauPreselectionCriteria = []
    allPatTauPreselectionCriteria.append(["Eta23",    "abs(eta) < 2.3"                     ])
    allPatTauPreselectionCriteria.append(["Pt20",     "pt > 20."                           ])
    allPatTauPreselectionCriteria.extend(patTauPreselectionCriteria)
    allPatTauPreselectionCriteria.append(["PFRelIso", "userFloat('pfLooseIsoPt06') < 1000."])

    lastModuleName = None
    
    for patTauPreselectionCriterion in allPatTauPreselectionCriteria:

        if len(patTauPreselectionCriterion) != 2:
            raise ValueError("Invalid format of pat::Tau selection criteria !!")
                
        patTauSelectionModuleName = "selectedPat%ss%s%s" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1], patTauPreselectionCriterion[0])
        patTauSelectionModule = cms.EDFilter("PATTauSelector",
            cut = cms.string(patTauPreselectionCriterion[1]),
            filter = cms.bool(False)                                 
        )
        setattr(process, patTauSelectionModuleName, patTauSelectionModule)
        patTauSelectionModules.append(patTauSelectionModule)

        lastModuleName = patTauSelectionModuleName

    patTauSelConfigurator = objSelConfigurator(
        patTauSelectionModules,
        src = patTauCollectionName,
        pyModuleName = __name__,
        doSelIndividual = False
    )
    setattr(patTauSelConfigurator, "systematics", tauSystematics)   
    patTauSelectionSequenceName = "selectPat%ss%sForTauIdEff" % (tauIdAlgorithmName[0], tauIdAlgorithmName[1])
    patTauSelectionSequence = patTauSelConfigurator.configure(process = process)
    setattr(process, patTauSelectionSequenceName, patTauSelectionSequence)
    sequence += patTauSelectionSequence

    selTauCollectionName = lastModuleName

    #--------------------------------------------------------------------------------
    # produce collections of Jets shifted Up/Down in energy
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    sequence += process.prodSmearedJets 

    #--------------------------------------------------------------------------------
    # produce collections of Muons, Tau-jet candidates and Jets 
    # shifted Up/Down in energy/momentum
    #--------------------------------------------------------------------------------

    patJetSelectionModules = []

    process.load("TauAnalysis.RecoTools.patJetSelection_cff")
    selectedPatJetsAntiOverlapWithLeptonsVetoModuleName = \
        composeModuleName(["selectedPatJetsForTauIdEffAntiOverlapWithLeptonsVeto", "".join(tauIdAlgorithmName)])
    selectedPatJetsAntiOverlapWithLeptonsVetoModule = process.selectedPatJetsAntiOverlapWithLeptonsVeto.clone(
        srcNotToBeFiltered = cms.VInputTag(
            "selectedPatElectronsTrkIPcumulative",
            composeModuleName([patMuonCollectionName, "cumulative"]),
            composeModuleName([selTauCollectionName,  "cumulative"])
        )
    )
    setattr(process, selectedPatJetsAntiOverlapWithLeptonsVetoModuleName, selectedPatJetsAntiOverlapWithLeptonsVetoModule)
    patJetSelectionModules.append(selectedPatJetsAntiOverlapWithLeptonsVetoModule)

    selectedPatJetsEta21ModuleName = composeModuleName(["selectedPatJetsEta21", "".join(tauIdAlgorithmName)])
    selectedPatJetsEta21Module = copy.deepcopy(process.selectedPatJetsEta21)
    setattr(process, selectedPatJetsEta21ModuleName, selectedPatJetsEta21Module)
    patJetSelectionModules.append(selectedPatJetsEta21Module)

    selectedPatJetsEt20ModuleName = composeModuleName(["selectedPatJetsEt20", "".join(tauIdAlgorithmName)])
    selectedPatJetsEt20Module = copy.deepcopy(process.selectedPatJetsEt20)
    setattr(process, selectedPatJetsEt20ModuleName, selectedPatJetsEt20Module)
    patJetSelectionModules.append(selectedPatJetsEt20Module)

    patJetSelConfigurator = objSelConfigurator(
        patJetSelectionModules,
        src = "patJets",
        pyModuleName = __name__,
        doSelIndividual = False
    )
    setattr(patJetSelConfigurator, "systematics", jetSystematics)
    patJetSelectionSequenceName = composeModuleName(["selectPatJetsForTauIdEff", "".join(tauIdAlgorithmName)])
    patJetSelectionSequence = patJetSelConfigurator.configure(process = process)
    setattr(process, patJetSelectionSequenceName, patJetSelectionSequence)
    sequence += patJetSelectionSequence

    #--------------------------------------------------------------------------------
    # produce Up/Down shifted MET collections
    #--------------------------------------------------------------------------------

    # produce collections of pat::(PF)Jets
    #  o of Pt > 10 GeV for MET Type I correction uncertainties
    #  o of Pt < 10 GeV for MET Type II ("unclustered energy") correction uncertainties
    #
    # NOTE: the "unclustered energy" is smeared by an additional resolution uncertainty of 10% * sum(Et) of jets with Pt < 10 GeV
    #      ( as defined in TauAnalysis/Configuration/python/tools/sysUncertaintyTools.py )
    #
    selectedPatJetsForMEtTypeIcorrModuleName = composeModuleName(["selectedPatJetsForMEtTypeIcorrModule", "".join(tauIdAlgorithmName)])
    selectedPatJetsForMEtTypeIcorrModule = cms.EDFilter("PATJetSelector",
        src = cms.InputTag(composeModuleName([selectedPatJetsAntiOverlapWithLeptonsVetoModuleName, "cumulative"])),
        cut = cms.string('pt > 10.'), 
        filter = cms.bool(False)
    )
    setattr(process, selectedPatJetsForMEtTypeIcorrModuleName, selectedPatJetsForMEtTypeIcorrModule)
    sequence += selectedPatJetsForMEtTypeIcorrModule

    selectedPatJetsForMEtTypeIcorrSysJetEnUpModuleName = \
        composeModuleName(["selectedPatJetsForMEtTypeIcorr", "".join(tauIdAlgorithmName), "sysJetEnUp"])
    selectedPatJetsForMEtTypeIcorrSysJetEnUpModule = process.patJetsJECshiftUp.clone(
        src = cms.InputTag(selectedPatJetsForMEtTypeIcorrModuleName)
    )
    setattr(process, selectedPatJetsForMEtTypeIcorrSysJetEnUpModuleName, selectedPatJetsForMEtTypeIcorrSysJetEnUpModule)
    sequence += selectedPatJetsForMEtTypeIcorrSysJetEnUpModule

    selectedPatJetsForMEtTypeIcorrSysJetEnDownModuleName = \
        composeModuleName(["selectedPatJetsForMEtTypeIcorr", "".join(tauIdAlgorithmName), "sysJetEnDown"])
    selectedPatJetsForMEtTypeIcorrSysJetEnDownModule = process.patJetsJECshiftDown.clone(
        src = cms.InputTag(selectedPatJetsForMEtTypeIcorrModuleName)
    )
    setattr(process, selectedPatJetsForMEtTypeIcorrSysJetEnDownModuleName, selectedPatJetsForMEtTypeIcorrSysJetEnDownModule)
    sequence += selectedPatJetsForMEtTypeIcorrSysJetEnDownModule

    selectedPatJetsForMEtTypeIIcorrModuleName = composeModuleName(["selectedPatJetsForMEtTypeIIcorr", "".join(tauIdAlgorithmName)])
    selectedPatJetsForMEtTypeIIcorrModule = selectedPatJetsForMEtTypeIcorrModule.clone(
        cut = cms.string('pt < 10.')
    )
    setattr(process, selectedPatJetsForMEtTypeIIcorrModuleName, selectedPatJetsForMEtTypeIIcorrModule)
    sequence += selectedPatJetsForMEtTypeIIcorrModule

    selectedPatJetsForMEtTypeIIcorrSysJetEnUpModuleName = \
        composeModuleName(["selectedPatJetsForMEtTypeIIcorr", "".join(tauIdAlgorithmName), "sysJetEnUp"])
    selectedPatJetsForMEtTypeIIcorrSysJetEnUpModule = process.patJetsJECshiftUp.clone(
        src = cms.InputTag(selectedPatJetsForMEtTypeIIcorrModuleName)
    )
    selectedPatJetsForMEtTypeIIcorrSysJetEnUpModule.jecUncertaintyValue = cms.double(0.10)
    setattr(process, selectedPatJetsForMEtTypeIIcorrSysJetEnUpModuleName, selectedPatJetsForMEtTypeIIcorrSysJetEnUpModule)
    sequence += selectedPatJetsForMEtTypeIIcorrSysJetEnUpModule

    selectedPatJetsForMEtTypeIIcorrSysJetEnDownModuleName = \
        composeModuleName(["selectedPatJetsForMEtTypeIIcorr", "".join(tauIdAlgorithmName), "sysJetEnDown"])
    selectedPatJetsForMEtTypeIIcorrSysJetEnDownModule = process.patJetsJECshiftDown.clone(
        src = cms.InputTag(selectedPatJetsForMEtTypeIIcorrModuleName)
    )
    selectedPatJetsForMEtTypeIIcorrSysJetEnDownModule.jecUncertaintyValue = cms.double(0.10)
    setattr(process, selectedPatJetsForMEtTypeIIcorrSysJetEnDownModuleName, selectedPatJetsForMEtTypeIIcorrSysJetEnDownModule)
    sequence += selectedPatJetsForMEtTypeIIcorrSysJetEnDownModule

    smearMEtUnclustedEnergyResolution = cms.double(0.10)

    smearedMETmoduleName = composeModuleName(["smearedMET", "".join(tauIdAlgorithmName)])
    smearedMETmodule = cms.EDProducer("SmearedMETProducer",
        src = cms.InputTag(patMEtCollectionName),
        smearedParticles = cms.PSet()
    )
    setattr(process, smearedMETmoduleName, smearedMETmodule)
    smearedMETconfigurator = objProdConfigurator(
        smearedMETmodule,
        pyModuleName = __name__,
        systematics = {
            "sysMuonPtUp" : {
                "smearedParticles.muons.srcOriginal"       : \
                    cms.InputTag(composeModuleName([patMuonCollectionName, "cumulative"])),
                "smearedParticles.muons.srcSmeared"        : \
                    cms.InputTag(composeModuleName([patMuonCollectionName, "sysMuonPtUp", "cumulative"]))
            },
            "sysMuonPtDown" : {
                "smearedParticles.muons.srcOriginal"       : \
                    cms.InputTag(composeModuleName([patMuonCollectionName, "cumulative"])),
                "smearedParticles.muons.srcSmeared"        : \
                    cms.InputTag(composeModuleName([patMuonCollectionName, "sysMuonPtDown", "cumulative"]))
            },
            "sysTauJetEnUp" : {
                "smearedParticles.taus.srcOriginal"       : \
                    cms.InputTag(composeModuleName([selTauCollectionName, "cumulative"])),
                "smearedParticles.taus.srcSmeared"        : \
                    cms.InputTag(composeModuleName([selTauCollectionName, "sysTauJetEnUp", "cumulative"]))
            },
            "sysTauJetEnDown" : {
                "smearedParticles.taus.srcOriginal"       : \
                    cms.InputTag(composeModuleName([selTauCollectionName, "cumulative"])),
                "smearedParticles.taus.srcSmeared"        : \
                    cms.InputTag(composeModuleName([selTauCollectionName, "sysTauJetEnDown", "cumulative"]))
            },
            "sysJetEnUp" : {
                "smearedParticles.jetsTypeI.srcOriginal"  : \
                    cms.InputTag(selectedPatJetsForMEtTypeIcorrModuleName),
                "smearedParticles.jetsTypeI.srcSmeared"   : \
                    cms.InputTag(selectedPatJetsForMEtTypeIcorrSysJetEnUpModuleName),
                "smearedParticles.jetsTypeII.srcOriginal" : \
                    cms.InputTag(selectedPatJetsForMEtTypeIIcorrModuleName),
                "smearedParticles.jetsTypeII.srcSmeared"  : \
                    cms.InputTag(selectedPatJetsForMEtTypeIIcorrSysJetEnUpModuleName),
                "smearedParticles.jetsTypeII.smearByResolutionUncertainty" : smearMEtUnclustedEnergyResolution
            },
            "sysJetEnDown" : {
                "smearedParticles.jetsTypeI.srcOriginal"  : \
                    cms.InputTag(selectedPatJetsForMEtTypeIcorrModuleName),
                "smearedParticles.jetsTypeI.srcSmeared"   : \
                    cms.InputTag(selectedPatJetsForMEtTypeIcorrSysJetEnDownModuleName),
                "smearedParticles.jetsTypeII.srcOriginal" : \
                    cms.InputTag(selectedPatJetsForMEtTypeIIcorrModuleName),
                "smearedParticles.jetsTypeII.srcSmeared"  : \
                    cms.InputTag(selectedPatJetsForMEtTypeIIcorrSysJetEnDownModuleName),
                "smearedParticles.jetsTypeII.smearByResolutionUncertainty" : smearMEtUnclustedEnergyResolution
            }
        }
    )
    prodSmearedMETsequenceName = composeModuleName(["prodSmearedMET", "".join(tauIdAlgorithmName)])
    prodSmearedMETsequence = smearedMETconfigurator.configure(process = process)
    setattr(process, prodSmearedMETsequenceName, prodSmearedMETsequence)
    sequence += prodSmearedMETsequence

    #--------------------------------------------------------------------------------
    # produce collections of smeared Muon + Tau-jet candidate pairs
    # **without** Z-recoil corrections applied
    #--------------------------------------------------------------------------------

     #--------------------------------------------------------------------------------
    # define selection of Muon + Tau-jet candidate pairs
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.CandidateTools.muTauPairProduction_cff")
    allMuTauPairsModuleName = composeModuleName(["allMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    allMuTauPairsModule = process.allMuTauPairs.clone(
        srcLeg1 = cms.InputTag(composeModuleName([patMuonCollectionName, "cumulative"])),
        srcLeg2 = cms.InputTag(composeModuleName([selTauCollectionName,  "cumulative"])),
        srcMET = cms.InputTag(patMEtCollectionName),
        doSVreco = cms.bool(runSVfit),
        doPFMEtSign = cms.bool(runSVfit)
    )
    if not runSVfit:
        if hasattr(allMuTauPairsModule, "nSVfit"):
            delattr(allMuTauPairsModule, "nSVfit")
        if hasattr(allMuTauPairsModule, "pfMEtSign"):
            delattr(allMuTauPairsModule, "pfMEtSign")
    
    if addGenInfo:
        setattr(allMuTauPairsModule, "srcGenParticles", cms.InputTag('genParticles'))
    else:
        setattr(allMuTauPairsModule, "srcGenParticles", cms.InputTag(''))
    setattr(process, allMuTauPairsModuleName, allMuTauPairsModule)
    sequence += allMuTauPairsModule
    muTauPairProdConfigurator = objProdConfigurator(
        allMuTauPairsModule,
        pyModuleName = __name__
    )
    setattr(muTauPairProdConfigurator, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag(composeModuleName([patMuonCollectionName, "sysMuonPtUp",     "cumulative"])),
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysMuonPtUp"]))
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag(composeModuleName([patMuonCollectionName, "sysMuonPtDown",   "cumulative"])),
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysMuonPtDown"]))
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag(composeModuleName([selTauCollectionName,  "sysTauJetEnUp",   "cumulative"])),
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysTauJetEnUp"]))
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag(composeModuleName([selTauCollectionName,  "sysTauJetEnDown", "cumulative"])),
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysTauJetEnDown"]))
        },
        "sysJetEnUp" : {
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysJetEnUp"]))
        },
        "sysJetEnDown" : {
            "srcMET"  : cms.InputTag(composeModuleName([smearedMETmoduleName,  "sysJetEnDown"]))
        }
    })
    prodMuTauPairSequenceName = composeModuleName(["prodMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    prodMuTauPairSequence = muTauPairProdConfigurator.configure(process = process)
    setattr(process, prodMuTauPairSequenceName, prodMuTauPairSequence)
    sequence += prodMuTauPairSequence

    muTauPairSystematicsForTauIdEff = {
        "sysMuonPtUp"     : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysMuonPtUp"])),
        "sysMuonPtDown"   : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysMuonPtDown"])),
        "sysTauJetEnUp"   : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysTauJetEnUp"])),
        "sysTauJetEnDown" : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysTauJetEnDown"])),
        "sysJetEnUp"      : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysJetEnUp"])),
        "sysJetEnDown"    : cms.InputTag(composeModuleName([allMuTauPairsModuleName, "sysJetEnDown"]))
    }

    #print("muTauPairSystematicsForTauIdEff:", muTauPairSystematicsForTauIdEff)

    process.load("TauAnalysis.CandidateTools.muTauPairSelection_cfi")
    selectedMuTauPairsModuleName = composeModuleName(["selectedMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    selectedMuTauPairsModule = process.selectedMuTauPairsAntiOverlapVeto.clone(
        cut = cms.string('dR12 > 0.7')
    )
    setattr(process, selectedMuTauPairsModuleName, selectedMuTauPairsModule)        
    muTauPairSelConfigurator = objSelConfigurator(
        [ selectedMuTauPairsModule ],
        src = allMuTauPairsModuleName,
        pyModuleName = __name__,
        doSelIndividual = True
    )
    setattr(muTauPairSelConfigurator, "systematics", muTauPairSystematicsForTauIdEff)
    muTauPairSelectionSequenceName = composeModuleName(["selectMu", "".join(tauIdAlgorithmName), "PairsForTauIdEff"])
    muTauPairSelectionSequence = muTauPairSelConfigurator.configure(process = process)
    setattr(process, muTauPairSelectionSequenceName, muTauPairSelectionSequence)
    sequence += muTauPairSelectionSequence

    #--------------------------------------------------------------------------------
    # produce collections of smeared Muon + Tau-jet candidate pairs
    # **with** Z-recoil corrections applied
    #--------------------------------------------------------------------------------

    if applyZrecoilCorrection:
        # produce Muon + Tau-jet candidate pairs with Z-recoil corrections applied
        configZllRecoilCorrection = \
          configureZllRecoilCorrection(process, allMuTauPairsModuleName, "ZllRecoilCorrectionMuTauPair")
        sequence += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']

        # apply Z-recoil corrections to shifted/smeared diTau objects
        muTauPairSystematicsForTauIdEffZllRecoilCorrected = {}
        for sysName, src in muTauPairSystematicsForTauIdEff.items():
            if sysName.find("ZllRecoilCorrection") != -1:
                continue;

            configZllRecoilCorrection = configureZllRecoilCorrection(process, src.getModuleLabel(), "ZllRecoilCorrectionMuTauPair")
            sequence += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']        
            muTauPairSystematicsForTauIdEffZllRecoilCorrected[sysName] = \
                configZllRecoilCorrection['diTauProducerModuleZllRecoilCorrectedName']

            # add uncertainties on Z-recoil correction
            configZllRecoilCorrection = configureZllRecoilCorrection(process, allMuTauPairsModuleName,
                                                                     "ZllRecoilCorrectionMuTauPair", +1., "SysUp")
            sequence += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
            muTauPairSystematicsForTauIdEffZllRecoilCorrected["sysZllRecoilCorrectionUp"] = \
                cms.InputTag(composeModuleName([allMuTauPairsModuleName, "ZllRecoilCorrectedSysUp"]))
            configZllRecoilCorrection = configureZllRecoilCorrection(process, allMuTauPairsModuleName,
                                                                     "ZllRecoilCorrectionMuTauPair", -1., "SysDown")
            sequence += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
            muTauPairSystematicsForTauIdEffZllRecoilCorrected["sysZllRecoilCorrectionDown"] = \
                cms.InputTag(composeModuleName([allMuTauPairsModuleName, "ZllRecoilCorrectedSysDown"]))

        #print("muTauPairSystematicsForTauIdEffZllRecoilCorrected:", muTauPairSystematicsForTauIdEffZllRecoilCorrected)

        selectedMuTauPairsZllRecoilCorrectedModuleName = \
            composeModuleName(["selectedMu", "".join(tauIdAlgorithmName), "PairsForTauIdEffZllRecoilCorrected"])
        selectedMuTauPairsZllRecoilCorrectedModule = copy.deepcopy(selectedMuTauPairsModule)
        setattr(process, selectedMuTauPairsZllRecoilCorrectedModuleName, selectedMuTauPairsZllRecoilCorrectedModule)
        muTauPairSelConfiguratorZllRecoilCorrected = objSelConfigurator(
            [ selectedMuTauPairsZllRecoilCorrectedModule ],
            src = configZllRecoilCorrection['diTauProducerModuleZllRecoilCorrectedName'],
            pyModuleName = __name__,
            doSelIndividual = True
        )    
        setattr(muTauPairSelConfiguratorZllRecoilCorrected, "systematics", muTauPairSystematicsForTauIdEffZllRecoilCorrected)
        muTauPairSelectionSequenceZllRecoilCorrectedName = \
            composeModuleName(["selectMu", "".join(tauIdAlgorithmName), "PairsForTauIdEffZllRecoilCorrected"])
        muTauPairSelectionSequenceZllRecoilCorrected = muTauPairSelConfiguratorZllRecoilCorrected.configure(process = process)
        setattr(process, muTauPairSelectionSequenceZllRecoilCorrectedName, muTauPairSelectionSequenceZllRecoilCorrected)
        sequence += muTauPairSelectionSequenceZllRecoilCorrected

    retVal = {}
    retVal["sequence"] = sequence
    retVal["patTauCollections"] = [
        composeModuleName([selTauCollectionName,                                                             "cumulative"]),
        composeModuleName([selTauCollectionName,                               "sysTauJetEnUp",              "cumulative"]),
        composeModuleName([selTauCollectionName,                               "sysTauJetEnDown",            "cumulative"])
    ]        
    retVal["muTauPairCollections"] = [
        composeModuleName([selectedMuTauPairsModuleName,                                                     "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysMuonPtUp",                "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysMuonPtDown",              "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysTauJetEnUp",              "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysTauJetEnDown",            "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysJetEnUp",                 "cumulative"]),
        composeModuleName([selectedMuTauPairsModuleName,                       "sysJetEnDown",               "cumulative"])
    ]
    if applyZrecoilCorrection:
        retVal["muTauPairCollections"].extend([
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName,                               "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysMuonPtUp",                "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysMuonPtDown",              "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysTauJetEnUp",              "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysTauJetEnDown",            "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysJetEnUp",                 "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysJetEnDown",               "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysZllRecoilCorrectionUp",   "cumulative"]),
            composeModuleName([selectedMuTauPairsZllRecoilCorrectedModuleName, "sysZllRecoilCorrectionDown", "cumulative"])
        ])
    retVal["patJetCollections"] = [
        composeModuleName([selectedPatJetsEt20ModuleName,                                                    "cumulative"]),
        composeModuleName([selectedPatJetsEt20ModuleName,                      "sysJetEnUp",                 "cumulative"]),
        composeModuleName([selectedPatJetsEt20ModuleName,                      "sysJetEnDown",               "cumulative"])
    ]    
    return retVal
    
