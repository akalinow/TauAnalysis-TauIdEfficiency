import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.patLeptonSelection_cff import *
from TauAnalysis.RecoTools.patLeptonSystematics_cff import *
from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.diTauPairSelectionAllKinds_cff import *
from TauAnalysis.RecoTools.patJetSelection_cff import *
from TauAnalysis.RecoTools.tools.configureZllRecoilCorrection import configureZllRecoilCorrection
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *

def configurePatTupleProductionTauIdEffMeasSpecific(process):

    process.producePatTupleTauIdEffMeasSpecific = cms.Sequence(
        process.selectPrimaryVertex
       + process.selectPatMuons
       + process.patTausLoosePFIsoEmbedded + process.selectPatTausForMuTau
       + process.muTauPairsForTauIdEff + process.selectMuTauPairs 
    )

    #--------------------------------------------------------------------------------
    # produce collections of Muons, Tau-jet candidates and Jets 
    # shifted Up/Down in energy/momentum
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patLeptonSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatMuons, process.prodSmearedMuons + process.selectPatMuons)
    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatTausForMuTau, process.prodSmearedTaus + process.selectPatTausForMuTau)

    setattr(patMuonSelConfigurator, "systematics", muonSystematics)
    process.selectPatMuons = patMuonSelConfigurator.configure(process = process)

    setattr(patTauSelConfiguratorForMuTau, "systematics", tauSystematics)
    process.selectPatTausForMuTau = patTauSelConfiguratorForMuTau.configure(process = process)

    process.load("TauAnalysis.RecoTools.patJetSelection_cff")
    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatJets
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedJets

    setattr(patJetSelConfigurator, "systematics", jetSystematics)
    process.selectPatJets = patJetSelConfigurator.configure(process = process)

    #--------------------------------------------------------------------------------
    # produce Muon + Tau-jet candidate pairs with Z-recoil corrections applied
    #--------------------------------------------------------------------------------

    configZllRecoilCorrection = \
      configureZllRecoilCorrection(process, "allMuTauPairs", "ZllRecoilCorrectionMuTauPair")
    process.produceMuTauPairsAll += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
    
    #--------------------------------------------------------------------------------
    # produce Up/Down shifted MET collections
    #--------------------------------------------------------------------------------

    process.smearedMET = cms.EDProducer("SmearedMETProducer",
        src = cms.InputTag('patPFMETs'),
        smearedParticles = cms.PSet()
    )
    smearedMETconfigurator = objProdConfigurator(
        process.smearedMET,
        pyModuleName = __name__,
        systematics = {
            "sysMuonPtUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatMuonsTrkIPcumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtUpCumulative')
            },
            "sysMuonPtDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatMuonsTrkIPcumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtDownCumulative')
            },
            "sysTauJetEnUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauElectronVetoSysTauJetEnUpCumulative')
            },
            "sysTauJetEnDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauElectronVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauElectronVetoSysTauJetEnDownCumulative')
            },
            "sysJetEnUp" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoSysJetEnUpCumulative')
            },
            "sysJetEnDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoSysJetEnDownCumulative')
            }
        }
    )
    process.prodSmearedMET = smearedMETconfigurator.configure(process = process)
    process.produceMuTauPairsAll += process.prodSmearedMET 

    setattr(muTauPairProdConfigurator, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtUp')
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtDown')
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauElectronVetoSysTauJetEnUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnUp')
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauElectronVetoSysTauJetEnDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnDown')
        },
        "sysJetEnUp" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnUp')
        },
        "sysJetEnDown" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnDown')
        }
    })
    process.produceMuTauPairs = muTauPairProdConfigurator.configure(process = process)
    process.produceMuTauPairsAll.replace(process.produceMuTauPairs, process.allMuTauPairs)
    process.produceMuTauPairsAll += process.produceMuTauPairs

    setattr(patMuTauPairSelConfiguratorOS, "systematics", muTauPairSystematics)
    setattr(patMuTauPairSelConfiguratorOS, "src", "allMuTauPairsZllRecoilCorrected")
    process.selectMuTauPairsOS = patMuTauPairSelConfiguratorOS.configure(process = process)
    setattr(patMuTauPairSelConfiguratorSS, "systematics", muTauPairSystematics)
    setattr(patMuTauPairSelConfiguratorSS, "src", "allMuTauPairsZllRecoilCorrected")
    process.selectMuTauPairsSS = patMuTauPairSelConfiguratorSS.configure(process = process)

    # apply Z-recoil corrections to shifted/smeared diTau objects
    for sysName, src in muTauPairSystematics.items():
        if sysName.find("ZllRecoilCorrection") != -1:
            continue;

        configZllRecoilCorrection = configureZllRecoilCorrection(process, src.getModuleLabel(), "ZllRecoilCorrectionMuTauPair")
        process.produceMuTauPairsAll += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
        muTauPairSystematics[sysName] = configZllRecoilCorrection['diTauProducerModuleZllRecoilCorrectedName']

    # add uncertainties on Z-recoil correction
    configZllRecoilCorrection = configureZllRecoilCorrection(process, "allMuTauPairs",
                                                             "ZllRecoilCorrectionMuTauPair", +1., "SysUp")
    process.produceMuTauPairsAll += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
    configZllRecoilCorrection = configureZllRecoilCorrection(process, "allMuTauPairs",
                                                             "ZllRecoilCorrectionMuTauPair", -1., "SysDown")
    process.produceMuTauPairsAll += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
