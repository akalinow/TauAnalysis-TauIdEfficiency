import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.patLeptonSelection_cff import *
from TauAnalysis.RecoTools.patLeptonSystematics_cff import *
from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.diTauPairSelectionAllKinds_cff import *
from TauAnalysis.RecoTools.patJetSelection_cff import *
from TauAnalysis.RecoTools.tools.configureZllRecoilCorrection import configureZllRecoilCorrection
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *
from TauAnalysis.TauIdEfficiency.filterTauIdEffSample_cfi import \
       patMuonSelConfiguratorForTauIdEff, patTauSelConfiguratorForTauIdEff, patMuTauPairSelConfiguratorForTauIdEff

def configurePatTupleProductionTauIdEffMeasSpecific(process):

    #--------------------------------------------------------------------------------
    # produce collections for "central values" of Muons and Tau-jets
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patMuonMomentumCorrection_cfi")
    
    process.patMuonsLoosePFIsoEmbedded04.src = cms.InputTag('patMuonsMuScleFitCorrectedMomentum')
    
    process.producePatTupleTauIdEffMeasSpecific = cms.Sequence(
        process.selectPrimaryVertex
       + process.patMuonsMuScleFitCorrectedMomentum + process.patMuonsLoosePFIsoEmbedded + process.selectPatMuonsForTauIdEff
       + process.patTausLoosePFIsoEmbedded + process.selectPatTausForTauIdEff
       + process.muTauPairsForTauIdEff
    )

    #--------------------------------------------------------------------------------
    # produce collections of Muons, Tau-jet candidates and Jets 
    # shifted Up/Down in energy/momentum
    #--------------------------------------------------------------------------------

    process.load("TauAnalysis.RecoTools.patLeptonSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatMuonsForTauIdEff,
                                                        process.prodSmearedMuons + process.selectPatMuonsForTauIdEff)
    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatTausForTauIdEff,
                                                        process.prodSmearedTaus + process.selectPatTausForTauIdEff)

    setattr(patMuonSelConfiguratorForTauIdEff, "systematics", muonSystematics)
    process.selectPatMuonsForTauIdEff = patMuonSelConfiguratorForTauIdEff.configure(process = process)

    setattr(patTauSelConfiguratorForTauIdEff, "systematics", tauSystematics)
    process.selectPatTausForTauIdEff = patTauSelConfiguratorForTauIdEff.configure(process = process)

    # CV: pat::Electron and pat::Tau collections needed for pat::Jet overlap removal
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatElectrons
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatTaus

    process.load("TauAnalysis.RecoTools.patJetSelection_cff")
    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedJets 
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatJets

    setattr(patJetSelConfigurator, "systematics", jetSystematics)    
    process.selectPatJets = patJetSelConfigurator.configure(process = process)

    #--------------------------------------------------------------------------------
    # produce Muon + Tau-jet candidate pairs with Z-recoil corrections applied
    #--------------------------------------------------------------------------------

    configZllRecoilCorrection = \
      configureZllRecoilCorrection(process, "muTauPairsForTauIdEff", "ZllRecoilCorrectionMuTauPair")
    process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
    
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
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoSysTauJetEnUpCumulative')
            },
            "sysTauJetEnDown" : {
                "smearedParticles.srcOriginal" : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoCumulative'),
                "smearedParticles.srcSmeared"  : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoSysTauJetEnDownCumulative')
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
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedMET 

    muTauPairProdConfiguratorForTauIdEff = objProdConfigurator(
        process.muTauPairsForTauIdEff,
        pyModuleName = __name__
    )

    setattr(muTauPairProdConfiguratorForTauIdEff, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtUp')
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsTrkIPsysMuonPtDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtDown')
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoSysTauJetEnUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnUp')
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForMuTauEcalCrackVetoSysTauJetEnDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnDown')
        },
        "sysJetEnUp" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnUp')
        },
        "sysJetEnDown" : {
            "srcMET"  : cms.InputTag('smearedMETsysJetEnDown')
        }
    })
    process.produceMuTauPairsForTauIdEff = muTauPairProdConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.produceMuTauPairsForTauIdEff

    muTauPairSystematicsForTauIdEff = {
        "sysMuonPtUp"                : cms.InputTag('muTauPairsForTauIdEffSysMuonPtUp'),
        "sysMuonPtDown"              : cms.InputTag('muTauPairsForTauIdEffSysMuonPtDown'),
        "sysTauJetEnUp"              : cms.InputTag('muTauPairsForTauIdEffSysTauJetEnUp'),
        "sysTauJetEnDown"            : cms.InputTag('muTauPairsForTauIdEffSysTauJetEnDown'),
        "sysJetEnUp"                 : cms.InputTag('muTauPairsForTauIdEffSysJetEnUp'),
        "sysJetEnDown"               : cms.InputTag('muTauPairsForTauIdEffSysJetEnDown'),
        "sysZllRecoilCorrectionUp"   : cms.InputTag('muTauPairsForTauIdEffZllRecoilCorrectedSysUp'),
        "sysZllRecoilCorrectionDown" : cms.InputTag('muTauPairsForTauIdEffZllRecoilCorrectedSysDown')
    }

    # apply Z-recoil corrections to shifted/smeared diTau objects
    for sysName, src in muTauPairSystematicsForTauIdEff.items():
        if sysName.find("ZllRecoilCorrection") != -1:
            continue;

        configZllRecoilCorrection = configureZllRecoilCorrection(process, src.getModuleLabel(), "ZllRecoilCorrectionMuTauPair")
        process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
        muTauPairSystematicsForTauIdEff[sysName] = configZllRecoilCorrection['diTauProducerModuleZllRecoilCorrectedName']

    # add uncertainties on Z-recoil correction
    configZllRecoilCorrection = configureZllRecoilCorrection(process, "muTauPairsForTauIdEff",
                                                             "ZllRecoilCorrectionMuTauPair", +1., "SysUp")
    process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
    configZllRecoilCorrection = configureZllRecoilCorrection(process, "muTauPairsForTauIdEff",
                                                             "ZllRecoilCorrectionMuTauPair", -1., "SysDown")
    process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']

    setattr(patMuTauPairSelConfiguratorForTauIdEff, "systematics", muTauPairSystematicsForTauIdEff)
    setattr(patMuTauPairSelConfiguratorForTauIdEff, "src", "muTauPairsForTauIdEffZllRecoilCorrected")
    process.selectMuTauPairsForTauIdEff = patMuTauPairSelConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.selectMuTauPairsForTauIdEff

    
