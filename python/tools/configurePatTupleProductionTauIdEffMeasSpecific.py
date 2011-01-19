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

def configurePatTupleProductionTauIdEffMeasSpecific(process, applyZrecoilCorrection = False):

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

    # set InputTags of JEC shifted pat::Tau collections
    # to pat::Taus with pfLooseIsolation embedded
    process.patTausJECshiftUp.src = cms.InputTag("patTausLoosePFIsoEmbedded06")
    process.patTausJECshiftDown.src = cms.InputTag("patTausLoosePFIsoEmbedded06")

    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatMuonsForTauIdEff,
                                                        process.prodSmearedMuons + process.selectPatMuonsForTauIdEff)
    process.producePatTupleTauIdEffMeasSpecific.replace(process.selectPatTausForTauIdEff,
                                                        process.prodSmearedTaus + process.selectPatTausForTauIdEff)

    setattr(patMuonSelConfiguratorForTauIdEff, "systematics", muonSystematics)
    process.selectPatMuonsForTauIdEff = patMuonSelConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatMuonsForTauIdEff
    
    setattr(patTauSelConfiguratorForTauIdEff, "systematics", tauSystematics)
    process.selectPatTausForTauIdEff = patTauSelConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatTausForTauIdEff

    # CV: pat::Electron collection needed for pat::Jet overlap removal
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatElectrons

    process.load("TauAnalysis.RecoTools.patJetSelection_cff")
    process.selectedPatJetsForTauIdEffAntiOverlapWithLeptonsVeto = copy.deepcopy(process.selectedPatJetsAntiOverlapWithLeptonsVeto)
    process.selectedPatJetsForTauIdEffAntiOverlapWithLeptonsVeto.srcNotToBeFiltered = cms.VInputTag(
        "selectedPatElectronsTrkIPcumulative",
        "selectedPatMuonsForTauIdEffTrkIPcumulative",
        "selectedPatTausForTauIdEffEcalCrackVetoCumulative"
    )
    process.selectedPatJetsForTauIdEffEta21 = copy.deepcopy(process.selectedPatJetsEta21)
    process.selectedPatJetsForTauIdEffEt20 = copy.deepcopy(process.selectedPatJetsEt20)
    patJetSelConfiguratorForTauIdEff = objSelConfigurator(
        [ process.selectedPatJetsForTauIdEffAntiOverlapWithLeptonsVeto,
          process.selectedPatJetsForTauIdEffEta21,
          process.selectedPatJetsForTauIdEffEt20 ],
        src = "patJets",
        pyModuleName = __name__,
        doSelIndividual = False
    )
    setattr(patJetSelConfiguratorForTauIdEff, "systematics", jetSystematics)   
    process.selectPatJetsForTauIdEff = patJetSelConfiguratorForTauIdEff.configure(process = process)

    process.load("TauAnalysis.RecoTools.patJetSystematics_cff")
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedJets 
    process.producePatTupleTauIdEffMeasSpecific += process.selectPatJetsForTauIdEff

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
    process.selectedPatJetsForMEtTypeIcorr = cms.EDFilter("PATJetSelector",
        src = cms.InputTag('selectedPatJetsForTauIdEffAntiOverlapWithLeptonsVetoCumulative'),                                    
        cut = cms.string('pt > 10.'), 
        filter = cms.bool(False)
    )

    process.selectedPatJetsForMEtTypeIcorrSysJetEnUp = process.patJetsJECshiftUp.clone(
        src = cms.InputTag('selectedPatJetsForMEtTypeIcorr')
    )
    process.selectedPatJetsForMEtTypeIcorrSysJetEnDown = process.patJetsJECshiftDown.clone(
        src = cms.InputTag('selectedPatJetsForMEtTypeIcorr')
    )

    process.selectedPatJetsForMEtTypeIIcorr = process.selectedPatJetsForMEtTypeIcorr.clone(
        cut = cms.string('pt < 10.')
    )

    process.selectedPatJetsForMEtTypeIIcorrSysJetEnUp = process.patJetsJECshiftUp.clone(
        src = cms.InputTag('selectedPatJetsForMEtTypeIIcorr')
    )
    process.selectedPatJetsForMEtTypeIIcorrSysJetEnUp.jecUncertaintyValue = cms.double(0.10)
    process.selectedPatJetsForMEtTypeIIcorrSysJetEnDown = process.patJetsJECshiftDown.clone(
        src = cms.InputTag('selectedPatJetsForMEtTypeIIcorr')
    )
    process.selectedPatJetsForMEtTypeIIcorrSysJetEnDown.jecUncertaintyValue = cms.double(0.10)

    process.producePatJetsForMEtCorrUncertainty = cms.Sequence(
        process.selectedPatJetsForMEtTypeIcorr
       * process.selectedPatJetsForMEtTypeIcorrSysJetEnUp * process.selectedPatJetsForMEtTypeIcorrSysJetEnDown
       * process.selectedPatJetsForMEtTypeIIcorr
       * process.selectedPatJetsForMEtTypeIIcorrSysJetEnUp * process.selectedPatJetsForMEtTypeIIcorrSysJetEnDown
    )
    process.producePatTupleTauIdEffMeasSpecific += process.producePatJetsForMEtCorrUncertainty

    smearMEtUnclustedEnergyResolution = cms.double(0.10)

    process.smearedMET = cms.EDProducer("SmearedMETProducer",
        src = cms.InputTag('patPFMETs'),
        smearedParticles = cms.PSet()
    )
    smearedMETconfigurator = objProdConfigurator(
        process.smearedMET,
        pyModuleName = __name__,
        systematics = {
            "sysMuonPtUp" : {
                "smearedParticles.muons.srcOriginal" : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
                "smearedParticles.muons.srcSmeared"  : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtUpCumulative')
            },
            "sysMuonPtDown" : {
                "smearedParticles.muons.srcOriginal" : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
                "smearedParticles.muons.srcSmeared"  : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtDownCumulative')
            },
            "sysTauJetEnUp" : {
                "smearedParticles.taus.srcOriginal" : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoCumulative'),
                "smearedParticles.taus.srcSmeared"  : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnUpCumulative')
            },
            "sysTauJetEnDown" : {
                "smearedParticles.taus.srcOriginal" : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoCumulative'),
                "smearedParticles.taus.srcSmeared"  : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnDownCumulative')
            },
            "sysJetEnUp" : {
                "smearedParticles.jetsTypeI.srcOriginal" : cms.InputTag('selectedPatJetsForMEtTypeIcorr'),
                "smearedParticles.jetsTypeI.srcSmeared"  : cms.InputTag('selectedPatJetsForMEtTypeIcorrSysJetEnUp'),
                "smearedParticles.jetsTypeII.srcOriginal" : cms.InputTag('selectedPatJetsForMEtTypeIIcorr'),
                "smearedParticles.jetsTypeII.srcSmeared"  : cms.InputTag('selectedPatJetsForMEtTypeIIcorrSysJetEnUp'),
                "smearedParticles.jetsTypeII.smearByResolutionUncertainty" : smearMEtUnclustedEnergyResolution
            },
            "sysJetEnDown" : {
                "smearedParticles.jetsTypeI.srcOriginal" : cms.InputTag('selectedPatJetsForMEtTypeIcorr'),
                "smearedParticles.jetsTypeI.srcSmeared"  : cms.InputTag('selectedPatJetsForMEtTypeIcorrSysJetEnDown'),
                "smearedParticles.jetsTypeII.srcOriginal" : cms.InputTag('selectedPatJetsForMEtTypeIIcorr'),
                "smearedParticles.jetsTypeII.srcSmeared"  : cms.InputTag('selectedPatJetsForMEtTypeIIcorrSysJetEnDown'),
                "smearedParticles.jetsTypeII.smearByResolutionUncertainty" : smearMEtUnclustedEnergyResolution
            }
        }
    )
    process.prodSmearedMET = smearedMETconfigurator.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.prodSmearedMET 

    #--------------------------------------------------------------------------------
    # produce collections of smeared Muon + Tau-jet candidate pairs
    # **without** Z-recoil corrections applied
    #--------------------------------------------------------------------------------

    muTauPairProdConfiguratorForTauIdEff = objProdConfigurator(
        process.muTauPairsForTauIdEff,
        pyModuleName = __name__
    )

    setattr(muTauPairProdConfiguratorForTauIdEff, "systematics", {
        "sysMuonPtUp" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtUp')
        },
        "sysMuonPtDown" : {
            "srcLeg1" : cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtDownCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysMuonPtDown')
        },
        "sysTauJetEnUp" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnUpCumulative'),
            "srcMET"  : cms.InputTag('smearedMETsysTauJetEnUp')
        },
        "sysTauJetEnDown" : {
            "srcLeg2" : cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnDownCumulative'),
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
        "sysJetEnDown"               : cms.InputTag('muTauPairsForTauIdEffSysJetEnDown')
    }

    print("muTauPairSystematicsForTauIdEff:", muTauPairSystematicsForTauIdEff)
    
    setattr(patMuTauPairSelConfiguratorForTauIdEff, "systematics", muTauPairSystematicsForTauIdEff)    
    process.selectMuTauPairsForTauIdEff = \
      patMuTauPairSelConfiguratorForTauIdEff.configure(process = process)
    process.producePatTupleTauIdEffMeasSpecific += process.selectMuTauPairsForTauIdEff

    #--------------------------------------------------------------------------------
    # produce collections of smeared Muon + Tau-jet candidate pairs
    # **with** Z-recoil corrections applied
    #--------------------------------------------------------------------------------

    if applyZrecoilCorrection:
        # produce Muon + Tau-jet candidate pairs with Z-recoil corrections applied
        configZllRecoilCorrection = \
          configureZllRecoilCorrection(process, "muTauPairsForTauIdEff", "ZllRecoilCorrectionMuTauPair")
        process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']

        # apply Z-recoil corrections to shifted/smeared diTau objects
        muTauPairSystematicsForTauIdEffZllRecoilCorrected = {}
        for sysName, src in muTauPairSystematicsForTauIdEff.items():
            if sysName.find("ZllRecoilCorrection") != -1:
                continue;

            configZllRecoilCorrection = configureZllRecoilCorrection(process, src.getModuleLabel(), "ZllRecoilCorrectionMuTauPair")
            process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
            muTauPairSystematicsForTauIdEffZllRecoilCorrected[sysName] = \
              configZllRecoilCorrection['diTauProducerModuleZllRecoilCorrectedName']

            # add uncertainties on Z-recoil correction
            configZllRecoilCorrection = configureZllRecoilCorrection(process, "muTauPairsForTauIdEff",
                                                                 "ZllRecoilCorrectionMuTauPair", +1., "SysUp")
            process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
            muTauPairSystematicsForTauIdEffZllRecoilCorrected["sysZllRecoilCorrectionUp"] = \
              cms.InputTag('muTauPairsForTauIdEffZllRecoilCorrectedSysUp')
            configZllRecoilCorrection = configureZllRecoilCorrection(process, "muTauPairsForTauIdEff",
                                                                 "ZllRecoilCorrectionMuTauPair", -1., "SysDown")
            process.producePatTupleTauIdEffMeasSpecific += configZllRecoilCorrection['patPFMETsZllRecoilCorrectionSequence']
            muTauPairSystematicsForTauIdEffZllRecoilCorrected["sysZllRecoilCorrectionDown"] = \
              cms.InputTag('muTauPairsForTauIdEffZllRecoilCorrectedSysDown')

        print("muTauPairSystematicsForTauIdEffZllRecoilCorrected:", muTauPairSystematicsForTauIdEffZllRecoilCorrected)

        process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVeto = \
          copy.deepcopy(process.selectedMuTauPairsForTauIdEffAntiOverlapVeto)
        process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedMt1MET = \
          copy.deepcopy(process.selectedMuTauPairsForTauIdEffMt1MET)
        process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedPzetaDiff = \
          copy.deepcopy(process.selectedMuTauPairsForTauIdEffPzetaDiff)
        process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedZeroCharge = \
          copy.deepcopy(process.selectedMuTauPairsForTauIdEffZeroCharge)

        patMuTauPairSelConfiguratorForTauIdEffZllRecoilCorrected = objSelConfigurator(
            [ process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVeto,
              process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedMt1MET,
              process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedPzetaDiff,
              process.selectedMuTauPairsForTauIdEffZllRecoilCorrectedZeroCharge ],
            src = "muTauPairsForTauIdEffZllRecoilCorrected",
            pyModuleName = __name__,
            doSelIndividual = True
        )
    
        setattr(patMuTauPairSelConfiguratorForTauIdEffZllRecoilCorrected, "systematics", muTauPairSystematicsForTauIdEffZllRecoilCorrected)
        process.selectMuTauPairsForTauIdEffZllRecoilCorrected = \
          patMuTauPairSelConfiguratorForTauIdEffZllRecoilCorrected.configure(process = process)
        process.producePatTupleTauIdEffMeasSpecific += process.selectMuTauPairsForTauIdEffZllRecoilCorrected 
    
