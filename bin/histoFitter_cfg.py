import FWCore.ParameterSet.Config as cms
import TauAnalysis.TauIdEfficiency.plotMacroConfigVariables as variables
import TauAnalysis.TauIdEfficiency.plotMacroConfigTauIDs as tauIds

process = cms.PSet()


process.histoFitter = cms.PSet(
    channels = cms.VPSet(
        
        ),
    lumi = cms.double(),
    regions = cms.VPSet( #it may happen that some regions are needed for costraint but not for fitting
        cms.PSet(name = cms.string("ABCD"), toFit=cms.bool(False)),
        cms.PSet(name = cms.string("A"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("B"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("C"), toFit=cms.bool(False)),
        cms.PSet(name = cms.string("C1"), toFit=cms.bool(False)),
        cms.PSet(name = cms.string("C1p"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("C1f"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("C2"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("C2p"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("C2f"), toFit=cms.bool(True)),
        cms.PSet(name = cms.string("D"), toFit=cms.bool(True)),   # generic background control region (SS, tight muon isolation)
        ),
    tauIds = cms.vstring(
        tauIds.oldHPSLoose.name,
        tauIds.oldHPSMedium.name,
        tauIds.oldHPSTight.name,
        tauIds.combinedTancOnePercent.name,
        tauIds.combinedTancHalfPercent.name,
        tauIds.combinedTancQuarterPercent.name),
    sysShifts = cms.vstring("CENTRAL_VALUE",
                            "SysTauJetEnUp",
                            "SysTauJetEnDown",
                            "SysJetEnUp",
                            "SysJetEnDown",
                            "SysZllRecoilCorrectionUp",
                            "SysZllRecoilCorrectionDown"),#vstring
    fitVariables = cms.VPSet(variables.diTauMt),
    runClosureTest = cms.bool(False),
    outFile = cms.string(),
    fitSettings = cms.PSet(
        calculateStartingOSSSratio = cms.bool(True),
        calculateStartingMuonIsoRatio = cms.bool(True),
        calculateStartingTauKineBkSigRatio = cms.bool(True),
        calculateStartingTauIdPasFailRatio = cms.bool(True),
        fitTauIdEffC2 = cms.bool(True),
        fixPDiTauCharge_OS_SS = cms.vstring("Ztautau","TTplusJets"),
        fixpMuonIso_tight_loose = cms.string("Ztautau","Zmumu","TTplusJets"),
        fixpDiTauKine_Sig_Bgr = cms.string("Ztautau","Zmumu","TTplusJets"),
        fixpTauId_passed_failed = cms.vstring(), #This has some sense only for developping/stuging systematics, maybe...
        fitConstraintsC1 = cms.VPSet(
            cms.PSet( channel = cms.string("Zmumu"), nsigma = cms.double(2.)),
            cms.PSet( channel = cms.string("QCD"), nsigma = cms.double(1.)),
            cms.PSet( channel = cms.string("WplusJets"), nsigma = cms.double(1.)),
            cms.PSet( channel = cms.string("TTplusJets"), nsigma = cms.double(1.)),
            ),
        )
    )
