import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTupleSpecific_template_cfg import *

##for processAttrName in dir(process):
##    processAttr = getattr(process, processAttrName)
##    if isinstance(processAttr, cms.EDFilter):
##        if processAttr.type_() == "PATTauSelector":
##            print "--> Disabling cut %s" % processAttrName
##            setattr(processAttr, "cut", cms.string('pt > 10. & abs(eta) < 2.5'))

for processAttrName in dir(process):
    processAttr = getattr(process, processAttrName)
    if isinstance(processAttr, cms.EDFilter):
        if processAttr.type_() == "PATPFTauSelectorForTauIdEff":
            print "--> Disabling cut %s" % processAttrName
            setattr(processAttr, "cut", cms.string('pt > 0.'))
            setattr(processAttr, "produceAll", cms.bool(True))

if hasattr(process.patTupleOutputModule, "SelectEvents"):
    delattr(process.patTupleOutputModule, "SelectEvents")

process.patPFTausCleanedHPSpTaNC.finalCut = cms.string('') #Disable cut 'pfJetRef().pt() > 5 & abs(pfJetRef().eta()) < 2.5'
process.patPFTausVertexMatchedHPSpTaNC.dzMax = cms.double(1000) #Disable vertex cut
process.selectedPatPFTausHPSpTaNCAntiOverlapWithMuonsVetoCumulative.dRmin = cms.double(0.) #Disable dRmin wrt muon
process.selectedPatPFTausHPSpTaNCForTauIdEffCumulative.produceAll = cms.bool(True) #Turns off all the cuts but produces the userFloats

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
#    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_10_1_h65.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_11_0_erV.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_12_0_r8u.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_13_0_EDP.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_14_0_kJO.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_15_0_n6K.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_16_0_MlT.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_17_0_NqG.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_18_0_YBQ.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_19_0_s5y.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_1_1_5dZ.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_20_0_zj0.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_21_0_JLp.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_22_0_HwW.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_23_0_W47.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_24_0_cl7.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_25_1_UNY.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_2_1_DLd.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_3_1_Xzs.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_4_1_l2x.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_5_1_cFw.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_6_1_TSd.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_7_1_QSe.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_8_1_PVt.root',
    'rfio:/castor/cern.ch/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/skimGenZtoMuTauWithinAcc_Ztautau_2011Jun30v2_9_1_rJ0.root'                                ),
  skipEvents = cms.untracked.uint32(0)            
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents.input = cms.untracked.int32(-1) 

process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
    ntupleName = cms.string("tauIdEffNtuple"),
    sources = cms.PSet(
              vertexMultReweight = process.vertexMultReweight_template
              )
)

process.totalEventsProcessed = cms.EDProducer("EventCountProducer")
process.addenda = cms.Path(process.totalEventsProcessed+process.ntupleProducer)

import TauAnalysis.Configuration.pathModifiers as pathModifiers
#Removes Sys*Up collections, case insensitive
pathModifiers.PathRemover(process.p,['sys','up'],True)
pathModifiers.PathRemover(process.p,['sys','down'],True)

process.patTupleOutputModule = cms.OutputModule("PoolOutputModule",
                                              cms.PSet(
                                                  outputCommands = cms.untracked.vstring(
                                                      'drop *',
                                                      'keep edmMergeableCounter_*_*_*',
                                                      'keep *TriggerEvent_*_*_*',
                                                      'keep patMuons_selectedPatMuonsForTauIdEffTrkIPcumulative_*_*',
                                                      'keep *_patMuonsStandAlone_*_*',
                                                      'keep *_patMuons_*_*',
                                                      'keep *_patMuonsLoosePFIsoEmbedded06_*_*',
                                                      'keep patTaus_selectedPatPFTausShrinkingConePFRelIsoCumulative_*_*',
                                                      'keep *_selectedPatPFTausHPSpTaNCForTauIdEffCumulative_*_*',
                                                      'keep *Vertex*_*_*_*',
                                                      'keep *_ak5PFJets_*_*',
                                                      'keep *_ak5CaloJets_*_*',
                                                      'keep *PFCandidate*_*_*_*',
                                                      'keep *GenJet*_*_*_*',
#                                                      'keep *recoGenMET_*_*_*',
                                                      'keep *_*Met*_*_*',
                                                      'keep *_*MET*_*_*',
                                                      'keep *_genParticles_*_*',
                                                      'keep *_generator_*_*', 
#                                                      'keep *_*GenJets_*_*',
                                                      'keep *PFCandidate*_*_*_*',
                                                      'keep *_addPileupInfo_*_*',
                                                      'keep *_tauGenJets_*_*'
                                                      )               
                                                  ),
#                                              process.tauIdEffSampleEventSelection,
                                              fileName = cms.untracked.string("/nfs/data4/verzetti/tagprobe/NTuples/Ztautau_muhadr_noselection_V2/tauIdEffMeasPatTuple_noSelection.root"),
                                              maxSize = cms.untracked.int32(10000) #Limit the file size to 1GB then opens another one (or at least should)
                                              )

process.schedule = cms.Schedule(process.addenda,process.p,process.muonPFTauFixedConeSkimPath,process.muonPFTauShrinkingConeSkimPath,process.muonPFTauHPSskimPath,process.muonPFTauHPSpTaNCskimPath,process.o)

processDumpFile = open('produceTauIdEffMeasPATTuple_noTauPresel.dump' , 'w')
print >> processDumpFile, process.dumpPython()



