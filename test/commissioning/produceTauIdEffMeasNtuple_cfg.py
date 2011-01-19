import FWCore.ParameterSet.Config as cms

process = cms.Process("prodTauIdEffMeasNtuple")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#--------------------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        ##'file:/data1/veelken/CMSSW_3_6_x/skims/pseudoData_Ztautau.root'
        ##'file:/data1/veelken/CMSSW_3_6_x/skims/Ztautau_1_1_sXK.root'
        ##'file:/data1/veelken/CMSSW_3_8_x/skims/test/mcDYttPU156bx_GEN_SIM_RECO_1_1_1VV.root'
        'file:/data1/veelken/CMSSW_3_8_x/skims/AHtoMuTau/selEvents_AHtoMuTau_HPSloose_friis_RECO.root'
    ),
    skipEvents = cms.untracked.uint32(0)            
)

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# print gen. particle information
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(100)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

##isMC = True # use for MC
isMC = False # use for Data
HLTprocessName = "HLT" # use for non-reprocessed MC samples and Data
##HLTprocessName = "REDIGI36X" # use for Spring'10 reprocessed MC
##HLTprocessName = "REDIGI38XPU" # use for Fall'10 reprocessed MC with pile-up
##HLTprocessName = "REDIGI38X" # use for Fall'10 reprocessed MC without pile-up
pfCandidateCollection = "particleFlow" # pile-up removal disabled
##pfCandidateCollection = "pfNoPileUp" # pile-up removal enabled
applyZrecoilCorrection = False
#applyZrecoilCorrection = True
#---------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__isMC = #isMC#
#__HLTprocessName = #HLTprocessName#
#__pfCandidateCollection = #pfCandidateCollection#
#__applyZrecoilCorrection = #applyZrecoilCorrection#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define GlobalTag to be used for event reconstruction
# (only relevant for HPS tau reconstruction algorithm)
if isMC:
    process.GlobalTag.globaltag = cms.string('START38_V14::All')
else:   
    process.GlobalTag.globaltag = cms.string('GR_R_38X_V15::All')
#--------------------------------------------------------------------------------    

#--------------------------------------------------------------------------------
# define skimming criteria
# (in order to be able to produce Tau Ntuple directly from unskimmed Monte Carlo/datasets;
#  HLT single jet trigger passed && either two CaloJets or two PFJets of Pt > 10 GeV within |eta| < 2.5)
process.load('TauAnalysis.TauIdEfficiency.filterTauIdEffSample_cfi')

process.hltMu.selector.src = cms.InputTag('TriggerResults::%s' % HLTprocessName)

if isMC:
    process.dataQualityFilters.remove(process.hltPhysicsDeclared)
    process.dataQualityFilters.remove(process.dcsstatus)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce collections of objects needed as input for PAT-tuple production
# (e.g. rerun reco::Tau identification algorithms with latest tags)
#
from TauAnalysis.TauIdEfficiency.tools.configurePrePatProduction import configurePrePatProduction

configurePrePatProduction(process, pfCandidateCollection = pfCandidateCollection, addGenInfo = isMC)

process.prePatProductionSequence.remove(process.tautagging)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce PAT objects
#
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# configure PAT trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process, hltProcess = HLTprocessName, outputModule = '')
process.patTrigger.addL1Algos = cms.bool(True)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import *

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
#switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)

# comment-out to take new HPS + TaNC combined tau id. algorithm
switchToPFTauHPSpTaNC(process)

# disable preselection on of pat::Taus
# (disabled also in TauAnalysis/RecoTools/python/patPFTauConfig_cfi.py ,
#  but re-enabled after switching tau collection)
process.cleanPatTaus.preselection = cms.string('')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment to replace caloJets by pfJets
switchJetCollection(process, jetCollection = cms.InputTag("ak5PFJets"), outputModule = '')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# set Boolean swich to true in order to apply type-1 corrections
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('patMETs'), cms.InputTag('patPFMETs'))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTauIdEffMeas_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPhaseSpaceEventInfo_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPileUpEventInfo_cfi")

process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,

        # variables specifying x,y,z coordinates of primary event vertices
        vertex = process.vertex_template,

        # number of reconstructed primary event vertices
        # with sum(trackPt) exceeding different thresholds
        vertexMultiplicity = process.vertexMultiplicity_template,

        # global variables describing the underlying event/
        # amount of hadronic activity                                            
        jets = process.jets_template.clone(
            src = cms.InputTag("selectedPatJetsForTauIdEffEt20Cumulative"),
            srcNotToBeFiltered = cms.VInputTag()
        ),                                                
        caloMet_rec = process.caloMet_template,
        pfMet_rec = process.pfMet_template,

        # variables specific to tau id. efficiency measurement
        tauIdEffMeas01 = process.tauIdEffMeas_template01.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        ),
        tauIdEffMeas02 = process.tauIdEffMeas_template02.clone(
            src = cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoCumulative')
        ),
        tauIdEffMeas03woZllRecoilCorr = process.tauIdEffMeas_template03.clone(
            src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoCumulative')
        )                                        
    )
)

# add branches for collections of Muons, Tau-jet candidates and Jets shifted Up/Down in energy/momentum 
# and branches for MET with Z-recoil corrections applied
from TauAnalysis.TauIdEfficiency.tools.configurePatTupleProductionTauIdEffMeasSpecific \
       import configurePatTupleProductionTauIdEffMeasSpecific
configurePatTupleProductionTauIdEffMeasSpecific(process, applyZrecoilCorrection)

setattr(process.ntupleProducer.sources, "tauIdEffMeas01MuonPtUp", process.tauIdEffMeas_template01.clone(
    src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtUpCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrMuonPtUp", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysMuonPtUpCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas01MuonPtDown", process.tauIdEffMeas_template01.clone(
    src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPsysMuonPtDownCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrMuonPtDown", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysMuonPtDownCumulative')
))

setattr(process.ntupleProducer.sources, "tauIdEffMeas02TauJetEnUp", process.tauIdEffMeas_template02.clone(
    src = cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnUpCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrTauJetEnUp", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysTauJetEnUpCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas02TauJetEnDown", process.tauIdEffMeas_template02.clone(
    src = cms.InputTag('selectedPatTausForTauIdEffEcalCrackVetoSysTauJetEnDownCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrTauJetEnDown", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysTauJetEnDownCumulative')
))

setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrJetEnUp", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysJetEnUpCumulative')
))
setattr(process.ntupleProducer.sources, "tauIdEffMeas03woZllRecoilCorrJetEnDown", process.tauIdEffMeas_template03.clone(
    src = cms.InputTag('selectedMuTauPairsForTauIdEffAntiOverlapVetoSysJetEnDownCumulative')
))

if applyZrecoilCorrection:
    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorr", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoCumulative')
    ))

    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrMuonPtUp", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysMuonPtUpCumulative')
    ))
    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrMuonPtDown", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysMuonPtDownCumulative')
    ))

    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrTauJetEnUp", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysTauJetEnUpCumulative')
    ))
    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrTauJetEnDown", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysTauJetEnDownCumulative')
    ))

    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrJetEnUp", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysJetEnUpCumulative')
    ))
    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrJetEnDown", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysJetEnDownCumulative')
    ))

    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrZllRecoilCorrectionUp", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysZllRecoilCorrectionUpCumulative')
    ))
    setattr(process.ntupleProducer.sources, "tauIdEffMeas03wZllRecoilCorrZllRecoilCorrectionDown", process.tauIdEffMeas_template03.clone(
        src = cms.InputTag('selectedMuTauPairsForTauIdEffZllRecoilCorrectedAntiOverlapVetoSysZllRecoilCorrectionDownCumulative')
    ))

if isMC:
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
    setattr(process.ntupleProducer.sources, "genPileUpEventInfo", process.genPileUpEventInfo_template)

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save ntuple
#
process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            "drop *",
            "keep *_*ntupleProducer*_*_*"
        )                               
    ),
    process.tauIdEffSampleEventSelection,
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("tauIdEffMeasEDNtuple.root")      
)
#--------------------------------------------------------------------------------

process.p = cms.Path(
    process.prePatProductionSequence
   + process.patDefaultSequence
   + process.producePatTupleTauIdEffMeasSpecific
   #+ process.printEventContent
   #+ process.printGenParticleList
   + process.ntupleProducer
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.ntupleOutputModule)

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    process.muonPFTauSkimPath,
    process.o
)

# print-out all python configuration parameter information
#print process.dumpPython()

