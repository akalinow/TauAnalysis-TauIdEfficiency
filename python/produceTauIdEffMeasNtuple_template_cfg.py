import FWCore.ParameterSet.Config as cms

from TauAnalysis.TauIdEfficiency.produceTauIdEffMeasPATTuple_template_cfg import *

#--------------------------------------------------------------------------------
#
# produce Ntuple
#
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigRunLumiSectionEventNumber_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTrigger_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigVertex_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauFixedCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauShrinkingCone_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPS_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigPFTauHPSpTaNC_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigTauIdEffMeas_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGlobalVariables_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenJets_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPhaseSpaceEventInfo_cfi")
process.load("TauAnalysis.TauIdEfficiency.ntupleConfigGenPileUpEventInfo_cfi")

process.ntupleProducer = cms.EDProducer("ObjValEDNtupleProducer",
                                        
    ntupleName = cms.string("tauIdEffNtuple"),
                                        
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc

        # run number, luminosity section, event number variables
        runLumiSectionEventNumber = process.runLumiSectionEventNumber_template,

        # variables indicating decision of HLT trigger paths
        trigger = process.trigger_template,

        # variables specifying x,y,z coordinates of primary event vertices
        vertex = process.vertex_template,

        # number of reconstructed primary event vertices
        # with sum(trackPt) exceeding different thresholds
        vertexMultiplicity = process.vertexMultiplicity_template,

        # number of reconstructed global (global || stand-alone || tracker) muons in the event
        numGlobalMuons = process.numGlobalMuons_template,
        numStandAloneMuons = process.numGoodMuons_template,                                    

        # global variables describing the underlying event/
        # amount of hadronic activity                                                                                            
        caloMet_rec = process.caloMet_template,
        pfMet_rec = process.pfMet_template,

        # variables specific to tau id. efficiency measurement
        tauIdEffMeas01 = process.tauIdEffMeas_template01.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        ),                                                   
        muonTriggerEff = process.tauIdEffMeas_muonTriggerEff.clone(
            src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative')
        )                                                                               
    )
)

patTauTemplates = {
    "pfTauFixedCone"     : process.pfTausFixedCone_recInfo,
    "pfTauShrinkingCone" : process.pfTausShrinkingCone_recInfo,
    "pfTauHPS"           : process.pfTausHPS_recInfo,
    "pfTauHPSpTaNC"      : process.pfTausHPSpTaNC_recInfo
}

muTauTemplates = {
    "pfTauFixedCone"     : [ process.tauIdEffMeas_template03pfTau ],
    "pfTauShrinkingCone" : [ process.tauIdEffMeas_template03pfTau ],
    "pfTauHPS"           : [ process.tauIdEffMeas_template03pfTau ],
    "pfTauHPSpTaNC"      : [ process.tauIdEffMeas_template03pfTau ]
}

# add branches for collections of Tau-jet candidates and Muon + Tau-Jet candidate pairs
# (including collecions shifted Up/Down in energy/momentum and with Z-recoil corrections applied)
for algorithm in patTupleConfig["algorithms"]:

    for patJetCollection in patTupleConfig[algorithm]["patJetCollections"]:
        attrJetsName = "tauIdEffMeasJets_%s" % patJetCollection
        attrJets = process.jets_template.clone(
            src = cms.InputTag(patJetCollection),
            srcNotToBeFiltered = cms.VInputTag()
        )
        setattr(process.ntupleProducer.sources, attrJetsName, attrJets)

    for patTauCollection in patTupleConfig[algorithm]["patTauCollections"]:                
        attr02Name = "tauIdEffMeas02_%s" % patTauCollection
        attr02 = patTauTemplates[algorithm].clone(
            src = cms.InputTag(patTauCollection),
            vector = cms.bool(True),
            indices = cms.vuint32([0])
        )
        setattr(process.ntupleProducer.sources, attr02Name, attr02)

    for muTauPairCollection in patTupleConfig[algorithm]["muTauPairCollections"]:
        attr03Name = "tauIdEffMeas03_%s" % muTauPairCollection
        attr03 = muTauTemplates[algorithm][0].clone(
            src = cms.InputTag(muTauPairCollection)
        )
        setattr(process.ntupleProducer.sources, attr03Name, attr03)

if isMC:
    # add in information about generator level visible taus and all generator level jets
    setattr(process.ntupleProducer.sources, "tauGenJets", process.tauGenJets_genInfo)
    setattr(process.ntupleProducer.sources, "genJets", process.genJets_genInfo)
    setattr(process.ntupleProducer.sources, "genPhaseSpaceEventInfo", process.genPhaseSpaceEventInfo_template)
    setattr(process.ntupleProducer.sources, "genPileUpEventInfo", process.genPileUpEventInfo_template)
    # add reweighting factors to be applied to Monte Carlo simulated events
    # in order to match vertex multiplicity distribution in Data                                             
    setattr(process.ntupleProducer.sources, "vertexMultReweight", process.vertexMultReweight_template)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
# Save Ntuple
#
process.ntupleOutputModule = cms.OutputModule("PoolOutputModule",
    cms.PSet(
        outputCommands = cms.untracked.vstring(
            "drop *",
            "keep *_*ntupleProducer*_*_*",
            "keep edmMergeableCounter_*_*_*"
        )                               
    ),
    process.tauIdEffSampleEventSelection,
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

#--------------------------------------------------------------------------------
#                   Modify the content of the extractors
#--------------------------------------------------------------------------------
import TauAnalysis.Configuration.pathModifiers as pathModifiers
additionalNtupleVariables = [
    [ 'PATTauVectorValExtractor', 'productionVertexZ',  "vertex().z()"                    ],
    [ 'PATTauVectorValExtractor', 'preselLoosePFIsoPt', "userFloat('preselLoosePFIsoPt')" ],
    [ 'PATTauVectorValExtractor', 'origTauPt',          "userFloat('origTauPt')"          ],
    [ 'PATTauVectorValExtractor', 'origTauEta',         "userFloat('origTauEta')"         ],
    [ 'PATTauVectorValExtractor', 'origTauPhi',         "userFloat('origTauPhi')"         ],
    [ 'PATTauVectorValExtractor', 'origTauMass',        "userFloat('origTauMass')"        ],
    [ 'PATMuTauPairValExtractor', 'tauVertexZ',         "leg2().vertex().z()"             ],
    [ 'PATMuTauPairValExtractor', 'muVertexZ',          "leg1().vertex().z()"             ]    
]
for additionalNtupleVariable in additionalNtupleVariables:
    if len(additionalNtupleVariable) != 3:
        raise ValueError("Invalid format of 'additionalNtupleVariable' !!")

    pathModifiers.ExtractorAddColumn(process.ntupleProducer.sources,
      additionalNtupleVariable[0], additionalNtupleVariable[1], cms.string(additionalNtupleVariable[2]), True)
#--------------------------------------------------------------------------------

# define order in which different paths are run
process.schedule = cms.Schedule(
    process.p,
    process.muonPFTauFixedConeSkimPath,
    process.muonPFTauShrinkingConeSkimPath,
    process.muonPFTauHPSskimPath,
    process.muonPFTauHPSpTaNCskimPath,
    process.o
)
