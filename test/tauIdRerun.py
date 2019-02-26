import FWCore.ParameterSet.Config as cms

def addMVA_WPs_run2_2017(process):

    print "Setup new TauIDs to be evaluated at MiniAOD"
    process.load("PhysicsTools.NanoAOD.taus_updatedMVAIds_cff")
    process.slimmedTausNewMVAIDs = cms.EDProducer("PATTauIDEmbedder",
                                                  src = cms.InputTag('slimmedTaus'),
                                                  tauIDSources = cms.PSet(
                                                      #oldDM
            byIsolationMVArun2v1DBoldDMwLTraw2017v2 = cms.InputTag('patTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw'),
                                                      byVVLooseIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVVLooseIsolationMVArun2v1DBoldDMwLT'),
                                                      byVLooseIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT'),
                                                      byLooseIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT'),
                                                      byMediumIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT'),
                                                      byTightIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT'),
                                                      byVTightIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT'),
                                                      byVVTightIsolationMVArun2v1DBoldDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT'),
                                                      #newDM
            byIsolationMVArun2v1DBnewDMwLTraw2017v2 = cms.InputTag('patTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw'),
                                                      byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVVLooseIsolationMVArun2v1DBnewDMwLT'),
                                                      byVLooseIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT'),
                                                      byLooseIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT'),
                                                      byMediumIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT'),
                                                      byTightIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT'),
                                                      byVTightIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT'),
                                                      byVVTightIsolationMVArun2v1DBnewDMwLT2017v2 = cms.InputTag('patTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT'),
                                                  )
    )
    process.newTauMVAIDsSeq = cms.Sequence(
        process.patTauDiscriminationByIsolationMVArun2v1DBoldDMwLTSeq+
        process.patTauDiscriminationByIsolationMVArun2v1DBnewDMwLTSeq+
        process.slimmedTausNewMVAIDs
    )

    '''
    #DNN-based tau-IDs
    # deepTau
    from RecoTauTag.RecoTau.DeepTauId_cfi import deepTauIdraw
    process.deepTauIdraw = deepTauIdraw.clone(
        electrons = cms.InputTag('slimmedElectrons'),
        muons = cms.InputTag('slimmedMuons'),
        taus = cms.InputTag('slimmedTaus'),
        graph_file = cms.string('RecoTauTag/TrainingFiles/data/DeepTauId/deepTau_2017v1_20L1024N.pb')
    )
    process.newTauMVAIDsSeq.replace(process.slimmedTausNewMVAIDs,
                                    process.deepTauIdraw+process.slimmedTausNewMVAIDs)
    deepTau2017v1Sources_ = cms.PSet(        
        deepTau2017v1tauVSe = cms.InputTag('deepTauIdraw', 'tauVSe'),
        deepTau2017v1tauVSmu = cms.InputTag('deepTauIdraw', 'tauVSmu'),
        deepTau2017v1tauVSjet = cms.InputTag('deepTauIdraw', 'tauVSjet'),
        deepTau2017v1tauVSall = cms.InputTag('deepTauIdraw', 'tauVSall'),
    )
    process.slimmedTausNewMVAIDs.tauIDSources = cms.PSet(
        process.slimmedTausNewMVAIDs.tauIDSources,
        deepTau2017v1Sources_
    )
    ## DPFTau2016v0/v1
    from RecoTauTag.RecoTau.DPFIsolation_cfi import DPFIsolation
    process.DPFIsolationv0 = DPFIsolation.clone(
        electrons = cms.InputTag('slimmedElectrons'),
        muons = cms.InputTag('slimmedMuons'),
        taus = cms.InputTag('slimmedTaus'),
        graph_file = cms.string('RecoTauTag/TrainingFiles/data/DPFTauId/DPFIsolation_2017v0.pb')
    )
    process.DPFIsolationv1 = DPFIsolation.clone(
        electrons = cms.InputTag('slimmedElectrons'),
        muons = cms.InputTag('slimmedMuons'),
        taus = cms.InputTag('slimmedTaus'),
        graph_file = cms.string('RecoTauTag/TrainingFiles/data/DPFTauId/DPFIsolation_2017v1.pb')
    )
    process.dpfIsoSeq = cms.Sequence(process.DPFIsolationv0+process.DPFIsolationv1)
    process.newTauMVAIDsSeq.replace(process.slimmedTausNewMVAIDs,
                                    process.dpfIsoSeq+process.slimmedTausNewMVAIDs)
    dpfTau2016Sources_ = cms.PSet(        
        DPFTau_2016_v0tauVSall = cms.InputTag('DPFIsolationv0', 'tauVSall'),

        DPFTau_2016_v1tauVSall = cms.InputTag('DPFIsolationv1', 'tauVSall'),
    )
    process.slimmedTausNewMVAIDs.tauIDSources = cms.PSet(
        process.slimmedTausNewMVAIDs.tauIDSources,
        dpfTau2016Sources_
    )
    '''
   




