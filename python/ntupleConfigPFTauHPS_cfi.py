import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to PFTaus
# reconstructed by hadron + strips (HPS) algorithm
#--------------------------------------------------------------------------------

pfTausHPS_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeHPS"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for PFTau
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),

        # kinematic variables for PFJet associated to PFTau
        jetPt = cms.string("pfTauTagInfoRef().pfjetRef().pt()"),
        jetEta = cms.string("pfTauTagInfoRef().pfjetRef().eta()"),
        jetPhi = cms.string("pfTauTagInfoRef().pfjetRef().phi()"),

        # flags for tag/probe and Pt_index for distinguishing between
        # highest Pt, second highest Pt and third highest Pt jet
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
        ptIndex = cms.string("userFloat('pt_index')"),
        
        # tau id. discriminators
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byIsolationLoose = cms.string("tauID('byLooseIsolation')"),
        byIsolationMedium = cms.string("tauID('byMediumIsolation')"),
        byIsolationTight = cms.string("tauID('byTightIsolation')"),
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')"),                             

        mass = cms.string("mass"),
        prongs = cms.string("signalPFChargedHadrCands().size()"),
        gammas = cms.string("signalPFGammaCands().size()"),
        nIsoCharged = cms.string("isolationPFChargedHadrCands().size()"),
        nIsoGamma = cms.string("isolationPFGammaCands().size()"),
        nIsoNeutral = cms.string("isolationPFNeutrHadrCands().size()"),
        isoChargedEtSum = cms.string("isolationPFChargedHadrCandsPtSum()"),
        isoGammaEtSum = cms.string("isolationPFGammaCandsEtSum()")
    )
)                

pfTausHPS_genInfo = pfTausHPS_recInfo.clone(
    pluginType = cms.string("PATTauVectorGenJetValExtractor"),

    columns = cms.PSet(
        # generator level information
        genMatch = cms.string("genMatch"),
        
        genPt = cms.string("genPt"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi"),
        
        genDecayMode = cms.string("genDecayMode")
    )
)    
