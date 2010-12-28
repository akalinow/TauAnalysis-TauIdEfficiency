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
        mass = cms.string("mass()"),

        # kinematic variables for PFJet associated to PFTau
        jetPt = cms.string("pfTauTagInfoRef().pfjetRef().pt()"),
        jetEta = cms.string("pfTauTagInfoRef().pfjetRef().eta()"),
        jetPhi = cms.string("pfTauTagInfoRef().pfjetRef().phi()"),

        # loose PFIsolation Pt sum
        # (computed by summing PFChargedHadrons of Pt > 1.0 GeV + PFGammas of Pt > 1.5 GeV
        #  in pfNoPileUp collection within region 0.15 < dR < 0.4/0.6 centered on PFTau direction)
        ptSumLooseIsolation04 = cms.string("userFloat('pfLooseIsoPt04')"),
        ptSumLooseIsolation06 = cms.string("userFloat('pfLooseIsoPt06')"),

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
        pfElectronMVA = cms.string("? leadPFCand().isNonnull() ? leadPFCand().mva_e_pi() : +1."),
        againstMuon = cms.string("tauID('againstMuon')"),                             

        # multiplicity and Pt sum of charged/neutral particles within signal and isolation cones
        numChargedParticlesSignalCone = cms.string("signalPFChargedHadrCands().size()"),
        numNeutralHadronsSignalCone = cms.string("signalPFNeutrHadrCands().size()"),
        numPhotonsSignalCone = cms.string("signalPFGammaCands().size()"),
        numParticlesSignalCone = cms.string("signalPFCands().size()"),

        numChargedParticlesIsoCone = cms.string("isolationPFChargedHadrCands().size()"),
        numNeutralHadronsIsoCone = cms.string("isolationPFNeutrHadrCands().size()"),
        numPhotonsIsoCone = cms.string("isolationPFGammaCands().size()"),
        numParticlesIsoCone = cms.string("isolationPFCands().size()"),

        ptSumChargedParticlesIsoCone = cms.string("isolationPFChargedHadrCandsPtSum()"),
        ptSumPhotonsIsoCone = cms.string("isolationPFGammaCandsEtSum()"),
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
        
        genDecayMode = cms.string("genDecayMode"),
        
        genPdgId = cms.string("genPdgId")
    ),
    
    srcGenParticles = cms.InputTag('genParticles'),
    skipPdgIdsGenParticleMatch = cms.vint32( 12, 14, 16 )
)

pfTausHPS_mcEmbeddingInfo = pfTausHPS_recInfo.clone(
    pluginType = cms.string("PATTauVectorValExtractor"),

    columns = cms.PSet(
        # flags for tag/probe
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
    )
)    

