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

        # charge of PFTau
        # (sum of charges of charged particles within signal cone)
        charge = cms.string("charge()"),

        # transverse momentum of leading charged/neutral particle within signal cone
        leadChargedParticlePt = cms.string("? leadPFChargedHadrCand().isNonnull() ? leadPFChargedHadrCand().pt() : 0."),
        leadNeutralParticlePt = cms.string("? leadPFNeutralCand().isNonnull() ? leadPFNeutralCand().pt() : 0."),
        leadParticlePt = cms.string("? leadPFCand().isNonnull() ? leadPFCand().pt() : 0."),

        # multiplicity and Pt sum of charged/neutral particles within signal and isolation cones
        numChargedParticlesSignalCone = cms.string("signalPFChargedHadrCands().size()"),
        numNeutralHadronsSignalCone = cms.string("signalPFNeutrHadrCands().size()"),
        numPhotonsSignalCone = cms.string("signalPFGammaCands().size()"),
        numParticlesSignalCone = cms.string("signalPFCands().size()"),

        numChargedParticlesIsoCone = cms.string("isolationPFChargedHadrCands().size()"),
        numNeutralHadronsIsoCone = cms.string("isolationPFNeutrHadrCands().size()"),
        numPhotonsIsoCone = cms.string("isolationPFGammaCands().size()"),
        numParticlesIsoCone = cms.string("isolationPFCands().size()"),

        numChargedParticles = cms.string("signalPFChargedHadrCands().size() + isolationPFChargedHadrCands().size()"),
        numNeutralHadrons = cms.string("signalPFNeutrHadrCands().size() + isolationPFNeutrHadrCands().size()"),
        numPhotons = cms.string("signalPFGammaCands().size() + isolationPFGammaCands().size()"),
        numParticles = cms.string("signalPFCands().size() + isolationPFCands().size()"),

        ptSumChargedParticlesIsoCone = cms.string("isolationPFChargedHadrCandsPtSum()"),
        ptSumPhotonsIsoCone = cms.string("isolationPFGammaCandsEtSum()"),

        # kinematic variables for PFJet associated to PFTau
        jetPt = cms.string("pfJetRef().pt()"),
        jetEta = cms.string("pfJetRef().eta()"),
        jetPhi = cms.string("pfJetRef().phi()"),

	# jet width
        jetWidth = cms.string("? (etaetaMoment() + phiphiMoment()) > 0. ? sqrt(etaetaMoment() + phiphiMoment()) : 0."),

        # loose PFIsolation Pt sum
        # (computed by summing PFChargedHadrons of Pt > 1.0 GeV + PFGammas of Pt > 1.5 GeV
        #  in pfNoPileUp collection within region 0.15 < dR < 0.4/0.6 centered on PFTau direction)
        ptSumLooseIsolation04 = cms.string("userFloat('pfLooseIsoPt04')"),
        ptSumLooseIsolation06 = cms.string("userFloat('pfLooseIsoPt06')"),

        # flags for tag/probe and Pt_index for distinguishing between
        # highest Pt, second highest Pt and third highest Pt jet
        tagJet30v1 = cms.string("userFloat('tagJet30v1')"),
        probeJet30v1 = cms.string("userFloat('probeJet30v1')"),
        tagJet60v1 = cms.string("userFloat('tagJet60v1')"),
        probeJet60v1 = cms.string("userFloat('probeJet60v1')"),
        tagJet80v1 = cms.string("userFloat('tagJet80v1')"),
        probeJet80v1 = cms.string("userFloat('probeJet80v1')"),
        tagJet110v1 = cms.string("userFloat('tagJet110v1')"),
        probeJet110v1 = cms.string("userFloat('probeJet110v1')"),
        tagJet150v1 = cms.string("userFloat('tagJet150v1')"),
        probeJet150v1 = cms.string("userFloat('probeJet150v1')"),
        
        # tau id. discriminators        
        byDecayModeFinding = cms.string("tauID('decayModeFinding')"),
        byIsolationLoose = cms.string("tauID('byLooseIsolation')"),
        byIsolationMedium = cms.string("tauID('byMediumIsolation')"),
        byIsolationTight = cms.string("tauID('byTightIsolation')"),
        
        # discriminators against electrons/muons
        againstElectronLoose = cms.string("tauID('againstElectronLoose')"),
        againstElectronMedium = cms.string("tauID('againstElectronMedium')"),
        againstElectronTight = cms.string("tauID('againstElectronTight')"),
        pfElectronMVA = cms.string("? leadPFCand().isNonnull() ? leadPFCand().mva_e_pi() : +1."),
        againstMuonLoose = cms.string("tauID('againstMuonLoose')"),
        againstMuonTight = cms.string("tauID('againstMuonTight')")
    )
)

pfTausHPS_recJetIdInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorPFJetIdValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeHPS"),
    srcJet = cms.InputTag("patJetsAK5PF"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # jetId bits
        jetIdLoose = cms.string("loose"),
        jetIdTight = cms.string("tight")
    )
)

pfTausHPS_recJetCorrMomInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorJetCorrMomValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeHPS"),
    srcJet = cms.InputTag("patJetsAK5PF"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for PFJet associated to PFTau,
        # with L2+L3 (MC) / L2+L3+Residual (Data) jet energy corrections applied
        jetCorrPt = cms.string("pt()"),
        jetCorrEta = cms.string("eta()"),
        jetCorrPhi = cms.string("phi()")
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

