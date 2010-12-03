import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to PFTaus
# reconstructed by shrinking signal cone algorithm using ellipse for photon isolation
#--------------------------------------------------------------------------------

pfTausShrinkingConeEllipticPhotonIso_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingConeEllipticPhotonIso"),
    
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

        ptSumChargedParticlesIsoCone = cms.string("isolationPFChargedHadrCandsPtSum()"),
        ptSumPhotonsIsoCone = cms.string("isolationPFGammaCandsEtSum()"),

        # kinematic variables for PFJet associated to PFTau
        jetPt = cms.string("pfTauTagInfoRef().pfjetRef().pt()"),
        jetEta = cms.string("pfTauTagInfoRef().pfjetRef().eta()"),
        jetPhi = cms.string("pfTauTagInfoRef().pfjetRef().phi()"),

        # jet width
        jetWidth = cms.string("sqrt(etaetaMoment() + phiphiMoment())"),

        # loose PFIsolation Pt sum
        # (computed by summing PFChargedHadrons of Pt > 1.0 GeV + PFGammas of Pt > 1.5 GeV
        #  in pfNoPileUp collection within region 0.15 < dR < 0.6 centered on PFTau direction)
        ptSumLooseIsolation = cms.string("userFloat('pfLooseIsoPt')"),

        # flags for tag/probe and Pt_index for distinguishing between
        # highest Pt, second highest Pt and third highest Pt jet
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
        ptIndex = cms.string("userFloat('pt_index')"),

        # tau id. discriminators based on leading track/PFChargedHadron                                                 
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byLeadTrackPtCut = cms.string("tauID('leadingTrackPtCut')"),
        byTrackIsolation = cms.string("tauID('trackIsolation')"),
        byEcalIsolation = cms.string("tauID('ecalIsolation')"),
        byIsolation = cms.string("tauID('byIsolation')"),
        
        # tau id. discriminators based on leading pion (PFChargedHadron or PFGamma)                                    
        byLeadPionPtCut = cms.string("tauID('leadingPionPtCut')"),
        byTrackIsolationUsingLeadingPion = cms.string("tauID('trackIsolationUsingLeadingPion')"),
        byEcalIsolationUsingLeadingPion = cms.string("tauID('ecalIsolationUsingLeadingPion')"),
        byIsolationUsingLeadingPion = cms.string("tauID('byIsolationUsingLeadingPion')"),
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')")
    )
)

pfTausShrinkingConeEllipticPhotonIso_genInfo = pfTausShrinkingConeEllipticPhotonIso_recInfo.clone(
    pluginType = cms.string("PATTauVectorGenJetValExtractor"),

    columns = cms.PSet(
        # generator level information
        genMatch = cms.string("genMatch"),
        
        genPt = cms.string("genPt"),
        genMass = cms.string("genMass"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi"),
        
        genDecayMode = cms.string("genDecayMode")
    )
)

pfTausShrinkingConeEllipticPhotonIso_mcEmbeddingInfo = pfTausShrinkingConeEllipticPhotonIso_recInfo.clone(
    pluginType = cms.string("PATTauVectorValExtractor"),

    columns = cms.PSet(
        # flags for tag/probe
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
    )
)    

