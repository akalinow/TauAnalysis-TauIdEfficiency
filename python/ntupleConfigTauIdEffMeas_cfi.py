import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# Muon quantities
#--------------------------------------------------------------------------------

tauIdEffMeas_template01 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False),
    indices = cms.vuint32([0]), # Store values for first object only
        
    # Extractor plugin
    pluginType = cms.string("PATMuonVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),

    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for Muon
        pt = cms.string("leg1.pt()"),
        eta = cms.string("leg1.eta()"),
        phi = cms.string("leg1.phi()"),

        # loose isolation
        ptSumLooseIsolation = cms.string("leg2.userFloat('pfLooseIsoPt')")
    )
)

#--------------------------------------------------------------------------------
# PFTau quantities
#--------------------------------------------------------------------------------

tauIdEffMeas_template02 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False),
    indices = cms.vuint32([0]), # Store values for first object only
        
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for PFTau
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),

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
        jetPt = cms.string("pfTauTagInfoRef().pfjetRef().pt()"),
        jetEta = cms.string("pfTauTagInfoRef().pfjetRef().eta()"),
        jetPhi = cms.string("pfTauTagInfoRef().pfjetRef().phi()"),

        # jet width
        jetWidth = cms.string("sqrt(etaetaMoment() + phiphiMoment())"),

        # loose PFIsolation Pt sum
        # (computed by summing PFChargedHadrons of Pt > 1.0 GeV + PFGammas of Pt > 1.5 GeV
        #  in pfNoPileUp collection within region 0.15 < dR < 0.6 centered on PFTau direction)
        ptSumLooseIsolation = cms.string("userFloat('pfLooseIsoPt')"),

        # tau id. discriminators based on leading track/PFChargedHadron                                                 
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byLeadTrackPtCut = cms.string("tauID('leadingTrackPtCut')"),
        
        # tau id. discriminators based on leading pion (PFChargedHadron or PFGamma)                                    
        byLeadPionPtCut = cms.string("tauID('leadingPionPtCut')"),

        # tau decay mode
        decayMode = cms.string('decayMode()'),
        
        # additional TaNC (neural network) based discriminators
        byTaNCraw = cms.string("tauID('byTaNCraw')"),
        byTaNC = cms.string("tauID('byTaNC')"),
        byTaNCloose = cms.string("tauID('byTaNCloose')"),
        byTaNCmedium = cms.string("tauID('byTaNCmedium')"),
        byTaNCtight = cms.string("tauID('byTaNCtight')"),        

        # HPS discriminators
        byDecayMode = cms.string("tauID('byDecayMode')"),
        byHPSloose = cms.string("tauID('byHPSloose')"),
        byHPSmedium = cms.string("tauID('byHPSmedium')"),
        byHPStight = cms.string("tauID('byHPStight')"),
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')"),
        againstCaloMuon = cms.string("tauID('againstCaloMuon')")
    )
)

#--------------------------------------------------------------------------------
# Muon + PFTau quantities (extracted from diTau object)
#--------------------------------------------------------------------------------

tauIdEffMeas_template03 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False),
    indices = cms.vuint32([0]), # Store values for first object only
    
    # Extractor plugin
    pluginType = cms.string("PATMuTauPairVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),

    # Variables to compute for this source
    columns = cms.PSet(
        charge = cms.string("leg1.charge() + leg2.leadPFChargedHadrCand.charge()"),
        
        visMass = cms.string("p4Vis.mass()"),
        SVfitMass1 = cms.string("svFitSolution('psKine_MEt').mass()"),
        SVfitMass2 = cms.string("svFitSolution('psKine_MEt_ptBalance').mass()"),
        Mt = cms.string("mt1MET()"),
        Ht = cms.string("leg1.pt() + leg2.pt() + met.pt()")
    )
)
