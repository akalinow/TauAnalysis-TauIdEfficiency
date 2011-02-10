import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to CaloTaus
#--------------------------------------------------------------------------------

caloTaus_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patCaloTausDijetTagAndProbe"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for CaloTau
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),
        mass = cms.string("mass()"),

        # charge of CaloTau
        # (sum of charges of tracks within signal cone)
        charge = cms.string("charge()"),

        # kinematic variables for "Jet + Tracks" (JPT) Jet associated to CaloTau
        # (NOTE: energy of JPTJet corresponds to corrected jet energy)
        jetPt = cms.string("caloTauTagInfoRef().jetRef().pt()"),
        jetEta = cms.string("caloTauTagInfoRef().jetRef().eta()"),
        jetPhi = cms.string("caloTauTagInfoRef().jetRef().phi()"),

        # kinematic variables for CaloJet associated to CaloTau
        # (NOTE: energy of CaloJet corresponds to uncorrected/"raw" calorimeter energy)
        caloJetPt = cms.string("caloTauTagInfoRef().calojetRef().pt()"),
        caloJetEta = cms.string("caloTauTagInfoRef().calojetRef().eta()"),
        caloJetPhi = cms.string("caloTauTagInfoRef().calojetRef().phi()"),

        # loose PFIsolation Pt sum
        # (computed by summing PFChargedHadrons of Pt > 1.0 GeV + PFGammas of Pt > 1.5 GeV
        #  in pfNoPileUp collection within region 0.15 < dR < 0.6 centered on PFTau direction)
        ptSumLooseIsolation04 = cms.string("userFloat('pfLooseIsoPt04')"),
        ptSumLooseIsolation06 = cms.string("userFloat('pfLooseIsoPt06')"),

        # momentum of leading charged/neutral particle within signal cone
        leadTrackPt = cms.string("? leadTrack().isNonnull() ? leadTrack().pt() : 0."),

        # invariant mass of tracks within signal cone
        invariantMassSignalTracks = cms.string("signalTracksInvariantMass()"),

        # multiplicity of tracks within signal/isolations cones
        numSignalTracks = cms.string("signalTracks().size()"),
        numIsolationTracks = cms.string("isolationTracks().size()"),
        numTracks = cms.string("signalTracks().size() + isolationTracks().size()"),

        # Pt sum of tracks/Et sum of ECAL energy deposits within isolation cone
        ptSumIsolationTracks = cms.string("isolationTracksPtSum()"),
        etSumIsolationECAL = cms.string("isolationECALhitsEtSum()"),

        # flags for tag/probe and Pt_index for distinguishing between
        # highest Pt, second highest Pt and third highest Pt jet
        tagJet15U = cms.string("userFloat('tagJet15U')"),
        probeJet15U = cms.string("userFloat('probeJet15U')"),
        tagJet30U = cms.string("userFloat('tagJet30U')"),
        probeJet30U = cms.string("userFloat('probeJet30U')"),
        tagJet50U = cms.string("userFloat('tagJet50U')"),
        probeJet50U = cms.string("userFloat('probeJet50U')"),
        tagJet70U = cms.string("userFloat('tagJet70U')"),
        probeJet70U = cms.string("userFloat('probeJet70U')"),
        tagJet100U = cms.string("userFloat('tagJet100U')"),
        probeJet100U = cms.string("userFloat('probeJet100U')"),
        ptIndex = cms.string("userFloat('pt_index')"),
        
        # tau id. discriminators based on leading track
        byLeadTrackFinding = cms.string("tauID('leadingTrackFinding')"),
        byLeadTrackPtCut = cms.string("tauID('leadingTrackPtCut')"),
        byTrackIsolation = cms.string("tauID('trackIsolation')"),
        byEcalIsolation = cms.string("tauID('ecalIsolation')"),
        byIsolation = cms.string("tauID('byIsolation')"),
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')")
    )
)

caloTaus_recJetIdInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorCaloJetIdValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patCaloTausDijetTagAndProbe"),
    srcJet = cms.InputTag("patJetsAK5Calo"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # jetId bits
        jetIdMinimal = cms.string("minimal"),
        jetIdLooseAOD = cms.string("loose_AOD"),
        jetIdLoose = cms.string("loose"),
        jetIdTight = cms.string("tight")
    )
)   

caloTaus_genInfo = caloTaus_recInfo.clone(
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

caloTaus_mcEmbeddingInfo = caloTaus_recInfo.clone(
    pluginType = cms.string("PATTauVectorValExtractor"),

    columns = cms.PSet(
        # flags for tag/probe
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
    )
)    

