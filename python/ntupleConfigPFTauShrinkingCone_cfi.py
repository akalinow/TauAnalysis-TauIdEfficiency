import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# define variables specific to PFTaus
# reconstructed by shrinking signal cone algorithm
#--------------------------------------------------------------------------------

from RecoTauTag.TauTagTools.PFTauMVAInputDiscriminatorTranslator_cfi import \
        produceTancMVAInputDiscriminators

pfTausShrinkingCone_recInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingCone"),
    
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
        
        # additional TaNC (neural network) based discriminators
        decayMode = cms.string('decayMode()'),

        # Load decay mode kinematic information (hacked into the discriminators)
        decayModePt = cms.string("tauID('decayModePt')"),
        decayModeEta = cms.string("tauID('decayModeEta')"),
        decayModePhi = cms.string("tauID('decayModePhi')"),
        decayModeMass = cms.string("tauID('decayModeMass')"),

        byTaNC = cms.string("tauID('byTaNC')"),
        byTaNCfrOnePercent = cms.string("tauID('byTaNCfrOnePercent')"),
        byTaNCfrHalfPercent = cms.string("tauID('byTaNCfrHalfPercent')"),
        byTaNCfrQuarterPercent = cms.string("tauID('byTaNCfrQuarterPercent')"),        
        byTaNCfrTenthPercent = cms.string("tauID('byTaNCfrTenthPercent')"),
        ##transformedTaNCoutput = cms.string("tauID('transformedTaNCoutput')")
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        againstMuon = cms.string("tauID('againstMuon')")
    )
)

# Insert TaNC inputs into ntuple
print "Embedding TaNC inputs into shrinking cone ntuple"
for tancInputInfo in produceTancMVAInputDiscriminators.discriminants:
    name = "TaNC" + tancInputInfo.name.value()
    if hasattr(tancInputInfo, "indices"): 
        # multiple input
        for index in tancInputInfo.indices:
            collectionName = name + str(index)
            setattr(pfTausShrinkingCone_recInfo.columns, collectionName, cms.string(
                "tauID('%s')" % collectionName))
    else: 
        # single input
        setattr(pfTausShrinkingCone_recInfo.columns, name, cms.string(
                "tauID('%s')" % name ))

pfTausShrinkingCone_genInfo = pfTausShrinkingCone_recInfo.clone(
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
