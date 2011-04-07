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

        # tau decay mode
        decayMode = cms.string('decayMode()'),
        
        # additional TaNC (neural network) based discriminators
        byTaNC = cms.string("tauID('byTaNC')"),
        byTaNCfrOnePercent = cms.string("tauID('byTaNCfrOnePercent')"),
        byTaNCfrHalfPercent = cms.string("tauID('byTaNCfrHalfPercent')"),
        byTaNCfrQuarterPercent = cms.string("tauID('byTaNCfrQuarterPercent')"),        
        byTaNCfrTenthPercent = cms.string("tauID('byTaNCfrTenthPercent')"),
        
        # discriminators against electrons/muons
        againstElectron = cms.string("tauID('againstElectron')"),
        pfElectronMVA = cms.string("? leadPFCand().isNonnull() ? leadPFCand().mva_e_pi() : +1."),
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

pfTausShrinkingCone_recJetIdInfo = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(True), # Store a value for all objects in this collection
    #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects
    
    # Extractor plugin
    pluginType = cms.string("PATTauVectorPFJetIdValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag("patPFTausDijetTagAndProbeShrinkingCone"),
    srcJet = cms.InputTag("patJetsAK5PF"),
    
    # Variables to compute for this source
    columns = cms.PSet(
        # jetId bits
        jetIdLoose = cms.string("loose"),
        jetIdTight = cms.string("tight")
    )
)   

pfTausShrinkingCone_genInfo = pfTausShrinkingCone_recInfo.clone(
    pluginType = cms.string("PATTauVectorGenJetValExtractor"),

    columns = cms.PSet(
        # generator level information
        genMatch = cms.string("genMatch"),
        
        genPt = cms.string("genPt"),
        genMass = cms.string("genMass"),
        genEta = cms.string("genEta"),
        genPhi = cms.string("genPhi"),

        genDecayMode = cms.string("genDecayMode"),
        
        genPdgId = cms.string("genPdgId")
    ),
    
    srcGenParticles = cms.InputTag('genParticles'),
    skipPdgIdsGenParticleMatch = cms.vint32( 12, 14, 16 )
)    

pfTausShrinkingCone_mcEmbeddingInfo = pfTausShrinkingCone_recInfo.clone(
    pluginType = cms.string("PATTauVectorValExtractor"),

    columns = cms.PSet(
        # flags for tag/probe
        tag = cms.string("userFloat('tag')"),
        probe = cms.string("userFloat('probe')"),
    )
)    
