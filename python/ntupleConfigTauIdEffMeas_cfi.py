import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# Muon quantities
#--------------------------------------------------------------------------------

tauIdEffMeas_template01 = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False),
    indices = cms.vuint32([0]), # Store values for first object only
        
    # Extractor plugin
    pluginType = cms.string("PATMuonValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),

    # Variables to compute for this source
    columns = cms.PSet(
        # kinematic variables for Muon
        pt = cms.string("pt()"),
        eta = cms.string("eta()"),
        phi = cms.string("phi()"),

        # loose isolation
        ptSumLooseIsolation03 = cms.string("userFloat('pfLooseIsoPt03')"),
        ptSumLooseIsolation04 = cms.string("userFloat('pfLooseIsoPt04')"),
        ptSumLooseIsolation06 = cms.string("userFloat('pfLooseIsoPt06')")
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
    pluginType = cms.string("PATMuTauPairValExtractor"),
    
    # Collection to extract from
    src = cms.InputTag(""),

    # Variables to compute for this source
    columns = cms.PSet(
        chargeTauLeadTrack = cms.string(""),
        charge = cms.string(""),

	Mt = cms.string("mt1MET()"),  
        pZeta = cms.string("pZeta()"),
        pZetaVis = cms.string("pZetaVis()"),

        visMass = cms.string("p4Vis.mass()"),
        Ht = cms.string("leg1.pt() + leg2.pt() + met.pt()")
    )
)

tauIdEffMeas_template03pfTau = copy.deepcopy(tauIdEffMeas_template03)
tauIdEffMeas_template03pfTau.columns.chargeTauLeadTrack = \
  cms.string("? leg2.leadPFChargedHadrCand().isNonnull() ? (leg1.charge() + leg2.leadPFChargedHadrCand.charge()) : -1000.")
tauIdEffMeas_template03pfTau.columns.charge = \
  cms.string("leg1.charge() + leg2.charge()")

#--------------------------------------------------------------------------------
# Single muon Trigger efficiency
#--------------------------------------------------------------------------------

tauIdEffMeas_muonTriggerEff = cms.PSet(
    # Select multiplicy of object(s) to store
    vector = cms.bool(False),
    indices = cms.vuint32([0]), # Store values for first object only
    
    # Extractor plugin
    pluginType = cms.string("KineEventReweightExtractor"),

    # Collection to extract from
    src = cms.InputTag(""),

    # Variables to compute for this source
    columns = cms.PSet(
        muonHLTeff = cms.PSet(    
            weightLookupTable = cms.PSet(
                inputFileName = cms.FileInPath('TauAnalysis/TauIdEfficiency/data/singleMuHLTeff.root'),
                lutName = cms.string("hEff")
            ),
    
            variables = cms.PSet(
                x = cms.PSet(
                    pluginType = cms.string("PATMuonValExtractor"),
                    src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
                    value = cms.string("pt()"),
                    vector = cms.bool(False),
                    indices = cms.vuint32([0])
                ),
                y = cms.PSet(
                    pluginType = cms.string("PATMuonValExtractor"),
                    src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
                    value = cms.string("eta()"),
                    vector = cms.bool(False),
                    indices = cms.vuint32([0])
                ),
                z = cms.PSet(
                    pluginType = cms.string("PATMuonValExtractor"),
                    src = cms.InputTag('selectedPatMuonsForTauIdEffTrkIPcumulative'),
                    value = cms.string("? pt() > 0. ? userFloat('pfLooseIsoPt03')/pt() : -1."),
                    vector = cms.bool(False),
                    indices = cms.vuint32([0])
                )
            )
        )
    )
)
