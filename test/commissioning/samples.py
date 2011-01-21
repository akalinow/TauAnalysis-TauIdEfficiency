from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: samples.py,v 1.36 2011/01/21 10:33:02 veelken Exp $

'''

# Map the complicated PFTau EDProducer names to something more managable.  These
# are the strings that can/should now be used to retrieve the ntuple from the
# ntuple manager.  Other samples (W+jets etc) can rename their tau collections
# to match these.
dijetSampleAliasMap = {
    'patPFTausDijetTagAndProbeHPS': 'hps',
    'patPFTausDijetTagAndProbeShrinkingCone': 'shrinking',
    'patPFTausDijetTagAndProbeFixedCone': 'fixed',
    'patPFTausDijetTagAndProbeHPSpTaNC': 'hpstanc',
    'patCaloTausDijetTagAndProbe': 'calo'
}

muEnrichedSampleAliasMap = {
    'patPFTausCleanedHPS': 'hps',
    'patPFTausCleanedHPSpTaNC': 'hpstanc',
    'patPFTausCleanedShrinkingCone': 'shrinking',
    'patPFTausCleanedFixedCone': 'fixed',
    'patCaloTausCleaned': 'calo'
}

#wJetsSampleAliasMap = {
    #'patPFTausCleanedHPS': 'hps',
    #'patPFTausCleanedHPSpTaNC': 'hpstanc',
    #'patPFTausCleanedShrinkingCone': 'shrinking',
    #'patPFTausCleanedFixedCone': 'fixed',
    #'patCaloTausCleaned': 'calo'
#}

wJetsSampleAliasMap = {
    'patPFTausLoosePFIsoEmbeddedHPS': 'hps',
    'patPFTausLoosePFIsoEmbedded06HPSpTaNC': 'hpstanc',
    'patPFTausLoosePFIsoEmbeddedShrinkingCone': 'shrinking',
    'patPFTausLoosePFIsoEmbeddedFixedCone': 'fixed',
    'patCaloTausLoosePFIsoEmbedded': 'calo'
}

muEnrichedSampleAliasMap = wJetsSampleAliasMap

# Locations of the luminosity maps
_DATA_LUMI_MAP_FILE = os.path.join(
    os.environ['CMSSW_BASE'], 'src/TauAnalysis/TauIdEfficiency/test/'
    'commissioning', 'dataLumiMap.json')

_MC_LUMI_MAP_FILE = os.path.join(
    os.environ['CMSSW_BASE'], 'src/TauAnalysis/TauIdEfficiency/test/'
    'commissioning', 'mcLumiMap.json')

# Arugments: lumi map file, name of output collection, merge/add, list of samples
# to take from the JSON file
print "loading definition of Ztautau signal Monte Carlo samples..."
ztautau_mc = build_sample(
    _MC_LUMI_MAP_FILE, "mc_ztautau",
    "merge", datasets = ["Ztautau"],
    alias_map=dijetSampleAliasMap)
#ztautau_mc = None
zttPU156bx_mc = build_sample(
    _MC_LUMI_MAP_FILE, "mc_zttPU156bx",
    "merge", datasets = ["ZtautauPU156bx"],
    alias_map=dijetSampleAliasMap)

print "loading definition of QCD background Monte Carlo samples..."
qcddijet_mc = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcddijet", "merge",
    take_every=1, datasets = [
        "mcQCDdiJetPtHat5to15",
        "mcQCDdiJetPtHat15to30",
        "mcQCDdiJetPtHat30to50",
        "mcQCDdiJetPtHat50to80"
    ],
    alias_map = dijetSampleAliasMap)
qcddijet_mc_noCuts = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcddijet_noCuts", "merge",
    take_every=1, datasets = [
        "mcQCDdiJetPtHat5to15_noCuts",
        "mcQCDdiJetPtHat15to30_noCuts",
        "mcQCDdiJetPtHat30to50_noCuts",
        "mcQCDdiJetPtHat50to80_noCuts"
    ],
    alias_map = dijetSampleAliasMap)
qcddijetPU156bx_mc = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcddijetPU156bx", "merge",
    take_every=1, datasets = [
        "mcQCDdiJetPtHat5to15PU156bx",
        "mcQCDdiJetPtHat15to30PU156bx",
        "mcQCDdiJetPtHat30to50PU156bx",
        "mcQCDdiJetPtHat50to80PU156bx"
    ],
    alias_map = dijetSampleAliasMap)
qcddijetPU156bx_mc_noCuts = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcddijetPU156bx_noCuts", "merge",
    take_every=1, datasets = [
        "mcQCDdiJetPtHat5to15PU156bx_noCuts",
        "mcQCDdiJetPtHat15to30PU156bx_noCuts",
        "mcQCDdiJetPtHat30to50PU156bx_noCuts",
        "mcQCDdiJetPtHat50to80PU156bx_noCuts"
    ],
    alias_map = dijetSampleAliasMap)

ppmux10_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux10", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt10"],
   alias_map = muEnrichedSampleAliasMap)
ppmux10_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux10_noCuts", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt10_noCuts"],
   alias_map = muEnrichedSampleAliasMap)
ppmux10PU156bx_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux10PU156bx", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt10PU156bx"],
   alias_map = muEnrichedSampleAliasMap)
ppmux10PU156bx_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux10PU156bx_noCuts", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt10PU156bx_noCuts"],
   alias_map = muEnrichedSampleAliasMap)
ppmux15_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux15", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt15"],
   alias_map = muEnrichedSampleAliasMap)
ppmux15_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux15_noCuts", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt15_noCuts"],
   alias_map = muEnrichedSampleAliasMap)
ppmux15PU156bx_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux15PU156bx", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt15PU156bx"],
   alias_map = muEnrichedSampleAliasMap)
ppmux15PU156bx_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux15PU156bx_noCuts", "merge",
   take_every=1, datasets = ["mcQCDppMuXptHatGt20PtMuGt15PU156bx_noCuts"],
   alias_map = muEnrichedSampleAliasMap)

print "loading definition of W + jets background Monte Carlo samples..."
wjets_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_wjets", "merge",
   take_every=1, datasets = ["mcWplusJets"],
   alias_map = wJetsSampleAliasMap)
wjets_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_wjets_noCuts", "merge",
   take_every=1, datasets = ["mcWplusJets_noCuts"],
   alias_map = wJetsSampleAliasMap)
wjetsPU156bx_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_wjetsPU156bx", "merge",
   take_every=1, datasets = ["mcWplusJetsPU156bx"],
   alias_map = wJetsSampleAliasMap)
wjetsPU156bx_mc_noCuts = build_sample(
   _MC_LUMI_MAP_FILE, "mc_wjetsPU156bx_noCuts", "merge",
   take_every=1, datasets = ["mcWplusJetsPU156bx_noCuts"],
   alias_map = wJetsSampleAliasMap)

# For data, we use the add mode, to concatenate data
print "loading definition of Data samples..."
data_dijet = build_sample(
    _DATA_LUMI_MAP_FILE, "qcdDiJet_data", "add",
    take_every=1, datasets = [
        "qcdDiJet_data_MinBiasCommissioning_Jun14ReReco",
        #"qcdDiJet_data_JetMETTau_Run2010Ai_Sep17ReReco",
        #"qcdDiJet_data_JetMET_Run2010Aii_Sep17ReReco",
        #"qcdDiJet_data_Jet_Run2010B_Nov4ReReco"
    ],
    alias_map = dijetSampleAliasMap)

data_ppmux = build_sample(
    _DATA_LUMI_MAP_FILE, "qcdMuEnriched_data", "add",
    take_every=1, datasets = [
        "qcdMuEnriched_data_Mu_Run2010A_Nov4ReReco",
        "qcdMuEnriched_data_Mu_Run2010B_Nov4ReReco"
    ],
    alias_map = muEnrichedSampleAliasMap)

data_wjets = build_sample(
    _DATA_LUMI_MAP_FILE, "wJets_data", "add",
    take_every=1, datasets = [
        "wJets_data_Mu_Run2010A_Nov4ReReco",
        "wJets_data_Mu_Run2010B_Nov4ReReco"
    ],
    alias_map = wJetsSampleAliasMap)
