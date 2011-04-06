from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: samples.py,v 1.41 2011/03/10 11:27:41 veelken Exp $

'''

# Map the complicated PFTau EDProducer names to something more managable.  These
# are the strings that can/should now be used to retrieve the ntuple from the
# ntuple manager.  Other samples (W+jets etc) can rename their tau collections
# to match these.
dijetSampleAliasMap = {
    'patPFTausDijetTagAndProbeHPS':               'hps',
    'patPFTausDijetTagAndProbeShrinkingCone':     'shrinking',
    'patPFTausDijetTagAndProbeFixedCone':         'fixed',
    'patPFTausDijetTagAndProbeHPSpTaNC':          'hpstanc',
    'patCaloTausDijetTagAndProbe':                'calo'
}

wJetsSampleAliasMap = {
    'patPFTausLoosePFIsoEmbedded06HPS':           'hps',
    'patPFTausLoosePFIsoEmbedded06HPSpTaNC':      'hpstanc',
    'patPFTausLoosePFIsoEmbedded06ShrinkingCone': 'shrinking',
    'patPFTausLoosePFIsoEmbedded06FixedCone':     'fixed',
    'patCaloTausLoosePFIsoEmbedded06':            'calo'
}

muEnrichedSampleAliasMap = wJetsSampleAliasMap

#take_every = 1
take_every = 100 # NOTE: to be used for debugging purposes only !!

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

print "loading definition of QCD background Monte Carlo samples..."
qcddijet_mc = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcddijet", "merge",
    take_every = take_every, datasets = [
        "mcQCDdiJetPtHat15to30",
        "mcQCDdiJetPtHat30to50",
        "mcQCDdiJetPtHat50to80",
        "mcQCDdiJetPtHat80to120",
        "mcQCDdiJetPtHat120to170",
        "mcQCDdiJetPtHat170to300"
    ],
    alias_map = dijetSampleAliasMap)

ppmux15_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_ppmux15", "merge",
   take_every = take_every, datasets = ["mcQCDppMuXptHatGt20PtMuGt15"],
   alias_map = muEnrichedSampleAliasMap)

print "loading definition of W + jets background Monte Carlo samples..."
wjets_mc = build_sample(
   _MC_LUMI_MAP_FILE, "mc_wjets", "merge",
   take_every = take_every, datasets = ["mcWplusJets"],
   alias_map = wJetsSampleAliasMap)

# For data, we use the add mode, to concatenate data
print "loading definition of Data samples..."
data_dijet = build_sample(
    _DATA_LUMI_MAP_FILE, "qcdDiJet_data", "add",
    take_every = take_every, datasets = [
        "qcdDiJet_data_Jet_Run2011A_PromptReco"
    ],
    alias_map = dijetSampleAliasMap)
data_ppmux = build_sample(
    _DATA_LUMI_MAP_FILE, "qcdMuEnriched_data", "add",
    take_every = take_every, datasets = [
        "qcdMuEnriched_data_Mu_Run2011A_PromptReco"
    ],
    alias_map = muEnrichedSampleAliasMap)
data_wjets = build_sample(
    _DATA_LUMI_MAP_FILE, "wJets_data", "add",
    take_every = take_every, datasets = [
        "wJets_data_Mu_Run2011A_PromptReco",
    ],
    alias_map = wJetsSampleAliasMap)
