from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: samples.py,v 1.17 2010/11/09 16:51:18 friis Exp $

'''

# Map the complicated PFTau EDProducer names to something more managable.  These
# are the strings that can/should now be used to retrieve the ntuple from the
# ntuple manager.  Other samples (W+jets etc) can rename their tau collections
# to match these.
dijetSampleAliasMap = {
    'patPFTausDijetTagAndProbeHPS': 'hps',
    'patPFTausDijetTagAndProbeShrinkingCone': 'shrinking',
    'patPFTausDijetTagAndProbeFixedCone': 'fixed',
}

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
ztautau_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_ztt",
                          "merge", datasets = ["Ztautau"],
                          take_every=40,
                          alias_map=dijetSampleAliasMap)

zllPU156bx_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_zllPU156bx",
                             "merge", datasets = ["ZllPU156bx"],
                             alias_map=dijetSampleAliasMap)
zttPU156bx_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_zttPU156bx",
                             "merge", datasets = ["ZttPU156bx"],
                             alias_map=dijetSampleAliasMap)

# Merge multiple pt hat bins
#
# CV: restrict analysis to first pt hat bin for now...
#
print "loading definition of QCD (pythia 6) background Monte Carlo samples..."
##qcd_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd", "merge", "QCD_Pt15", "QCD_Pt30", "QCD_Pt80", "QCD_Pt170")
qcd_mc_pythia6 = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcd_pythia6", "merge",
    take_every=1, datasets = ["QCD_Pt15_pythia6"],
    alias_map = dijetSampleAliasMap
)
qcd_mc_pythia6_recoTrackDowngrade = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcd_pythia6_recoTrackDowngrade",
    "merge", take_every=1, datasets = ["QCD_Pt15_pythia6_recoTrackDowngrade"],
    alias_map = dijetSampleAliasMap
)

print "loading definition of QCD (pythia 8) background Monte Carlo samples..."
qcd_mc_pythia8 = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcd_pythia8", "merge",
    take_every=1, datasets = ["QCD_Pt15_pythia8"],
    alias_map = dijetSampleAliasMap
)

qcd_mc_pythia8_recoTrackDowngrade = build_sample(
    _MC_LUMI_MAP_FILE, "mc_qcd_pythia8_recoTrackDowngrade",
    "merge", take_every=1, datasets = ["QCD_Pt15_pythia8_recoTrackDowngrade"],
    alias_map = dijetSampleAliasMap
)

print "loading definition of min. Bias (pythia 6) background Monte Carlo samples..."
minbias_mc_pythia6 = build_sample(
    _MC_LUMI_MAP_FILE, "mc_minbias_pythia6", "merge",
    take_every=1, datasets = ["minBias_pythia6"],
    alias_map = dijetSampleAliasMap
)

print "loading definition of min. Bias (pythia 8) background Monte Carlo samples..."
minbias_mc_pythia8 = build_sample(
    _MC_LUMI_MAP_FILE, "mc_minbias_pythia8",
    "merge", take_every=1, datasets = ["minBias_pythia8"])

# For data, we use the add mode, to concatenate data
print "loading definition of Data samples..."
data_runs132440to133802 = build_sample(
    _DATA_LUMI_MAP_FILE, "data_runs132440to133802", "add",
    take_every=1, datasets = ["data_runs132440_133802"],
    alias_map = dijetSampleAliasMap)
data_runs135821to141887 = build_sample(
    _DATA_LUMI_MAP_FILE, "data_runs135821to141887", "add",
    take_every=1, datasets = ["data_runs135821_141887"],
    alias_map = dijetSampleAliasMap)
data_runs141950to144114 = build_sample(
    _DATA_LUMI_MAP_FILE, "data_runs141950to144114", "add",
    take_every=1, datasets = ["data_runs141950_144114"],
    alias_map = dijetSampleAliasMap)
data = build_sample(
    _DATA_LUMI_MAP_FILE, "data", "add",
    take_every=100, datasets = ["data_runs132440_133802", "data_runs135821_141887", "data_runs141950_144114"],
    alias_map = dijetSampleAliasMap)
