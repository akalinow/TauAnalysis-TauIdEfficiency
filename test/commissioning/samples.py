from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: samples.py,v 1.9 2010/06/26 01:15:44 friis Exp $

'''

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
ztautau_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_ztt", "merge", datasets = ["Ztautau"])

# Merge multiple pt hat bins
#
# CV: restrict analysis to first pt hat bin for now...
#
print "loading definition of QCD (pythia 6) background Monte Carlo samples..."
##qcd_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd", "merge", "QCD_Pt15", "QCD_Pt30", "QCD_Pt80", "QCD_Pt170")
qcd_mc_pythia6 = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd_pythia6", "merge", take_every=1, datasets = ["QCD_Pt15_pythia6"])

print "loading definition of QCD (pythia 8) background Monte Carlo samples..."
qcd_mc_pythia8 = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd_pythia8", "merge", take_every=1, datasets = ["QCD_Pt15_pythia8"])

print "loading definition of min. Bias (pythia 6) background Monte Carlo samples..."
minbias_mc_pythia6 = build_sample(_MC_LUMI_MAP_FILE, "mc_minbias_pythia6", "merge", take_every=1, datasets = ["minBias_pythia6"])

print "loading definition of min. Bias (pythia 8) background Monte Carlo samples..."
minbias_mc_pythia8 = build_sample(_MC_LUMI_MAP_FILE, "mc_minbias_pythia8", "merge", take_every=1, datasets = ["minBias_pythia8"])

# For data, we use the add mode, to concatenate data
print "loading definition of Data samples..."
data = build_sample(_DATA_LUMI_MAP_FILE, "data", "add", take_every=1, datasets = ["Data_rerecoMay27th"])
