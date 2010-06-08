from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: samples.py,v 1.1 2010/05/27 21:10:23 friis Exp $

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
ztautau_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_ztt", "merge", "Ztautau")

# Merge multiple pt hat bins
qcd_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd", "merge", "QCD_Pt15", "QCD_Pt30", "QCD_Pt80", "QCD_Pt170")

minbias_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_minbias", "merge", "minBias")

# For data, we use the add mode, to concatenate data
data = build_sample(_DATA_LUMI_MAP_FILE, "data", "add", "Data_part01")
