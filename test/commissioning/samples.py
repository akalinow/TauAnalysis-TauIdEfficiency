from TauAnalysis.TauIdEfficiency.ntauples.sample_builder import build_sample
import os

'''

samples.py

Central defintion of data sources for commissioning.

$Id: $

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
ztautau_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_ztt", "merge", "ztautau")

# Merge multiple pt hat bins
qcd_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_qcd", "merge", "QCD_20_30", 
                      "QCD_30_50")

minbias_mc = build_sample(_MC_LUMI_MAP_FILE, "mc_minbias", "merge", "minbias")

# For data, we use the add mode, to concatenate data
data_first_two_weeks = build_sample(_DATA_LUMI_MAP_FILE, "data_two_weeks", 
                                    "add", "data_week_1", "data_week_2")

data_first_three_weeks = build_sample(_DATA_LUMI_MAP_FILE, "data_three_weeks", 
                                      "add", "data_week_1", "data_week_2", 
                                      "data_week_3")
