#!/usr/bin/env python

import os

hltPaths = [
    #'HLT_L1SingleMuOpen_v*',
    #'HLT_L1SingleMuOpen_DT_v*',
    #'HLT_L1SingleMu10_v*',
    #'HLT_IsoMu15_v*',
    #'HLT_IsoMu12_v*',
    'HLT_Mu15_v*',    # for QCD muon enriched 
    #'HLT_IsoMu17_v*', # for W --> mu nu, Z --> mu+ mu- and Tau id. efficiency
    'HLT_IsoMu24_v*',
    #'HLT_DoubleMu7_v*',
    #'HLT_Mu13_Mu8_v*',
    #'HLT_Mu17_Mu8_v*',
    'HLT_Jet30_v*',   # for QCD multi-jet
    'HLT_Jet60_v*',
    'HLT_Jet80_v*',
    'HLT_Jet110_v*',
    'HLT_Jet150_v*'
    #'HLT_Mu8_Ele17_CaloIdL_v*',
    #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
    #'HLT_Mu17_Ele8_CaloIdL_v*',
    #'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*'
    #'HLT_IsoMu15_L1ETM20_v*', # for Tau id. efficiency in RunB (NOTE: events are in MET dataset)
    #'HLT_Mu15_L1ETM20_v*'     #   -""-
    #'HLT_IsoMu15_LooseIsoPFTau15_v*',
    #'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v*'
    #'HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
    #'HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
    #'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*'
    #'HLT_Ele27_WP80_PFMT50_v*'
]

jsonFiles = {
    #'2011RunA'          : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                     + 'Cert_2011RunA_7TeV_PromptReco_Collisions11_JSON.txt',
    #'2011RunB'          : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                     + 'Cert_2011RunB_7TeV_PromptReco_Collisions11_JSON.txt',
    'runs160404-180252' : '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/' \
                         + 'Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt',
    #'May10ReReco-v1'     : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                      + 'Cert_May10ReReco-v1_JSON.txt',
    #'PromptReco-v4'      : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                      + 'Cert_PromptReco-v4_JSON.txt',
    #'05Aug2011-v1'       : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                      + 'Cert_05Aug2011-v1_JSON.txt',
    #'PromptReco-v6'      : '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/TauAnalysis/TauIdEfficiency/test/commissioning/' \
    #                      + 'Cert_PromptReco-v6_JSON.txt'
}

executable_lumiCalc = '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/RecoLuminosity/LumiDB/scripts/lumiCalc2.py'
executable_analyzeLumiCalcOutput = './analyzeLumiCalcOutput.py'

MakeFile_dict = {}
for hltPath in hltPaths:
    MakeFile_dict[hltPath] = {}
    for jsonFileName, jsonFile in jsonFiles.items():
        MakeFile_dict[hltPath][jsonFileName] = {}
        hltPath_abbr = hltPath
        hltPath_abbr = hltPath_abbr.replace("_v*", "")
        hltPath_abbr = hltPath_abbr.replace("*", "")    
        MakeFile_dict[hltPath][jsonFileName]['outputFileName'] = 'lumiCalc_%s_%s.out' % (hltPath_abbr, jsonFileName)
        MakeFile_dict[hltPath][jsonFileName]['lumiCalc_command'] = \
          '%s recorded -hltpath "%s" -i %s' % (executable_lumiCalc, hltPath, jsonFile)

def make_MakeFile_vstring(list_of_strings):
    retVal = ""
    for i, string_i in enumerate(list_of_strings):
        if i > 0:
            retVal += " "
        retVal += string_i
    return retVal
    
makeFileName = "Makefile_runLumiCalc"
makeFile = open(makeFileName, "w")
makeFile.write("\n")
outputFileNames = []
for hltPath in hltPaths:
    for jsonFileName in jsonFiles.keys():
        outputFileNames.append(MakeFile_dict[hltPath][jsonFileName]['outputFileName'])
makeFile.write("all: parseLumiCalcOutput\n")
makeFile.write("\n")
for hltPath in hltPaths:
    for jsonFileName in jsonFiles.keys():
        makeFile.write("%s:\n" % MakeFile_dict[hltPath][jsonFileName]['outputFileName'])
        makeFile.write("\t%s &> %s\n" %
          (MakeFile_dict[hltPath][jsonFileName]['lumiCalc_command'],
           MakeFile_dict[hltPath][jsonFileName]['outputFileName']))
makeFile.write("\n")    
makeFile.write("parseLumiCalcOutput: %s\n" % make_MakeFile_vstring(outputFileNames))
for hltPath in hltPaths:
    for jsonFileName in jsonFiles.keys():
        makeFile.write("\t%s %s\n" %
          (executable_analyzeLumiCalcOutput,
           MakeFile_dict[hltPath][jsonFileName]['outputFileName']))
makeFile.write("\techo 'Finished running lumiCalc.'\n")    
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
for hltPath in hltPaths:
    for jsonFileName in jsonFiles.keys():
        makeFile.write("\trm -f %s\n" % MakeFile_dict[hltPath][jsonFileName]['outputFileName'])
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)

