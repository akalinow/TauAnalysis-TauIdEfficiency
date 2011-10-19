#!/usr/bin/env python

import os

hltPaths = [
    #'HLT_L1SingleMuOpen_v*',
    #'HLT_L1SingleMuOpen_DT_v*',
    #'HLT_L1SingleMu10_v*',
    #'HLT_IsoMu15_v*',
    #'HLT_IsoMu12_v*',
    #'HLT_Mu15_v*',    # for QCD muon enriched 
    #'HLT_IsoMu17_v*', # for W --> mu nu, Z --> mu+ mu-
    #'HLT_Jet30_v*',   # for QCD multi-jet
    #'HLT_Jet60_v*',
    #'HLT_Jet80_v*',
    #'HLT_Jet110_v*',
    #'HLT_Jet150_v*'
    'HLT_Mu17_Ele8_CaloIdL_v*',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',
    'HLT_Mu8_Ele17_CaloIdL_v*',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*'
]

jsonFile = \
  '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-177053_7TeV_PromptReco_Collisions11_JSON.txt'
#'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt'
#'/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-172802_7TeV_PromptReco_Collisions11_JSON_v3.txt'

executable_lumiCalc = '/afs/cern.ch/user/v/veelken/scratch0/CMSSW_4_2_4_patch1/src/RecoLuminosity/LumiDB/scripts/lumiCalc2.py'
executable_analyzeLumiCalcOutput = './analyzeLumiCalcOutput.py'

MakeFile_dict = {}
for hltPath in hltPaths:
    MakeFile_dict[hltPath] = {}
    hltPath_abbr = hltPath
    hltPath_abbr = hltPath_abbr.replace("_v*", "")
    hltPath_abbr = hltPath_abbr.replace("*", "")    
    MakeFile_dict[hltPath]['outputFileName'] = 'lumiCalc_%s.out' % hltPath_abbr
    MakeFile_dict[hltPath]['lumiCalc_command'] = \
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
    outputFileNames.append(MakeFile_dict[hltPath]['outputFileName'])
makeFile.write("all: parseLumiCalcOutput\n")
makeFile.write("\n")
for hltPath in hltPaths:
    makeFile.write("%s:\n" % MakeFile_dict[hltPath]['outputFileName'])
    makeFile.write("\t%s &> %s\n" % (MakeFile_dict[hltPath]['lumiCalc_command'], MakeFile_dict[hltPath]['outputFileName']))
makeFile.write("\n")    
makeFile.write("parseLumiCalcOutput: %s\n" % make_MakeFile_vstring(outputFileNames))
for hltPath in hltPaths:
    makeFile.write("\t%s %s\n" % (executable_analyzeLumiCalcOutput, MakeFile_dict[hltPath]['outputFileName']))
makeFile.write("\techo 'Finished running lumiCalc.'\n")    
makeFile.write("\n")
makeFile.write(".PHONY: clean\n")
makeFile.write("clean:\n")
for hltPath in hltPaths:
    makeFile.write("\trm -f %s\n" % MakeFile_dict[hltPath]['outputFileName'])
makeFile.write("\techo 'Finished deleting old files.'\n")
makeFile.write("\n")
makeFile.close()

print("Finished building Makefile. Now execute 'make -j 8 -f %s'." % makeFileName)

