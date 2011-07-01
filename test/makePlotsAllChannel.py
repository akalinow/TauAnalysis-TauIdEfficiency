#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms
import os
import sys
import TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi as samples
import TauAnalysis.TauIdEfficiency.templates as templates
import glob

samplesToAnalyze = [
    'data_SingleMu_Run2011A_PromptReco_v1',
    'data_SingleMu_Run2011A_PromptReco_v2',
#    'DYtautauM10to20_powheg',
    'Ztautau_powheg',
#    'qqZll',
#    'DYmumuM10to20_pythia',
    'Zmumu_powheg',
    'PPmuXptGt20Mu15',
    'WplusJets_madgraph',
    'TTplusJets_madgraph'
]

outputFilePath = "/nfs/data4/verzetti/scratch/CMSSW_4_1_3/src/TauAnalysis/TauIdEfficiency/test/"
inputFilePath = "/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/"
run = False

if inputFilePath[-1] != '/':
    inputFilePath += '/'

if not os.path.isdir(inputFilePath):
    print inputFilePath + ' is not a directory, exiting...'
    sys.exit(0)

if not os.path.isdir(outputFilePath):
    print outputFilePath + ' is not a directory, creating it...'
    os.mkdir(outputFilePath)

if outputFilePath[-1] != '/':
    outputFilePath += '/'

if len(samplesToAnalyze) == 0:
    samplesToAnalyze = samples.recoSampleDefinitionsZtoMuTau_7TeV['SAMPLES_TO_ANALYZE']

for currentSample in samplesToAnalyze:
    infiles = glob.glob(inputFilePath+"tauIdEffMeasPATTuple*"+currentSample+'*.root')
    if len(infiles) == 0:
        print 'WARNING! No files found for sample: ' +currentSample
        continue
    infilesPrint = ""
    for infile in infiles[:-1]:
        infilesPrint += "'file:"+infile + "', "
    infilesPrint += "'file:"+infiles[-1] + "'"
    outfileName = outputFilePath + "regionPlots_" + currentSample + ".root"
    channelName = ''
    if currentSample.find('data') != -1:
        channelName = 'data'
    else:
        channelName = currentSample
    ## if not os.path.isdir('channelsConfig'):
    ##     print 'channelsConfig directory not found, making it...'
    ##     os.mkdir('channelsConfig')
    cfgfileName = 'plotChannel_'+currentSample+'_cfi.py'
    cfgfile = file(cfgfileName,'w')
    cfgfile.write(templates.plotMacroTemplate % (infilesPrint,outfileName,channelName))
    command = "./make_plots_channel_cfg.py 'load="+cfgfileName.replace('.py','')+"' &"#.replace('/','.').replace('.py','')+"' &"
    if run:
        os.system(command)
    else:
        print command

