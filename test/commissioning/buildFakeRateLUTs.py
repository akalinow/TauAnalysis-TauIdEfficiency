#!/usr/bin/env python

#--------------------------------------------------------------------------------
# Build Makefile snippets needed to run
#   ls config/fakerate*cfg | xargs -P 5 -I % make -f Makefile.FakeRateMethod -j 2 CFG=%
# command for building k-NearestNeighbour trees used as look-up tables
# for tau id. efficiencies/fake-rates
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

print("<buildFakeRateLUTs>:")

import os

import TauAnalysis.TauIdEfficiency.fakeRateDef as fakerates

skipExisting = True

sample_filterfunc = None
# Uncomment to skip loose iso
sample_filterfunc = lambda x: x.find('closure') != -1

type_filter = None
type_filter = lambda x: x.find('WplusJets') != -1

for fakeRateType, fakeRateConfig in fakerates.fakeRateDef.items():
    for sampleLabel, sampleConfig in fakeRateConfig['samples'].items():
        for fake_rate_label, fake_rate_impl in fakeRateConfig['fake_rates'].items():

            if fakerates.fake_rates_to_make and (
                fake_rate_label not in fakerates.fake_rates_to_make):
                print "Skipping", fake_rate_label
                continue

            if type_filter:
                if not type_filter(fakeRateType):
                    print "Skipping", fakeRateType
                    continue

            if sample_filterfunc:
                if not sample_filterfunc(sampleLabel):
                    print "Skipping", sampleLabel
                    continue


            denominator = fake_rate_impl['denominator']
            numerator = fake_rate_impl['numerator']

            if not os.path.exists("./config"):
                os.mkdir("./config")

            cfgFileName = './config/fakerate_%s_%s_%s.cfg' % (
                fake_rate_label, fakeRateType, sampleLabel)

            if os.path.exists(cfgFileName) and skipExisting:
                print "Skipping ", cfgFileName
                continue

            cfgFile = open(cfgFileName, "w")
            cfgFile.write("[fake_rate]\n")
            cfgFile.write("sample = %s\n" % sampleConfig['name'])
            cfgFile.write("denominator = %s\n" % denominator)
            cfgFile.write("numerator = %s\n" % numerator)
            if fakeRateConfig.get('hlt') is not None:
                cfgFile.write("hlt = %s\n" % fakeRateConfig['hlt'])
            cfgFile.close()

print("Configuration files build:")
print("execute 'ls config/fakerate*cfg | xargs -P 5 -I % make -f Makefile.FakeRateMethod -j 2 CFG=%' in order to start make process.")
