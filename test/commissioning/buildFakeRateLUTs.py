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

from TauAnalysis.TauIdEfficiency.fakeRateDef import fakeRateDef, \
        numerators_to_make

skipExisting = True

for fakeRateLabel, fakeRateConfig in fakeRateDef.items():
    for sampleLabel, sampleConfig in fakeRateConfig['samples'].items():
        for numeratorLabel, numerator in fakeRateConfig['numerators'].items():
            if numeratorLabel not in numerators_to_make:
                print "Skipping", numeratorLabel
                continue
            denominator = fakeRateConfig['denominator']

            if not os.path.exists("./config"):
                os.mkdir("./config")

            cfgFileName = './config/fakerate_%s_%s_%s.cfg' % (numeratorLabel, fakeRateLabel, sampleLabel)
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
