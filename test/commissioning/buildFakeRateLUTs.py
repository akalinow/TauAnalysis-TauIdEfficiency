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

import subprocess

from TauAnalysis.TauIdEfficiency.fakeRateDef import fakeRateDef

for fakeRateLabel, fakeRateConfig in fakeRateDef.items():
    for sampleLabel, sampleConfig in fakeRateConfig['samples'].items():
        for numeratorLabel, numerator in fakeRateConfig['numerators'].items():
            denominator = fakeRateConfig['denominator']
            
            cfgFileName = './config/fakerate_%s_%s_%s.cfg' % (numeratorLabel, fakeRateLabel, sampleLabel)
            cfgFile = open(cfgFileName, "w")
            cfgFile.write("[fake_rate]")
            cfgFile.write("sample = %s" % sampleConfig['name'])
            cfgFile.write("denominator = %s" % denominator)
            cfgFile.write("numerator = %s" % numerator)
            if fakeRateConfig.get('hlt') is not None:
                cfgFile.write("hlt = %s" % fakeRateConfig['hlt'])
            cfgFile.close()

print("Configuration files build:")
print("execute 'ls config/fakerate*cfg | xargs -P 5 -I % make -f Makefile.FakeRateMethod -j 2 CFG=%' in order to start make process.")

        
