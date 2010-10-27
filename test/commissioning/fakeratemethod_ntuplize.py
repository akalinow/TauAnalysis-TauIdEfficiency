#!/usr/bin/env python

'''

Select events and train a TMVA k-Nearest Neighboor classier using the fake rate
method.

Authors: Evan K. Friis, Christian Veelken (UC Davis)

'''

# cheating
from RecoLuminosity.LumiDB import argparse

parser = argparse.ArgumentParser(
    description = "Build an output ntuple of jetPt, eta and width"
    " given an algorithm and numerator/denominator selection."
)
parser.add_argument('--ntuple', help="Name of ntuple")
parser.add_argument('--num', help="Numerator")
parser.add_argument('--den', help="Denominator")
parser.add_argument('-skipn', action='store_true', default=False,
                    help="Don't make numerator ntuple")
parser.add_argument('-skipd', action='store_true', default=False,
                    help="Don't make denominator ntuple")
parser.add_argument('--output', help="Output file")

options=parser.parse_args()

import sys
# Otherwise ROOT will try to parse the options as well Jeeze Louise
sys.argv = []
import ROOT
ROOT.gROOT.SetBatch(True)

# Definition of input files.
import samples_cache as samples

ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

# Get our specific ntuple
ntuple = ntuple_manager.get_ntuple(options.ntuple)

# Get the HLT ntuple
hlt = ntuple_manager.get_ntuple("TriggerResults")

# Build the denominator query
denominator = hlt.expr('$hltJet15U > 0.5') & \
        ntuple.expr(options.den)

# Build the numerator query
numerator = ntuple.expr(options.num) & denominator

draw_string = ntuple.expr('$jetPt:$jetEta:$jetWidth')

output_file = ROOT.TFile(options.output, "RECREATE")
output_file.cd()

events_list = list(samples.data.events_and_weights())
if len(events_list) > 1:
    print "More than one TChain at once not supported!"
    sys.exit(1)

events = events_list[0][0]

print "Ntuple has entries: ", events.GetEntries()
# Prevent root from cutting us off at 10^6
events.SetEstimate(100000000L)

if not options.skipn:
    # Build a TPolyMaker3D with our points
    events.Draw(str(draw_string), str(numerator))
    numerator_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
    numerator_tuple.SetName("numerator")

    print "Built nuemrator pre-ntuple with %i entries" % \
            numerator_tuple.GetN()
    numerator_tuple.Write()

if not options.skipd:
    # Build a TPolyMaker3D with our points
    drawn = events.Draw(str(draw_string), str(denominator))
    denominator_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
    denominator_tuple.SetName("denominator")
    print "Built denominator pre-ntuple with %i entries" % \
            denominator_tuple.GetN()
    denominator_tuple.Write()

