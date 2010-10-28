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
parser.add_argument('--hlt', help="HLT ntuple selection")
parser.add_argument('-passing', action='store_true', default=False,
                    help="Build ntuple for events that pass")
parser.add_argument('-failing', action='store_true', default=False,
                    help="Build ntuple for events that fail")
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
denominator = hlt.expr(options.hlt) & \
        ntuple.expr(options.den)

# Build the queries for passing and failing
passing = ntuple.expr(options.num) & denominator
failing = ntuple.expr(options.num).false() & denominator

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

if options.passing:
    # Build a TPolyMaker3D with our points
    events.Draw(str(draw_string), str(passing))
    numerator_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
    numerator_tuple.SetName("pass")

    print "Built 'passing' pre-ntuple with %i entries" % \
            numerator_tuple.GetN()

    if numerator_tuple.GetN() == events.GetEstimate():
        print "Error: Number of entries produced is at the upper limit " \
                "of TTree::Draw.  You need to set events.SetEstimate to " \
                "be larger than the total number of selected rows."
        sys.exit(1)
    numerator_tuple.Write()

if options.failing:
    # Build a TPolyMaker3D with our points
    drawn = events.Draw(str(draw_string), str(failing))
    denominator_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
    denominator_tuple.SetName("failing")
    print "Built 'failining' pre-ntuple with %i entries" % \
            denominator_tuple.GetN()
    if denominator_tuple.GetN() == events.GetEstimate():
        print "Error: Number of entries produced is at the upper limit " \
                "of TTree::Draw.  You need to set events.SetEstimate to " \
                "be larger than the total number of selected rows."
        sys.exit(1)
    denominator_tuple.Write()

