#!/usr/bin/env python

'''

Select events and train a TMVA k-Nearest Neighboor classier using the fake rate
method.

Authors: Evan K. Friis, Christian Veelken (UC Davis)

'''

# cheating
from RecoLuminosity.LumiDB import argparse
import ConfigParser

parser = argparse.ArgumentParser(
    description = "Build an output ntuple of jetPt, eta and width"
    " given an algorithm and numerator/denominator selection."
)
parser.add_argument('--ntuple', help="Name of ntuple")
parser.add_argument('-passing', action='store_true', default=False,
                    help="Build ntuple for events that pass")
parser.add_argument('-failing', action='store_true', default=False,
                    help="Build ntuple for events that fail")
parser.add_argument('--output', help="Output file")
parser.add_argument('--config', help="Configuration file")
parser.add_argument('--index', help="Sample index")
parser.add_argument(
    '--kin', help='[tau] or [jet], default = tau', default='jet')

options=parser.parse_args()

import sys
# Otherwise ROOT will try to parse the options as well Jeeze Louise
sys.argv = []
import ROOT
ROOT.gROOT.SetBatch(True)

# Load configuration file
config = ConfigParser.ConfigParser()
config.read(options.config)

options.sample = config.get('fake_rate', 'sample')
options.num = config.get('fake_rate', 'numerator')
options.den = config.get('fake_rate', 'denominator')
# Allow missing HLT value - no cut will be applied
options.hlt = (config.has_option('fake_rate', 'hlt') and
               config.get('fake_rate', 'hlt') or None)

print "Building fake rate mini-ntuple with options:"
print " sample =", options.sample
print " numerator =", options.num
print " denominator =", options.den
print " hlt bit =", (options.hlt is not None and options.hlt or "NO SELECTION")

# Definition of input files.
import samples_cache as samples

ntuple_manager = getattr(
    samples, options.sample).build_ntuple_manager("tauIdEffNtuple")

# Get our specific ntuple
try:
    ntuple = ntuple_manager.get_ntuple(options.ntuple)
except:
    print ntuple_manager
    raise

# Build the denominator query
denominator = None
if options.hlt:
    # Get the HLT ntuple
    hlt = ntuple_manager.get_ntuple("patTriggerEvent")
    denominator = hlt.expr(options.hlt) & ntuple.expr(options.den)
else:
    denominator = ntuple.expr(options.den)

# Build the queries for passing and failing
passing = ntuple.expr(options.num) & denominator
failing = ntuple.expr(options.num).false() & denominator

# Make sure no spurious quotes have entered the strings
def check_for_quotes(my_string):
    my_string = str(my_string)
    " Make sure there aren't any quotes in a string "
    if my_string.find('"') != -1 or my_string.find("'") != -1:
        print "Either the numerator, denominator or HLT bit contains" \
                " extra quotation marks, please remove them."
        print " Offending string:", my_string
        sys.exit(1)

check_for_quotes(passing)
check_for_quotes(failing)


draw_string = "bad"
if options.kin.find('tau') != -1:
    draw_string = ntuple.expr('$jetWidth:abs($eta):$pt')
elif options.kin.find('jet') != -1:
    draw_string = ntuple.expr('$jetWidth:$jetEta:$jetPt')
else:
    raise ValueError("bad kinematic var option - use jet or tau!")
print "Parameterization string:", draw_string

output_file = ROOT.TFile(options.output, "RECREATE")
output_file.cd()

events_list = list(getattr(samples, options.sample).events_and_weights())

# WARNING: current method of filling k-NN tree using TPolyMarker3D objects
#          does not support event weights
#         --> need to make sure that all events have the same weights,
#             by requiring sample to either contain one set of files
#             or all sets to have the same weights in case sample contains multiple set of files
#         (Note that this means it is not possible to fill k-NN tree
#          with QCD Monte Carlo generated in different PtHat bins)
#import TauAnalysis.TauIdEfficiency.ntauples.helpers as helpers

#events = None
#if len(events_list) == 1:
    #events = events_list[0][0]
#elif len(events_list) > 1:
    #events = ROOT.TChain("Events")
    #isFirst = True
    #weight = 0.
    #for entry in events_list:
        #if isFirst:
            #helpers.copy_aliases_from(entry[0], events)
            #weight = entry[1]
            #isFirst = False
        #else:
            #if weight != entry[1]:
                #raise ValueError("Events with different weights not supported yet !!")
        #print("--> adding %s events to TChain..." % entry[0].GetEntries())
        #print entry[0]
        #events.Add(entry[0])

#if events is None:
    #raise ValueError("Failed to access TChain !!")

weights = ROOT.TArrayD(len(events_list))

# Normalize the sample weights
weight_norm = sum(weight for events, weight in events_list)

to_write = []

max_entries = sum(events.GetEntries() for events, weight in events_list)

# Increase max events if necessary
map(lambda x: max(10000, x.SetEstimate(max_entries*3)),
    (x for x, weight in events_list))

ROOT.SetMemoryPolicy( ROOT.kMemoryStrict )

# Write a TPolyMarker containing indices/evts/runs of passed & failed jets in
# events.
write_indices = False

for index, (events, weight) in enumerate(events_list):
    #if not index == 0:
    #    helpers.copy_aliases_from(events_list[0][0], events)
    print "Ntuple has entries: ", events.GetEntries(), "and weight:", weight
    # Prevent root from cutting us off at 10^6
    #if index != 1:
    #    continue
    weights[index] = weight/weight_norm
    if options.passing:
        # Build a TPolyMaker3D with our points
        selected_events = events.Draw(str(draw_string), str(passing))
        #selected_events = events.Draw(str(draw_string))
        print "TTree::Drew ntuple with %i events" % selected_events
        my_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
        #my_tuple = numerator_tuple.Clone("passing_%i" % index)
        my_tuple.SetName("passing_%i" % index)

        print "Built 'passing' pre-ntuple with %i entries" % \
               my_tuple.GetN()

        if my_tuple.GetN() == events.GetEstimate():
           raise ValueError("Number of entries produced is at the upper limit " \
                            "of TTree::Draw.  You need to set events.SetEstimate to " \
                            "be larger than the total number of selected rows.")
        bytes = my_tuple.Write()
        print "Wrote %i bytes" % bytes
        print "Deleting..."
        my_tuple.IsA().Destructor(my_tuple)
        print "Done."

        # If desired, write indexes to log
        if write_indices:
            index_draw_str = str(ntuple.expr("$index:event:run"))
            selected_events = events.Draw(index_draw_str, str(passing))
            print " Drew index ntuple with %i events" % selected_events
            my_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
            #my_tuple = numerator_tuple.Clone("passing_%i" % index)
            my_tuple.SetName("passing_indices_%i" % index)
            bytes = my_tuple.Write()
            print "Wrote %i bytes" % bytes
            print "Deleting..."
            my_tuple.IsA().Destructor(my_tuple)
            print "Done."

    if options.failing:
        # Build a TPolyMaker3D with our points
        drawn = events.Draw(str(draw_string), str(failing))
        print "TTree::Drew ntuple with %i events" % drawn
        my_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
        #my_tuple = denominator_tuple.Clone("failing_%i" % index)
        my_tuple.SetName("failing_%i" % index)
        print "Built 'failing' pre-ntuple with %i entries" % \
                my_tuple.GetN()

        if my_tuple.GetN() == events.GetEstimate():
            raise ValueError("Number of entries produced is at the upper limit " \
                             "of TTree::Draw.  You need to set events.SetEstimate to " \
                             "be larger than the total number of selected rows.")
        my_tuple.Write()
        ROOT.gPad.Clear()

        if write_indices:
            index_draw_str = str(ntuple.expr("$index:event:run"))
            selected_events = events.Draw(index_draw_str, str(failing))
            print " Drew index ntuple with %i events" % selected_events
            my_tuple = ROOT.gPad.GetPrimitive("TPolyMarker3D")
            #my_tuple = numerator_tuple.Clone("passing_%i" % index)
            my_tuple.SetName("failing_indices_%i" % index)
            bytes = my_tuple.Write()
            print "Wrote %i bytes" % bytes
            print "Deleting..."
            my_tuple.IsA().Destructor(my_tuple)
            print "Done."

output_file.WriteObject(weights, "weights")
