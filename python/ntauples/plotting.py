import ROOT
import array

'''
Functions to produce plots using TauNtuples

Author: Evan K. Friis, UC Davis
'''

def draw(events, expression, selection="", output_name="", binning=(), options="goff"):
    ''' Plot a histogram of [expression] given [events].

    If [events] is given as list of tuples of type (TTree, weight), the output
    histogram will be the weighted combination of the events in each tree.
    
    Can optionally pass an event selection.  Output histogram is stored in 
    output_name.  Specific binning is optional, and of the form (nBins, x_low, x_high).
    The options are passed to TTree::Draw.
    '''
    # Put into common format
    if not isinstance(events, list):
        events = [ (events, 1.0) ]
    if not output_name:
        output_name = "htemp"

    options = "e, %s" % options

    # Sum the histogram from each group of events
    for event_source_index, (event_source, weight) in enumerate(events):
        # Build the draw string
        expression_str = str(expression.value)
        selection_with_weight = str(selection*weight)

        # Check if this is the first event
        if event_source_index == 0:
            # Determine how to build the histogram
            if not binning:
                # If no binning information is provided
                event_source.Draw(expression_str + ">>" + output_name, 
                                  selection_with_weight, options)
            elif len(binning) == 3: 
                # Fixed bins
                event_source.Draw(expression_str + ">>" + output_name + "(" + 
                                  ','.join(str(x) for x in binning) + ")",
                                  selection_with_weight, options)
            else: 
                # Variable bin sizes, must pre-build histo
                histogram = ROOT.TH1F(output_name, output_name, "", 
                                      len(binning)-1, array.array('d', binning))
                # Append to histogram
                event_source.Draw(expression_str + ">>+" + output_name, 
                                  selection_with_weight, options)
        else:
            # Otherwise we are appending to an existing histogram
            event_source.Draw(expression_str+">>+">>output_name, 
                              selection_with_weight, options)
            
    return ROOT.gDirectory.Get(output_name)

def efficiency(events, expression, numerator="", denominator="", output_name="", **kwargs):
    ''' Compute the efficiency versus expression

    Returns a tuple containing a background TH1F (used to draw the axis)
    and a TGraphAsymmErrors that gives the efficiency of the numerator selection
    w.r.t the denominator selection, parameterized by the given expression.  Numerator
    and denomiantor are computed using draw(...).  The kwargs are passed to draw(...).
    '''
    if not output_name:
        output_name = "eff_temp"
    numerator_h = draw(events, expression, numerator, "numerator_temp", **kwargs)
    denominator_h = draw(events, expression, denominator, "denominator_temp", **kwargs)
    #print numerator_h, denominator_h
    # Build a blank histogram w/ correct x & y axes to draw the efficiency on
    histogram_background = numerator_h.Clone("%s_bkg" % output_name)
    histogram_background.GetXaxis().SetTitle(str(expression))
    histogram_background.GetYaxis().SetTitle("Efficiency")
    histogram_background.GetYaxis().SetRangeUser(1e-4, 1)
    histogram_background.Reset()
    histogram_background.SetTitle("%s/%s" % (numerator, denominator))
    histogram_background.SetStats(False)
    # Compute efficiency
    efficiency = ROOT.TGraphAsymmErrors(numerator_h, denominator_h)
    efficiency.SetName(output_name)
    return (histogram_background, efficiency)

