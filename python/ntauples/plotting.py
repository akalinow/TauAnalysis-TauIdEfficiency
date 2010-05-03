import ROOT

'''
Functions to produce plots using TauNtuples

Author: Evan K. Friis, UC Davis
'''

def draw(expression, selection="", output_name="", binning=(), options="goff"):
    " Plot a histogram of expression "
    if not output_name:
        output_name = "htemp"
    events_to_draw_from = expression.events
    # Make sure selection and expression are using the same events
    # just as a sanity check
    if hasattr(selection, "events"):
        if events_to_draw_from is not selection.events:
            raise ValueError, "expression event source is not the same as selection!"
    # Build binning string if desired
    binning_str = ""
    if binning:
        binning_str = "(%s)" % ','.join(str(x) for x in binning)
    draw_command = expression.value + ">>" + output_name + binning_str
    events_to_draw_from.Draw(str(draw_command), str(selection), options)
    return ROOT.gDirectory.Get(output_name)

def efficiency(expression, numerator="", output_name="", denominator="", **kwargs):
    if not output_name:
        output_name = "eff_temp"
    numerator_h = draw(expression, numerator, "numerator_temp", **kwargs)
    denominator_h = draw(expression, denominator, "denominator_temp", **kwargs)
    #print numerator_h, denominator_h
    # Build a blank histogram w/ correct x & y axes to draw the efficiency on
    histogram_background = numerator_h.Clone("%s_bkg" % output_name)
    histogram_background.GetXaxis().SetTitle(str(expression))
    histogram_background.GetYaxis().SetTitle("Efficiency")
    histogram_background.GetYaxis().SetRangeUser(1e-4, 1)
    histogram_background.Reset()
    histogram_background.SetTitle("%s/%s" % (numerator, denominator))
    # Compute efficiency
    efficiency = ROOT.TGraphAsymmErrors(numerator_h, denominator_h)
    efficiency.SetName(output_name)
    return (histogram_background, efficiency)

