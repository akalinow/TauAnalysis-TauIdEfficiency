import ROOT
import array

'''
Functions to produce plots using TauNtuples

Author: Evan K. Friis, UC Davis
'''

def draw(events, expression, selection = "", output_name = "", binning = (), options = "goff", maxNumEntries = 1000000000):
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

        # Check if this is the first event source
        if event_source_index == 0:
            # Determine how to build the histogram
            #
            # CV: always pre-build histogram before calling TTree::Draw,
            #     to make sure that TH1::Sumw2 gets called and binErrors are computed correctly.
            #
            if not binning:
                raise ValueError("Undefined binning Parameter --> Please update your code !!")
            histogram = None
            if len(binning) == 3:
                # Fixed bin width
                histogram = ROOT.TH1F(output_name, output_name,
                                      binning[0], binning[1], binning[2])
            else:
                # Variable bin widhts
                histogram = ROOT.TH1F(output_name, output_name, 
                                      len(binning) - 1, array.array('d', binning))
            histogram.Sumw2()
            # Fill histogram
            event_source.Draw(expression_str + ">>+" + output_name, 
                              selection_with_weight, options, maxNumEntries)
        else:
            # Otherwise we are appending to an existing histogram
            event_source.Draw(expression_str + ">>+" + output_name, 
                              selection_with_weight, options, maxNumEntries)
            
    return ROOT.gDirectory.Get(output_name)

def getEfficiencyGraph(numerator, denominator):
    ''' Compute efficiency graph from numerator and denominator histograms.

        Make numerator and denominator histograms integer
        before computing the tau id. efficiency/fake-rate,
        in order to work-around the problem
        that the TGraphAsymmErrors(TH1*, TH1*) constructor does not work in ROOT 5.27/06b (default in CMSSW_4_1_3),
        in case the histograms contain weighted entries.

        The binContents of numerator and denominator histograms are scaled
        to the "effective" number of (denominator) entries:
          effNum = (sum(weights)/sqrt(sum(weights^2)))^2 = binContent/(binError^2)
    '''

    numerator_cloned = numerator.Clone()
    denominator_cloned = denominator.Clone()
    
    for iBin in range(denominator_cloned.GetNbinsX() + 1):
        
        scaleFactor = 0.
        if denominator.GetBinError(iBin) > 0:
            scaleFactor = denominator.GetBinContent(iBin)/ROOT.TMath.Power(denominator.GetBinError(iBin), 2.)
            
        numerator_cloned.SetBinContent(iBin, ROOT.TMath.Nint(scaleFactor*numerator.GetBinContent(iBin)))
        numerator_cloned.SetBinError(iBin, ROOT.TMath.Sqrt(numerator_cloned.GetBinContent(iBin)));
        
        denominator_cloned.SetBinContent(iBin, ROOT.TMath.Nint(scaleFactor*denominator.GetBinContent(iBin)));
        denominator_cloned.SetBinError(iBin, ROOT.TMath.Sqrt(denominator_cloned.GetBinContent(iBin)));

    return ROOT.TGraphAsymmErrors(numerator_cloned, denominator_cloned)

def efficiency(events, expression, numerator = "", denominator = "", binning = (),
               output_name = "", maxNumEntries = 1000000000, **kwargs):
    ''' Compute the efficiency versus expression

    Returns a tuple containing a background TH1F (used to draw the axis)
    and a TGraphAsymmErrors that gives the efficiency of the numerator selection
    w.r.t the denominator selection, parameterized by the given expression.  Numerator
    and denomiantor are computed using draw(...).  The kwargs are passed to draw(...).
    '''

    #print("<efficiency>:")
    #print " kwargs = ", kwargs 
    
    if not output_name:
        output_name = "eff_temp"        
    numerator_h   = draw(events, expression, numerator, "numerator_temp",
                         binning = binning, maxNumEntries = maxNumEntries, **kwargs)
    print("--> numerator = %.0f" % numerator_h.Integral())
    denominator_h = draw(events, expression, denominator, "denominator_temp",
                         binning = binning, maxNumEntries = maxNumEntries, **kwargs)
    print("--> denominator = %.0f" % denominator_h.Integral())
    #FIXME clean this up
    from math import sqrt
    nNum = float(numerator_h.Integral())
    nDenom = float(denominator_h.Integral())
    if nDenom == 0:
        print("Error in <efficiency>: nDenom = 0 --> skipping !!")
        return 
    elif nNum > nDenom:
        print("Error in <efficiency>: nNum = %e exceeds nDenom = %e --> skipping !!" % (nNum, nDenom))
        return 
    err = 1/nDenom*sqrt(nNum*(1-nNum/nDenom) )
    efficiencyLogHack("%e / %e = %e +- %e\n" % (nNum, nDenom, nNum/nDenom, err))
    efficiencyLogHack("", timestamp=True)
    #end cleanup
    
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
    print("CV: enabling temporary work-around for problem with TGraphAsymmErrors(TH1*, TH1*)")
    print("    in case of histograms containing weighted entries")
    print "numerator_h:"
    print "numBins = ", numerator_h.GetNbinsX()
    for i in range(numerator_h.GetNbinsX() + 1):
        print " bin i = ", i, ":", numerator_h.GetBinContent(i), " +/- ", numerator_h.GetBinError(i)
    print "denominator_h:"
    print "numBins = ", denominator_h.GetNbinsX()
    for i in range(denominator_h.GetNbinsX() + 1):
        print " bin i = ", i, ":", denominator_h.GetBinContent(i), " +/- ", denominator_h.GetBinError(i)
    # CV: temporary workaround
    #    - TGraphAsymmErrors(TH1*, TH1*) constructor does not work in ROOT 5.27/06b (default in CMSSW_4_1_3),
    #      in case histograms contain weighted entries,
    #      see posting in RootTalk mailing list http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=12534
    efficiency = getEfficiencyGraph(numerator_h, denominator_h)
    efficiency.SetName(output_name)
    return (histogram_background, efficiency)

def efficiencyLogHack(output, timestamp = False):
    from time import strftime
    import sys
    ts = ""
    if timestamp:
        ts = strftime("%D - %T")
    f = open("plots/efficiency_%s.log"%("_".join(sys.argv[1:])),"a")
    f.write("%s: %s"%(ts, output))
    f.close()
