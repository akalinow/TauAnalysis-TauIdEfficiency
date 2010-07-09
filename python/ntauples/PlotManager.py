import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import math

from ROOT import gROOT
# Prevent complaining about X server
gROOT.SetBatch(True)

import ROOT
import copy
import helpers

class LegendMaker(object):
    ''' Accumulate objects and produces TLegends of them 

    Add each time using [add_object].  Once all objects
    have been added a TLegend can be retrieved via
    the [make_legend] method.  The LegendMaker retains
    ownership/a reference of the TLegend

    '''
    def __init__(self):
        self.entries = []
        self.output = None
    def add_object(self, object, name, option="lpf"):
        # Correctly deal with filled histograms
        if option.find("hist") != -1:
            option = "lf"
        self.entries.append((object, name, option))
    def make_legend(self, x_low=0.64, y_low=0.74, x_high=0.89, y_high=0.89):
        self.output = ROOT.TLegend(x_low, y_low, x_high, y_high, "", "NDC")
        # Turn off all the stupid defaults
        self.output.SetFillColor(0)
        for entry in self.entries:
            self.output.AddEntry(*entry)
        return self.output


class PlotManager(object):
    def __init__(self):
        ROOT.TGaxis.SetMaxDigits(3)
        self.samples = {}
        # Keep track of sample ordering
        self.sample_order = []
        self.int_lumi = None

    def add_sample(self, sample_mananger, nice_name, **style_options):
        ''' Add a sample to the plot manager

        Add a sample (MC_QCD, Data, etc) to the PlotManager. The style options
        to use for plots of this sample can be sepecified using [style_options]
        keyword arguments.  The [sample_manager] input should be an instance
        of the NtupleSampleCollection class.
        
        The available style options are:
            draw_option 
            marker_color 
            marker_size 
            marker_style 
            fill_color
            fill_style
            line_width
            line_color
            x_axis_title 
            y_axis_title 

        '''
        
        if sample_mananger.name in self.samples:
            raise KeyError, "attempt to add a sample to the plot manager that already exits"
        # Create a new dictionary
        sample_dict = self.samples.setdefault(sample_mananger.name, {})
        sample_dict['sample'] = sample_mananger
        sample_dict['int_lumi'] = sample_mananger.effective_luminosity()
        sample_dict['nice_name'] = nice_name
        self.sample_order.append(sample_mananger.name)

        # Get any desired style overrides for this sample
        sample_style = copy.copy(style.DEFAULT_STYLE)
        sample_style.update(style_options)
        sample_dict['style'] = sample_style

    def update_style(self, sample, **style_options):
        " Update a samples style options "
        self.samples[sample]['style'].update(style_options)

    def set_integrated_lumi(self, target_int_lumi):
        """ Set the desired integrated luminsotity for all samples
        
        In the case target_int_lumi is None, no reweighting will be 
        applied.
        """
        self.int_lumi = target_int_lumi

    def distribution(self, expression, selection, binning=(), normalize = None, verbose=True, **options):
        " Compare a distribution from the different samples "
        if verbose: print "Plotting", str(expression), str(selection)
        # Get a unique name for the ROOT name
        unique_name = helpers.make_unique_name(expression, selection)
        result_dict = {}
        result_dict['result'] = ROOT.THStack("thStack"+unique_name, "")
        result_dict['keep'] = []
        # Keep track of all the plots we make
        result_dict['samples'] = {}
        result_dict['legend'] = LegendMaker()
        #keep track of all the plots made (yeah one could extract that from the TStack with some gymnastics...
        all_plots = {}
        for sample_name in self.sample_order:
            if verbose: print " * Sample: ", sample_name
            sample_info = self.samples[sample_name]
            # Keep track of relevant info for each sample
            plot_dict = {}
            # List to prevent garbage collection
            plot_dict['keep'] = []
            result_dict['samples'][sample_name] = plot_dict
            my_plot = plot.draw(
                # events for this sample
                list(sample_info['sample'].events_and_weights(self.int_lumi)),
                expression, selection, 
                output_name="%s%s" % (sample_name,unique_name),
                binning=binning
            )
            if my_plot == None: raise StandardError, "Error drawing '%s' for '%s' with selection '%s' no histogram was created!"%( expression, sample_name, selection)
            if verbose: print " * * - got plot", my_plot, "integral:", my_plot.Integral()
            all_plots[sample_name] = my_plot
            plot_dict['plot'] = my_plot
            # Get some handy stats
            plot_dict['rms'] = my_plot.GetRMS()
            plot_dict['mean'] = my_plot.GetMean()
            plot_dict['max'] = my_plot.GetMaximum()
            plot_dict['entries'] = my_plot.GetEntries()
            plot_dict['integral'] = my_plot.Integral()
            plot_dict['underflow'] = my_plot.GetBinContent(0)
            plot_dict['overflow'] = my_plot.GetBinContent(
                my_plot.GetNbinsX()+1)
            # Update the styles for this histogram
            plot_dict['keep'].append(
                style.update_histo_style(my_plot, sample_info['style']))
            #if we want to normalize to this distribution it will be by its integral
            if sample_name == normalize:
                normalize = plot_dict['integral']
            # Add this histogram to the stack (will be plotted using 
            # nostack option)
            result_dict['result'].Add(my_plot, sample_info['style']['draw_option'])
            # Add an entry in the legend
            result_dict['legend'].add_object(my_plot, sample_info['nice_name'], 
                                             sample_info['style']['draw_option'])
        #take care of normalization
        if not normalize == None:
            for my_plot_name, my_plot in all_plots.iteritems():
                my_plot.Scale( float(normalize) / my_plot.Integral())
            
        # Now draw the total stack to current pad
        result_dict['result'].Draw("nostack")
        # If the upper limit for y is not set, make a nice guess
        # Give us some breathing room on top of the plot to put our labels and stuff
        if 'y_max' not in options:
            if 'logy' in options and options['logy']:
                result_dict['result'].SetMaximum(result_dict['result'].GetMaximum()*20)
            else:
                result_dict['result'].SetMaximum(result_dict['result'].GetMaximum()*1.5)
        # Apply any canvas style options
        canvas_style = copy.deepcopy(style.DEFAULT_STYLE)
        canvas_style.update(options)
        result_dict['keep'].append(
            style.update_canvas_style(ROOT.gPad, canvas_style))
        # Apply any style overrides
        result_dict['keep'].append(
            style.update_histo_style(result_dict['result'], options))
        # and return all the crap we have generated in case the user is interested
        return result_dict

    def efficiency(self, expression, numerator, denominator,
                   binning=(), verbose=True, x_error_bars=True, **options):
        ''' Compare efficiencies for different samples 

        Compute the efficiency of the [numerator] selection with respect
        to the [denominator] selection, parameterized by [expression]. The efficiencies
        are computed with binomial errors using TGraphAsymmErrors.

        The horizontal error bars can be kept by setting [x_error_bars] to True.

        '''
        # Get a unique name
        unique_name = helpers.make_unique_name(expression, numerator, denominator)
        # store the background histogram
        result_dict = {}
        result_dict['background'] = None
        result_dict['keep'] = []
        result_dict['samples'] = {}
        result_dict['legend'] = LegendMaker()
        # Loop over each sample
        if verbose: print "Plotting", str(expression), str(numerator), str(denominator)
        for sample_name in self.sample_order:
            if verbose: print "  *  Sample:", sample_name
            sample_info = self.samples[sample_name]
            # Build a background histogram and efficiency curve for this sample
            # Nota Bene: TGraphAsymmErrors cannot handle non-unit weights.  For
            # non-merged samples (i.e. no pt hat bins), we can just set the target
            # luminosity to None, and each sample will get unit scaling.  Since we 
            # only care about the ratio, this is fine.  However, pt-hat binned samples
            # the errors will be overestimated in regions of phase space where high-lumi
            # samples are the only contribution.  (i.e. it will think that there are a few
            # events for QCD jets of PT = 100, when actually the QCD_80 sample contributes
            # a lot.
            plot.efficiencyLogHack(" for %s "%sample_name)
            my_bkg_hist, my_eff = plot.efficiency(
                #list(sample_info['sample'].events_and_weights(self.int_lumi)),
                list(sample_info['sample'].events_and_weights(None)),
                expression, numerator, denominator, binning=binning,
                output_name="%s%s"%(sample_name, unique_name))

            # If this is the first histogram drawn, add the background
            if result_dict['background'] is None:
                # Apply style to background histo
                background_style = copy.copy(style.DEFAULT_STYLE)
                background_style.update(options)
                style.update_histo_style(my_bkg_hist, background_style)
                # Draw and store background
                result_dict['background'] = my_bkg_hist
                my_bkg_hist.Draw()
            # Update style off efficiency curve
            result_dict['keep'].append(
                style.update_histo_style(my_eff, sample_info['style']))
            # Check if we want to remove the horizontal error bars
            if not x_error_bars:
                style.remove_x_error_bars(my_eff)
            # Store and draw efficiency
            result_dict['samples'][sample_name] = my_eff
            my_eff.Draw('e1p, same')
            # Build legend
            result_dict['legend'].add_object(
                my_eff, sample_info['nice_name'], 'p')

        # Apply any canvas style updates
        canvas_style = copy.deepcopy(style.DEFAULT_STYLE)
        canvas_style.update(options)
        result_dict['keep'].append(
            style.update_canvas_style(ROOT.gPad, canvas_style))
        return result_dict

    def multi_efficiency(self, expression, denominator, numerators, binning,
                         style_map={"mc_qcd_pythia8": style.MC_STYLES, "data":style.DATA_STYLES},
                         **options):
        if len(numerators) < 1:
            raise ValueError, "you have to pass at least one numerator (got '%s')!"%numerators
        if len(binning) < 3:
            raise ValueError, "the binning tuple has to be at leal of length 3 (numbins, first, last) or more (list of bin edges) (got '%s')!"%binning
        output = {}
        my_legend = LegendMaker()
        output['legend'] = my_legend
        output['numerators'] = {}
        output['keep'] = []
        
        cumulative_numerator = None
        for numerator_info in numerators:
            if cumulative_numerator is None:
                cumulative_numerator = numerator_info['expr']
            else:
                cumulative_numerator = cumulative_numerator & \
                        numerator_info['expr']

            numerator_name = numerator_info['style_name']
            plot.efficiencyLogHack("####### %s ######\n"%numerator_name)
            plot.efficiencyLogHack("",timestamp=True)
            #set style for the given numerator
            for sample_name in self.samples:
                style_to_use = style_map[sample_name][numerator_info["style_name"]]
                self.update_style(sample_name, **style_to_use)
            
            # Build the single efficiency
            print "Building fake rates for", numerator_name
            single_eff = self.efficiency(expression, cumulative_numerator,
                                         denominator, binning, options)

            # No GC please
            output['keep'].append(single_eff)

            # Build background
            if 'background' not in output:
                output['background'] = single_eff['background']
                background_style = copy.copy(style.DEFAULT_STYLE)
                background_style.update(options)
                style.update_histo_style(output['background'], background_style)

            # Set the correct marker for each and store it
            output['numerators'][numerator_name] = {}
            for sample_name, efficiency in single_eff['samples'].iteritems():
                print "Getting efficiency for", sample_name
                output['numerators'][numerator_name][sample_name] = efficiency
                my_legend.add_object(
                    efficiency, " ".join(
                        [numerator_info['nice_name'], 
                         self.samples[sample_name]['nice_name']]), "p")

        # Now draw everything
        print "Drawing background", output['background']
        output['background'].Draw()
        for numerator, numerator_info in output['numerators'].iteritems():
            for sample_name, efficiency in numerator_info.iteritems():
                print "Drawing", numerator, sample_name, efficiency
                efficiency.Draw("e1p,same")
                
        # Update the canvas
        canvas_style = copy.deepcopy(style.DEFAULT_STYLE)
        canvas_style.update(options)
        output['keep'].append(
            style.update_canvas_style(ROOT.gPad, canvas_style))
        return output

    def plot_dist_deviations(self, plot_result, base_sample, probe_samples = [], **options):
        output = {}
        # Construct background
        background = plot_result['samples'][base_sample]['plot'].Clone(
            plot_result['samples'][base_sample]['plot'].GetName() + "_diff")

        background.Reset()
        background_style_options = {
            'y_max' : 1.0,
            'y_min': -1.0,
        }
        ROOT.gPad.SetLogy(False)
        background_style_options.update(options)
        style.update_histo_style(background, options)
        # CV: SetRangeUser overwrites y-axis range;
        #     I don't think this is what we want (?)
        #background.GetYaxis().SetRangeUser(-2, 2)
        background.SetStats(False)
        background.Draw()
        output['background'] = background
        output['samples'] = {}
        output['legend'] = LegendMaker()

        base_histo = plot_result['samples'][base_sample]['plot']
        for probe_name in probe_samples:
            print "Comparing", probe_name
            probe_histo = plot_result['samples'][probe_name]['plot']
            difference_histo = probe_histo.Clone(probe_histo.GetName()+"_diff")
            difference_histo.Reset()
            output['samples'][probe_name] = difference_histo
            # Set each bin
            for bin in range(difference_histo.GetNbinsX()+1):
                probe_content = probe_histo.GetBinContent(bin)
                base_content = base_histo.GetBinContent(bin)
                probe_error = probe_histo.GetBinError(bin)
                base_error = base_histo.GetBinError(bin)
                print bin, base_content, probe_content
                if probe_content > 0:
                    diff = (base_content - probe_content)/probe_content

                    base_error_norm = 0
                    if base_content > 0:
                        base_error_norm = base_error/base_content

                    err = (base_content/probe_content)*math.sqrt(
                        (probe_error/probe_content)**2 + base_error_norm**2)

                    difference_histo.SetBinContent(bin, diff)
                    difference_histo.SetBinError(bin, err)
                    
            style.update_histo_style(difference_histo, self.samples[probe_name]['style'])
            output['legend'].add_object(
                difference_histo, self.samples[probe_name]['nice_name'], 'p')
            difference_histo.Draw("same,pe")

        return output






