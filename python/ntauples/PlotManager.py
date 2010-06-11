import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

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
        self.entries.append((object, name, option))
    def make_legend(self, x_low=0.70, y_low=0.8, x_high=0.95, y_high=0.95):
        self.output = ROOT.TLegend(x_low, y_low, x_high, y_high, "", "NDC")
        # Turn off all the stupid defaults
        self.output.SetFillColor(0)
        for entry in self.entries:
            self.output.AddEntry(*entry)
        return self.output


class PlotManager(object):
    def __init__(self):
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

    def distribution(self, expression, selection, binning=(), 
                     title="CMS Preliminary", verbose=True, **options):
        " Compare a distribution from the different samples "
        if verbose:
            print "Plotting", str(expression), str(selection)
        # Get a unique name for the ROOT name
        unique_name = helpers.make_unique_name(expression, selection)
        result_dict = {}
        result_dict['result'] = ROOT.THStack(unique_name, title)
        # Keep track of all the plots we make
        result_dict['samples'] = {}
        result_dict['legend'] = LegendMaker()
        for sample_name in self.sample_order:
            if verbose: print " * Sample: ", sample_name
            sample_info = self.samples[sample_name]
            # Keep track of relevant info for each sample
            plot_dict = {}
            result_dict['samples'][sample_name] = plot_dict
            my_plot = plot.draw(
                # events for this sample
                list(sample_info['sample'].events_and_weights(self.int_lumi)),
                expression, selection, 
                output_name="%s%s" % (sample_name,unique_name),
                binning=binning
            )
            if verbose: print " * * - got plot", my_plot, "integral:", my_plot.Integral()
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
            style.update_histo_style(my_plot, sample_info['style'])
            # Add this histogram to the stack (will be plotted using 
            # nostack option)
            result_dict['result'].Add(my_plot, sample_info['style']['draw_option'])
            # Add an entry in the legend
            result_dict['legend'].add_object(my_plot, sample_info['nice_name'], 
                                             sample_info['style']['draw_option'])

        # Now draw the total stack to current pad
        result_dict['result'].Draw("nostack")
        # Apply any canvas style options
        canvas_style = copy.deepcopy(style.DEFAULT_STYLE)
        canvas_style.update(options)
        style.update_canvas_style(ROOT.gPad, canvas_style)
        # Apply any style overrides
        style.update_histo_style(result_dict['result'], options)
        # and return all the crap we have generated in case the user is interested
        return result_dict

    def efficiency(self, expression, numerator, denominator,
                   binning=(), x_error_bars=False, **options):
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
        result_dict['samples'] = {}
        result_dict['legend'] = LegendMaker()
        # Loop over each sample
        for sample_name in self.sample_order:
            sample_info = self.samples[sample_name]
            # Build a background histogram and efficiency curve for this sample
            my_bkg_hist, my_eff = plot.efficiency(
                list(sample_info['sample'].events_and_weights(self.int_lumi)),
                expression, numerator, denominator, binning=binning,
                output_name="%s%s"%(sample_name, unique_name))

            if result_dict['background'] is None:
                # Apply style to background histo
                background_style = copy.copy(style.DEFAULT_STYLE)
                background_style.update(options)
                style.update_histo_style(my_bkg_hist, background_style)
                # Draw and store background
                result_dict['background'] = my_bkg_hist
                my_bkg_hist.Draw()
            # Update style off efficiency curve
            style.update_histo_style(my_eff, sample_info['style'])
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
        style.update_canvas_style(ROOT.gPad, canvas_style)
        return result_dict
