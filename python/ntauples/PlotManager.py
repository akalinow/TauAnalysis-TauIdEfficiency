import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot

from ROOT import gROOT
# Prevent complaining about X server
gROOT.SetBatch(True)

import ROOT
import copy
import hashlib

_GOOD_MARKERS = []
_GOOD_COLORS = []

_HISTO_METHOD_MAP = {
    'draw_option' : lambda histo: histo.SetDrawOption ,
    'marker_color' : lambda histo: histo.SetMarkerColor,
    'marker_size' : lambda histo: histo.SetMarkerSize,
    'marker_style' : lambda histo: histo.SetMarkerStyle,
    'fill_color': lambda histo: histo.SetFillColor,
    'fill_style': lambda histo: histo.SetFillStyle,
    'line_width': lambda histo: histo.SetLineWidth,
    'line_color': lambda histo: histo.SetLineColor,
    'x_axis_title' : lambda histo: histo.GetXaxis().SetTitle,
    'y_axis_title' : lambda histo: histo.GetYaxis().SetTitle,
}

_DEFAULT_STYLE = {
    'draw_option' : "",
    'marker_color' : ROOT.EColor.kBlack,
    'marker_size' : 1.0,
    'marker_style': 22,
    'fill_color': ROOT.EColor.kWhite,
    'fill_style': 0,
    'line_width': 1,
    'line_color': ROOT.EColor.kBlack,
}

def _update_histo_style(histo, style_dict):
    " Update a histograms style, given style dict "
    for style_item, value in style_dict.iteritems():
        _HISTO_METHOD_MAP[style_item](histo)(value)

class LegendMaker(object):
    " Accumulate objects and produces TLegends of them "
    def __init__(self):
        self.entries = []
    def add_object(self, object, name, option="lpf"):
        self.entries.append((object, name, option))

class PlotManager(object):
    def __init__(self):
        self.samples = {}
        self.int_lumi = None

    def add_sample(self, sample_mananger, nice_name, **style_options):
        if sample_mananger.name in self.samples:
            raise KeyError, "attempt to add a sample to the plot manager that already exits"
        # Create a new dictionary
        sample_dict = self.samples.setdefault(sample_mananger.name, {})
        sample_dict['sample'] = sample_mananger
        sample_dict['int_lumi'] = sample_mananger.effective_luminosity()
        sample_dict['nice_name'] = nice_name

        # Get any desired style overrides for this sample
        sample_style = copy.copy(_DEFAULT_STYLE)
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

    def compare_distributions(self, expression, selection, binning=(), 
                              **options):
        " Compare a distribution from the different samples "
        # Get a unique name for the ROOT name
        unique_name = hashlib.md5(str(expression) + str(selection)).hexdigest()
        result_dict = {}
        result_dict['result'] = ROOT.THStack(unique_name, "blah")
        # Keep track of all the plots we make
        result_dict['samples'] = {}
        result_dict['legend'] = LegendMaker()
        for sample_name, sample_info in self.samples.iteritems():
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
            _update_histo_style(my_plot, sample_info['style'])
            # Add this histogram to the stack (will be plotted using 
            # nostack option)
            result_dict['result'].Add(my_plot, sample_info['style']['draw_option'])
            # Add an entry in the legend
            result_dict['legend'].add_object(my_plot, sample_info['nice_name'], 
                                             sample_info['style']['draw_option'])

        # Now draw the total stack to current pad
        result_dict['result'].Draw("nostack")
        # and return all the crap we have generated in case the user is interested
        return result_dict
