from ROOT import gROOT
# Prevent complaining about X server
gROOT.SetBatch(True)
import ROOT

'''

Defintion of style objects used by 
classes in the ntauples package.

Author: Evan K. Friis (UC Davis)

'''

DEFAULT_STYLE = {
    'draw_option' : "",
    'marker_color' : ROOT.EColor.kBlack,
    'marker_size' : 1.0,
    'marker_style': 22,
    'fill_color': ROOT.EColor.kWhite,
    'fill_style': 0,
    'line_width': 1,
    'line_color': ROOT.EColor.kBlack,
}

HISTO_METHOD_MAP = {
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


def update_histo_style(histo, style_dict):
    " Update a histograms style, given style dict "
    for style_item, value in style_dict.iteritems():
        HISTO_METHOD_MAP[style_item](histo)(value)

