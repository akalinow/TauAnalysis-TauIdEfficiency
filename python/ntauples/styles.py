from ROOT import gROOT
# Prevent complaining about X server
gROOT.SetBatch(True)
import ROOT

# Setup common styles
ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetTitleBorderSize(0)

'''

Defintion of style objects used by 
classes in the ntauples package.

Author: Evan K. Friis (UC Davis)

'''

# Define the CMS preliminary labels
CMS_PRELIMINARY_UPPER_LEFT = ROOT.TPaveText(0.12, 0.87, 0.45, 0.92, "NDC")
CMS_PRELIMINARY_UPPER_LEFT.AddText("CMS Preliminary 7 TeV")
CMS_PRELIMINARY_UPPER_LEFT.SetTextAlign(13)
CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.04)
CMS_PRELIMINARY_UPPER_LEFT.SetFillStyle(0)
CMS_PRELIMINARY_UPPER_LEFT.SetBorderSize(0)

# Define the luminosity labels
LUMI_LABEL_UPPER_LEFT = ROOT.TPaveText(0.12, 0.80, 0.45, 0.85, "NDC")
LUMI_LABEL_UPPER_LEFT.AddText("#int L = 8.19nb^{-1}")
LUMI_LABEL_UPPER_LEFT.SetTextAlign(13)
LUMI_LABEL_UPPER_LEFT.SetTextSize(0.04)
LUMI_LABEL_UPPER_LEFT.SetFillStyle(0)
LUMI_LABEL_UPPER_LEFT.SetBorderSize(0)

# Define the kinematic cut labels
PT_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.12, 0.73, 0.45, 0.78, "NDC")
PT_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 10 GeV/c")
PT_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.04)
PT_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

ETA_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.12, 0.73, 0.45, 0.78, "NDC")
ETA_CUT_LABEL_UPPER_LEFT.AddText("|#eta| < 2.5")
ETA_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.04)
ETA_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
ETA_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

PT_ETA_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.12, 0.73, 0.45, 0.78, "NDC")
PT_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 10 GeV/c, |#eta| < 2.5")
PT_ETA_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.04)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

# Labels for different pt cuts
PT_20_ETA_CUT_LABEL_UPPER_LEFT = PT_ETA_CUT_LABEL_UPPER_LEFT.Clone()
PT_20_ETA_CUT_LABEL_UPPER_LEFT.Clear()
PT_20_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c, |#eta| < 2.5")

PT_40_ETA_CUT_LABEL_UPPER_LEFT = PT_ETA_CUT_LABEL_UPPER_LEFT.Clone()
PT_40_ETA_CUT_LABEL_UPPER_LEFT.Clear()
PT_40_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 40 GeV/c, |#eta| < 2.5")

DEFAULT_STYLE = {
    'draw_option' : "",
    'marker_color' : ROOT.EColor.kBlack,
    'marker_size' : 1.0,
    'marker_style': 20,
    'fill_color': ROOT.EColor.kWhite,
    'fill_style': 0,
    #'line_width': 1.0,
    'line_color': ROOT.EColor.kBlack,
    'title' : "",
    'y_axis_title' : "Fake Rate",
    'horizontal_grid' : 0,
    'vertical_grid' : 1,
    'labels' : [CMS_PRELIMINARY_UPPER_LEFT, LUMI_LABEL_UPPER_LEFT],
    # Easy way to take on more labels
    'extra_labels' : []
}

DATA_STYLE = {
    'draw_option' : "e1p",      # draw error bars
    'marker_style' : 20,        # dot
    'marker_size' : 1.0,
    #'line_width' : 1.5,
    'marker_color' : ROOT.EColor.kBlack,
    'line_color' :  ROOT.EColor.kBlack,
    'fill_style' : 0,
}

QCD_MC_STYLE_HIST = {
    'draw_option' : "hist",
    'line_color' :  ROOT.EColor.kAzure-2,
    #'line_width' : 1.5,
    'fill_color' :  ROOT.EColor.kAzure-4,
    'fill_style' : 3003,
    'marker_style' : 24,
    'marker_color' : ROOT.EColor.kAzure,
}

QCD_MC_STYLE_DOTS = {
    'draw_option' : "e1p",
    'marker_style' : 24,        # open dot
    'marker_size' : 1.0,
    'marker_color' : ROOT.EColor.kAzure,
    'line_color' :  ROOT.EColor.kBlack,
    #'line_width' : 1.5,
}

MINBIAS_MC_STYLE = {
    'draw_option' : "e1p",
    'marker_style' : 21,        # square dot
    'marker_size' : 1.0,
    'marker_color' : ROOT.EColor.kRed,
    'line_color' :  ROOT.EColor.kBlack,
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
    'y_min' : lambda histo: histo.SetMinimum,
    'y_max' : lambda histo: histo.SetMaximum,
    'title' : lambda histo: histo.SetTitle,
}

def draw_labels(pad):
    ''' Returns a function that takes a list of labels 
    and draws them on the current pad
    '''
    def apply_labels(labels):
        pad.cd()
        to_keep = []
        for label in labels:
            label.Draw()
            # Prevent garbage collection
            to_keep.append(label)
        return to_keep
    return apply_labels

CANVAS_METHOD_MAP = {
    'horizontal_grid' : lambda pad: pad.SetGridx,
    'vertical_grid' : lambda pad: pad.SetGridy,
    'logy' : lambda pad: pad.SetLogy,
    'labels' : draw_labels,
    'extra_labels' : draw_labels
}

def update_canvas_style(pad, style_dict):
    " Update any canvas level options "
    to_keep = [] # prevent GC
    for style_item, value in style_dict.iteritems():
        if style_item in CANVAS_METHOD_MAP:
            # Keep the item if necessary
            item = CANVAS_METHOD_MAP[style_item](pad)(value)
            if item is not None:
                to_keep.append(item)
    return to_keep 

def update_histo_style(histo, style_dict):
    " Update a histograms style, given style dict "
    to_keep = []
    for style_item, value in style_dict.iteritems():
        if style_item in HISTO_METHOD_MAP:
            item = HISTO_METHOD_MAP[style_item](histo)(value)
            if item is not None:
                to_keep.append(item)
        elif style_item in CANVAS_METHOD_MAP:
            pass # ignore
        else:
            print "Warning: Unrecognized style option %s" % style_item
    return to_keep

def remove_x_error_bars(tgraph):
    " Remove horizontal error bars from graph "
    for point in range(tgraph.GetN()):
        tgraph.SetPointEXhigh(point, 0)
        tgraph.SetPointEXlow(point, 0)




