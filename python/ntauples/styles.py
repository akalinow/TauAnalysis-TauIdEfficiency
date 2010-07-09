from copy import deepcopy
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
CMS_PRELIMINARY_UPPER_LEFT = ROOT.TPaveText(0.125, 0.865, 0.45, 0.905, "NDC")
CMS_PRELIMINARY_UPPER_LEFT.AddText("CMS Preliminary")
CMS_PRELIMINARY_UPPER_LEFT.SetTextAlign(13)
CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.040)
CMS_PRELIMINARY_UPPER_LEFT.SetFillStyle(0)
CMS_PRELIMINARY_UPPER_LEFT.SetBorderSize(0)

# Define the luminosity labels
LUMI_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.810, 0.45, 0.850, "NDC")
LUMI_LABEL_UPPER_LEFT.AddText("L = 8.4nb^{-1}")
LUMI_LABEL_UPPER_LEFT.SetTextAlign(13)
LUMI_LABEL_UPPER_LEFT.SetTextSize(0.035)
LUMI_LABEL_UPPER_LEFT.SetFillStyle(0)
LUMI_LABEL_UPPER_LEFT.SetBorderSize(0)

#define Z--> tau tau label
ZTAUTAU_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.810, 0.45, 0.850, "NDC")
ZTAUTAU_LABEL_UPPER_LEFT.AddText("Z #rightarrow #tau #tau Simulation")
ZTAUTAU_LABEL_UPPER_LEFT.SetTextAlign(13)
ZTAUTAU_LABEL_UPPER_LEFT.SetTextSize(0.035)
ZTAUTAU_LABEL_UPPER_LEFT.SetFillStyle(0)
ZTAUTAU_LABEL_UPPER_LEFT.SetBorderSize(0)

#define center-of-mass energy label
SQRTS_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.765, 0.45, 0.805, "NDC")
SQRTS_LABEL_UPPER_LEFT.AddText("#sqrt{s} = 7 TeV")
SQRTS_LABEL_UPPER_LEFT.SetTextAlign(13)
SQRTS_LABEL_UPPER_LEFT.SetTextSize(0.035)
SQRTS_LABEL_UPPER_LEFT.SetFillStyle(0)
SQRTS_LABEL_UPPER_LEFT.SetBorderSize(0)

# Define the kinematic cut labels
PT_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.715, 0.45, 0.755, "NDC")
PT_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 10 GeV/c")
PT_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.035)
PT_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

PTVIS_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.715, 0.45, 0.755, "NDC")
PTVIS_CUT_LABEL_UPPER_LEFT.AddText("P_{T}^{vis} > 10 GeV/c")
PTVIS_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
PTVIS_CUT_LABEL_UPPER_LEFT.SetTextSize(0.035)
PTVIS_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
PTVIS_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

ETA_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.715, 0.45, 0.755, "NDC")
ETA_CUT_LABEL_UPPER_LEFT.AddText("|#eta| < 2.5")
ETA_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.035)
ETA_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
ETA_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

ETAVIS_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.715, 0.45, 0.755, "NDC")
ETAVIS_CUT_LABEL_UPPER_LEFT.AddText("|#eta^{vis}| < 2.5")
ETAVIS_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
ETAVIS_CUT_LABEL_UPPER_LEFT.SetTextSize(0.035)
ETAVIS_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
ETAVIS_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

PT_ETA_CUT_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.715, 0.45, 0.755, "NDC")
PT_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 10 GeV/c, |#eta| < 2.5")
PT_ETA_CUT_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.035)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_ETA_CUT_LABEL_UPPER_LEFT.SetBorderSize(0)

PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.665, 0.45, 0.755, "NDC")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("P_{T} > 10 GeV/c,\n")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta| < 2.5")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.035)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetBorderSize(0)

PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.650, 0.45, 0.750, "NDC")
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("P_{T}^{vis}(gen) > 10 GeV/c,\n")
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta|(gen) < 2.5")
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.040)
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetBorderSize(0)

# Labels for different pt cuts
PT_20_ETA_CUT_LABEL_UPPER_LEFT = PT_ETA_CUT_LABEL_UPPER_LEFT.Clone()
PT_20_ETA_CUT_LABEL_UPPER_LEFT.Clear()
PT_20_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c, |#eta| < 2.5")

PT_40_ETA_CUT_LABEL_UPPER_LEFT = PT_ETA_CUT_LABEL_UPPER_LEFT.Clone()
PT_40_ETA_CUT_LABEL_UPPER_LEFT.Clear()
PT_40_ETA_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 40 GeV/c, |#eta| < 2.5")

EFFICIENCY_STYLES = {
    'matching' : {
        'marker_color' : ROOT.EColor.kBlack,
        'marker_style' : 29, # closed star
        'marker_size': 2,
    },
    'byLeadTrackFinding' : {
        'marker_color' : ROOT.EColor.kRed,
        'marker_style' : 21, # closed square
    },
    'byLeadTrackPtCut' : {
        'marker_color' : ROOT.EColor.kBlue,
        'marker_style' : 22, # closed upward-facing triangle
    },
    'byTrackIsolation' : {
        'marker_color' : ROOT.EColor.kGreen + 2,
        'marker_style' : 24, # open circle
    },
    'byEcalIsolation' : {
        'marker_color' : ROOT.EColor.kOrange + 7,
        'marker_style' : 25 # open square
    },
    'OneOrThreeProng' : {
        'marker_color' : ROOT.EColor.kMagenta + 2,
        'marker_style' : 26, # open upward-facing triangle
    },
    # For HPS
    'byDecayModeFinding' : {
        'marker_color' : 28,
        'marker_style' : 23, # closed downward-facing triangle
    },
}

openToCloseMap = { 24:20, 25:21, 26:22 }
MC_STYLES = {"OneOrThreeProng": {
    'marker_color' : ROOT.EColor.kBlack,
    'marker_style' : 24, # open dot
    }
             }
DATA_STYLES = deepcopy( MC_STYLES)
for dataStyle in DATA_STYLES:
    DATA_STYLES[dataStyle]["marker_style"] = openToCloseMap[DATA_STYLES[dataStyle]["marker_style"]]

# TaNC, HPS and calo tau colors follow the same pattern
for iso_type, tanc_type in zip(
    [ 'byTrackIsolation', 'byEcalIsolation', 'OneOrThreeProng' ],
    [ 'byTaNCfrOnePercent', 'byTaNCfrHalfPercent', 'byTaNCfrQuarterPercent' ]):
    EFFICIENCY_STYLES[tanc_type] = EFFICIENCY_STYLES[iso_type]
    MC_STYLES[tanc_type] = EFFICIENCY_STYLES[iso_type]
    DATA_STYLES[tanc_type] = deepcopy(EFFICIENCY_STYLES[iso_type])
    DATA_STYLES[tanc_type]["marker_style"] = openToCloseMap[EFFICIENCY_STYLES[iso_type]["marker_style"]]

for iso_type, hps_type in zip(
    [ 'byTrackIsolation', 'byEcalIsolation', 'OneOrThreeProng' ],
    [ 'byIsolationLoose', 'byIsolationMedium', 'byIsolationTight' ]):
    EFFICIENCY_STYLES[hps_type] = EFFICIENCY_STYLES[iso_type]
    MC_STYLES[hps_type] = EFFICIENCY_STYLES[iso_type]
    DATA_STYLES[hps_type] = deepcopy(EFFICIENCY_STYLES[iso_type])
    DATA_STYLES[hps_type]["marker_style"] = openToCloseMap[EFFICIENCY_STYLES[iso_type]["marker_style"]]

for iso_type, calo_type in zip(
    [ 'byTrackIsolation', 'byEcalIsolation', 'OneOrThreeProng' ],
    [ 'byIsolation', 'byEcalIsolation_calo', 'OneOrThreeProng_calo' ]):
    EFFICIENCY_STYLES[calo_type] = EFFICIENCY_STYLES[iso_type]
    MC_STYLES[calo_type] = EFFICIENCY_STYLES[iso_type]
    DATA_STYLES[calo_type] = deepcopy(EFFICIENCY_STYLES[iso_type])
    DATA_STYLES[calo_type]["marker_style"] = openToCloseMap[EFFICIENCY_STYLES[iso_type]["marker_style"]]
    
# make sure that all errorbars have the right colors
for effStyle in EFFICIENCY_STYLES:
    EFFICIENCY_STYLES[effStyle]["line_color"] = EFFICIENCY_STYLES[effStyle]["marker_color"]
    EFFICIENCY_STYLES[effStyle]["fill_color"] = EFFICIENCY_STYLES[effStyle]["marker_color"]
for mcStyle in MC_STYLES:
    MC_STYLES[mcStyle]["line_color"] = MC_STYLES[mcStyle]["marker_color"]
    MC_STYLES[mcStyle]["fill_color"] = MC_STYLES[mcStyle]["marker_color"]
for dataStyle in DATA_STYLES:
    DATA_STYLES[dataStyle]["line_color"] = DATA_STYLES[dataStyle]["marker_color"]
    DATA_STYLES[dataStyle]["fill_color"] = DATA_STYLES[dataStyle]["marker_color"]

def make_nice_numbers(x, sig_figs = 2):
    format_string = "%0." + str(sig_figs) + "e"
    result = format_string%x
    if result.endswith("e+00"):
        result = result[:-4]
    else:
        result = result.replace("e", " 10^{")
        result = result.replace("{-0", "{-")
        result = result.replace("{+0", "{+")
        result+="}"
        
    return result

# Build a pave with the mean and RMS
def make_mean_rms_pave(plot, sig_figs = 2, x_low=0.55, y_low=0.77, x_high=0.89, y_high=0.89):
    mean = plot.GetMean()
    rms = plot.GetRMS()
    output = ROOT.TPaveText(x_low, y_low, x_high, y_high, "brNDC")
    output.SetTextAlign(33) # top right alignment
    output.SetBorderSize(0)
    output.SetTextSize(0.045)
    output.SetFillStyle(0)
    
    output.AddText("mean = "+ make_nice_numbers(mean) )
    output.AddText("rms = "+ make_nice_numbers(rms))
    return output

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
    'labels' : [CMS_PRELIMINARY_UPPER_LEFT, LUMI_LABEL_UPPER_LEFT, SQRTS_LABEL_UPPER_LEFT],
    # Easy way to take on more labels
    'extra_labels' : []
}

DATA_STYLE = {
    'draw_option' : "e1p",      # draw error bars
    'marker_style' : 20,        # solid circle
    'marker_size' : 1.0,
    #'line_width' : 1.5,
    'marker_color' : ROOT.EColor.kBlack,
    'line_color' :  ROOT.EColor.kBlack,
    'fill_style' : 0,
}

QCD_MC_PYTHIA6_STYLE_HIST = {
    'draw_option' : "hist",
    'line_color' :  ROOT.EColor.kGreen-2,
    #'line_width' : 1.5,
    'fill_color' :  ROOT.EColor.kGreen-4,
    'fill_style' : 3005,
    'marker_style' : 25,
    'marker_color' : ROOT.EColor.kGreen,
}

QCD_MC_PYTHIA6_STYLE_DOTS = {
    'draw_option' : "e1p",
    'marker_style' : 25,        # open square
    'marker_size' : 1.0,
    'marker_color' : ROOT.EColor.kGreen,
    'line_color' :  ROOT.EColor.kBlack,
    #'line_width' : 1.5,
}

QCD_MC_PYTHIA8_STYLE_HIST = {
    'draw_option' : "hist",
    'line_color' :  ROOT.EColor.kAzure-2,
    #'line_width' : 1.5,
    'fill_color' :  ROOT.EColor.kAzure-4,
    'fill_style' : 3004,
    'marker_style' : 24,
    'marker_color' : ROOT.EColor.kAzure,
}

QCD_MC_PYTHIA8_STYLE_DOTS = {
    'draw_option' : "e1p",
    'marker_style' : 24,        # open circle
    'marker_size' : 1.0,
    'marker_color' : ROOT.EColor.kAzure,
    'line_color' :  ROOT.EColor.kBlack,
    #'line_width' : 1.5,
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
    'marker_style' : 24,        # open circle
    'marker_size' : 1.0,
    'marker_color' : ROOT.EColor.kAzure,
    'line_color' :  ROOT.EColor.kBlack,
    #'line_width' : 1.5,
}

MINBIAS_MC_STYLE = {
    'draw_option' : "e1p",
    'marker_style' : 21,        # solid square 
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
    'y_axis_label_size' : lambda histo: histo.GetYaxis().SetLabelSize,
    'x_axis_label_size' : lambda histo: histo.GetXaxis().SetLabelSize,
    'n_divisions' : lambda histo: histo.GetXaxis().SetNdivisions,
    'y_min' : lambda histo: histo.SetMinimum,
    'y_max' : lambda histo: histo.SetMaximum,
    'title' : lambda histo: histo.SetTitle,
}

def adjust_coordinates(pad, x1, y1, x2, y2):

    result_dict = {}

    # Adjust position and size of labels/legends to margins of pad
    pad_topMargin = pad.GetTopMargin()
    pad_rightMargin = pad.GetRightMargin()
    pad_bottomMargin = pad.GetBottomMargin()
    pad_leftMargin = pad.GetLeftMargin()

    pad_widthX = 1.0 - (pad_leftMargin + pad_rightMargin)
    pad_widthY = 1.0 - (pad_topMargin + pad_bottomMargin)

    pad_defaultTopMargin = 0.1
    pad_defaultRightMargin = 0.1
    pad_defaultBottomMargin = 0.1
    pad_defaultLeftMargin = 0.1
    
    pad_defaultWidthX = 1.0 - (pad_defaultLeftMargin + pad_defaultRightMargin)
    pad_defaultWidthY = 1.0 - (pad_defaultTopMargin + pad_defaultBottomMargin)
    
    result_dict['x1'] = pad_leftMargin + (x1 - pad_defaultLeftMargin)*(pad_widthX)/pad_defaultWidthX
    result_dict['y1'] = pad_bottomMargin + (y1 - pad_defaultBottomMargin)*(pad_widthY)/pad_defaultWidthY
    result_dict['x2'] = pad_leftMargin + (x2 - pad_defaultLeftMargin)*(pad_widthX)/pad_defaultWidthX
    result_dict['y2'] = pad_bottomMargin + (y2 - pad_defaultBottomMargin)*(pad_widthY)/pad_defaultWidthY

    return result_dict

def draw_labels(pad):
    ''' Returns a function that takes a list of labels 
    and draws them on the current pad
    '''
    def apply_labels(labels):
        pad.cd()
        to_keep = []
        for label in labels:
            #label.Draw()

            # Adjust position and size of label to margins of pad
            newLabel = label.Clone()
            adjusted_coordinates = adjust_coordinates(pad, label.GetX1(), label.GetY1(), label.GetX2(), label.GetY2())
            newLabel.SetX1(adjusted_coordinates['x1'])
            newLabel.SetY1(adjusted_coordinates['y1'])
            newLabel.SetX2(adjusted_coordinates['x2'])
            newLabel.SetY2(adjusted_coordinates['y2'])
            newLabel.Draw()
            
            # Prevent garbage collection
            #to_keep.append(label)
            to_keep.append(newLabel)
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




