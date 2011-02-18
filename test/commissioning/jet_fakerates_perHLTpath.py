#!/usr/bin/env python

'''
Plot fake-rates of different tau id algorithms for
QCD multi-jet events triggered by different HLT paths

Author: Christian Veelken

'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as colors
import copy

# Define of input files.
import samples_cache as samples

# Define sample names
samples_dijet_data = {}
samples_dijet_data['HLT_Jet15U']  = samples.data_dijet_wPU
samples_dijet_data['HLT_Jet30U']  = samples.data_dijet_wPU
samples_dijet_data['HLT_Jet50U']  = samples.data_dijet_wPU
samples_dijet_data['HLT_Jet70U']  = samples.data_dijet_wPU
samples_dijet_data['HLT_Jet100U'] = samples.data_dijet_wPU

sample_dijet_mc                 = samples.qcddijetPU156bx_mc

# Define integrated luminosity of analyzed dataset
intLumiData = 36100 # nb^-1

# Define denominator (phase-space) for which fake-rates are to be determined
denominator_phase_space = "$jetPt > 20.0 & abs($jetEta) < 2.3"

# Define maximum number of Ntuple entries to process
# (NOTE: use values different from 1000000000 for debugging purposes only !!)
maxNumEntries = 1000000000
#maxNumEntries = 10000

# Define marker styles
custom_styles_dijet_data = {}
custom_styles_dijet_data['HLT_Jet15U'] = copy.deepcopy(style.DATA_STYLE)
custom_styles_dijet_data['HLT_Jet15U']['marker_color']  = colors.color_black.value()
custom_styles_dijet_data['HLT_Jet15U']['line_color']    = colors.color_black.value()
custom_styles_dijet_data['HLT_Jet15U']['marker_style']  = 20 # closed circle
custom_styles_dijet_data['HLT_Jet30U'] = copy.deepcopy(style.DATA_STYLE)
custom_styles_dijet_data['HLT_Jet30U']['marker_color']  = colors.color_violett.value()
custom_styles_dijet_data['HLT_Jet30U']['line_color']    = colors.color_violett.value()
custom_styles_dijet_data['HLT_Jet30U']['marker_style']  = 21 # closed square
custom_styles_dijet_data['HLT_Jet50U'] = copy.deepcopy(style.DATA_STYLE)
custom_styles_dijet_data['HLT_Jet50U']['marker_color']  = colors.color_green.value()
custom_styles_dijet_data['HLT_Jet50U']['line_color']    = colors.color_green.value()
custom_styles_dijet_data['HLT_Jet50U']['marker_style']  = 22 # closed upward-facing triangle
custom_styles_dijet_data['HLT_Jet70U'] = copy.deepcopy(style.DATA_STYLE)
custom_styles_dijet_data['HLT_Jet70U']['marker_color']  = colors.color_red.value()
custom_styles_dijet_data['HLT_Jet70U']['line_color']    = colors.color_red.value()
custom_styles_dijet_data['HLT_Jet70U']['marker_style']  = 23 # closed downward-facing triangle
custom_styles_dijet_data['HLT_Jet100U'] = copy.deepcopy(style.DATA_STYLE)
custom_styles_dijet_data['HLT_Jet100U']['marker_color'] = colors.color_darkBlue.value()
custom_styles_dijet_data['HLT_Jet100U']['line_color']   = colors.color_darkBlue.value()
custom_styles_dijet_data['HLT_Jet100U']['marker_style'] = 29 # closed star

dijet_dataHLTpaths = [ 'HLT_Jet15U', 'HLT_Jet30U', 'HLT_Jet50U', 'HLT_Jet70U', 'HLT_Jet100U' ]

custom_style_dijet_mc = copy.deepcopy(style.QCD_MC_STYLE_DOTS)
custom_style_dijet_mc['marker_color']                  = colors.color_black.value()
custom_style_dijet_mc['line_color']                    = colors.color_black.value()
custom_style_dijet_mc['marker_style']                  = 24 # open circle

# Adjust position of labels
custom_CMS_PRELIMINARY_UPPER_LEFT = copy.deepcopy(style.CMS_PRELIMINARY_UPPER_LEFT)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetX1(style.CMS_PRELIMINARY_UPPER_LEFT.GetX1() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetX2(style.CMS_PRELIMINARY_UPPER_LEFT.GetX2() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetY1(style.CMS_PRELIMINARY_UPPER_LEFT.GetY1() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetY2(style.CMS_PRELIMINARY_UPPER_LEFT.GetY2() + 0.050)

custom_LUMI_LABEL_UPPER_LEFT = copy.deepcopy(style.LUMI_LABEL_UPPER_LEFT)
custom_LUMI_LABEL_UPPER_LEFT.SetX1(style.LUMI_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetX2(style.LUMI_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetY1(style.LUMI_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetY2(style.LUMI_LABEL_UPPER_LEFT.GetY2() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.Clear()
custom_LUMI_LABEL_UPPER_LEFT.AddText("Data, L = 36.1pb^{-1}")
custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.0375)

custom_SQRTS_LABEL_UPPER_LEFT = copy.deepcopy(style.SQRTS_LABEL_UPPER_LEFT)
custom_SQRTS_LABEL_UPPER_LEFT.SetX1(style.SQRTS_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_SQRTS_LABEL_UPPER_LEFT.SetX2(style.SQRTS_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_SQRTS_LABEL_UPPER_LEFT.SetY1(style.SQRTS_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_SQRTS_LABEL_UPPER_LEFT.SetY2(style.SQRTS_LABEL_UPPER_LEFT.GetY2() + 0.050)

custom_PT_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_CUT_LABEL_UPPER_LEFT)
custom_PT_CUT_LABEL_UPPER_LEFT.SetX1(style.PT_CUT_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_PT_CUT_LABEL_UPPER_LEFT.SetX2(style.PT_CUT_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY1(style.PT_CUT_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY2(style.PT_CUT_LABEL_UPPER_LEFT.GetY2() + 0.050)
custom_PT_CUT_LABEL_UPPER_LEFT.Clear()
custom_PT_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c")

custom_ETA_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.ETA_CUT_LABEL_UPPER_LEFT)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetX1(style.ETA_CUT_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetX2(style.ETA_CUT_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY1(style.ETA_CUT_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY2(style.ETA_CUT_LABEL_UPPER_LEFT.GetY2() + 0.050)
custom_ETA_CUT_LABEL_UPPER_LEFT.Clear()
custom_ETA_CUT_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")

custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX1(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX2(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY1(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY2(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetY2() + 0.050)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.Clear()
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c,\n")
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")

# Define labels for different algorithms
algo_label_template = ROOT.TPaveText(0.50, 0.615, 0.88, 0.655, "NDC")
algo_label_template.SetTextAlign(33)
algo_label_template.SetTextSize(0.035)
algo_label_template.SetTextColor(2)
algo_label_template.SetFillStyle(0)
algo_label_template.SetBorderSize(0)

algo_label_fixed = copy.deepcopy(algo_label_template)
algo_label_fixed.AddText("Fixed Cone")

algo_label_shrinking = copy.deepcopy(algo_label_template)
algo_label_shrinking.AddText("Shrinking Cone")

algo_label_tanc_loose = copy.deepcopy(algo_label_template)
algo_label_tanc_loose.AddText("TaNC loose")
algo_label_tanc_medium = copy.deepcopy(algo_label_template)
algo_label_tanc_medium.AddText("TaNC medium")
algo_label_tanc_tight = copy.deepcopy(algo_label_template)
algo_label_tanc_tight.AddText("TaNC tight")

algo_label_hps_loose = copy.deepcopy(algo_label_template)
algo_label_hps_loose.AddText("HPS loose")
algo_label_hps_medium = copy.deepcopy(algo_label_template)
algo_label_hps_medium.AddText("HPS medium")
algo_label_hps_tight = copy.deepcopy(algo_label_template)
algo_label_hps_tight.AddText("HPS tight")

algo_label_calo = copy.deepcopy(algo_label_template)
algo_label_calo.AddText("TCTau")

def makeFakeratePlots(algorithm, numerator,
                      denominators_dijet_data, plotters_dijet_data, dijet_dataHLTpaths,
                      denominator_dijet_mc, plotter_dijet_mc,
                      expression, binning,
                      x_axis_title, y_min = 1e-3, y_max = 0.9, extra_labels = [], extra_labels_diff = [],
                      maxNumEntries = 1000000000,
                      outputFileName = "jet_fakerates_perHLTpath.pdf"):

    #print("<makeFakeratePlots>:")

    result_algorithm = {}
    
    #----------------------------------------------------------------------------
    # Make "regular" fake-rate plots
    #----------------------------------------------------------------------------

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    canvas.cd()

    for dijet_dataHLTpath in dijet_dataHLTpaths:
        eff_dijet_data = plotters_dijet_data[dijet_dataHLTpath].efficiency(
            expression,
            denominators_dijet_data[dijet_dataHLTpath] & numerator,
            denominators_dijet_data[dijet_dataHLTpath],
            binning = binning,
            maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
            x_axis_title = x_axis_title,
            y_axis_title = "Fake Rate",
            y_min = y_min, y_max = y_max, logy = True)
        print "eff_dijet_data", eff_dijet_data
        result_algorithm[dijet_dataHLTpath] = eff_dijet_data

    eff_dijet_mc = plotter_dijet_mc.efficiency(
        expression,
        denominator_dijet_mc & numerator,
        denominator_dijet_mc,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Fake Rate",
        y_min = y_min, y_max = y_max, logy = True)
    print "eff_dijet_mc", eff_dijet_mc
    result_algorithm['dijet_mc'] = eff_dijet_mc

    # Add Monte Carlo expectation to all data plots
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        result_algorithm[dijet_dataHLTpath]['samples'][sample_dijet_mc.name] = \
          result_algorithm['dijet_mc']['samples'][sample_dijet_mc.name]

    # Draw fake-rate plots
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        result_algorithm[dijet_dataHLTpath]['samples'][samples_dijet_data[dijet_dataHLTpath].name].Draw('epsame')
    result_algorithm['dijet_mc']['samples'][sample_dijet_mc.name].Draw('epsame')

    # Draw legend
    legend = ROOT.TLegend(0.45, 0.68, 0.88, 0.88, "", "brNDC")
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(1)
    legend.SetBorderSize(1)
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        legend.AddEntry(result_algorithm[dijet_dataHLTpath]['samples'][samples_dijet_data[dijet_dataHLTpath].name],
                        "Data: %s" % dijet_dataHLTpath, "p")
    legend.AddEntry(result_algorithm['dijet_mc']['samples'][sample_dijet_mc.name],
                    "MC", "p")
    legend.Draw()

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
    style.LUMI_LABEL_UPPER_LEFT.Draw()
    style.SQRTS_LABEL_UPPER_LEFT.Draw()

    for extra_label in extra_labels:
        extra_label.Draw()

    canvas.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.png" % algorithm)
    canvas.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.pdf" % algorithm)

    #----------------------------------------------------------------------------
    # Make fake-rate plots with normalized differences added below plots
    #----------------------------------------------------------------------------
    
    canvas_diff = ROOT.TCanvas("canvas_diff", "canvas_diff", 500, 625)
    canvas_diff.cd()
    
    topPad_diff = ROOT.TPad("topPad_diff", "topPad_diff", 0.00, 0.35, 1.00, 1.00)
    topPad_diff.SetFillColor(10)
    topPad_diff.SetLogy(True)
    topPad_diff.SetTopMargin(0.05)
    topPad_diff.SetLeftMargin(0.18)
    topPad_diff.SetBottomMargin(0.00)
    topPad_diff.SetRightMargin(0.05)
    topPad_diff.SetGridy()
    topPad_diff.Draw()
    topPad_diff.cd()

    # Adjust offset of y-axis title
    eff_dijet_mc['background'].GetYaxis().SetNdivisions(505)
    eff_dijet_mc['background'].GetYaxis().SetTickLength(0.04)
    eff_dijet_mc['background'].GetYaxis().SetTitleOffset(1.2)
    eff_dijet_mc['background'].GetYaxis().SetTitleSize(0.05) 
    eff_dijet_mc['background'].GetYaxis().SetLabelSize(0.05)
    eff_dijet_mc['background'].Draw()

    # Draw fake-rate plots
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        result_algorithm[dijet_dataHLTpath]['samples'][samples_dijet_data[dijet_dataHLTpath].name].Draw('epsame')
    result_algorithm['dijet_mc']['samples'][sample_dijet_mc.name].Draw('epsame')

    # Draw legend
    legend.SetX1NDC(legend.GetX1NDC() + 0.050)
    legend.SetX2NDC(legend.GetX2NDC() + 0.050)
    legend.SetY1NDC(legend.GetY1NDC() + 0.050)
    legend.SetY2NDC(legend.GetY2NDC() + 0.050)
    legend.Draw()

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    custom_CMS_PRELIMINARY_UPPER_LEFT.Draw()
    custom_LUMI_LABEL_UPPER_LEFT.Draw()
    custom_SQRTS_LABEL_UPPER_LEFT.Draw()
    
    for extra_label_diff in extra_labels_diff:
        extra_label_diff.Draw()

    canvas_diff.cd()
    
    bottomPad_diff = ROOT.TPad("bottomPad_diff", "bottomPad_diff", 0.00, 0.00, 1.00, 0.35)
    bottomPad_diff.SetFillColor(10)
    bottomPad_diff.SetLogy(False)
    bottomPad_diff.SetTopMargin(0.00)
    bottomPad_diff.SetLeftMargin(0.18)
    bottomPad_diff.SetBottomMargin(0.20)
    bottomPad_diff.SetRightMargin(0.05)
    bottomPad_diff.SetGridy()
    bottomPad_diff.Draw()
    bottomPad_diff.cd()

    # Show normalized difference (Data - MC)/MC
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        diff_dijet = plotters_dijet_data[dijet_dataHLTpath].plot_eff_deviations(
            result_algorithm[dijet_dataHLTpath],
            sample_dijet_mc.name,
            [ samples_dijet_data[dijet_dataHLTpath].name ],
            x_axis_title = x_axis_title,
            y_axis_title = "#frac{Data - Simulation}{Simulation}",
            y_min = -0.34, y_max = +0.34, logy = False
        )
        result_algorithm['%s_diff' % dijet_dataHLTpath] = diff_dijet

    diff_dijet = result_algorithm['%s_diff' % 'HLT_Jet15U']
    diff_dijet['background'].GetYaxis().SetNdivisions(505)
    diff_dijet['background'].GetXaxis().SetTitleOffset(1.1)
    diff_dijet['background'].GetXaxis().SetTitleSize(0.08)
    diff_dijet['background'].GetXaxis().SetLabelSize(0.08)
    diff_dijet['background'].GetXaxis().SetTickLength(0.055)
    diff_dijet['background'].GetYaxis().CenterTitle()
    diff_dijet['background'].GetYaxis().SetTitleOffset(0.9)
    diff_dijet['background'].GetYaxis().SetTitleSize(0.08)
    diff_dijet['background'].GetYaxis().SetLabelSize(0.08)
    diff_dijet['background'].GetYaxis().SetTickLength(0.04)
    diff_dijet['background'].Draw()

    for dijet_dataHLTpath in dijet_dataHLTpaths:
        result_algorithm['%s_diff' % dijet_dataHLTpath]['samples'][samples_dijet_data[dijet_dataHLTpath].name].Draw("epsame")

    canvas_diff.cd()
    canvas_diff.Update()

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.png" % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.pdf" % algorithm)

    return result_algorithm
       
if __name__ == "__main__":

    # define PlotManagers for QCD Monte Carlo and data samples
    plotter_dijet_mc = PlotManager()
    plotter_dijet_mc.add_sample(sample_dijet_mc, "QCDj Simulation",
                                **custom_style_dijet_mc)
    plotters_dijet_data = {}
    for dijet_dataHLTpath in dijet_dataHLTpaths:
        plotter_dijet_data = PlotManager()
        plotter_dijet_data.add_sample(samples_dijet_data[dijet_dataHLTpath], "QCDj Data: %s" % dijet_dataHLTpath,
                                      **custom_styles_dijet_data[dijet_dataHLTpath])
        plotters_dijet_data[dijet_dataHLTpath] = plotter_dijet_data
    
    # Build the ntuple manager
    ntuple_manager = sample_dijet_mc.build_ntuple_manager("tauIdEffNtuple")

    nTuples = {
        "shrinkingCone" : ntuple_manager.get_ntuple("shrinking"),
        "fixedCone"     : ntuple_manager.get_ntuple("fixed"),
        ##"TaNC"          : ntuple_manager.get_ntuple("shrinking"), # "old" TaNC
        "TaNC"          : ntuple_manager.get_ntuple("hpstanc"),     # "new" TaNC
        "hps"           : ntuple_manager.get_ntuple("hps"),
        "calo"          : ntuple_manager.get_ntuple("calo")
    }
    
    hlt = ntuple_manager.get_ntuple("patTriggerEvent")

    # Define binning options
    binning_pt  = (0, 10,  15,  20,  25, 30, 40, 60, 100, 200)
    binning_eta = (15, -2.5,  +2.5)
    binning_phi = (15, -3.14, +3.14)

    # Define the numerators to plot
    pfString    = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTrackIsolation > 0.5 " \
                 + " & $byEcalIsolation > 0.5 & ($numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3)"
    caloString  = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byIsolation > 0.5 " \
                 + " & $etSumIsolationECAL < 5 & ($numSignalTracks ==1 || $numSignalTracks ==3)"
    
    numerators = {
        "shrinkingCone" : {
            'expr'       : nTuples["shrinkingCone"].expr(pfString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        },
        "fixedCone" : {
            'expr'       : nTuples["fixedCone"].expr(pfString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        },
        ##"TaNC_loose" : { # "old" TaNC
        ##    'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrOnePercent > 0.5 & $pt > 15.0'),
        ##    'style_name' : "byTaNCfrOnePercent",
        ##    'nice_name'  : "TaNC loose"
        ##},
        ##"TaNC_medium" : {
        ##    'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrHalfPercent > 0.5 & $pt > 15.0'),
        ##    'style_name' :"byTaNCfrHalfPercent",
        ##    'nice_name'  : "TaNC medium"
        ##},
        ##"TaNC_tight" : {
        ##    'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrQuarterPercent > 0.5 & $pt > 15.0'),
        ##    'style_name' : "byTaNCfrQuarterPercent",
        ##    'nice_name'  : "TaNC tight"
        ##},
        "TaNC_loose" : { # "new" TaNC
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCloose > 0.5 & $pt > 15.0'),
            'style_name' : "byTaNCloose",
            'nice_name'  : "TaNC loose"
        },
        "TaNC_medium" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCmedium > 0.5 & $pt > 15.0'),
            'style_name' :"byTaNCmedium",
            'nice_name'  : "TaNC medium"
        },
        "TaNC_tight" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCtight > 0.5 & $pt > 15.0'),
            'style_name' : "byTaNCtight",
            'nice_name'  : "TaNC tight"
        },
        "hps_loose" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationLoose > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationLoose",
            'nice_name'  : "HPS loose"
        },
        "hps_medium" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationMedium > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationMedium",
            'nice_name'  : "HPS medium"
        },
        "hps_tight" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationTight > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationTight",
            'nice_name'  : "HPS tight"
        },
        "calo" : {
            'expr'       : nTuples["calo"].expr(caloString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        }
    }

    extra_labels = {}
    extra_labels['fixedCone']     = [ algo_label_fixed       ]
    extra_labels['shrinkingCone'] = [ algo_label_shrinking   ]
    extra_labels['TaNC_loose']    = [ algo_label_tanc_loose  ]
    extra_labels['TaNC_medium']   = [ algo_label_tanc_medium ]
    extra_labels['TaNC_tight']    = [ algo_label_tanc_tight  ]
    extra_labels['hps_loose']     = [ algo_label_hps_loose   ]
    extra_labels['hps_medium']    = [ algo_label_hps_medium  ]
    extra_labels['hps_tight']     = [ algo_label_hps_tight   ]
    extra_labels['calo']          = [ algo_label_calo        ]
        
    #for algorithm_discriminator in [ 'TaNC_loose', 'TaNC_medium', 'TaNC_tight',
    #                                 'hps_loose',  'hps_medium',  'hps_tight'  ]:
    for algorithm_discriminator in [ 'TaNC_loose' ]:

        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue

        algorithm     = algorithm_discriminator
        discriminator = ""
        if algorithm_discriminator.find("_") != -1:
            algorithm     = algorithm_discriminator[:algorithm_discriminator.find("_")]
            discriminator = algorithm_discriminator[algorithm_discriminator.find("_") + 1:]
        print("algorithm = %s, discriminator = %s" % (algorithm, discriminator))
        
        # Define the denominators
        denominator_jetId     = nTuples[algorithm].expr("$jetIdLoose > 0.5")
        denominator_e_mu_veto = nTuples[algorithm].expr("$pfElectronMVA < 0.6 & $againstMuon > 0.5")
        denominator_dijet     = nTuples[algorithm].expr(denominator_phase_space) & denominator_jetId & denominator_e_mu_veto
        denominator_dijet_mc  = denominator_dijet & hlt.expr('$hltJet15Ubit > 0.5') & nTuples[algorithm].expr("$probeJet15U > 0.5")
        denominators_dijet_data = {}
        denominators_dijet_data['HLT_Jet15U'] = \
          denominator_dijet & hlt.expr('$hltJet15Ubit > 0.5') & nTuples[algorithm].expr("$probeJet15U > 0.5")
        denominators_dijet_data['HLT_Jet30U'] = \
          denominator_dijet & hlt.expr('$hltJet30Ubit > 0.5') & nTuples[algorithm].expr("$probeJet30U > 0.5")
        denominators_dijet_data['HLT_Jet50U'] = \
          denominator_dijet & hlt.expr('$hltJet50Ubit > 0.5') & nTuples[algorithm].expr("$probeJet50U > 0.5")
        denominators_dijet_data['HLT_Jet70U'] = \
          denominator_dijet & hlt.expr('$hltJet70Ubit > 0.5') & nTuples[algorithm].expr("$probeJet70U > 0.5")
        denominators_dijet_data['HLT_Jet100U'] = \
          denominator_dijet & hlt.expr('$hltJet100Ubit > 0.5') & nTuples[algorithm].expr("$probeJet100U > 0.5")

        extra_labels_pt  = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_pt.append(style.ETA_CUT_LABEL_UPPER_LEFT)
        extra_labels_diff_pt  = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_pt.append(custom_ETA_CUT_LABEL_UPPER_LEFT)

        makeFakeratePlots(algorithm_discriminator, numerators[algorithm_discriminator]['expr'],
                          denominators_dijet_data, plotters_dijet_data, dijet_dataHLTpaths,
                          denominator_dijet_mc, plotter_dijet_mc,
                          nTuples[algorithm].expr('$jetPt'), binning_pt,
                          'Jet P_{T} [GeV/c]', 6.5e-4, 5.9, extra_labels_pt, extra_labels_diff_pt,
                          maxNumEntries,
                          "plots/jet_fakerates_perHLTpath_jetPt.pdf")

        extra_labels_eta = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_eta.append(style.PT_CUT_LABEL_UPPER_LEFT)
        extra_labels_diff_eta = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_eta.append(custom_PT_CUT_LABEL_UPPER_LEFT)

        ##makeFakeratePlots(algorithm_discriminator, numerators[algorithm_discriminator]['expr'],
        ##                  denominators_dijet_data, plotters_dijet_data, dijet_dataHLTpaths,
        ##                  denominator_dijet_mc, plotter_dijet_mc,
        ##                  nTuples[algorithm].expr('$jetEta'), binning_eta,
        ##                  'Jet #eta', 6.5e-4, 5.9, extra_labels_eta, extra_labels_diff_eta,
        ##                  maxNumEntries,
        ##                  "plots/jet_fakerates_perHLTpath_jetEta.pdf")

        extra_labels_phi = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_phi.append(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
        extra_labels_diff_phi = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_phi.append(custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)

        ##makeFakeratePlots(algorithm_discriminator, numerators[algorithm_discriminator]['expr'],
        ##                  denominators_dijet_data, plotters_dijet_data, dijet_dataHLTpaths,
        ##                  denominator_dijet_mc, plotter_dijet_mc,
        ##                  nTuples[algorithm].expr('$jetPhi'), binning_phi,
        ##                  'Jet #phi', 6.5e-4, 5.9, extra_labels_phi, extra_labels_diff_phi,
        ##                  maxNumEntries,
        ##                  "plots/jet_fakerates_perHLTpath_jetPhi.pdf")







