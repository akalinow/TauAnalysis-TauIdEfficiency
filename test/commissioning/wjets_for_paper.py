#!/usr/bin/env python

import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples_cache as samples
import copy
import os
import sys

# Define sample names
sample_wjets_data = samples.data_wjets
sample_wjets_mc   = samples.wjetsPU156bx_mc

# style Definitions
custom_style_wjets_data = copy.deepcopy(style.DATA_STYLE)
custom_style_wjets_mc = copy.deepcopy(style.QCD_MC_STYLE_HIST)

custom_CMS_PRELIMINARY_UPPER_LEFT = copy.deepcopy(style.CMS_PRELIMINARY_UPPER_LEFT)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetX2(0.60)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetY1(0.850)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.055)

custom_LUMI_LABEL_UPPER_LEFT = copy.deepcopy(style.LUMI_LABEL_UPPER_LEFT)
custom_LUMI_LABEL_UPPER_LEFT.SetX2(0.60)
custom_LUMI_LABEL_UPPER_LEFT.SetY2(0.840)
custom_LUMI_LABEL_UPPER_LEFT.SetY1(0.790)
custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.045)
custom_LUMI_LABEL_UPPER_LEFT.Clear()
custom_LUMI_LABEL_UPPER_LEFT.AddText("Data, L = 36.1pb^{-1}")
custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.0375)

custom_SQRTS_LABEL_UPPER_LEFT = copy.deepcopy(style.SQRTS_LABEL_UPPER_LEFT)
custom_SQRTS_LABEL_UPPER_LEFT.SetX2(0.60)
custom_SQRTS_LABEL_UPPER_LEFT.SetY2(0.785)
custom_SQRTS_LABEL_UPPER_LEFT.SetY1(0.735)
custom_SQRTS_LABEL_UPPER_LEFT.SetTextSize(0.045)

custom_PT_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_CUT_LABEL_UPPER_LEFT)
custom_PT_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
custom_PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)
custom_PT_CUT_LABEL_UPPER_LEFT.Clear()
custom_PT_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c")

custom_ETA_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.ETA_CUT_LABEL_UPPER_LEFT)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)
custom_ETA_CUT_LABEL_UPPER_LEFT.Clear()
custom_ETA_CUT_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")

custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX2(0.60)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY2(0.740)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY1(0.640)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.045)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.Clear()
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c,\n")
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")

custom_DEFAULT_LABELS = [
    custom_CMS_PRELIMINARY_UPPER_LEFT,
    custom_LUMI_LABEL_UPPER_LEFT,
    custom_SQRTS_LABEL_UPPER_LEFT
]

def makeJetPlot(expression, selection, extra_labels,
                binning, x_axis_title, y_min, y_max, logy,
                filename):

    #----------------------------------------------------------------------------
    # Plot distribution of Data vs. Monte Carlo predictions
    #----------------------------------------------------------------------------

    # Create canvas
    canvas_distribution = ROOT.TCanvas("canvas_distribution", "canvas_distribution", 500, 500)     
    canvas_distribution.SetLeftMargin(0.12)
    canvas_distribution.SetBottomMargin(0.12)

    canvas_distribution.cd()

    # Create plot
    distribution_result = plotter.distribution(
        expression = expression,
        selection = selection,
        extra_labels = extra_labels,
        binning = binning,
        normalize = sample_wjets_data.name,
        x_axis_title = x_axis_title,
        y_axis_title = "Number of Events",
        y_min = y_min, y_max = y_max, logy = logy
    )

    # Adjust offset of y-axis title
    distribution_result['result'].GetYaxis().SetTitleOffset(1.4)
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coordinates to make_legend(...)
    legend_x1 = 0.60
    legend_y1 = 0.69
    legend_x2 = 0.88
    legend_y2 = 0.88
    distribution_result['legend'].make_legend(legend_x1, legend_y1, legend_x2, legend_y2).Draw()

    canvas_distribution.cd()
    canvas_distribution.Update()
    canvas_distribution.SaveAs("plots/" + filename + ".png")
    canvas_distribution.SaveAs("plots/" + filename + ".pdf")

    canvas_distribution.IsA().Destructor(canvas_distribution)

    #----------------------------------------------------------------------------
    # Add plot of (normalized) difference (Data - MC)/MC
    #----------------------------------------------------------------------------

    # Create new canvas
    canvas_diff = ROOT.TCanvas("canvas_diff", "canvas_diff", 500, 625)

    canvas_diff.cd()
    
    topPad_diff = ROOT.TPad("topPad_diff", "topPad_diff", 0.00, 0.35, 1.00, 1.00)
    topPad_diff.SetFillColor(10)
    topPad_diff.SetLogy(logy)
    topPad_diff.SetTopMargin(0.05)
    topPad_diff.SetLeftMargin(0.18)
    topPad_diff.SetBottomMargin(0.00)
    topPad_diff.SetRightMargin(0.05)
    topPad_diff.SetGridy()
    topPad_diff.Draw()
    topPad_diff.cd()

    distribution_result['result'].Draw("nostack")

    # Ommit x-axis title
    distribution_result['result'].GetXaxis().SetTitle("")
    
    distribution_result['result'].GetYaxis().SetTitleSize(0.05) 
    distribution_result['result'].GetYaxis().SetLabelSize(0.05)
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    legend_adjusted_coordinates = style.adjust_coordinates(topPad_diff, legend_x1, legend_y1, legend_x2, legend_y2)
    distribution_result['legend'].make_legend(
        legend_adjusted_coordinates['x1'], legend_adjusted_coordinates['y1'],
        legend_adjusted_coordinates['x2'], legend_adjusted_coordinates['y2']).Draw()

    # Define temporary disctionary to prevent automatic garbage collection
    to_keep = {}

    draw_labels_topPad = style.draw_labels(topPad_diff)
    to_keep["1"] = draw_labels_topPad(custom_DEFAULT_LABELS)
    to_keep["2"] = draw_labels_topPad(extra_labels)

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

    # Make plot of (normalized) difference (Data - MC)/MC
    diff_result = plotter.plot_dist_deviations(
        distribution_result,
        sample_wjets_data.name,
        [ sample_wjets_mc.name ],
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.28, y_max = +0.28, logy = False
    )

    diff_result['background'].GetYaxis().SetNdivisions(505)

    diff_result['background'].GetXaxis().SetTitleOffset(1.1)
    diff_result['background'].GetXaxis().SetTitleSize(0.08) 
    diff_result['background'].GetXaxis().SetLabelSize(0.08)
    diff_result['background'].GetXaxis().SetTickLength(0.055)

    diff_result['background'].GetYaxis().SetTitleOffset(0.9)
    diff_result['background'].GetYaxis().SetTitleSize(0.08) 
    diff_result['background'].GetYaxis().SetLabelSize(0.08)
    diff_result['background'].GetYaxis().SetTickLength(0.04)

    canvas_diff.cd() 
    canvas_diff.Update()    
    canvas_diff.SaveAs("plots/" + filename + "_diff.png")
    canvas_diff.SaveAs("plots/" + filename + "_diff.pdf")

    canvas_diff.IsA().Destructor(canvas_diff)

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    if not os.path.isdir('plots'):
        os.mkdir('plots')

    # Build the plot manager.  The plot manager keeps track of all the samples
    # and ensures they are correctly normalized w.r.t. luminosity.  See 
    # samples.py for available samples.
    plotter = PlotManager()

    # Add each sample we want to plot/compare
    plotter.add_sample(sample_wjets_data, "W #rightarrow #mu #nu Data",       **custom_style_wjets_data)
    plotter.add_sample(sample_wjets_mc,   "W #rightarrow #mu #nu Simulation", **custom_style_wjets_mc)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(sample_wjets_data.effective_luminosity())

    # Build the ntuple manager
    ntuple_manager = sample_wjets_data.build_ntuple_manager("tauIdEffNtuple")

    # Get HLT trigger decisions
    hlt_ntuple = ntuple_manager.get_ntuple("patTriggerEvent")
 
    ######################################################
    ####      Plot PFJet distributions                ####
    ######################################################

    # Get the (shrinking cone) PFTau ntuple
    pfTau_ntuple = ntuple_manager.get_ntuple("shrinking")
    muTau_ntuple = ntuple_manager.get_ntuple("patMuonPFTauPairsShrinkingCone")

    # Basic requirement HLT + Probe object
    pfTau_selection_phase_space = pfTau_ntuple.expr("$jetPt > 20.0 & abs($jetEta) < 2.3")
    pfTau_selection_e_mu_veto   = pfTau_ntuple.expr("$pfElectronMVA < 0.6 & $againstMuon > 0.5")
    pfTau_selection_wjets       = pfTau_selection_phase_space & pfTau_selection_e_mu_veto
    
    # Make plots of jetPt, jetEta and jetPhi
    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
        selection = pfTau_selection_wjets,
        extra_labels = [ custom_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 6.5e-1, y_max = 1.e+8, logy = True,
        filename = "wjets_pfJetPt"
    )

    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
        selection = pfTau_selection_wjets,
        extra_labels = [ custom_PT_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "wjets_pfJetEta"
    )

    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetPhi'),
        selection = pfTau_selection_wjets,
        extra_labels = [ custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (50, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "wjets_pfJetPhi"
    )

    makeJetPlot(
        expression = muTau_ntuple.expr('$Mt'),
        selection = pfTau_selection_wjets,
        extra_labels = [ custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (20, 0, 100),
        x_axis_title = "E_{T}^{miss}",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "wjets_Mt"
    )
    
    makeJetPlot(
        expression = muTau_ntuple.expr('$ptW'),
        selection = pfTau_selection_wjets,
        extra_labels = [ custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (20, 0, 100),
        x_axis_title = "P_{T}^{W}",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "wjets_ptW"
    )
    

    
