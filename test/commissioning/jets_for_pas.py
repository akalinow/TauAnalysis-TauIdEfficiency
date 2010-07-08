#!/usr/bin/env python

import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples_cache as samples
import os
import sys

# style Definitions
jets_CMS_PRELIMINARY_UPPER_LEFT = style.CMS_PRELIMINARY_UPPER_LEFT.Clone()
jets_CMS_PRELIMINARY_UPPER_LEFT.SetX2(0.60)
jets_CMS_PRELIMINARY_UPPER_LEFT.SetY1(0.850)
jets_CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.055)

jets_LUMI_LABEL_UPPER_LEFT = style.LUMI_LABEL_UPPER_LEFT.Clone()
jets_LUMI_LABEL_UPPER_LEFT.SetX2(0.60)
jets_LUMI_LABEL_UPPER_LEFT.SetY2(0.840)
jets_LUMI_LABEL_UPPER_LEFT.SetY1(0.790)
jets_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_SQRTS_LABEL_UPPER_LEFT = style.SQRTS_LABEL_UPPER_LEFT.Clone()
jets_SQRTS_LABEL_UPPER_LEFT.SetX2(0.60)
jets_SQRTS_LABEL_UPPER_LEFT.SetY2(0.785)
jets_SQRTS_LABEL_UPPER_LEFT.SetY1(0.735)
jets_SQRTS_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_PT_CUT_LABEL_UPPER_LEFT = style.PT_CUT_LABEL_UPPER_LEFT.Clone()
jets_PT_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
jets_PT_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
jets_PT_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
jets_PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_ETA_CUT_LABEL_UPPER_LEFT = style.ETA_CUT_LABEL_UPPER_LEFT.Clone()
jets_ETA_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.Clone()
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX2(0.60)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY2(0.740)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY1(0.640)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_DEFAULT_LABELS = [
    jets_CMS_PRELIMINARY_UPPER_LEFT,
    jets_LUMI_LABEL_UPPER_LEFT,
    jets_SQRTS_LABEL_UPPER_LEFT
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
        normalize = "data",
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
    to_keep["1"] = draw_labels_topPad(jets_DEFAULT_LABELS)
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
        "data",
        [ "mc_qcd_pythia8" ],
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
    plotter.add_sample(samples.qcd_mc_pythia8, "Simulation", **style.QCD_MC_STYLE_HIST)
    plotter.add_sample(samples.data, "Data", **style.DATA_STYLE)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple manager
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get HLT trigger decisions
    hlt = ntuple_manager.get_ntuple("TriggerResults")
 
    ######################################################
    ####      Plot PFJet distributions                ####
    ######################################################

    # Get the (shrinking cone) PFTau ntuple
    pfTau_ntuple = ntuple_manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone")

    # Basic requirement HLT + Probe object
    pfTau_base_selection = hlt.expr('$hltJet15U > 0.5') & pfTau_ntuple.expr('$probe > 0.5')

    # Make plots of jetPt, jetEta and jetPhi
    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
        selection = pfTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & pfTau_base_selection,
        extra_labels = [ jets_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 6.5e-1, y_max = 1.e+8, logy = True,
        filename = "pfJetPt"
    )

    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
        selection = pfTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & pfTau_base_selection,
        extra_labels = [ jets_PT_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "pfJetEta"
    )

    makeJetPlot(
        expression = pfTau_ntuple.expr('$jetPhi'),
        selection = pfTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (50, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "pfJetPhi"
    )

    ######################################################
    ####      Plot Calo/JPTJet distributions          ####
    ######################################################

    # Get the Calo/TCTau ntuple
    caloTau_ntuple = ntuple_manager.get_ntuple("patCaloTausDijetTagAndProbe")

    # Basic requirement HLT + Probe object
    caloTau_base_selection = hlt.expr('$hltJet15U > 0.5') & caloTau_ntuple.expr('$probe > 0.5')

    # Make plots of jetPt, jetEta and jetPhi
    makeJetPlot(
        expression = caloTau_ntuple.expr('$jetPt'),
        selection = caloTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & caloTau_base_selection,
        extra_labels = [ jets_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 6.5e-1, y_max = 1.e+8, logy = True,
        filename = "caloJetPt"
    )

    makeJetPlot(
        expression = caloTau_ntuple.expr('$jetEta'),
        selection = caloTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & caloTau_base_selection,
        extra_labels = [ jets_PT_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "caloJetEta"
    )

    makeJetPlot(
        expression = caloTau_ntuple.expr('$jetPhi'),
        selection = caloTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & caloTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = -4000, y_max = 50000, logy = False,
        filename = "caloJetPhi"
    )
