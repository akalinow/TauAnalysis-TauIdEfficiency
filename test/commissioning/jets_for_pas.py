#!/usr/bin/env python

import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples as samples
import os
import sys

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
 
    # Make some plots
    canvas = ROOT.TCanvas("canvas", "canvas", 500, 500)
    pad = ROOT.TPad("pad", "pad", 0.01, 0.06, 0.99, 0.99)
    pad.SetFillColor(10)
    pad.Draw()
    pad.Divide(1,1)
    pad.cd(1)

    ######################################################
    ####      Plot PFJet distributions                ####
    ######################################################

    # Get the (shrinking cone) PFTau ntuple
    pfTau_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")

    # Basic requirement HLT + Probe object
    pfTau_base_selection = hlt.expr('$hltJet15U > 0.5') & pfTau_ntuple.expr('$probe > 0.5')

    #Compare basic distributions
    pfJetPt_result = plotter.distribution(
        expression = pfTau_ntuple.expr('$jetPt'),
        selection = pfTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & pfTau_base_selection,
        extra_labels = [ style.ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, 0, 100),
        normalize = "data",
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_axis_title = "Number of Events",
        y_min = 1e0, y_max = 1e10, logy = True
    )
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    pfJetPt_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()
    pfJetPt_result['result'].GetYaxis().SetTitleOffset(1.4)

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/pfJetPt.png")
    canvas.SaveAs("plots/pfJetPt.pdf")
  
    pfJetEta_result = plotter.distribution(
        expression = pfTau_ntuple.expr('$jetEta'),
        selection = pfTau_ntuple.expr('$jetPt > 10') & pfTau_base_selection,
        extra_labels = [ style.PT_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -2.5, 2.5),
        normalize = "data",
        x_axis_title = "Jet #eta",
        y_axis_title = "Number of Events",
        y_min = 0, y_max = 60000, logy = False
    )

    pfJetEta_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/pfJetEta.png")
    canvas.SaveAs("plots/pfJetEta.pdf")
  
    pfJetPhi_result = plotter.distribution(
        expression = pfTau_ntuple.expr('$jetPhi'),
        selection = pfTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & pfTau_base_selection,
        extra_labels = [ style.PT_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -3.14, 3.14),
        normalize = "data",
        x_axis_title = "Jet #phi",
        y_axis_title = "Number of Events",
        y_min = 0, y_max = 60000, logy = False
    )
    
    pfJetPhi_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/pfJetPhi.png")
    canvas.SaveAs("plots/pfJetPhi.pdf")

    ######################################################
    ####      Plot Calo/JPTJet distributions          ####
    ######################################################

    # Get the Calo/TCTau ntuple
    caloTau_ntuple = ntuple_manager.get_ntuple(
        "patCaloTausDijetTagAndProbe")

    # Basic requirement HLT + Probe object
    caloTau_base_selection = hlt.expr('$hltJet15U > 0.5') & caloTau_ntuple.expr('$probe > 0.5')

    #Compare basic distributions
    caloJetPt_result = plotter.distribution(
        expression = caloTau_ntuple.expr('$jetPt'),
        selection = caloTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & caloTau_base_selection,
        extra_labels = [ style.ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, 0, 100),
        normalize = "data",
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_axis_title = "Number of Events",
        y_min = 1e0, y_max = 1e10, logy = True
    )
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    caloJetPt_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()
    caloJetPt_result['result'].GetYaxis().SetTitleOffset(1.4)

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/caloJetPt.png")
    canvas.SaveAs("plots/caloJetPt.pdf")
  
    caloJetEta_result = plotter.distribution(
        expression = caloTau_ntuple.expr('$jetEta'),
        selection = caloTau_ntuple.expr('$jetPt > 10') & caloTau_base_selection,
        extra_labels = [ style.PT_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -2.5, 2.5),
        normalize = "data",
        x_axis_title = "Jet #eta",
        y_axis_title = "Number of Events",
        y_min = 0, y_max = 60000, logy = False
    )

    caloJetEta_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/caloJetEta.png")
    canvas.SaveAs("plots/caloJetEta.pdf")
  
    caloJetPhi_result = plotter.distribution(
        expression = caloTau_ntuple.expr('$jetPhi'),
        selection = caloTau_ntuple.expr('$jetPt > 10 && abs($jetEta) < 2.5') & caloTau_base_selection,
        extra_labels = [ style.PT_ETA_CUT_LABEL_UPPER_LEFT ],
        binning = (50, -3.14, 3.14),
        normalize = "data",
        x_axis_title = "Jet #phi",
        y_axis_title = "Number of Events",
        y_min = 0, y_max = 60000, logy = False
    )
    
    caloJetPhi_result['legend'].make_legend(0.60, 0.69, 0.88, 0.88).Draw()

    pad.RedrawAxis()
    canvas.Update()
    pad.Update()
  
    canvas.SaveAs("plots/caloJetPhi.png")
    canvas.SaveAs("plots/caloJetPhi.pdf")
