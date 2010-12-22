#!/usr/bin/env python

'''
Make control plots comparing PFLooseIsolation Pt sums for
 o tau --> e
 o tau --> mu
 o tau --> had.
decays
Authors: Christian Veelken

'''

import ROOT
import samples_cache as samples
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
from TauAnalysis.TauIdEfficiency.ntauples.plotting import draw
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import sys
import copy

maxNumEntries = 1000000000
#maxNumEntries = 10000

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)

    # Get ZTT samples
    ##ztt = samples.ztautau_mc
    ztt = samples.zttPU156bx_mc

    # Get the ntuple we produced
    ntuple_manager = ztt.build_ntuple_manager("tauIdEffNtuple")

    genTaus = ntuple_manager.get_ntuple("tauGenJets")
    genSelection_e = genTaus.expr('$genDecayMode == 0 && $genPt > 10 && abs($genEta) < 2.1')
    genSelection_mu = genTaus.expr('$genDecayMode == 1 && $genPt > 10 && abs($genEta) < 2.1')
    genSelection_had = genTaus.expr('$genDecayMode > 1.5 && $genPt > 10 && abs($genEta) < 2.1')

    recTaus = ntuple_manager.get_ntuple("hpstanc")
    recSelection = recTaus.expr('$byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5')

    recExpression = recTaus.expr('$ptSumLooseIsolation')

    ztt_events = list(ztt.events_and_weights())[0][0]

    canvas = ROOT.TCanvas("canvas", "canvas", 500, 500)
    canvas.SetLeftMargin(0.11)
    canvas.SetGridy(1)

    histogram_e = draw(
        ztt_events,
        expression = recExpression,
        selection = genSelection_e & recSelection,
        binning = (20, 0, 10),
        output_name = 'e',
        maxNumEntries = maxNumEntries
    )
    histogram_e.Sumw2()
    histogram_e.Scale(1./histogram_e.Integral())
    histogram_e.SetMarkerSize(2)
    histogram_e.SetMarkerStyle(20)
    histogram_e.SetMarkerColor(3)
    histogram_e.SetLineColor(3)
    histogram_e.Draw("e1p")

    histogram_mu = draw(
        ztt_events,
        expression = recExpression,
        selection = genSelection_mu & recSelection,
        binning = (20, 0, 10),
        output_name = 'mu',
        maxNumEntries = maxNumEntries
    )
    histogram_mu.Sumw2()
    histogram_mu.Scale(1./histogram_mu.Integral())
    histogram_mu.SetMarkerSize(2)
    histogram_mu.SetMarkerStyle(20)
    histogram_mu.SetMarkerColor(4)
    histogram_mu.SetLineColor(4)
    histogram_mu.Draw("e1psame")

    histogram_had = draw(
        ztt_events,
        expression = recExpression,
        selection = genSelection_had & recSelection,
        binning = (20, 0, 10),
        output_name = 'had',
        maxNumEntries = maxNumEntries
    )
    histogram_had.Sumw2()
    histogram_had.Scale(1./histogram_had.Integral())
    histogram_had.SetMarkerSize(2)
    histogram_had.SetMarkerStyle(20)
    histogram_had.SetMarkerColor(2)
    histogram_had.SetLineColor(2)
    histogram_had.Draw("e1psame")

    # Draw the legend - you can pass NDC xl, yl, xh, yh coordinates to make_legend(...)
    legend = ROOT.TLegend(0.60, 0.69, 0.88, 0.88)
    legend.AddEntry(histogram_e, "#tau #rightarrow e", "p")
    legend.AddEntry(histogram_mu, "#tau #rightarrow #mu", "p")
    legend.AddEntry(histogram_had, "#tau #rightarrow had.", "p")
    legend.Draw()
    
    canvas.SaveAs("./plots/pfLooseIsolation.png")
    canvas.SaveAs("./plots/pfLooseIsolation.pdf")
        
