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
    # Uncomment to add QCD
    plotter.add_sample(samples.ztautau_mc, "Z->#tau#tau MC", **style.QCD_MC_STYLE_HIST)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple maanger
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get the shrinking ntuple
    hps_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeHPS")


    canvas = ROOT.TCanvas("example", "example", 500, 500)
    
    pt_resol = plotter.distribution(
        expression=hps_ntuple.expr(('$jetPt-$genPt)/$genPt'),
        selection=hps_ntuple.expr('abs($jetEta) < 2.5') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & hps_ntuple.expr('$genMatch > 0.5'),
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(50, -1, 1),
        x_axis_title = "#tau p_{T} Resolution",
        logy=False
    )
    
    canvas.SaveAs("plots/pt_resol.png")
    canvas.SaveAs("plots/pt_resol.pdf")


 

