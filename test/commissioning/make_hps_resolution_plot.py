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
    ntuple_manager = samples.ztautau_mc.build_ntuple_manager("tauIdEffNtuple")

    # Get the shrinking ntuple
    hps_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeHPS")
        
    theSelection = hps_ntuple.expr('abs($jetEta) < 2.5') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & hps_ntuple.expr('$genMatch > 0.5') & hps_ntuple.expr('$byIsolationMedium > 0.5')


    canvas = ROOT.TCanvas("example", "example", 500, 500)
    
    pt_resol = plotter.distribution(
        expression=hps_ntuple.expr('($pt-$genPt)/$genPt'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(200, -0.5, 0.5),
        x_axis_title = "#tau p_{T} Resolution",
        logy=False
    )
    
    canvas.SaveAs("plots/pt_resol.png")
    canvas.SaveAs("plots/pt_resol.pdf")
    
    eta_resol = plotter.distribution(
        expression=hps_ntuple.expr('$eta-$genEta'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(100, -0.01, 0.01),
        x_axis_title = "#tau #eta Resolution",
        logy=False
    )
    
    canvas.SaveAs("plots/eta_resol.png")
    canvas.SaveAs("plots/eta_resol.pdf")
    
    phi_resol = plotter.distribution(
        expression=hps_ntuple.expr('$phi-$genPhi'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(100, -0.01, 0.01),
        x_axis_title = "#tau #phi Resolution",
        logy=False
    )
    
    canvas.SaveAs("plots/phi_resol.png")
    canvas.SaveAs("plots/phi_resol.pdf")



 

