#!/usr/bin/env python


import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style



# Defintion of input files.
import samples_cache as samples
import os


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
    shrinking_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")
        
    theSelection = shrinking_ntuple.expr('abs($jetEta) < 2.5') & \
            shrinking_ntuple.expr('$byLeadTrackFinding > 0.5') & \
            shrinking_ntuple.expr('$genMatch > 0.5') & \
            shrinking_ntuple.expr('$genDecayMode > 1.5') & \
            shrinking_ntuple.expr('$byTaNCfrQuarterPercent > 0.5')


    canvas = ROOT.TCanvas("example", "example", 500, 500)
    
    pt_resol = plotter.distribution(
        expression=shrinking_ntuple.expr('($pt-$genPt)/$genPt'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(200, -0.5, 0.5),
        x_axis_title = "#tau p_{T} Resolution",
        logy=False
    )

    # Make a pave text w/ mean rms
    stat_label = style.make_mean_rms_pave(pt_resol['samples']['mc_ztt']['plot'])
    stat_label.Draw()
    
    
    canvas.SaveAs("plots/tanc_pt_resolution.png")
    canvas.SaveAs("plots/tanc_pt_resolution.pdf")
    
#    mass_resol = plotter.distribution(
#        expression=shrinking_ntuple.expr('($mass-$genMass)'),
#        selection=theSelection,
#        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
#        binning=(200, -0.5, 0.5),
#        x_axis_title = "#tau mass Resolution",
#        logy=False
#    )
#
#    # Make a pave text w/ mean rms
#    stat_label = style.make_mean_rms_pave(mass_resol['samples']['mc_ztt']['plot'])
#    stat_label.Draw()
#    
#    
#    canvas.SaveAs("plots/tanc_mass_resolution.png")
#    canvas.SaveAs("plots/tanc_mass_resolution.pdf")
    
    eta_resol = plotter.distribution(
        expression=shrinking_ntuple.expr('$eta-$genEta'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(100, -0.01, 0.01),
        x_axis_title = "#tau #eta Resolution",
        n_divisions = 404,
        logy=False
    )
    
    stat_label = style.make_mean_rms_pave(eta_resol['samples']['mc_ztt']['plot'])
    stat_label.Draw()
    
    canvas.SaveAs("plots/tanc_eta_resolution.png")
    canvas.SaveAs("plots/tanc_eta_resolution.pdf")
    
    phi_resol = plotter.distribution(
        expression=shrinking_ntuple.expr('$phi-$genPhi'),
        selection=theSelection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning=(100, -0.01, 0.01),
        x_axis_title = "#tau #phi Resolution",
        n_divisions = 404,
        logy=False
    )
    
    stat_label = style.make_mean_rms_pave(phi_resol['samples']['mc_ztt']['plot'])
    stat_label.Draw()
    
    canvas.SaveAs("plots/tanc_phi_resolution.png")
    canvas.SaveAs("plots/tanc_phi_resolution.pdf")

 

