import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples_cache as samples
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
   
    plotter.add_sample(samples.qcd_mc, "QCD MC", **style.QCD_MC_STYLE_HIST)

#    plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)

    plotter.add_sample(samples.data, "Data (7 TeV)", **style.DATA_STYLE)


    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple maanger
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get the caloTau ntuple
    calo_ntuple = ntuple_manager.get_ntuple(
        "patCaloTausDijetTagAndProbe")

    hlt = ntuple_manager.get_ntuple("TriggerResults")

    # Make some plots_Eff
    canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    # Plot # different triggers
    trigger_results_expr =  hlt.expr('$hltJet15U')

    # Basic requirement HLT + Probe object
    # N.B. currently disabled, no HLT info in ntuples!
    binnin_pt=(20, 0, 100)
    binnin_eta=(25, -2.5, 2.5)
    binnin_phi=(50, -3.14, 3.14)

    base_selection = calo_ntuple.expr('$probe > 0.5') & hlt.expr('$hltJet15U > 0.5')


    denominator = calo_ntuple.expr(
        'abs($jetEta) < 2.5 & $jetPt > 10') & base_selection


    ######Efficiency for Discrimination byLeadTrackFinding

    #Pt#
    pt_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackFinding') & denominator,
        binning = binnin_pt,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )
    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetPt.png")
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetPt.pdf")

    #Eta#
    eta_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackFinding') & denominator,
        binning=binnin_eta,
        x_axis_title = "Jet Eta",
         y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetEta.png")
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetEta.pdf")


    #Phi#
    Phi_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackFinding') & denominator,
        binning=binnin_phi,
        x_axis_title = "Jet Phi",
         y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    Phi_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetPhi.png")
    canvas.SaveAs("plots/caloTau_LeadTrackFinding_eff_jetPhi.pdf")



    ######Efficiency for Discrimination byLeadTrackPtCut

    #Pt#
    pt_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackPtCut && $byLeadTrackFinding') & denominator,
        binning = binnin_pt,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )
    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_eff_jetPt.png")
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_jetPt.pdf")

    #Eta#
    eta_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackPtCut && $byLeadTrackFinding') & denominator,
        binning=binnin_eta,
        x_axis_title = "Jet Eta",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_eff_jetEta.png")
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_jetEta.pdf")


    #Phi#
    Phi_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byLeadTrackPtCut && $byLeadTrackFinding') & denominator,
        binning=binnin_phi,
        x_axis_title = "Jet Phi",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    Phi_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_eff_jetPhi.png")
    canvas.SaveAs("plots/caloTau_LeadTrackPtCut_jetPhi.pdf")


    ######Efficiency for Discrimination byIsolation

    #Pt#
    pt_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator,
        binning = binnin_pt,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )
    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetPt.png")
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetPt.pdf")

    #Eta#
    eta_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator,
        binning=binnin_eta,
        x_axis_title = "Jet Eta",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetEta.png")
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetEta.pdf")


    #Phi#
    Phi_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator,
        binning=binnin_phi,
        x_axis_title = "Jet Phi",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    Phi_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetPhi.png")
    canvas.SaveAs("plots/caloTau_TrkIsolation_eff_jetPhi.pdf")



    ######Efficiency for Discrimination by cut on etSumIsolationECAL

    #Pt#
    pt_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5'),
        binning = binnin_pt,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )
    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_ecalIsolation_eff_jetPt.png")
    canvas.SaveAs("plots/caloTau_ecalIsolation_jetPt.pdf")

    #Eta#
    eta_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5'),
        binning=binnin_eta,
        x_axis_title = "Jet Eta",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_ecalIsolation_eff_jetEta.png")
    canvas.SaveAs("plots/caloTau_ecalIsolation_jetEta.pdf")


    #Phi#
    Phi_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5'),
        binning=binnin_phi,
        x_axis_title = "Jet Phi",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    Phi_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_ecalIsolation_eff_jetPhi.png")
    canvas.SaveAs("plots/caloTau_ecalIsolation_jetPhi.pdf")


    ######Efficiency for Discrimination by number of Tracks in SignalCone

   #Pt#
    pt_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5') & calo_ntuple.expr('$numSignalTracks ==1 | $numSignalTracks ==3') ,
        binning = binnin_pt,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )
    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_prong_eff_jetPt.png")
    canvas.SaveAs("plots/caloTau_prong_eff_jetPt.pdf")

    #Eta#
    eta_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5') & calo_ntuple.expr('$numSignalTracks ==1 | $numSignalTracks ==3') ,
        binning=binnin_eta,
        x_axis_title = "Jet Eta",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_prong_eff_jetEta.png")
    canvas.SaveAs("plots/caloTau_prong_eff_jetEta.pdf")

    #Phi#
    Phi_eff_result = plotter.efficiency(
        expression=calo_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = calo_ntuple.expr('$byIsolation && $byLeadTrackPtCut &&$byLeadTrackFinding') & denominator
        & calo_ntuple.expr('$etSumIsolationECAL < 5') & calo_ntuple.expr('$numSignalTracks ==1 | $numSignalTracks ==3') ,
        binning=binnin_phi,
        x_axis_title = "Jet Phi",
        y_min = 1e-4, y_max = 10, logy = True
    )
    # Add a legend
    Phi_eff_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots/caloTau_prong_eff_jetPhi.png")
    canvas.SaveAs("plots/caloTau_prong_eff_jetPhi.pdf")
    
