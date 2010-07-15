import sys

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import os
import samples_cache as samples

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    if not os.path.isdir('plots_Control'):
        os.mkdir('plots_Control')

    # Build the plot manager.  The plot manager keeps track of all the samples
    # and ensures they are correctly normalized w.r.t. luminosity.  See 
    # samples.py for available samples.
    plotter = PlotManager()

    plotter.add_sample(samples.qcd_mc_pythia6, "PYTHIA 6", ** style.QCD_MC_PYTHIA6_STYLE_DOTS)

    plotter.add_sample(samples.qcd_mc_pythia8, "PYTHIA 8", ** style.QCD_MC_PYTHIA8_STYLE_HIST)

    plotter.add_sample(samples.data, "Data", ** style.DATA_STYLE)


    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple maanger
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get the shrinking ntuple
    calo_ntuple = ntuple_manager.get_ntuple("patCaloTausDijetTagAndProbe")

    hlt = ntuple_manager.get_ntuple("TriggerResults")

    # Make some plots_Control
    canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    # Plot # different triggers
    trigger_results_expr = hlt.expr('$hltJet15U')

    trigger_results = plotter.distribution(
                                           expression=hlt.expr('$hltJet15U'),
                                           selection=hlt.expr('1'), # no selection
                                           binning=(5, -2.5, 2.5),
                                           #        binning = (0, 5,10,15,25,40,60,80,110),
                                           x_axis_title="HLT_Jet15U Result",
                                           y_min=1, logy=True,
                                           # Can pass a list of TPaveTexts to draw on the plot
                                           # CMS preliminary in the upper left is the default
                                           #labels = [styles.CMS_PRELIMINARY_UPPER_LEFT]
                                           )

    trigger_results['legend'].make_legend().Draw()

    canvas.SaveAs("plots_Control/hltJet15U_result.png")
    canvas.SaveAs("plots_Control/hltJet15U_result.pdf")

    # Basic requirement HLT + Probe object
    # N.B. currently disabled, no HLT info in ntuples!
    base_selection = calo_ntuple.expr('$probe > 0.5') & hlt.expr('$hltJet15U > 0.5')


    ##########################################################################################

    jetpt_result = plotter.distribution(
                                        expression=calo_ntuple.expr('$jetPt'),
                                        selection=calo_ntuple.expr('abs($jetEta) < 2.5') & base_selection,
                                        binning=(50, 0, 100),
                                        x_axis_title="Jet P_{T} [GeV/c]",
                                        y_min=1e-1, logy=True
                                        )

    jetpt_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_jetPt.png")
    canvas.SaveAs("plots_Control/caloTau_jetPt.pdf")
    ##########################################################################################

    jetEEta_result = plotter.distribution(
                                          expression=calo_ntuple.expr('$jetEta'),
                                          selection=calo_ntuple.expr('$jetPt > 10') & base_selection,
                                          binning=(50, -2.5, 2.5),
                                          x_axis_title="Jet #eta",
                                          y_min=1e-4, logy=False
                                          )
    jetEEta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_jetEta.png")
    canvas.SaveAs("plots_Control/caloTau_jetEta.pdf")
    ##########################################################################################

    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$jetPhi'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(50, -3.14, 3.14),
                                         x_axis_title="Jet #phi",
                                         y_min=1e-1, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_jetPhi.png")
    canvas.SaveAs("plots_Control/caloTau_jetPhi.pdf")
    ##########################################################################################


    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$invariantMassSignalTracks'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(50, 0, 5),
                                         x_axis_title="invariantMassSignalTracks [GeV]",
                                         y_min=1e-1, y_max=1e8, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_invariantMassSignalTracks.png")
    canvas.SaveAs("plots_Control/caloTau_invariantMassSignalTracks.pdf")
    ##########################################################################################

    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$numSignalTracks'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(15, -0.5, 14.5),
                                         x_axis_title="numSignalTracks",
                                         y_min=1e-1, y_max=1e8, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_numSignalTracks.png")
    canvas.SaveAs("plots_Control/caloTau_numSignalTracks.pdf")

    ##########################################################################################
    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$numIsolationTracks'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(25, -0.5, 24.5),
                                         x_axis_title="numIsolationTracks",
                                         y_min=1e-1, y_max=1e8, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_numIsolationTracks.png")
    canvas.SaveAs("plots_Control/caloTau_numIsolationTracks.pdf")

##########################################################################################

    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$etSumIsolationECAL'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(80, 0, 40),
                                         x_axis_title="etSumIsolationECAL [GeV]",
                                         y_min=1e-1, y_max=1e8, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()
    canvas.SaveAs("plots_Control/caloTau_etSumIsolationECAL.png")
    canvas.SaveAs("plots_Control/caloTau_etSumIsolationECAL.pdf")
##########################################################################################

    jeteta_result = plotter.distribution(
                                         expression=calo_ntuple.expr('$ptSumIsolationTracks'),
                                         selection=calo_ntuple.expr('$jetPt > 10 & abs($jetEta) < 2.5 ') & base_selection,
                                         binning=(80, 0, 40),
                                         x_axis_title="ptSumIsolationTracks [GeV/c]",
                                         y_min=1e-1, y_max=1e8, logy=True
                                         )
    jeteta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots_Control/caloTau_ptSumIsolationTracks.png")
    canvas.SaveAs("plots_Control/caloTau_ptSumIsolationTracks.pdf")
    ##########################################################################################
