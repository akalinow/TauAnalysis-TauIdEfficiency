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
    plotter.add_sample(samples.qcd_mc, "QCD MC", **style.QCD_MC_STYLE_HIST)

    #plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)

    plotter.add_sample(samples.data, "Data (7 TeV)", **style.DATA_STYLE)


    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple maanger
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get the shrinking ntuple
    hps_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")
        
    print hps_ntuple

    hlt = ntuple_manager.get_ntuple("TriggerResults")

    # Make some plots
    canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    #Plot different triggers
    trigger_results_expr =  hlt.expr('$hltJet15U')

    trigger_results = plotter.distribution(
        expression=hlt.expr('$hltJet15U'),
        selection=hlt.expr('1'), # no selection
        binning = (2, 0, 2),
        x_axis_title = "HLT_Jet15U Result",
        y_min = 1
        # Can pass a list of TPaveTexts to draw on the plot
        # CMS preliminary in the upper left is the default
        #labels = [styles.CMS_PRELIMINARY_UPPER_LEFT]  
    )

    trigger_results['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hltJet15U_result.png")
    canvas.SaveAs("plots/hltJet15U_result.pdf")

    # Basic requirement HLT + Probe object
    # N.B. currently disabled, no HLT info in ntuples!
    base_selection = hlt.expr('$hltJet15U > 0.5') & hps_ntuple.expr('$probe > 0.5')
    #base_selection = shrinking_ntuple.expr('1')

    # Compare basic distributions
    jetpt_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetPt'),
        selection=hps_ntuple.expr('abs($jetEta) < 2.5') & base_selection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e2,  y_max = 4e8, logy=True
    )

    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    jetpt_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/jetPt.png")
    canvas.SaveAs("plots/jetPt.pdf")

    jeteta_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetEta'),
        selection=hps_ntuple.expr('$jetPt > 10') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1, y_max = 80000, logy=False
    )
    jeteta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/jetEta.png")
    canvas.SaveAs("plots/jetEta.pdf")

    jetphi_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetPhi'),
        selection=hps_ntuple.expr('$jetPt > 10') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (50, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1, y_max = 80000
    )
    jeteta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/jetPhi.png")
    canvas.SaveAs("plots/jetPhi.pdf")

    # Compare basic distributions
    taupt_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetPt'),
        selection=hps_ntuple.expr('abs($jetEta) < 2.5') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.ETA_CUT_LABEL_UPPER_LEFT],
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1, y_max = 12000, logy=False
    )

    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    taupt_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsPt.png")
    canvas.SaveAs("plots/hpsPt.pdf")

    taueta_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetEta'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1, y_max = 5000
    )
    taueta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsEta.png")
    canvas.SaveAs("plots/hpsEta.pdf")

    tauphi_result = plotter.distribution(
        expression=hps_ntuple.expr('$jetPhi'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (50, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1, y_max = 5000
    )
    taueta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsPhi.png")
    canvas.SaveAs("plots/hpsPhi.pdf")
    
    tauIsoCands_result = plotter.distribution(
        expression=hps_ntuple.expr('$numParticlesIsoCone'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (50, 0, 50),
        x_axis_title = "Num Particles Iso Cone",
        y_min = 1, y_max = 20000
    )
    tauIsoCands_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/tauIsoCands.png")
    canvas.SaveAs("plots/tauIsoCands.pdf")

    tauSigCands_result = plotter.distribution(
        expression=hps_ntuple.expr('$numChargedParticlesSignalCone'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (6, 0, 5),
        x_axis_title = "Num Charged Particles Signal Cone",
        y_min = 1, y_max = 120000
    )
    tauIsoCands_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/tauSigCands.png")
    canvas.SaveAs("plots/tauSigCands.pdf")

    mass1prongWpi0_result = plotter.distribution(
        expression=hps_ntuple.expr('$mass'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$numChargedParticlesSignalCone == 1') & hps_ntuple.expr('$numPhotonsSignalCone > 0') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (15, 0, 1.4),
        x_axis_title = "One Prong + #pi0 Mass",
        y_min = 1, y_max = 16000
    )
    mass1prongWpi0_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/mass1prongWpi0.png")
    canvas.SaveAs("plots/mass1prongWpi0.pdf")

    mass3prong_result = plotter.distribution(
        expression=hps_ntuple.expr('$mass'),
        selection=hps_ntuple.expr('$jetPt > 10') & hps_ntuple.expr('$numChargedParticlesSignalCone == 3') & hps_ntuple.expr('$byLeadTrackFinding > 0.5') & base_selection,
        extra_labels = [style.PT_CUT_LABEL_UPPER_LEFT],
        binning = (11, 0.6, 1.6),
        x_axis_title = "Three Prong Mass",
        y_min = 1, y_max = 12000
    )
    mass3prong_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/mass3prong.png")
    canvas.SaveAs("plots/mass3prong.pdf")

    ######################################################
    ####      Plot efficiencies                       ####
    ######################################################

    # Change the style of the QCD from filled histogram to dots
    # name mc_qcd is defined in samples.py
    plotter.update_style("mc_qcd", **style.QCD_MC_STYLE_DOTS)
    pt_binning_fine = (0, 5, 10, 15, 20, 25, 35, 45, 60, 80, 100)
    

    denominator = hps_ntuple.expr('abs($jetEta) < 2.5 & $jetPt > 10') & base_selection
    numerator = hps_ntuple.expr('$byLeadTrackFinding > 0.5') & denominator

#     eta_eff_result = plotter.efficiency(
#         expression=hps_ntuple.expr('abs($jetEta)'),
#         denominator = denominator,
#         numerator = numerator,
#         binning = (25, 0, 2.5),
#         x_axis_title = "Jet |#eta|",
#         y_min = 1e-4, y_max = 5, logy = True,
#     )
# 
#     # Add a legend
#     eta_eff_result['legend'].make_legend().Draw()
# 
#     canvas.SaveAs("plots/hpsCone_TaNCHalf_eff_jetEta.png")
#     canvas.SaveAs("plots/hpsCone_TaNCHalf_eff_jetEta.pdf")

    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = pt_binning_fine,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 10, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetPt.png")
    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetPt.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (25, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1e-3, y_max = 10, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetEta.png")
    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetEta.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (31, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1e-3, y_max = 10, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetPhi.png")
    canvas.SaveAs("plots/hpsCone_decayMode_eff_jetPhi.pdf")


    numerator = hps_ntuple.expr('$byLeadTrackFinding > 0.5') & hps_ntuple.expr('$byIsolationLoose > 0.5') & denominator
 
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = pt_binning_fine,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetPt.png")
    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetPt.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (25, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetEta.png")
    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetEta.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (31, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetPhi.png")
    canvas.SaveAs("plots/hpsCone_decayModeWLooseIso_eff_jetPhi.pdf")


    numerator = hps_ntuple.expr('$byLeadTrackFinding > 0.5') & hps_ntuple.expr('$byIsolationMedium > 0.5') & denominator
 
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = pt_binning_fine,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetPt.png")
    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetPt.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (25, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetEta.png")
    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetEta.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (31, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetPhi.png")
    canvas.SaveAs("plots/hpsCone_decayModeWMediumIso_eff_jetPhi.pdf")
    
    numerator = hps_ntuple.expr('$byLeadTrackFinding > 0.5') & hps_ntuple.expr('$byIsolationTight > 0.5') & denominator
 
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = pt_binning_fine,
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-5, y_max = 1, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetPt.png")
    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetPt.pdf")

    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetEta'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (25, -2.5, 2.5),
        x_axis_title = "Jet #eta",
        y_min = 1e-5, y_max = 1, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetEta.png")
    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetEta.pdf")
    
    pt_eff_result = plotter.efficiency(
        expression=hps_ntuple.expr('$jetPhi'),
        denominator = denominator,
        numerator = numerator,
        extra_labels = [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        binning = (31, -3.14, 3.14),
        x_axis_title = "Jet #phi",
        y_min = 1e-5, y_max = 1, logy = True,
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetPhi.png")
    canvas.SaveAs("plots/hpsCone_decayModeWTightIso_eff_jetPhi.pdf")