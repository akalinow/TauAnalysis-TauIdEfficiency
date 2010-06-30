#!/usr/bin/env python

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Definition of input files.
import samples_cache as samples

plotter = PlotManager()

# Add each sample we want to plot/compare
# Uncomment to add QCD
plotter.add_sample(samples.qcd_mc_pythia8, "Simulation", **style.QCD_MC_PYTHIA8_STYLE_HIST)

#plotter.add_sample(samples.qcd_mc_pythia6, "QCD (Pythia 6)", **style.QCD_MC_PYTHIA6_STYLE_HIST)

#plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)

plotter.add_sample(samples.data, "Data", **style.DATA_STYLE)

# Normalize everything to the data luminosity
plotter.set_integrated_lumi(samples.data.effective_luminosity())

# Build the ntuple manager
ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

shrinking_ntuple = ntuple_manager.get_ntuple(
    "patPFTausDijetTagAndProbeShrinkingCone")

hlt = ntuple_manager.get_ntuple("TriggerResults")

# Make some plots
canvas = ROOT.TCanvas("pas", "pas", 500, 500)
canvas.cd()

from shrinkingConePlots import base_selection, basic_kinematic_cut, pt_binning_fine, \
        eta_binning_fine, phi_binning_fine, lead_pion_selection

# Define the numerators to plot
numerators = [
    {
        'expr': shrinking_ntuple.expr('$byTaNCfrOnePercent') & lead_pion_selection,
        "style_name":"byTaNCfrOnePercent",
        'nice_name': "TaNC 1.00%",
    },
    {
        'expr': shrinking_ntuple.expr('$byTaNCfrHalfPercent') & lead_pion_selection,
        "style_name":"byTaNCfrHalfPercent",
        'nice_name': "TaNC 0.50%"
    },
    {
        'expr': shrinking_ntuple.expr('$byTaNCfrQuarterPercent') & lead_pion_selection,
        "style_name":"byTaNCfrQuarterPercent",
        'nice_name': "TaNC 0.25%"
    },
]

denominator = base_selection & basic_kinematic_cut

# Update the numerator for each efficiency to ensure it is a subset of the
# denominator
for eff_info in numerators: 
    eff_info['expr'] = denominator & eff_info['expr'] 

pt_effs = plotter.multi_efficiency(
    shrinking_ntuple.expr('$jetPt'), 
    denominator,
    numerators, 
    binning=pt_binning_fine, 
    y_min = 1e-4,
    y_max = 10,
    x_axis_title='Jet P_{T}',
    labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
              style.LUMI_LABEL_UPPER_LEFT,
              style.ETA_CUT_LABEL_UPPER_LEFT
              ],
    logy = True)

pt_effs['legend'].make_legend().Draw()

canvas.SaveAs("plots/TaNC_multieff_vs_pt_for_pas.png")
canvas.SaveAs("plots/TaNC_multieff_vs_pt_for_pas.pdf")


eta_effs = plotter.multi_efficiency(
    shrinking_ntuple.expr('$jetEta'), 
    denominator,
    numerators, 
    binning=eta_binning_fine, 
    y_min = 1e-4,
    y_max = 10,
    x_axis_title='Jet #eta',
    labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
              style.LUMI_LABEL_UPPER_LEFT,
              style.PT_CUT_LABEL_UPPER_LEFT
              ],
    logy = True)

eta_effs['legend'].make_legend().Draw()

canvas.SaveAs("plots/TaNC_multieff_vs_eta_for_pas.png")
canvas.SaveAs("plots/TaNC_multieff_vs_eta_for_pas.pdf")

phi_effs = plotter.multi_efficiency(
    shrinking_ntuple.expr('$jetPhi'), 
    denominator,
    numerators, 
    binning=phi_binning_fine, 
    y_min = 1e-4,
    y_max = 10,
    x_axis_title='Jet #phi',
    labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
              style.LUMI_LABEL_UPPER_LEFT,
              style.PT_ETA_CUT_LABEL_UPPER_LEFT
              ],
    logy = True)

phi_effs['legend'].make_legend().Draw()

canvas.SaveAs("plots/TaNC_multieff_vs_phi_for_pas.png")
canvas.SaveAs("plots/TaNC_multieff_vs_phi_for_pas.pdf")









