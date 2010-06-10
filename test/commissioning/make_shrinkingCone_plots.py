import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager \
        import TauNtupleManager
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager \
        import PlotManager

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
    #plotter.add_sample(samples.qcd_mc, "QCD MC", 
    #                   fill_color=ROOT.EColor.kBlue-5, draw_option="hist",
    #                   line_color=ROOT.EColor.kBlue, fill_style=1)

    plotter.add_sample(samples.minbias_mc, "Minbias MC", 
                       fill_color=ROOT.EColor.kGreen-5, draw_option="pe",
                       marker_size=2, line_color=ROOT.EColor.kBlue, fill_style=0)

    plotter.add_sample(samples.data, "Data (7 TeV)", marker_size=2,
                       marker_color=ROOT.EColor.kRed, draw_option="pe")


    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())

    # Build the ntuple maanger
    ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

    # Get the shrinking ntuple
    shrinking_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")

    # Make some plots
    canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    # Compare basic distributions
    jetpt_result = plotter.distribution(
        expression=shrinking_ntuple.expr('$jetPt'),
        selection=shrinking_ntuple.expr('abs($jetEta) < 2.5'),
        binning = (50, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-2, logy=True
    )

    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    jetpt_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/shrinkingCone_jetPt.png")
    canvas.SaveAs("plots/shrinkingCone_jetPt.pdf")

    jeteta_result = plotter.distribution(
        expression=shrinking_ntuple.expr('$jetEta'),
        selection=shrinking_ntuple.expr('abs($jetPt) > 5'),
        binning = (50, -2.5, 2.5),
        x_axis_title = "Jet #eta"
    )
    jeteta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/shrinkingCone_jetEta.png")
    canvas.SaveAs("plots/shrinkingCone_jetEta.pdf")

    jeteta_result = plotter.distribution(
        expression=shrinking_ntuple.expr('$jetPhi'),
        selection=shrinking_ntuple.expr('abs($jetPt) > 5'),
        binning = (50, 0, 2*3.14),
        x_axis_title = "Jet #phi"
    )
    jeteta_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/shrinkingCone_jetPhi.png")
    canvas.SaveAs("plots/shrinkingCone_jetPhi.pdf")

    denominator = shrinking_ntuple.expr('abs($jetEta) < 2.5 & $jetPt > 5')
    numerator = shrinking_ntuple.expr('$byTaNCfrHalfPercent') & denominator

    eta_eff_result = plotter.efficiency(
        expression=shrinking_ntuple.expr('abs($jetEta)'),
        denominator = denominator,
        numerator = numerator,
        binning = (50, 0, 2.5),
        x_axis_title = "Jet |#eta|",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("plots/shrinkingCone_TaNCHalf_eff_jetPt.png")
    canvas.SaveAs("plots/shrinkingCone_TaNCHalf_eff_jetPt.pdf")

    eta_eff_result = plotter.efficiency(
        expression=shrinking_ntuple.expr('$jetPt'),
        denominator = denominator,
        numerator = numerator,
        binning = (100, 0, 100),
        x_axis_title = "Jet P_{T} [GeV/c]",
        y_min = 1e-4, y_max = 5, logy = True,
    )

    # Add a legend
    eta_eff_result['legend'].make_legend().Draw()

