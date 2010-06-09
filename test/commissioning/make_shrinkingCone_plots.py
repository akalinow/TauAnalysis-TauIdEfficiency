import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager \
        import TauNtupleManager
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager \
        import PlotManager

# Defintion of input files.
import samples

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")

    # Pull out the defintion of the ntuple (this should be improved)
    # Don't ever access this file directly, just use it to parse event content...
    dummy_file = ROOT.TFile.Open("tauIdEff_ntuple.root", "READ")
    ntuple_manager = TauNtupleManager(dummy_file.Get("Events"), "tauIdEffNtuple")

    # Get the shrinking ntuple
    shrinking_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")

    # Build the plot manager.  The plot manager keeps track of all the samples
    # and ensures they are correctly normalized w.r.t. luminosity.  See 
    # samples.py for available samples.
    data_sample = samples.data
    mc_sample = samples.qcd_mc
   
    plotter = PlotManager()
    # Add each sample we want to plot/compare
    plotter.add_sample(mc_sample, "MC", 
                       fill_color=ROOT.EColor.kBlue-5, draw_option="hist",
                       line_color=ROOT.EColor.kBlue, fill_style=1)

    plotter.add_sample(data_sample, "Data", marker_size=2,
                       marker_color=ROOT.EColor.kRed, draw_option="pe")

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(data_sample.effective_luminosity())

    # Make some plots
    canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    pt_result = plotter.distribution(
        expression=shrinking_ntuple.expr('$pt'),
        selection=shrinking_ntuple.expr('abs($eta) < 2.5'),
        binning = (10, 0, 50),
        x_axis_title = "P_{T}"
    )
    # this should draw a comparison on the canvas, but pt_result
    # now contains some helpful stuff.
    print "MC average pt: %f" % \
            pt_result['samples']['mc_qcd']['mean']

    canvas.SaveAs("test_pt_compare.png")

    denominator = shrinking_ntuple.expr('abs($eta) < 2.5 & $pt > 5')
    numerator = shrinking_ntuple.expr('$byTaNCfrHalfPercent') & denominator

    pt_eff_result = plotter.efficiency(
        expression=shrinking_ntuple.expr('$pt'),
        denominator = denominator,
        numerator = numerator,
        binning = (10, 0, 50),
        x_axis_title = "P_{T}"
    )

    # Add a legend
    pt_eff_result['legend'].make_legend().Draw()

    canvas.SaveAs("test_pt_eff.png")
