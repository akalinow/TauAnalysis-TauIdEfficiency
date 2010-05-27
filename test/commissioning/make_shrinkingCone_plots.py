import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager \
        import TauNtupleManager
from TauAnalysis.TauIdEfficiency.ntauples.SampleManager \
        import NtupleSample, NtupleSampleCollection
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager \
        import PlotManager

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")

    # Define some pseudo samples for now
    # Multple weeks of data
    fakeData1 = NtupleSample("data_week1", int_lumi=1.0, prescale=2.0, 
                            files=['tauIdEff_ntuple.root'])

    fakeData2 = NtupleSample("data_week2", int_lumi=2.0, prescale=4.0, 
                            files=['tauIdEff_ntuple.root'])

    # MC binned in pt hat
    fakeMC1 = NtupleSample("mc_pt30", int_lumi=1.0, prescale=1.0,
                          files=['tauIdEff_ntuple.root'])
    fakeMC2 = NtupleSample("mc_pt50", int_lumi=1.0, prescale=1.0,
                          files=['tauIdEff_ntuple.root'])

    # Define the data collections.  Samples are concatenated.
    fakeDataCollection = NtupleSampleCollection(
        name = "data_first_weeks", subsamples = [fakeData1, fakeData2],
        mode = "add")
    # total effective luminosity = sum = 1/2 + 2/4 = 1
    print "Total data int. luminosity", fakeDataCollection.effective_luminosity()

    # Define the MC collection.  Samples are merged (i.e. Pt hat bins)
    fakeMCCollection = NtupleSampleCollection(
        name = "mc_binned", 
        subsamples=[
            fakeMC1, 
            fakeMC2
        ], mode = "merge")

    # Pull out the defintion of the ntuple (this should be improved)
    # Don't ever access this file directly, just use it to parse event content...
    dummy_file = ROOT.TFile("tauIdEff_ntuple.root")
    ntuple_manager = TauNtupleManager(dummy_file.Get("Events"), "tauIdEffNtuple")

    # Get the shrinking ntuple
    shrinking_ntuple = ntuple_manager.get_ntuple(
        "patPFTausDijetTagAndProbeShrinkingCone")

    # Build the plot manager for this selection
    plotter = PlotManager()
    plotter.add_sample(fakeMCCollection, "MC", 
                       fill_color=ROOT.EColor.kBlue-5, draw_option="hist",
                       line_color=ROOT.EColor.kBlue, fill_style=1)

    plotter.add_sample(fakeDataCollection, "Data", marker_size=2,
                       marker_color=ROOT.EColor.kRed, draw_option="pe")

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(fakeDataCollection.effective_luminosity())

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
            pt_result['samples']['mc_binned']['mean']

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
