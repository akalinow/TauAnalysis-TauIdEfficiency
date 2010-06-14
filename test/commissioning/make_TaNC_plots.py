import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager \
        import TauNtupleManager
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager \
        import PlotManager

# Defintion of input files.
import samples_cache as samples
import os

#do not like the code replication but other options are just to involved
def makePtPlots( canvas, name = "Charged", objects ="Outlier" , count = 4, selection = "abs($jetEta) < 2.5 & $jetPt > 5", formats = ["png","pdf"], subdir="plots" ):
    result = {}
    for number in range(count):
        # Compare basic distributions
        try:
            result["%s%sPt%s"%(name,objects,number)] = plotter.distribution(
                expression=shrinking_ntuple.expr("$TaNC%s%sPt%s"%(name,objects,number)),
                selection=shrinking_ntuple.expr(selection),
                binning = (100, 0, 16),
                x_axis_title = "p_t of %s %s No. %s [GeV/c]"%(name,objects,number),
                y_min = 1e-2, logy=True
                )
            # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
            result["%s%sPt%s"%(name,objects,number)]['legend'].make_legend().Draw()
            
            for format in formats:
                canvas.SaveAs("plots/shrinkingCone_%sOutlierPt%s.%s"%(name,number,format))
        except StandardError, err:
            print "TaNC%s%sPt%s not found: %s"%(name,objects,number,err)
    
    return result

def makeAnglePlots( canvas, name = "Charged", objects ="Outlier" , count = 4, selection = "abs($jetEta) < 2.5 & $jetPt > 5", formats = ["png","pdf"], subdir="plots"):
    result = {}    
    for number in range(count):
        # Compare basic distributions
        try:
            result["%s%sAngle%s"%(name,objects,number)] = plotter.distribution(
                expression=shrinking_ntuple.expr("$TaNC%s%sAngle%s"%(name,objects,number)),
                selection=shrinking_ntuple.expr(selection),
                binning = (100, 0, 1.5),
                x_axis_title = "#DeltaR %s %s No. %s"%(name,objects,number),
                y_min = 1e-2, logy=True
                )
            
            # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
            result["%s%sAngle%s"%(name,objects,number)]['legend'].make_legend().Draw()
            
            for format in formats:
                canvas.SaveAs("%s/shrinkingCone_%s%sAngle%s.%s"%(subdir,name,objects,number,format))
        except StandardError, err:
            print "TaNC%s%sAngle%s not found: %s"%(name,objects,number,err)
    return result

def makeKiniticPlots(canvas, selection = "abs($jetEta) < 2.5 & $jetPt > 5", formats = ["png","pdf"], subdir="plots"):
    result["Pt"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCPt'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 70),
        x_axis_title = "#tau p_{T} [GeV/c]",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["Pt"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_Pt.%s"%(subdir,format))
    
    result["Eta"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCEta'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 2.5),
        x_axis_title = "|#eta|",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["Eta"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_Eta.%s"%(subdir,format))

    result["MainTrackAngle"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCMainTrackAngle'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 1),
        x_axis_title = "#DeltaR of main Track",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["MainTrackAngle"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_MainTrackAngle.%s"%(subdir,format))

    result["MainTrackPt"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCMainTrackPt'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 70),
        x_axis_title = "p_{T} of main Track",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["MainTrackPt"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_MainTrackPt.%s"%(subdir,format))

    result["InvariantMassOfSignal"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCInvariantMassOfSignal'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 2.5),
        x_axis_title = "inv. mass of signal particles",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["InvariantMassOfSignal"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_InvariantMassOfSignal.%s"%(subdir,format))

        

    for dalitz in range(2):
        result["Dalitz%s"%dalitz] = plotter.distribution(
            expression=shrinking_ntuple.expr("$TaNCDalitz%s"%dalitz),
            selection=shrinking_ntuple.expr(selection),
            binning = (100, 0, 4),
            x_axis_title = "Dalitz%s"%dalitz,
            y_min = 1e-2, logy=True
            )
        
        # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
        result["Dalitz%s"%dalitz]['legend'].make_legend().Draw()
    
        for format in formats:
            canvas.SaveAs("%s/shrinkingCone_Dalitz%s.%s"%(subdir,dalitz,format))


    return result

def makeOutlierPlots(canvas, selection = "abs($jetEta) < 2.5 & $jetPt > 5", formats = ["png","pdf"], subdir="plots"):
    result={}
    for charge in ["","Charged","Neutral"]:
        result.update( makePtPlots(canvas, charge,"Outlier", 4, selection, formats, subdir) )
        result.update( makeAnglePlots(canvas, charge,"Outlier", 4, selection, formats, subdir) )

        # Compare basic distributions
        result["%sOutlierSumPt"%charge] = plotter.distribution(
            expression=shrinking_ntuple.expr('$TaNC%sOutlierSumPt'%(charge)),
            selection=shrinking_ntuple.expr(selection),
            binning = (100, 0, 70),
            x_axis_title = "Sum %s Outlier P_{T} [GeV/c]"%(charge),
            y_min = 1e-2, logy=True
            )
        
        # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
        result["%sOutlierSumPt"%charge]['legend'].make_legend().Draw()

        for format in formats:
            canvas.SaveAs("%s/shrinkingCone_%sOutlierSumPt.%s"%(subdir,charge,format))

    for charge in ["","Charged"]:
        # Compare basic distributions
        result["OutlierN%s"%charge] = plotter.distribution(
            expression=shrinking_ntuple.expr('$TaNCOutlierN%s'%(charge)),
            selection=shrinking_ntuple.expr(selection),
            binning = (101, -0.5, 100.5),
            x_axis_title = "%s Outlier Count"%(charge),
            y_min = 1e-2, logy=True
            )
        
        # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
        result["OutlierN%s"%charge]['legend'].make_legend().Draw()

        for format in formats:
            canvas.SaveAs("%s/shrinkingCone_OutlierN%s.%s"%(subdir,charge,format))
    result["OutlierMass"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCOutlierMass'),
        selection=shrinking_ntuple.expr(selection),
        binning = (100, 0, 70),
        x_axis_title = "Inv. Mass of Outliers [GeV/c]",
        y_min = 1e-2, logy=True
        )
    
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    result["OutlierMass"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_OutlierMass.%s"%(subdir,format))
        
    return result

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

    result = {}

    selection = "abs($jetEta) < 2.5 & $jetPt > 5"
    formats = ["png","pdf"]
    subdir="plots/combined"

    result["DecayMode"] = plotter.distribution(
        expression=shrinking_ntuple.expr('$TaNCPt'),
        selection=shrinking_ntuple.expr(selection),
        binning = (21, -0.5, 20.5),
        x_axis_title = "#tau decay mode",
        y_min = 1e-2, logy=True
        )
    result["DecayMode"]['legend'].make_legend().Draw()
    
    for format in formats:
        canvas.SaveAs("%s/shrinkingCone_DecayMode.%s"%(subdir,format))

    result.update(makeKiniticPlots(canvas, selection=selection, formats=formats, subdir=subdir))
    result.update(makeOutlierPlots(canvas, selection=selection, formats=formats, subdir=subdir))
    result.update( makePtPlots(canvas, "","PiZero",2, selection=selection, formats=formats, subdir=subdir) )
    result.update( makePtPlots(canvas, "","Track",2, selection=selection, formats=formats, subdir=subdir) )
    result.update( makeAnglePlots(canvas, "","PiZero",2, selection=selection, formats=formats, subdir=subdir) )
    result.update( makeAnglePlots(canvas, "","Track",2, selection=selection, formats=formats, subdir=subdir) )

    
