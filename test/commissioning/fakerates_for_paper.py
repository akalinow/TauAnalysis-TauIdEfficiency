#!/usr/bin/env python

'''
Plot fake-rates of different tau id algorithms for
 o QCD multi-jet
 o QCD muon enriched
 o W --> mu nu + jets
Data compared to Monte Carlo predictions

Used to make plots for TAU-11-001

Author: Christian Veelken

'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import TauAnalysis.DQMTools.plotterStyleDefinitions_cfi as colors
import copy

# Define of input files.
import samples_cache as samples

# Define sample names
sample_dijet_data  = samples.data_dijet
sample_dijet_mc    = samples.qcddijet_mc
sample_ppmux_data  = samples.data_ppmux
sample_ppmux_mc    = samples.ppmux15_mc
sample_wjets_data  = samples.data_wjets
sample_wjets_mc    = samples.wjets_mc

# Define integrated luminosity of analyzed dataset
intLumiData = 5100 # nb^-1

# Define denominator (phase-space) for which fake-rates are to be determined
denominator_phase_space = "$jetPt > 20.0 & abs($jetEta) < 2.3"

# Define maximum number of Ntuple entries to process
# (NOTE: use values different from 1000000000 for debugging purposes only !!)
maxNumEntries = 1000000000
#maxNumEntries = 100000

# Define marker styles
custom_style_dijet_data = copy.deepcopy(style.DATA_STYLE)
custom_style_dijet_data['marker_color']  = colors.color_darkBlue.value()
custom_style_dijet_data['line_color']    = colors.color_darkBlue.value()
custom_style_dijet_data['marker_style']  = 22 # closed upward-facing triangle
custom_style_dijet_mc = copy.deepcopy(style.QCD_MC_STYLE_DOTS)
custom_style_dijet_mc['marker_color']    = colors.color_darkBlue.value()
custom_style_dijet_mc['line_color']      = colors.color_darkBlue.value()
custom_style_dijet_mc['marker_style']    = 26 # open upward-facing triangle

custom_style_ppmux_data = copy.deepcopy(style.DATA_STYLE)
custom_style_ppmux_data['marker_color']      = colors.color_violett.value()
custom_style_ppmux_data['line_color']        = colors.color_violett.value()
custom_style_ppmux_data['marker_style']      = 21 # closed square
custom_style_ppmux_mc = copy.deepcopy(style.QCD_MC_STYLE_DOTS)
custom_style_ppmux_mc['marker_color']        = colors.color_violett.value()
custom_style_ppmux_mc['line_color']          = colors.color_violett.value()
custom_style_ppmux_mc['marker_style']        = 25 # open square

custom_style_wjets_data = copy.deepcopy(style.DATA_STYLE)
custom_style_wjets_data['marker_color']      = colors.color_red.value()
custom_style_wjets_data['line_color']        = colors.color_red.value()
custom_style_wjets_data['marker_style']      = 20 # closed circle
custom_style_wjets_mc = copy.deepcopy(style.QCD_MC_STYLE_DOTS)
custom_style_wjets_mc['marker_color']        = colors.color_red.value()
custom_style_wjets_mc['line_color']          = colors.color_red.value()
custom_style_wjets_mc['marker_style']        = 24 # open circle

# Adjust position of labels
custom_CMS_PRELIMINARY_UPPER_LEFT = copy.deepcopy(style.CMS_PRELIMINARY_UPPER_LEFT)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetX1(style.CMS_PRELIMINARY_UPPER_LEFT.GetX1() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetX2(style.CMS_PRELIMINARY_UPPER_LEFT.GetX2() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetY1(style.CMS_PRELIMINARY_UPPER_LEFT.GetY1() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetY2(style.CMS_PRELIMINARY_UPPER_LEFT.GetY2() + 0.050)
custom_CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.045)

custom_LUMI_LABEL_UPPER_LEFT = copy.deepcopy(style.LUMI_LABEL_UPPER_LEFT)
custom_LUMI_LABEL_UPPER_LEFT.SetX1(style.LUMI_LABEL_UPPER_LEFT.GetX1() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetX2(style.LUMI_LABEL_UPPER_LEFT.GetX2() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetY1(style.LUMI_LABEL_UPPER_LEFT.GetY1() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.SetY2(style.LUMI_LABEL_UPPER_LEFT.GetY2() + 0.050)
custom_LUMI_LABEL_UPPER_LEFT.Clear()
custom_LUMI_LABEL_UPPER_LEFT.AddText("#sqrt{s} = 7 TeV, L = %.1f pb^{-1}" % (intLumiData/1000))
custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.045)

custom_PT_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_CUT_LABEL_UPPER_LEFT)
custom_PT_CUT_LABEL_UPPER_LEFT.SetX1(style.PT_CUT_LABEL_UPPER_LEFT.GetX1() + 0.100)
custom_PT_CUT_LABEL_UPPER_LEFT.SetX2(style.PT_CUT_LABEL_UPPER_LEFT.GetX2() + 0.100)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY1(style.PT_CUT_LABEL_UPPER_LEFT.GetY1() + 0.100)
custom_PT_CUT_LABEL_UPPER_LEFT.SetY2(style.PT_CUT_LABEL_UPPER_LEFT.GetY2() + 0.100)
custom_PT_CUT_LABEL_UPPER_LEFT.Clear()
custom_PT_CUT_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c")
custom_PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

custom_ETA_CUT_LABEL_UPPER_LEFT = copy.deepcopy(style.ETA_CUT_LABEL_UPPER_LEFT)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetX1(style.ETA_CUT_LABEL_UPPER_LEFT.GetX1()+ 0.100)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetX2(style.ETA_CUT_LABEL_UPPER_LEFT.GetX2() + 0.100)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY1(style.ETA_CUT_LABEL_UPPER_LEFT.GetY1() + 0.100)
custom_ETA_CUT_LABEL_UPPER_LEFT.SetY2(style.ETA_CUT_LABEL_UPPER_LEFT.GetY2() + 0.100)
custom_ETA_CUT_LABEL_UPPER_LEFT.Clear()
custom_ETA_CUT_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")
custom_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = copy.deepcopy(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX1(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetX1() + 0.100)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX2(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetX2() + 0.100)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY1(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetY1() + 0.100)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY2(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.GetY2() + 0.100)
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.Clear()
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("P_{T} > 20 GeV/c,\n")
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta| < 2.3")
custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.045)

# Define white background "patch" to clear grid lines in area used for labels and legend
white_patch = ROOT.TPave(0.10, 0.65, 0.8975, 0.8950, 0, "NDC")
white_patch.SetFillColor(10)

# Define labels for different algorithms
algo_label_template = ROOT.TPaveText(0.50, 0.615, 0.88, 0.655, "NDC")
algo_label_template.SetTextAlign(33)
algo_label_template.SetTextSize(0.035)
algo_label_template.SetTextColor(2)
algo_label_template.SetFillStyle(0)
algo_label_template.SetBorderSize(0)

algo_label_fixed = copy.deepcopy(algo_label_template)
algo_label_fixed.AddText("Fixed Cone")

algo_label_shrinking = copy.deepcopy(algo_label_template)
algo_label_shrinking.AddText("Shrinking Cone")

algo_label_tanc_loose = copy.deepcopy(algo_label_template)
algo_label_tanc_loose.AddText("TaNC loose")
algo_label_tanc_medium = copy.deepcopy(algo_label_template)
algo_label_tanc_medium.AddText("TaNC medium")
algo_label_tanc_tight = copy.deepcopy(algo_label_template)
algo_label_tanc_tight.AddText("TaNC tight")

algo_label_hps_loose = copy.deepcopy(algo_label_template)
algo_label_hps_loose.AddText("HPS loose")
algo_label_hps_medium = copy.deepcopy(algo_label_template)
algo_label_hps_medium.AddText("HPS medium")
algo_label_hps_tight = copy.deepcopy(algo_label_template)
algo_label_hps_tight.AddText("HPS tight")

algo_label_calo = copy.deepcopy(algo_label_template)
algo_label_calo.AddText("TCTau")

def makeFakeratePlots(algorithm, numerator,
                      denominator_dijet, plotter_dijet,
                      denominator_ppmux, plotter_ppmux,
                      denominator_wjets, plotter_wjets,
                      expression, binning,
                      x_axis_title, y_min = 1e-3, y_max = 0.9, extra_labels = [], extra_labels_diff = [],
                      maxNumEntries = 1000000000,
                      outputFileName = "fakerates_for_paper.pdf"):

    #print("<makeFakeratePlots>:")

    result_algorithm = {}
    
    #----------------------------------------------------------------------------
    # Make "regular" fake-rate plots
    #----------------------------------------------------------------------------

    canvas = ROOT.TCanvas("canvas_diff", "canvas_diff", 600, 600)
    canvas.cd()

    # Show difference in fake-rate between QCD di-jet, QCD muon enriched and W + jets events
    eff_dijet = plotter_dijet.efficiency(
        expression,
        denominator_dijet & numerator,
        denominator_dijet,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Fake Rate",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['dijet'] = eff_dijet

    eff_ppmux = plotter_ppmux.efficiency(
        expression,
        denominator_ppmux & numerator,
        denominator_ppmux,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Fake Rate",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['ppmux'] = eff_ppmux

    eff_wjets = plotter_wjets.efficiency(
        expression,
        denominator_wjets & numerator,
        denominator_wjets,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Fake Rate",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['wjets'] = eff_wjets

    # Draw fake-rate plots
    eff_dijet['samples'][sample_dijet_data.name].Draw('epsame')
    eff_dijet['samples'][sample_dijet_mc.name].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_data.name].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_mc.name].Draw('epsame')
    eff_wjets['samples'][sample_wjets_data.name].Draw('epsame')
    eff_wjets['samples'][sample_wjets_mc.name].Draw('epsame')

    # Clear grid lines in area used for labels and legend
    #white_patch.Draw()
    #eff_dijet['samples'][sample_dijet_data_wPU.name].Draw('axissame')

    # Draw legend
    legend = ROOT.TLegend(0.45, 0.61, 0.88, 0.88, "", "brNDC")
    legend.SetTextSize(0.03)
    legend.SetFillColor(0)
    legend.SetLineColor(1)
    legend.SetBorderSize(1)
    legend.AddEntry(eff_wjets['samples'][sample_wjets_data.name], "W #rightarrow #mu #nu Data",       "p")
    legend.AddEntry(eff_wjets['samples'][sample_wjets_mc.name],   "W #rightarrow #mu #nu Simulation", "p")
    legend.AddEntry(eff_dijet['samples'][sample_dijet_data.name], "QCDj Data",                        "p")
    legend.AddEntry(eff_dijet['samples'][sample_dijet_mc.name],   "QCDj Simulation",                  "p")
    legend.AddEntry(eff_ppmux['samples'][sample_ppmux_data.name], "QCD#mu Data",                      "p")
    legend.AddEntry(eff_ppmux['samples'][sample_ppmux_mc.name],   "QCD#mu Simulation",                "p")
    legend.Draw()

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
    style.LUMI_LABEL_UPPER_LEFT.Draw()
    style.SQRTS_LABEL_UPPER_LEFT.Draw()

    for extra_label in extra_labels:
        extra_label.Draw()

    canvas.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.png"  % algorithm)
    canvas.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.pdf"  % algorithm)
    canvas.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.root" % algorithm)

    #----------------------------------------------------------------------------
    # Make fake-rate plots with normalized differences added below plots
    #----------------------------------------------------------------------------
    
    canvas_diff = ROOT.TCanvas("canvas_diff", "canvas_diff", 500, 625)
    canvas_diff.cd()
    
    topPad_diff = ROOT.TPad("topPad_diff", "topPad_diff", 0.00, 0.35, 1.00, 1.00)
    topPad_diff.SetFillColor(10)
    topPad_diff.SetLogy(True)
    topPad_diff.SetTopMargin(0.05)
    topPad_diff.SetLeftMargin(0.18)
    topPad_diff.SetBottomMargin(0.00)
    topPad_diff.SetRightMargin(0.05)
    #topPad_diff.SetGridy()
    topPad_diff.Draw()
    topPad_diff.cd()

    # Adjust offset of y-axis title
    eff_dijet['background'].GetYaxis().SetNdivisions(505)
    eff_dijet['background'].GetYaxis().SetTickLength(0.04)
    eff_dijet['background'].GetYaxis().SetTitleOffset(1.2)
    eff_dijet['background'].GetYaxis().SetTitleSize(0.05) 
    eff_dijet['background'].GetYaxis().SetLabelSize(0.05)
    eff_dijet['background'].Draw()

    # Draw fake-rate plots
    eff_dijet['samples'][sample_dijet_data.name].Draw('epsame')
    eff_dijet['samples'][sample_dijet_mc.name].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_data.name].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_mc.name].Draw('epsame')
    eff_wjets['samples'][sample_wjets_data.name].Draw('epsame')
    eff_wjets['samples'][sample_wjets_mc.name].Draw('epsame')

    # Clear grid lines in area used for labels and legend
    #white_patch.Draw()
    #eff_dijet['background'].Draw('axissame')
    
    # Draw legend
    legend.SetX1NDC(legend.GetX1NDC() + 0.050)
    legend.SetX2NDC(legend.GetX2NDC() + 0.050)
    legend.SetY1NDC(legend.GetY1NDC() + 0.050)
    legend.SetY2NDC(legend.GetY2NDC() + 0.050)
    legend.Draw()

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    custom_CMS_PRELIMINARY_UPPER_LEFT.Draw()
    custom_LUMI_LABEL_UPPER_LEFT.Draw()
    
    for extra_label_diff in extra_labels_diff:
        extra_label_diff.Draw()

    canvas_diff.cd()
    
    bottomPad_diff = ROOT.TPad("bottomPad_diff", "bottomPad_diff", 0.00, 0.00, 1.00, 0.35)
    bottomPad_diff.SetFillColor(10)
    bottomPad_diff.SetLogy(False)
    bottomPad_diff.SetTopMargin(0.00)
    bottomPad_diff.SetLeftMargin(0.18)
    bottomPad_diff.SetBottomMargin(0.20)
    bottomPad_diff.SetRightMargin(0.05)
    #bottomPad_diff.SetGridy()
    bottomPad_diff.Draw()
    bottomPad_diff.cd()

    # Show normalized difference (Data - MC)/MC
    diff_dijet = plotter_dijet.plot_eff_deviations(
        eff_dijet,
        sample_dijet_mc.name,
        [ sample_dijet_data.name ],
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.34, y_max = +0.34, logy = False
    )
    result_algorithm['dijet_DataToMC_diff'] = diff_dijet
    
    diff_ppmux = plotter_ppmux.plot_eff_deviations(
        eff_ppmux,
        sample_ppmux_mc.name,
        [ sample_ppmux_data.name ],
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.34, y_max = +0.34, logy = False
    )
    result_algorithm['ppmux_DataToMC_diff'] = diff_ppmux
    
    diff_wjets = plotter_wjets.plot_eff_deviations(
        eff_wjets,
        sample_wjets_mc.name,
        [ sample_wjets_data.name ],
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.34, y_max = +0.34, logy = False
    )
    result_algorithm['wjets_DataToMC_diff'] = diff_wjets

    diff_dijet['background'].GetYaxis().SetNdivisions(505)
    diff_dijet['background'].GetXaxis().SetTitleOffset(1.1)
    diff_dijet['background'].GetXaxis().SetTitleSize(0.08)
    diff_dijet['background'].GetXaxis().SetLabelSize(0.08)
    diff_dijet['background'].GetXaxis().SetTickLength(0.055)
    diff_dijet['background'].GetYaxis().CenterTitle()
    diff_dijet['background'].GetYaxis().SetTitleOffset(0.9)
    diff_dijet['background'].GetYaxis().SetTitleSize(0.08)
    diff_dijet['background'].GetYaxis().SetLabelSize(0.08)
    diff_dijet['background'].GetYaxis().SetTickLength(0.04)
    diff_dijet['background'].Draw()

    diff_dijet['samples'][sample_dijet_data.name].Draw("epsame")
    diff_ppmux['samples'][sample_ppmux_data.name].Draw("epsame")
    diff_wjets['samples'][sample_wjets_data.name].Draw("epsame")

    canvas_diff.cd()
    canvas_diff.Update()

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.png"  % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.pdf"  % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.root" % algorithm)

    return result_algorithm
       
if __name__ == "__main__":

    # define PlotManagers for QCD multi-jet, QCD muon enriched and W + jets samples
    plotter_dijet = PlotManager()
    plotter_dijet.add_sample(sample_dijet_data, "QCDj Data",                        **custom_style_dijet_data)
    plotter_dijet.add_sample(sample_dijet_mc,   "QCDj Simulation",                  **custom_style_dijet_mc)
    plotter_dijet.set_integrated_lumi(intLumiData)
    
    plotter_ppmux = PlotManager()
    plotter_ppmux.add_sample(sample_ppmux_data, "QCD#mu Data",                      **custom_style_ppmux_data)
    plotter_ppmux.add_sample(sample_ppmux_mc,   "QCD#mu Simulation",                **custom_style_ppmux_mc)
    plotter_ppmux.set_integrated_lumi(intLumiData)
        
    plotter_wjets = PlotManager()
    plotter_wjets.add_sample(sample_wjets_data, "W #rightarrow #mu #nu Data",       **custom_style_wjets_data)
    plotter_wjets.add_sample(sample_wjets_mc,   "W #rightarrow #mu #nu Simulation", **custom_style_wjets_mc)
    plotter_wjets.set_integrated_lumi(intLumiData)
    
    # Build the ntuple manager
    ntuple_manager = sample_wjets_data.build_ntuple_manager("tauIdEffNtuple")

    nTuples = {
        "shrinkingCone" : ntuple_manager.get_ntuple("shrinking"),
        "fixedCone"     : ntuple_manager.get_ntuple("fixed"),
        "TaNC"          : ntuple_manager.get_ntuple("hpstanc"),
        "hps"           : ntuple_manager.get_ntuple("hps"),
        "calo"          : ntuple_manager.get_ntuple("calo")
    }
    
    hlt = ntuple_manager.get_ntuple("patTriggerEvent")

    # Define binning options
    binning_pt  = (0, 10,  15,  20,  25, 30, 40, 60, 100, 200)
    binning_eta = (15, -2.5,  +2.5)
    binning_phi = (15, -3.14, +3.14)

    # Define the numerators to plot
    pfString    = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTrackIsolation > 0.5 " \
                 + " & $byEcalIsolation > 0.5 & ($numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3)"
    caloString  = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byIsolation > 0.5 " \
                 + " & $etSumIsolationECAL < 5 & ($numSignalTracks ==1 || $numSignalTracks ==3)"
    
    numerators = {
        "shrinkingCone" : {
            'expr'       : nTuples["shrinkingCone"].expr(pfString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        },
        "fixedCone" : {
            'expr'       : nTuples["fixedCone"].expr(pfString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        },
        "TaNC_loose" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCloose > 0.5 & $pt > 15.0'),
            'style_name' : "byTaNCloose",
            'nice_name'  : "TaNC loose"
        },
        "TaNC_medium" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCmedium > 0.5 & $pt > 15.0'),
            'style_name' :"byTaNCmedium",
            'nice_name'  : "TaNC medium"
        },
        "TaNC_tight" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCtight > 0.5 & $pt > 15.0'),
            'style_name' : "byTaNCtight",
            'nice_name'  : "TaNC tight"
        },
        "hps_loose" : {
            'expr'       : nTuples["hps"].expr('$bgDecayModeFinding > 0.5 & $byIsolationLoose > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationLoose",
            'nice_name'  : "HPS loose"
        },
        "hps_medium" : {
            'expr'       : nTuples["hps"].expr('$bgDecayModeFinding > 0.5 & $byIsolationMedium > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationMedium",
            'nice_name'  : "HPS medium"
        },
        "hps_tight" : {
            'expr'       : nTuples["hps"].expr('$bgDecayModeFinding > 0.5 & $byIsolationTight > 0.5 & $pt > 15.0'),
            'style_name' : "byIsolationTight",
            'nice_name'  : "HPS tight"
        },
        "calo" : {
            'expr'       : nTuples["calo"].expr(caloString),
            'style_name' : "OneOrThreeProng",
            'nice_name'  : ""
        }
    }

    extra_labels = {}
    extra_labels['fixedCone']     = [ algo_label_fixed       ]
    extra_labels['shrinkingCone'] = [ algo_label_shrinking   ]
    extra_labels['TaNC_loose']    = [ algo_label_tanc_loose  ]
    extra_labels['TaNC_medium']   = [ algo_label_tanc_medium ]
    extra_labels['TaNC_tight']    = [ algo_label_tanc_tight  ]
    extra_labels['hps_loose']     = [ algo_label_hps_loose   ]
    extra_labels['hps_medium']    = [ algo_label_hps_medium  ]
    extra_labels['hps_tight']     = [ algo_label_hps_tight   ]
    extra_labels['calo']          = [ algo_label_calo        ]
        
    ##for algorithm_discriminator in [ 'TaNC_loose', 'TaNC_medium', 'TaNC_tight',
    ##                                 'hps_loose',  'hps_medium',  'hps_tight'  ]:
    for algorithm_discriminator in [ 'TaNC_loose', 'hps_loose' ] :

        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue

        algorithm     = algorithm_discriminator
        discriminator = ""
        if algorithm_discriminator.find("_") != -1:
            algorithm     = algorithm_discriminator[:algorithm_discriminator.find("_")]
            discriminator = algorithm_discriminator[algorithm_discriminator.find("_") + 1:]
        print("algorithm = %s, discriminator = %s" % (algorithm, discriminator))
        
        # Define the denominators
        denominator_jetId     = nTuples[algorithm].expr("$jetIdLoose > 0.5")
        ##denominator_jetId     = nTuples[algorithm].expr("$jetPt > 20.0")
        denominator_e_mu_veto = nTuples[algorithm].expr("$againstElectronLoose < 0.6 & $againstMuonTight > 0.5")
        denominator_dijet     = nTuples[algorithm].expr(denominator_phase_space) & denominator_jetId & denominator_e_mu_veto 
        ##                     & hlt.expr('$hltJet30v1bit > 0.5') & nTuples[algorithm].expr("$probeJet30v1 > 0.5")
        denominator_ppmux     = nTuples[algorithm].expr(denominator_phase_space) & denominator_jetId & denominator_e_mu_veto 
        ##                     & hlt.expr('$hltMu15v2bit > 0.5') 
        denominator_wjets     = nTuples[algorithm].expr(denominator_phase_space) & denominator_jetId & denominator_e_mu_veto
        ##                     & hlt.expr('$hltMu15v2bit > 0.5')

        extra_labels_pt  = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_pt.append(style.ETA_CUT_LABEL_UPPER_LEFT)
        extra_labels_diff_pt  = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_pt.append(custom_ETA_CUT_LABEL_UPPER_LEFT)
        
        makeFakeratePlots(algorithm_discriminator, numerators[algorithm_discriminator]['expr'],
                          denominator_dijet, plotter_dijet, 
                          denominator_ppmux, plotter_ppmux, 
                          denominator_wjets, plotter_wjets, 
                          nTuples[algorithm].expr('$jetPt'), binning_pt,
                          'Jet P_{T} [GeV/c]', 6.5e-4, 0.59, extra_labels_pt, extra_labels_diff_pt,
                          maxNumEntries,
                          "plots/fakerates_for_paper_jetPt.pdf")

        extra_labels_eta = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_eta.append(style.PT_CUT_LABEL_UPPER_LEFT)
        extra_labels_diff_eta = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_eta.append(custom_PT_CUT_LABEL_UPPER_LEFT)
        
        ##makeFakeratePlots(algorithm_discriminator, numerator[algorithm]['expr'],
        ##                  denominator_dijet, plotter_dijet, 
        ##                  denominator_ppmux, plotter_ppmux, 
        ##                  denominator_wjets, plotter_wjets,
        ##                  nTuples[algorithm].expr('$jetEta'), binning_eta,
        ##                  'Jet #eta', 6.5e-4, 0.59, extra_labels_eta, extra_labels_diff_eta,
        ##                  maxNumEntries,
        ##                  "plots/fakerates_for_paper_jetEta.pdf")

        extra_labels_phi = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_phi.append(style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
        extra_labels_diff_phi = copy.deepcopy(extra_labels[algorithm_discriminator])
        extra_labels_diff_phi.append(custom_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT)
        
        ##makeFakeratePlots(algorithm_discriminator, numerator[algorithm]['expr'],
        ##                  denominator_dijet, plotter_dijet, 
        ##                  denominator_ppmux, plotter_ppmux, 
        ##                  denominator_wjets, plotter_wjets,
        ##                  nTuples[algorithm].expr('$jetPhi'), binning_phi,
        ##                  'Jet #phi', 6.5e-4, 0.59, extra_labels_phi, extra_labels_diff_phi,
        ##                  maxNumEntries,
        ##                  "plots/fakerates_for_paper_jetPhi.pdf")







