#!/usr/bin/env python

#
# This Macro assumes:
# mySamples.py exists and
# data2, ppMuX, minBias, qcd dataset  
# are defined there
#

import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import  TauAnalysis.TauIdEfficiency.ntauples.plotting as plot2
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import mySamples

import os
import sys
import copy

# style Definitions
jets_CMS_PRELIMINARY_UPPER_LEFT = style.CMS_PRELIMINARY_UPPER_LEFT.Clone()
jets_CMS_PRELIMINARY_UPPER_LEFT.SetX2(0.60)
jets_CMS_PRELIMINARY_UPPER_LEFT.SetY1(0.850)
jets_CMS_PRELIMINARY_UPPER_LEFT.SetTextSize(0.055)

jets_LUMI_LABEL_UPPER_LEFT = style.LUMI_LABEL_UPPER_LEFT.Clone()
jets_LUMI_LABEL_UPPER_LEFT.SetX2(0.60)
jets_LUMI_LABEL_UPPER_LEFT.SetY2(0.840)
jets_LUMI_LABEL_UPPER_LEFT.SetY1(0.790)
jets_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_SQRTS_LABEL_UPPER_LEFT = style.SQRTS_LABEL_UPPER_LEFT.Clone()
jets_SQRTS_LABEL_UPPER_LEFT.SetX2(0.60)
jets_SQRTS_LABEL_UPPER_LEFT.SetY2(0.785)
jets_SQRTS_LABEL_UPPER_LEFT.SetY1(0.735)
jets_SQRTS_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_PT_CUT_LABEL_UPPER_LEFT = style.PT_CUT_LABEL_UPPER_LEFT.Clone()
jets_PT_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
jets_PT_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
jets_PT_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
jets_PT_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_ETA_CUT_LABEL_UPPER_LEFT = style.ETA_CUT_LABEL_UPPER_LEFT.Clone()
jets_ETA_CUT_LABEL_UPPER_LEFT.SetX2(0.60)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetY2(0.740)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetY1(0.690)
jets_ETA_CUT_LABEL_UPPER_LEFT.SetTextSize(0.045)

PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = ROOT.TPaveText(0.125, 0.665, 0.45, 0.755, "NDC")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("E_{T} > 10 GeV,\n")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.AddText("|#eta| < 2.5")
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextAlign(13)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.035)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetFillStyle(0)
PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetBorderSize(0)

jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT = PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.Clone()
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetX2(0.60)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY2(0.740)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetY1(0.640)
jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT.SetTextSize(0.045)

jets_DEFAULT_LABELS = [
    jets_CMS_PRELIMINARY_UPPER_LEFT,
#    jets_LUMI_LABEL_UPPER_LEFT,
    jets_SQRTS_LABEL_UPPER_LEFT
]

def plot_ratio_deviations(h_data, h_mc):
	output = {}
	difference_histo = h_mc.Clone(h_mc.GetName()+"_diff")
	#difference_histo.Reset()
	for bin in range(difference_histo.GetNbinsX()+1):
		probe_content = h_mc.GetBinContent(bin)
		base_content = h_data.GetBinContent(bin)
		probe_error = h_mc.GetBinError(bin)
		base_error = h_data.GetBinError(bin)
		if probe_content > 0:
			diff = (base_content - probe_content)/probe_content
			base_error_norm = 0
			if base_content > 0:
				base_error_norm = base_error/base_content
			err = (base_content/probe_content)*ROOT.TMath.sqrt(
		        	(probe_error/probe_content)**2 + base_error_norm**2)
                	difference_histo.SetBinContent(bin, diff)
                	difference_histo.SetBinError(bin, err)

	difference_histo.GetXaxis().SetTitleOffset(1.5)
    	difference_histo.GetXaxis().SetTitleSize(0.08)
    	difference_histo.GetXaxis().SetLabelSize(0.08)
    	difference_histo.GetXaxis().SetTickLength(0.055)
    	difference_histo.GetYaxis().SetNdivisions(505)
    	difference_histo.GetYaxis().SetTitleOffset(0.9)
    	difference_histo.GetYaxis().SetTitleSize(0.06)
    	difference_histo.GetYaxis().SetLabelSize(0.08)
    	difference_histo.GetYaxis().SetTickLength(0.04)
    	difference_histo.SetStats(0)
    	difference_histo.SetMaximum(+0.80)
    	difference_histo.SetMinimum(-0.80)
    	#difference_histo.Draw("")
    	output['diff'] = difference_histo
    	return output

def makeJetEffPlot(expression, denominator,numerator, 
		extra_labels,
                binning, x_axis_title, y_min, y_max, logy,
                filename):

    #----------------------------------------------------------------------------
    # Plot distribution of Data vs. Monte Carlo predictions
    #----------------------------------------------------------------------------

    canvas_diff = ROOT.TCanvas("canvas_diff", "canvas_diff", 500, 625)
    canvas_diff.cd()
    topPad_diff = ROOT.TPad("topPad_diff", "topPad_diff", 0.00, 0.35, 1.00, 1.00)
    topPad_diff.SetFillColor(10)
    topPad_diff.SetLogy(logy)
    topPad_diff.SetTopMargin(0.05)
    topPad_diff.SetLeftMargin(0.15)
    topPad_diff.SetBottomMargin(0.00)
    topPad_diff.SetRightMargin(0.05)
    topPad_diff.SetGridy()
    topPad_diff.Draw()
    topPad_diff.cd()

    eff_result = plotter.efficiency(
        expression = expression,
        denominator = denominator,
        numerator = numerator,
        binning = binning,
	maxNumEntries = 1000000000, verbose = True, x_error_bars = True, 
        extra_labels = extra_labels,
        x_axis_title = x_axis_title,
        y_axis_title = "Efficiency",
        y_min = y_min, y_max = y_max, logy = logy
    )
    # Adjust offset of y-axis title
    eff_result['background'].GetYaxis().SetNdivisions(505)
    eff_result['background'].GetYaxis().SetTickLength(0.04)
    eff_result['background'].GetYaxis().SetTitleOffset(1.2)
    eff_result['background'].GetYaxis().SetTitleSize(0.05) 
    eff_result['background'].GetYaxis().SetLabelSize(0.05)
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coordinates to make_legend(...)
    legend_x1 = 0.58
    legend_y1 = 0.69
    legend_x2 = 0.88
    legend_y2 = 0.88
  
    # Draw the legend - you can pass NDC xl, yl, xh, yh coords to make_legend(...)
    legend_adjusted_coordinates = style.adjust_coordinates(topPad_diff, legend_x1, legend_y1, legend_x2, legend_y2)
    eff_result['legend'].make_legend(
        legend_adjusted_coordinates['x1'], legend_adjusted_coordinates['y1'],
        legend_adjusted_coordinates['x2'], legend_adjusted_coordinates['y2']).Draw()

    canvas_diff.cd()
    canvas_diff.Update()
    
    bottomPad_diff = ROOT.TPad("bottomPad_diff", "bottomPad_diff", 0.00, 0.00, 1.00, 0.35)
    bottomPad_diff.SetFillColor(10)
    bottomPad_diff.SetLogy(False)
    bottomPad_diff.SetTopMargin(0.00)
    bottomPad_diff.SetLeftMargin(0.15)
    bottomPad_diff.SetBottomMargin(0.30)
    bottomPad_diff.SetRightMargin(0.05)
    bottomPad_diff.SetGridy()
    bottomPad_diff.Draw()
    bottomPad_diff.cd()

    h_data = plot2.draw(list(plotter.samples['data2']['sample'].events_and_weights(None)),expression, numerator, "data2Num", binning=binning ) # muon stream data
    h_data.Divide(plot2.draw(list(plotter.samples['data2']['sample'].events_and_weights(None)),expression, denominator, "data2Den", binning=binning ) )
    h_mc = plot2.draw(list(plotter.samples['ppMuX']['sample'].events_and_weights(None)),expression, numerator, "ppMuXNum", binning=binning )
    h_mc.Divide(plot2.draw(list(plotter.samples['ppMuX']['sample'].events_and_weights(None)),expression, denominator, "ppMuXDen", binning=binning ) )
    h_mc0 = plot2.draw(list(plotter.samples['minBias']['sample'].events_and_weights(None)),expression, numerator, "minBiasNum", binning=binning )
    h_mc0.Divide(plot2.draw(list(plotter.samples['minBias']['sample'].events_and_weights(None)),expression, denominator, "minBiasDen", binning=binning ) )
    h_mc2 = plot2.draw(list(plotter.samples['qcd']['sample'].events_and_weights(None)),expression, numerator, "qcdNum", binning=binning )
    h_mc2.Divide(plot2.draw(list(plotter.samples['qcd']['sample'].events_and_weights(None)),expression, denominator, "qcdDen", binning=binning ) )

# muon stream and ppmux
    diff_ppMuX = plot_ratio_deviations(h_data, h_mc)
    style.update_histo_style(diff_ppMuX['diff'], style.QCD_MC_PYTHIA8_STYLE_HIST)
    diff_ppMuX['diff'].SetTitle(";%s;#frac{#Delta(Data,MC)  }{MC};"%(x_axis_title))
    diff_ppMuX['diff'].Draw("")

    diff_minBias = plot_ratio_deviations(h_mc0, h_mc2)
    style.update_histo_style(diff_minBias['diff'], style.QCD_MC_PYTHIA6_STYLE_HIST)
    diff_minBias['diff'].Draw("same")

    leg2 = ROOT.TLegend(0.8,0.38,0.94,0.5);
    leg2.AddEntry(diff_ppMuX['diff'],"Muon,ppMuX  ","lep");
    leg2.AddEntry(diff_minBias['diff'],"JetMETTau,QCD  ","lep");
    leg2.SetFillStyle(0);
    leg2.Draw();
    canvas_diff.cd()
    canvas_diff.Update()
    canvas_diff.SaveAs("plots/" + filename  + "_dataMC_diff.pdf")

# diffs mc data
    bottomPad_diff.cd()
    diff_data = plot_ratio_deviations(h_mc0, h_data)
    style.update_histo_style(diff_data['diff'], style.DATA_STYLE)
    diff_data['diff'].SetTitle(";%s;#frac{#Delta(JetMETTau,Muon)   }{Muon};"%(x_axis_title))
    diff_data['diff'].Draw("")
    diff_mc = plot_ratio_deviations(h_mc2, h_mc)
    style.update_histo_style(diff_mc['diff'], style.QCD_MC_STYLE_HIST)
    diff_mc['diff'].Draw("same")
#
    leg = ROOT.TLegend(0.8,0.38,0.94,0.5);
    leg.AddEntry(diff_data['diff'],"#Delta(Data)  ","lep");
    leg.AddEntry(diff_mc['diff'],"#Delta(MC)  ","lep");
    leg.SetFillStyle(0);
    leg.Draw();
#
    canvas_diff.cd()
    canvas_diff.Update()
    canvas_diff.SaveAs("plots/" + filename  + "_muon_jetmettau_diff.pdf")

    canvas_diff.IsA().Destructor(canvas_diff)

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    if not os.path.isdir('plots'):
        os.mkdir('plots')

    # Build the plot manager.  The plot manager keeps track of all the samples
    # and ensures they are correctly normalized w.r.t. luminosity.  See 
    # samples.py for available samples.
    plotter = PlotManager()

    # Add each sample we want to plot/compare
    plotter.add_sample(mySamples.data2, "Muon Data(L=1.08pb^{-1})   ", **style.DATA_STYLE)
    plotter.add_sample(mySamples.ppMuX_mc, "PPMuX MC", **style.QCD_MC_PYTHIA6_STYLE_HIST)
    plotter.add_sample(mySamples.minBias_mc, "JetMETTau(Data L=12.8nb^{-1})   ", **style.MINBIAS_MC_STYLE )
    plotter.add_sample(mySamples.qcd, "QCD MC", **style.QCD_MC_PYTHIA8_STYLE_HIST)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(mySamples.data2.effective_luminosity())

    # Build the ntuple manager
    ntuple_manager = mySamples.data2.build_ntuple_manager("tauIdEffNtuple")

    # Get HLT trigger decisions
#    hlt = ntuple_manager.get_ntuple("TriggerResults")
 
    ######################################################
    ####      Plot PFJet distributions                ####
    ######################################################
    #pt_binning_fine = (0, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100)
    pt_binning_fine = (0, 10, 20, 40, 60, 70, 100)
    #pt_binning_fine = (20,0,100)
    eta_binning_fine = (20, -2.5, 2.5)
#phi_binning_fine = (25, -3.14, 3.14)

    pfTau_ntuple = ntuple_manager.get_ntuple(
        "patPFTausCleanedShrinkingCone")
    pfTau_base_selection = pfTau_ntuple.expr('$jetPt > 10 & abs( $jetEta) < 2.5') #& ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') )# & hlt.expr('$hltMu5 == 1')
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPhi'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byIsolationUsingLeadingPion') & pfTau_ntuple.expr('$againstMuon') & pfTau_ntuple.expr('$againstElectron') & ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') ) & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (20,-3.2,3.2),
        x_axis_title = "Jet #phi",
        y_min = 1.1e-5, y_max = 100., logy =True,
        filename = "pfJetPhiEff"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPhi'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byEcalIsolationUsingLeadingPion') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (20,-3.2,3.2),
        x_axis_title = "Jet #phi",
        y_min = 0, y_max = 1.2, logy =False,
        filename = "pfJetPhiEcalEff"
    )
    pfTau_base_selection = pfTau_ntuple.expr('$jetPt > 10 & abs( $jetEta) < 2.5 ')
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadTrackFinding') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = pt_binning_fine,
        x_axis_title = "Jet E_{T} [GeV]",
        y_min = 0.001, y_max = 1.5, logy =False,
        filename = "pfJetEtLeadTrack"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = pt_binning_fine,
        x_axis_title = "Jet E_{T} [GeV]",
        y_min = 0.001, y_max = 1.5, logy =False,
        filename = "pfJetEtLeadPt"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byEcalIsolationUsingLeadingPion') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = pt_binning_fine,
        x_axis_title = "Jet E_{T} [GeV]",
        y_min = 0.001, y_max = 1.2, logy =False,
        filename = "pfJetEtECalIso"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byEcalIsolationUsingLeadingPion')  & pfTau_ntuple.expr('$byTrackIsolationUsingLeadingPion') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = pt_binning_fine,
        x_axis_title = "Jet E_{T} [GeV]",
        y_min = 0.001, y_max = 0.20, logy =False,
        filename = "pfJetEtTrackerIso"
    )

    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPt'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byIsolationUsingLeadingPion') & pfTau_ntuple.expr('$againstMuon') & pfTau_ntuple.expr('$againstElectron') & ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') ) & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = pt_binning_fine,
        x_axis_title = "Jet E_{T} [GeV]",
        y_min = 1.1e-5, y_max = 10., logy =True,
        filename = "pfJetEtEff"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byIsolationUsingLeadingPion') & pfTau_ntuple.expr('$againstMuon') & pfTau_ntuple.expr('$againstElectron') & ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') ) & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = eta_binning_fine,
        x_axis_title = "Jet #eta",
        y_min = 1.1e-5, y_max = 10., logy =True,
        filename = "pfJetEtaEff"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetWidth'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byIsolationUsingLeadingPion') & pfTau_ntuple.expr('$againstMuon') & pfTau_ntuple.expr('$againstElectron') & ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') ) & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (30, 0, 0.3),
        x_axis_title = "Jet Width",
        y_min = 1.1e-5, y_max = 1000., logy =True,
        filename = "pfJetWidthEff"
    )
####

    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadTrackFinding') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = eta_binning_fine,
        x_axis_title = "Jet #eta",
        y_min = 0.001, y_max = 1.5, logy =False,
        filename = "pfJetEtaLeadTrack"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = eta_binning_fine,
        x_axis_title = "Jet #eta",
        y_min = 0.001, y_max = 1.0, logy =False,
        filename = "pfJetEtaLeadPt"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byEcalIsolationUsingLeadingPion') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = eta_binning_fine,
        x_axis_title = "Jet #eta",
        y_min = 0.001, y_max = 1.0, logy =False,
        filename = "pfJetEtaECalIso"
    )
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetEta'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byEcalIsolationUsingLeadingPion')  & pfTau_ntuple.expr('$byTrackIsolationUsingLeadingPion') & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = eta_binning_fine,
        x_axis_title = "Jet #eta",
        y_min = 0.001, y_max = 0.20, logy =False,
        filename = "pfJetEtaTrackerIso"
    )

    pfTau_base_selection = pfTau_base_selection + pfTau_ntuple.expr('$jetPt@.at(0)')
    makeJetEffPlot(
        expression = pfTau_ntuple.expr('$jetPhi'),
	denominator =  pfTau_base_selection,
        numerator = pfTau_ntuple.expr('$byLeadPionPtCut') & pfTau_ntuple.expr('$byIsolationUsingLeadingPion') & pfTau_ntuple.expr('$againstMuon') & pfTau_ntuple.expr('$againstElectron') & ( pfTau_ntuple.expr('$numChargedParticlesSignalCone == 1') | pfTau_ntuple.expr('$numChargedParticlesSignalCone == 3') ) & pfTau_base_selection,
        extra_labels = [ jets_PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT ],
        binning = (20,-3.2,3.2),
        x_axis_title = "Jet #phi",
        y_min = 1.1e-5, y_max = 100., logy =True,
        filename = "pfLeadJetPhiEff"
    )
