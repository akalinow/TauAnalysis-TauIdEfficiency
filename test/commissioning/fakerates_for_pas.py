#!/usr/bin/env python

'''
Plot fake-rates of different tau id algorithms for Data compared to QCD multi-jet Monte Carlo predictions

Authors: Matthias Edelhoff, Christian Veelken

'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import copy

# Definition of input files.
import samples_cache as samples

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

algo_label_tanc = copy.deepcopy(algo_label_template)
algo_label_tanc.AddText("TaNC")

algo_label_hps = copy.deepcopy(algo_label_template)
algo_label_hps.AddText("HPS")

algo_label_calo = copy.deepcopy(algo_label_template)
algo_label_calo.AddText("TCTau")

def makeFakeratePlots(algorithm, y_min = 1e-4, y_max = 9.9, extra_labels = []):
    
    denominator_section = "$probe > 0.5 & $jetPt > 10.0 & abs($jetEta) < 2.5"
    
    denominator = hlt.expr('$hltJet15U > 0.5') & nTuples[algorithm].expr(denominator_section)
    
    for eff_info in numerators[algorithm]: 
        eff_info['expr'] = denominator & eff_info['expr']

    result_algorithm = {}    
        
    pt_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetPt'), 
        denominator,
        numerators[algorithm], 
        binning = pt_binning_fine, 
        y_min = y_min,
        y_max = y_max,
        x_axis_title = 'Jet P_{T} [GeV/c]',
        labels = [
            style.CMS_PRELIMINARY_UPPER_LEFT,
            style.LUMI_LABEL_UPPER_LEFT,
            style.SQRTS_LABEL_UPPER_LEFT,
            style.ETA_CUT_LABEL_UPPER_LEFT
        ],
        extra_labels = extra_labels,
        logy = True
    )
    result_algorithm['Pt'] = pt_effs
    result_algorithm['Pt']['extra_labels'] = [ style.ETA_CUT_LABEL_UPPER_LEFT, ]
    
    pt_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_pt_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_pt_for_pas.pdf"%(algorithm))
    
    eta_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetEta'), 
        denominator,
        numerators[algorithm], 
        binning = eta_binning_fine, 
        y_min = y_min,
        y_max = y_max,
        x_axis_title = 'Jet #eta',
        labels = [
            style.CMS_PRELIMINARY_UPPER_LEFT,
            style.LUMI_LABEL_UPPER_LEFT,
            style.SQRTS_LABEL_UPPER_LEFT,
            style.PT_CUT_LABEL_UPPER_LEFT
        ],
        extra_labels = extra_labels,
        logy = True
    )
    result_algorithm['eta'] = eta_effs
    result_algorithm['eta']['extra_labels'] = [ style.PT_CUT_LABEL_UPPER_LEFT, ]
    
    eta_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_eta_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_eta_for_pas.pdf"%(algorithm))
    
    phi_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetPhi'), 
        denominator,
        numerators[algorithm], 
        binning = phi_binning_fine, 
        y_min = y_min,
        y_max = y_max,
        x_axis_title = 'Jet #phi',
        labels = [
            style.CMS_PRELIMINARY_UPPER_LEFT,
            style.LUMI_LABEL_UPPER_LEFT,
            style.SQRTS_LABEL_UPPER_LEFT,
            style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT
        ],
        extra_labels = extra_labels,
        logy = True)
    result_algorithm['phi'] = phi_effs
    result_algorithm['phi']['extra_labels'] = [ style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT, ]
    
    phi_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_phi_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_phi_for_pas.pdf"%(algorithm))

    return result_algorithm

def makeFakerateComparisonPlots(numerators):

    # Make comparison plot of fake-rate after all tau id. criteria are applied
    # for different algorithms
    for x_var in [ "Pt", "eta", "phi" ]:

        canvas.Clear()

        # build legend
        legend = ROOT.TLegend(0.45, 0.68, 0.88, 0.88, "","brNDC")
        legend.SetTextSize(0.03)
        legend.SetFillColor(0);
        legend.SetLineColor(1);
        legend.SetBorderSize(1);

        fakerate_results['fixedCone'][x_var]['background'].Draw("")

        if numerators['fixedCone'] != "":
            numerator_fixed = numerators['fixedCone']
            fakerate_fixed = fakerate_results['fixedCone'][x_var]['numerators'][numerator_fixed]['data']
            fakerate_fixed.SetMarkerStyle(20)
            fakerate_fixed.SetMarkerColor(ROOT.EColor.kGreen + 2)
            fakerate_fixed.SetLineColor(ROOT.EColor.kGreen + 2)
            fakerate_fixed.Draw("e1p, same")
            legend.AddEntry(fakerate_fixed, "Fixed signal cone", "P")

        if numerators['shrinkingCone'] != "":
            numerator_shrinking = numerators['shrinkingCone']
            fakerate_shrinking = fakerate_results['shrinkingCone'][x_var]['numerators'][numerator_shrinking]['data']
            fakerate_shrinking.SetMarkerStyle(21)
            fakerate_shrinking.SetMarkerColor(ROOT.EColor.kRed)
            fakerate_shrinking.SetLineColor(ROOT.EColor.kRed)
            fakerate_shrinking.Draw("e1p,same")
            legend.AddEntry(fakerate_shrinking, "Shrinking signal cone", "P")

        if numerators['TaNC'] != "":
            numerator_tanc = numerators['TaNC']
            fakerate_tanc = fakerate_results['TaNC'][x_var]['numerators']['byTaNCfrHalfPercent']['data']
            fakerate_tanc.SetMarkerStyle(22)
            fakerate_tanc.SetMarkerColor(ROOT.EColor.kBlue)
            fakerate_tanc.SetLineColor(ROOT.EColor.kBlue)
            fakerate_tanc.Draw("e1p,same")
            legend.AddEntry(fakerate_tanc, "TaNC 0.50%", "P")

        fakerate_hps = fakerate_results['hps'][x_var]['numerators']['byIsolationMedium']['data']
        fakerate_hps.SetMarkerStyle(23)
        fakerate_hps.SetMarkerColor(28)
        fakerate_hps.SetLineColor(28)
        fakerate_hps.Draw("e1p,same")
        legend.AddEntry(fakerate_hps, "HPS medium isolation", "P")

        fakerate_tctau = fakerate_results['calo'][x_var]['numerators']['OneOrThreeProng']['data']
        fakerate_tctau.SetMarkerStyle(29)
        fakerate_tctau.SetMarkerColor(ROOT.EColor.kBlack)
        fakerate_tctau.SetLineColor(ROOT.EColor.kBlack)
        fakerate_tctau.Draw("e1p,same")
        legend.AddEntry(fakerate_tctau, "TCTau", "P")

        # Draw legend
        legend.Draw()
            
        # Draw the preliminary label
        style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
        style.LUMI_LABEL_UPPER_LEFT.Draw()
        style.SQRTS_LABEL_UPPER_LEFT.Draw()
        for extra_label in fakerate_results['calo'][x_var]['extra_labels']:
            extra_label.Draw()
        
        # Save the plot
        canvas.SaveAs("plots/%s.png" % '_'.join(['fakerate_algo_comparison', 'vs', x_var]))
        canvas.SaveAs("plots/%s.pdf" % '_'.join(['fakerate_algo_comparison', 'vs', x_var]))
                    
if __name__ == "__main__":

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

    nTuples = {
        "shrinkingCone": ntuple_manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone"),
        "fixedCone": ntuple_manager.get_ntuple("patPFTausDijetTagAndProbeFixedCone"),
        "TaNC": ntuple_manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone"),
        "hps": ntuple_manager.get_ntuple("patPFTausDijetTagAndProbeHPS"),
        "calo": ntuple_manager.get_ntuple("patCaloTausDijetTagAndProbe")
    }
    
    hlt = ntuple_manager.get_ntuple("TriggerResults")
    
    # Make some plots
    canvas = ROOT.TCanvas("pas", "pas", 500, 500)
    canvas.cd()
    
    pt_binning_fine = (0, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100)
    eta_binning_fine = (25, -2.5, 2.5)
    phi_binning_fine = (25, -3.14, 3.14)

    # Define the numerators to plot
    lead_pion_selection = '$byLeadPionPtCut > 0.5 '
    pfString = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTrackIsolation > 0.5 "
    pfString += " & $byEcalIsolation > 0.5 & ($numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3)"
    caloString = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byIsolation > 0.5 "
    caloString += " & $etSumIsolationECAL < 5 & ($numSignalTracks ==1 || $numSignalTracks ==3)"
    
    numerators = {
        "shrinkingCone":[
            {
                'expr': nTuples["shrinkingCone"].expr(pfString),
                'style_name':"OneOrThreeProng",
                'nice_name': ""
            }
        ],
        "fixedCone":[
            {
                'expr': nTuples["fixedCone"].expr(pfString),
                'style_name':"OneOrThreeProng",
                'nice_name': ""
            }
        ],
        "TaNC":[
            {
                'expr': nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrOnePercent > 0.5'),
                'style_name':"byTaNCfrOnePercent",
                'nice_name': "TaNC 1.00%"
            },
            {
                'expr': nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrHalfPercent > 0.5'),
                'style_name':"byTaNCfrHalfPercent",
                'nice_name': "TaNC 0.50%"
            },
            {
                'expr': nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrQuarterPercent > 0.5'),
                'style_name':"byTaNCfrQuarterPercent",
                'nice_name': "TaNC 0.25%"
            }
         ],
         "hps":[
            {
                'expr': nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationLoose > 0.5'),
                'style_name':"byIsolationLoose",
                'nice_name': "Loose Isolation"
            },
            {
                'expr': nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationMedium > 0.5'),
                'style_name':"byIsolationMedium",
                'nice_name': "Medium Isolation"
            },
            {
                'expr': nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationTight > 0.5'),
                'style_name':"byIsolationTight",
                'nice_name': "Tight Isolation"
            }
        ],
        "calo":[
            {
                'expr': nTuples["calo"].expr(caloString) ,
                'style_name':"OneOrThreeProng",
                'nice_name': ""
            }
        ]
    }

    fakerate_results = {}

    extra_labels = {}
    extra_labels['fixedCone'] = [ algo_label_fixed, ]
    extra_labels['shrinkingCone'] = [ algo_label_shrinking, ]
    extra_labels['TaNC'] = [ algo_label_tanc, ]
    extra_labels['hps'] = [ algo_label_hps, ]
    extra_labels['calo'] = [ algo_label_calo, ]
        
    for algorithm in numerators:
        # Update the numerator for each efficiency to ensure it is a subset of the
        # denominator
        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue
        
        result_algorithm = makeFakeratePlots(algorithm, extra_labels = extra_labels[algorithm])
        fakerate_results[algorithm] = result_algorithm

    # Make comparison plot of fake-rate after all tau id. criteria are applied
    # for different algorithms
    numerators = {}
    numerators['fixedCone'] = ""
    numerators['shrinkingCone'] = 'OneOrThreeProng'
    numerators['TaNC'] = 'byTaNCfrHalfPercent'
    numerators['hps'] = 'byIsolationMedium'
    numerators['calo'] = 'OneOrThreeProng'
    makeFakerateComparisonPlots(numerators)
            
        
        







