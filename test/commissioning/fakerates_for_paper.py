#!/usr/bin/env python

'''
Plot fake-rates of different tau id algorithms for
 o QCD multi-jet
 o QCD muon enriched
 o W --> mu nu + jets
Data compared to Monte Carlo predictions

Author: Christian Veelken

'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import copy

# Define of input files.
import samples_cache as samples

# Define sample names
sample_dijet_data_woPU = samples.data_dijet_woPU
sample_dijet_mc_woPU   = samples.qcddijet_mc
sample_dijet_data_wPU  = samples.data_dijet_wPU
sample_dijet_mc_wPU    = samples.qcddijetPU156bx_mc
sample_ppmux_data      = samples.data_ppmux
sample_ppmux_mc        = samples.ppmux10PU156bx_mc
sample_wjets_data      = samples.data_wjets
sample_wjets_mc        = samples.wjetsPU156bx_mc

# Define integrated luminosity of analyzed dataset
intLumiData = 36100 # nb^-1

# Define denominator (phase-space) for which fake-rates are to be determined
denominator_phase_space = "$jetPt > 10.0 & abs($jetEta) < 2.5"

# Define maximum number of Ntuple entries to process
# (NOTE: use values different from 1000000000 for debugging purposes only !!)
#maxNumEntries = 1000000000
maxNumEntries = 10000

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
                      x_axis_title, y_min = 1e-3, y_max = 0.9, extra_labels = [],
                      maxNumEntries = 1000000000,
                      outputFileName = "fakerates_for_paper.pdf"):

    print("<makeFakeratePlots>:")

    result_algorithm = {}   

    canvas_diff = ROOT.TCanvas("canvas_diff", "canvas_diff", 500, 625)
    canvas_diff.cd()
    
    topPad_diff = ROOT.TPad("topPad_diff", "topPad_diff", 0.00, 0.35, 1.00, 1.00)
    topPad_diff.SetFillColor(10)
    topPad_diff.SetLogy(True)
    topPad_diff.SetTopMargin(0.05)
    topPad_diff.SetLeftMargin(0.15)
    topPad_diff.SetBottomMargin(0.00)
    topPad_diff.SetRightMargin(0.05)
    topPad_diff.SetGridy()
    topPad_diff.Draw()
    topPad_diff.cd()

    # Make fake-rate plots
    eff_dijet = plotter_dijet.efficiency(
        expression, 
        denominator_dijet,
        denominator_dijet & numerator,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Efficiency",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['dijet'] = eff_dijet

    eff_ppmux = plotter.efficiency(
        expression, 
        denominator_ppmux,
        denominator_ppmux & numerator,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Efficiency",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['ppmux'] = eff_ppmux

    eff_wjets = plotter.efficiency(
        expression, 
        denominator_wjets,
        denominator_wjets & numerator,
        binning = binning,
        maxNumEntries = maxNumEntries, verbose = True, x_error_bars = True,
        x_axis_title = x_axis_title,
        y_axis_title = "Efficiency",
        y_min = y_min, y_max = y_max, logy = True)
    result_algorithm['wjets'] = eff_wjets

    # Adjust offset of y-axis title
    #eff_dijet['background'].GetYaxis().SetTitle("Efficiency" )
    eff_dijet['background'].GetYaxis().SetNdivisions(505)
    eff_dijet['background'].GetYaxis().SetTickLength(0.04)
    eff_dijet['background'].GetYaxis().SetTitleOffset(1.2)
    eff_dijet['background'].GetYaxis().SetTitleSize(0.05) 
    eff_dijet['background'].GetYaxis().SetLabelSize(0.05)
    eff_dijet['background'].Draw()

    # Draw fake-rate plots
    eff_dijet['samples'][sample_dijet_data_wPU].Draw('epsame')
    eff_dijet['samples'][sample_dijet_mc_wPU].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_data].Draw('epsame')
    eff_ppmux['samples'][sample_ppmux_mc].Draw('epsame')
    eff_wjets['samples'][sample_wjets_data].Draw('epsame')
    eff_wjets['samples'][sample_wjets_mc].Draw('epsame')

    # Draw legend
    legend = ROOT.TLegend(0.45, 0.68, 0.88, 0.88, "","brNDC")
    legend.SetTextSize(0.03)
    legend.SetFillColor(0);
    legend.SetLineColor(1);
    legend.SetBorderSize(1);
    legend.AddEntry(eff_wjets['samples'][sample_wjets_data.name],         "W #rightarrow #mu #nu Data",       "p")
    legend.AddEntry(eff_wjets['samples'][sample_wjets_mc.name],           "W #rightarrow #mu #nu Simulation", "p")
    legend.AddEntry(eff_dijet_wPU['samples'][sample_dijet_data_wPU.name], "QCDj Data",                        "p")
    legend.AddEntry(eff_dijet_wPU['samples'][sample_dijet_mc_wPU.name],   "QCDj Simulation",                  "p")
    legend.AddEntry(eff_ppmux['samples'][sample_ppmux_data.name],         "QCD#mu Data",                      "p")
    legend.AddEntry(eff_ppmux['samples'][sample_ppmux_mc.name],           "QCD#mu Simulation",                "p")
    legend.Draw();

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
    custom_LUMI_LABEL_UPPER_LEFT = copy.deepcopy(style.LUMI_LABEL_UPPER_LEFT)
    custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.0375)
    custom_LUMI_LABEL_UPPER_LEFT.Clear()
    custom_LUMI_LABEL_UPPER_LEFT.AddText("Data, L = 36.1pb^{-1}")
    custom_LUMI_LABEL_UPPER_LEFT.Draw()
    style.SQRTS_LABEL_UPPER_LEFT.Draw()

    for extra_label in extra_labels:
        extra_label.Draw()

    algo_label.Draw()

    print "outputFileName", outputFileName[:outputFileName.find(".")] + "_%s_diff.png" % algorithm

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.png" % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s.pdf" % algorithm)

    canvas_diff.cd()

    bottomPad_diff = ROOT.TPad("bottomPad_diff", "bottomPad_diff", 0.00, 0.00, 1.00, 0.35)
    bottomPad_diff.SetFillColor(10)
    bottomPad_diff.SetLogy(False)
    bottomPad_diff.SetTopMargin(0.00)
    bottomPad_diff.SetLeftMargin(0.18)
    bottomPad_diff.SetBottomMargin(0.20)
    bottomPad_diff.SetRightMargin(0.05)
    bottomPad_diff.Draw()
    bottomPad_diff.cd()

    # Show (normalized) difference (Data - MC)/MC
    diff_dijet_wPU = plotter_dijet.plot_dist_deviations(
        eff_dijet,
        sample_dijet_mc_wPU.name,
        sample_dijet_data_wPU.name,
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.28, y_max = +0.28, logy = False
    )
    result_algorithm['dijet_DataToMC_diff'] = diff_dijet_wPU
    
    diff_ppmux = plotter_ppmux.plot_dist_deviations(
        eff_ppmux,
        sample_ppmux_mc.name,
        sample_ppmux_data.name,
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.28, y_max = +0.28, logy = False
    )
    result_algorithm['ppmux_DataToMC_diff'] = diff_ppmux
    
    diff_wjets = plotter_wjets.plot_dist_deviations(
        eff_wjets,
        sample_wjets_mc.name,
        sample_wjets_data.name,
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{Data - Simulation}{Simulation}",
        y_min = -0.28, y_max = +0.28, logy = False
    )
    result_algorithm['wjets_DataToMC_diff'] = diff_wjets

    diff_dijet_wPU['background'].GetYaxis().SetNdivisions(505)
    diff_dijet_wPU['background'].GetXaxis().SetTitleOffset(1.1)
    diff_dijet_wPU['background'].GetXaxis().SetTitleSize(0.08)
    diff_dijet_wPU['background'].GetXaxis().SetLabelSize(0.08)
    diff_dijet_wPU['background'].GetXaxis().SetTickLength(0.055)
    diff_dijet_wPU['background'].GetYaxis().SetTitleOffset(0.9)
    diff_dijet_wPU['background'].GetYaxis().SetTitleSize(0.08)
    diff_dijet_wPU['background'].GetYaxis().SetLabelSize(0.08)
    diff_dijet_wPU['background'].GetYaxis().SetTickLength(0.04)

    canvas_diff.cd()
    canvas_diff.Update()

    print "outputFileName", outputFileName[:outputFileName.find(".")] + "_%s_diff.png" % algorithm

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.png" % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_diff.pdf" % algorithm)

    # Show difference in fake-rate between QCD di-jet with/without pile-up
    topPad_diff.cd()

    # Draw fake-rate plots
    eff_dijet['samples'][sample_dijet_data_woPU.name].Draw('epsame')
    eff_dijet['samples'][sample_dijet_mc_woPU.name].Draw('epsame')
    eff_dijet['samples'][sample_dijet_data_wPU.name].Draw('epsame')
    eff_dijet['samples'][sample_dijet_mc_wPU.name].Draw('epsame')

    # Draw legend
    legend = ROOT.TLegend(0.45, 0.68, 0.88, 0.88, "","brNDC")
    legend.SetTextSize(0.03)
    legend.SetFillColor(0);
    legend.SetLineColor(1);
    legend.SetBorderSize(1);
    legend.AddEntry(eff_dijet_woPU['samples'][sample_dijet_data_woPU.name], "QCDj Data, wo. PU", "p")
    legend.AddEntry(eff_dijet_woPU['samples'][sample_dijet_mc_woPU.name], "QCDj Simulation, wo. PU", "p")
    legend.AddEntry(eff_dijet_wPU['samples'][sample_dijet_data_wPU.name], "QCDj Data, w. PU", "p")
    legend.AddEntry(eff_dijet_wPU['samples'][sample_dijet_mc_wPU.name], "QCDj Simulation, w. PU", "p")
    legend.Draw();

    # Draw CMS preliminary, luminosity and center-of-mass energy labels
    style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
    custom_LUMI_LABEL_UPPER_LEFT = copy.deepcopy(style.LUMI_LABEL_UPPER_LEFT)
    custom_LUMI_LABEL_UPPER_LEFT.SetTextSize(0.0375)
    custom_LUMI_LABEL_UPPER_LEFT.Clear()
    custom_LUMI_LABEL_UPPER_LEFT.AddText("Data, L = 36.1pb^{-1}")
    custom_LUMI_LABEL_UPPER_LEFT.Draw()
    style.SQRTS_LABEL_UPPER_LEFT.Draw()

    for extra_label in extra_labels:
        extra_label.Draw()

    algo_label.Draw()

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_pileup.png" % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_pileup.pdf" % algorithm)

    canvas_diff.cd()

    bottomPad_diff.cd()

    # Show (normalized) difference (Data - MC)/MC
    diff_dijet_pileup_data = plotter_dijet.plot_dist_deviations(
        eff_dijet,
        sample_dijet_data_woPU.name,
        sample_dijet_data_wPU.name,
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{w. PU - wo. PU}{wo. PU}",
        y_min = -0.28, y_max = +0.28, logy = False
    )
    result_algorithm['dijet_Data_pileup_diff'] = diff_dijet_pileup_data
        
    diff_dijet_pileup_mc = plotter_ppmux.plot_dist_deviations(
        eff_dijet,
        sample_dijet_mc_woPU.name,
        sample_dijet_mc_wPU.name,
        x_axis_title = x_axis_title,
        y_axis_title = "#frac{w. PU - wo. PU}{wo. PU}",
        y_min = -0.28, y_max = +0.28, logy = False
    )
    result_algorithm['dijet_MC_pileup_diff'] = diff_ppmux

    diff_dijet_pileup_data['background'].GetYaxis().SetNdivisions(505)
    diff_dijet_pileup_data['background'].GetXaxis().SetTitleOffset(1.1)
    diff_dijet_pileup_data['background'].GetXaxis().SetTitleSize(0.08)
    diff_dijet_pileup_data['background'].GetXaxis().SetLabelSize(0.08)
    diff_dijet_pileup_data['background'].GetXaxis().SetTickLength(0.055)
    diff_dijet_pileup_data['background'].GetYaxis().SetTitleOffset(0.9)
    diff_dijet_pileup_data['background'].GetYaxis().SetTitleSize(0.08)
    diff_dijet_pileup_data['background'].GetYaxis().SetLabelSize(0.08)
    diff_dijet_pileup_data['background'].GetYaxis().SetTickLength(0.04)

    canvas_diff.cd()
    canvas_diff.Update()

    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_pileup_diff.png" % algorithm)
    canvas_diff.SaveAs(outputFileName[:outputFileName.find(".")] + "_%s_pileup_diff.pdf" % algorithm)

    return result_algorithm
       
if __name__ == "__main__":

    # define PlotManagers for QCD multi-jet, QCD muon enriched and W + jets samples
    plotter_dijet = PlotManager()
    plotter_dijet.add_sample(sample_dijet_data_woPU, "QCDj Data, wo. PU",                **style.DATA_STYLE)
    plotter_dijet.add_sample(sample_dijet_mc_woPU,   "QCDj Simulation, wo. PU",          **style.MINBIAS_MC_STYLE)
    plotter_dijet.add_sample(sample_dijet_data_wPU,  "QCDj Data, w. PU",                 **style.DATA_STYLE)
    plotter_dijet.add_sample(sample_dijet_mc_wPU,    "QCDj Simulation, w. PU",           **style.QCD_MC_PYTHIA6_STYLE_HIST)
    plotter_dijet.set_integrated_lumi(intLumiData)
    
    plotter_ppmux = PlotManager()
    plotter_ppmux.add_sample(sample_ppmux_data,      "QCD#mu Data",                      **style.DATA_STYLE)
    plotter_ppmux.add_sample(sample_ppmux_mc, "       QCD#mu Simulation",                **style.QCD_MC_STYLE_HIST)
    plotter_ppmux.set_integrated_lumi(intLumiData)
        
    plotter_wjets = PlotManager()
    plotter_wjets.add_sample(sample_wjets_data,      "W #rightarrow #mu #nu Data",       **style.DATA_STYLE)
    plotter_wjets.add_sample(sample_wjets_mc,        "W #rightarrow #mu #nu Simulation", **style.WJETS_MC_STYLE_HIST)
    plotter_wjets.set_integrated_lumi(intLumiData)
    
    # Build the ntuple manager
    ntuple_manager = sample_wjets_data.build_ntuple_manager("tauIdEffNtuple")

    nTuples = {
        "shrinkingCone" : ntuple_manager.get_ntuple("shrinking"),
        "fixedCone"     : ntuple_manager.get_ntuple("fixed"),
        "TaNC"          : ntuple_manager.get_ntuple("shrinking"),
        "hps"           : ntuple_manager.get_ntuple("hps"),
        "calo"          : ntuple_manager.get_ntuple("calo")
    }
    
    hlt = ntuple_manager.get_ntuple("patTriggerEvent")

    # Define binning options
    binning_pt  = (0, 10,  15,  20,  25, 30, 40, 50, 65, 80, 100)
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
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrOnePercent > 0.5'),
            'style_name' : "byTaNCfrOnePercent",
            'nice_name'  : "TaNC loose"
        },
        "TaNC_medium" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrHalfPercent > 0.5'),
            'style_name' :"byTaNCfrHalfPercent",
            'nice_name'  : "TaNC medium"
        },
        "TaNC_tight" : {
            'expr'       : nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTaNCfrQuarterPercent > 0.5'),
            'style_name' : "byTaNCfrQuarterPercent",
            'nice_name'  : "TaNC tight"
        },
        "hps_loose" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationLoose > 0.5'),
            'style_name' : "byIsolationLoose",
            'nice_name'  : "HPS loose"
        },
        "hps_medium" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationMedium > 0.5'),
            'style_name' : "byIsolationMedium",
            'nice_name'  : "HPS medium"
        },
        "hps_tight" : {
            'expr'       : nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationTight > 0.5'),
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
        
    #for algorithm in numerators:
    for algorithm_discriminator in [ 'TaNC_loose', 'hps_loose' ]:

        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue
        
        algorithm     = algorithm_discriminator[:algorithm_discriminator.find("_")]
        discriminator = algorithm_discriminator[algorithm_discriminator.find("_") + 1:]
        print("algorithm = %s, discriminator = %s" % (algorithm, discriminator))
        
        # Define the denominators
        denominator_e_mu_veto = nTuples[algorithm].expr("$pfElectronMVA < 0.6 & $againstMuon > 0.5")
        denominator_dijet     = nTuples[algorithm].expr(denominator_phase_space) & denominator_e_mu_veto \
                               & hlt.expr('$hltJet15Ubit > 0.5') & nTuples[algorithm].expr("$probe > 0.5")
        denominator_ppmux     = nTuples[algorithm].expr(denominator_phase_space) & denominator_e_mu_veto
        denominator_wjets     = nTuples[algorithm].expr(denominator_phase_space) & denominator_e_mu_veto

        makeFakeratePlots(algorithm, numerators[algorithm_discriminator]['expr'],
                          denominator_dijet, plotter_dijet, 
                          denominator_ppmux, plotter_ppmux, 
                          denominator_wjets, plotter_wjets, 
                          nTuples[algorithm].expr('$jetPt'), binning_pt,
                          'Jet P_{T} [GeV/c]', 1e-4, 9.9, extra_labels[algorithm_discriminator],
                          maxNumEntries,
                          "plots/fakerates_for_pas_jetPt.pdf")
        ##makeFakeratePlots(algorithm, numerator[algorithm]['expr'],
        ##                  denominator_dijet, plotter_dijet, 
        ##                  denominator_ppmux, plotter_ppmux, 
        ##                  denominator_wjets, plotter_wjets,
        ##                  nTuples[algorithm].expr('$jetEta'), binning_eta,
        ##                  'Jet #eta', 1e-4, 9.9, extra_labels[algorithm],
        ##                  maxNumEntries,
        ##                  "plots/fakerates_for_pas_jetEta.pdf")
        ##makeFakeratePlots(algorithm, numerator[algorithm]['expr'],
        ##                  denominator_dijet, plotter_dijet, 
        ##                  denominator_ppmux, plotter_ppmux, 
        ##                  denominator_wjets, plotter_wjets,
        ##                  nTuples[algorithm].expr('$jetPhi'), binning_phi,
        ##                  'Jet #phi', 1e-4, 9.9, extra_labels[algorithm],
        ##                  maxNumEntries,
        ##                  "plots/fakerates_for_pas_jetPhi.pdf")







