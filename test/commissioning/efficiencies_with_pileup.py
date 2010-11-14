#!/usr/bin/env python

'''
Plot efficiency of different tau id algorithms for ZTT hadronic taus
as function of the number of vertices reconstructed in each event

Authors: Christian Veelken

'''

import ROOT
import samples_cache as samples
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
from TauAnalysis.TauIdEfficiency.ntauples.plotting import draw, efficiencyLogHack
import TauAnalysis.TauIdEfficiency.ntauples.styles as style
import sys
import copy

# Define labels for different algorithms
algo_label_template = ROOT.TPaveText(0.50, 0.14, 0.88, 0.18, "NDC")
algo_label_template.SetTextAlign(33)
algo_label_template.SetTextSize(0.035)
algo_label_template.SetTextColor(2)
algo_label_template.SetFillStyle(0)
algo_label_template.SetBorderSize(0)

algo_label_shrinking = copy.deepcopy(algo_label_template)
algo_label_shrinking.AddText("Shrinking Cone")

algo_label_tanc = copy.deepcopy(algo_label_template)
algo_label_tanc.AddText("TaNC")

algo_label_hps = copy.deepcopy(algo_label_template)
algo_label_hps.AddText("HPS")

algo_label_calo = copy.deepcopy(algo_label_template)
algo_label_calo.AddText("TCTau")

algo_labels = {
    'shrinkingCone' : algo_label_shrinking,
    'TaNC'          : algo_label_tanc,
    'hps'           : algo_label_hps,
    'calo'          : algo_label_calo
}    

# Define colors for different vertex multiplicities
marker_styles = {
    1 : {
        'color' : 1,
        'style' : 20
    },
    2 : {
        'color' : 596,
        'style' : 21
    },
    3 : {
        'color' : 877,
        'style' : 23
    },
    4 : {
        'color' : 628,
        'style' : 24
    },
    5 : {
        'color' : 797,
        'style' : 25
    },
    6 : {
        'color' : 396,
        'style' : 26
    }
}

# Define numerators for different tau id. algorithms
numerator_selections_by_algorithm = {
    'shrinkingCone' : '$byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                     + ' && $byTrackIsolation > 0.5 && $byEcalIsolation > 0.5' \
                     + ' && ($numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3)',
    'TaNC'          : '$byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                     + ' && $byTaNCfrHalfPercent > 0.5',
    'hps'           : '$byLeadTrackFinding > 0.5' \
                     + ' && $byIsolationMedium',
    'calo'          : '$byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                     + ' && $byIsolation > 0.5 && $etSumIsolationECAL < 5' \
                     + ' && ($numSignalTracks == 1 || $numSignalTracks == 3)'
}

algorithms = [ 'shrinkingCone', 'TaNC', 'hps', 'calo' ]

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)

    # Get ZTT samples
    #ztt = samples.zttPU156bx_mc # particleFlow
    ztt = samples.zttPU156bxPFnoPileUp_mc # pfNoPileUp

    # Build the plot manager.  The plot manager keeps track of all the samples
    # and ensures they are correctly normalized w.r.t. luminosity.  See 
    # samples.py for available samples.

    plotter = PlotManager()

    # Add each sample we want to plot/compare
    plotter.add_sample(ztt, "Z #rightarrow #ell #ell, BX156", **style.QCD_MC_STYLE_HIST)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data_wjets.effective_luminosity())

    # Get the ntuple we produced
    ntuple_manager = ztt.build_ntuple_manager("tauIdEffNtuple")

    # Get generator level tau ntuple
    genTaus = ntuple_manager.get_ntuple("tauGenJets")

    # Get ntuple with vertex information
    vertex = ntuple_manager.get_ntuple("offlinePrimaryVertices")

    canvas = ROOT.TCanvas("canvas", "canvas", 500, 500)
    canvas.SetLeftMargin(0.11)
    canvas.SetGridy(1)
    
    # Make plots of vertex multiplicity
    vertex_multiplicity = plotter.distribution(
        expression = vertex.expr('$numVerticesPtGt10'),
        selection = vertex.expr('$numVerticesPtGt10 > -1'),
        extra_labels = [],
        binning = (10, -0.5, 9.5),
        x_axis_title = "Num. Vertices with sum(track P_{T}) > 10 GeV",
        y_axis_title = "Number of Events",
        y_min = 6.5e-1, y_max = 1.e+3, logy = True,
    )
    vertex_multiplicity['result'].GetYaxis().SetTitleOffset(1.4)
    vertex_multiplicity['result'].Draw("hist")
    canvas.SaveAs("plots/numVerticesPtGt10_efficiency_with_pileup.png")
    canvas.SaveAs("plots/numVerticesPtGt10_efficiency_with_pileup.pdf")

    canvas.SetLogy(False)
           
    # Binning for PT
    pt_bins = ( 0, 10, 15, 20, 25, 30, 40, 60, 80, 100, 150 )

    parameterizations = {
        'pt' : {
            'expr_str': '$genPt',
            'binning': pt_bins,
            'label': 'Generated #tau visible P_{T} [GeV/c]',
            'cutLabel': style.ETAVIS_CUT_LABEL_UPPER_LEFT,
        }, 
        'eta' : {
            'expr_str': '$genEta',
            'binning': (10, -2.5, 2.5),
            'label': 'Generated #tau visible #eta',
            'cutLabel': style.PTVIS_CUT_LABEL_UPPER_LEFT,
        }
    }

    nTuples = {
        "shrinkingCone" : "patPFTausDijetTagAndProbeShrinkingCone",
        "TaNC"          : "patPFTausDijetTagAndProbeShrinkingCone",
        "hps"           : "patPFTausDijetTagAndProbeHPS",
        "calo"          : "patCaloTausDijetTagAndProbe"
    }

    
    denominator_selection_base = '$genDecayMode > 1.5 && $genPt > 10 && abs($genEta) < 2.5'

    numerator_selection_base = '$genMatch > 0.5 && $genDecayMode > 1.5 && $genPt > 10 && abs($genEta) < 2.5'

    ztt_events = list(ztt.events_and_weights())[0][0]

    efficiency_results = {}

    extra_labels = {}
    extra_labels['shrinkingCone'] = [ algo_label_shrinking, ]
    extra_labels['TaNC'] = [ algo_label_tanc, ]
    extra_labels['hps'] = [ algo_label_hps, ]
    extra_labels['calo'] = [ algo_label_calo, ]
    
    for algorithm in algorithms:
        print "Plotting", algorithm
        # Get the ntuple
        ntuple = ntuple_manager.get_ntuple(nTuples[algorithm])
        # Loop over different horizontal axes
        for x_var, x_var_info in parameterizations.iteritems():

            # build legend
            legend = ROOT.TLegend(0.45, 0.68, 0.88, 0.88, "","brNDC")
            legend.SetTextSize(0.03)
            legend.SetFillColor(0);
            legend.SetLineColor(1);
            legend.SetBorderSize(1);

            numerator_effs = []
            for numVertices in [ 1, 2, 3, 4, -5 ]:

                vertex_multiplicity_selection = None
                if numVertices > 0:
                    vertex_multiplicity_selection = '$numVerticesPtGt10 == %i' % numVertices
                else:
                    vertex_multiplicity_selection = '$numVerticesPtGt10 >= %i' % abs(numVertices)

                # Build denominator histogram
                print "Building denominator for", x_var + ",", "num. vertices = %i" % numVertices
                denominator_selection = genTaus.expr(denominator_selection_base) \
                                       & vertex.expr(vertex_multiplicity_selection)
                denominator = draw(
                    ztt_events,
                    expression = genTaus.expr(x_var_info['expr_str']),
                    selection = denominator_selection,
                    binning = x_var_info['binning'],
                    output_name = '_'.join(["denominator%i" % abs(numVertices), x_var, algorithm])
                )
                print "denominator: entries = %i" % denominator.GetEntries()

                # Build numerator histogram
                print "Building numerator for", x_var + ",", algorithm + ",", "num. vertices = %i" % numVertices
                numerator_selection = ntuple.expr(numerator_selection_base) & ntuple.expr(numerator_selections_by_algorithm[algorithm]) \
                                     & vertex.expr(vertex_multiplicity_selection)
                numerator = draw(
                    ztt_events,
                    expression = ntuple.expr(x_var_info['expr_str']),
                    selection = numerator_selection,
                    binning = x_var_info['binning'],
                    output_name = '_'.join(["numerator%i" % abs(numVertices), x_var, algorithm])
                )
                print "numerator: entries = %i" % numerator.GetEntries()
                
                # FIXME: clean this up
                from math import sqrt
                nNum = float(numerator.Integral())
                nDenom = denominator.Integral()
                err = 1/nDenom*sqrt(nNum*(1 - nNum/nDenom) )
                efficiencyLogHack("%s -> %i: " % (nTuples[algorithm], abs(numVertices)), timestamp = True)
                efficiencyLogHack("%e / %e = %e +- %e\n" % (nNum, nDenom, nNum/nDenom, err))
                # end cleanup

                # Create efficiency graph = numerator/denominator histogram
                my_eff = ROOT.TGraphAsymmErrors(numerator, denominator)
                                    
                # Overwrite the additional spacing by root
                if len(x_var_info['binning']) == 3:
                    my_eff.GetXaxis().SetRangeUser(x_var_info['binning'][1],
                                                   x_var_info['binning'][2])
                else:
                    my_eff.GetXaxis().SetRangeUser(min(x_var_info['binning']),
                                                   max(x_var_info['binning']))
                    
                # Update style                
                my_eff.SetMarkerColor(marker_styles[abs(numVertices)]['color'])
                my_eff.SetMarkerStyle(marker_styles[abs(numVertices)]['style'])
                my_eff.SetLineColor(marker_styles[abs(numVertices)]['color'])
                
                numerator_effs.append(my_eff)

                label = ""
                if numVertices <= 0:
                    label += "#geq "
                label += "%i " % abs(numVertices)
                if numVertices == 1:
                    label += "Vertex"
                else:
                    label += "Vertices"
                legend.AddEntry(my_eff, label, "P")
                
            # Make the actual plots
            print "Building plots"
            for index, numerator_eff in enumerate(numerator_effs):
                # Plot the background on the first go-round
                if index == 0:
                    numerator_eff.Draw("Ap")
                    numerator_eff.GetHistogram().GetXaxis().SetTitle(x_var_info['label'])
                    numerator_eff.GetHistogram().GetYaxis().SetTitle("Efficiency %s" % algorithm)
                    numerator_eff.GetHistogram().GetYaxis().SetTitleOffset(1.2)
                    numerator_eff.GetHistogram().GetYaxis().SetRangeUser(0, 1.5)
                else:
                    numerator_eff.Draw("p, same")

            # Draw legend
            legend.Draw()
            
            # Draw the preliminary label
            style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
            style.ZTAUTAU_LABEL_UPPER_LEFT.Draw()
            style.SQRTS_LABEL_UPPER_LEFT.Draw()
            x_var_info["cutLabel"].Draw()
            for extra_label in extra_labels[algorithm]:
                extra_label.Draw()

            # Save the plot
            canvas.SaveAs("plots/%s.png" % '_'.join([algorithm, "efficiency_with_pileup", 'vs', x_var]))
            canvas.SaveAs("plots/%s.pdf" % '_'.join([algorithm, "efficiency_with_pileup", 'vs', x_var]))


