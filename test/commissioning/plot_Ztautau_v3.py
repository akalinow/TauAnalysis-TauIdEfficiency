#!/usr/bin/env python

'''
Plot efficiency of different tau id algorithms for ZTT hadronic taus

Authors: Aruna Nayak, Evan K. Friis

'''

import ROOT
import samples_cache as samples
from TauAnalysis.TauIdEfficiency.ntauples.plotting import draw
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)

    # Get ZTT samples
    ztt = samples.ztautau_mc

    # Get the ntuple we produced
    ntuple_manager = ztt.build_ntuple_manager("tauIdEffNtuple")

    # Get generator level tau ntuple
    genTaus = ntuple_manager.get_ntuple("tauGenJets")

    # Binning for PT
    pt_bins = (0,2.5,5.,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,
               35,37.5,40,42.5,45,47.5,50,52.5,55,57.5,60,63,66,69,72,
               76,80,85,90,95,100,110,120,130,140,150)

    # Define styles of all of the numerators
    numerators = {
        'matching' : {
            'expr_str': '1',
            'label' : "Matched",
        },
        'byLeadTrackFinding': {
            'expr_str': '$byLeadTrackFinding',
            'label' : "Lead Track Finding",
        },
        'byLeadTrackPtCut': {
            'expr_str': '$byLeadTrackPtCut',
            'label' : "Lead Track P_{T} Cut",
        },
        'byTrackIsolation': {
            'expr_str': '$byTrackIsolation',
            'label' : "Charged Hadron Isolation",
        },
        'byEcalIsolation': {
            'expr_str': '$byEcalIsolation',
            'label' : "Gamma Isolation",
        },
        # For calotau
        'byIsolation' : {
            'expr_str': '$byIsolation',
            'label' : 'Isolation'
        },
        'OneOrThreeProng': {
            'expr_str': '$numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3',
            'label' : "1 or 3 Prong",
        },
        'byEcalIsolation_calo' : {
            'expr_str': '$etSumIsolationECAL < 5',
            'label' : 'EcalIsolation'
        },
        'OneOrThreeProng_calo': {
            'expr_str': '$numSignalTracks ==1 || $numSignalTracks ==3',
            'label' : "1 or 3 Prong",
        },
        'byTaNCfrOnePercent': {
            'expr_str': '$byTaNCfrOnePercent',
            'label' : "TaNC 1.00%",
        },
        'byTaNCfrHalfPercent': {
            'expr_str': '$byTaNCfrHalfPercent',
            'label' : "TaNC 0.50%",
        },
        'byTaNCfrQuarterPercent': {
            'expr_str': '$byTaNCfrQuarterPercent',
            'label' : "TaNC 0.25%",
        },
        'byIsolationLoose' : {
            'expr_str': '$byIsolationLoose',
            'label': "Medium Isolation",
        },
        'byIsolationMedium' : {
            'expr_str': '$byIsolationMedium',
            'label': "Medium Isolation",
        },
        'byIsolationTight' : {
            'expr_str': '$byIsolationTight',
            'label': "Tight Isolation",
        },
    }

    parameterizations = {
        'pt' : {
            'expr_str': '$genPt',
            'binning': pt_bins,
            'label': 'Generated #tau visible P_{T}',
        }, 
        'eta' : {
            'expr_str': '$genEta',
            'binning': (60, -3, 3),
            'label': 'Generated #tau visible #eta',
        },
    }

    # Define 'orders' for tau discriminators
    # traditional tau ID sequence
    standard_sequence = [
        'matching',
        'byLeadTrackFinding',
        'byLeadTrackPtCut',
        'byTrackIsolation',
        'byEcalIsolation',
        'OneOrThreeProng'
    ]

    tanc_sequence = [
        'matching',
        'byLeadTrackPtCut',
        'byTaNCfrOnePercent',
        'byTaNCfrHalfPercent',
        'byTaNCfrQuarterPercent',
    ]

    hps_sequence = [
        'matching',
        'byLeadTrackFinding',
        'byIsolationLoose',
        'byIsolationMedium',
        'byIsolationTight'
    ]

    calo_sequence = [
        'matching',
        'byLeadTrackFinding',
        'byLeadTrackPtCut',
        'byIsolation',
        'byEcalIsolation_calo',
        'OneOrThreeProng_calo'
   ]

    # Match up sequences to tau algos
    sequences_and_algos = [
        (standard_sequence, "iso", "patPFTausDijetTagAndProbeShrinkingCone"),
        (standard_sequence, "iso", "patPFTausDijetTagAndProbeFixedCone"),
        (tanc_sequence, "tanc", "patPFTausDijetTagAndProbeShrinkingCone"),
        (hps_sequence, "hps", "patPFTausDijetTagAndProbeHPS"),
        (calo_sequence, "calo", "patCaloTausDijetTagAndProbe"),
    ]

    #denom_selection = genTaus.expr(
    #    '$genPt > 5 && abs($genEta) < 2.5 && $genDecayMode > 1.5')
    denom_selection = genTaus.expr(
        '$genPt > 10 && abs($genEta) < 2.5 && $genDecayMode > 1.5')
    
    #denom_selection_from_reco_str = '$jetPt > 15 && abs($jetEta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5 && $genPt > 5 && abs($genEta) < 2.5'
    denom_selection_from_reco_str = \
            '$pt > 10 && abs($eta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5 && $genPt > 10 && abs($genEta) < 2.5'
    

    ztt_events = list(ztt.events_and_weights())[0][0]

    c1 = ROOT.TCanvas("c1","",0,0,500,500)
    c1.SetGridy(1)
    
    for sequence, sequence_name, algo in sequences_and_algos:
        print "Plotting", sequence_name, "for", algo
        # Get the ntuple
        ntuple = ntuple_manager.get_ntuple(algo)
        # Loop over different horizontal axes
        for x_var, x_var_info in parameterizations.iteritems():
            #t1 = ROOT.TLegend(0.6, 0.10, 0.7, 0.35, "","brNDC")
            t1 = ROOT.TLegend(0.45, 0.65, 0.88, 0.85, "","brNDC")
            t1.SetTextSize(0.03)
            t1.SetFillColor(0);
            t1.SetLineColor(1);
            t1.SetBorderSize(1);
            # build denominator
            print "Building denominator for", x_var
            denominator = draw(
                ztt_events,
                expression = genTaus.expr(x_var_info['expr_str']),
                selection = denom_selection,
                binning = x_var_info['binning'],
                output_name = '_'.join(["denom", x_var, sequence_name, algo]),
            )
            # Loop over sequence of numerators
            numerator_effs = []
            running_cut = ntuple.expr(denom_selection_from_reco_str)
            for numerator_name in sequence:
                # Get the description of this numerator
                numerator_info = numerators[numerator_name]
                print "Building numerator for", numerator_name
                running_cut = ntuple.expr(numerator_info['expr_str']) & running_cut
                # Draw numerator
                print x_var_info['expr_str']
                print running_cut
                numerator = draw(
                    ztt_events,
                    expression = ntuple.expr(x_var_info['expr_str']),
                    selection = running_cut,
                    binning = x_var_info['binning'],
                    output_name = '_'.join([numerator_name, x_var, sequence_name, algo]),
                )

                my_eff = ROOT.TGraphAsymmErrors(numerator, denominator)
                # Update style
                style.update_histo_style(
                    my_eff, style.EFFICIENCY_STYLES[numerator_name])
                numerator_effs.append(my_eff)
                t1.AddEntry(my_eff, numerator_info['label'],"P")

            # Make the actual plots
            print "Building plots"
            for index, numerator_eff in enumerate(numerator_effs):
                # Plot the background on the first go-round
                if index == 0:
                    numerator_eff.Draw("Ap")
                    numerator_eff.GetHistogram().GetXaxis().SetTitle(x_var_info['label'])
                    numerator_eff.GetHistogram().GetYaxis().SetTitle("Efficiency")
                    numerator_eff.GetHistogram().GetYaxis().SetRangeUser(0, 1.5)
                else:
                    numerator_eff.Draw("p,same")
            # Draw legend
            t1.Draw()
            # Draw the preliminary label
            style.CMS_PRELIMINARY_UPPER_LEFT.Draw()
            # Save the plot
            c1.SaveAs("plots/%s.png" % '_'.join(['ztt', algo, sequence_name, 'vs', x_var]))
            c1.SaveAs("plots/%s.pdf" % '_'.join(['ztt', algo, sequence_name, 'vs', x_var]))

                    


                



            



