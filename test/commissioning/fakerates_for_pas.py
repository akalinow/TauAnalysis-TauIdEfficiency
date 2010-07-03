#!/usr/bin/env python

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Definition of input files.
import samples_cache as samples

def makeFakeratePlots( algorithm ):
    denominator_section = "$probe > 0.5 & $jetPt > 10.0 & abs($jetEta) < 2.5"
    denominator = hlt.expr('$hltJet15U > 0.5') & nTuples[algorithm].expr(denominator_section)
    for eff_info in numerators[algorithm]: 
        eff_info['expr'] = denominator & eff_info['expr']
    pt_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetPt'), 
        denominator,
        numerators[algorithm], 
        binning=pt_binning_fine, 
        y_min = 1e-4,
        y_max = 10,
        x_axis_title='Jet P_{T} [GeV/c]',
        labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                  style.LUMI_LABEL_UPPER_LEFT,
                  style.SQRTS_LABEL_UPPER_LEFT,
                  style.ETA_CUT_LABEL_UPPER_LEFT
                  ],
        logy = True)
    
    pt_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_pt_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_pt_for_pas.pdf"%(algorithm))
    
    eta_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetEta'), 
        denominator,
        numerators[algorithm], 
        binning=eta_binning_fine, 
        y_min = 1e-4,
        y_max = 10,
        x_axis_title='Jet #eta',
        labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                  style.LUMI_LABEL_UPPER_LEFT,
                  style.SQRTS_LABEL_UPPER_LEFT,
                  style.PT_CUT_LABEL_UPPER_LEFT
                  ],
        logy = True)
    
    eta_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_eta_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_eta_for_pas.pdf"%(algorithm))
    
    phi_effs = plotter.multi_efficiency(
        nTuples[algorithm].expr('$jetPhi'), 
        denominator,
        numerators[algorithm], 
        binning=phi_binning_fine, 
        y_min = 1e-4,
        y_max = 10,
        x_axis_title='Jet #phi',
        labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                  style.LUMI_LABEL_UPPER_LEFT,
                  style.SQRTS_LABEL_UPPER_LEFT,
                  style.PT_ETA_CUT_TWO_LINE_LABEL_UPPER_LEFT,
                  ],
        logy = True)
    
    phi_effs['legend'].make_legend(0.45, 0.68, 0.88, 0.88).Draw()
    
    canvas.SaveAs("plots/%s_multifake_vs_phi_for_pas.png"%(algorithm))
    canvas.SaveAs("plots/%s_multifake_vs_phi_for_pas.pdf"%(algorithm))

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
        "calo": ntuple_manager.get_ntuple("patCaloTausDijetTagAndProbe"),
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
                'expr':  nTuples["shrinkingCone"].expr(pfString) ,
                "style_name":"OneOrThreeProng",
                ##'nice_name': "all shrinkingCone",
                'nice_name': "",
                },
            ],
        "fixedCone":[
            {
                'expr':  nTuples["fixedCone"].expr(pfString) ,
                "style_name":"OneOrThreeProng",
                ##'nice_name': "all fixedCone",
                'nice_name': "",
                },
            ],
        "TaNC":[
            {
                'expr':  nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrOnePercent > 0.5'),
                "style_name":"byTaNCfrOnePercent",
                'nice_name': "TaNC 1.00%",
                },
            {
                'expr':  nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrHalfPercent > 0.5'),
                "style_name":"byTaNCfrHalfPercent",
                'nice_name': "TaNC 0.50%"
                },
            {
                'expr':  nTuples["TaNC"].expr('$byLeadTrackFinding > 0.5 & $byTaNCfrQuarterPercent > 0.5'),
                "style_name":"byTaNCfrQuarterPercent",
                'nice_name': "TaNC 0.25%"
                },
            ],
         "hps":[
            {
                'expr':  nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationLoose > 0.5'),
                "style_name":"byIsolationLoose",
                'nice_name': "Loose Isolation",
                },
            {
                'expr':  nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationMedium > 0.5'),
                "style_name":"byIsolationMedium",
                'nice_name': "Medium Isolation",
                },
            {
                'expr':  nTuples["hps"].expr('$byLeadTrackFinding > 0.5 & $byIsolationTight > 0.5'),
                "style_name":"byIsolationTight",
                'nice_name': "Tight Isolation",
                },
            ],
        "calo":[
            {
                'expr':  nTuples["calo"].expr(caloString) ,
                "style_name":"OneOrThreeProng",
                'nice_name': "",
                },
            ],
        }
        
    for algorithm in numerators:
        # Update the numerator for each efficiency to ensure it is a subset of the
        # denominator
        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue
        makeFakeratePlots( algorithm )
            
        
        







