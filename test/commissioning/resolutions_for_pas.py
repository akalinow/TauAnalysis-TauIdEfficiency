#!/usr/bin/env python

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Definition of input files.
import samples_cache as samples

def getResolutionExpression(algorithm, variable):
    result = None
    resolutions ={
        "TaNC": {
            "pt": nTuples["TaNC"].expr('($decayModePt-$genPt)/$genPt'),
            "eta": nTuples["TaNC"].expr('$decayModeEta-$genEta'),
            "phi": nTuples["TaNC"].expr('$decayModePhi-$genPhi'),
            },
        "other":{
            "pt": '($pt-$genPt)/$genPt',
            "eta": '$eta-$genEta',
            "phi": '$phi-$genPhi',
            }
        }
    if algorithm in resolutions.keys() and not algorithm == "other":
        result = resolutions[algorithm][variable]
    else:
        result = nTuples[algorithm].expr(resolutions["other"][variable])
    return result

def makeResolutionPlots( algorithm ):
    for selection in selections[algorithm]:
        plotter.update_style("mc_ztt", **style.MC_STYLES[selection["style_name"]])
        pt_resol = plotter.distribution(
            expression= getResolutionExpression(algorithm, "pt"),
            selection=selection["expr"],
            labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                      style.ZTAUTAU_LABEL_UPPER_LEFT,
                      style.SQRTS_LABEL_UPPER_LEFT,
                      style.PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT],
            binning=(100, -0.5, 0.5),
            x_axis_title = "(P_{T}^{vis}(rec) - P_{T}^{vis}(gen))/P_{T}^{vis}(gen)",
            y_min = 0., y_max=0.175,
            normalize = 1.,
            logy=False
            )

        # Make a pave text w/ mean rms
        stat_label = style.make_mean_rms_pave(pt_resol['samples']['mc_ztt']['plot'])
        stat_label.Draw()
    
        canvas.SaveAs("plots/%s%s_pt_resolution.png"%(algorithm,selection["style_name"]))
        canvas.SaveAs("plots/%s%s_pt_resolution.pdf"%(algorithm,selection["style_name"]))        

        plotter.update_style("mc_ztt", **style.MC_STYLES[selection["style_name"]])
        eta_resol = plotter.distribution(
            expression= getResolutionExpression(algorithm, "eta"),
            selection=selection["expr"],
            labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                      style.ZTAUTAU_LABEL_UPPER_LEFT,
                      style.SQRTS_LABEL_UPPER_LEFT,
                      style.PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT],
            binning=(100, -0.01, 0.01),
            n_divisions = 404,
            x_axis_title = "#eta(rec) - #eta(gen)",
            y_min = 0., y_max=0.175,
            normalize = 1.,
            logy=False
            )

        # Make a pave text w/ mean rms
        stat_label = style.make_mean_rms_pave(eta_resol['samples']['mc_ztt']['plot'])
        stat_label.Draw()
    
        canvas.SaveAs("plots/%s%s_eta_resolution.png"%(algorithm,selection["style_name"]))
        canvas.SaveAs("plots/%s%s_eta_resolution.pdf"%(algorithm,selection["style_name"]))

        phi_resol = plotter.distribution(
            expression=getResolutionExpression(algorithm, "phi"),
            selection=selection["expr"],
            labels = [style.CMS_PRELIMINARY_UPPER_LEFT,
                      style.ZTAUTAU_LABEL_UPPER_LEFT,
                      style.SQRTS_LABEL_UPPER_LEFT,
                      style.PT_ETA_GEN_CUT_TWO_LINE_LABEL_UPPER_LEFT],
            binning=(100, -0.01, 0.01),
            x_axis_title = "#phi(rec) - #phi(gen)",
            y_min = 0., y_max=0.175,
            n_divisions = 404,
            normalize = 1.,
            logy=False
            )

        # Make a pave text w/ mean rms
        stat_label = style.make_mean_rms_pave(phi_resol['samples']['mc_ztt']['plot'])
        stat_label.Draw()
    
        canvas.SaveAs("plots/%s%s_phi_resolution.png"%(algorithm,selection["style_name"]))
        canvas.SaveAs("plots/%s%s_phi_resolution.pdf"%(algorithm,selection["style_name"]))
    

if __name__ == "__main__":

    plotter = PlotManager()

    # Add each sample we want to plot/compare
    # Uncomment to add QCD
    plotter.add_sample(samples.ztautau_mc, "Z->#tau#tau MC", **style.QCD_MC_STYLE_HIST)

    # Normalize everything to the data luminosity
    plotter.set_integrated_lumi(samples.data.effective_luminosity())
    
    # Build the ntuple manager
    ntuple_manager = samples.ztautau_mc.build_ntuple_manager("tauIdEffNtuple")

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
    
    from shrinkingConePlots import pt_binning_fine, eta_binning_fine, phi_binning_fine
    from shrinkingConePlots import lead_pion_selection
    
    # Define the selections to plot
    pfString = "$byLeadTrackFinding > 0.5 & $byLeadTrackPtCut > 0.5 & $byTrackIsolation > 0.5 "
    pfString += " & $byEcalIsolation > 0.5 & ($numChargedParticlesSignalCone == 1 || $numChargedParticlesSignalCone == 3)"
    caloString = "$byLeadTrackFinding > 0.5 &  & $byLeadTrackPtCut > 0.5 & $byIsolation > 0.5 "
    caloString += " & $etSumIsolationECAL < 5 & ($numSignalTracks ==1 || $numSignalTracks ==3)"
    
    selections = {
        "shrinkingCone":[
            {
                'expr':  nTuples["shrinkingCone"].expr(pfString) ,
                "style_name":"OneOrThreeProng",
                'nice_name': "all shrinkingCone",
                },
            ],
        "fixedCone":[
            {
                'expr':  nTuples["fixedCone"].expr(pfString) ,
                "style_name":"OneOrThreeProng",
                'nice_name': "all fixedCone",
                },
            ],
        "TaNC":[
            {
                'expr':  nTuples["TaNC"].expr('$byTaNCfrOnePercent > 0.5') & lead_pion_selection,
                "style_name":"byTaNCfrOnePercent",
                'nice_name': "TaNC 1.00%",
                },
            {
                'expr':  nTuples["TaNC"].expr('$byTaNCfrHalfPercent > 0.5') & lead_pion_selection,
                "style_name":"byTaNCfrHalfPercent",
                'nice_name': "TaNC 0.50%"
                },
            {
                'expr':  nTuples["TaNC"].expr('$byTaNCfrQuarterPercent > 0.5') & lead_pion_selection,
                "style_name":"byTaNCfrQuarterPercent",
                'nice_name': "TaNC 0.25%"
                },
            ],
        "hps":[
            {
                'expr':  nTuples["hps"].expr('$byIsolationLoose > 0.5') & lead_pion_selection,
                "style_name":"byIsolationLoose",
                'nice_name': "Loose Isolation",
                },
            {
                'expr':  nTuples["hps"].expr('$byIsolationMedium > 0.5') & lead_pion_selection,
                "style_name":"byIsolationMedium",
                'nice_name': "Medium Isolation",
                },
            {
                'expr':  nTuples["hps"].expr('$byIsolationTight > 0.5') & lead_pion_selection,
                "style_name":"byIsolationTight",
                'nice_name': "Tight Isolation",
                },
            ],
        "calo":[
            {
                'expr':  nTuples["calo"].expr(caloString) ,
                "style_name":"OneOrThreeProng",
                'nice_name': "all caloTau",
                },
            ],
        }    
    for algorithm in selections:
        import sys
        if sys.argv[1:] != [] and (not algorithm in sys.argv[1:]):
            continue
        #the basic selection is the same for all nTuples.
        base_selection = nTuples[algorithm].expr("$genMatch > 0.5 && $genDecayMode > 1.5 && $genPt > 10. && abs($genEta) < 2.5 ")
    
        for res_info in selections[algorithm]: 
            res_info['expr'] = base_selection & res_info['expr']
            makeResolutionPlots( algorithm )
            
        
        







