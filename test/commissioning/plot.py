#!/usr/bin/env python

'''
Example of use of TauNtuple software

Author: Evan K. Friis, UC Davis
'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager import TauNtupleManager
import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    file = ROOT.TFile.Open("tauIdEffEDNtuple_qcdDiJet.root", "READ")
    
    # Get the events tree (this can also be a TChain)
    events = file.Get("Events")

    # Get the ntuple we produced
    manager = TauNtupleManager(events, "tauIdEffNtuple")

    # Get the list of collections availabe for our ntuple
    print manager

    # Lets use selectedPatTaus
    selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone")

    # All this helper funciton does it makes it easy for use
    # to interface to TTree::Draw, etc.

    # Get list of variables available
    print selectedPatTaus
    

    my_selection = selectedPatTaus.expr('$pt > 5')

    # The following operators are available
    # + - * / & | < > <= >=

    # There are also tools to make plots.  Here is an example to draw the Pt spectrum
    # for taus in the barrel 
    canvas = ROOT.TCanvas("example", "example", 500, 500)

    pt_hist = plot.draw(events, expression=selectedPatTaus.expr('$pt'), 
                        selection=selectedPatTaus.expr('abs($eta) < 2.5'),
                        binning=(10, 0, 50))

    # pt_hist is a TH1F
    pt_hist.Draw()
    canvas.SaveAs("pt_spectrum.png")

    # Plot track angle 0 TaNC discriminant
    track_angle_hist = plot.draw(events, expression=selectedPatTaus.expr('$TaNCTrackAngle0'),
                        selection=(selectedPatTaus.expr('abs($eta) < 2.5') & selectedPatTaus.expr('$byLeadPionPtCut')),
                        binning=(40, 0, 0.4))

    track_angle_hist.Draw()
    canvas.SaveAs("track_angle_0.png")

    # Plot decay mode
    decay_mode_hist = plot.draw(events, expression=selectedPatTaus.expr('$decayMode'),
                        selection=(selectedPatTaus.expr('abs($eta) < 2.5') & selectedPatTaus.expr('$byLeadPionPtCut')),
                        binning=(14, -0.5, 13.5))

    decay_mode_hist.Draw()
    canvas.SaveAs("decay_mode_hist.png")

    # We also have access to the generator level (w/ no reco matching required) in a seperate ntuple
    true_taus = manager.get_ntuple("tauGenJets")

    true_tau_pt = plot.draw(events, expression=true_taus.expr('$genPt'),
                            selection=(true_taus.expr('abs($genEta) < 2.5')),
                            binning=(20, 0, 100)
                           )
    true_tau_pt.Draw()
    canvas.SaveAs("true_tau_pt.png")


    # We can easily make an efficiency plot as well
    # Let's compute the TaNC eff w.r.t Leading Piont for taus with Pt > 5 
    # with eta < 2.5
    denom_selection = selectedPatTaus.expr('$pt > 5') & \
            selectedPatTaus.expr('abs($eta) < 2.5')  

    numerator_selection = denom_selection & selectedPatTaus.expr('$byTaNCfrOnePercent') &\
            selectedPatTaus.expr('$byLeadPionPtCut > 0.5')

    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    bkg_histo, efficiency = plot.efficiency(
        events,
        expression=selectedPatTaus.expr('$pt'),
        numerator=numerator_selection,
        denominator=denom_selection,
        binning=(10, 0, 50))

    bkg_histo.GetYaxis().SetRangeUser(1e-4, 0.5)
    bkg_histo.SetTitle("")
    bkg_histo.GetXaxis().SetTitle("PT")
    bkg_histo.Draw()
    ROOT.gPad.SetLogy(True)
    efficiency.Draw("s, p")

    canvas.SaveAs("tanc_pt_eff.pdf")

    numerator_selection = denom_selection & selectedPatTaus.expr('$byIsolation') &\
            selectedPatTaus.expr('$byLeadPionPtCut > 0.5')
    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    bkg_histo, efficiency = plot.efficiency(
        events,
        expression=selectedPatTaus.expr('$pt'),
        numerator=numerator_selection,
        denominator=denom_selection,
        binning=(10, 0, 50))

    bkg_histo.GetYaxis().SetRangeUser(1e-4, 0.5)
    bkg_histo.SetTitle("")
    bkg_histo.GetXaxis().SetTitle("PT")
    bkg_histo.Draw()
    efficiency.Draw("s, p")

    canvas.SaveAs("iso_pt_eff.pdf")







