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
    file = ROOT.TFile("example_ntuple.root", "READ")
    
    # Get the events tree (this can also be a TChain)
    events = file.Get("Events")

    # Get the ntuple we produced
    manager = TauNtupleManager(events, "exampleNtuple")

    # Get the list of collections availabe for our ntuple
    print manager
    # returns
    #== TauNtupleManager
    #===  Available ntuples for label exampleNtuple
    #==== selectedPatTaus

    # So lets use selectedPatTaus
    selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone")

    # All this helper funciton does it makes it easy for use
    # to interface to TTree::Draw, etc.

    # Get list of variables available
    print selectedPatTaus
    # returns
    # Tau Ntuple - collection selectedPatTaus
    #  pt
    #  eta
    #  byLeadPionPtCut
    #  byIsolation
    #  byLeadPionPtCut
    #  byTaNCfrOnePercent
    
    # Now it easy to create expressions using the expr method
    print selectedPatTaus.expr('$pt')
    # returns
    # exampleNtuple#selectedPatTaus#pt

    # You can also build selection string, with operator overloading
    # or without.  The following are logically equivalent
    my_selection = selectedPatTaus.expr('$pt > 20 && $pt < 50')
    print my_selection

    # or
    my_first_selection = selectedPatTaus.expr('$pt') > 20 
    my_second_selection = selectedPatTaus.expr('$pt') < 50 
    # NB the notation for logical AND
    my_selection = my_first_selection & my_second_selection
    print my_selection

    # The following operators are available
    # + - * / & | < > <= >=

    # There are also tools to make plots.  Here is an example to draw the Pt spectrum
    # for taus in the barrel 
    canvas = ROOT.TCanvas("example", "example", 500, 500)

    pt_hist = plot.draw(events, expression=selectedPatTaus.expr('$pt'), 
                        selection=selectedPatTaus.expr('$eta > -2.1 && $eta < +2.1'),
                        binning=(10, 0, 50))

    # pt_hist is a TH1F
    pt_hist.Draw()
    canvas.SaveAs("pt_spectrum.png")

    # We can easily make an efficiency plot as well
    # Let's compute the TaNC eff w.r.t Leading Piont for taus with Pt > 5 
    # with eta < 2.5
    denom_selection = selectedPatTaus.expr('$pt > 5') & \
            selectedPatTaus.expr('$eta > -2.5 && $eta < +2.5') & \
            selectedPatTaus.expr('$byLeadPionPtCut > 0.5')

    numerator_selection = denom_selection & selectedPatTaus.expr('$byTaNCfrOnePercent')

    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    bkg_histo, efficiency = plot.efficiency(
        events,
        expression=selectedPatTaus.expr('$pt'),
        numerator=numerator_selection,
        denominator=denom_selection,
        binning=(10, 0, 50))

    bkg_histo.Draw()
    efficiency.Draw("s, p")

    canvas.SaveAs("pt_eff.png")







