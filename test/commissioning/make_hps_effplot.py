#!/usr/bin/env python

'''
Example of use of TauNtuple software

Author: Evan K. Friis, UC Davis
'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager import TauNtupleManager
import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    #file = ROOT.TFile("tauIdEff_ntuple.root", "READ")
    
    #Get the events tree (this can also be a TChain)
    chain = ROOT.TChain("Events")
    chain.Add("/scratch/swanson/TauIDPAS_Data/v2_1_2/tau*.root")
    events = chain
    
    #events = file.Get("Events")

    # Get the ntuple we produced
    manager = TauNtupleManager(events, "tauIdEffNtuple")

    # Get the list of collections availabe for our ntuple
    print manager

    # So lets use selectedPatTaus
    selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeHPS")

    # All this helper funciton does it makes it easy for use
    # to interface to TTree::Draw, etc.


    # The following operators are available
    # + - * / & | < > <= >=

    # There are also tools to make plots.  Here is an example to draw the Pt spectrum
    # for taus in the barrel 
    canvas = ROOT.TCanvas("example", "example", 500, 500)

 

    # We can easily make an efficiency plot as well
   
    denom_selection = selectedPatTaus.expr('$jetPt > 10') & selectedPatTaus.expr('abs($jetEta) < 2.5') & selectedPatTaus.expr('$genMatch>0.5') 

    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') #&\
    #        selectedPatTaus.expr('$againstElectron > 0.5')

    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    bkg_histo, efficiency1 = plot.efficiency(
        events,
        expression=selectedPatTaus.expr('$jetPt'),
        numerator=numerator_selection,
        denominator=denom_selection,
        binning=(25, 0, 50))
    
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationLoose > 0.5')
   
    bkg_histo, efficiency2 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetPt'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,0,50))
        
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationMedium > 0.5')
            
    bkg_histo, efficiency3 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetPt'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,0,50))
    	
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationTight > 0.5')
            
    bkg_histo, efficiency4 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetPt'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,0,50))

    bkg_histo.GetYaxis().SetRangeUser(0, 1.0)
    bkg_histo.SetTitle("HPS #tau Efficiency for Z->#tau+#tau MC")
    bkg_histo.GetXaxis().SetTitle("Jet P_{T} (GeV)")
    bkg_histo.Draw()
    ROOT.gPad.SetLogy(False)
    efficiency1.SetLineColor(2)
    efficiency1.SetMarkerStyle(20)
    efficiency1.SetMarkerColor(2)
    efficiency1.Draw("s, p")
    efficiency2.SetLineColor(4)
    efficiency2.SetMarkerStyle(20)
    efficiency2.SetMarkerColor(4)
    efficiency2.Draw("s, p")
    efficiency3.SetMarkerStyle(20)
    efficiency3.Draw("s, p")
    efficiency4.SetLineColor(8)
    efficiency4.SetMarkerStyle(20)
    efficiency4.SetMarkerColor(8)
    efficiency4.Draw("s, p")
    Legend = ROOT.TLegend(0.65,0.3,0.89,0.14,"")
    Legend.AddEntry(efficiency1,"Lead Track Finding","p")
    Legend.AddEntry(efficiency2,"Loose Isolation","p")
    Legend.AddEntry(efficiency3,"Medium Isolation","p")
    Legend.AddEntry(efficiency4,"Tight Isolation","p")
    Legend.SetFillColor(0)
    Legend.Draw()

   
    canvas.SaveAs("hps_pt_eff.pdf")


    # We can easily make an efficiency plot as well
 
    denom_selection = selectedPatTaus.expr('$jetPt > 5') & selectedPatTaus.expr('abs($jetEta) < 2.5') & selectedPatTaus.expr('$genMatch>0.5') 

    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') #&\
    #        selectedPatTaus.expr('$againstElectron > 0.5')

    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    bkg_histo, efficiency1 = plot.efficiency(
        events,
        expression=selectedPatTaus.expr('$jetEta'),
        numerator=numerator_selection,
        denominator=denom_selection,
        binning=(25, -2.5, 2.5))
    
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationLoose > 0.5')
   
    bkg_histo, efficiency2 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetEta'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,-2.5,2.5))
        
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationMedium > 0.5')
            
    bkg_histo, efficiency3 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetEta'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,-2.5,2.5))
    	
    numerator_selection = denom_selection & selectedPatTaus.expr('$byLeadTrackFinding > 0.5') & selectedPatTaus.expr('$byIsolationTight > 0.5')
            
    bkg_histo, efficiency4 = plot.efficiency(
    	events,
    	expression=selectedPatTaus.expr('$jetEta'),
    	numerator=numerator_selection,
    	denominator=denom_selection,
    	binning=(25,-2.5,2.5))

    bkg_histo.GetYaxis().SetRangeUser(0, 1.0)
    bkg_histo.SetTitle("HPS #tau Efficiency for Z->#tau+#tau MC")
    bkg_histo.GetXaxis().SetTitle("#eta")
    bkg_histo.Draw()
    ROOT.gPad.SetLogy(False)
    efficiency1.SetLineColor(2)
    efficiency1.SetMarkerStyle(20)
    efficiency1.SetMarkerColor(2)
    efficiency1.Draw("s, p")
    efficiency2.SetLineColor(4)
    efficiency2.SetMarkerStyle(20)
    efficiency2.SetMarkerColor(4)
    efficiency2.Draw("s, p")
    efficiency3.SetMarkerStyle(20)
    efficiency3.Draw("s, p")
    efficiency4.SetLineColor(8)
    efficiency4.SetMarkerStyle(20)
    efficiency4.SetMarkerColor(8)
    efficiency4.Draw("s, p")
    Legend.Draw()


    canvas.SaveAs("hps_eta_eff.pdf")
    
   