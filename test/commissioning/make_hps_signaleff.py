#!/usr/bin/env python

'''
Example of use of TauNtuple software

Author: Evan K. Friis, UC Davis
'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager import TauNtupleManager
from TauAnalysis.TauIdEfficiency.ntauples.plotting import *

# Load FWLite libraries (prevents warnings)
from PhysicsTools.PythonAnalysis import *
from ROOT import *
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()
gROOT.SetBatch(True)
gStyle.SetOptStat(0)

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    #file = ROOT.TFile("commissioning_ntuple.root", "READ")
    #file = ROOT.TFile("/tmp/anayak/my_tauIdEff_ntuple.root", "READ")
    #multiple files
    events = TChain("Events")
    events.Add("/scratch/swanson/TauIDPAS_Data/v2_1_2/tauId*.root")
    #events.Add("/tmp/anayak/Z2TauTau/tauIdEffEDNtuple_qcdDiJet_*.root")
        
    # Get the events tree (this can also be a TChain)
    #events = file.Get("Events")

    # Get the ntuple we produced
    manager = TauNtupleManager(events, "tauIdEffNtuple")

    # Get the list of collections availabe for our ntuple
    print manager

    # So lets use selectedPatTaus
    selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeHPS")

    genTaus = manager.get_ntuple("tauGenJets")
    

    # The following operators are available
    # + - * / & | < > <= >=

    gStyle.SetErrorX(0)
    
    def plotPt(histo1, histo2, histo3, histo4, histo5, histo6):
        bool_=True
        xaxis="pt"
        if bool_:
            histo1.GetYaxis().SetRangeUser(5E-3,1.3)
            histo1.SetTitle("")
            histo1.GetYaxis().SetTitle("Efficiency")
            if xaxis == "eta":
                histo1.GetXaxis().SetTitle("gen Jet #eta")
            if xaxis == "pt":
                histo1.GetXaxis().SetTitle("gen Jet p_{T} [GeV/c]")
            histo1.Draw()
        histo2.SetMarkerSize(1.5)
        histo2.SetMarkerStyle(20)
        histo2.SetMarkerColor(1)
        histo3.SetMarkerStyle(21)
        histo3.SetMarkerColor(2)
        histo3.SetMarkerSize(1.5)
        histo4.SetMarkerStyle(22)
        histo4.SetMarkerColor(4)
        histo4.SetMarkerSize(1.5)
        histo5.SetMarkerColor(6)
        histo5.SetMarkerStyle(23)
        histo5.SetMarkerSize(1.5)
        histo6.SetMarkerColor(7)
        histo6.SetMarkerStyle(24)
        histo6.SetMarkerSize(1.5)
        histo2.Draw("pZsame")
        histo3.Draw("pZsame")
        histo4.Draw("pZsame")
        histo5.Draw("pZsame")
        histo6.Draw("pZsame")
        
    #Make efficiency plots for Pt

    denom_selection = genTaus.expr('$genPt > 10 && abs($genEta) < 2.5 && $genDecayMode > 1.5')
    
    #base_selection = selectedPatTaus.expr('$genPt > 5 && abs($genEta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5')
    base_selection = selectedPatTaus.expr('$pt > 10 && abs($eta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5')
    my_expression = genTaus.expr('$genPt')
    myPtBins=(50,0.,100.)

    bkg_dummy, my_eff_dummy = efficiency(
        events,
        expression=my_expression,
        numerator=genTaus.expr('$genPt > 5000 && abs($genEta) < 2.5'),
        denominator=denom_selection,
        output_name= "dummy",
        binning=(50, 0., 100.)
        )
    
    # plot distribution
    pt_denom = draw(
        events,
        expression=my_expression,
        selection=denom_selection,
        output_name="pt_denom",
        binning=myPtBins
        )

    pt_num_match = draw(
        events,
        expression=selectedPatTaus.expr('$genPt'),
        selection=base_selection,
        output_name="pt_num_match",
        binning=myPtBins
        )

    my_eff_match = ROOT.TGraphAsymmErrors(pt_num_match, pt_denom)

    pt_num_leadTrack = draw(
        events,
        expression=selectedPatTaus.expr('$genPt'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding'),
        output_name="pt_num_leadTrack",
        binning=myPtBins
        )

    my_eff_leadTrack = ROOT.TGraphAsymmErrors(pt_num_leadTrack, pt_denom)

    pt_num_looseIso = draw(
        events,
        expression=selectedPatTaus.expr('$genPt'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationLoose'),
        output_name="pt_num_looseIso",
        binning=myPtBins
        )

    my_eff_looseIso = ROOT.TGraphAsymmErrors(pt_num_looseIso, pt_denom)

    pt_num_mediumIso = draw(
        events,
        expression=selectedPatTaus.expr('$genPt'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationMedium'),
        output_name="pt_num_mediumIso",
        binning=myPtBins
        )
    my_eff_mediumIso = ROOT.TGraphAsymmErrors(pt_num_mediumIso, pt_denom)

    pt_num_tightIso = draw(
        events,
        expression=selectedPatTaus.expr('$genPt'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationTight'),
        output_name="pt_num_tightIso",
        binning=myPtBins
        )

    my_eff_tightIso = ROOT.TGraphAsymmErrors(pt_num_tightIso, pt_denom)
    
    c1 = TCanvas("c1","",0,0,500,500)
    plotPt(bkg_dummy,my_eff_match,my_eff_leadTrack,my_eff_looseIso,my_eff_mediumIso,my_eff_tightIso)
    c1.SetGridy(1)
    gStyle.SetErrorX(0);
    
    simpleStyle = TStyle("simpleStyle","");
    simpleStyle.SetCanvasColor(0);
    simpleStyle.SetPadColor(0);
    simpleStyle.SetFrameFillColor(0);
    simpleStyle.SetStatColor(0);
    simpleStyle.SetOptStat(0);
    simpleStyle.SetTitleFillColor(0);
    simpleStyle.SetCanvasBorderMode(0);
    simpleStyle.SetPadBorderMode(0);
    simpleStyle.SetFrameBorderMode(0);
    simpleStyle.SetLegendBorderSize(1);
    simpleStyle.SetErrorX(0);
    simpleStyle.cd();
    
    t1 = TLegend(0.64, 0.74, 0.89, 0.89, "","brNDC")
    t1.AddEntry(my_eff_match,"Matched","P");
    t1.AddEntry(my_eff_leadTrack,"Decay Mode Finding","P");
    t1.AddEntry(my_eff_looseIso,"Loose Isolation","P");
    t1.AddEntry(my_eff_mediumIso,"Medium Isolation","P");
    t1.AddEntry(my_eff_tightIso,"Tight Isolation","P");
    t1.SetFillColor(0);
    t1.Draw();
    
    t2 = ROOT.TPaveText(0.12, 0.87, 0.45, 0.92, "NDC")
    t2.AddText("CMS Preliminary 7 TeV")
    t2.SetTextAlign(13)
    t2.SetTextSize(0.04)
    t2.SetFillStyle(0)
    t2.SetBorderSize(0)
    t2.Draw();
    
    t3 = ROOT.TPaveText(0.12, 0.8, 0.45, 0.85, "NDC")
    t3.AddText("P_{T} > 10 GeV/c, |#eta| < 2.5")
    t3.SetTextAlign(13)
    t3.SetTextSize(0.04)
    t3.SetFillStyle(0)
    t3.SetBorderSize(0)
    t3.Draw()

    
    c1.SaveAs("hps_pt_eff.png")
    c1.SaveAs("hps_pt_eff.pdf")
    
    
    gStyle.SetErrorX(0)
    
    def plotEta(histo1, histo2, histo3, histo4, histo5, histo6):
        bool_=True
        xaxis="eta"
        if bool_:
            histo1.GetYaxis().SetRangeUser(5E-3,1.3)
            histo1.SetTitle("")
            histo1.GetYaxis().SetTitle("Efficiency")
            if xaxis == "eta":
                histo1.GetXaxis().SetTitle("gen Jet #eta")
            if xaxis == "pt":
                histo1.GetXaxis().SetTitle("gen Jet p_{T} [GeV/c]")
            histo1.Draw()
        histo2.SetMarkerSize(1.5)
        histo2.SetMarkerStyle(20)
        histo2.SetMarkerColor(1)
        histo3.SetMarkerStyle(21)
        histo3.SetMarkerColor(2)
        histo3.SetMarkerSize(1.5)
        histo4.SetMarkerStyle(22)
        histo4.SetMarkerColor(4)
        histo4.SetMarkerSize(1.5)
        histo5.SetMarkerColor(6)
        histo5.SetMarkerStyle(23)
        histo5.SetMarkerSize(1.5)
        histo6.SetMarkerColor(7)
        histo6.SetMarkerStyle(24)
        histo6.SetMarkerSize(1.5)
        histo2.Draw("pZsame")
        histo3.Draw("pZsame")
        histo4.Draw("pZsame")
        histo5.Draw("pZsame")
        histo6.Draw("pZsame")
        
    #Make efficiency plots for Eta
    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    
    denom_selection = genTaus.expr('$genPt > 10 && abs($genEta) < 2.5 && $genDecayMode > 1.5')
    
    base_selection = selectedPatTaus.expr('$pt > 15 && abs($eta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5')
    my_expression = genTaus.expr('$genEta')
    myEtaBins=(50,-2.5,2.5)

    bkg_dummy, my_eff_dummy = efficiency(
        events,
        expression=my_expression,
        numerator=genTaus.expr('$genPt > 5000 && abs($genEta) < 2.5'),
        denominator=denom_selection,
        output_name= "dummy",
        binning=myEtaBins
        )
    
    # plot distribution
    pt_denom = draw(
        events,
        expression=my_expression,
        selection=denom_selection,
        output_name="eta_denom",
        binning=myEtaBins
        )

    pt_num_match = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection,
        output_name="eta_num_match",
        binning=myEtaBins
        )

    my_eff_match = ROOT.TGraphAsymmErrors(eta_num_match, eta_denom)

    eta_num_leadTrack = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding'),
        output_name="eta_num_leadTrack",
        binning=myEtaBins
        )

    my_eff_leadTrack = ROOT.TGraphAsymmErrors(eta_num_leadTrack, eta_denom)

    eta_num_looseIso = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationLoose'),
        output_name="eta_num_looseIso",
        binning=myEtaBins
        )

    my_eff_looseIso = ROOT.TGraphAsymmErrors(eta_num_looseIso, eta_denom)

    eta_num_mediumIso = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationMedium'),
        output_name="eta_num_mediumIso",
        binning=myEtaBins
        )
    my_eff_mediumIso = ROOT.TGraphAsymmErrors(eta_num_mediumIso, eta_denom)

    eta_num_tightIso = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byIsolationTight'),
        output_name="eta_num_tightIso",
        binning=myEtaBins
        )

    my_eff_tightIso = ROOT.TGraphAsymmErrors(eta_num_tightIso, eta_denom)
    
    c2 = TCanvas("c1","",0,0,500,500)
    plotEta(bkg_dummy,my_eff_match,my_eff_leadTrack,my_eff_looseIso,my_eff_mediumIso,my_eff_tightIso)
    c2.SetGridy(1)
    gStyle.SetErrorX(0);
    
    simpleStyle.cd();
    t1.Draw();
    t2.Draw();
    t3.Draw()
    
    c2.SaveAs("hps_eta_eff.png")
    c2.SaveAs("hps_eta_eff.pdf")
