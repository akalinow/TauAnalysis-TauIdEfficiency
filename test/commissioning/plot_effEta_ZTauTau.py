#!/usr/bin/env python

'''
Example of use of TauNtuple software

Author: Evan K. Friis, UC Davis
'''

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager import TauNtupleManager
#import TauAnalysis.TauIdEfficiency.ntauples.plotting as plot
from TauAnalysis.TauIdEfficiency.ntauples.plotting import *

# Load FWLite libraries (prevents warnings)
from PhysicsTools.PythonAnalysis import *
from ROOT import *
gSystem.Load("libFWCoreFWLite.so")
AutoLibraryLoader.enable()
gROOT.SetBatch(True)
#gStyle.SetOptStat(0)

#gPad.SetLogy(True)

#MC = False

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    #file = ROOT.TFile("commissioning_ntuple.root", "READ")
    #file = ROOT.TFile("/tmp/anayak/my_tauIdEff_ntuple.root", "READ")
    #file = ROOT.TFile("rfio:/castor/cern.ch/user/a/abdollah/TauIdCommissioning/qcdDiJet/v1_1/tauIdEff_ntuple_10_2.root", "READ")
    # Get the events tree (this can also be a TChain)
    #events = file.Get("Events")

    #multiple files
    events = TChain("Events")
    events.Add("rfio:/castor/cern.ch/user/a/abdollah/TauIdCommissioning/qcdDiJet/v2_1_2/tauIdEffEDNtuple_qcdDiJet_*.root")
    #events.Add("/tmp/anayak/Z2TauTau/tauIdEffEDNtuple_qcdDiJet_*.root")
    
    # Get the ntuple we produced
    manager = TauNtupleManager(events, "tauIdEffNtuple")

    # Get the list of collections availabe for our ntuple
    print manager

    # So lets use selectedPatTaus
    #selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone")
    selectedPatTaus = manager.get_ntuple("patPFTausDijetTagAndProbeFixedCone")

    # All this helper funciton does it makes it easy for use
    # to interface to TTree::Draw, etc.

    # Get list of variables available
    print selectedPatTaus
    # returns
    # Tau Ntuple - collection selectedPatTaus
    #  byIsolation
    #  byLeadPionPt
    #  byTaNCfrOne
    #  absEta
    #  pt

    genTaus = manager.get_ntuple("tauGenJets")
    print genTaus
    
    # Now it easy to create expressions using the expr method
    print selectedPatTaus.expr('$pt')
    # returns
    # exampleNtuple#selectedPatTaus#pt

    # The following operators are available
    # + - * / & | < > <= >=

    gStyle.SetErrorX(0)
    
    def plot(histo1, histo2, histo3, histo4, histo5, histo6, bool_, xaxis):
        if bool_:
            histo1.GetYaxis().SetRangeUser(5E-3,1.2)
            histo1.SetTitle("")
            histo1.GetYaxis().SetTitle("Efficiency")
            if xaxis == "eta":
                #histo1.GetXaxis().SetTitle("Jet #eta ")
                histo1.GetXaxis().SetTitle("gen Jet #eta ")
            if xaxis == "pt":
                histo1.GetXaxis().SetTitle("gen Jet p_{T} (GeV/c)")
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
        histo5.SetMarkerStyle(23)
        histo5.SetMarkerColor(6)
        histo5.SetMarkerSize(1.5)
        histo6.SetMarkerStyle(24)
        histo6.SetMarkerColor(7)
        histo6.SetMarkerSize(1.5)
        gStyle.SetErrorX(0)
        histo2.Draw("pZsame")
        histo3.Draw("pZsame")
        histo4.Draw("pZsame")
        histo5.Draw("pZsame")
        histo6.Draw("pZsame")

    #Make efficiency plots
    # The efficiency function returns a tuple with
    # a histo background + a TGraph asymmerrors
    
    denom_selection = genTaus.expr('$genPt > 5 && abs($genEta) < 2.5 && $genDecayMode > 1.5')
    #base_selection = selectedPatTaus.expr('$genPt > 5 && abs($genEta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5')
    base_selection = selectedPatTaus.expr('$jetPt > 15 && abs($jetEta) < 2.5 && $genMatch > 0.5 && $genDecayMode > 1.5')
    my_expression = genTaus.expr('$genEta')
    myEtaBins=(60,-3,3)
    
    bkg_dummy, my_eff_dummy = efficiency(
        events,
        expression=my_expression,
        numerator=genTaus.expr('$genPt > 5000 && abs($genEta) < 2.5'),
        denominator=denom_selection,
        output_name= "match",
        binning=myEtaBins
        )
    
    # plot distribution
    eta_denom = draw(
        events,
        expression=my_expression,
        selection=denom_selection,
        output_name="eta_denom",
        binning=myEtaBins
        )

    eta_num_match = draw(
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

    eta_num_leadTrackPt = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byLeadTrackPtCut'),
        output_name="eta_num_leadTrackPt",
        binning=myEtaBins
        )
    my_eff_leadTrackPt = ROOT.TGraphAsymmErrors(eta_num_leadTrackPt, eta_denom)

    eta_num_TrackIso = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byLeadTrackPtCut') &\
        selectedPatTaus.expr('$byTrackIsolation'),
        output_name="eta_num_TrackIso",
        binning=myEtaBins
        )
    my_eff_TrackIso = ROOT.TGraphAsymmErrors(eta_num_TrackIso, eta_denom)

    eta_num_EcalIso = draw(
        events,
        expression=selectedPatTaus.expr('$genEta'),
        selection=base_selection & selectedPatTaus.expr('$byLeadTrackFinding && $byLeadTrackPtCut') &\
        selectedPatTaus.expr('$byTrackIsolation && $byEcalIsolation'),
        output_name="eta_num_EcalIso",
        binning=myEtaBins
        )
    my_eff_EcalIso = ROOT.TGraphAsymmErrors(eta_num_EcalIso, eta_denom)
    
    c1 = TCanvas("c1","",0,0,800,600)
    plot(bkg_dummy,my_eff_match,my_eff_leadTrack,my_eff_leadTrackPt,my_eff_TrackIso,my_eff_EcalIso,True,"eta")
    #c1.SetLogy(1)
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
    simpleStyle.SetErrorX(0);
    simpleStyle.cd();
    
    t1 = TLegend(0.6, 0.10, 0.7, 0.35, "","brNDC")
    t1.SetTextSize(0.04)
    t1.AddEntry(my_eff_match,"Matched","P");
    t1.AddEntry(my_eff_leadTrack,"Lead Track Finding","P");
    t1.AddEntry(my_eff_leadTrackPt,"Lead Track p_{T} cut","P");
    t1.AddEntry(my_eff_TrackIso,"Track Isolation","P");
    t1.AddEntry(my_eff_EcalIso,"Gamma Isolation","P");
    t1.SetFillColor(0);
    t1.SetLineColor(0);
    t1.SetBorderSize(0);
    t1.Draw();
    
    t2= TLatex(0.20,0.85,"CMS Preliminary 2010 7 TeV");
    t2.SetNDC();
    t2.SetTextSize(0.04);
    t2.Draw();
    c1.SaveAs("Eff_Z2TT_Iso_Eta_genMatch.png")
    
