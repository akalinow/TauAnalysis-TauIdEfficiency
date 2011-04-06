makeVertexMultiplicityPlots()
{
  Double_t binEdges[22] = { -0.5,  0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5, 9.5, 
			    10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5 };

  TH1* histoData = new TH1F("histoData", "histoData", 21, binEdges);
  TH1* histoMC   = new TH1F("histoMC",   "histoMC",   21, binEdges);

  TString cacheFilePathData = "/data2/friis/tau_fakerate_ntuples/user/v/veelken/TauIdCommissioning/WplusJets/v5_0b/";

  TString prefix = "double_ntupleProducer_tauIdEffNtuple";
  TString suffix = "_prodCommissioningWplusJetsEnrichedNtuple.obj";

  TString branchNumVerticesPtGt5 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt5").Append(suffix);
  TString branchNumVerticesPtGt10 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt10").Append(suffix);
  TString branchNumVerticesPtGt15 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt15").Append(suffix);
  TString branchNumVerticesPtGt20 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt20").Append(suffix);
  
  TString branchGenPileUp = TString(prefix).Append("#addPileupInfo#numPileUpInteractions").Append(suffix);

  TChain* chainData = new TChain("Events");
  chainData->Add(TString(cacheFilePathData).Append("data_Mu_Run2011A_PromptReco/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainData->Draw(TString(branchNumVerticesPtGt10).Append(">>+histoData"));

  std::cout << histoData->GetName() << " created." << std::endl;
  if ( !histoData->GetSumw2N()     ) histoData->Sumw2();
  if (  histoData->Integral() > 0. ) histoData->Scale(1./histoData->Integral());

  delete chainData;

  TString cacheFilePathMC = "/data2/friis/tau_fakerate_ntuples/user/v/veelken/TauIdCommissioning/WplusJets/v5_0b/";

  TChain* chainMC = new TChain("Events");
  chainMC->Add(TString(cacheFilePathMC).Append("WplusJets/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainMC->Draw(TString(branchNumVerticesPtGt10).Append(">>+histoMC"));

  std::cout << histoMC->GetName() << " created." << std::endl;
  if ( !histoMC->GetSumw2N()     ) histoMC->Sumw2();
  if (  histoMC->Integral() > 0. ) histoMC->Scale(1./histoMC->Integral());

  TH2* histoPileUpCorrelationNumVerticesPtGt5  = new TH2F("histoPileUpCorrelationNumVerticesPtGt5", 
							  "histoPileUpCorrelationNumVerticesPtGt5",  10, -0.5, +9.5, 10, -0.5, +9.5);
  TH2* histoPileUpCorrelationNumVerticesPtGt10 = new TH2F("histoPileUpCorrelationNumVerticesPtGt10", 
							  "histoPileUpCorrelationNumVerticesPtGt10", 10, -0.5, +9.5, 10, -0.5, +9.5);
  TH2* histoPileUpCorrelationNumVerticesPtGt15 = new TH2F("histoPileUpCorrelationNumVerticesPtGt15", 
							  "histoPileUpCorrelationNumVerticesPtGt15", 10, -0.5, +9.5, 10, -0.5, +9.5);
  TH2* histoPileUpCorrelationNumVerticesPtGt20 = new TH2F("histoPileUpCorrelationNumVerticesPtGt20", 
							  "histoPileUpCorrelationNumVerticesPtGt20", 10, -0.5, +9.5, 10, -0.5, +9.5);

  chainMC->Draw(TString(branchNumVerticesPtGt5).Append(":").Append(branchGenPileUp).Append(">>+histoPileUpCorrelationNumVerticesPtGt5"));
  chainMC->Draw(TString(branchNumVerticesPtGt10).Append(":").Append(branchGenPileUp).Append(">>+histoPileUpCorrelationNumVerticesPtGt10"));
  chainMC->Draw(TString(branchNumVerticesPtGt15).Append(":").Append(branchGenPileUp).Append(">>+histoPileUpCorrelationNumVerticesPtGt15"));
  chainMC->Draw(TString(branchNumVerticesPtGt20).Append(":").Append(branchGenPileUp).Append(">>+histoPileUpCorrelationNumVerticesPtGt20"));
  
  delete chainMC;

  TH1* histoReweight = (TH1*)histoData->Clone("histoReweight");
  histoReweight->Divide(histoData, histoMC);

  TFile* outputFile = new TFile("vertexMultiplicityPlots.root", "RECREATE");
  histoData->Write();
  histoMC->Write();
  histoReweight->Write();
  histoPileUpCorrelationNumVerticesPtGt5->Write();
  histoPileUpCorrelationNumVerticesPtGt10->Write();
  histoPileUpCorrelationNumVerticesPtGt15->Write();
  histoPileUpCorrelationNumVerticesPtGt20->Write();
  delete outputFile;

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  histoMC->SetTitle("Num. Vertices #Sigma trackPt > 10 GeV");
  histoMC->SetStats(false);
  //histoMC->SetMaximum(1.);
  histoMC->SetMaximum(0.4);
  histoMC->SetXTitle("Num. Vertices");
  histoMC->SetYTitle("a.u.");
  histoMC->SetMarkerStyle(24);
  histoMC->SetMarkerColor(kRed);
  histoMC->SetLineColor(kRed);
  histoMC->Draw("e1p");

  histoData->SetMarkerStyle(20);
  histoData->SetMarkerColor(kBlack);
  histoData->SetLineColor(kBlack);
  histoData->Draw("e1psame");

  TLegend legend(0.64, 0.69, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);

  legend.AddEntry(histoMC, "Spring'11 MC", "p");
  legend.AddEntry(histoData, "2011 Data", "p");
  legend.Draw();

  canvas->Update();
  canvas->Print("vertexMultiplicity.pdf");

  canvas->Clear();

  histoReweight->SetTitle("Vertex Multiplicity Reweight factors");
  histoReweight->SetStats(false);
  histoReweight->SetXTitle("Num. Vertices");
  histoReweight->SetYTitle("Reweight Factor");
  histoReweight->SetMarkerStyle(20);
  histoReweight->SetMarkerColor(kBlack);
  histoReweight->SetLineColor(kBlack);
  histoReweight->Draw("e1p");

  canvas->Update();
  canvas->Print("vertexMultiplicityReweights.pdf");
}
