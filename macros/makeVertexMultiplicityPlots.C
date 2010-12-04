makeVertexMultiplicityPlots()
{
  Double_t binEdges[9] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5 };

  TH1* histoData = new TH1F("histoData", "histoData", 8, binEdges);
  TH1* histoMC   = new TH1F("histoMC",   "histoMC",   8, binEdges);

  TString cacheFilePathData = "/data1/friis/tau_fakerate_ntuples/user/v/veelken/TauIdCommissioning/WplusJets/v3_3/";

  TString prefix = "double_ntupleProducer_tauIdEffNtuple";
  TString suffix = "_prodCommissioningWplusJetsEnrichedNtuple.obj";

  TString branchNumVerticesPtGt5 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt5").Append(suffix);
  TString branchNumVerticesPtGt10 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt10").Append(suffix);
  TString branchNumVerticesPtGt15 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt15").Append(suffix);
  TString branchNumVerticesPtGt20 = TString(prefix).Append("#offlinePrimaryVertices#numVerticesPtGt20").Append(suffix);
  
  TString branchGenPileUp = TString(prefix).Append("#addPileupInfo#genPtHat").Append(suffix);

  TChain* chainData = new TChain("Events");
  chainData->Add(TString(cacheFilePathData).Append("runs132440_145761/tauIdEffEDNtuple_wPlusJetsEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs145762_147116/tauIdEffEDNtuple_wPlusJetsEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs147117_149442/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainData->Draw(TString(branchNumVerticesPtGt10).Append(">>+histoData"));

  std::cout << histoData->GetName() << " created." << std::endl;
  if ( !histoData->GetSumw2N()     ) histoData->Sumw2();
  if (  histoData->Integral() > 0. ) histoData->Scale(1./histoData->Integral());

  delete chainData;

  TString cacheFilePathMC = "/data1/friis/tau_fakerate_ntuples/user/v/veelken/TauIdCommissioning/WplusJets/v3_3/";

  TChain* chainMC = new TChain("Events");
  chainMC->Add(TString(cacheFilePathMC).Append("mcWplusJetsPU156bx/tauIdEffEDNtuple_wPlusJetsEnriched*"));

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
  histoMC->SetMaximum(1.);
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

  legend.AddEntry(histoMC, "MC 156bxPU", "p");
  legend.AddEntry(histoData, "Data", "p");
  legend.Draw();

  canvas->Update();
  canvas->Print("vertexMultiplicity.png");

  canvas->Clear();

  histoReweight->SetTitle("Vertex Multiplicity Reweight factors");
  histoReweight->SetStats(false);
  histoReweight->SetMarkerStyle(20);
  histoReweight->SetMarkerColor(kBlack);
  histoReweight->SetLineColor(kBlack);
  histoReweight->Draw("e1p");

  canvas->Update();
  canvas->Print("vertexMultiplicityReweights.png");
}
