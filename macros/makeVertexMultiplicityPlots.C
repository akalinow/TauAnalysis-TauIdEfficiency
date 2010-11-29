makeVertexMultiplicityPlots()
{
  Double_t binEdges[9] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5 };

  TH1* histoData = new TH1F("histoData", "histoData", 8, binEdges);
  TH1* histoMC   = new TH1F("histoMC",   "histoMC",   8, binEdges);

  TString cacheFilePathData = "/tmp/tau_fakerate_cache/user/v/veelken/TauIdCommissioning/WplusJets/v3_2/";

  TChain* chainData = new TChain("Events");
  chainData->Add(TString(cacheFilePathData).Append("runs132440_145761/tauIdEffEDNtuple_wPlusJetsEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs145762_147116/tauIdEffEDNtuple_wPlusJetsEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs147117_149442/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainData->Draw("double_ntupleProducer_tauIdEffNtuple#offlinePrimaryVertices#numVerticesPtGt10_prodCommissioningWplusJetsEnrichedNtuple.obj>>+histoData");

  std::cout << histoData->GetName() << " created." << std::endl;
  if ( !histoData->GetSumw2N() ) histoData->Sumw2();
  histoData->Scale(1./histoData->Integral());

  delete chainData;

  TString cacheFilePathMC = "/tmp/tau_fakerate_cache/user/v/veelken/TauIdCommissioning/WplusJets/v3_2/";

  TChain* chainMC = new TChain("Events");
  chainMC->Add(TString(cacheFilePathMC).Append("mcWplusJetsPU156bx/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainMC->Draw("double_ntupleProducer_tauIdEffNtuple#offlinePrimaryVertices#numVerticesPtGt10_prodCommissioningWplusJetsEnrichedNtuple.obj>>+histoMC");

  std::cout << histoMC->GetName() << " created." << std::endl;
  if ( !histoMC->GetSumw2N() ) histoMC->Sumw2();
  histoMC->Scale(1./histoMC->Integral());

  delete chainMC;

  TH1* histoReweight = (TH1*)histoData->Clone("histoReweight");
  histoReweight->Divide(histoData, histoMC);

  TFile* outputFile = new TFile("vertexMultiplicityPlots.root", "RECREATE");
  histoData->Write();
  histoMC->Write();
  histoReweight->Write();
  delete outputFile;

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
/*
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
 */
  histoReweight->SetTitle("Vertex Multiplicity Reweight factors");
  histoReweight->SetStats(false);
  histoReweight->SetMarkerStyle(20);
  histoReweight->SetMarkerColor(kBlack);
  histoReweight->SetLineColor(kBlack);
  histoReweight->Draw("e1p");

  canvas->Update();
  canvas->Print("vertexMultiplicityReweights.png");
}
