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
/*
  TString cacheFilePathData = "/tmp/tau_fakerate_cache/user/v/veelken/TauIdCommissioning/qcdMuEnriched/v3_2/";

  TChain* chainData = new TChain("Events");
  chainData->Add(TString(cacheFilePathData).Append("runs132440_145761/tauIdEffEDNtuple_qcdMuEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs145762_147116/tauIdEffEDNtuple_qcdMuEnriched*"));
  chainData->Add(TString(cacheFilePathData).Append("runs147117_149442/tauIdEffEDNtuple_qcdMuEnriched*"));
  stc::cout << "chainData: entries = " << chainData->GetEntries() << std::endl;

  chainData->Draw("double_ntupleProducer_tauIdEffNtuple#offlinePrimaryVertices#numVerticesPtGt10_prodCommissioningQDCmuEnrichedNtuple.obj");
 */
  //TH1* histoData = (TH1*)gROOT->FindObject("htemp")->Clone("histoData");
  std::cout << histoData->GetName() << " created." << std::endl;
  if ( !histoData->GetSumw2N() ) histoData->Sumw2();
  histoData->Scale(1./histoData->Integral());

  delete chainData;
/*
  TString cacheFilePathMC = "/tmp/tau_fakerate_cache/user/v/veelken/TauIdCommissioning/qcdDiJet/v3_0/";

  TChain* chainMC = new TChain("Events");
  chainMC->Add(TString(cacheFilePathMC).Append("mcDYttPU156bx/tauIdEffEDNtuple_qcdDiJet*"));
  stc::cout << "chainMC: entries = " << chainMC->GetEntries() << std::endl;

  chainMC->Draw("double_ntupleProducer_tauIdEffNtuple#offlinePrimaryVertices#numVerticesPtGt10_prodCommissioningQDCdiJetNtuple.obj");
 */
  TString cacheFilePathMC = "/tmp/tau_fakerate_cache/user/v/veelken/TauIdCommissioning/WplusJets/v3_2/";

  TChain* chainMC = new TChain("Events");
  chainMC->Add(TString(cacheFilePathMC).Append("mcWplusJetsPU156bx/tauIdEffEDNtuple_wPlusJetsEnriched*"));

  chainMC->Draw("double_ntupleProducer_tauIdEffNtuple#offlinePrimaryVertices#numVerticesPtGt10_prodCommissioningWplusJetsEnrichedNtuple.obj>>+histoMC");

  //TH1* histoMC = (TH1*)gROOT->FindObject("htemp")->Clone("histoMC");
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
}
