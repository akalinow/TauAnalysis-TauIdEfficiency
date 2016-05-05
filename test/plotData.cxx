/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::string fName = "TnP_MuonToTau_MisID_Data.root";
std::string fNameMC = "TnP_MuonToTau_MisID_MC.root";
std::string topDirectory = "tpTree/";
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFitCanvas(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;
  std::string binnedVars = "__byLooseCombinedIsolationDeltaBetaCorr3Hits_bin0__decayModeFinding_bin0__pair_deltaR_bin0__pair_dz_bin0__pair_probeMultiplicity_bin0__vpvPlusExpo/";
  binnedVars = "__byLooseCombinedIsolationDeltaBetaCorr3Hits_bin0__decayModeFinding_bin0__pair_deltaR_bin0__pair_dz_bin0__pair_probeMultiplicity_bin0__tag_triggerMatch_bin0__vpvPlusExpo/"; 
  
  std::string endPattern = "__vpvPlusExpo/";

  dirName = topDirectory+category+"/abseta_bin0"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;
  TCanvas *fit_canvasEta0 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta0->SetName("fit_canvasEta0");
  fit_canvasEta0->Draw();
  fit_canvasEta0->Print(("./fig_png/"+category+"_fit_canvasEta0_Data.png").c_str());

  TPad *aPad = (TPad*)fit_canvasEta0->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing0_Data.png").c_str());

  dirName = topDirectory+category+"/abseta_bin1"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;
  TCanvas *fit_canvasEta1 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta1->SetName("fit_canvasEta1");
  fit_canvasEta1->Draw();
  fit_canvasEta1->Print(("./fig_png/"+category+"_fit_canvasEta1_Data.png").c_str());

  aPad = (TPad*)fit_canvasEta1->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing1_Data.png").c_str());

  dirName = topDirectory+category+"/abseta_bin2"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;
  TCanvas *fit_canvasEta2 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta2->SetName("fit_canvasEta2");
  fit_canvasEta2->Draw();
  fit_canvasEta2->Print(("./fig_png/"+category+"_fit_canvasEta2_Data.png").c_str());

  aPad = (TPad*)fit_canvasEta2->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing2_Data.png").c_str());

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphAsymmErrors *getMCResult(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fNameMC.c_str());

  std::string aPattern = "abseta";
  category.replace(category.find(aPattern),aPattern.size(),"abseta0");

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;

  TCanvas *aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());
  TGraphAsymmErrors *hxy_fit_eff_eta0 = (TGraphAsymmErrors*)aEffCanvas->FindObject("hxy_fit_eff");
  hxy_fit_eff_eta0->SetName("hxy_fit_eff_eta0");

  aPattern = "abseta0";
  category.replace(category.find(aPattern),aPattern.size(),"abseta12");
  dirName = topDirectory+category+"/fit_eff_plots/";
  objName = "abseta_PLOT";
  objectPath = dirName+objName;
  
  aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());
  TGraphAsymmErrors *hxy_fit_eff_eta12 = (TGraphAsymmErrors*)aEffCanvas->FindObject("hxy_fit_eff");
  hxy_fit_eff_eta12->SetName("hxy_fit_eff_eta12");

  hxy_fit_eff_eta0->Print("all");

  Double_t x, y;
  hxy_fit_eff_eta0->GetPoint(0,x,y);

  hxy_fit_eff_eta12->Expand(3);
  hxy_fit_eff_eta12->SetPoint(2,x,y);
  hxy_fit_eff_eta12->SetPointError(2,
				   hxy_fit_eff_eta0->GetErrorXlow(0), hxy_fit_eff_eta0->GetErrorXhigh(0),
				   hxy_fit_eff_eta0->GetErrorYlow(0), hxy_fit_eff_eta0->GetErrorYhigh(0));
  
  hxy_fit_eff_eta12->Print("all");
  
  return hxy_fit_eff_eta12;
  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotMistagRate(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;

  TCanvas *aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());
  aEffCanvas->SetName("aEffCanvas");

  TGraphAsymmErrors *hxy_fit_eff = (TGraphAsymmErrors*)aEffCanvas->FindObject("hxy_fit_eff");
  TH1F *hFrame = (TH1F*)aEffCanvas->FindObject("frame");
  hFrame->SetXTitle("|#eta|");
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(3E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(5E-6);
    hFrame->SetMaximum(1.5E-3);
  }

  category = category+"_mcTrue";
  TGraphAsymmErrors *aGraphMCTrue  = getMCResult(category);
  aGraphMCTrue->SetName("aGraphMCTrue");
  aGraphMCTrue->SetLineColor(2);
  aGraphMCTrue->SetMarkerColor(2);
  
  TLegend *aLegend = new TLegend(0.12,0.2,0.3,0.35,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetFillStyle(4000);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(hxy_fit_eff,"Data","lp");
  aLegend->AddEntry(aGraphMCTrue,"#splitline{MC Z#rightarrow #mu#mu with}{tag&probe matched to #mu}","lp");
  
  aEffCanvas->Draw();
  aGraphMCTrue->Draw("p");
  hxy_fit_eff->Draw("p");
  aLegend->Draw();
  
  aEffCanvas->Print(("fig_png/"+category+"_misTagRateDATAvsMatchedMC.png").c_str());

  std::cout<<"DATA: "<<std::endl;
  hxy_fit_eff->Print("all");
  std::cout<<"MC: "<<std::endl;
  aGraphMCTrue->Print("all");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void test(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;

  TCanvas *aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());
  aEffCanvas->SetName("aEffCanvas");
  aEffCanvas->SetLogy();
  TGraphAsymmErrors *hxy_fit_eff = (TGraphAsymmErrors*)aEffCanvas->FindObject("hxy_fit_eff");

  hxy_fit_eff->Print("all");
  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotAll(){

  plotFitCanvas("againstMuonLoose3_pt_abseta");
  plotFitCanvas("againstMuonTight3_pt_abseta");

  plotMistagRate("againstMuonLoose3_pt_abseta");
  plotMistagRate("againstMuonTight3_pt_abseta");
  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
