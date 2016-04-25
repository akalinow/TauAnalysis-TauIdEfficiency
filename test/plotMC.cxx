/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::string fName = "TnP_MuonToTau_MisID_MC.root";
  std::string topDirectory = "tpTree/";
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFitCanvas(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;
  std::string binnedVars = "__byLooseCombinedIsolationDeltaBetaCorr3Hits_bin0__decayModeFinding_bin0__pair_deltaR_bin0__pair_dz_bin0__pair_probeMultiplicity_bin0__vpvPlusExpo/";

  std::string aPattern = "pair_deltaR_bin0";
  if(category.find("mcTrue")!=std::string::npos) binnedVars.replace(binnedVars.find(aPattern),
								    aPattern.size(),
								    "mcTrue_bin0__"+aPattern);

  aPattern = "abseta";
  category.replace(category.find(aPattern),aPattern.size(),"abseta0");
  
  dirName = topDirectory+category+"/abseta_bin0"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;

  TCanvas *fit_canvasEta0 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta0->SetName("fit_canvasEta0");
  fit_canvasEta0->Draw();
  fit_canvasEta0->Print(("./fig_png/"+category+"_fit_canvasEta0_MC.png").c_str());

  TPad *aPad = (TPad*)fit_canvasEta0->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing0_MC.png").c_str());

  aPattern = "abseta0";
  category.replace(category.find(aPattern),aPattern.size(),"abseta12");
  
  dirName = topDirectory+category+"/abseta_bin0"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;
  TCanvas *fit_canvasEta1 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta1->SetName("fit_canvasEta1");
  fit_canvasEta1->Draw();
  fit_canvasEta1->Print(("fig_png/"+category+"_fit_canvasEta1_MC.png").c_str());

  aPad = (TPad*)fit_canvasEta1->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing1_MC.png").c_str());

  dirName = topDirectory+category+"/abseta_bin1"+binnedVars;
  objName = "fit_canvas";
  objectPath = dirName+objName;
  TCanvas *fit_canvasEta2 = (TCanvas*)file.Get(objectPath.c_str());
  fit_canvasEta2->SetName("fit_canvasEta2");
  fit_canvasEta2->Draw();
  fit_canvasEta2->Print(("fig_png/"+category+"_fit_canvasEta2_MC.png").c_str());

  aPad = (TPad*)fit_canvasEta2->FindObject("fit_canvas_1");
  aPad->SetCanvasSize(1000,1000);
  aPad->Print(("./fig_png/"+category+"_fit_passing2_MC.png").c_str());
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphAsymmErrors *getMCResult(std::string category = "againstMuonLoose3_pt_abseta"){

  TFile file(fName.c_str());

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

  TGraphAsymmErrors *aGraph  = getMCResult(category);

  category = category+"_mcTrue";
  TGraphAsymmErrors *aGraphMCTrue  = getMCResult(category);
  
  aGraphMCTrue->SetName("aGraphMCTrue");
  aGraphMCTrue->SetLineColor(2);
  aGraphMCTrue->SetMarkerColor(2);

  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(3E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(5E-6);
    hFrame->SetMaximum(1.5E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("|#eta|");   

  TLegend *aLegend = new TLegend(0.15,0.2,0.3,0.4,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetFillStyle(4000);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(aGraph,"MC Z#rightarrow ll","lp");
  aLegend->AddEntry(aGraphMCTrue,"#splitline{MC Z#rightarrow #mu#mu with}{tag&probe matched to #mu}","lp");

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  aCanvas->SetHighLightColor(2);
  aCanvas->Range(0,0,1,1);
  aCanvas->SetFillColor(0);
  aCanvas->SetBorderMode(0);
  aCanvas->SetBorderSize(2);
  aCanvas->SetFrameBorderMode(0);
  
  hFrame->Draw();
  aGraph->Draw("p");
  aGraphMCTrue->Draw("p");
  aLegend->Draw();
  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedMC.png").c_str());
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotAll(){

  plotFitCanvas("againstMuonLoose3_pt_abseta");
  plotFitCanvas("againstMuonTight3_pt_abseta");

  plotFitCanvas("againstMuonLoose3_pt_abseta_mcTrue");
  plotFitCanvas("againstMuonTight3_pt_abseta_mcTrue");
  
  plotMistagRate("againstMuonLoose3_pt_abseta");
  plotMistagRate("againstMuonTight3_pt_abseta");
  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
