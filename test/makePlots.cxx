/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::string fNameMC =   "TnP_MuonToTau_MisID_MC.root";
std::string fNameData = "TnP_MuonToTau_MisID_Data.root";

std::string topDirectory = "tpTree/";
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFitCanvas(std::string category = "againstMuonLoose3_Zmumu", bool isData=false){

  std::string fName = fNameMC;
  if(isData) fName = fNameData;
  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "fit_canvas";
  std::string objectPath = dirName+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";
  fitModelName = "Zll_Model";
  std::string binnedVars =   "__mcTrue_bin0__";  
  if(isData){
    fitModelName = "Data_Model";
    binnedVars =   "__";
  }

  for(unsigned int iEta=0;iEta<3;++iEta){

    std::string etaBinNumber = std::to_string(iEta);
    std::string nameSuffix = "Eta"+etaBinNumber;
    if(isData) nameSuffix+="Data";
    dirName = topDirectory+category+"/abseta_bin"+etaBinNumber+binnedVars+fitModelName+"_Eta"+etaBinNumber+"/";
    objectPath = dirName+objName;
       
    TCanvas *fit_canvas = (TCanvas*)file.Get(objectPath.c_str());

    std::cout<<objectPath<<" fit_canvas: "<<fit_canvas<<std::endl;
    if(!fit_canvas) return;
    
    fit_canvas->Draw();
    fit_canvas->Print(("./fig_png/"+category+"_fit_canvas_"+nameSuffix+".png").c_str());
    TPad *aPad = (TPad*)fit_canvas->FindObject("fit_canvas_1");
    aPad->SetCanvasSize(1000,1000);
    aPad->Print(("./fig_png/"+category+"_fit_passing"+nameSuffix+".png").c_str());    
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphAsymmErrors *getResultGraph(std::string category = "againstMuonLoose3", bool count=false, bool isData=false){

  std::string fName = fNameMC;
  if(isData) fName = fNameData;
  TFile file(fName.c_str());

  std::string fitType = "fit";
  if(count) fitType = "cnt";

  std::string dirName = topDirectory+category+"/"+fitType+"_eff_plots/";

  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;

  TCanvas *aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());

  std::cout<<"objectPath: "<<objectPath<<" aEffCanvas: "<<aEffCanvas<<std::endl;
  if(!aEffCanvas) return 0;

  TGraphAsymmErrors *hxy_fit_eff = (TGraphAsymmErrors*)aEffCanvas->FindObject(("hxy_"+fitType+"_eff").c_str());

  return hxy_fit_eff;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphAsymmErrors* getRatioGraph(TGraphAsymmErrors *graph1,
				 TGraphAsymmErrors *graph2){

  TGraphAsymmErrors *grRatio = new TGraphAsymmErrors(3);
  Double_t x1, y1, x2, y2; 
  Double_t eX1, eX2, eY1, eY2;
  for(int iPoint=0;iPoint<3;++iPoint){
    graph1->GetPoint(iPoint,x1,y1);
    graph2->GetPoint(iPoint,x2,y2);
    eX1 = graph1->GetErrorX(iPoint);
    eX2 = graph2->GetErrorX(iPoint);
    eY1 = graph1->GetErrorY(iPoint);
    eY2 = graph2->GetErrorY(iPoint);

    float ratio = y1/y2;
    float error = ratio*sqrt(pow(eY1/y1,2) +  pow(eY2/y2,2));
    grRatio->SetPoint(iPoint, x1,ratio);
    grRatio->SetPointError(iPoint,eX1,eX1,error,error);
  }

  grRatio->SetMarkerColor(2);
  grRatio->SetMarkerStyle(21);
  grRatio->SetMarkerSize(1.5);

  return grRatio;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotMistagRateMC(std::string category = "againstMuonLoose3"){

  TGraphAsymmErrors *aGraph  = getResultGraph(category+"_Zll");
  TGraphAsymmErrors *aGraphMCTrue  = getResultGraph(category+"_Zmumu");
  TGraphAsymmErrors *aGraphMCTrueCount  = getResultGraph(category+"_Zmumu", true);
  
  if(!aGraph || !aGraphMCTrue || !aGraphMCTrueCount) return;

  std::cout<<"Zll: "<<std::endl;
  aGraph->Print();
  std::cout<<"Zmumu fit: "<<std::endl;
  aGraphMCTrue->Print();
  std::cout<<"Zmumu count: "<<std::endl;
  aGraphMCTrueCount->Print();
  
  aGraphMCTrue->SetName("aGraphMCTrue");
  aGraphMCTrue->SetLineColor(2);
  aGraphMCTrue->SetMarkerColor(2);

  aGraphMCTrueCount->SetName("aGraphMCTrueCount");
  aGraphMCTrueCount->SetLineColor(4);
  aGraphMCTrueCount->SetMarkerColor(4);

  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(5E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(5E-6);
    hFrame->SetMaximum(1.5E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("|#eta|");
  hFrame->SetYTitle("");

  TLegend *aLegend = new TLegend(0.12,0.1,0.3,0.45,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetFillStyle(4000);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(aGraph,"MC Z#rightarrow ll","lp");
  aLegend->AddEntry(aGraphMCTrue,"#splitline{MC Z#rightarrow #mu#mu with}{tag&probe matched to #mu}","lp");
  aLegend->AddEntry(aGraphMCTrueCount,"#splitline{tag&probe matched to #mu}{event count}","lp");

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
  aGraphMCTrueCount->Draw("p");
  aLegend->Draw();

  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedMC.png").c_str());

  aCanvas->SetLeftMargin(0.1);
  aCanvas->SetWindowSize(700,200);
  TGraphAsymmErrors *grRatio = getRatioGraph(aGraph, aGraphMCTrueCount);
  hFrame->SetYTitle("#frac{DY #rightarrow ll with fit}{DY #rightarrow #mu #mu with count}");
  hFrame->SetMaximum(1.1);
  hFrame->SetMinimum(0.9);
  hFrame->GetYaxis()->SetTitleOffset(0.6);
  hFrame->GetYaxis()->SetLabelSize(0.08);
  hFrame->GetYaxis()->SetTitleSize(0.08);
  hFrame->GetXaxis()->SetLabelSize(0.08);
  hFrame->GetXaxis()->SetTitleSize(0.08);
  hFrame->Draw();  
  grRatio->Draw("p");

  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedMC_Ratio.png").c_str());
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotMistagRateData(std::string category = "againstMuonLoose3"){

  bool useCount = false;
  bool isData = true;
  
  TGraphAsymmErrors *aGraphData  = getResultGraph(category, useCount, isData);

  useCount = true;
  isData = false;
  TGraphAsymmErrors *aGraphMCTrueCount  = getResultGraph(category+"_Zmumu", useCount, isData);

  if(!aGraphData || !aGraphMCTrueCount) return;

  std::cout<<"DATA: "<<std::endl;
  aGraphData->Print("all");
  std::cout<<"MC: "<<std::endl;
  aGraphMCTrueCount->Print("all");
  
  aGraphData->SetName("aGraphData");
  aGraphData->SetLineColor(2);
  aGraphData->SetMarkerColor(2);

  aGraphMCTrueCount->SetName("aGraphMCTrueCount");
  aGraphMCTrueCount->SetLineColor(4);
  aGraphMCTrueCount->SetMarkerColor(4);

  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(4E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(5E-6);
    hFrame->SetMaximum(1.5E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("|#eta|");
  hFrame->SetYTitle("");

  TLegend *aLegend = new TLegend(0.12,0.1,0.3,0.45,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetFillStyle(4000);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(aGraphData,"Data","lp");
  aLegend->AddEntry(aGraphMCTrueCount,"#splitline{tag&probe matched to #mu}{event count}","lp");

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  aCanvas->SetHighLightColor(2);
  aCanvas->Range(0,0,1,1);
  aCanvas->SetFillColor(0);
  aCanvas->SetBorderMode(0);
  aCanvas->SetBorderSize(2);
  aCanvas->SetFrameBorderMode(0);
  
  hFrame->Draw();
  aGraphData->Draw("p");
  aGraphMCTrueCount->Draw("p");
  aLegend->Draw();

  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedData.png").c_str());

  aCanvas->SetLeftMargin(0.1);
  aCanvas->SetWindowSize(700,200);
  TGraphAsymmErrors *grRatio = getRatioGraph(aGraphData, aGraphMCTrueCount);
  hFrame->SetYTitle("#frac{DATA}{DY #rightarrow #mu #mu with count}");
  hFrame->SetMaximum(3.0);
  hFrame->SetMinimum(0.0);
  hFrame->GetYaxis()->SetTitleOffset(0.6);
  hFrame->GetYaxis()->SetLabelSize(0.08);
  hFrame->GetYaxis()->SetTitleSize(0.08);
  hFrame->GetXaxis()->SetLabelSize(0.08);
  hFrame->GetXaxis()->SetTitleSize(0.08);
  hFrame->Draw();  
  grRatio->Draw("p");

  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedData_Ratio.png").c_str());

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotMistagRate(std::string category = "againstMuonLoose3", bool isData=false){

  if(isData) plotMistagRateData(category);
  else plotMistagRateMC(category);

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotMisIdVsMass(){

  /* 
60-70

x[0]=0.577131, y[0]=0.0822561, exl[0]=0.577131, exh[0]=0.622869, eyl[0]=0.000700314, eyh[0]=0.000705739

70-80

x[0]=0.572356, y[0]=0.0224239, exl[0]=0.572356, exh[0]=0.627644, eyl[0]=0.000274601, eyh[0]=0.000277894

80-85

x[0]=0.574981, y[0]=0.00611134, exl[0]=0.574981, exh[0]=0.625019, eyl[0]=0.000123982, eyh[0]=0.000126497

85-95

x[0]=0.572822, y[0]=0.00256906, exl[0]=0.572822, exh[0]=0.627178, eyl[0]=2.86164e-05, eyh[0]=2.8935e-05

95-100

x[0]=0.581184, y[0]=0.00258284, exl[0]=0.581184, exh[0]=0.618816, eyl[0]=7.47759e-05, eyh[0]=7.69588e-05

100-110

x[0]=0.578683, y[0]=0.00298973, exl[0]=0.578683, exh[0]=0.621317, eyl[0]=0.000126367, eyh[0]=0.000131787

110-120

x[0]=0.583116, y[0]=0.00396451, exl[0]=0.583116, exh[0]=0.616884, eyl[0]=0.000261604, eyh[0]=0.000279323
  */
  
  const Int_t n = 7;
  Double_t x[n]   = {65, 75, 82.5, 90, 97.5, 105, 115};
  Double_t y[n]   = {0.0822561, 0.0224239, 0.00611134, 0.00256906, 0.00258284, 0.00298973, 0.00396451}; 
  Double_t exl[n] = {5, 5, 2.5, 5, 2.5, 5, 5};
  Double_t exh[n] = {5, 5, 2.5, 5, 2.5, 5, 5};
  Double_t eyl[n] = {0.000700314, 0.000274601, 0.000123982, 2.86164e-05, 7.47759e-05, 0.000126367, 0.000261604};
  Double_t eyh[n] = {0.000700314, 0.000274601, 0.000123982, 2.86164e-05, 7.47759e-05, 0.000126367, 0.000261604};
 
 TGraphAsymmErrors *aGraph = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
 aGraph->SetMarkerColor(4);
 aGraph->SetMarkerStyle(21);
 aGraph->Print();


 TH1F *hFrame = new TH1F("hFrame","",20,60,120);
 hFrame->SetStats(kFALSE);
 hFrame->SetMinimum(0.001);
 hFrame->SetMaximum(0.09);
 hFrame->GetYaxis()->SetTitleOffset(1.3);
 hFrame->SetXTitle("m_{#mu #mu} [GeV]");
 hFrame->SetYTitle("#mu #rightarrow #mu mis. ID.");

 TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
 aCanvas->SetFillColor(0);
 aCanvas->SetBorderMode(0);
 aCanvas->SetBorderSize(2);
 aCanvas->SetFrameBorderMode(0);
 
 hFrame->Draw();
 aGraph->Draw("P");

 aCanvas->Print("fig_png/MisIdVsMass.png");
 
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void fixParamsForPdf(RooAbsPdf *aPdf){

  RooArgSet *aPdfVariables = aPdf->getVariables();
  aPdfVariables->Print();
  RooFIter aIterator = aPdfVariables->fwdIterator();
  RooRealVar *aVar = (RooRealVar*)aIterator.next();
  while(aVar){
    if(std::string(aVar->GetName()).find("mass")==std::string::npos) aVar->setConstant();
    aVar = (RooRealVar*)aIterator.next();
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void getParamsMC(std::string category = "againstMuonLoose3_Zmumu"){

  TFile file(fNameMC.c_str(),"UPDATE");

  std::string dirName = topDirectory+category;
  std::string objName = "w";
  std::string objectPath = dirName+"/"+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";  
  std::string aPattern = "";
  std::string binnedVars =   "__mcTrue_bin0__";
    
  for(unsigned int iEta=0;iEta<3;++iEta){

    std::string etaBinNumber = std::to_string(iEta);
  
    dirName = topDirectory+category+"/abseta_bin"+etaBinNumber+binnedVars+fitModelName+"_Eta"+etaBinNumber+"/";
    objectPath = dirName+objName;
 
    RooWorkspace *aWorkspace = (RooWorkspace*)file.Get(objectPath.c_str());

    std::cout<<objectPath<<" aWorkspace: "<<aWorkspace<<std::endl;
    if(!aWorkspace) return;
    
    aWorkspace->SetName("workspaceConst");

    std::string pdfNames[5] = {"signalFail", "signalPass", "backgroundFail", "backgroundPass", "backgroundZtautauPass"};
    
    for(auto aName : pdfNames){
      RooAbsPdf *aPdf = (RooAbsPdf*)aWorkspace->obj(aName.c_str());
      if(aPdf) fixParamsForPdf(aPdf);
    }
    file.cd(dirName.c_str());
    aWorkspace->Write("workspaceFixedParams");
  }   
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotAll(){

  bool isData = false;
  
  plotFitCanvas("againstMuonLoose3_Zmumu",isData);
  plotFitCanvas("againstMuonLoose3_Ztautau",isData);
  plotFitCanvas("againstMuonLoose3_Zll",isData);
  plotMistagRate("againstMuonLoose3",isData);

  plotMistagRate("againstMuonTight3",isData);

  isData = true;
  plotFitCanvas("againstMuonLoose3",isData);
  plotMistagRate("againstMuonLoose3",isData);
    
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void fixModelParameters(){

  fNameMC = "TnP_MuonToTau_MisID_MC_Templates.root";

  getParamsMC("againstMuonLoose3_Zmumu");
  getParamsMC("againstMuonLoose3_Ztautau");

  getParamsMC("againstMuonTight3_Zmumu");
  getParamsMC("againstMuonTight3_Ztautau");

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
