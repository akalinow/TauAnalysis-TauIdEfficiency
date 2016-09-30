/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::string fNameMC =   "TnP_MuonToTau_MisID_MC.root";
std::string fNameData = "TnP_MuonToTau_MisID_Data.root";

std::string topDirectory = "tpTree/";
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "tdrstyle.C"
#include "CMS_lumi.C"
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFitCanvas(std::string category = "againstMuonLoose3_Zmumu", bool isData=false, bool allProbes=false){

  std::string fName = fNameMC;
  if(isData) fName = fNameData;
  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "fit_canvas";
  std::string objectPath = dirName+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";
  std::string binnedVars =   "";
  std::string workingPointName = "Loose";
  if(category.find("Tight")!=string::npos) workingPointName = "Tight"; 
  fitModelName = "Zll_Model";
  
  if(category.find("Zmumu")!=string::npos){
    binnedVars =   "mcTrue_bin0__";  
  }
  if(isData){
    fitModelName = "Data_Model";
    binnedVars =   "";
  }

  for(unsigned int iPt=0;iPt<5;++iPt){

    std::string etaBinNumber = std::to_string(iPt);
    std::string nameSuffix = "Pt"+etaBinNumber;
    if(isData) nameSuffix+="Data";
    dirName = topDirectory+category+"/"+binnedVars+"pt_bin"+etaBinNumber+"__"+fitModelName+"_"+workingPointName+"Pt"+etaBinNumber+"/";
    objectPath = dirName+objName;
       
    TCanvas *fit_canvas = (TCanvas*)file.Get(objectPath.c_str());

    std::cout<<objectPath<<" fit_canvas: "<<fit_canvas<<std::endl;
    if(!fit_canvas) return;
    
    fit_canvas->Draw();
    fit_canvas->Print(("./fig_png/"+category+"_fit_canvas_"+nameSuffix+".png").c_str());
    TPad *aPad = 0;    
    if(allProbes) aPad = (TPad*)fit_canvas->FindObject("fit_canvas_3");
    else aPad = (TPad*)fit_canvas->FindObject("fit_canvas_1");
    
    aPad->SetLeftMargin(0.15);
    aPad->SetTopMargin(0.05);
    aPad->SetRightMargin(0.05);
    RooHist *aDataHist = (RooHist*)aPad->FindObject("h_data_binned");
    RooCurve *aSumModel = 0;
    if(allProbes) aSumModel = (RooCurve*)aPad->FindObject("simPdf_Norm[mass]");    
    else aSumModel = (RooCurve*)aPad->FindObject("pdfPass_Norm[mass]");        
    RooCurve *aBkgModel = 0;
    if(allProbes)aBkgModel = (RooCurve*)aPad->FindObject("simPdf_Norm[mass]_Comp[backgroundPass,backgroundFail]");
    else aBkgModel = (RooCurve*)aPad->FindObject("pdfPass_Norm[mass]_Comp[backgroundPass]");    
    TFrame *aFrame = (TFrame*)aPad->FindObject("TFrame");
    TH1D *hFrame;
    TAxis *xAxis = (TAxis*)aDataHist->FindObject("xaxis");
    TAxis *yAxis = (TAxis*)aPad->FindObject("yaxis");
    TList *aList = aPad->GetListOfPrimitives();
    TIter next(aList);
    while (TObject *obj = next()){
      std::string objName(obj->GetName());
      if(objName.find("frame_")!=std::string::npos) hFrame = (TH1D*)obj;
    }

    hFrame->GetXaxis()->SetTitle("Tag-Probe mass [GeV/c^{2}]");
    hFrame->GetYaxis()->SetTitleOffset(1.8);
    hFrame->GetYaxis()->SetTitle("Events");
    
    aBkgModel->SetLineColor(2);
    aSumModel->SetLineColor(4);
    
    TLegend *aLegend;
    if(iPt<3) aLegend = new TLegend(0.2,0.75,0.45,0.9); 
    else aLegend = new TLegend(0.55,0.75,0.9,0.9);
    aLegend->SetTextSize(0.05);
    aLegend->SetBorderSize(0);
    aLegend->AddEntry(aDataHist,"Observed","lp");
    aLegend->AddEntry(aSumModel,"S+B fit","lp");
    aLegend->AddEntry(aBkgModel,"B component","lp");

    std::string label = "muon veto}{passing probes}";
    if(allProbes) label = "muon veto}{all probes}";
    if(category.find("Loose")!=std::string::npos) label = "#splitline{Loose " + label;
    if(category.find("Tight")!=std::string::npos) label = "#splitline{Tight " + label;
    TLatex * aLatex = new TLatex(0,0,"");
      
    int iPeriod = 4;//Luminosity period
    int iPosX = 0;
    lumiTextSize = 1.0;
    cmsTextSize = 0.85;
    if(allProbes){
      cmsText = "        CMS";
      extraText = "           Preliminary";
    }
    CMS_lumi(aPad, iPeriod, iPosX);
    aLegend->Draw();

    aPad->SetCanvasSize(1000,1000);
    
    if(allProbes){
      //aPad->SetLogy();
      aLatex->DrawLatexNDC(0.18,0.2,label.c_str());
      aPad->Print(("./fig_png/"+category+"_fit_all"+nameSuffix+".png").c_str());
    }
    else{
      aLatex->DrawLatexNDC(0.18,0.15,label.c_str());
      aPad->Print(("./fig_png/"+category+"_fit_passing"+nameSuffix+".png").c_str());
    }
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

  std::string objName = "pt_PLOT";
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

  TGraphAsymmErrors *grRatio = new TGraphAsymmErrors(5);
  Double_t x1, y1, x2, y2; 
  Double_t eX1, eX2, eY1, eY2;
  for(int iPoint=0;iPoint<5;++iPoint){
    graph1->GetPoint(iPoint,x1,y1);
    graph2->GetPoint(iPoint,x2,y2);
    eX1 = graph1->GetErrorXlow(iPoint);
    eX2 = graph2->GetErrorXhigh(iPoint);
    eY1 = graph1->GetErrorYlow(iPoint);
    eY2 = graph2->GetErrorYhigh(iPoint);

    float ratio = y1/y2;
    float error = ratio*sqrt(pow(eY1/y1,2) +  pow(eY2/y2,2));
    grRatio->SetPoint(iPoint, x1,ratio);
    grRatio->SetPointError(iPoint,eX1,eX2,error,error);
  }

  grRatio->SetMarkerColor(1);
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

  TH1F *hFrame = new TH1F("hFrame","",3,10,100);
  hFrame->SetMinimum(1E-3);
  hFrame->SetMaximum(1.5E-2);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(5E-6);
    hFrame->SetMaximum(1.5E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("p_{T} [GeV/c]");
  hFrame->SetYTitle("");

  TLegend *aLegend = new TLegend(0.6,0.65,0.9,0.9,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(aGraph,"MC Z#rightarrow ll","lp");
  aLegend->AddEntry(aGraphMCTrue,"MC Z#rightarrow #mu#mu, fit","lp");
  aLegend->AddEntry(aGraphMCTrueCount,"MC Z#rightarrow #mu#mu, count","lp");

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  aCanvas->Divide(1,2);
  TPad *pad1 = (TPad*)aCanvas->GetPad(1);
  TPad *pad2 = (TPad*)aCanvas->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.38);
  pad1->SetTopMargin(0.07);
  pad2->SetBottomMargin(0.25);

  aCanvas->cd(1);
  hFrame->GetXaxis()->SetLabelColor(10);
  hFrame->DrawCopy();
  aGraph->Draw("p");
  aGraphMCTrue->Draw("p");
  aGraphMCTrueCount->Draw("p");
  aLegend->Draw();
  
  aCanvas->cd(2);
  TGraphAsymmErrors *grRatio = getRatioGraph(aGraph, aGraphMCTrueCount);
  hFrame->SetYTitle("#frac{DY #rightarrow ll with fit}{DY #rightarrow #mu #mu with count}");
  hFrame->SetMaximum(1.2);
  hFrame->SetMinimum(0.8);
  hFrame->GetXaxis()->SetLabelColor(1);
  hFrame->GetYaxis()->SetTitleOffset(0.6);
  hFrame->GetYaxis()->SetLabelSize(0.1);
  hFrame->GetYaxis()->SetTitleSize(0.1);
  hFrame->GetXaxis()->SetLabelSize(0.1);
  hFrame->GetXaxis()->SetTitleSize(0.12);
  hFrame->Draw();  
  grRatio->Draw("p");

  std::cout<<"Ratio Zll fit/Zmumu count: "<<std::endl;
  grRatio->Print();

  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedMC.png").c_str());

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
  aGraphData->SetLineColor(1);
  aGraphData->SetLineStyle(2);
  aGraphData->SetMarkerColor(1);
  aGraphData->SetMarkerStyle(21);

  aGraphMCTrueCount->SetName("aGraphMCTrueCount");
  aGraphMCTrueCount->SetLineColor(2);
  aGraphMCTrueCount->SetMarkerColor(2);

  TH1F *hFrame = new TH1F("hFrame","",5,10,100);
  hFrame->SetMinimum(2E-3);
  hFrame->SetMaximum(1.5E-2);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(0);
    hFrame->SetMaximum(4.5E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("p_{T}^{#mu} [GeV/c]");
  hFrame->SetYTitle("#mu #rightarrow #tau misidentification rate");
  TH1F *hFrameClone = (TH1F*)hFrame->Clone("hFrameClone");
  
  TLegend *aLegend = new TLegend(0.65,0.75,0.9,0.9);
  aLegend->SetTextSize(0.07);
  aLegend->SetBorderSize(0);
  aLegend->AddEntry(aGraphData,"Observed","lp");
  aLegend->AddEntry(aGraphMCTrueCount,"Z #rightarrow #mu #mu simul.","lp");

  std::string label = "}{muon veto}";
  if(category.find("Loose")!=std::string::npos) label = "#splitline{Loose " + label;
  if(category.find("Tight")!=std::string::npos) label = "#splitline{Tight " + label;
  TLatex * aLatex = new TLatex(0,0,"");
  aLatex->SetTextSize(0.07);

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  aCanvas->Divide(1,2);
  TPad *pad1 = (TPad*)aCanvas->GetPad(1);
  TPad *pad2 = (TPad*)aCanvas->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.38);
  pad1->SetTopMargin(0.07);
  pad2->SetBottomMargin(0.25);

  aCanvas->cd(1);
  hFrame->GetXaxis()->SetLabelColor(10);
  hFrame->DrawCopy();
  aGraphMCTrueCount->Draw("p");
  aGraphData->Draw("p");
  aLegend->Draw();
  aLatex->DrawLatexNDC(0.7,0.6,label.c_str());
  int iPeriod = 4;//Luminosity period
  int iPosX = 0;
  lumiTextSize = 0.9;
  cmsTextSize = 1.0;
  CMS_lumi(pad1, iPeriod, iPosX);
   
  aCanvas->cd(2);
  TGraphAsymmErrors *grRatio = getRatioGraph(aGraphData, aGraphMCTrueCount);
  hFrame->SetYTitle("#frac{DATA}{Simulation}");
  hFrame->SetMaximum(1.4);
  hFrame->SetMinimum(0.6);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMaximum(3.2);
    hFrame->SetMinimum(0.9);
  }
  hFrame->GetXaxis()->SetLabelColor(1);
  hFrame->GetYaxis()->SetTitleOffset(0.6);
  hFrame->GetYaxis()->SetLabelSize(0.1);
  hFrame->GetYaxis()->SetTitleSize(0.1);
  hFrame->GetXaxis()->SetLabelSize(0.1);
  hFrame->GetXaxis()->SetTitleSize(0.12);
  hFrame->Draw();  
  grRatio->Draw("p");
  grRatio->Print();
  
  aCanvas->Print(("fig_png/"+category+"_misTagRateMCvsMatchedData.png").c_str());
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
  std::string binnedVars = "__mcTrue_bin0__";
    
  for(unsigned int iEta=0;iEta<5;++iEta){

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
void getDataHistPdf(std::string category = "againstMuonLoose3_Zmumu"){

  TFile *file = new TFile(fNameMC.c_str(),"UPDATE");

  std::string dirName = topDirectory+category;
  std::string objName = "w";
  std::string objectPath = dirName+"/"+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";  
  std::string aPattern = "";
  std::string binnedVars =   "mcTrue_bin0__";
    
  for(unsigned int iPt=0;iPt<5;++iPt){

    std::string ptBinNumber = std::to_string(iPt);
  
    dirName = topDirectory+category+"/"+binnedVars+"pt_bin"+ptBinNumber+"__"+fitModelName+"_Pt"+ptBinNumber+"/";
    objectPath = dirName+objName;
 
    RooWorkspace *aWorkspace = (RooWorkspace*)file->Get(objectPath.c_str());

    std::cout<<objectPath<<" aWorkspace: "<<aWorkspace<<std::endl;
    if(!aWorkspace) return;
    
    RooRealVar *mass = (RooRealVar*)aWorkspace->obj("mass");
    RooCategory *aCategory = (RooCategory*)aWorkspace->obj("_efficiencyCategory_");
    RooAbsData* dataAll = (RooAbsData*)aWorkspace->obj("data_binned");
    RooDataHist* dataPass = (RooDataHist*)dataAll->reduce(RooArgSet(*mass),"_efficiencyCategory_==_efficiencyCategory_::Passed"); 
    RooDataHist* dataFail = (RooDataHist*)dataAll->reduce(RooArgSet(*mass),"_efficiencyCategory_==_efficiencyCategory_::Failed");

    std::string namePrefix = "signal";
    if(category.find("Zmumu")!=std::string::npos) namePrefix = "signal";
    else if(category.find("Ztautau")!=std::string::npos) namePrefix = "backgroundZtautau";

    RooHistPdf *aHistPassPdf = new RooHistPdf((namePrefix+"Pass").c_str(),"",RooArgSet(*mass),*dataPass,2);
    RooHistPdf *aHistFailPdf = new RooHistPdf((namePrefix+"Fail").c_str(),"",RooArgSet(*mass),*dataFail,2);
    /*
    TCanvas *canvas = new TCanvas("fit_canvas");
    RooPlot* frame = mass->frame();
    //dataFail->plotOn(frame);
    //aHistFailPdf->plotOn(frame);
    dataPass->plotOn(frame);
    aHistPassPdf->plotOn(frame);
    frame->Draw();
    canvas->Draw();
    */
    file->cd(dirName.c_str());

    RooWorkspace *aNewWorkspace = new RooWorkspace("workspaceHistPdf");
    aNewWorkspace->import(*aHistPassPdf);
    aNewWorkspace->import(*aHistFailPdf);
    if(category.find("Ztautau")!=std::string::npos){
      namePrefix = "background";
      aHistPassPdf->SetName((namePrefix+"Pass").c_str());
      aHistFailPdf->SetName((namePrefix+"Fail").c_str());
      aNewWorkspace->import(*aHistPassPdf);
      aNewWorkspace->import(*aHistFailPdf);
    }    
    aNewWorkspace->Write("workspaceHistPdf");
  }
  delete file;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotDitributions(std::string variable){
  
  TFile *mcFile = new TFile("/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/16_06_2016/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v17_ext4/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_M_50_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v17_ext4/160531_115349/0000/tnpZ_MC.root");
  TFile *dataFile = new TFile("/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP/SingleMuon_Run2016B_PromptReco_v2_v30/SingleMuon/SingleMuon_Run2016B_PromptReco_v2_v30/160707_133020/0000/tnpZ_Data.root");
  
  TTree *dataTree = (TTree*)dataFile->Get("tpTree/fitter_tree");
  TTree *mcTree = (TTree*)mcFile->Get("tpTree/fitter_tree");

  TCut baseCut = "tag_triggerMatch==1 && tag_dB<0.004 && pair_probeMultiplicity==1 && pair_MET<25 &&  pair_MTtag<40 && abs(pair_dz)<0.05 && pair_deltaR>0.5 && decayModeFinding==1 && byLooseCombinedIsolationDeltaBetaCorr3Hits==1 && mass>85 && mass<95";
  TCut mcCut = "mcTrue==1";
  TCut etaCut = "1";

  
  TH1F *hMC = 0;  
  if(variable.find("pt")!=std::string::npos){
    hMC = new TH1F("hMC","",20,0,100);
    hMC->SetXTitle("p_{T}^{probe} [GeV/c]");
  }
  if(variable.find("eta")!=std::string::npos){
    hMC = new TH1F("hMC","",10,-2.4,2.4);
    hMC->SetXTitle("#eta^{probe}");
  }
  if(variable=="mass"){
    hMC = new TH1F("hMC","",10,80,100);
    hMC->SetXTitle("m_{tag-probe} [GeV/c^{2}");
  }
  if(variable.find("MT")!=std::string::npos){
    hMC = new TH1F("hMC","",10,00,100);
    hMC->SetXTitle("m_{T}^{tag-MET} [GeV/c^{2}");
  }

  TH1F *hData = (TH1F*)hMC->Clone("hData");

  mcTree->Draw((variable+">>hMC").c_str(),baseCut && mcCut && etaCut,"goff");
  dataTree->Draw((variable+">>hData").c_str(),baseCut && etaCut,"goff");

  hMC->Scale(1.0/hMC->Integral(0,hMC->GetNbinsX()+1));
  hData->Scale(1.0/hData->Integral(0,hData->GetNbinsX()+1));

  if(hData->GetMaximum()>hMC->GetMaximum()){
    hMC->SetMaximum(1.05*hData->GetMaximum());
  }

  hMC->SetLineWidth(3);
  hData->SetLineWidth(3);

  hMC->SetLineColor(2);
  hData->SetLineColor(1);

  hMC->SetStats(kFALSE);

  TCanvas *aCanvas = new TCanvas("aCanvas","",460,500);
  //aCanvas->SetLogy();
  hMC->Draw();
  hData->Draw("same");

  TLegend *aLegend = new TLegend(0.7,0.7,0.9,0.85,NULL,"brNDC");
  aLegend->SetTextSize(0.05);
  aLegend->SetFillStyle(4000);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(hMC,"DY MC","l");
  aLegend->AddEntry(hData,"DATA","l");
  aLegend->Draw();
  aCanvas->Print(TString::Format("fig_png/%s.png",variable.c_str()).Data());
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotAll(){

  setTDRStyle();
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_13TeV  = "15.9 fb^{-1}";  // default is "20.1 fb^{-1}"

  
  plotFitCanvas("againstMuonLoose3_Zll",false);
  plotMistagRate("againstMuonLoose3",false);

  plotFitCanvas("againstMuonLoose3",true);
  plotMistagRate("againstMuonLoose3",true);

  plotFitCanvas("againstMuonTight3",true);
  plotMistagRate("againstMuonTight3",true);
  return;


  
  bool isData = false;

  plotFitCanvas("againstMuonLoose3_Zmumu",isData);
  plotFitCanvas("againstMuonLoose3_Ztautau",isData);

  plotFitCanvas("againstMuonLoose3_Zll",isData);
  plotFitCanvas("againstMuonTight3_Zll",isData);
  
  plotMistagRate("againstMuonLoose3",isData);
  plotMistagRate("againstMuonTight3",isData);

  isData = true;
  plotFitCanvas("againstMuonLoose3",isData);
  plotFitCanvas("againstMuonLoose3",isData,true);
  plotMistagRate("againstMuonLoose3",isData);

  plotFitCanvas("againstMuonTight3",isData);
  plotFitCanvas("againstMuonTight3",isData,true);
  plotMistagRate("againstMuonTight3",isData);


  return;
  plotDitributions("pt");
  plotDitributions("eta");
  plotDitributions("mass");
  plotDitributions("tag_pt");
  plotDitributions("tag_eta");
  plotDitributions("pair_MTtag");
      
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void fixModelParameters(){

  fNameMC = "TnP_MuonToTau_MisID_MC_Templates.root";

  getDataHistPdf("againstMuonLoose3_Zmumu");
  getDataHistPdf("againstMuonLoose3_Ztautau");

  getDataHistPdf("againstMuonTight3_Zmumu");
  getDataHistPdf("againstMuonTight3_Ztautau");
  return;

  getParamsMC("againstMuonLoose3_Zmumu");
  getParamsMC("againstMuonLoose3_Ztautau");

  getParamsMC("againstMuonTight3_Zmumu");
  getParamsMC("againstMuonTight3_Ztautau");

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
