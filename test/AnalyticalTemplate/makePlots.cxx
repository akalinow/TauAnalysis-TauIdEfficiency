/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
std::string fNameMC =   "TnP_MuonToTau_MisID_MC.root";
std::string fNameMCTemplate =   "TnP_MuonToTau_MisID_MC_Templates.root";
std::string fNameData = "TnP_MuonToTau_MisID_Data.root";

std::string topDirectory = "tpTree/";
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "tdrstyle.C"
#include "CMS_lumi.C"
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFitCanvas(std::string category = "againstMuonLoose3_Zmumu",
		   bool isData=false, 
		   bool allProbes=false,
		   bool isTemplate=false
		   ){

  std::string fName = fNameMC;
  if(isData) fName = fNameData;
  if(isTemplate) fName = fNameMCTemplate;
  TFile file(fName.c_str());

  std::string dirName = topDirectory+category+"/fit_eff_plots/";
  std::string objName = "fit_canvas";
  std::string objectPath = dirName+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";
  std::string binnedVars =   "__";
  std::string workingPointName = "Loose";
  if(category.find("Tight")!=string::npos) workingPointName = "Tight"; 
  fitModelName = "Zll_Model";
  
  if(category.find("Zmumu")!=string::npos || category.find("Ztautau")!=string::npos){
    binnedVars =   "__mcTrue_bin0__";  
  }
  if(isData){
    fitModelName = "Data_Model";
    binnedVars =   "__";
  }
  if(isTemplate){
    std::string type = category.substr(category.find("_")+1);
    std::cout<<"type: "<<type<<std::endl;
    fitModelName = type+"_Model";
    workingPointName = "";
  }

  for(unsigned int iEta=0;iEta<5;++iEta){

    std::string etaBinNumber = std::to_string(iEta);
    std::string nameSuffix = "Eta"+etaBinNumber;
    if(isData) nameSuffix+="Data";
    dirName = topDirectory+category+"/abseta_bin"+etaBinNumber+binnedVars+fitModelName+"_"+workingPointName+"Eta"+etaBinNumber+"/";
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

    float max = hFrame->GetMaximum();
    hFrame->SetMaximum(1.3*max);
    hFrame->GetXaxis()->SetTitle("Tag-Probe mass [GeV/c^{2}]");
    hFrame->GetYaxis()->SetTitleOffset(1.8);
    hFrame->GetYaxis()->SetTitle("Events");
    
    aBkgModel->SetLineColor(2);
    aSumModel->SetLineColor(4);
    
    TLegend *aLegend;
    //if(iEta<3) aLegend = new TLegend(0.2,0.75,0.45,0.9); 
    aLegend = new TLegend(0.55,0.75,0.9,0.9);
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
      
    int iPeriod = 5;//Luminosity period
    int iPosX = 0;
    lumiTextSize = 1.0;
    cmsTextSize = 0.85;
    if(allProbes){
      cmsText = "        CMS";
      extraText = "           Preliminary";
    }
    if(!isTemplate) CMS_lumi(aPad, iPeriod, iPosX);
    if(isTemplate) aPad->SetLogy();
    aLegend->Draw();

    aPad->SetCanvasSize(1000,1000);
    
    if(allProbes){
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

  std::string objName = "abseta_PLOT";
  std::string objectPath = dirName+objName;

  TCanvas *aEffCanvas = (TCanvas*)file.Get(objectPath.c_str());

  std::cout<<"objectPath: "<<objectPath<<" aEffCanvas: "<<aEffCanvas<<std::endl;
  if(!aEffCanvas) return 0;

  TGraphAsymmErrors *hxy_fit_eff = (TGraphAsymmErrors*)aEffCanvas->FindObject(("hxy_"+fitType+"_eff").c_str());

  ///Fix against weighted eta values
  int nPoints = hxy_fit_eff->GetN();
  Double_t x,y;
  std::cout<<"nPoints: "<<nPoints<<std::endl;
  float edges[6] = {0.0, 0.4, 0.8, 1.2, 1.7, 2.3};

  for(unsigned int iPoint=0;iPoint<nPoints;iPoint++){
    float width = edges[iPoint+1] - edges[iPoint];    
    hxy_fit_eff->SetPointEXlow(iPoint,width/2.0);
    hxy_fit_eff->SetPointEXhigh(iPoint,width/2.0);
    hxy_fit_eff->GetPoint(iPoint,x,y);

    float eyHigh = hxy_fit_eff->GetErrorYhigh(iPoint);
    float eyLow = hxy_fit_eff->GetErrorYlow(iPoint);
    float ey = std::min(eyHigh, eyLow);
    hxy_fit_eff->SetPointEYlow(iPoint,ey);
    hxy_fit_eff->SetPointEYhigh(iPoint,ey);
    
    x = edges[iPoint] + (edges[iPoint+1] - edges[iPoint])/2.0;    
    hxy_fit_eff->SetPoint(iPoint,x,y);
  }
  /////////////////////////////////

  return hxy_fit_eff;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphAsymmErrors* getRatioGraph(TGraphAsymmErrors *graph1,
				 TGraphAsymmErrors *graph2){

  int nPoints = graph1->GetN();
  
  TGraphAsymmErrors *grRatio = new TGraphAsymmErrors(nPoints);
  Double_t x1, y1, x2, y2; 
  Double_t eX1, eX2, eY1, eY2;
  for(int iPoint=0;iPoint<nPoints;++iPoint){
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

  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(7E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(0);
    hFrame->SetMaximum(7E-3);
  }

  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("|#eta|");
  hFrame->SetYTitle("#mu #rightarrow #tau misidentification rate");

  TLegend *aLegend = new TLegend(0.55,0.55,0.9,0.9);
  aLegend->SetTextSize(0.05);
  aLegend->SetBorderSize(0);
  aLegend->SetFillColor(10);
  aLegend->AddEntry(aGraph,"MC Z#rightarrow ll","lp");
  aLegend->AddEntry(aGraphMCTrue,"#splitline{MC Z#rightarrow #mu#mu with}{tag&probe matched to #mu}","lp");
  aLegend->AddEntry(aGraphMCTrueCount,"#splitline{tag&probe matched to #mu}{event count}","lp");

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  aCanvas->Divide(1,2);
  TPad *pad1 = (TPad*)aCanvas->GetPad(1);
  TPad *pad2 = (TPad*)aCanvas->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.38);
  pad1->SetTopMargin(0.07);
  pad1->SetBottomMargin(0.13);
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
  hFrame->SetMaximum(1.15);
  hFrame->SetMinimum(0.85);
  hFrame->GetXaxis()->SetLabelColor(1);
  hFrame->GetYaxis()->SetTitleOffset(0.6);
  hFrame->GetYaxis()->SetLabelSize(0.1);
  hFrame->GetYaxis()->SetTitleSize(0.1);
  hFrame->GetXaxis()->SetLabelSize(0.1);
  hFrame->GetXaxis()->SetTitleSize(0.12);
  hFrame->Draw();  
  grRatio->Draw("p");
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

  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(1E-4);
  hFrame->SetMaximum(7E-3);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMinimum(0);
    hFrame->SetMaximum(7E-3);
  }
  hFrame->SetStats(kFALSE);
  hFrame->GetYaxis()->SetTitleOffset(1.2);
  hFrame->SetXTitle("|#eta|");
  hFrame->SetYTitle("#mu #rightarrow #tau misidentification rate");
  TH1F *hFrameClone = (TH1F*)hFrame->Clone("hFrameClone");
  
  TLegend *aLegend = new TLegend(0.6,0.75,0.85,0.9);
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
  pad1->SetBottomMargin(0.13);
  pad2->SetBottomMargin(0.25);

  aCanvas->cd(1);
  hFrame->GetXaxis()->SetLabelColor(10);
  hFrame->DrawCopy();
  aGraphMCTrueCount->Draw("p");
  aGraphData->Draw("p");
  aLegend->Draw();
  aLatex->DrawLatexNDC(0.7,0.6,label.c_str());
  int iPeriod = 5;//Luminosity period
  int iPosX = 0;
  lumiTextSize = 0.9;
  cmsTextSize = 1.0;
  CMS_lumi(pad1, iPeriod, iPosX);
   
  aCanvas->cd(2);
  TGraphAsymmErrors *grRatio = getRatioGraph(aGraphData, aGraphMCTrueCount);
  hFrame->SetYTitle("#frac{DATA}{Simulation}");
  hFrame->SetMaximum(2.0);
  hFrame->SetMinimum(0.5);
  if(category.find("Tight")!=std::string::npos){
    hFrame->SetMaximum(2.0);
    hFrame->SetMinimum(0.5);
  }
  hFrame->GetXaxis()->SetLabelColor(1);
  hFrame->GetYaxis()->SetTitleOffset(0.4);
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
string getExpoEff(float number, float error=0){

  char text[200];
  int base;
  if(number>1) base = (unsigned int)(fabs(log(number)/log(10.0)) + 0.0);
  else base = (unsigned int)(fabs(log(number)/log(10.0)) + 0.99);
  
  if(base<2){
    if(error<0.001) sprintf(text," ${\\mathrm %4.4f \\pm %4.4f}$ ",number,error);
    else if(error<0.01) sprintf(text," ${\\mathrm %3.3f \\pm %3.3f}$ ",number,error);
    else if(error<0.1) sprintf(text," ${\\mathrm %3.2f \\pm %3.2f}$ ",number,error);
    else if(error<10) sprintf(text," ${\\mathrm %3.1f \\pm %3.1f}$ ",number,error);
    else if(error<100) sprintf(text," ${\\mathrm %3.0f \\pm %3.0f}$ ",number,error);
    else if(error<1000) sprintf(text," ${\\mathrm %2.0f \\pm %2.0f}$ ",number,error);
    else if(error<10000) sprintf(text," ${\\mathrm %2.0f \\pm %2.0f}$ ",number,error);
    else std::cout<<"expoEff: error: "<<error<<std::endl;
    string result;
    result.append(text);
    return result;
  }
  int sgn = (unsigned int)(log(number)/fabs(log(number)));
  float tmp = number*pow(10.0,-sgn*base);
  if(base!=0) sprintf(text,"${\\mathrm (%3.2f \\pm",tmp);
  else sprintf(text,"${\\mathrm  %3.2f $",tmp);
  string result;
  result.append(text);

  sgn = (unsigned int)(log(error)/fabs(log(error)));
  tmp = error*pow(10.0,-sgn*base);
  if(base!=0) sprintf(text," %3.2f) \\cdot 10^{%d}}$",tmp,sgn*base);
  else sprintf(text," %3.2f}$",tmp);
  result.append(text);

  return result;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void makeMistagRateTable(){

  Double_t valueMCX, valueMCY, errorMCY;
  Double_t valueDataX, valueDataY, errorDataY;
  Double_t valueRatioX, valueRatioY, errorRatioY;

  float etaRanges[6] = {0.0, 0.4, 0.8, 1.2, 1.7, 2.3};
  
  ofstream out("table.tex");
  out<<"\\documentclass{article}"<<std::endl;
  out<<std::endl;
  out<<"\\begin{document} "<<std::endl;
  out<<"\\begin{center} "<<std::endl;
  out<<"\\begin{tabular}{|c|c|c|c|} \\hline "<<std::endl;
  out<<"  WP & Simulation & Data & Data/Simulation \\\\  "<<std::endl;

  int nPoints = 5;
  
    for(unsigned int iEta=0;iEta<nPoints;++iEta){

      if(iEta==0) out<<"  \\hline \\multicolumn{4}{|c|}{$|\\eta|<"<<etaRanges[iEta+1]<<"$} \\\\ \\hline "<<std::endl;
      else out<<"  \\hline \\multicolumn{4}{|c|}{$"<<etaRanges[iEta]<<"<|\\eta|<"<<etaRanges[iEta+1]<<"$} \\\\ \\hline "<<std::endl;
      
      for(unsigned int iWP=0;iWP<2;++iWP){

	std::string category = "againstMuonLoose3";
	if(iWP==1) category = "againstMuonTight3";
	
	TGraphAsymmErrors *aGraphZllFit  = getResultGraph(category+"_Zll");
	TGraphAsymmErrors *aGraphData  = getResultGraph(category, false, true);
	TGraphAsymmErrors *aGraphZmumuCount  = getResultGraph(category+"_Zmumu", true, false);
	TGraphAsymmErrors *aGraphRatio = getRatioGraph(aGraphData, aGraphZmumuCount);
	if(!aGraphData || !aGraphZmumuCount || !aGraphZllFit) return;

	if(iEta==0){
	std::cout<<category<<std::endl;
	std::cout<<"DATA: "<<std::endl;
	aGraphData->Print("all");
	std::cout<<"Z->mumu count: "<<std::endl;
	aGraphZmumuCount->Print("all");
	std::cout<<"Z->ll fit: "<<std::endl;
	aGraphZllFit->Print("all");
	std::cout<<"DATA/Z->mumu count: "<<std::endl;
	aGraphRatio->Print("all");
	}
	
      
	aGraphZmumuCount->GetPoint(iEta, valueMCX, valueMCY);
	aGraphData->GetPoint(iEta, valueDataX, valueDataY);
	aGraphRatio->GetPoint(iEta, valueRatioX, valueRatioY);
	
	errorDataY = aGraphData->GetErrorY(iEta);
	errorMCY = aGraphZmumuCount->GetErrorY(iEta);
	errorRatioY = aGraphRatio->GetErrorY(iEta);
	
	out<<category<<" "
	   <<" & "<<getExpoEff(valueMCY,errorMCY)
	   <<" & "<<getExpoEff(valueDataY,errorDataY)
	   <<" & "<<getExpoEff(valueRatioY,errorRatioY)
	   <<"\\\\"<<std::endl;
      }
    }
    out<<"\\hline"<<std::endl;
    out<<"\\end{tabular}"<<std::endl;
    out<<"\\end{center}"<<std::endl;
    out<<"\\end{document}"<<std::endl;
    out.close();
  
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
TGraphErrors * getFittedParamGraph(std::string paramName,
				   bool isData = false,
				   std::string category = "againstMuonLoose3",
				   std::string pdfName = "signalPass"){
  
  std::string dirName = topDirectory+category;
  std::string objName = "w";
  std::string objectPath = dirName+"/"+objName;
  std::string fitModelName = "Zll_Model";
  std::string aPattern = "";
  std::string binnedVars =   "__mcTrue_bin0__";
  std::string workingPointName = category.substr(category.find("Muon")+4,5);
  std::string fName = fNameMC;

  if(isData){
    binnedVars = "__Data_";
    fName = fNameData;
    fitModelName = "Model";
  }
  else{
    category+="_Zmumu";
  }
  
  TFile file(fName.c_str());

  float y[5] = {0,0,0,0,0};
  float yError[5] = {0,0,0,0,0};
  
  for(unsigned int iEta=0;iEta<5;++iEta){

    std::string etaBinNumber = std::to_string(iEta);  
    dirName = topDirectory+category+"/abseta_bin"+etaBinNumber+binnedVars+fitModelName+"_"+workingPointName+"Eta"+etaBinNumber+"/";
    objectPath = dirName+objName;
 
    RooWorkspace *aWorkspace = (RooWorkspace*)file.Get(objectPath.c_str());
    std::cout<<objectPath<<" aWorkspace: "<<aWorkspace<<std::endl;
    if(!aWorkspace) return 0;
    
    RooAbsPdf *aPdf = (RooAbsPdf*)aWorkspace->obj(pdfName.c_str());
    RooArgSet *aPdfVariables = aPdf->getVariables();
    RooRealVar *aVar = (RooRealVar*)aPdfVariables->find(paramName.c_str());
    if(aVar){
      float value = aVar->getValV();
      float error = aVar->getError();
      y[iEta] = value;
      yError[iEta] = error;
    }   
  }

  float x[5] = {0.199823, 0.597364, 0.996381, 1.44261, 1.97812};
  float xError[5] = {0.2, 0.2, 0.2, 0.25, 0.32};

  TGraphErrors *gr = new TGraphErrors(5, x, y, xError, yError);

  return gr;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void plotFittedParams(std::string paramName,
		      std::string category = "againstMuonLoose3_Zmumu",
		      std::string pdfName = "signalPass"){
  
  bool isData = false;
  TGraphErrors *grMC = getFittedParamGraph(paramName, isData, category, pdfName);
  isData = true;
  TGraphErrors *grData = getFittedParamGraph(paramName, isData, category, pdfName);

  grMC->SetMarkerColor(2);
  grMC->SetLineColor(2);

  grData->SetLineColor(1);
  grData->SetLineStyle(2);
  grData->SetMarkerColor(1);
  grData->SetMarkerStyle(21);

  grMC->Draw("AC");
  grData->Draw("AC");
  float maxY = max(grMC->GetYaxis()->GetXmax(), grData->GetYaxis()->GetXmax());
  float minY = min(grMC->GetYaxis()->GetXmin(), grData->GetYaxis()->GetXmin());
  maxY+=(maxY-minY)*0.4;
  
  TH1F *hFrame = new TH1F("hFrame","",3,0,2.3);
  hFrame->SetMinimum(minY);
  hFrame->SetMaximum(maxY);
  hFrame->SetStats(kFALSE);
  hFrame->SetXTitle("|#eta|");
  hFrame->SetYTitle(paramName.c_str());

  TLegend *aLegend = new TLegend(0.2,0.75,0.6,0.9);
  aLegend->SetTextSize(0.07);
  aLegend->SetBorderSize(0);
  aLegend->AddEntry(grData,"Observed","lp");
  aLegend->AddEntry(grMC,"Z #rightarrow #mu #mu simul.","lp");

  TCanvas *aCanvas = new TCanvas("aCanvas", "",4,29,700,500);
  hFrame->Draw();
  grMC->Draw("p");
  grData->Draw("p");
  aLegend->Draw();
  std::string workingPointName = category.substr(category.find("Muon")+4,5);
  std::string figName = "fig_png/"+paramName+"_"+workingPointName+".png";
  aCanvas->Print(figName.c_str());   
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void fixParamsForPdf(RooAbsPdf *aPdf){

  RooArgSet *aPdfVariables = aPdf->getVariables();
  aPdfVariables->Print();
  RooFIter aIterator = aPdfVariables->fwdIterator();
  RooRealVar *aVar = (RooRealVar*)aIterator.next();
  while(aVar){
    if(std::string(aVar->GetName()).find("mass")==std::string::npos &&
       std::string(aVar->GetName()).find("sigmaFail")==std::string::npos &&
       std::string(aVar->GetName()).find("sigmaPass")==std::string::npos &&
       std::string(aVar->GetName()).find("meanFail")==std::string::npos &&
       std::string(aVar->GetName()).find("meanPass")==std::string::npos
       ) aVar->setConstant();
    if(std::string(aVar->GetName()).find("meanFail")!=std::string::npos ||
       std::string(aVar->GetName()).find("meanPass")!=std::string::npos){
      aVar->setRange(88.0,92.0);
    }
    
    aVar = (RooRealVar*)aIterator.next();
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void fixParamsMC(std::string category = "againstMuonLoose3_Zmumu"){

  TFile file(fNameMC.c_str(),"UPDATE");

  std::string dirName = topDirectory+category;
  std::string objName = "w";
  std::string objectPath = dirName+"/"+objName;
  std::string fitModelName = category.substr(category.find("3_")+2,category.size())+"_Model";  
  std::string aPattern = "";
  std::string binnedVars =   "__mcTrue_bin0__";
    
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
  //lumi_13TeV  = "36.8 fb^{-1}";  // default is "20.1 fb^{-1}" Run2016
  lumi_13TeV  = "41.4 fb^{-1}";  // default is "20.1 fb^{-1}" Run2017


  plotFittedParams("sigmaPass", "againstMuonLoose3", "signalPass");
  plotFittedParams("sigmaFail", "againstMuonLoose3", "signalFail");
  plotFittedParams("meanPass", "againstMuonLoose3", "signalPass");
  plotFittedParams("meanFail", "againstMuonLoose3", "signalFail");

  plotFittedParams("sigmaPass", "againstMuonTight3", "signalPass");
  plotFittedParams("sigmaFail", "againstMuonTight3", "signalFail");
  plotFittedParams("meanPass", "againstMuonTight3", "signalPass");
  plotFittedParams("meanFail", "againstMuonTight3", "signalFail");
  
  bool isData = false;
  bool allProbes = false;
  bool isTemplate = true;

  plotFitCanvas("againstMuonLoose3_Zmumu",isData, allProbes, isTemplate);
  plotFitCanvas("againstMuonLoose3_Ztautau",isData, allProbes, isTemplate);

  plotFitCanvas("againstMuonLoose3_Zmumu",isData);
  plotFitCanvas("againstMuonTight3_Zmumu",isData);
  
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

  makeMistagRateTable();

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

  fixParamsMC("againstMuonLoose3_Zmumu");
  fixParamsMC("againstMuonLoose3_Ztautau");

  fixParamsMC("againstMuonTight3_Zmumu");
  fixParamsMC("againstMuonTight3_Ztautau");

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
