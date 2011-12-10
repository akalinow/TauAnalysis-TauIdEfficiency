
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <limits>

TH1* getHistogram(TFile* inputFile, const TString& dqmDirectory, const TString& meName)
{  
  TString histogramName = dqmDirectory;
  if ( histogramName.Length() > 0 && !histogramName.EndsWith("/") ) histogramName.Append("/");
  histogramName.Append(meName);

  TH1* histogram = (TH1*)inputFile->Get(histogramName.Data());
  std::cout << "histogramName = " << histogramName.Data() << ": histogram = " << histogram;
  if ( histogram ) std::cout << ", integral = " << histogram->Integral();
  std::cout << std::endl; 

  if ( !histogram->GetSumw2N() ) histogram->Sumw2();

  if ( histogram->GetDimension() == 1 ) histogram->Rebin(5);

  return histogram;
}

TH1* getHistogram(TFile* inputFile, const TString& meName, const TString& process, const TString& region)
{
  TString suffix = "all";
  if      ( region.Contains("p") ) suffix = "passed";
  else if ( region.Contains("f") ) suffix = "failed";
  TString meName_full = Form(meName.Data(), process.Data(), region.Data(), suffix.Data());
  return getHistogram(inputFile, "", meName_full);
}

TH1* getHistogram(TFile* inputFile, const TString& meName, TObjArray& processes, const TString& region)
{
  TH1* histogram_sum = 0;
  
  int numProcesses = processes.GetEntries();
  for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
    TObjString* process = dynamic_cast<TObjString*>(processes.At(iProcess));
    assert(process);
    TH1* histogram = getHistogram(inputFile, meName, process->GetString(), region);
    if ( !histogram_sum ) {
      TString histogramName_sum = histogram->GetName();
      histogramName_sum = histogramName_sum.ReplaceAll(process->GetString(), "sum");
      histogram_sum = (TH1*)histogram->Clone(histogramName_sum.Data());
      if ( !histogram_sum->GetSumw2N() ) histogram_sum->Sumw2();
    } else {
      histogram_sum->Add(histogram);
    }
  }

  return histogram_sum;
}

//-------------------------------------------------------------------------------
//
// Integral of Crystal Ball function for fitting trigger efficiency turn-on curves
// (code from Pascal Paganini)
//

double integralCrystalBall(double m, double m0, double sigma, double alpha, double n, double norm) 
{
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  
  double sig = fabs((double)sigma);
  
  double t = (m - m0)/sig;
  
  if (alpha < 0) t = -t;
  
  double absAlpha = fabs(alpha / sig);
  double a = TMath::Power(n/absAlpha, n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  
  if ( a >= std::numeric_limits<double>::max() ) return -1.;
  
  double approxErf;
  double arg = absAlpha / sqrt2;
  if      ( arg >  5. ) approxErf =  1;
  else if ( arg < -5. ) approxErf = -1;
  else                  approxErf = erf(arg);
  
  double leftArea = (1 + approxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  
  if ( t <= absAlpha ) {
    arg = t / sqrt2;
    if      ( arg >  5.) approxErf =  1;
    else if ( arg < -5.) approxErf = -1;
    else                 approxErf = erf(arg);
    return norm * (1 + approxErf) * sqrtPiOver2 / area;
  } else {
    return norm * (leftArea +  a * (1/TMath::Power(t-b, n - 1) - 1/TMath::Power(absAlpha - b, n - 1)) / (1 - n)) / area;
  }
}
//-------------------------------------------------------------------------------

Double_t integralCrystalBall_f(Double_t* x, Double_t* par) 
{
  return integralCrystalBall(x[0], par[0], par[1], par[2], par[3], par[4]);
}

Double_t integralCrystalBall_f_div_f(Double_t* x, Double_t* par) 
{
  double test = integralCrystalBall(x[0], par[0], par[1], par[2], par[3], par[4]);
  double ref  = integralCrystalBall(x[0], par[5], par[6], par[7], par[8], par[9]);
  return ( ref > 0. ) ? (test/ref - 1.) : 0.;
}

double square(double x)
{
  return x*x;
}

TGraphAsymmErrors* makeGraph_data_div_mc(const TGraph* graph_data, const TGraph* graph_mc)
{
  TGraphAsymmErrors* graph_data_div_mc = new TGraphAsymmErrors(graph_data->GetN());
  
  for ( int iPoint = 0; iPoint < graph_data->GetN(); ++iPoint ) {
    double x_data, y_data;
    graph_data->GetPoint(iPoint, x_data, y_data);
    double yErrUp_data = graph_data->GetErrorYhigh(iPoint);
    double yErrDown_data = graph_data->GetErrorYlow(iPoint);
    
    double x_mc, y_mc;
    graph_mc->GetPoint(iPoint, x_mc, y_mc);
    double yErrUp_mc = graph_mc->GetErrorYhigh(iPoint);
    double yErrDown_mc = graph_mc->GetErrorYlow(iPoint);
    
    assert(x_data == x_mc);
    
    if ( !(y_mc > 0.) ) continue;

    double yDiv = (y_data - y_mc)/y_mc;
    double yDivErrUp = 0.;
    if ( y_data > 0. ) yDivErrUp += square(yErrUp_data/y_data);
    if ( y_mc   > 0. ) yDivErrUp += square(yErrDown_mc/y_mc);
    yDivErrUp *= square(y_data/y_mc);
    yDivErrUp = TMath::Sqrt(yDivErrUp);
    double yDivErrDown = 0.;
    if ( y_data > 0. ) yDivErrDown += square(yErrDown_data/y_data);
    if ( y_mc   > 0. ) yDivErrDown += square(yErrUp_mc/y_mc);
    yDivErrDown *= square(y_data/y_mc);
    yDivErrDown = TMath::Sqrt(yDivErrDown);
    
    //std::cout << "x = " << x_data << ": y = " << yDiv << " + " << yDivErrUp << " - " << yDivErrDown << std::endl;
    
    graph_data_div_mc->SetPoint(iPoint, x_data, yDiv);
    graph_data_div_mc->SetPointError(iPoint, 0., 0., yDivErrDown, yDivErrUp);
  }

  return graph_data_div_mc;
}

void makePlot(const TString& title, 
	      TGraphAsymmErrors* graph_Data_passed, TF1* fit_Data_passed, const TString& legendEntry_Data_passed,
	      TGraphAsymmErrors* graph_mcSum_passed, TF1* fit_mcSum_passed, const TString& legendEntry_mcSum_passed,
	      TGraphAsymmErrors* graph_Data_failed, TF1* fit_Data_failed, const TString& legendEntry_Data_failed,
	      TGraphAsymmErrors* graph_mcSum_failed, TF1* fit_mcSum_failed, const TString& legendEntry_mcSum_failed,
	      const TString& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 900);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", 10, 0., 100.);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  dummyHistogram_top->SetMaximum(1.2);
  dummyHistogram_top->SetMinimum(0.);
  
  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetTitle("calo-E_{T}^{miss} / GeV");
  xAxis_top->SetTitleOffset(1.15);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle("#varepsilon");
  yAxis_top->SetTitleOffset(1.2);

  dummyHistogram_top->Draw("axis");

  graph_Data_passed->SetLineColor(4);
  graph_Data_passed->SetMarkerColor(4);
  graph_Data_passed->SetMarkerStyle(20);
  graph_Data_passed->Draw("p");

  fit_Data_passed->SetLineColor(graph_Data_passed->GetLineColor());
  fit_Data_passed->SetLineWidth(2);
  fit_Data_passed->Draw("same");

  graph_mcSum_passed->SetLineColor(7);
  graph_mcSum_passed->SetMarkerColor(7);
  graph_mcSum_passed->SetMarkerStyle(24);
  graph_mcSum_passed->Draw("p");

  fit_mcSum_passed->SetLineColor(graph_mcSum_passed->GetLineColor());
  fit_mcSum_passed->SetLineWidth(2);
  fit_mcSum_passed->Draw("same");

  graph_Data_failed->SetLineColor(2);
  graph_Data_failed->SetMarkerColor(2);
  graph_Data_failed->SetMarkerStyle(21);
  graph_Data_failed->Draw("p");

  fit_Data_failed->SetLineColor(graph_Data_failed->GetLineColor());
  fit_Data_failed->SetLineWidth(2);
  fit_Data_failed->Draw("same");
 
  graph_mcSum_failed->SetLineColor(6);
  graph_mcSum_failed->SetMarkerColor(6);
  graph_mcSum_failed->SetMarkerStyle(25);
  graph_mcSum_failed->Draw("p");

  fit_mcSum_failed->SetLineColor(graph_mcSum_failed->GetLineColor());
  fit_mcSum_failed->SetLineWidth(2);
  fit_mcSum_failed->Draw("same");
 
  TLegend* legend = new TLegend(0.61, 0.16, 0.89, 0.47, "", "brNDC"); 
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(graph_Data_passed,  legendEntry_Data_passed.Data(),  "p");
  legend->AddEntry(graph_mcSum_passed, legendEntry_mcSum_passed.Data(), "p");
  legend->AddEntry(graph_Data_failed,  legendEntry_Data_failed.Data(),  "p");
  legend->AddEntry(graph_mcSum_failed, legendEntry_mcSum_failed.Data(), "p");
  legend->Draw();

  TPaveText* label = 0;
  if ( title.Length() > 0 ) {
    label = new TPaveText(0.175, 0.89, 0.48, 0.94, "NDC");
    label->AddText(title.Data());
    label->SetTextAlign(13);
    label->SetTextSize(0.045);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->Draw();
  }
  
  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();
  
  TH1* dummyHistogram_bottom = new TH1D("dummyHistogram_bottom", "dummyHistogram_bottom", 10, 0., 100.);
  
  dummyHistogram_bottom->SetMinimum(-1.0);
  dummyHistogram_bottom->SetMaximum(+1.0);

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle("calo-E_{T}^{miss} / GeV");
  xAxis_bottom->SetTitleOffset(1.20);
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);

  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Data-Simulation}{Simulation}");
  yAxis_bottom->SetTitleOffset(0.85);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04);

  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);
  dummyHistogram_bottom->Draw("axis");
 
  TGraphAsymmErrors* graph_Data_div_mc_passed = makeGraph_data_div_mc(graph_Data_passed, graph_mcSum_passed);
  graph_Data_div_mc_passed->SetLineColor(graph_Data_passed->GetLineColor());
  graph_Data_div_mc_passed->SetMarkerColor(graph_Data_passed->GetMarkerColor());
  graph_Data_div_mc_passed->SetMarkerStyle(graph_Data_passed->GetMarkerStyle());
  graph_Data_div_mc_passed->Draw("p");
  
  TF1* fit_Data_div_mc_passed = 
    new TF1("fit_Data_div_mc_passed", &integralCrystalBall_f_div_f, 
	    fit_mcSum_passed->GetMinimumX(), fit_mcSum_passed->GetMaximumX(), 2*fit_mcSum_passed->GetNpar());
  for ( int iPar = 0; iPar < fit_mcSum_passed->GetNpar(); ++iPar ) {
    fit_Data_div_mc_passed->SetParameter(iPar, fit_Data_passed->GetParameter(iPar));
    fit_Data_div_mc_passed->SetParameter(iPar + fit_mcSum_passed->GetNpar(), fit_mcSum_passed->GetParameter(iPar));
  }
  fit_Data_div_mc_passed->SetLineColor(graph_Data_div_mc_passed->GetLineColor());
  fit_Data_div_mc_passed->SetLineWidth(2);
  fit_Data_div_mc_passed->Draw("same");

  TGraphAsymmErrors* graph_Data_div_mc_failed = makeGraph_data_div_mc(graph_Data_failed, graph_mcSum_failed);
  graph_Data_div_mc_failed->SetLineColor(graph_Data_failed->GetLineColor());
  graph_Data_div_mc_failed->SetMarkerColor(graph_Data_failed->GetMarkerColor());
  graph_Data_div_mc_failed->SetMarkerStyle(graph_Data_failed->GetMarkerStyle());
  graph_Data_div_mc_failed->Draw("p");
  
  TF1* fit_Data_div_mc_failed = 
    new TF1("fit_Data_div_mc_failed", &integralCrystalBall_f_div_f, 
	    fit_mcSum_failed->GetMinimumX(), fit_mcSum_failed->GetMaximumX(), 2*fit_mcSum_failed->GetNpar());
  for ( int iPar = 0; iPar < fit_mcSum_failed->GetNpar(); ++iPar ) {
    fit_Data_div_mc_failed->SetParameter(iPar, fit_Data_failed->GetParameter(iPar));
    fit_Data_div_mc_failed->SetParameter(iPar + fit_mcSum_failed->GetNpar(), fit_mcSum_failed->GetParameter(iPar));
  }
  fit_Data_div_mc_failed->SetLineColor(graph_Data_div_mc_failed->GetLineColor());
  fit_Data_div_mc_failed->SetLineWidth(2);
  fit_Data_div_mc_failed->Draw("same");
 
  canvas->Update();
  canvas->Print(outputFileName.Data());

  delete legend;
  delete label;
  delete dummyHistogram_top;
  delete topPad;
  delete dummyHistogram_bottom;
  delete bottomPad;
  delete canvas;
}

void getBinomialBounds(Int_t n, Int_t r, Float_t& rMin, Float_t& rMax)
{
  rMin = 0.;
  rMax = 0.;

  if ( n == 0 ){
    return;
  }
  if ( r < 0 ){
    std::cerr << "Error in <getBinomialBounds>: n = " << n << ", r = " << r << std::endl;
    return;
  }
  
  if ( ((Double_t)r*(n - r)) > (9.*n) ){
    rMin = r - TMath::Sqrt((Double_t)r*(n - r)/((Double_t)n));
    rMax = r + TMath::Sqrt((Double_t)r*(n - r)/((Double_t)n));
    return;
  }

  Double_t binomialCoefficient = 1.;

  Double_t rMinLeft       = 0.;
  Double_t rMinMiddle     = TMath::Max(0.5*r, n - 1.5*r);
  Double_t rMinRight      = n;
  Double_t rMinLeftProb   = 0.;
  Double_t rMinMiddleProb = 0.5;
  Double_t rMinRightProb  = 1.;
  while ( (rMinRight - rMinLeft) > (0.001*n) ){

    rMinMiddleProb = 0;
    for ( Int_t i = r; i <= n; i++ ){
      binomialCoefficient = 1;

      for ( Int_t j = n; j > i; j-- ){
        binomialCoefficient *= j/((Double_t)(j - i));
      }

      rMinMiddleProb += binomialCoefficient*TMath::Power(rMinMiddle/((Double_t)(n)), i)
                       *TMath::Power((n - rMinMiddle)/((Double_t)(n)), n - i);
    }

    if ( rMinMiddleProb > 0.16 ){
      rMinRight     = rMinMiddle;
      rMinRightProb = rMinMiddleProb;
    } else if ( rMinMiddleProb < 0.16 ){
      rMinLeft      = rMinMiddle;
      rMinLeftProb  = rMinMiddleProb;
    } else {
      rMinLeft      = rMinRight     = rMinMiddle;
      rMinLeftProb  = rMinRightProb = rMinMiddleProb;
    }

    rMinMiddle = 0.5*(rMinLeft + rMinRight);

    if ( rMinLeft > r ){
      rMinMiddle = rMinLeft = rMinRight = 0;
    }
  }

  Double_t rMaxLeft       = 0.;
  Double_t rMaxMiddle     = TMath::Min(1.5*r, n - 0.5*r);
  Double_t rMaxRight      = n;
  Double_t rMaxLeftProb   = 1.;
  Double_t rMaxMiddleProb = 0.5;
  Double_t rMaxRightProb  = 0.;
  while ( (rMaxRight - rMaxLeft) > (0.001*n) ){

    rMaxMiddleProb = 0;
    for ( Int_t i = 0; i <= r; i++ ){
      binomialCoefficient = 1;
      
      for ( Int_t j = n; j > (n - i); j-- ){
        binomialCoefficient *= j/((Double_t)(i - (n - j)));
      }

      rMaxMiddleProb += binomialCoefficient*TMath::Power(rMaxMiddle/((Double_t)(n)), i)
                       *TMath::Power((n - rMaxMiddle)/((Double_t)(n)), n - i);
    }

    if ( rMaxMiddleProb > 0.16 ){
      rMaxLeft      = rMaxMiddle;
      rMaxLeftProb  = rMaxMiddleProb;
    } else if ( rMaxMiddleProb < 0.16 ){
      rMaxRight     = rMaxMiddle;
      rMaxRightProb = rMaxMiddleProb;
    } else {
      rMaxLeft      = rMaxRight     = rMaxMiddle;
      rMaxLeftProb  = rMaxRightProb = rMaxMiddleProb;
    }

    rMaxMiddle = 0.5*(rMaxLeft + rMaxRight);

    if ( rMaxRight < r ){
      rMaxMiddle = rMaxLeft = rMaxRight = n;
    }
  }

  rMin = rMinMiddle;
  rMax = rMaxMiddle;
}

TGraphAsymmErrors* getEfficiency(const TH1* histogram_numerator, const TH1* histogram_denominator)
{
  Int_t error = 0;
  if ( !(histogram_numerator->GetNbinsX()           == histogram_denominator->GetNbinsX())           ) error = 1;
  if ( !(histogram_numerator->GetXaxis()->GetXmin() == histogram_denominator->GetXaxis()->GetXmin()) ) error = 1;
  if ( !(histogram_numerator->GetXaxis()->GetXmax() == histogram_denominator->GetXaxis()->GetXmax()) ) error = 1;
  
  if ( error ){
    std::cerr << "Error in <getEfficiency>: Dimensionality of histograms does not match !!" << std::endl;
    return 0;
  }
  
  TAxis* xAxis = histogram_numerator->GetXaxis();

  Int_t nBins = xAxis->GetNbins();
  TArrayF x(nBins);
  TArrayF dxUp(nBins);
  TArrayF dxDown(nBins);
  TArrayF y(nBins);
  TArrayF dyUp(nBins);
  TArrayF dyDown(nBins);

  for ( Int_t ibin = 1; ibin <= nBins; ibin++ ){
    Int_t nObs = TMath::Nint(histogram_denominator->GetBinContent(ibin));
    Int_t rObs = TMath::Nint(histogram_numerator->GetBinContent(ibin));

    Float_t xCenter = histogram_denominator->GetBinCenter(ibin);
    Float_t xWidth  = histogram_denominator->GetBinWidth(ibin);

    x[ibin - 1]      = xCenter;
    dxUp[ibin - 1]   = 0.5*xWidth;
    dxDown[ibin - 1] = 0.5*xWidth;
    
    if ( nObs > 0 ){
      Float_t rMin = 0.;
      Float_t rMax = 0.;
      
      getBinomialBounds(nObs, rObs, rMin, rMax);

      y[ibin - 1]      = rObs/((Float_t)nObs);
      dyUp[ibin - 1]   = (rMax - rObs)/((Float_t)nObs);
      dyDown[ibin - 1] = (rObs - rMin)/((Float_t)nObs);
    } else{
      y[ibin - 1]      = 0.;
      dyUp[ibin - 1]   = 0.;
      dyDown[ibin - 1] = 0.;
    }
  }
  
  TString name  = TString(histogram_numerator->GetName()).Append("Graph");
  TString title = histogram_numerator->GetTitle();

  TGraphAsymmErrors* graph = 
    new TGraphAsymmErrors(nBins, x.GetArray(), y.GetArray(), 
			  dxDown.GetArray(), dxUp.GetArray(), dyDown.GetArray(), dyUp.GetArray());

  graph->SetName(name);
  graph->SetTitle(title);

  return graph;
}

void initializeCrystalBall(TF1* fit)
{
  fit->SetParameter(0, 3.e+1);
  fit->SetParameter(1, 2.e+1);
  fit->SetParameter(2, 1.e+1);
  fit->SetParameter(3, 1.e+2);
  fit->SetParameter(4, 1.e+0);
}

void makePlot_wrapper(const TString& title, 
		      const TH1* histogram_numerator_Data_passed, const TH1* histogram_denominator_Data_passed, 
		      const TH1* histogram_numerator_mcSum_passed, const TH1* histogram_denominator_mcSum_passed, 
		      const TH1* histogram_numerator_Data_failed, const TH1* histogram_denominator_Data_failed, 
		      const TH1* histogram_numerator_mcSum_failed, const TH1* histogram_denominator_mcSum_failed, 
		      const TString& legendEntry_Data, const TString& legendEntry_mcSum, 
		      const TString& region_passed, const TString& region_failed,
		      const TString& outputFileName)
{
  TAxis* xAxis = histogram_numerator_Data_passed->GetXaxis();
  Double_t xMin = xAxis->GetXmin();
  Double_t xMax = xAxis->GetXmax();

  std::cout << "processing Data, region " << region_passed.Data() << "..." << std::endl;
  TGraphAsymmErrors* graph_efficiency_Data_passed = 
    getEfficiency(histogram_numerator_Data_passed, histogram_denominator_Data_passed);   
  TF1* fit_efficiency_Data_passed = new TF1("fit_efficiency_Data_passed", &integralCrystalBall_f, xMin, xMax, 5);
  initializeCrystalBall(fit_efficiency_Data_passed);
  graph_efficiency_Data_passed->Fit(fit_efficiency_Data_passed, "0");
  TString legendEntry_Data_passed = Form("%s, %s", legendEntry_Data.Data(), region_passed.Data());

  std::cout << "processing mcSum, region " << region_passed.Data() << "..." << std::endl;
  TGraphAsymmErrors* graph_efficiency_mcSum_passed = 
    getEfficiency(histogram_numerator_mcSum_passed, histogram_denominator_mcSum_passed);   
  TF1* fit_efficiency_mcSum_passed = new TF1("fit_efficiency_mcSum_passed", &integralCrystalBall_f, xMin, xMax, 5);
  initializeCrystalBall(fit_efficiency_mcSum_passed);
  graph_efficiency_mcSum_passed->Fit(fit_efficiency_mcSum_passed, "0");
  TString legendEntry_mcSum_passed = Form("%s, %s", legendEntry_mcSum.Data(), region_passed.Data());

  std::cout << "processing Data, region " << region_failed.Data() << "..." << std::endl;
  TGraphAsymmErrors* graph_efficiency_Data_failed = 
    getEfficiency(histogram_numerator_Data_failed, histogram_denominator_Data_failed);   
  TF1* fit_efficiency_Data_failed = new TF1("fit_efficiency_Data_failed", &integralCrystalBall_f, xMin, xMax, 5);
  initializeCrystalBall(fit_efficiency_Data_failed);
  graph_efficiency_Data_failed->Fit(fit_efficiency_Data_failed, "0");
  TString legendEntry_Data_failed = Form("%s, %s", legendEntry_Data.Data(), region_failed.Data());

  std::cout << "processing mcSum, region " << region_failed.Data() << "..." << std::endl;
  TGraphAsymmErrors* graph_efficiency_mcSum_failed = 
    getEfficiency(histogram_numerator_mcSum_failed, histogram_denominator_mcSum_failed);   
  TF1* fit_efficiency_mcSum_failed = new TF1("fit_efficiency_mcSum_failed", &integralCrystalBall_f, xMin, xMax, 5);
  initializeCrystalBall(fit_efficiency_mcSum_failed);
  graph_efficiency_mcSum_failed->Fit(fit_efficiency_mcSum_failed, "0");
  TString legendEntry_mcSum_failed = Form("%s, %s", legendEntry_mcSum.Data(), region_failed.Data());

  makePlot(title,
	   graph_efficiency_Data_passed,  fit_efficiency_Data_passed,  legendEntry_Data_passed,
	   graph_efficiency_mcSum_passed, fit_efficiency_mcSum_passed, legendEntry_mcSum_passed,
	   graph_efficiency_Data_failed,  fit_efficiency_Data_failed,  legendEntry_Data_failed,
	   graph_efficiency_mcSum_failed, fit_efficiency_mcSum_failed, legendEntry_mcSum_failed,
	   outputFileName);

  delete graph_efficiency_Data_passed;
  delete fit_efficiency_Data_passed;
  delete graph_efficiency_mcSum_passed;
  delete fit_efficiency_mcSum_passed;
  delete graph_efficiency_Data_failed;
  delete fit_efficiency_Data_failed;
  delete graph_efficiency_mcSum_failed;
  delete fit_efficiency_mcSum_failed;
}

void makeCaloMEtTriggerEffPlots_Ztautau()
{
//--- stop ROOT from keeping references to all histograms
  TH1::AddDirectory(false);

//--- suppress the output canvas 
  gROOT->SetBatch(true);

  TString inputFilePath = "/data1/veelken/tmp/muonPtGt17/V10_5tauEnRecovery_L1ETM20Eff/2011RunB/tauIdEfficiency/";
  
  TString inputFileName = "analyzeTauIdEffHistograms_all_2011Oct30V10_5tauEnRecovery.root";

  TString dqmDirectory = "";

  TString meName_numerator_data = "%s_%s_numCaloMEt_HLT_IsoMu15_L1ETM20_tauDiscrHPScombLooseDBcorr_%s";
  TString meName_numerator_mc   = "%s_%s_numCaloMEt_L1_ETM20_tauDiscrHPScombLooseDBcorr_%s"; 
  TString meName_denominator    = "%s_%s_denomCaloMEt_tauDiscrHPScombLooseDBcorr_%s"; 

  TObjArray processes;
  processes.Add(new TObjString("Ztautau"));
  processes.Add(new TObjString("Zmumu"));
  processes.Add(new TObjString("WplusJets"));
  processes.Add(new TObjString("QCD"));
  processes.Add(new TObjString("TTplusJets"));

  TString region_passed = "C1p";
  TString region_failed = "C1f";

  TString inputFileName_full = inputFilePath;
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName);
  std::cout << "opening inputFile = " << inputFileName_full.Data() << std::endl;
  TFile* inputFile = new TFile(inputFileName_full.Data());

  TH1* histogram_numerator_Data_passed    = getHistogram(inputFile, meName_numerator_data, "Data",    region_passed);
  TH1* histogram_denominator_Data_passed  = getHistogram(inputFile, meName_denominator,    "Data",    region_passed);
  TH1* histogram_numerator_mcSum_passed   = getHistogram(inputFile, meName_numerator_mc,   processes, region_passed);
  TH1* histogram_denominator_mcSum_passed = getHistogram(inputFile, meName_denominator,    processes, region_passed);

  TH1* histogram_numerator_Data_failed    = getHistogram(inputFile, meName_numerator_data, "Data",    region_failed);
  TH1* histogram_denominator_Data_failed  = getHistogram(inputFile, meName_denominator,    "Data",    region_failed);
  TH1* histogram_numerator_mcSum_failed   = getHistogram(inputFile, meName_numerator_mc,   processes, region_failed);
  TH1* histogram_denominator_mcSum_failed = getHistogram(inputFile, meName_denominator,    processes, region_failed);

  TString outputFileName = Form("makeCaloMEtTriggerEffPlots_Ztautau.png");
  makePlot_wrapper("",
		   histogram_numerator_Data_passed, histogram_denominator_Data_passed, 
		   histogram_numerator_mcSum_passed, histogram_denominator_mcSum_passed,
		   histogram_numerator_Data_failed, histogram_denominator_Data_failed,  
		   histogram_numerator_mcSum_failed, histogram_denominator_mcSum_failed,
		   "Data", "Simulation", region_passed, region_failed,
		   outputFileName);

  delete inputFile;
}
