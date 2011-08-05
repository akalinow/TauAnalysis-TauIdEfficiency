#ifndef TauAnalysis_GenSimTools_tauIdEffAuxFunctions_h
#define TauAnalysis_GenSimTools_tauIdEffAuxFunctions_h

#include "TauAnalysis/DQMTools/interface/histogramAuxFunctions.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooIntegralMorph.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TPolyMarker3D.h>
#include <TPaveText.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TMatrixD.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <time.h>

bool isSystematicShift(const std::string& key)
{
  if ( TString(key.data()).CountChar('_') >= 3 ) return true;
  else return false;
}

std::string getKey(const std::string& observable, const std::string& tauId, const std::string tauIdValue = "all",
		   const std::string sysShift = "CENTRAL_VALUE")
{
  //std::cout << "<getKey>:" << std::endl;
  
  std::string key = std::string(observable).append("_").append(tauId).append("_").append(tauIdValue);
  if ( sysShift != "CENTRAL_VALUE" ) key.append("_").append(sysShift);
  return key;
}

std::string getObservable(const std::string& key)
{
  //std::cout << "<getObservable>:" << std::endl;
  //std::cout << " key = " << key << std::endl;

  std::string observable;

  size_t idx = key.find("_");
  if ( idx != std::string::npos ) observable = std::string(key, 0, idx);
  else                            observable = key; 

  //std::cout << "--> observable = " << observable << std::endl;

  return observable;
}

std::vector<std::string> getTauIdValues(const std::string& region)
{
  std::vector<std::string> retVal;

  if      ( region.find("p") != std::string::npos ) retVal.push_back(std::string("passed"));
  else if ( region.find("f") != std::string::npos ) retVal.push_back(std::string("failed"));
  else                                              retVal.push_back(std::string("all"));

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

double getIntegral(const TH1* histogram, bool inclUnderflowBin, bool inclOverflowBin)
{
//--------------------------------------------------------------------------------
// Compute integral of histogram including/excluding underflow and overflow bins
//--------------------------------------------------------------------------------

  //std::cout << "<getIntegral>:" << std::endl;
  //std::cout << " histogram        = " << histogram << std::endl;
  //std::cout << " name             = " << histogram->GetName() << std::endl;
  //std::cout << " inclUnderflowBin = " << inclUnderflowBin << std::endl;
  //std::cout << " inclOverflowBin  = " << inclOverflowBin << std::endl;

  int firstBin = ( inclUnderflowBin ) ?                            0 :                      1;
  int lastBin  = ( inclOverflowBin  ) ? (histogram->GetNbinsX() + 1) : histogram->GetNbinsX();
  
  //std::cout << " firstBin = " << firstBin << ", lastBin = " << lastBin << std::endl;

  double integral = 0;
  for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    //std::cout << " binContent(" << iBin << ") = " << binContent << std::endl;
    integral += binContent;
  }

  //std::cout << "--> integral = " << integral << std::endl;

  return integral;
}

TH1* normalize(const TH1* histogram, double norm = 1.)
{
//--------------------------------------------------------------------------------
// Normalize histogram passed as function argument to unit area
// (for use as shape template)
//--------------------------------------------------------------------------------

  //std::cout << "<normalize>:" << std::endl;

  TH1* retVal = (TH1*)histogram->Clone();

  if ( !retVal->GetSumw2N() ) retVal->Sumw2();

  double integral = getIntegral(retVal, true, true);
  if ( integral != 0. ) retVal->Scale(norm/integral);
    
  return retVal;
}

//
//-------------------------------------------------------------------------------
//

double getTemplateNorm_fitted(
         const std::string& process, const std::string& region, const std::string& key, 
	 const std::string& tauId, const std::string& fitVariable,
	 std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > >& normFactorsAll_fitted,
	 std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions)
{
  //std::cout << "<getTemplateNorm_fitted>:" << std::endl;
  //std::cout << " process = " << process << std::endl;
  //std::cout << " region = " << region << std::endl;
  //std::cout << " observable = " << observable << std::endl;
  //std::cout << " tauId = " << tauId << std::endl;
  //std::cout << " fitVariable = " << fitVariable << std::endl;

  //std::cout << " normFactorsAll_fitted = " << normFactorsAll_fitted[process][region][tauId][fitVariable] << std::endl;
  //std::cout << " fittedFractions[process] = " << fittedFractions[process][region][key] << std::endl;

  double retVal = normFactorsAll_fitted[process][region][tauId][fitVariable]*fittedFractions[process][region][key];
  //std::cout << "--> retVal = " << retVal << std::endl; 

  return retVal;
}

TH1* compFittedTemplateShape(const TH1* histogram, const RooAbsPdf* pdf)
{
  const RooIntegralMorph* pdf_morph = dynamic_cast<const RooIntegralMorph*>(pdf);
  if ( pdf_morph ) {
    std::string histogramFittedShapeName = std::string(histogram->GetName()).append("_fittedShape");
    TH1* histogramFittedShape = (TH1*)histogram->Clone(histogramFittedShapeName.data());
    TAxis* xAxis = histogram->GetXaxis();
    RooRealVar x("x", "x", xAxis->GetXmin(), xAxis->GetXmax());
    for ( int iBin = 1; iBin <= (histogramFittedShape->GetNbinsX() + 1); ++iBin ) {
      double xMin = xAxis->GetBinLowEdge(iBin);
      double xMax = xAxis->GetBinUpEdge(iBin);
      TString binLabel = Form("bin%i", iBin);
      x.setRange(binLabel.Data(), xMin, xMax);
      RooAbsReal* pdfIntegral_bin = pdf->createIntegral(x, RooFit::NormSet(x), RooFit::Range(binLabel.Data()));
      histogramFittedShape->SetBinContent(iBin, pdfIntegral_bin->getVal());
      delete pdfIntegral_bin;
    }
    return histogramFittedShape;
  } else throw cms::Exception("getNumber")  
      << "PDF object passed as function argument is not of type RooIntegralMorph !!\n";
}

//
//-------------------------------------------------------------------------------
//

template <typename T>
void applyStyleOption(T* histogram, const std::string& histogramTitle,
		      const std::string& xAxisTitle, const std::string& yAxisTitle = "Events")
{
  //std::cout << "<applyStyleOption>:" << std::endl;
  //std::cout << " histogramName = " << histogram->GetName() << std::endl;

  histogram->SetTitle(histogramTitle.data());

  if ( histogram->GetXaxis() ) {
    histogram->GetXaxis()->SetTitle(xAxisTitle.data());
    histogram->GetXaxis()->SetTitleOffset(1.15);
    //histogram->GetXaxis()->SetTitleSize(0.05); 
    //histogram->GetXaxis()->SetLabelSize(0.05);
  } //else std::cerr << "Histogram = " << histogram->GetName() << " has no valid x-Axis !!" << std::endl;

  if ( histogram->GetYaxis() ) {
    histogram->GetYaxis()->SetTitle(yAxisTitle.data());
    histogram->GetYaxis()->SetTitleOffset(1.65);
    //histogram->GetYaxis()->SetTitleSize(0.05); 
    //histogram->GetYaxis()->SetLabelSize(0.05);
  } //else std::cerr << "Histogram = " << histogram->GetName() << " has no valid y-Axis !!" << std::endl;
}

void drawCMSprelimaryLabels(double xOffset = 0.150, double yOffset = 0.8075)
{
  static TPaveText* cmsPreliminaryLabel = 0;
  if ( !cmsPreliminaryLabel ) {
    cmsPreliminaryLabel = new TPaveText(xOffset, yOffset + 0.0525, xOffset + 0.32, yOffset + 0.0925, "NDC");
    cmsPreliminaryLabel->AddText("CMS Preliminary 2011");
    cmsPreliminaryLabel->SetTextAlign(13);
    cmsPreliminaryLabel->SetTextSize(0.045);
    cmsPreliminaryLabel->SetFillStyle(0);
    cmsPreliminaryLabel->SetBorderSize(0);
  }
  cmsPreliminaryLabel->Draw();

  static TPaveText* cmsLuminosityLabel = 0;
  if ( !cmsLuminosityLabel ) {
    cmsLuminosityLabel = new TPaveText(xOffset + 0.005, yOffset, xOffset + 0.32, yOffset + 0.0400, "NDC");
    cmsLuminosityLabel->AddText("#sqrt{s} = 7 TeV, L = 879.6 pb^{-1}");
    cmsLuminosityLabel->SetTextAlign(13);
    cmsLuminosityLabel->SetTextSize(0.045);
    cmsLuminosityLabel->SetFillStyle(0);
    cmsLuminosityLabel->SetBorderSize(0);
  }
  cmsLuminosityLabel->Draw();
}

void drawHistograms(TH1* histogramZtautau, double normZtautau,
		    TH1* histogramZmumu, double normZmumu,
		    TH1* histogramQCD, double normQCD,
		    TH1* histogramWplusJets, double normWplusJets,
		    TH1* histogramTTplusJets, double normTTplusJets,
		    TH1* histogramData,
		    const std::string& histogramTitle, const std::string& xAxisTitle, 
		    const std::string& outputFileName, 
		    bool runStatTest = false)
{
//--------------------------------------------------------------------------------
// Make control plots of sum(MC) versus Data.
// If normalization factors are passed as function argument,
// normalize all MC distributions accordingly;
// else assume MC distributions passed as function arguments
// are already properly normalized (by cross-section)
//--------------------------------------------------------------------------------

  std::cout << "<drawHistograms>:" << std::endl;
  std::cout << " Ztautau:       histogram = " << histogramZtautau       << ", norm = " << normZtautau       << std::endl;
  std::cout << " Zmumu:         histogram = " << histogramZmumu         << ", norm = " << normZmumu         << std::endl;
  std::cout << " QCD:           histogram = " << histogramQCD           << ", norm = " << normQCD           << std::endl;
  std::cout << " WplusJets:     histogram = " << histogramWplusJets     << ", norm = " << normWplusJets     << std::endl;
  std::cout << " TTplusJets:    histogram = " << histogramTTplusJets    << ", norm = " << normTTplusJets    << std::endl;
  std::cout << " Data:          histogram = " << histogramData          << std::endl;
  std::cout << " Options: " << std::endl;
  std::cout << "  histogramTitle = " << histogramTitle << std::endl;
  std::cout << "  xAxisTitle     = " << xAxisTitle << std::endl;
  std::cout << "  outputFileName = " << outputFileName << std::endl;

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);

//--- scale template histograms to given normalization factors
  TH1* templateZtautau    = ( normZtautau    > 0. ) ? normalize(histogramZtautau,    normZtautau)    : histogramZtautau;
  applyStyleOption(templateZtautau, histogramTitle, xAxisTitle);
  templateZtautau->SetFillStyle(1001);
  templateZtautau->SetFillColor(628);
  
//--- sum Z/gamma* --> mu+ mu- contributions
  TH1* templateZmumu      = ( normZmumu      > 0. ) ? normalize(histogramZmumu,      normZmumu)      : histogramZmumu;
  applyStyleOption(templateZmumu, histogramTitle, xAxisTitle);
  templateZmumu->SetFillStyle(1001);
  templateZmumu->SetFillColor(596);

  TH1* templateQCD        = ( normQCD        > 0. ) ? normalize(histogramQCD,        normQCD)        : histogramQCD;
  applyStyleOption(templateQCD, histogramTitle, xAxisTitle);
  templateQCD->SetFillStyle(1001);
  templateQCD->SetFillColor(797);

  TH1* templateWplusJets  = ( normWplusJets  > 0. ) ? normalize(histogramWplusJets,  normWplusJets)  : histogramWplusJets;
  applyStyleOption(templateWplusJets, histogramTitle, xAxisTitle);
  templateWplusJets->SetFillStyle(1001);
  templateWplusJets->SetFillColor(856);

  TH1* templateTTplusJets = ( normTTplusJets > 0. ) ? normalize(histogramTTplusJets, normTTplusJets) : histogramTTplusJets;
  applyStyleOption(templateTTplusJets, histogramTitle, xAxisTitle);
  templateTTplusJets->SetFillStyle(1001);
  templateTTplusJets->SetFillColor(618);
  
  applyStyleOption(histogramData, histogramTitle, xAxisTitle);
  histogramData->SetLineColor(1);
  histogramData->SetMarkerColor(1);
  histogramData->SetMarkerStyle(20);

//--- draw histograms for individual processes
  THStack smSum("smSum", "smSum");
  smSum.Add(templateTTplusJets);
  smSum.Add(templateZmumu);  
  smSum.Add(templateWplusJets);
  smSum.Add(templateQCD);
  smSum.Add(templateZtautau);

  smSum.SetMaximum(1.5*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));
  smSum.Draw("hist");
  applyStyleOption(&smSum, histogramTitle, xAxisTitle);
  
  histogramData->SetStats(false);
  histogramData->Draw("ep1same");

  TLegend legend(0.64, 0.59, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  legend.AddEntry(histogramData,      "Data",                                       "p");
  legend.AddEntry(templateZtautau,    "Z/#gamma^{*} #rightarrow #tau^{+} #tau^{-}", "f");
  legend.AddEntry(templateQCD,        "QCD",                                        "f");
  legend.AddEntry(templateWplusJets,  "W + jets",                                   "f");
  legend.AddEntry(templateZmumu,      "Z/#gamma^{*} #rightarrow #mu^{+} #mu^{-}",   "f");
  legend.AddEntry(templateTTplusJets, "t#bar{t} + jets",                            "f");
  legend.Draw();

  drawCMSprelimaryLabels();

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  gSystem->mkdir(outputFilePath.data(), true);
  canvas->Print(outputFilePath.append(outputFileName).data());

//--- draw histograms for individual processes except Zmumu and W + jets:
//    the visible mass shapes of Zmumu and W + jets are so similar 
//    that the fit mostly constrains the sum of Zmumu and W + jets backgrounds, 
//    not their individual contributions
  TH1* templateEWKbgSum = (TH1*)templateZmumu->Clone();
  templateEWKbgSum->Add(templateWplusJets);
  applyStyleOption(templateEWKbgSum, histogramTitle, xAxisTitle);
  templateEWKbgSum->SetFillStyle(1001);
  templateEWKbgSum->SetFillColor(856);

  THStack smSum2("smSum2", "smSum2");
  smSum2.Add(templateTTplusJets);
  smSum2.Add(templateEWKbgSum);  
  smSum2.Add(templateQCD);
  smSum2.Add(templateZtautau);

  smSum2.SetMaximum(1.5*TMath::Max(smSum2.GetMaximum(), histogramData->GetMaximum()));
  smSum2.Draw("hist");
  applyStyleOption(&smSum2, histogramTitle, xAxisTitle);
  
  histogramData->Draw("ep1same");

  TLegend legend2(0.64, 0.64, 0.89, 0.89, "", "brNDC"); 
  legend2.SetBorderSize(0);
  legend2.SetFillColor(0);
  
  legend2.AddEntry(histogramData,      "Data",                                       "p");
  legend2.AddEntry(templateZtautau,    "Z/#gamma^{*} #rightarrow #tau^{+} #tau^{-}", "f");
  legend2.AddEntry(templateQCD,        "QCD",                                        "f");
  legend2.AddEntry(templateEWKbgSum,   "EWK Bgr.",                                   "f");
  legend2.AddEntry(templateTTplusJets, "t#bar{t} + jets",                            "f");
  legend2.Draw();

  drawCMSprelimaryLabels();

  canvas->Update();
  std::string outputFilePath2 = std::string("./plots/");
  gSystem->mkdir(outputFilePath2.data(), true);
  TString outputFileName2 = outputFileName;
  outputFileName2.ReplaceAll(".", "_ewkBgSum.");
  canvas->Print(outputFilePath2.append(outputFileName2).data());

//--- draw histograms for all background processes summed plus Ztautau signal 
  TH1* templateSMbgSum = (TH1*)templateZmumu->Clone();
  templateSMbgSum->Add(templateQCD);
  templateSMbgSum->Add(templateWplusJets);
  templateSMbgSum->Add(templateTTplusJets);
  applyStyleOption(templateSMbgSum, histogramTitle, xAxisTitle);
  templateSMbgSum->SetFillStyle(1001);
  templateSMbgSum->SetFillColor(42);

  templateZtautau->SetFillStyle(1001);
  templateZtautau->SetFillColor(46);

  THStack smSum3("smSum3", "smSum3");
  smSum3.Add(templateSMbgSum);
  smSum3.Add(templateZtautau);

  smSum3.SetMaximum(1.5*TMath::Max(smSum3.GetMaximum(), histogramData->GetMaximum()));	
  smSum3.Draw("hist");
  applyStyleOption(&smSum3, histogramTitle, xAxisTitle);

  histogramData->Draw("ep1same");

  TLegend legend3(0.64, 0.71, 0.89, 0.89, "", "brNDC"); 
  legend3.SetBorderSize(0);
  legend3.SetFillColor(0);
  
  legend3.AddEntry(histogramData,   "Data",                            "p");
  legend3.AddEntry(templateZtautau, "Z #rightarrow #tau^{+} #tau^{-}", "f");
  legend3.AddEntry(templateSMbgSum, "#Sigma Backgrounds",              "f");
  legend3.Draw();

  drawCMSprelimaryLabels();

  canvas->Update();
  std::string outputFilePath3 = std::string("./plots/");
  gSystem->mkdir(outputFilePath3.data(), true);
  TString outputFileName3 = outputFileName;
  outputFileName3.ReplaceAll(".", "_smBgSum.");
  canvas->Print(outputFilePath3.append(outputFileName3).data());

//--- draw histograms for all signal and background processes summed
  TH1* templateSMsum = (TH1*)templateSMbgSum->Clone();
  templateSMsum->Add(templateZtautau);  
  templateSMsum->SetFillStyle(1001);
  templateSMsum->SetFillColor(10);
  templateSMsum->SetLineColor(1);
  templateSMsum->SetLineWidth(2);

  templateSMsum->SetMaximum(1.5*TMath::Max(templateSMsum->GetMaximum(), histogramData->GetMaximum()));
  applyStyleOption(templateSMsum, histogramTitle, xAxisTitle);
  templateSMsum->Draw("hist");
  templateSMsum->SetStats(false);

  histogramData->Draw("ep1same");

  TLegend legend4(0.74, 0.76, 0.89, 0.89, "", "brNDC"); 
  legend4.SetBorderSize(0);
  legend4.SetFillColor(0);
  
  legend4.AddEntry(histogramData, "Data", "p");
  legend4.AddEntry(templateSMsum, "MC",   "l");
  legend4.Draw();
  
  drawCMSprelimaryLabels();

  canvas->Update();
  std::string outputFilePath4 = std::string("./plots/");
  gSystem->mkdir(outputFilePath4.data(), true);
  TString outputFileName4 = outputFileName;
  outputFileName4.ReplaceAll(".", "_smSum.");
  canvas->Print(outputFilePath4.append(outputFileName4).data());  

  if ( runStatTest ) {
    std::cout << "<runStatTest>:" << std::endl;
    std::cout << " histogram = " << histogramData->GetName() << std::endl;
    std::cout << " Data: entries = " << getIntegral(histogramData, false, false)
	      << " (" << getIntegral(histogramData, true, true) << " incl. underflow/overflow bins)" << std::endl;
    std::cout << " MC: entries  = " << getIntegral(templateSMsum, false, false)
	      << " (" << getIntegral(templateSMsum, true, true) << " incl. underflow/overflow bins)" << std::endl;
    std::cout << " p(chi^2) = " << histogramData->Chi2Test(templateSMsum, "UW") << std::endl;
    std::cout << " p(KS) = " << histogramData->KolmogorovTest(templateSMsum, "") << std::endl;
  }

  if ( templateZtautau       != histogramZtautau       ) delete templateZtautau;
  if ( templateZmumu         != histogramZmumu         ) delete templateZmumu;
  if ( templateQCD           != histogramQCD           ) delete templateQCD;
  if ( templateWplusJets     != histogramWplusJets     ) delete templateWplusJets;
  if ( templateTTplusJets    != histogramTTplusJets    ) delete templateTTplusJets;
  delete templateEWKbgSum;
  delete templateSMbgSum;
  delete templateSMsum;

  delete canvas;
}

//
//-------------------------------------------------------------------------------
//

void addToFormula(std::string& formula, const std::string& expression, TObjArray& arguments, RooAbsReal* p)
{
//-------------------------------------------------------------------------------
// Build formula and argument list for RooFormulaVar
//
// NOTE: auxiliary function for makeRooFormulaVar
//
//-------------------------------------------------------------------------------

  //std::cout << "<addToFormula>:" << std::endl;

  if ( expression != "" ) {
    if ( formula != "" ) formula.append("*");

    if      ( expression == "regular"  ) formula.append(p->GetName());
    else if ( expression == "inverted" ) formula.append("(1 - ").append(p->GetName()).append(")");
    else assert(0);
    
    arguments.Add(p);
  }
}

std::vector<std::string> getObservables(const std::string& region, const std::vector<std::string>& fitVariables)
{
  std::vector<std::string> retVal;

  if        ( region           == "ABCD"            ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
  } else if ( region.find("A") != std::string::npos ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else if ( region.find("B") != std::string::npos ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else if ( region.find("C") != std::string::npos ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else if ( region.find("D") != std::string::npos ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else {
    std::cout << "Error in <getObservables>: undefined region = " << region << " !!" << std::endl;
  }

  return retVal;
}

RooHistPdf* makeRooHistPdf(TH1* templateHistogram, RooAbsReal* fitVar)
{
  //std::cout << "<makeRooHistPdf>:" << std::endl;
  
  std::string templateDataHistName = std::string(templateHistogram->GetName()).append("_dataHist");
  RooDataHist* templateDataHist = new RooDataHist(templateDataHistName.data(), templateDataHistName.data(), *fitVar, templateHistogram);
  std::string templatePdfName = std::string(templateHistogram->GetName()).append("_histPdf");
  RooHistPdf* templatePdf = new RooHistPdf(templatePdfName.data(), templatePdfName.data(), *fitVar, *templateDataHist);
  return templatePdf;
}

RooGaussian* makeFitConstraint(RooAbsReal* p, double value, double error)
{
  //std::cout << "<makeFitConstraint>:" << std::endl;

  std::string pValueName = std::string(p->GetName()).append("_constValue");
  RooConstVar* pValue = new RooConstVar(pValueName.data(), pValueName.data(), value);
  std::string pErrorName = std::string(p->GetName()).append("_constError");
  RooConstVar* pError = new RooConstVar(pErrorName.data(), pErrorName.data(), error);

  std::string constraintName = std::string(p->GetName()).append("_constraint");
  RooGaussian* constraint = new RooGaussian(constraintName.data(), constraintName.data(), *p, *pValue, *pError);
  return constraint;
}

void loadHistograms(
  std::map<std::string, std::map<std::string, TH1*> >& histogramMap,
  TDirectory* inputDirectory, const std::string& process, const std::vector<std::string>& regions,
  const std::vector<std::string>& tauIds, const std::vector<std::string>& fitVariables, const std::string& sysShift)
{
//--------------------------------------------------------------------------------
// Load template histograms/distributions observed in data from ROOT file
//--------------------------------------------------------------------------------

  //std::cout << "<loadHistograms>:" << std::endl;

  for ( std::vector<std::string>::const_iterator region = regions.begin();
	region != regions.end(); ++region ) {
    
    std::vector<std::string> observables = getObservables(*region, fitVariables);
    std::vector<std::string> tauIdValues = getTauIdValues(*region);
    
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {
	for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {

	  std::string histogramName = std::string(process).append("_").append(*region).append("_").append(*observable);
	  histogramName.append("_").append(*tauId).append("_").append(*tauIdValue);
	  if ( sysShift != "CENTRAL_VALUE" ) histogramName.append("_").append(sysShift);
 
	  TH1* histogram = dynamic_cast<TH1*>(inputDirectory->Get(histogramName.data()));
	  if ( !histogram ) {
	    std::cout << "Error in <loadHistograms>: failed to load histogram = " << histogramName 
		      << " from file/directory = " << inputDirectory->GetName() << " --> aborting !!";
	    assert(0);
	  }

	  //std::cout << " histogram = " << histogram << std::endl;
	  //std::cout << " name      = " << histogram->GetName() << std::endl;
	  //std::cout << " integral  = " << histogram->Integral() << std::endl;
	  
	  int numBins = histogram->GetNbinsX();
	  if      ( (numBins % 3) == 0                  ) histogram->Rebin(3);
	  else if ( (numBins % 4) == 0 && numBins >= 36 ) histogram->Rebin(4);
	  else                                            histogram->Rebin(2);

	  // CV: scale MC histograms by 0.80 to account for crab jobs lost when processing Data
	  //    (temporary fix, 2011/07/11)
	  if ( process != "Data" ) {
	    if ( !histogram->GetSumw2N() ) histogram->Sumw2();	    
	    histogram->Scale(0.80);
	  }

	  std::string key = getKey(*observable, *tauId, *tauIdValue, sysShift);	
	  //std::cout << "--> key = " << key << std::endl;
	  if ( histogram != 0 ) histogramMap[*region][key] = histogram;
	}
      }
    }
  }
}

struct sysUncertaintyEntry
{
  sysUncertaintyEntry(const std::string& sysNameUp, const std::string& sysNameDown, const std::string& sysNameDiff)
    : sysNameUp_(sysNameUp),
      sysNameDown_(sysNameDown),
      sysNameDiff_(sysNameDiff)
  {}
  ~sysUncertaintyEntry() {}
  std::string sysNameUp_;
  std::string sysNameDown_;
  std::string sysNameDiff_;
};

void compSysHistograms(std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,
		       const sysUncertaintyEntry& sysUncertainty)
{
  for ( std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >::iterator process = templatesAll.begin();
	process != templatesAll.end(); ++process ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::iterator region = process->second.begin();
	  region != process->second.end(); ++region ) {
      for ( std::map<std::string, TH1*>::iterator keyUp = region->second.begin();
	    keyUp != region->second.end(); ++keyUp ) {
	if ( keyUp->first.find(sysUncertainty.sysNameUp_) != std::string::npos ) {
	  //std::cout << "keyUp = " << keyUp->first << std::endl;
	  TH1* sysHistogramUp = keyUp->second;
	  std::string keyDown = TString(keyUp->first).ReplaceAll(sysUncertainty.sysNameUp_.data(), 
								 sysUncertainty.sysNameDown_.data()).Data();	  
	  //std::cout << "keyDown = " << keyDown << std::endl;
	  TH1* sysHistogramDown = region->second[keyDown];
	  assert(sysHistogramDown);

	  assert(isCompatibleBinning(sysHistogramUp, sysHistogramDown));
  
	  if ( !sysHistogramUp->GetSumw2N()   ) sysHistogramUp->Sumw2();
	  if ( !sysHistogramDown->GetSumw2N() ) sysHistogramDown->Sumw2();

	  TH1* sysHistogramDiff = (TH1*)sysHistogramUp->Clone();

	  unsigned numBins = sysHistogramUp->GetNbinsX();
	  for ( unsigned iBin = 0; iBin <= (numBins + 1); ++iBin ) {
	    double binContentDiff = 0.5*(sysHistogramUp->GetBinContent(iBin) - sysHistogramDown->GetBinContent(iBin));
	    sysHistogramDiff->SetBinContent(iBin, binContentDiff);
	    
	    double binErrorUp   = sysHistogramUp->GetBinError(iBin);
	    double binErrorDown = sysHistogramDown->GetBinError(iBin);
	    double binErrorDiff = TMath::Sqrt(binErrorUp*binErrorUp + binErrorDown*binErrorDown);
	    sysHistogramDiff->SetBinError(iBin, binErrorDiff);
	  }

	  std::string keyDiff = TString(keyUp->first).ReplaceAll(sysUncertainty.sysNameUp_.data(), 
								 sysUncertainty.sysNameDiff_.data()).Data();
	  //std::cout << "keyDiff = " << keyDiff << std::endl;
	  region->second[keyDiff] = sysHistogramDiff;
	}
      }	    
    }
  }
} 

void sumHistograms(std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,
		   const std::vector<std::string>& processesToSum, const std::string& processNameSum,
		   double mcToDataScaleFactor = 1.0)
{
//--------------------------------------------------------------------------------
// Compute total Standard Model expectation by summing all signal and background contributions
//--------------------------------------------------------------------------------

  for ( std::vector<std::string>::const_iterator process = processesToSum.begin();
	process != processesToSum.end(); ++process ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = templatesAll["Ztautau"].begin();
	  region != templatesAll["Ztautau"].end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	TH1* histogram = templatesAll[*process][region->first][key->first];

	TH1* histogramSum = templatesAll[processNameSum][region->first][key->first];
	if ( !histogramSum ) {
	  std::string histogramName = histogram->GetName();
	  std::string histogramSumName = std::string(processNameSum).append(std::string(histogramName, histogramName.find("_")));
	  histogramSum = (TH1*)histogram->Clone(histogramSumName.data());
	  histogramSum->Scale(mcToDataScaleFactor);
	  templatesAll[processNameSum][region->first][key->first] = histogramSum;
	} else {	    	      
	  histogramSum->Add(histogram, mcToDataScaleFactor);
	}
      }
    }
  }
}

void scaleBins(std::map<std::string, std::map<std::string, TH1*> >& histograms, const std::string& region, 
	       const std::string& observable, double xMin, double xMax, double scaleFactor,
	       const std::vector<std::string>& tauIds)
{
//--------------------------------------------------------------------------------
// Scale bins of distribution/template that are within the range xMin..xMax
// up/down by scale-factor given as function argument
//
// NOTE: to be used for systematic uncertainty studies only !!
//
//--------------------------------------------------------------------------------

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {

    std::vector<std::string> keys;
    keys.push_back(getKey(observable, *tauId));
    keys.push_back(getKey(observable, *tauId, "passed"));
    keys.push_back(getKey(observable, *tauId, "failed"));

    for ( std::vector<std::string>::const_iterator key = keys.begin();
	  key != keys.end(); ++key ) {
      if ( histograms[region].find(*key) == histograms[region].end() ) continue;
      
      TH1* histogram = histograms[region][*key];
      
      int numBins = histogram->GetNbinsX();
      int numBins_scaled = 0;
      for ( int iBin = 1; iBin <= numBins; ++iBin ) {
	double binCenter  = histogram->GetBinCenter(iBin);
	if ( binCenter > xMin && binCenter < xMax ) {
	  double binContent = histogram->GetBinContent(iBin);
	  histogram->SetBinContent(iBin, scaleFactor*binContent);
	  ++numBins_scaled;
	}
      }

      std::cout << "histogram: name = " << histogram->GetName() << ", #scaled bins = " << numBins_scaled << std::endl;
      
//--- in case scale-factor is zero, merge bins set to zero with adjacent (non-zero) bin
//    in order to avoid empty bins (which might confuse RooFit)
      if ( scaleFactor == 0. ) {
	TAxis* axis = histogram->GetXaxis();

	double leftEdge_mergedBin  = 0;
	double rightEdge_mergedBin = 0;
	const double epsilon = 1.e-3;
	if        ( axis->FindBin(xMax) < numBins ) {
	  leftEdge_mergedBin  = axis->GetBinLowEdge(axis->FindBin(xMin));
          rightEdge_mergedBin = axis->GetBinUpEdge(axis->FindBin(xMax - epsilon) + 1);
	} else if ( axis->FindBin(xMin) > 1       ) {
	  leftEdge_mergedBin  = axis->GetBinLowEdge(axis->FindBin(xMin) - 1);
          rightEdge_mergedBin = axis->GetBinUpEdge(axis->FindBin(xMax - epsilon));
	} else assert(0);

	std::cout << "leftEdge_mergedBin = " << leftEdge_mergedBin << std::endl;
	std::cout << "rightEdge_mergedBin = " << rightEdge_mergedBin << std::endl;

	int numBins_rebinned = numBins - (numBins_scaled + 1);
	double* binEdges_rebinned = new double[numBins_rebinned + 1];
	int iBin_rebinned = 0;
	for ( int iBin = 1; iBin <= numBins; ++iBin ) {
	  double leftEdge  = axis->GetBinLowEdge(iBin);
	  double rightEdge = axis->GetBinUpEdge(iBin);

	  std::cout << "iBin = " << iBin << ": leftEdge = " << leftEdge << ", rightEdge = " << rightEdge << std::endl;

	  if        ( rightEdge <= leftEdge_mergedBin || 
		      leftEdge  >= rightEdge_mergedBin ) {
	    ++iBin_rebinned;
	    binEdges_rebinned[iBin_rebinned - 1] = leftEdge;
	  } else if ( leftEdge  == leftEdge_mergedBin  ) {
	    ++iBin_rebinned;
	    binEdges_rebinned[iBin_rebinned - 1] = leftEdge;
	  } 
	}

	std::cout << "numBins_rebinned = " << numBins_rebinned << std::endl;
	std::cout << "iBin_rebinned = " << iBin_rebinned << std::endl;

	assert(iBin_rebinned == (numBins_rebinned + 1));
	binEdges_rebinned[numBins_rebinned] = axis->GetXmax();

	std::cout << "numBins = " << numBins << " --> numBins_rebinned = " << numBins_rebinned << std::endl;
	std::cout << "binEdges = { ";
	for ( int iBin = 1; iBin <= (numBins_rebinned + 1); ++iBin ) {
	  std::cout << binEdges_rebinned[iBin - 1];
	  if ( iBin <= numBins_rebinned ) std::cout << ", ";
	}
	std::cout << " }" << std::endl;
	
	std::string histogramName_rebinned = std::string(histogram->GetName()).append("_rebinned");
	histograms[region][*key] = histogram->Rebin(numBins_rebinned, histogramName_rebinned.data(), binEdges_rebinned);
	delete binEdges_rebinned;
      }
    }
  }
}

//
//-------------------------------------------------------------------------------
//

std::map<std::string, std::map<std::string, std::map<std::string, double> > > compNumEvents(
  std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,
  const std::vector<std::string>& processes, const std::map<std::string, std::map<std::string, TH1*> >& distributionsData)
{
  std::map<std::string, std::map<std::string, std::map<std::string, double> > > retVal; // key = (process/"sum", region, observable)

  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	retVal[*process][region->first][key->first] = getIntegral(templatesAll[*process][region->first][key->first], true, true);
	//std::cout << "numEvents[" << (*process) << "][" << region->first << "][" << key->first << "] = "
	//          << retVal[*process][region->first][key->first] << std::endl;
	retVal["sum"][region->first][key->first] += retVal[*process][region->first][key->first];
      }
    }
  }

  std::cout << std::endl;

  return retVal;
}

std::map<std::string, std::map<std::string, std::map<std::string, double> > > compFittedFractions(
  std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,
  std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll,
  const std::vector<std::string>& processes, const std::map<std::string, std::map<std::string, TH1*> >& distributionsData)
{
  std::map<std::string, std::map<std::string, std::map<std::string, double> > > retVal; // key = (process, region, observable)
  
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	retVal[*process][region->first][key->first] = 
	  getIntegral(templatesAll[*process][region->first][key->first], false, false)/numEventsAll[*process][region->first][key->first];
	//std::cout << "fittedFraction[" << (*process) << "][" << region->first << "][" << key->first << "] = "
	//	    << retVal[*process][region->first][key->first] << std::endl;
      }
    }
  }

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

void savePseudoExperimentHistograms(const std::map<std::string, std::map<std::string, TH1*> >& histograms, 
				    const std::string& xAxisLabel, const std::string& outputFileNameSuffix)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator tauId = histograms.begin();
	tauId != histograms.end(); ++tauId ) {
    for ( std::map<std::string, TH1*>::const_iterator fitVariable = tauId->second.begin();
	  fitVariable != tauId->second.end(); ++fitVariable ) {
      TH1* histogram = fitVariable->second;

      TH1* histogram_normalized = normalize(histogram);
      
      applyStyleOption(histogram_normalized, histogram->GetTitle(), xAxisLabel, "a.u");

      canvas->Update();
      std::string outputFileName = std::string(histogram->GetName()).append(outputFileNameSuffix);
      std::string outputFilePath = std::string("./plots/");
      gSystem->mkdir(outputFilePath.data(), true);
      canvas->Print(outputFilePath.append(outputFileName).data());
    }
  }

  delete canvas;
}

//
//-------------------------------------------------------------------------------
//

std::pair<double, double> getNumber(TDirectory* inputDirectory, const TString& auxHistogramName, 
				    int auxHistogramBin, int assertAuxHistogramNumBins = 1)
{
  //std::cout << "<getNumber>:" << std::endl;
  //std::cout << " auxHistogramName = " << auxHistogramName << std::endl;
  //std::cout << " auxHistogramBin = " << auxHistogramBin << std::endl;
  //std::cout << " assertAuxHistogramNumBins = " << assertAuxHistogramNumBins << std::endl;

  TH1* histogram = dynamic_cast<TH1*>(inputDirectory->Get(auxHistogramName.Data()));
  if ( !histogram ) 
    throw cms::Exception("getNumber")  
      << "Failed to find histogram = " << auxHistogramName << " in input file/directory = " << inputDirectory->GetName() << " !!\n";
  
  int numBins = histogram->GetNbinsX();
  // CV: check that histogram has the expected number of bins,
  //     if it hasn't, FWLiteTauIdEffPreselNumbers/TauIdEffCutFlowTable has probably changed
  //     and the auxHistogramBin parameter needs to be updated !!
  //    (in particular if auxHistogramBin is negative)
  if ( numBins != assertAuxHistogramNumBins ) 
    throw cms::Exception("getNumber")  
      << "Histogram = " << auxHistogramName << " has incompatible binning:" 
      << " found = " << numBins << ", expected = " << assertAuxHistogramNumBins << " !!\n";

  int x;
  if   ( auxHistogramBin >= 0 ) x = auxHistogramBin;
  else                          x = numBins - TMath::Abs(auxHistogramBin);
  if ( !(x >= 0 && x < numBins) ) 
    throw cms::Exception("getNumber") 
      << "Invalid auxHistogramBin = " << auxHistogramBin << ", expected range = " << (-numBins) << ".." << (numBins - 1) << " !!\n";

  int bin = histogram->FindBin(x);

  std::pair<double, double> retVal;
  retVal.first = histogram->GetBinContent(bin);
  retVal.second = histogram->GetBinError(bin);

  //std::cout << "--> returning " << retVal.first << " +/- " << retVal.second << std::endl;

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

void printTimeStamp()
{
  time_t now_raw;
  time(&now_raw);
  struct tm* now_local = localtime(&now_raw);
  std::cout << "it is NOW: " << asctime(now_local) << std::endl;
}

#endif
