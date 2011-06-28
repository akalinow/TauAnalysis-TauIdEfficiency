#ifndef TauAnalysis_GenSimTools_tauIdEffAuxFunctions_h
#define TauAnalysis_GenSimTools_tauIdEffAuxFunctions_h

#include "TauAnalysis/DQMTools/interface/histogramAuxFunctions.h"

#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
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

void applyStyleOption(TH1* histogram, const std::string& histogramTitle,
		      const std::string& xAxisTitle, const std::string& yAxisTitle = "Number of Events")
{
  //std::cout << "<applyStyleOption>:" << std::endl;

  histogram->SetStats(false);

  histogram->SetTitle(histogramTitle.data());
  histogram->SetTitle("");
    
  histogram->GetXaxis()->SetTitle(xAxisTitle.data());
  histogram->GetXaxis()->SetTitleOffset(1.15);
  //histogram->GetXaxis()->SetTitleSize(0.05); 
  //histogram->GetXaxis()->SetLabelSize(0.05);

  histogram->GetYaxis()->SetTitle(yAxisTitle.data());
  histogram->GetYaxis()->SetTitleOffset(1.20);
  //histogram->GetYaxis()->SetTitleSize(0.05); 
  //histogram->GetYaxis()->SetLabelSize(0.05);
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

  //std::cout << "<drawHistograms>:" << std::endl;
  //std::cout << " Ztautau:       histogram = " << histogramZtautau       << ", norm = " << normZtautau       << std::endl;
  //std::cout << " Zmumu:         histogram = " << histogramZmumu         << ", norm = " << normZmumu         << std::endl;
  //std::cout << " QCD:           histogram = " << histogramQCD           << ", norm = " << normQCD           << std::endl;
  //std::cout << " WplusJets:     histogram = " << histogramWplusJets     << ", norm = " << normWplusJets     << std::endl;
  //std::cout << " TTplusJets:    histogram = " << histogramTTplusJets    << ", norm = " << normTTplusJets    << std::endl;
  //std::cout << " Data:          histogram = " << histogramData          << std::endl;
    
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetLeftMargin(0.12);
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

  smSum.SetTitle(templateZtautau->GetTitle());
  smSum.SetMaximum(1.4*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));

  smSum.Draw("hist");
  histogramData->SetStats(false);
  histogramData->Draw("ep1same");
  //smSum.Draw("axissame");

  TLegend legend(0.64, 0.64, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  legend.AddEntry(histogramData,      "Data",                                       "p");
  legend.AddEntry(templateZtautau,    "Z/#gamma^{*} #rightarrow #tau^{+} #tau^{-}", "f");
  legend.AddEntry(templateQCD,        "QCD",                                        "f");
  legend.AddEntry(templateWplusJets,  "W + jets",                                   "f");
  legend.AddEntry(templateZmumu,      "Z/#gamma^{*} #rightarrow #mu^{+} #mu^{-}",   "f");
  legend.AddEntry(templateTTplusJets, "t#bar{t} + jets",                            "f");
  legend.Draw();

  TPaveText cmsPreliminaryLabel(0.140, 0.860, 0.46, 0.900, "NDC");
  cmsPreliminaryLabel.AddText("CMS Preliminary 2011");
  cmsPreliminaryLabel.SetTextAlign(13);
  cmsPreliminaryLabel.SetTextSize(0.045);
  cmsPreliminaryLabel.SetFillStyle(0);
  cmsPreliminaryLabel.SetBorderSize(0);
  cmsPreliminaryLabel.Draw();

  TPaveText cmsLuminosityLabel(0.145, 0.8075, 0.460, 0.8475, "NDC");
  cmsLuminosityLabel.AddText("#sqrt{s} = 7 TeV, L = 191 pb^{-1}");
  cmsLuminosityLabel.SetTextAlign(13);
  cmsLuminosityLabel.SetTextSize(0.045);
  cmsLuminosityLabel.SetFillStyle(0);
  cmsLuminosityLabel.SetBorderSize(0);
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  gSystem->mkdir(outputFilePath.data(), true);
  canvas->Print(outputFilePath.append(outputFileName).data());

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

  THStack smSum2("smSum", "smSum");
  smSum2.Add(templateSMbgSum);
  smSum2.Add(templateZtautau);

  smSum2.SetTitle(templateZtautau->GetTitle());
  smSum2.SetMaximum(1.4*TMath::Max(smSum2.GetMaximum(), histogramData->GetMaximum()));
	
  smSum2.Draw("hist");
  histogramData->Draw("ep1same");
  //smSum2.Draw("axissame");

  TLegend legend2(0.64, 0.74, 0.89, 0.89, "", "brNDC"); 
  legend2.SetBorderSize(0);
  legend2.SetFillColor(0);
  
  legend2.AddEntry(histogramData,   "Data",                            "p");
  legend2.AddEntry(templateZtautau, "Z #rightarrow #tau^{+} #tau^{-}", "f");
  legend2.AddEntry(templateSMbgSum, "#Sigma Backgrounds",              "f");
  legend2.Draw();

  cmsPreliminaryLabel.Draw();
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath2 = std::string("./plots/");
  gSystem->mkdir(outputFilePath2.data(), true);
  TString outputFileName2 = outputFileName;
  outputFileName2.ReplaceAll(".", "_smBgSum.");
  canvas->Print(outputFilePath2.append(outputFileName2).data());

//--- draw histograms for all signal and background processes summed
  TH1* templateSMsum = (TH1*)templateSMbgSum->Clone();
  templateSMsum->Add(templateZtautau);
  applyStyleOption(templateSMsum, histogramTitle, xAxisTitle);
  templateSMsum->SetFillStyle(1001);
  templateSMsum->SetFillColor(10);
  templateSMsum->SetLineColor(1);
  templateSMsum->SetLineWidth(2);

  templateSMsum->SetTitle(templateZtautau->GetTitle());
  templateSMsum->SetMaximum(1.4*TMath::Max(templateSMsum->GetMaximum(), histogramData->GetMaximum()));
	
  templateSMsum->Draw("hist");
  histogramData->Draw("ep1same");
  //templateSMsum->Draw("axissame");

  TLegend legend3(0.64, 0.79, 0.89, 0.89, "", "brNDC"); 
  legend3.SetBorderSize(0);
  legend3.SetFillColor(0);
  
  legend3.AddEntry(histogramData, "Data", "p");
  legend3.AddEntry(templateSMsum, "MC",   "l");
  legend3.Draw();
  
  cmsPreliminaryLabel.Draw();
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath3 = std::string("./plots/");
  gSystem->mkdir(outputFilePath3.data(), true);
  TString outputFileName3 = outputFileName;
  outputFileName3.ReplaceAll(".", "_smSum.");
  canvas->Print(outputFilePath3.append(outputFileName3).data());  

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
  delete templateSMbgSum;
  delete templateSMsum;

  delete canvas;
}

void drawHistograms(std::map<std::string, std::map<std::string, TH1*> >& distributionsData, 
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,    
		    std::map<std::string, double> normFactors,                                           
		    const std::string& region, const std::string& observable_key,
		    const std::string& histogramTitle, const std::string& xAxisTitle, 
		    const std::string& outputFileName,
		    bool runStatTest = false)
{
  //std::cout << "<drawHistograms (wrapper)>:" << std::endl;

  drawHistograms(templatesAll["Ztautau"][region][observable_key], normFactors["Ztautau"],
		 templatesAll["Zmumu"][region][observable_key], normFactors["Zmumu"],
		 templatesAll["QCD"][region][observable_key], normFactors["QCD"],
		 templatesAll["WplusJets"][region][observable_key], normFactors["WplusJets"],
		 templatesAll["TTplusJets"][region][observable_key], normFactors["TTplusJets"],
		 distributionsData[region][observable_key],
		 histogramTitle, xAxisTitle,
		 outputFileName,
		 runStatTest);
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

std::vector<std::string> getTauIdValues(const std::string& region)
{
  std::vector<std::string> retVal;

  if      ( region.find("p") != std::string::npos ) retVal.push_back(std::string("passed"));
  else if ( region.find("f") != std::string::npos ) retVal.push_back(std::string("failed"));
  else                                              retVal.push_back(std::string("all"));

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
  TFile* inputFile, const std::string& process, const std::vector<std::string>& regions,
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
 
	  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
	  if ( !histogram ) {
	    std::cout << "Error in <loadHistograms>: failed to load histogram = " << histogramName 
		      << " from file = " << inputFile->GetName() << " --> aborting !!";
	    assert(0);
	  }

	  //std::cout << " histogram = " << histogram << std::endl;
	  //std::cout << " name      = " << histogram->GetName() << std::endl;
	  //std::cout << " integral  = " << histogram->Integral() << std::endl;
	  
	  int numBins = histogram->GetNbinsX();
	  if      ( (numBins % 3) == 0                  ) histogram->Rebin(3);
	  else if ( (numBins % 4) == 0 && numBins >= 36 ) histogram->Rebin(4);
	  else                                            histogram->Rebin(2);

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
      TH1* sysHistogramUp   = 0;
      TH1* sysHistogramDown = 0;
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
	std::cout << "numEvents[" << (*process) << "][" << region->first << "][" << key->first << "] = "
		  << retVal[*process][region->first][key->first] << std::endl;
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

#endif
