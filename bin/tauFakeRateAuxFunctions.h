#ifndef TauAnalysis_TauIdEfficiency_tauFakeRateAuxFunctions_h
#define TauAnalysis_TauIdEfficiency_tauFakeRateAuxFunctions_h

#include "TauAnalysis/DQMTools/interface/histogramAuxFunctions.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
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

typedef std::map<std::string, TH1*>               histogramMapType1;
typedef std::map<std::string, histogramMapType1>  histogramMapType2;
typedef std::map<std::string, histogramMapType2>  histogramMapType3;
typedef std::map<std::string, histogramMapType3>  histogramMapType4;
typedef std::map<std::string, histogramMapType4>  histogramMapType5;

typedef std::map<std::string, TGraphAsymmErrors*> frMapType1;
typedef std::map<std::string, frMapType1>         frMapType2; 
typedef std::map<std::string, frMapType2>         frMapType3; 
typedef std::map<std::string, frMapType3>         frMapType4; 

typedef std::vector<std::string> vstring;

int defaultMarkerColor =    1;
int defaultMarkerStyle =   20;
int defaultMarkerSize  =    1.;

int defaultLineColor   =    1;
int defaultLineStyle   =    1;
int defaultLineSize    =    1;  

int defaultFillColor   =    5;
int defaultFillStyle   = 1001;

double square(double x)
{
  return x*x;
}

// define auxiliary class holding draw-options for control plots
struct histogramDrawOptionType 
{
  histogramDrawOptionType()
  {}
  histogramDrawOptionType(const edm::ParameterSet& cfg)
    : name_(cfg.getParameter<std::string>("name")),
      legendEntry_(cfg.getParameter<std::string>("legendEntry")),
      
      //drawOption_(cfg.getParameter<std::string>("drawOption")),
      drawOptionLegend_(cfg.getParameter<std::string>("drawOptionLegend"))
  {
    markerColor_ = ( cfg.exists("markerColor") ) ? cfg.getParameter<int>("markerColor")   : defaultMarkerColor;
    markerStyle_ = ( cfg.exists("markerStyle") ) ? cfg.getParameter<int>("markerStyle")   : defaultMarkerStyle;
    markerSize_  = ( cfg.exists("markerSize")  ) ? cfg.getParameter<double>("markerSize") : defaultMarkerSize;
      
    lineColor_   = ( cfg.exists("lineColor")   ) ? cfg.getParameter<int>("lineColor")     : defaultLineColor;
    lineStyle_   = ( cfg.exists("lineStyle")   ) ? cfg.getParameter<int>("lineStyle")     : defaultLineStyle;
    lineWidth_   = ( cfg.exists("lineWidth")   ) ? cfg.getParameter<int>("lineWidth")     : defaultLineSize;

    fillColor_   = ( cfg.exists("fillColor")   ) ? cfg.getParameter<int>("fillColor")     : defaultFillColor;
    fillStyle_   = ( cfg.exists("fillStyle")   ) ? cfg.getParameter<int>("fillStyle")     : defaultFillStyle;
  }
  ~histogramDrawOptionType() {}
  std::string name_;
  std::string legendEntry_;
  int markerColor_;
  int markerStyle_;
  double markerSize_;
  int lineColor_;
  int lineStyle_;
  int lineWidth_;
  int fillColor_;
  int fillStyle_;
  //std::string drawOption_;
  std::string drawOptionLegend_;
};

// define auxiliary class holding draw-options for fake-rate plots
struct graphDrawOptionType
{
  graphDrawOptionType()
  {}
  graphDrawOptionType(const edm::ParameterSet& cfg)
    : name_(cfg.getParameter<std::string>("name")),
      legendEntry_(cfg.getParameter<std::string>("legendEntry")),
      expMarkerStyle_(cfg.getParameter<unsigned>("markerStyleSim")),
      measMarkerStyle_(cfg.getParameter<unsigned>("markerStyleData")),
      markerSize_(cfg.getParameter<unsigned>("markerSize")),
      color_(cfg.getParameter<unsigned>("color"))
  {}
  ~graphDrawOptionType() {}
  std::string name_;
  std::string legendEntry_;
  unsigned expMarkerStyle_;
  unsigned measMarkerStyle_;
  unsigned markerSize_;
  unsigned color_;
};

void applyStyleOption_histogram(
  TH1* histogram, const std::string& histogramTitle, const histogramDrawOptionType& drawOptions, 
  const std::string& xAxisTitle, const std::string& yAxisTitle = "Events")
{
  //std::cout << "<applyStyleOption_histogram>:" << std::endl;
  //std::cout << " histogramName = " << histogram->GetName() << std::endl;

  histogram->SetTitle(histogramTitle.data());

  histogram->SetMarkerColor(drawOptions.markerColor_);
  histogram->SetMarkerStyle(drawOptions.markerStyle_);
  histogram->SetMarkerSize(drawOptions.markerSize_);

  histogram->SetLineColor(drawOptions.lineColor_);
  histogram->SetLineStyle(drawOptions.lineStyle_);
  histogram->SetLineWidth(drawOptions.lineWidth_);
  
  histogram->SetFillColor(drawOptions.fillColor_);
  histogram->SetFillStyle(drawOptions.fillStyle_);

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

void applyStyleOption_histogram(
  THStack* histogram, const std::string& histogramTitle,
  const std::string& xAxisTitle, const std::string& yAxisTitle = "Events")
{
  //std::cout << "<applyStyleOption_histogram>:" << std::endl;
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
  
void applyStyleOption_graph(
  TGraph* graph, const std::string& graphTitle,
  unsigned markerStyle, unsigned markerSize, unsigned color,
  const std::string& xAxisTitle, const std::string& yAxisTitle = "Fake-rate")
{
  //std::cout << "<applyStyleOption_graph>:" << std::endl;
  //std::cout << " graphName = " << graph->GetName() << std::endl;

  graph->SetTitle(graphTitle.data());

  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(markerSize);
  graph->SetMarkerColor(color);

  graph->SetLineColor(color);

  if ( graph->GetXaxis() ) {
    graph->GetXaxis()->SetTitle(xAxisTitle.data());
    graph->GetXaxis()->SetTitleOffset(1.15);
    //graph->GetXaxis()->SetTitleSize(0.05); 
    //graph->GetXaxis()->SetLabelSize(0.05);
  } //else std::cerr << "Graph = " << graph->GetName() << " has no valid x-Axis !!" << std::endl;

  if ( graph->GetYaxis() ) {
    graph->GetYaxis()->SetTitle(yAxisTitle.data());
    graph->GetYaxis()->SetTitleOffset(1.65);
    //graph->GetYaxis()->SetTitleSize(0.05); 
    //graph->GetYaxis()->SetLabelSize(0.05);
  } //else std::cerr << "Graph = " << graph->GetName() << " has no valid y-Axis !!" << std::endl;
}

//
//-------------------------------------------------------------------------------
//

void loadHistograms(
  histogramMapType5& histogramMap,
  TFile* inputFile, 
  const vstring& processes, const vstring& eventSelections, 
  const std::map<std::string, double>& avTriggerPrescales, 
  const vstring& tauIds, const vstring& observables, const vstring& regions)
{
//--------------------------------------------------------------------------------
// Load histograms from ROOT file
//--------------------------------------------------------------------------------

  //std::cout << "<loadHistograms>:" << std::endl;
  
  for ( vstring::const_iterator eventSelection = eventSelections.begin();
	eventSelection != eventSelections.end(); ++eventSelection ) {
    std::string inputDirectoryName = (*eventSelection);
    TDirectory* inputDirectory = dynamic_cast<TDirectory*>(inputFile->Get(inputDirectoryName.data()));
    if ( !inputDirectory ) 
      throw cms::Exception("loadHistograms") 
	<< "Failed to load histograms for event selection = " << (*eventSelection) << "," 
	<< " directory = " << inputDirectoryName << " does not exists in input file = " << inputFile->GetName() << " !!\n";
    
    if ( avTriggerPrescales.find(*eventSelection) == avTriggerPrescales.end() ) 
      throw cms::Exception("loadHistograms") 
	<< "No trigger prescale defined for event selection = " << (*eventSelection) << " !!\n";
    double avTriggerPrescale = avTriggerPrescales.find(*eventSelection)->second;
    assert(avTriggerPrescale > 0.);
    
    for ( vstring::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      for ( vstring::const_iterator tauId = tauIds.begin();
	    tauId != tauIds.end(); ++tauId ) {
	for ( vstring::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {
	  for ( vstring::const_iterator region = regions.begin();
		region != regions.end(); ++region ) {
	    std::string histogramName = std::string(*process).append("_").append(*region).append("_").append(*observable);
	    histogramName.append("_").append(*tauId);
	    
	    TH1* histogram = dynamic_cast<TH1*>(inputDirectory->Get(histogramName.data()));
	    if ( !histogram ) {
	      throw cms::Exception("loadHistograms") 
		<< "Failed to load histogram = " << histogramName << " from directory = " << inputDirectory->GetName() << " !!\n";
	      //std::cerr << "Failed to load histogram = " << histogramName 
	      //	  << " from directory = " << inputDirectory->GetName() << " !!\n";
	      //continue;
	    }
	    
	    if ( !histogram->GetSumw2N() ) histogram->Sumw2();

	    if ( process->find("Data") == std::string::npos ) histogram->Scale(1./avTriggerPrescale);

	    histogramMap[*eventSelection][*process][*tauId][*observable][*region] = histogram;
	  }
	}
      }
    }
  }
}

//
//-------------------------------------------------------------------------------
//

void getBinomialBounds(Int_t n, Int_t r, Float_t& rMin, Float_t& rMax)
{
  //std::cout << "<getBinomialBounds>:" << std::endl;
  //std::cout << " n = " << n << ", r = " << r << std::endl;

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

  Double_t rMinLeft       = 0.;
  Double_t rMinMiddle     = TMath::Max(0.5*r, n - 1.5*r);
  Double_t rMinRight      = n;
  Double_t rMinLeftProb   = 0.;
  Double_t rMinMiddleProb = 0.5;
  Double_t rMinRightProb  = 1.;
  while ( (rMinRight - rMinLeft) > (0.001*n) ){

    rMinMiddleProb = TMath::BinomialI(rMinMiddle/n, n, r);

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
      rMinMiddle = rMinLeft = rMinRight = 0.;
    }
  }

  Double_t rMaxLeft       = 0.;
  Double_t rMaxMiddle     = TMath::Min(1.5*r, n - 0.5*r);
  Double_t rMaxRight      = n;
  Double_t rMaxLeftProb   = 1.;
  Double_t rMaxMiddleProb = 0.5;
  Double_t rMaxRightProb  = 0.;
  while ( (rMaxRight - rMaxLeft) > (0.001*n) ){

    rMaxMiddleProb = 1. - TMath::BinomialI(rMaxMiddle/n, n, r);

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

TGraphAsymmErrors* getEfficiency(const TH1* histogram_numerator, const TH1* histogram_denominator,
				 bool applyEffStatCorrection = true)
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

  for ( Int_t iBin = 1; iBin <= nBins; iBin++ ){
    Float_t nObs_float = histogram_denominator->GetBinContent(iBin);
    Float_t rObs_float = histogram_numerator->GetBinContent(iBin);
    
    if ( applyEffStatCorrection && histogram_denominator->GetBinContent(iBin) > 0. && histogram_denominator->GetBinError(iBin) > 0. ) {
      Float_t effNum = square(histogram_denominator->GetBinContent(iBin)/histogram_denominator->GetBinError(iBin));
      Float_t scaleFactor = effNum/histogram_denominator->GetBinContent(iBin);
      nObs_float *= scaleFactor;
      rObs_float *= scaleFactor;
    }

    Int_t nObs = TMath::Nint(nObs_float);
    Int_t rObs = TMath::Nint(rObs_float);

    Float_t xCenter = histogram_denominator->GetBinCenter(iBin);
    Float_t xWidth  = histogram_denominator->GetBinWidth(iBin);

    x[iBin - 1]      = xCenter;
    dxUp[iBin - 1]   = 0.5*xWidth;
    dxDown[iBin - 1] = 0.5*xWidth;
    
    if ( nObs > 0 ){
      Float_t rMin = 0.;
      Float_t rMax = 0.;
      
      getBinomialBounds(nObs, rObs, rMin, rMax);

      y[iBin - 1]      = rObs/((Float_t)nObs);
      dyUp[iBin - 1]   = (rMax - rObs)/((Float_t)nObs);
      dyDown[iBin - 1] = (rObs - rMin)/((Float_t)nObs);
    } else{
      y[iBin - 1]      = 0.;
      dyUp[iBin - 1]   = 0.;
      dyDown[iBin - 1] = 0.;
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

void compFakeRates(
  histogramMapType5& histogramMap,
  frMapType4& frMap,
  const std::string& region_passed, const std::string& region_failed)
{
  for ( histogramMapType5::iterator eventSelection = histogramMap.begin();
	eventSelection != histogramMap.end(); ++eventSelection ) {
    for ( histogramMapType4::iterator process = eventSelection->second.begin();
	  process != eventSelection->second.end(); ++process ) {
      for ( histogramMapType3::iterator tauId = process->second.begin();
	    tauId != process->second.end(); ++tauId ) {
	for ( histogramMapType2::iterator observable = tauId->second.begin();
	      observable != tauId->second.end(); ++observable ) {
	  TH1* histogram_passed = observable->second[region_passed];
	  if ( !histogram_passed ) 
	    throw cms::Exception("compFakeRates") 
	      << "Failed to find histogram for 'passed' region !!\n";
	  
	  TH1* histogram_failed = observable->second[region_failed];
	  if ( !histogram_failed ) 
	    throw cms::Exception("compFakeRates") 
	      << "Failed to find histogram for 'failed' region !!\n";
	  
	  if ( !isCompatibleBinning(histogram_passed, histogram_failed) )
	    throw cms::Exception("compFakeRates") 
	      << "Incompatible binning of histograms for 'passed' = " << histogram_passed->GetName() 
	      << " and 'failed' = " << histogram_failed->GetName() << " region !!\n";

//--- compute efficiency graph from numerator and denominator histograms.
	  std::string numeratorName = std::string(process->first).append("_numerator_");
	  numeratorName.append(observable->first).append("_").append(tauId->first);
	  TH1* numerator = (TH1*)histogram_passed->Clone(numeratorName.data());
	  
	  std::string denominatorName = std::string(process->first).append("_denominator_");
	  denominatorName.append(observable->first).append("_").append(tauId->first);
	  TH1* denominator = (TH1*)histogram_failed->Clone(denominatorName.data());
	  denominator->Add(numerator);
//
//    CV: Make numerator and denominator histograms integer
//        before computing the tau id. efficiency/fake-rate,
//        in order to work-around the problem
//        that the TGraphAsymmErrors(TH1*, TH1*) constructor does not work in ROOT 5.27/06b (default in CMSSW_4_1_3),
//        in case the histograms contain weighted entries.
//
//        The binContents of numerator and denominator histograms are scaled
//        to the "effective" number of (denominator) entries:
//          effNum = (sum(weights)/sqrt(sum(weights^2)))^2 = (binContent/binError)^2	  
//
	  for ( int iBin = 1; iBin <= (numerator->GetNbinsX() + 1); ++iBin ) {
	    double scaleFactor = 0.;
	    if ( denominator->GetBinContent(iBin) > 0. && denominator->GetBinError(iBin) > 0. ) {
	      double effNum = square(denominator->GetBinContent(iBin)/denominator->GetBinError(iBin));
	      scaleFactor = effNum/denominator->GetBinContent(iBin);
	    }
	    
	    numerator->SetBinContent(iBin, TMath::Nint(scaleFactor*numerator->GetBinContent(iBin)));
	    numerator->SetBinError(iBin, TMath::Sqrt(numerator->GetBinContent(iBin)));
	    
	    denominator->SetBinContent(iBin, TMath::Nint(scaleFactor*denominator->GetBinContent(iBin)));
	    denominator->SetBinError(iBin, TMath::Sqrt(denominator->GetBinContent(iBin)));
	  }
	  
	  TGraphAsymmErrors* fr = new TGraphAsymmErrors(numerator, denominator);
	  //TGraphAsymmErrors* fr = getEfficiency(numerator, denominator);

	  frMap[eventSelection->first][process->first][tauId->first][observable->first] = fr;
	  
	  delete numerator;
	  delete denominator;
	}
      }
    }
  }
}

//
//-------------------------------------------------------------------------------
//

void sumHistograms(histogramMapType5& histogramMap,
		   const vstring& processesToSum, const std::string& processNameSum,
		   double mcToDataScaleFactor = 1.0)
{
//--------------------------------------------------------------------------------
// Compute total Standard Model expectation by summing all signal and background contributions
//--------------------------------------------------------------------------------

  for ( histogramMapType5::iterator eventSelection = histogramMap.begin();
	eventSelection != histogramMap.end(); ++eventSelection ) {
    for ( vstring::const_iterator process = processesToSum.begin();
	  process != processesToSum.end(); ++process ) {
      for ( histogramMapType3::iterator tauId = eventSelection->second[*process].begin();
	    tauId != eventSelection->second[*process].end(); ++tauId )
	for ( histogramMapType2::iterator observable = tauId->second.begin();
	      observable != tauId->second.end(); ++observable ) {
	  for ( histogramMapType1::iterator region = observable->second.begin();
		region != observable->second.end(); ++region ) {
	    TH1* histogram = region->second;

	    TH1* histogramSum = eventSelection->second[processNameSum][tauId->first][observable->first][region->first];
	    if ( !histogramSum ) {
	      std::string histogramName = histogram->GetName();
	      std::string histogramSumName = std::string(processNameSum).append(std::string(histogramName, histogramName.find("_")));
	      histogramSum = (TH1*)histogram->Clone(histogramSumName.data());
	      histogramSum->Scale(mcToDataScaleFactor);
	      histogramMap[eventSelection->first][processNameSum][tauId->first][observable->first][region->first] = histogramSum;
	    } else {	    	      
	      histogramSum->Add(histogram, mcToDataScaleFactor);
	    }
	  }
	}
    }
  }
}

//
//-------------------------------------------------------------------------------
//

std::vector<TPaveText*> drawLabels(
  const vstring& labels, double textSize = 0.0450, unsigned color = 1, double xOffset = 0.1500, double yOffset = 0.8075)
{
  std::vector<TPaveText*> retVal;

  size_t numLabels = labels.size();
  for ( size_t iLabel = 0; iLabel < numLabels; ++iLabel ) {
    double x0 = xOffset;
    double y0 = yOffset + (numLabels - (iLabel + 1))*(textSize + 0.0075);
    double x1 = x0 + 0.3200;
    double y1 = y0 + 0.0400;

    TPaveText* label = new TPaveText(x0, y0, x1, y1, "NDC");
    label->AddText(labels[iLabel].data());
    label->SetTextAlign(13);
    label->SetTextSize(textSize);
    label->SetTextColor(color);
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    
    label->Draw();
    
    retVal.push_back(label);
  }
  
  return retVal;
}

TPaveText* drawLabel(
  const std::string& label, double textSize = 0.045, unsigned color = 1, double xOffset = 0.150, double yOffset = 0.8075)
{
  vstring labels;
  labels.push_back(label);

  return drawLabels(labels, textSize, color, xOffset, yOffset)[0];
}

std::string getOutputFilePath(const std::string& outputFileName)
{
  size_t pos = outputFileName.find_last_of('/');
  std::string retVal;
  if ( pos != std::string::npos ) retVal = std::string(outputFileName, 0, pos + 1);
  else retVal = "./";
  retVal.append("plots/");
  return retVal;
}

std::string getOutputFileName(const std::string& outputFileName, const std::string& suffix)
{
  size_t pos1 = outputFileName.find_last_of('/');
  std::string outputFileName_base = ( pos1 != std::string::npos ) ?
    std::string(outputFileName, pos1 + 1) : outputFileName;
  
  size_t pos2 = outputFileName_base.find_last_of('.');
  if ( pos2 == std::string::npos )
    throw cms::Exception("drawGraphs") 
      << "Failed to find '.' in outputFileName_base = " << outputFileName << " passed as function argument !!\n";

  std::string retVal;
  retVal.append(std::string(outputFileName_base, 0, pos2));
  retVal.append(suffix);
  retVal.append(std::string(outputFileName_base, pos2));
  return retVal;
}

void drawHistograms(
  const std::string& observable, const std::string& xAxisTitle,
  TCanvas* canvas, histogramMapType5& histogramMap, 
  const vstring& processNamesSim, const std::string& processNameData, std::map<std::string, histogramDrawOptionType>& drawOptions,
  const vstring& eventSelectionsToPlot, const vstring& tauIdsToPlot, const vstring& regionsToPlot,
  const vstring& labels, const std::string& outputFileName)
{
  //std::cout << "<drawHistograms>:" << std::endl;
  //std::cout << " observable = " << observable << std::endl;
  //std::cout << " processNamesSim = " << format_vstring(processNamesSim) << std::endl;
  //std::cout << " processNameData = " << processNameData << std::endl;

  for ( vstring::const_iterator eventSelection = eventSelectionsToPlot.begin();
	eventSelection != eventSelectionsToPlot.end(); ++eventSelection ) {
    for ( vstring::const_iterator tauId = tauIdsToPlot.begin();
	  tauId != tauIdsToPlot.end(); ++tauId ) {
      for ( vstring::const_iterator region = regionsToPlot.begin();
	    region != regionsToPlot.end(); ++region ) {
	//std::cout << "--> drawing histograms for eventSelection = " << (*eventSelection) << ","
	//	    << " tauId = " << (*tauId) << ", region = " << (*region) << std::endl;

	canvas->SetLogy(true);
	canvas->Clear();
	canvas->SetLeftMargin(0.14);
	canvas->SetBottomMargin(0.12);
	
	TLegend legend(0.64, 0.64, 0.89, 0.89, "", "brNDC"); 
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	
	THStack smSum("smSum", "smSum");
	
	for ( vstring::const_iterator process = processNamesSim.begin();
	      process != processNamesSim.end(); ++process ) {
	  TH1* histogramSim = histogramMap[*eventSelection][*process][*tauId][observable][*region];
	  //std::cout << "process = " << (*process) << ": histogram = " << histogramSim 
	  //	      << " (integral = " << histogramSim->Integral() << ")" << std::endl;
	  const histogramDrawOptionType& drawOptionsSim = drawOptions[*process];
	  applyStyleOption_histogram(histogramSim, "", drawOptionsSim, xAxisTitle);
	  smSum.Add(histogramSim);
	  legend.AddEntry(histogramSim, drawOptionsSim.legendEntry_.data(), drawOptionsSim.drawOptionLegend_.data());
	}
	
	TH1* histogramData = histogramMap[*eventSelection][processNameData][*tauId][observable][*region];	  
	//std::cout << "Data: histogram = " << histogramData 
	//	    << " (integral = " << histogramData->Integral() << ")" << std::endl;
	const histogramDrawOptionType& drawOptionsData = drawOptions[processNameData];
	applyStyleOption_histogram(histogramData, "", drawOptionsData, xAxisTitle);
	legend.AddEntry(histogramData, drawOptionsData.legendEntry_.data(), drawOptionsData.drawOptionLegend_.data());
	
	double yMax = TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum());
	if ( canvas->GetLogy() ) smSum.SetMaximum(TMath::Exp(1.5*TMath::Log(yMax)));
	else smSum.SetMaximum(1.5*yMax);
	smSum.SetMinimum(1.);
	smSum.Draw("hist");
	applyStyleOption_histogram(&smSum, "", xAxisTitle);

	histogramData->SetStats(false);
	histogramData->Draw("ep1same");
	
	legend.Draw();
	std::vector<TPaveText*> labels_text = drawLabels(labels);
	
	canvas->Update();

	std::string outputFilePath = getOutputFilePath(outputFileName);
	gSystem->mkdir(outputFilePath.data(), true);
	std::string suffix = std::string("_").append(*eventSelection).append("_").append(*tauId).append("_").append(*region);
	suffix.append("_").append(observable);
	std::string outputFileName_full = std::string(outputFilePath).append(getOutputFileName(outputFileName, suffix));
	size_t idx = outputFileName_full.find_last_of('.');
	std::string outputFileName_plot = std::string(outputFileName_full, 0, idx);
	if ( idx != std::string::npos ) 
	  canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName_full, idx)).data());
	canvas->Print(std::string(outputFileName_plot).append(".png").data());
	canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

	for ( std::vector<TPaveText*>::iterator it = labels_text.begin();
	      it != labels_text.end(); ++it ) {
	  delete (*it);
	}	
      }
    }
  }
}

void drawGraphs(
  const std::string& observable, int numBinsX, double xMin, double xMax, const std::string& xAxisTitle,
  TCanvas* canvas, frMapType3& frMapToPlot, 
  const std::string& processNameSim, const std::string& processNameData, 
  const vstring& graphsToPlot, std::map<std::string, graphDrawOptionType>& drawOptions,
  const std::string& plotName, const vstring& labels, const std::string& addPlotLabel, const std::string& outputFileName)
{
  //std::cout << "<drawGraphs>:" << std::endl;

  canvas->SetLogy(false);
  canvas->Clear();
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  TPad* topPad = new TPad("topPad", "topPad", 0.00, 0.35, 1.00, 1.00);
  topPad->SetFillColor(10);
  topPad->SetTopMargin(0.04);
  topPad->SetLeftMargin(0.15);
  topPad->SetBottomMargin(0.03);
  topPad->SetRightMargin(0.05);
  topPad->SetLogy(true);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.24);
  bottomPad->SetRightMargin(0.05);
  
  double legendY0 = 0.95 - 2*graphsToPlot.size()*0.05;
  TLegend legend(0.53, legendY0, 0.94, 0.95, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);

  canvas->cd();
  topPad->Draw();
  topPad->cd();

  TH1* dummyHistogram_top = new TH1F("dummyHistogram_top", "dummyHistogram_top", numBinsX, xMin, xMax);

  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);
    
  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle("Fake-rate");
  yAxis_top->SetTitleOffset(0.90);
  yAxis_top->SetTitleSize(0.06);

  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  dummyHistogram_top->SetMaximum(1.e-0);
  dummyHistogram_top->SetMinimum(1.e-4);
  dummyHistogram_top->Draw("axis");

  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphSim = frMapToPlot[processNameSim][*graphName][observable];
    //std::cout << processNameSim << ", graphName = " << (*graphName) << ", observable = " << observable << ":"
    //	        << " graph = " << graphSim << std::endl;
    applyStyleOption_graph(graphSim, "", drawOptions[*graphName].expMarkerStyle_, drawOptions[*graphName].markerSize_, 
			   drawOptions[*graphName].color_, xAxisTitle);    

    TGraph* graphData = frMapToPlot[processNameData][*graphName][observable];
    //std::cout << processNameData << ", graphName = " << (*graphName) << ", observable = " << observable << ":"
    //	        << " graph = " << graphSim << std::endl;
    applyStyleOption_graph(graphData, "", drawOptions[*graphName].measMarkerStyle_, drawOptions[*graphName].markerSize_, 
			   drawOptions[*graphName].color_, xAxisTitle);

    graphData->Draw("P");
    legend.AddEntry(graphData, std::string(drawOptions[*graphName].legendEntry_).append(" Data").data(), "p");

    graphSim->Draw("P");      
    legend.AddEntry(graphSim, std::string(drawOptions[*graphName].legendEntry_).append(" Simulation").data(), "p");
  }

  legend.Draw();
  std::vector<TPaveText*> labels_text = drawLabels(labels, 0.0450, 1, 0.1650, 0.8575);
  TPaveText* addPlotLabel_text = drawLabel(addPlotLabel, 0.0450, 2, 0.1800, 0.1050);

  frMapType1 frMapDiff;

  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraphAsymmErrors* graphSim  = frMapToPlot[processNameSim][*graphName][observable];
    TGraphAsymmErrors* graphData = frMapToPlot[processNameData][*graphName][observable];
    // CV: unequal number of points may actually happen in case of empty bins
    //     in the numerator/denominator histograms (low event statistics) and is **not** an error
    //if ( !(graphSim->GetN() == graphData->GetN()) )
    //  throw cms::Exception("drawGraphs") 
    //	  << "Incompatible binning of graphs for Data = " << graphData->GetName() << " (" << graphData->GetN() << ")"
    //	  << " and Simulation = " << graphSim->GetName() << " (" << graphSim->GetN() << ")" << " !!\n";
    
    TGraphAsymmErrors* graphDiff = (TGraphAsymmErrors*)graphData->Clone();

    std::vector<bool> drawPoints(graphDiff->GetN()); // true = draw/false = remove point
    
    int numPoints = graphData->GetN();
    for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
      double xCenter, ySim, yData, dummy;
      graphSim->GetPoint(iPoint, dummy, ySim);
      graphData->GetPoint(iPoint, xCenter, yData);

      double dxLow     = graphData->GetErrorXlow(iPoint);
      double dxUp      = graphData->GetErrorXhigh(iPoint);
      double dyLowSim  = graphSim->GetErrorYlow(iPoint);
      double dyUpSim   = graphSim->GetErrorYhigh(iPoint);
      double dyLowData = graphData->GetErrorYlow(iPoint);
      double dyUpData  = graphData->GetErrorYhigh(iPoint);

      if ( ySim > 0. ) {
	graphDiff->SetPoint(iPoint, xCenter, (yData - ySim)/ySim);
	
	double dyLowDiff2 = square(dyLowData/ySim) + square((yData/ySim)*(dyUpSim/ySim));
        double dyUpDiff2 = square(dyUpData/ySim) + square((yData/ySim)*(dyLowSim/ySim));

	graphDiff->SetPointError(iPoint, dxLow, dxUp, TMath::Sqrt(dyLowDiff2), TMath::Sqrt(dyUpDiff2));
	drawPoints[iPoint] = true;
      } else { 
	drawPoints[iPoint] = false;
      }
    }
    
    for ( int iPoint = (numPoints - 1); iPoint >= 0; --iPoint ) {
      if ( !drawPoints[iPoint] ) graphDiff->RemovePoint(iPoint);
    }
    
    frMapDiff[*graphName] = graphDiff;
  }

  double maxDiff = 0.;    
  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphDiff = frMapDiff[*graphName];
    //std::cout << "graphName = " << (*graphName) << ", graphDiff = " << graphDiff << std::endl;

    int numPoints = graphDiff->GetN();
    for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {

      double x, diff;
      graphDiff->GetPoint(iPoint, x, diff);
      if ( diff > maxDiff ) maxDiff = diff;
      //double err = TMath::Max(graphDiff->GetErrorYlow(iPoint), graphDiff->GetErrorYhigh(iPoint));
      //if ( err  > maxDiff ) maxDiff = err;
    }
  }

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = (TH1*)dummyHistogram_top->Clone("dummyHistogram_bottom");

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(1.20);
  xAxis_bottom->SetTitleSize(0.08);
  xAxis_bottom->SetLabelOffset(0.02);
  xAxis_bottom->SetLabelSize(0.08);
  xAxis_bottom->SetTickLength(0.055);
  xAxis_bottom->SetLabelColor(1);
  xAxis_bottom->SetTitleColor(1);

  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("#frac{Data - Simulation}{Simulation}");
  yAxis_bottom->SetTitleOffset(0.85);
  yAxis_bottom->SetNdivisions(505);
  yAxis_bottom->CenterTitle();
  yAxis_bottom->SetTitleSize(0.08);
  yAxis_bottom->SetLabelSize(0.08);
  yAxis_bottom->SetTickLength(0.04);

  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);
  dummyHistogram_bottom->SetMaximum(+0.5);
  dummyHistogram_bottom->SetMinimum(-0.5);
  dummyHistogram_bottom->Draw("axis");

  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphDiff = frMapDiff[*graphName];
    applyStyleOption_graph(graphDiff, "", drawOptions[*graphName].measMarkerStyle_, 
			   drawOptions[*graphName].markerSize_, drawOptions[*graphName].color_, xAxisTitle);
    graphDiff->Draw("P");
  }

  canvas->Update();

  std::string outputFilePath = getOutputFilePath(outputFileName);
  gSystem->mkdir(outputFilePath.data(), true);
  std::string suffix = std::string("_fr").append(plotName).append("_").append(observable);
  std::string outputFileName_full = std::string(outputFilePath).append(getOutputFileName(outputFileName, suffix));
  size_t idx = outputFileName_full.find_last_of('.');
  std::string outputFileName_plot = std::string(outputFileName_full, 0, idx);
  if ( idx != std::string::npos ) 
    canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName_full, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

  for ( std::vector<TPaveText*>::iterator it = labels_text.begin();
	it != labels_text.end(); ++it ) {
    delete (*it);
  }

  delete addPlotLabel_text;

  for ( frMapType1::iterator it = frMapDiff.begin();
	it != frMapDiff.end(); ++it ) {
    delete it->second;
  }

  delete dummyHistogram_bottom;
  delete topPad;
  delete dummyHistogram_top;
  delete bottomPad;
}

#endif
