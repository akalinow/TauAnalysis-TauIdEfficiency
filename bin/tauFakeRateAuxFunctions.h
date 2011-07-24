#ifndef TauAnalysis_GenSimTools_tauFakeRateAuxFunctions_h
#define TauAnalysis_GenSimTools_tauFakeRateAuxFunctions_h

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

// define auxiliary class holding draw-options for control plots
struct histogramDrawOptionType 
{
  histogramDrawOptionType()
  {}
  histogramDrawOptionType(const edm::ParameterSet& cfg)
    : name_(cfg.getParameter<std::string>("name")),
      legendEntry_(cfg.getParameter<std::string>("legendEntry")),
      lineColor_(cfg.getParameter<int>("lineColor")),
      lineStyle_(cfg.getParameter<int>("lineStyle")),
      lineWidth_(cfg.getParameter<int>("lineWidth")),
      fillColor_(cfg.getParameter<int>("fillColor")),
      fillStyle_(cfg.getParameter<int>("fillStyle")),
      drawOption_(cfg.getParameter<std::string>("drawOption")),
      drawOptionLegend_(cfg.getParameter<std::string>("drawOptionLegend"))
  {}
  ~histogramDrawOptionType() {}
  std::string name_;
  std::string legendEntry_;
  int lineColor_;
  int lineStyle_;
  int lineWidth_;
  int fillColor_;
  int fillStyle_;
  std::string drawOption_;
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
      color_(cfg.getParameter<unsigned>("color"))
  {}
  ~graphDrawOptionType() {}
  std::string name_;
  std::string legendEntry_;
  unsigned expMarkerStyle_;
  unsigned measMarkerStyle_;
  unsigned color_;
};

void applyStyleOption_histogram(
  TH1* histogram, const std::string& histogramTitle, const histogramDrawOptionType& drawOptions, 
  const std::string& xAxisTitle, const std::string& yAxisTitle = "Events")
{
  //std::cout << "<applyStyleOption_histogram>:" << std::endl;
  //std::cout << " histogramName = " << histogram->GetName() << std::endl;

  histogram->SetTitle(histogramTitle.data());

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
  unsigned markerStyle, unsigned color,
  const std::string& xAxisTitle, const std::string& yAxisTitle = "Fake-rate")
{
  //std::cout << "<applyStyleOption_graph>:" << std::endl;
  //std::cout << " graphName = " << graph->GetName() << std::endl;

  graph->SetTitle(graphTitle.data());

  graph->SetMarkerStyle(markerStyle);
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
  const vstring& processes, const vstring& eventSelections, const vstring& tauIds,
  const vstring& observables, const vstring& regions)
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
	<< "Failed to load histograms for event selection," 
	<< " directory = " << inputDirectoryName << " does not exists in input file = " << inputFile->GetName() << " !!\n";
    
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
	    if ( !histogram ) 
	      throw cms::Exception("loadHistograms") 
		<< "Failed to load histogram = " << histogramName << " from directory = " << inputDirectory->GetName() << " !!\n";
	    
	    if ( !histogram->GetSumw2N() ) histogram->Sumw2();

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
//
//    CV: Make numerator and denominator histograms integer
//        before computing the tau id. efficiency/fake-rate,
//        in order to work-around the problem
//        that the TGraphAsymmErrors(TH1*, TH1*) constructor does not work in ROOT 5.27/06b (default in CMSSW_4_1_3),
//        in case the histograms contain weighted entries.
//
//        The binContents of numerator and denominator histograms are scaled
//        to the "effective" number of (denominator) entries:
//          effNum = (sum(weights)/sqrt(sum(weights^2)))^2 = binContent/(binError^2)

	  std::string numeratorName = std::string(process->first).append("_numerator_");
	  numeratorName.append(observable->first).append("_").append(tauId->first);
	  TH1* numerator = (TH1*)histogram_passed->Clone(numeratorName.data());
	  
	  std::string denominatorName = std::string(process->first).append("_denominator_");
	  denominatorName.append(observable->first).append("_").append(tauId->first);
	  TH1* denominator = (TH1*)histogram_failed->Clone(denominatorName.data());
	  denominator->Add(numerator);
	  
	  for ( int iBin = 1; iBin <= (numerator->GetNbinsX() + 1); ++iBin ) {
	    double scaleFactor = ( denominator->GetBinError(iBin) > 0. ) ?
	      denominator->GetBinContent(iBin)/TMath::Power(denominator->GetBinError(iBin), 2.) : 0.;
	    
	    numerator->SetBinContent(iBin, TMath::Nint(scaleFactor*numerator->GetBinContent(iBin)));
	    numerator->SetBinError(iBin, TMath::Sqrt(numerator->GetBinContent(iBin)));
	    
	    denominator->SetBinContent(iBin, TMath::Nint(scaleFactor*denominator->GetBinContent(iBin)));
	    denominator->SetBinError(iBin, TMath::Sqrt(denominator->GetBinContent(iBin)));
	  }
	  
	  TGraphAsymmErrors* fr = new TGraphAsymmErrors(numerator, denominator);
	  
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

void drawHistograms(
  const std::string& observable, const std::string& xAxisTitle,
  TCanvas* canvas, histogramMapType5& histogramMap, 
  const vstring& processNamesSim, const std::string& processNameData, std::map<std::string, histogramDrawOptionType>& drawOptions,
  const vstring& eventSelectionsToPlot, const vstring& tauIdsToPlot, const vstring& regionsToPlot,
  const vstring& labels, const std::string& outputFileName)
{
  for ( vstring::const_iterator eventSelection = eventSelectionsToPlot.begin();
	eventSelection != eventSelectionsToPlot.end(); ++eventSelection ) {
    for ( vstring::const_iterator tauId = tauIdsToPlot.begin();
	  tauId != tauIdsToPlot.end(); ++tauId ) {
      for ( vstring::const_iterator region = regionsToPlot.begin();
	    region != regionsToPlot.end(); ++region ) {

	canvas->SetLogy(true);
	canvas->Clear();
	canvas->SetLeftMargin(0.14);
	canvas->SetBottomMargin(0.12);
	
	TLegend legend(0.64, 0.59, 0.89, 0.89, "", "brNDC"); 
	legend.SetBorderSize(0);
	legend.SetFillColor(0);
	
	THStack smSum("smSum", "smSum");
	
	for ( vstring::const_iterator process = processNamesSim.begin();
	      process != processNamesSim.end(); ++process ) {
	  TH1* histogramSim = histogramMap[*eventSelection][*process][*tauId][observable][*region];
	  applyStyleOption_histogram(histogramSim, "", drawOptions[*process], xAxisTitle);
	  smSum.Add(histogramSim);
	  legend.AddEntry(histogramSim, drawOptions[*process].drawOptionLegend_.data(), "f");
	}
	
	TH1* histogramData = histogramMap[*eventSelection][processNameData][*tauId][observable][*region];	  
	applyStyleOption_histogram(histogramData, "", drawOptions[processNameData], xAxisTitle);
	legend.AddEntry(histogramData, drawOptions[processNameData].drawOptionLegend_.data(), "p");

	smSum.SetMaximum(1.5*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));
	smSum.Draw("hist");
	applyStyleOption_histogram(&smSum, "", xAxisTitle);

	histogramData->SetStats(false);
	histogramData->Draw("ep1same");
	
	legend.Draw();
	std::vector<TPaveText*> labels_text = drawLabels(labels);
	
	canvas->Update();
	std::string outputFilePath = std::string("./plots/");
	gSystem->mkdir(outputFilePath.data(), true);
	std::string suffix = std::string("_").append(*eventSelection).append("_").append(*tauId).append("_").append(*region);
	suffix.append("_").append(observable).append(".");
	canvas->Print(TString(outputFilePath.append(outputFileName).data()).ReplaceAll(".", suffix.data()));

	for ( std::vector<TPaveText*>::iterator it = labels_text.begin();
	      it != labels_text.end(); ++it ) {
	  delete (*it);
	}	
      }
    }
  }
}

double square(double x)
{
  return x*x;
}

void drawGraphs(
  const std::string& observable, const std::string& xAxisTitle,
  TCanvas* canvas, frMapType3& frMapToPlot, 
  const std::string& processNameSim, const std::string& processNameData, 
  const vstring& graphsToPlot, std::map<std::string, graphDrawOptionType>& drawOptions,
  const std::string& plotName, const vstring& labels, const std::string& addPlotLabel, const std::string& outputFileName)
{
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

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.20);
  bottomPad->SetRightMargin(0.05);
  
  TLegend legend(0.53, 0.63, 0.94, 0.95, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  canvas->cd();
  topPad->Draw();
  topPad->cd();
  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphSim = frMapToPlot[*graphName][processNameSim][observable];
    applyStyleOption_graph(graphSim, "", drawOptions[*graphName].expMarkerStyle_, drawOptions[*graphName].color_, xAxisTitle);

    TGraph* graphData = frMapToPlot[*graphName][processNameData][observable];
    applyStyleOption_graph(graphData, "", drawOptions[*graphName].measMarkerStyle_, drawOptions[*graphName].color_, xAxisTitle);

    TAxis* xAxis = graphData->GetXaxis();
    xAxis->SetLabelColor(10);
    xAxis->SetTitleColor(10);
    
    TAxis* yAxis = graphData->GetYaxis();
    yAxis->SetTitle("Fake-rate");
    yAxis->SetTitleOffset(1.15);
    yAxis->SetTitleSize(0.06);
    
    graphData->SetTitle("");
    graphData->SetMaximum(1.e0);
    graphData->SetMinimum(1.e-4);
    std::string drawOption_string = ( graphName == graphsToPlot.begin() ) ? "A" : "";
    graphData->Draw(drawOption_string.append("P").data());      
    legend.AddEntry(graphData, std::string(drawOptions[*graphName].legendEntry_).append(" Data").data(), "p");
    
    graphSim->Draw("P");      
    legend.AddEntry(graphSim, std::string(drawOptions[*graphName].legendEntry_).append(" Simulation").data(), "p");
  }

  legend.Draw();
  std::vector<TPaveText*> labels_text = drawLabels(labels);
  TPaveText* addPlotLabel_text = drawLabel(addPlotLabel, 0.035, 2, 0.5000, 0.6150);

  frMapType1 frMapDiff;

  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraphAsymmErrors* graphSim  = frMapToPlot[*graphName][processNameSim][observable];
    TGraphAsymmErrors* graphData = frMapToPlot[*graphName][processNameData][observable];
    
    if ( !(graphSim->GetN() == graphData->GetN()) )
      throw cms::Exception("drawGraphs") 
	<< "Incompatible binning of graphs for Data = " << graphData->GetName() 
	<< " and Simulation = " << graphSim->GetName() << " !!\n";
    
    TGraphAsymmErrors* graphDiff = (TGraphAsymmErrors*)graphData->Clone();
    
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
      }
    }

    frMapDiff[*graphName] = graphDiff;
  }

  double maxDiff = 0.;    
  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphDiff = frMapDiff[*graphName];

    int numPoints = graphDiff->GetN();
    for ( int iPoint = 0; iPoint < numPoints; ++iPoint ) {
      double x, diff;
      graphDiff->GetPoint(iPoint, x, diff);
      if ( diff > maxDiff ) maxDiff = diff;
      double err = TMath::Max(graphDiff->GetErrorYlow(iPoint), graphDiff->GetErrorYhigh(iPoint));
      if ( err  > maxDiff ) maxDiff = err;
    }
  }

  canvas->cd();
  bottomPad->Draw();
  bottomPad->cd();

  for ( vstring::const_iterator graphName = graphsToPlot.begin();
	graphName != graphsToPlot.end(); ++graphName ) {
    TGraph* graphDiff = frMapDiff[*graphName];
    applyStyleOption_graph(graphDiff, "", drawOptions[*graphName].measMarkerStyle_, drawOptions[*graphName].color_, xAxisTitle);

    TAxis* xAxis = graphDiff->GetXaxis();
    xAxis->SetTitle(xAxisTitle.data());
    xAxis->SetTitleOffset(1.20);
    xAxis->SetNdivisions(505);
    xAxis->SetTitleOffset(1.1);
    xAxis->SetTitleSize(0.08);
    xAxis->SetLabelOffset(0.02);
    xAxis->SetLabelSize(0.08);
    xAxis->SetTickLength(0.055);

    TAxis* yAxis = graphDiff->GetYaxis();
    yAxis->SetTitle("#frac{Data - Simulation}{Simulation}");
    yAxis->SetTitleOffset(1.15);
    yAxis->CenterTitle();
    yAxis->SetTitleOffset(0.9);
    yAxis->SetTitleSize(0.08);
    yAxis->SetLabelSize(0.08);
    yAxis->SetTickLength(0.04);

    graphDiff->SetTitle("");
    double maxDiff01 = 0.1*TMath::Ceil(1.2*maxDiff*10.);
    graphDiff->SetMaximum(+maxDiff01);
    graphDiff->SetMinimum(-maxDiff01);

    std::string drawOption_string = ( graphName == graphsToPlot.begin() ) ? "A" : "";
    graphDiff->Draw(drawOption_string.append("P").data());   
  }

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  gSystem->mkdir(outputFilePath.data(), true);
  std::string suffix = std::string("_fr").append(plotName).append("_").append(observable);
  canvas->Print(TString(outputFilePath.append(outputFileName).data()).ReplaceAll(".", suffix.data()));
  
  for ( std::vector<TPaveText*>::iterator it = labels_text.begin();
	it != labels_text.end(); ++it ) {
    delete (*it);
  }

  delete addPlotLabel_text;

  for ( frMapType1::iterator it = frMapDiff.begin();
	it != frMapDiff.end(); ++it ) {
    delete it->second;
  }

  delete topPad;
  delete bottomPad;
}

#endif
