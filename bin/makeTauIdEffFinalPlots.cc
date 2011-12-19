
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FWLite/interface/InputSource.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"

#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

struct tauIdEntryType
{
  tauIdEntryType(const std::string& name, const std::string& legendEntry, 
		 unsigned expMarkerStyle, unsigned measMarkerStyle, unsigned color)
    : name_(name),
      legendEntry_(legendEntry),
      expMarkerStyle_(expMarkerStyle),
      measMarkerStyle_(measMarkerStyle),
      color_(color)
  {}
  ~tauIdEntryType() {}
  std::string name_;
  std::string legendEntry_;
  unsigned expMarkerStyle_;
  unsigned measMarkerStyle_;
  unsigned color_;
};

template <typename T>
void clearMap(std::map<std::string, T*>& m)
{
  for ( typename std::map<std::string, T*>::iterator it = m.begin();
	it != m.end(); ++it ) {
    delete it->second;
  }
}

struct binEntryType
{
  binEntryType(const std::string& directory, double xBinCenter)
    : directory_(directory),
      xBinCenter_(xBinCenter)
  {}
  ~binEntryType() {}
  std::string directory_;
  double xBinCenter_;
};

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<makeTauIdEffFinalPlots>:" << std::endl;  

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeTauIdEffFinalPlots");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("makeTauIdEffFinalPlots") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgMakeTauIdEffPlots = cfg.getParameter<edm::ParameterSet>("makeTauIdEffFinalPlots");

  std::vector<tauIdEntryType> tauIds;

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIds = cfgMakeTauIdEffPlots.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauId = cfgTauIds.begin();
	cfgTauId != cfgTauIds.end(); ++cfgTauId ) {
    std::string name = cfgTauId->getParameter<std::string>("name");
    std::string legendEntry = cfgTauId->getParameter<std::string>("legendEntry");
    unsigned expMarkerStyle = cfgTauId->getParameter<unsigned>("markerStyleSim");
    unsigned measMarkerStyle = cfgTauId->getParameter<unsigned>("markerStyleData");
    unsigned color = cfgTauId->getParameter<unsigned>("color");
    tauIds.push_back(tauIdEntryType(name, legendEntry, expMarkerStyle, measMarkerStyle, color));
  }

  std::string expEff_label  = cfgMakeTauIdEffPlots.getParameter<std::string>("expEff_label");
  std::string measEff_label = cfgMakeTauIdEffPlots.getParameter<std::string>("measEff_label");

  double intLumiData = cfgMakeTauIdEffPlots.getParameter<double>("intLumiData"); // in units of fb^-1

  typedef std::vector<std::string> vstring;
  vstring fitVariables = cfgMakeTauIdEffPlots.getParameter<vstring>("fitVariables");

  typedef std::vector<double> vdouble;  
  vdouble xAxisBinning = cfgMakeTauIdEffPlots.getParameter<vdouble>("xAxisBinning");
  std::string xAxisTitle = cfgMakeTauIdEffPlots.getParameter<std::string>("xAxisTitle");

  int numBins = xAxisBinning.size() - 1;
  if ( !(numBins >= 1) ) 
    throw cms::Exception("makeTauIdEffFinalPlots") 
      << "Invalid configuration parameter 'xAxisBinning' = " << format_vdouble(xAxisBinning) << " !!\n";

  double xAxisMin = xAxisBinning[0];
  double xAxisMax = xAxisBinning[numBins];

  std::vector<binEntryType> tauIdValues;

  int binIdx = 0;

  vParameterSet cfgValues = cfgMakeTauIdEffPlots.getParameter<vParameterSet>("values");
  for ( vParameterSet::const_iterator cfgValue = cfgValues.begin();
	cfgValue != cfgValues.end(); ++cfgValue, ++binIdx ) {
    std::string directory = cfgValue->getParameter<std::string>("directory");
    double xBinCenter = cfgValue->getParameter<double>("xBinCenter");

//--- check that binCenter is within the range lowerBinEdge..upperBinEdge of 'xAxisBinning'
    if ( !(binIdx <= numBins && xAxisBinning[binIdx] < xBinCenter && xBinCenter < xAxisBinning[binIdx + 1]) )
      throw cms::Exception("makeTauIdEffFinalPlots") 
	<< "Configuration parameter 'xBinCenter' = " << xBinCenter << " is not compatible with 'xAxisBinning',"
	<< " expected range = " << xAxisBinning[binIdx] << ".." << xAxisBinning[binIdx + 1] << " !!\n";

    tauIdValues.push_back(binEntryType(directory, xBinCenter));
  }

  std::string outputFileName = cfgMakeTauIdEffPlots.getParameter<std::string>("outputFileName");

  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauIdEffFinalPlots") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string inputFileName = (*inputFiles.files().begin());
  
//--- open input file
  TFile* inputFile = TFile::Open(inputFileName.data());
  if ( !inputFile ) 
    throw cms::Exception("makeTauIdEffFinalPlots") 
      << "Failed to open inputFile = " << inputFileName << " !!\n";

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
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

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", 0.00, 0.00, 1.00, 0.35);
  bottomPad->SetFillColor(10);
  bottomPad->SetTopMargin(0.02);
  bottomPad->SetLeftMargin(0.15);
  bottomPad->SetBottomMargin(0.20);
  bottomPad->SetRightMargin(0.05);
  
  for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	fitVariable != fitVariables.end(); ++fitVariable ) {

    std::map<std::string, TGraphErrors*> expEffGraphs;  // key = tauId
    std::map<std::string, TGraphErrors*> measEffGraphs; // key = tauId
    std::map<std::string, TGraphErrors*> diffEffGraphs; // key = tauId

    for ( std::vector<tauIdEntryType>::iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {

      expEffGraphs[tauId->name_]  = new TGraphErrors(numBins);
      measEffGraphs[tauId->name_] = new TGraphErrors(numBins);
      diffEffGraphs[tauId->name_] = new TGraphErrors(numBins);

      int binIdx = 0;

      for ( std::vector<binEntryType>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue, ++binIdx ) {	
	double binLowerEdge = xAxisBinning[binIdx];
	double binUpperEdge = xAxisBinning[binIdx + 1];
	double binCenter  = 0.5*(binLowerEdge + binUpperEdge);

	TDirectory* inputDirectory = ( tauIdValue->directory_ != "" ) ?
	  dynamic_cast<TDirectory*>(inputFile->Get(tauIdValue->directory_.data())) : inputFile;
	if ( !inputDirectory ) 
	  throw cms::Exception("makeTauIdEffFinalPlots") 
	    << "Directory = " << tauIdValue->directory_ << " does not exists in input file = " << inputFileName << " !!\n";

	std::string expEffName = std::string(expEff_label).append("_").append(*fitVariable).append("_").append(tauId->name_);
	double expEff = getNumber(inputDirectory, expEffName.data(), 0).first;
        expEffGraphs[tauId->name_]->SetPoint(binIdx, binCenter, expEff);
	expEffGraphs[tauId->name_]->SetPointError(binIdx, 0.5*(binUpperEdge - binLowerEdge), 0.01);
	
	std::string measEffName = std::string(measEff_label).append("_").append(*fitVariable).append("_").append(tauId->name_);
	double measEff = getNumber(inputDirectory, measEffName.data(), 0).first;
	double measEffErr = getNumber(inputDirectory, measEffName.data(), 0).second;
	measEffGraphs[tauId->name_]->SetPoint(binIdx, binCenter, measEff);
	measEffGraphs[tauId->name_]->SetPointError(binIdx, 0.5*(binUpperEdge - binLowerEdge), measEffErr);

	if ( expEff > 0. ) {
	  diffEffGraphs[tauId->name_]->SetPoint(binIdx, binCenter, (measEff - expEff)/expEff);
	  diffEffGraphs[tauId->name_]->SetPointError(binIdx, 0.5*(binUpperEdge - binLowerEdge), measEffErr/expEff);
	}
      }
    }

    double maxDiff = 0.;    
    for ( int iBin = 0; iBin < numBins; ++iBin ) {
      for ( std::vector<tauIdEntryType>::const_iterator tauId = tauIds.begin();
	    tauId != tauIds.end(); ++tauId ) {
	double x, diff;
	diffEffGraphs[tauId->name_]->GetPoint(iBin, x, diff);
	double err = diffEffGraphs[tauId->name_]->GetErrorY(iBin);
	diff = TMath::Max(TMath::Abs(diff + err), TMath::Abs(diff - err));
	if ( diff > maxDiff ) maxDiff = diff;
      }
    }

    canvas->Clear();

    TLegend* legend = new TLegend(0.53, 0.63, 0.94, 0.95, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
  
    canvas->cd();
    topPad->Draw();
    topPad->cd();
    for ( std::vector<tauIdEntryType>::iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      TAxis* xAxis = measEffGraphs[tauId->name_]->GetXaxis();
      xAxis->SetLimits(xAxisMin, xAxisMax);
      xAxis->SetLabelColor(10);
      xAxis->SetTitleColor(10);

      TAxis* yAxis = measEffGraphs[tauId->name_]->GetYaxis();
      yAxis->SetTitle("Efficiency");
      yAxis->SetTitleOffset(1.15);
      yAxis->SetTitleSize(0.06);

      measEffGraphs[tauId->name_]->SetTitle("");
      measEffGraphs[tauId->name_]->SetMaximum(1.4);
      measEffGraphs[tauId->name_]->SetMinimum(0.0);

      measEffGraphs[tauId->name_]->SetMarkerStyle(tauId->measMarkerStyle_);
      measEffGraphs[tauId->name_]->SetMarkerColor(tauId->color_);
      measEffGraphs[tauId->name_]->SetLineColor(tauId->color_);
      std::string drawOption = ( tauId == tauIds.begin() ) ? "A" : "";
      measEffGraphs[tauId->name_]->Draw(drawOption.append("P").data());      
      legend->AddEntry(measEffGraphs[tauId->name_], std::string(tauId->legendEntry_).append(" Data").data(), "p");

      expEffGraphs[tauId->name_]->SetMarkerStyle(tauId->expMarkerStyle_);
      expEffGraphs[tauId->name_]->SetMarkerColor(tauId->color_);
      expEffGraphs[tauId->name_]->SetLineColor(tauId->color_);
      expEffGraphs[tauId->name_]->Draw("P");      
      legend->AddEntry(expEffGraphs[tauId->name_], std::string(tauId->legendEntry_).append(" Simulation").data(), "p");
    }
    legend->Draw();

    drawCMSprelimaryLabels(intLumiData, 0.175, 0.8375);

    canvas->cd();
    bottomPad->Draw();
    bottomPad->cd();

    for ( std::vector<tauIdEntryType>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      TAxis* xAxis = diffEffGraphs[tauId->name_]->GetXaxis();
      xAxis->SetLimits(xAxisMin, xAxisMax);
      xAxis->SetTitle(xAxisTitle.data());
      xAxis->SetTitleOffset(1.20);
      xAxis->SetNdivisions(505);
      xAxis->SetTitleOffset(1.1);
      xAxis->SetTitleSize(0.08);
      xAxis->SetLabelOffset(0.02);
      xAxis->SetLabelSize(0.08);
      xAxis->SetTickLength(0.055);

      TAxis* yAxis = diffEffGraphs[tauId->name_]->GetYaxis();
      yAxis->SetTitle("#frac{Data - Simulation}{Simulation}");
      yAxis->SetTitleOffset(1.15);
      yAxis->CenterTitle();
      yAxis->SetTitleOffset(0.9);
      yAxis->SetTitleSize(0.08);
      yAxis->SetLabelSize(0.08);
      yAxis->SetTickLength(0.04);

      diffEffGraphs[tauId->name_]->SetTitle("");
      double maxDiff01 = 0.1*TMath::Ceil(1.2*maxDiff*10.);
      diffEffGraphs[tauId->name_]->SetMaximum(+maxDiff01);
      diffEffGraphs[tauId->name_]->SetMinimum(-maxDiff01);

      diffEffGraphs[tauId->name_]->SetMarkerStyle(tauId->measMarkerStyle_);
      diffEffGraphs[tauId->name_]->SetMarkerColor(tauId->color_);
      diffEffGraphs[tauId->name_]->SetLineColor(tauId->color_);
      std::string drawOption = ( tauId == tauIds.begin() ) ? "A" : "";
      diffEffGraphs[tauId->name_]->Draw(drawOption.append("P").data());   
    }

    canvas->Update();
    size_t idx = outputFileName.find_last_of('.');
    std::string outputFileName_plot = std::string(outputFileName, 0, idx);
    outputFileName_plot.append("_").append(*fitVariable);
    if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
    canvas->Print(std::string(outputFileName_plot).append(".png").data());
    canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

    delete legend;

    clearMap(expEffGraphs);
    clearMap(measEffGraphs);
    clearMap(diffEffGraphs);
  }

  delete topPad;
  delete bottomPad;
  delete canvas;

  delete inputFile;

//--print time that it took macro to run
  std::cout << "finished executing makeTauIdEffFinalPlots macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("makeTauIdEffFinalPlots");

  return 0;
}
