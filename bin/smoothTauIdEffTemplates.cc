
/** \executable smoothTauIdEffTemplates
 *
 * Fit template histograms used in tau id. efficiency measurement by analytic functions,
 * in order to avoid bias on results arising from statistical fluctuations of template shapes
 *
 * \author Betty Calpas, RWTH Aachen
 *         Christian Veelken, LLR
 *
 * \version $Revision: 1.5 $
 *
 * $Id: smoothTauIdEffTemplates.cc,v 1.5 2012/06/15 20:48:53 calpas Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TBenchmark.h>
#include "TCanvas.h"
#include "TAxis.h"
#include "TDirectory.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"

using namespace RooFit;

void smoothHistogram(TH1* histogram, const std::string& fitFunctionType, 
		     TFileDirectory& histogramOutputDirectory,
		     bool makeControlPlots, const std::string& controlPlotFilePath)
{  
  TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();

  // general fit variable
  RooRealVar x("x", "x", xMin, xMax);

  RooAbsPdf* fitFunction = 0;
  std::vector<TObject*> objectsToDelete;
  if ( fitFunctionType == "LG1" ) {
    // create Landau function
    RooRealVar* meanl = new RooRealVar("meanl", "mean of Landau", 80., 60., 90.);
    RooRealVar* sigmal = new RooRealVar("sigmal", "sigma of Landau", 30., 20., 50.);  
    RooLandau* landau = new RooLandau("landau", "landau", x, *meanl, *sigmal);
    // create Gaussian
    RooRealVar* meang = new RooRealVar("meang", "mean of Gaussian", 1.); 
    RooRealVar* sigmag = new RooRealVar("sigmag","sigma of Gaussian", 1., 0.1, 10.);
    RooGaussian* gauss = new RooGaussian("gauss", "gauss", x, *meang, *sigmag);
    // create convolution of Landau with Gaussian
    fitFunction = new RooFFTConvPdf("LandauConvGauss", "LandauConvGauss", x, *landau, *gauss);
    objectsToDelete.push_back(meanl);
    objectsToDelete.push_back(sigmal);
    objectsToDelete.push_back(landau);
    objectsToDelete.push_back(meang);
    objectsToDelete.push_back(sigmag);
    objectsToDelete.push_back(gauss);
  } else if ( fitFunctionType == "CB1" ||
	      fitFunctionType == "CB2" ||
	      fitFunctionType == "CB3" ) {
    // create Crystal-ball function
    RooRealVar* cbmean  = 0;
    RooRealVar* cbsigma = 0;
    RooRealVar* n       = 0;
    RooRealVar* alpha   = 0;
    if        ( fitFunctionType == "CB1" ) {
      cbmean  = new RooRealVar("cbmean", "cbmean", 90., 20., 180.);
      cbsigma = new RooRealVar("cbsigma", "cbsigma", 1., 1., 40.); 
      n       = new RooRealVar("n", "", 0.2); 
      alpha   = new RooRealVar("alpha", "", 1.3);
    } else if ( fitFunctionType == "CB2" ) {
      cbmean  = new RooRealVar("cbmean", "cbmean", 70., 20., 180.);
      cbsigma = new RooRealVar("cbsigma", "cbsigma", 10., 1., 200.); 
      n       = new RooRealVar("n", "", 1.); 
      alpha   = new RooRealVar("alpha", "", 1.); 
    } else if ( fitFunctionType == "CB3" ) {
      cbmean  = new RooRealVar("cbmean", "cbmean" , 1., 20., 180.);
      cbsigma = new RooRealVar("cbsigma", "cbsigma" , 10., 1., 200.); 
      n       = new RooRealVar("n","", 1.); 
      alpha   = new RooRealVar("alpha", "", 1.);      
    }     
    fitFunction = new RooCBShape("CristalBall", "CristalBall", x, *cbmean, *cbsigma, *alpha, *n);
    objectsToDelete.push_back(cbmean);
    objectsToDelete.push_back(cbsigma);
    objectsToDelete.push_back(n);
    objectsToDelete.push_back(alpha);
  }

  if ( !fitFunction ) 
    throw cms::Exception("smoothTauIdEffTemplates") 
      << "Undefined fit-function type = " << fitFunctionType << " !!\n";
  objectsToDelete.push_back(fitFunction);
  
  // convert template histogram to RooDataHist object in order for it to be fitted
  RooDataHist data("", "", x, histogram, 1.0); 
  
  RooPlot* frame = x.frame(Title("Fit"));
  
  // fit shape template by analytic function
  fitFunction->fitTo(data);
  data.plotOn(frame);
  fitFunction->plotOn(frame, LineColor(kRed)); 
  
  // create smoothed output histogram
  std::string histogramName_smoothed = std::string(histogram->GetName()).append("_smoothed");
  int histogramBinning = xAxis->GetNbins();
  TH1* histogram_smoothed = fitFunction->createHistogram(histogramName_smoothed.data(), x, Binning(histogramBinning));

  // make control plot
  if ( makeControlPlots ) {
    TCanvas* canvas = new TCanvas("Fit", "Fit", 800, 600);
    gPad->SetLeftMargin(0.15); 
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->Draw();
    
    std::string outputFileName_plot = controlPlotFilePath;
    if ( outputFileName_plot.find_last_of("/") != (outputFileName_plot.length() - 1) ) outputFileName_plot.append("/");
    outputFileName_plot.append(Form("smoothTauIdEffTemplate_%s", histogram->GetName()));    
    canvas->Print(std::string(outputFileName_plot).append(".eps").data());
    canvas->Print(std::string(outputFileName_plot).append(".png").data());
    canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  }
  
  // scale output histogram to match integral of input histogram
  histogram_smoothed->Scale(histogram->Integral()/histogram_smoothed->Integral());

  // register output histogram with TFileService
  // (histogram will get saved in output file automatically;
  //  code copied from CommonTools/Utils/interface/TFileDirectory.h, version 1.9)
  TDirectory* dir = histogramOutputDirectory.getBareDirectory();
  dir->cd();
  ROOT::DirAutoAdd_t func = TH1::Class()->GetDirectoryAutoAdd();
  if ( func ) { 
    TH1AddDirectorySentry sentry; 
    func(histogram_smoothed, dir); 
  } else { 
    dir->Append(histogram_smoothed); 
  }

  for ( std::vector<TObject*>::iterator it = objectsToDelete.begin();
	it != objectsToDelete.end(); ++it ) {
    delete (*it);
  }
}

struct histogramEntryType
{
  std::string histogramName_;
  std::string fitFunctionType_;
};

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<smoothTauIdEffTemplates>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("smoothTauIdEffTemplates");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("smoothTauIdEffTemplates") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgSmoothTauIdEffTemplates = cfg.getParameter<edm::ParameterSet>("smoothTauIdEffTemplates"); 
 
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet histogramsToSmooth = cfgSmoothTauIdEffTemplates.getParameter<vParameterSet>("histogramsToSmooth");
  std::vector<histogramEntryType> inputHistogramEntries;
  for ( vParameterSet::const_iterator histogramToSmooth = histogramsToSmooth.begin();
	histogramToSmooth != histogramsToSmooth.end(); ++histogramToSmooth ) {
    histogramEntryType inputHistogramEntry;
    inputHistogramEntry.histogramName_ = histogramToSmooth->getParameter<std::string>("histogramName");
    inputHistogramEntry.fitFunctionType_ = histogramToSmooth->getParameter<std::string>("fitFunctionType");
    inputHistogramEntries.push_back(inputHistogramEntry);
  }
 
  bool makeControlPlots = cfgSmoothTauIdEffTemplates.getParameter<bool>("makeControlPlots");
  std::string controlPlotFilePath = cfgSmoothTauIdEffTemplates.getParameter<std::string>("controlPlotFilePath");

  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("smoothTauIdEffTemplates") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgSmoothTauIdEffTemplates.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("smoothTauIdEffTemplates") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory histogramOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;
  
  for ( std::vector<histogramEntryType>::const_iterator inputHistogramEntry = inputHistogramEntries.begin();
	inputHistogramEntry != inputHistogramEntries.end(); ++inputHistogramEntry ) {
    // retrieve template histogram to be fitted from input file
    TH1* inputHistogram = dynamic_cast<TH1*>(histogramInputFile->Get(inputHistogramEntry->histogramName_.data()));
    if ( !inputHistogram ) 
      throw cms::Exception("smoothTauIdEffTemplates") 
	<< "Failed to find histogram = " << inputHistogramEntry->histogramName_ << " in input file = " << histogramFileName << " !!\n";

    // fit template histogram by analytic function;
    // store smoothed shape template in outputFile
    smoothHistogram(inputHistogram, inputHistogramEntry->fitFunctionType_, 
		    histogramOutputDirectory, 
		    makeControlPlots, controlPlotFilePath);
  }
  
  return 0;
}
