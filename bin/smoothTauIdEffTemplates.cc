
/** \executable smoothTauIdEffTemplates
 *
 * Fit template histograms used in tau id. efficiency measurement by analytic functions,
 * in order to avoid bias on results arising from statistical fluctuations of template shapes
 *
 * \author Betty Calpas, RWTH Aachen
 *         Christian Veelken, LLR
 *
 * \version $Revision: 1.1 $
 *
 * $Id: smoothTauIdEffTemplates.cc,v 1.1 2012/06/16 17:21:40 veelken Exp $
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
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooSimWSTool.h"

#include "RooFitResult.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include <TVectorD.h>
#include "TIterator.h"
#include <math.h>

using namespace RooFit;
using namespace std;

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


////////////////////////////////////////////////////////////////////////

 cout << endl<< "XXXXXXXXXXXXXXXXXXXX Syst study XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl << endl;

  //Get FitFuncton result (mean, sigma...)
  RooFitResult* r = fitFunction->fitTo(data,Save()) ; 
  //r->Print();

  //Get the full covariance matrix
  const TMatrixDSym& cov = r->covarianceMatrix() ;
  //cov.Print() ;
  const TMatrixDSymEigen& eigencov = cov ;

  //Find the Eigenvectors and Eigen value of the covariance matrix
  const TVectorD eigenval = eigencov.GetEigenValues();
  const TMatrixD eigenvec = eigencov.GetEigenVectors();
  //cout << endl << "Eigen value associate to each eigenvector"<< endl;
  //eigenval.Print();
  //cout << endl << "Matrix of eigenvectors"<< endl;
  //eigenvec.Print();
  
   //Store function's parameter 
   RooArgSet* params  = fitFunction->getParameters(x);
   //Create iterator over the parameters
   TIterator* it = params->createIterator() ;

   //loop over eigenvector
   //http://root.cern.ch/root/html400/TVectorD.html#TVectorD:kSizeMax
   //cout << eigenval.GetNrows()  << endl; 
   //cout << eigenvec.GetNrows() << endl;
   for( int i = 0; i < eigenvec.GetNrows(); ++i ){
    //Get the Eigenvector[i]
    const TVectorD Eigenvec = eigenvec[i]; 
    //Get the Eigenvector[i] norm
    double Norm = sqrt(Eigenvec * Eigenvec); 
    //Get the Eigenvalue(i) associate to the eigenvector[i]
    double Eigenval = eigenval[i];
    //Get the error on Eigenval
    double Sigma = sqrt(Eigenval);
    //Get unit vector: projection of Eigenvector[i] over it's direction (i), divided by it's norm
    double Unit = Eigenvec(i)/Norm;

    //loop over fitFunction parameter
    //RooRealVar* var = NULL;
    RooRealVar* par;
    while(par = (RooRealVar*) it->Next()){

     //save parameter value because of parameter pointeur dependencies
     double par_value = par->getVal();
     //cout << endl << "par_value: " << par_value << endl; 

     //varie the fit parameter in direction of 1 Eigen Vector at a time 
     par->setVal(par_value + Sigma*Unit );
     //cout << endl << "par + err: " << par->getVal() << endl; 
     //cout << endl << "+ err: " <<  Eigenval * sqrt(Eigenval) / Norm  << endl; 

     //save Up syst. histogram
     std::string parName = par->GetName();
     std::string histogramName_smoothed_par_up = std::string(histogram->GetName()).append("_smoothed_").append(parName).append("_up");
     TH1* histogram_smoothed_par_up = fitFunction->createHistogram(histogramName_smoothed_par_up.data(), x, Binning(histogramBinning));

     //reset par value for down syst.
     par->setVal(par_value);

     //save Down syst. histogram
     par->setVal(par_value - Sigma*Unit );
     //cout << endl << "par - err: " << par->getVal() << endl; 
     //cout << endl << "- err: " <<  Eigenval * sqrt(Eigenval) / Norm << endl; 

     std::string histogramName_smoothed_par_down = std::string(histogram->GetName()).append("_smoothed_").append(parName).append("_down");
     TH1* histogram_smoothed_par_down = fitFunction->createHistogram(histogramName_smoothed_par_down.data(), x, Binning(histogramBinning));

     //retrive par value 
     par->setVal(par_value);

     }   
   }


  cout << endl <<  "XXXXXXXXXXXXXXXXXXXX End syst study XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxXXXXXXX"<<endl;


////////////////////////////////////////////////////////////////////////



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


/*

///////////////////////////////////////////////

  // clone fit up/down for each parameters for systematics
  RooAbsPdf* fitFunction_meanl_up     = (RooAbsPdf*) fitFunction->Clone("fitFunction_meanl_up");
  RooAbsPdf* fitFunction_meanl_down   = (RooAbsPdf*) fitFunction->Clone("fitFunction_meanl_down");
  RooAbsPdf* fitFunction_meang_up     = (RooAbsPdf*) fitFunction->Clone("fitFunction_meang_up");
  RooAbsPdf* fitFunction_meang_down   = (RooAbsPdf*) fitFunction->Clone("fitFunction_meang_down");
  RooAbsPdf* fitFunction_sigmal_up    = (RooAbsPdf*) fitFunction->Clone("fitFunction_sigmal_up");
  RooAbsPdf* fitFunction_sigmal_down  = (RooAbsPdf*) fitFunction->Clone("fitFunction_sigmal_down");
  RooAbsPdf* fitFunction_sigmag_up    = (RooAbsPdf*) fitFunction->Clone("fitFunction_sigmag_up");
  RooAbsPdf* fitFunction_sigmag_down  = (RooAbsPdf*) fitFunction->Clone("fitFunction_sigmag_down");
  //
  RooAbsPdf* fitFunction_cbmean_up    = (RooAbsPdf*) fitFunction->Clone("fitFunction_cbmean_up");
  RooAbsPdf* fitFunction_cbmean_down  = (RooAbsPdf*) fitFunction->Clone("fitFunction_cbmean_down");
  RooAbsPdf* fitFunction_cbsigma_up   = (RooAbsPdf*) fitFunction->Clone("fitFunction_cbsigma_up");
  RooAbsPdf* fitFunction_cbsigma_down = (RooAbsPdf*) fitFunction->Clone("fitFunction_cbsigma_down");
  RooAbsPdf* fitFunction_n_up         = (RooAbsPdf*) fitFunction->Clone("fitFunction_n_up");
  RooAbsPdf* fitFunction_n_down       = (RooAbsPdf*) fitFunction->Clone("fitFunction_n_down");
  RooAbsPdf* fitFunction_alpha_up     = (RooAbsPdf*) fitFunction->Clone("fitFunction_alpha_up");
  RooAbsPdf* fitFunction_alpha_down   = (RooAbsPdf*) fitFunction->Clone("fitFunction_alpha_down");


  // set new up/down parameters for systematics
  // Get all parameters of your model depending on the observable
  if ( fitFunctionType == "LG1" ) {


   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt2" << endl;

   RooArgSet* params = fitFunction->getParameters(x);

   RooRealVar* meanl = (RooRealVar*) params->find("meanl");
   RooRealVar* meang = (RooRealVar*) params->find("meang");
   RooRealVar* sigmal = (RooRealVar*) params->find("sigmal");
   RooRealVar* sigmag = (RooRealVar*) params->find("sigmag");

   double meanl_value = meanl->getVal();
   double meang_value = meang->getVal();
   double sigmal_value = sigmal->getVal();
   double sigmag_value = sigmag->getVal();
   //meanl->Print();  
   //cout << "meanl_value: " << meanl_value << endl;


   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt3" << endl;

   //inf. value??
    params = fitFunction_meang_up->getParameters(x);
   RooRealVar* meangUp = (RooRealVar*) params->find("meang");
   meangUp->setVal(meangUp->getVal() + meangUp->getError());   
   meangUp->setVal(meang_value);
   //
    params = fitFunction_meang_down->getParameters(x);
   RooRealVar* meangDown = (RooRealVar*) params->find("meang");
   meangDown->setVal(meangDown->getVal() - meangDown->getError());
   meangDown->setVal(meang_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt3" << endl;

    params = fitFunction_meanl_up->getParameters(x);
   RooRealVar* meanlUp = (RooRealVar*) params->find("meanl");
   meanlUp->Print();  
   meanlUp->setVal(meanlUp->getVal() + meanlUp->getError());
   meanlUp->Print();
   meanlUp->setVal(meanl_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt3" << endl;

    params = fitFunction_meanl_down->getParameters(x);
   RooRealVar* meanlDown = (RooRealVar*) params->find("meanl");
   meanlDown->Print();  
   meanlDown->setVal(meanlDown->getVal() - meanlDown->getError());
   meanlDown->Print();
   meanlDown->setVal(meanl_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt2" << endl;

    params = fitFunction_sigmag_up->getParameters(x);
   RooRealVar* sigmagUp = (RooRealVar*) params->find("sigmag");
   sigmagUp->Print();
   sigmagUp->setVal(sigmagUp->getVal() + sigmagUp->getError());
   sigmagUp->Print();
   sigmagUp->setVal(sigmag_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt4" << endl;

    params = fitFunction_sigmag_down->getParameters(x);
   RooRealVar* sigmagDown = (RooRealVar*) params->find("sigmag");
   sigmagDown->Print();
   sigmagDown->setVal(sigmagDown->getVal() - sigmagDown->getError());
   sigmagDown->Print();
   sigmagDown->setVal(sigmag_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt5" << endl;

    params = fitFunction_sigmal_up->getParameters(x);
   RooRealVar* sigmalUp = (RooRealVar*) params->find("sigmal");
   sigmalUp->Print();
   sigmalUp->setVal(sigmalUp->getVal() + sigmalUp->getError());
   sigmalUp->Print();
   sigmalUp->setVal(sigmal_value);

   cout<<"ooooooooooooooooooooooooooooooooooooottttttttttttttttttttttttttttttt5" << endl;

    params = fitFunction_sigmal_down->getParameters(x);
   RooRealVar* sigmalDown = (RooRealVar*) params->find("sigmal");
   sigmalDown->Print();
   sigmalDown->setVal(sigmalDown->getVal() - sigmalDown->getError());
   sigmalDown->Print();
   sigmalDown->setVal(sigmal_value);

   }

   else if (  fitFunctionType == "CB1" ||
	      fitFunctionType == "CB2" ||
	      fitFunctionType == "CB3" ) {
 

   RooArgSet* params = fitFunction->getParameters(x);

   RooRealVar* cbmean = (RooRealVar*) params->find("cbmean");
   RooRealVar* cbsigma = (RooRealVar*) params->find("cbsigma");
   RooRealVar* n = (RooRealVar*) params->find("n");
   RooRealVar* alpha = (RooRealVar*) params->find("alpha");

   double cbmean_value = cbmean->getVal();
   double cbsigma_value = cbsigma->getVal();
   double n_value = n->getVal();
   double alpha_value = alpha->getVal();
 
   params = fitFunction_cbmean_up->getParameters(x);
   RooRealVar* cbmeanUp = (RooRealVar*) params->find("cbmean");
   cbmeanUp->setVal(cbmeanUp->getVal() + cbmeanUp->getError());
   cbmeanUp->setVal(cbmean_value);


   params = fitFunction_cbmean_down->getParameters(x);
   RooRealVar* cbmeanDown = (RooRealVar*) params->find("cbmean");
   cbmeanDown->setVal(cbmeanDown->getVal() - cbmeanDown->getError());
   cbmeanDown->setVal(cbmean_value);


   params = fitFunction_cbsigma_up->getParameters(x);
   RooRealVar* cbsigmaUp = (RooRealVar*) params->find("cbsigma");
   cbsigmaUp->setVal(cbsigmaUp->getVal() - cbsigmaUp->getError());
   cbsigmaUp->setVal(cbsigma_value);


   params = fitFunction_cbsigma_down->getParameters(x);
   RooRealVar* cbsigmaDown = (RooRealVar*) params->find("cbsigma");
   cbsigmaDown->setVal(cbsigmaDown->getVal() - cbsigmaDown->getError());
   cbsigmaDown->setVal(cbsigma_value);


   params = fitFunction_n_up->getParameters(x);
   RooRealVar* nUp = (RooRealVar*) params->find("n");
   nUp->setVal(nUp->getVal() + nUp->getError());
   nUp->setVal(n_value);

   params = fitFunction_n_down->getParameters(x);
   RooRealVar* nDown = (RooRealVar*) params->find("n");
   nDown->setVal(nDown->getVal() - nDown->getError());
   nDown->setVal(n_value);

   params = fitFunction_alpha_up->getParameters(x);
   RooRealVar* alphaUp = (RooRealVar*) params->find("n");
   alphaUp->setVal(alphaUp->getVal() + alphaUp->getError());
   alphaUp->setVal(alpha_value);  
 
   params = fitFunction_alpha_down->getParameters(x);
   RooRealVar* alphaDown = (RooRealVar*) params->find("n");
   alphaDown->setVal(alphaDown->getVal() - alphaDown->getError());
   alphaDown->setVal(alpha_value);  

  }


  // create up/down histogram  
  std::string histogramName_smoothed_meanl_up   = std::string(histogram->GetName()).append("_smoothed_meanl_up");
  std::string histogramName_smoothed_meanl_down = std::string(histogram->GetName()).append("_smoothed_meanl_down");
  std::string histogramName_smoothed_meang_up   = std::string(histogram->GetName()).append("_smoothed_meang_up");
  std::string histogramName_smoothed_meang_down = std::string(histogram->GetName()).append("_smoothed_meang_down");
  //
  std::string histogramName_smoothed_sigmal_up   = std::string(histogram->GetName()).append("_smoothed_sigmal_up");
  std::string histogramName_smoothed_sigmal_down = std::string(histogram->GetName()).append("_smoothed_sigmal_down");
  std::string histogramName_smoothed_sigmag_up   = std::string(histogram->GetName()).append("_smoothed_sigmag_up");
  std::string histogramName_smoothed_sigmag_down = std::string(histogram->GetName()).append("_smoothed_sigmag_down");
  //
  std::string histogramName_smoothed_n_up       = std::string(histogram->GetName()).append("_smoothed_n_up");
  std::string histogramName_smoothed_n_down     = std::string(histogram->GetName()).append("_smoothed_n_down");
  std::string histogramName_smoothed_alpha_up   = std::string(histogram->GetName()).append("_smoothed_alpha_up");
  std::string histogramName_smoothed_alpha_down = std::string(histogram->GetName()).append("_smoothed_alpha_down");
  //
  TH1* histogram_smoothed_meanl_up = fitFunction_meanl_up->createHistogram(histogramName_smoothed_meanl_up.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_meanl_down = fitFunction_meanl_down->createHistogram(histogramName_smoothed_meanl_down.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_meang_up = fitFunction_meang_up->createHistogram(histogramName_smoothed_meang_up.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_meang_down = fitFunction_meang_down->createHistogram(histogramName_smoothed_meang_down.data(), x, Binning(histogramBinning)); 
  //
  TH1* histogram_smoothed_sigmal_up = fitFunction_sigmal_up->createHistogram(histogramName_smoothed_sigmal_up.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_sigmal_down = fitFunction_sigmal_down->createHistogram(histogramName_smoothed_sigmal_down.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_sigmag_up = fitFunction_sigmag_up->createHistogram(histogramName_smoothed_sigmag_up.data(), x, Binning(histogramBinning)); 
  TH1* histogram_smoothed_sigmag_down = fitFunction_sigmag_down->createHistogram(histogramName_smoothed_sigmag_down.data(), x, Binning(histogramBinning)); 
*/

/*
  // Loop over all bins 
  for( int ibin = 0; ibin < histogramBinning; ++ibin){

   double CentralValue = histogram_smoothed   ->GetBinContent(  xAxis->FindBin(ibin));

   //cout <<"CentralValue"<<CentralValue << endl;

  //3.1 you compute for each parameter of the analytic function the
    //           intUp_i : = integral(functionUp_i) from binEdgeLow to binEdgeHigh
    //           intDown_i = integral(functionDown_i) from binEdgeLow to binEdgeHigh
    //        functionUp = analytic function, but i-th fit parameter is set to best fit value + 1*uncertainty (uncertainty returned by RooFit/MINUIT)
    //        functionDown = analytic function, but i-th fit parameter is set to best fit value - 1*uncertainty (uncertainty returned by RooFit/MINUIT) 
   double intUp_meanl    = histogram_smoothed_meanl_up    ->GetBinContent(  xAxis->FindBin(ibin));
   double intDown_meanl  = histogram_smoothed_meanl_down  ->GetBinContent(  xAxis->FindBin(ibin));
   double intUp_meang    = histogram_smoothed_meang_up    ->GetBinContent(  xAxis->FindBin(ibin));
   double intDown_meang  = histogram_smoothed_meang_down  ->GetBinContent(  xAxis->FindBin(ibin));
   //
   double intUp_sigmal   = histogram_smoothed_sigmal_up   ->GetBinContent(  xAxis->FindBin(ibin));
   double intDown_sigmal = histogram_smoothed_sigmal_down ->GetBinContent(  xAxis->FindBin(ibin));
   double intUp_sigmag   = histogram_smoothed_sigmag_up   ->GetBinContent(  xAxis->FindBin(ibin));
   double intDown_sigmag = histogram_smoothed_sigmag_down ->GetBinContent(  xAxis->FindBin(ibin));

   //3.2 Then take maximum and miminum
                //intMax_i = max(intUp_i, intDown_i)
                //intMin_i = min(intUp_i, intDown_i) 
   double intMax_meanl = max(intUp_meanl, intDown_meanl);
   double intMin_meanl = min(intUp_meanl, intDown_meanl);
   double intMax_meang = max(intUp_meang, intDown_meang);
   double intMin_meang = min(intUp_meang, intDown_meang);
   //  
   double intMax_sigmal = max(intUp_sigmal, intDown_sigmal);
   double intMin_sigmal = min(intUp_sigmal, intDown_sigmal);
   double intMax_sigmag = max(intUp_sigmag, intDown_sigmag);
   double intMin_sigmag = min(intUp_sigmag, intDown_sigmag);

   //3.3 for each bin, you sum the differences of all max values to the central value in quadrature
                // intMax2Sum = (intMax-int)^2
                // and same for min
                // intMin2Sum = (intMin-int)^2 
   double intMax2Sum = (pow((intMax_meanl-CentralValue),2) +
                        pow((intMax_meang-CentralValue),2) +
                        pow((intMax_sigmal-CentralValue),2)+ 
                        pow((intMax_sigmag-CentralValue),2) )  ;

   double intMin2Sum = (pow((intMin_meanl-CentralValue),2) +
                        pow((intMin_meang-CentralValue),2) +
                        pow((intMin_sigmal-CentralValue),2)+ 
                        pow((intMin_sigmag-CentralValue),2) )  ;


    //3.4 the upper envelope of the smoother template is then computed bin-by-bin as
               //smoothedTemplateUp = smoothedTemplate + sqrt(intMax2Sum)
               //smoothedTemplateDown = smoothedTemplate - sqrt(intMin2Sum)
   double smoothedTemplateUp = CentralValue + sqrt(intMax2Sum);
   double smoothedTemplateDown = CentralValue - sqrt(intMin2Sum);

  }

*/

/////////////////////////////////////////////

