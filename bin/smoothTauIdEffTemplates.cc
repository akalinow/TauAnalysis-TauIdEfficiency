
/** \executable smoothTauIdEffTemplates
 *
 * Fit template histograms used in tau id. efficiency measurement by analytic functions,
 * in order to avoid bias on results arising from statistical fluctuations of template shapes
 *
 * \author Betty Calpas, RWTH Aachen
 *         Christian Veelken, LLR
 *
 * \version $Revision: 1.8 $
 *
 * $Id: smoothTauIdEffTemplates.cc,v 1.8 2012/10/12 09:24:34 kkaadze Exp $
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
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "RooExponential.h"
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
#include "TMatrix.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include <TVectorD.h>
#include <TVectorF.h>
#include <TVector.h>
#include <TVector3.h>
#include "TIterator.h"
#include <math.h>
#include <map>

using namespace RooFit;
using namespace std;

void addBinsNotFitted(TH1* histogram_fitted, const TH1* histogram, double xMin_fit, double xMax_fit)
{
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binCenter = xAxis->GetBinCenter(iBin);
    if ( binCenter >= xMin_fit && binCenter < xMax_fit ) {
      histogram_fitted->GetBinContent(iBin, histogram->GetBinContent(iBin));
      histogram_fitted->GetBinError(iBin, histogram->GetBinError(iBin));
    }
  }
}

double integral(const TH1* histogram, double xMin, double xMax)
{
  double integral = 0.;
  TAxis* xAxis = histogram->GetXaxis();
  int numBins = xAxis->GetNbins();
  for ( int iBin = 1; iBin <= numBins; ++iBin ) {
    double binCenter = xAxis->GetBinCenter(iBin);
    double binContent = histogram->GetBinContent(iBin);
    if ( binCenter >= xMin && binCenter < xMax ) integral += binContent;
  }
  return integral;
}

void smoothHistogram(TH1* histogram, const std::string& fitFunctionType, 
		     TFileDirectory& histogramOutputDirectory,
		     bool makeControlPlots, const std::string& controlPlotFilePath, 
		     bool limit_xMin_fit, double xMin_fit, 
		     bool limit_xMax_fit, double xMax_fit)
{  
  TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  if ( limit_xMin_fit ) xMin = TMath::Max(xMin, xMin_fit);
  double xMax = xAxis->GetXmax();
  if ( limit_xMax_fit ) xMax = TMath::Min(xMax, xMax_fit);
  
  // general fit variable
  RooRealVar x("x", "x", xMin, xMax);
  
  RooAbsPdf* fitFunction = 0;
  std::vector<TObject*> objectsToDelete;

  if ( fitFunctionType == "LG1" ) {
    // create Landau function
    RooRealVar*  meanl  = new RooRealVar ("meanl",  "mean of Landau",  80.,  60.   , 90.    );
    RooRealVar*  sigmal = new RooRealVar ("sigmal", "sigma of Landau", 30.,  20.   , 50.    );  
    RooLandau*   landau = new RooLandau  ("landau", "landau",            x,  *meanl, *sigmal);
    // create Gaussian
    RooRealVar*  meang  = new RooRealVar ("meang",  "mean of Gaussian",  1., 0.    , 15.    ); 
    RooRealVar*  sigmag = new RooRealVar ("sigmag", "sigma of Gaussian", 1., 0.1   , 10.    );
    RooGaussian* gauss  = new RooGaussian("gauss",  "gauss",             x,  *meang, *sigmag);
    // create convolution of Landau with Gaussian
    fitFunction = new RooFFTConvPdf("LandauConvGauss", "LandauConvGauss", x, *landau, *gauss);

    objectsToDelete.push_back(meanl);
    objectsToDelete.push_back(sigmal);
    objectsToDelete.push_back(landau);
    objectsToDelete.push_back(meang);
    objectsToDelete.push_back(sigmag);
    objectsToDelete.push_back(gauss);
  } else if (fitFunctionType == "EXP1" ){
    // create Expo
    RooRealVar*  lambda = new RooRealVar ("lambda", "slope", -1., -2. , 10.);  
    fitFunction = new RooExponential("expo", "exponential PDF", x, *lambda);
    objectsToDelete.push_back(lambda);
  } 

 else if (fitFunctionType == "EXP2") {
   //   RooAbsArg* par0= new RooRealVar ("par0","par0",0.5);
   RooAbsArg* par1= new RooRealVar ("par1","par1",-0.00001,-0.1,-0.0000000001);
   RooAbsArg* par2= new RooRealVar ("par2","par2",1,0.1,3);
   x.setRange(200,1500);
   fitFunction = new RooGenericPdf("g","TMath::Exp(par1*TMath::Power(x,par2))",RooArgList(x,*par1,*par2));
   objectsToDelete.push_back(par1);
   objectsToDelete.push_back(par2);
 }

  
 else if ( fitFunctionType == "CB1" ||
	   fitFunctionType == "CB2"   ){
   // create Crystal-ball function
   RooRealVar* cbmean  = 0;
   RooRealVar* cbsigma = 0;
   RooRealVar* n       = 0;
   RooRealVar* alpha   = 0;
   if        ( fitFunctionType == "CB1" ) {
      cbmean  = new RooRealVar("cbmean",  "cbmean",  90., 20., 180.);
      cbsigma = new RooRealVar("cbsigma", "cbsigma", 1.,  1.,  40. ); 
      n       = new RooRealVar("n",       "n",       0.2, 0.,  10. ); 
      alpha   = new RooRealVar("alpha",   "alpha",   1.3, 0.,  10. );
   } else if ( fitFunctionType == "CB2" ) {
     cbmean  = new RooRealVar("cbmean",  "cbmean",  50., 20., 180.);
     cbsigma = new RooRealVar("cbsigma", "cbsigma", 5., 1.,  200.); 
     n       = new RooRealVar("n",       "n",       1.,  0.,  10. ); 
     alpha   = new RooRealVar("alpha",   "alpha",   1.,  -10.,  10. ); 
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
  
  //Get FitFunction result (mean, sigma...)
  RooFitResult* r = fitFunction->fitTo(data,Save()) ; 
  cout << "Fit status " << endl;
  cout << r->status() << endl;
  cout<< "Fit results: " << endl;
  r->Print();
  
  //Get the full covariance matrix
  const TMatrixDSym& cov = r->covarianceMatrix() ;
  cout << endl << "Cov Matrix: "<< endl;
  cov.Print();

  int numFitParameter = cov.GetNrows();
 
  //Find the Eigenvectors and Eigen value of the covariance matrix
  const TMatrixDSymEigen& eigencov = cov;  

  //Get the matrix of Eigen vector
  const TMatrixD eigenvector_matrix = eigencov.GetEigenVectors();
  cout << endl << "Matrix of Eigenvectors"<< endl;
  eigenvector_matrix.Print();

  const TVectorD eigenvalues = eigencov.GetEigenValues();
  cout << endl << "Vector of Eigenvalues"<< endl;
  eigenvalues.Print();

  std::vector<TVectorD*> eigenvectors;
  for ( int i = 0; i < eigenvector_matrix.GetNrows(); ++i ) {
    TVectorD* eigenvector = new TVectorD(eigenvector_matrix.GetNrows());
    for ( int j = 0; j < eigenvector_matrix.GetNrows(); ++j ) {
      (*eigenvector)(j) =  eigenvector_matrix(j, i);
    }
    eigenvectors.push_back(eigenvector);
  }
  
  // check #1: M( eigenvectorInv_matrix * eigenvector_matrix ) = M(Identity) 
  TMatrixD eigenvectorInv_matrix(TMatrixD::kInverted, eigenvector_matrix);  
  const TMatrixD ID = eigenvectorInv_matrix*eigenvector_matrix;
  for ( int i = 0; i < numFitParameter; ++i ) {
    for ( int j = 0; j < numFitParameter; ++j ) {
      if ( i == j ) assert(TMath::Abs(ID(i, j) - 1.) < 1.e-3);
      else assert(TMath::Abs(ID(i, j)) < 1.e-3);
    }
  }

  // check #2: eigenvectorInv_matrix * cov * eigenvector_matrix = diagonal matrix of eigenvalues 
  // http://root.cern.ch/phpBB3/viewtopic.php?f=15&t=8663 
  const TMatrixD MDiag = eigenvectorInv_matrix * (cov * eigenvector_matrix); 
  for ( int i = 0; i < numFitParameter; ++i ) {
    for ( int j = 0; j < numFitParameter; ++j ) {
      if ( i == j ) assert(TMath::Abs(MDiag(i, j) - eigenvalues(i)) < (1.e-3*eigenvalues(i)));
      else assert(TMath::Abs(ID(i, j)) < 1.e-3);
    }
  }

  // check #3: cov * eigenvector = eigenvalue * eigenvector 
  for ( int i = 0; i < numFitParameter; ++i ) {
    const TVectorD& eigenvector = (*eigenvectors[i]);
    TVectorD cov_times_eigenvector = cov * eigenvector;
    for ( int j = 0; j < numFitParameter; ++j ) {
      assert(TMath::Abs(cov_times_eigenvector(j) - eigenvalues[i]*eigenvector(j)) < (1.e-3*eigenvalues[i]));
    }
  }

  //Store function's parameter 
  RooArgSet* params  = fitFunction->getParameters(x);
  vector<std::string> vfitFunction_parName;  
  TIterator* it_test = params->createIterator();
  while ( RooRealVar* par_test = (RooRealVar*) it_test->Next() ) {
    vfitFunction_parName.push_back(par_test->GetName());
  }
  
  //define mapping of functionType to fitParameterNames
  std::map<std::string, vstring> fitFunction_map;
  map<std::string, vstring>::iterator map_it;
  
  fitFunction_map.insert(std::pair<std::string, vstring>(fitFunctionType, vfitFunction_parName));
  
  //save parameter value to be retrive later
  //Converting a string to an array of characters with string.data(): http://msdn.microsoft.com/en-us/library/3372cxcy.aspx
  std::vector<double> param_values;
  for ( int i = 0; i < numFitParameter; ++i ) {
    RooRealVar* param = (RooRealVar*)params->find(vfitFunction_parName[i].data()); 
    double param_value = param->getVal(); 
    param_values.push_back(param_value);
  }
  
  for( int i = 0; i < numFitParameter; ++i ) {
    
    //Get the i-th Eigenvector
    const TVectorD& eigenvector = (*eigenvectors[i]);
    cout<< "Eigenvector[" << i << "]:" << endl;
    eigenvector.Print();
    //Compute norm of i-th Eigenvector
    double norm = eigenvector.Norm2Sqr();
    
    //Get the Eigenvalue associate to i-th Eigenvector
    double eigenvalue = eigenvalues[i];
    cout<< "Eigenvalue[i] = " << eigenvalue << std::endl;
    //Compute uncertainty on fit parameter
    double sigma = sqrt(eigenvalue);
    //Compute unit vector in direction of i-th Eigenvector
    const TVectorD direction = (1./norm)*eigenvector;
    //cout << "direction: " << direction[i] << endl;
    
    //save i eigenvec index for histo name
    std::ostringstream strsI;
    strsI << i;
    std::string Iindex = strsI.str();
    
    //Vary fit parameters in "Up" direction
    TIterator* itUp = params->createIterator();  
    int idxUp = 0;
    while ( RooRealVar* param = (RooRealVar*)itUp->Next() ) { 
      cout << endl << "param_value Up: " << param->getVal() << endl; 
      param->setVal(param->getVal() + sigma*direction[idxUp]);
      cout << endl << "param_value Up + err: " << param->getVal() << endl; 
      ++idxUp;
    }
    
    //save Up syst. histogram
    std::string histogramName_smoothed_Eigenvec_up = std::string(histogram->GetName()).append("_smoothed_").append("_EigenVec").append(Iindex).append("_up");
    TH1* histogram_smoothed_Eigenvec_up = fitFunction->createHistogram(histogramName_smoothed_Eigenvec_up.data(), x, Binning(histogramBinning));
    addBinsNotFitted(histogram_smoothed_Eigenvec_up, histogram, xMin, xMax);
    
    // Restore central values of fit parameters
    itUp = params->createIterator();  
    idxUp = 0;
    while ( RooRealVar* param = (RooRealVar*)itUp->Next() ) { 
      param->setVal(param_values[idxUp]);
      ++idxUp;
    }
    
    //Vary fit parameters in "Down" direction
    TIterator* itDown = params->createIterator();  
    int idxDown = 0;
    while ( RooRealVar* param = (RooRealVar*) itDown->Next() ) { 
      cout << endl << "param_value Down: " << param->getVal() << endl; 
      param->setVal(param->getVal() - sigma*direction[idxDown] );
      cout << endl << "param_value Down - err: " << param->getVal() << endl; 
      ++idxDown;
    }
    
    //save Down syst. histogram
    std::string histogramName_smoothed_Eigenvec_down = std::string(histogram->GetName()).append("_smoothed_").append("_EigenVec").append(Iindex).append("_down");
    TH1* histogram_smoothed_Eigenvec_down = fitFunction->createHistogram(histogramName_smoothed_Eigenvec_down.data(), x, Binning(histogramBinning));
    addBinsNotFitted(histogram_smoothed_Eigenvec_down, histogram, xMin, xMax);
    
    // Restore central values of fit parameters
    itDown = params->createIterator();  
    idxDown = 0;
    while ( RooRealVar* param = (RooRealVar*)itDown->Next() ) { 
      param->setVal(param_values[idxDown]);
      ++idxDown;
    }
  }
  
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
  histogram_smoothed->Scale(integral(histogram, xMin, xMax)/integral(histogram_smoothed, xMin, xMax));
  
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
  bool xMin_limited_;
  double xMin_;
  bool xMax_limited_;
  double xMax_;
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
    if ( histogramToSmooth->exists("xMin") ) {
      inputHistogramEntry.xMin_limited_ = true;
      inputHistogramEntry.xMin_ = histogramToSmooth->getParameter<double>("xMin");
    } else {
      inputHistogramEntry.xMin_limited_ = false;
      inputHistogramEntry.xMin_ = 0.;
    }
    if ( histogramToSmooth->exists("xMax") ) {
      inputHistogramEntry.xMax_limited_ = true;
      inputHistogramEntry.xMax_ = histogramToSmooth->getParameter<double>("xMax");
    } else {
      inputHistogramEntry.xMax_limited_ = false;
      inputHistogramEntry.xMax_ = 0.;
    }
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
		    makeControlPlots, controlPlotFilePath,
		    inputHistogramEntry->xMin_limited_, inputHistogramEntry->xMin_,
		    inputHistogramEntry->xMax_limited_, inputHistogramEntry->xMax_);
  }
  
  return 0;
}






