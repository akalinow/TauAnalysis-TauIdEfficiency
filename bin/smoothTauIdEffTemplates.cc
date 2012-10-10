
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
 * $Id: smoothTauIdEffTemplates.cc,v 1.5 2012/09/05 06:33:44 calpas Exp $
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

  } 


 else if (fitFunctionType == "EXP1" ){
    // create Expo
    RooRealVar*  lambda = new RooRealVar ("lambda", "slope", -1., -2. , 10.);  
    fitFunction = new RooExponential("expo", "exponential PDF", x, *lambda);
    objectsToDelete.push_back(lambda);
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
  
  
  cout << endl<< "XXXXXXXXXXXXXXXXXXXX Syst study XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl << endl;

  //Get FitFunction result (mean, sigma...)
  RooFitResult* r = fitFunction->fitTo(data,Save()) ; 
  cout<< "Fit results: " << endl;
  r->Print();
  
  //Get the full covariance matrix
  const TMatrixDSym& cov = r->covarianceMatrix() ;
  cout << endl << "Cov Matrix: "<< endl;
  cov.Print();
 
  //Find the Eigenvectors and Eigen value of the covariance matrix
  const TMatrixDSymEigen& eigencov = cov;  

  //Get the matrix of Eigen vector
  const TMatrixD eigenvec = eigencov.GetEigenVectors();
  cout << endl << "Matrix of eigenvectors"<< endl;
  eigenvec.Print();

  const TVectorD eigenval = eigencov.GetEigenValues();
  cout << endl << "Vector of Eigen value"<< endl;
  eigenval.Print();

  //The eigenvectors correspond to the column of the eigenvec Matrix  !!OK!! 
  cout << "Eigenvec Matrix column check: each eigenvector should match a eigenvec Matrix column!!!! "<< endl; 
  TVectorD EV0(4); // need to specify the vector size!!
  TVectorD EV1(4); 
  TVectorD EV2(4); 
  TVectorD EV3(4); 

  if(eigenvec.GetNrows() < 2){
    for(int i = 0 ; i < eigenvec.GetNrows(); ++i){
     EV0[i] = eigenvec(i,0);
     cout << "Column EVO: " << i << "  " << EV0[i] << endl;
    }
  }else if(eigenvec.GetNrows() < 5){
    for(int i = 0 ; i < eigenvec.GetNrows(); ++i){
     EV0[i] = eigenvec(i,0);
     EV1[i] = eigenvec(i,1);
     EV2[i] = eigenvec(i,2);
     EV3[i] = eigenvec(i,3);
     cout << "Column EVO: " << i << "  " << EV0[i] << endl;
     cout << "Column EV1: " << i << "  " << EV1[i] << endl;
     cout << "Column EV2: " << i << "  " << EV2[i] << endl;
     cout << "Column EV3: " << i << "  " << EV3[i] << endl;
    }
  }

  //check n1: M( invert eigenvec * eigenvec ) = M(Identity) !!OK!!
  TMatrixD inveigenvec(TMatrixD::kInverted,eigenvec);  
  const TMatrixD ID = inveigenvec * eigenvec;
  cout<< endl << "Matrice identite "<< endl;
  ID.Print(); // should print Identity Matrix 

  //Check n2: inveigenvec * cov * eigenvec = diag. of eigenval !!OK!!
  // http://root.cern.ch/phpBB3/viewtopic.php?f=15&t=8663 
  //TMatrixD Minv(TMatrixD::kInverted,eigenvec);  
  const TMatrixD MDiag = inveigenvec * cov * eigenvec; 
  cout << endl << " inveigenvec * cov * eigenvec =? diag. of eigenval" << endl;
  MDiag.Print();  //should print a Matrix with the eigenvalue on the diagonal 


  //Store function's parameter 
  RooArgSet* params  = fitFunction->getParameters(x);
  //Create iterator
  TIterator* it_test = params->createIterator() ;
  RooRealVar* par_test;
    
  vector<std::string> vfitFunction_parName;
  
  while(par_test = (RooRealVar*) it_test->Next()){
    vfitFunction_parName.push_back(par_test->GetName());
  }
  
  int size = vfitFunction_parName.size();
  //for (int i = 0; i < size ; ++i ){ 
  // cout << "parName: "<< vfitFunction_parName[i] << endl;
  //}

   //define mapping of functionType to fitParameterNames
   std::map<std::string, vstring> fitFunction_map;
   map<std::string, vstring>::iterator map_it;

   fitFunction_map.insert( pair< std::string, vstring > (fitFunctionType, vfitFunction_parName ));

   //for ( map_it = fitFunction_map.begin() ; map_it != fitFunction_map.end(); map_it++ ){
   //   cout << endl << (*map_it).first << " MAP " << endl; 
   //   cout << (*map_it).second << endl;  // !!error message!! 
   //}

   
   //save parameter value to be retrive later
   //Converting a string to an array of characters with string.data(): http://msdn.microsoft.com/en-us/library/3372cxcy.aspx

   double par0_value = 0;
   double par1_value = 0;
   double par2_value = 0;
   double par3_value = 0;
   RooRealVar* par0 ;
   RooRealVar* par1 ;
   RooRealVar* par2 ;
   RooRealVar* par3 ;

   if(size < 2){ 
   par0 = (RooRealVar*) params->find(vfitFunction_parName[0].data()); 
   par0_value = par0->getVal(); 
   //cout << "par0: " << par0_value << endl;
   }else if(size > 0 && size < 5){
   par0 = (RooRealVar*) params->find(vfitFunction_parName[0].data()); 
   par0_value = par0->getVal(); 
   //cout << "par0: " << par0_value << endl;
   par1 = (RooRealVar*) params->find(vfitFunction_parName[1].data());
   par1_value = par1->getVal(); 
   //cout << "par1: " << par1_value << endl;
   par2 = (RooRealVar*) params->find(vfitFunction_parName[2].data());
   par2_value = par2->getVal(); 
   //cout << "par2: " << par2_value << endl;
   par3 = (RooRealVar*) params->find(vfitFunction_parName[3].data());
   par3_value = par3->getVal(); 
   //cout << "par3: " << par3_value << endl;
   }
   

   //loop over eigenvector
   //http://root.cern.ch/root/html400/TVectorD.html#TVectorD:kSizeMax
   //cout << eigenval.GetNrows()  << endl; 
   //cout << eigenvec.GetNrows() << endl;

   for( int i = 0; i < eigenvec.GetNrows(); ++i ){
     cout<< "Eigen Vector: "<< i << endl;
     
     //Get the Eigenvector[i]
     const TVectorD Eigenvec = eigenvec[i]; 
     //Get the Eigenvector[i] norm
     double Norm = Eigenvec.Norm2Sqr();
     //double Norm1 = sqrt(Eigenvec * Eigenvec); //produce the same result then above 
     //cout << "Norm: "<< Norm << endl; 
     //cout << "Norm1: "<< Norm1 << endl; 
     
     //Get the Eigenvalue(i) associate to the eigenvector[i]
     double Eigenval = eigenval[i];
     //Get the error on Eigenval
     double Sigma = sqrt(Eigenval);
     //Get unit vector: projection of Eigenvector[i] over it's direction (i), divided by it's norm
     const TVectorD Direction = (1./Norm) * Eigenvec;
     //cout << "Direction: " << Direction[i] << endl;
     
     
     //save i eigenvec index for histo name
     std::ostringstream strsi;
     strsi << i;
     std::string Iindex = strsi.str();
     //cout << endl << "Eigenvect index: " << Iindex << endl;
  
   
     for( int j = 0; j < eigenval.GetNrows() ; ++j ){
       cout << "Eigen Value: " << j << endl << endl;	
       
       //Create iterator over the parameters
       TIterator* itUp    = params->createIterator() ;  // !!do not work if only one!!??
       TIterator* itDown  = params->createIterator() ;  
       RooRealVar* par;

      
       cout <<  "##################### Up ##########################"<< endl;

       //Varie all the fit parameters in on time in the direction of 1 Eigen Vector at a time 
       while(par = (RooRealVar*) itUp->Next()){ // !! no good in for Loop !!??
	 //Varie par Up
	 //par->Print();
	 cout << endl << "par_value Up: " << par->getVal() << endl; 
	 par->setVal(par->getVal() + Sigma*Direction[j] );
	 //cout << "+ err: " <<  Sigma*Direction[j] << endl; 
	 cout << endl << "par_value Up + err: " << par->getVal() << endl; 
	 //cout << "sigma: " << Sigma << endl;
	 //cout << "direction: " << Direction[j] << endl << endl;
       }
       
       
       //save j direction index for histo name
       std::ostringstream strsj;
       strsj << j;
       std::string Jindex = strsj.str();
       //cout <<"Direction index: " << Jindex << endl;
       
       //save Up syst. histogram
       std::string histogramName_smoothed_Eigenvec_up = std::string(histogram->GetName()).append("_smoothed_").append("_EigenVec").append(Iindex).append("_Direction").append(Jindex).append("_up");
       TH1* histogram_smoothed_Eigenvec_up = fitFunction->createHistogram(histogramName_smoothed_Eigenvec_up.data(), x, Binning(histogramBinning));
  
     
       //Retrive parameter value for down syst
       if(size < 2){ 
       par0->setVal(par0_value);
       }
       else if(size > 1 && size < 5){ 
       par0->setVal(par0_value);
       par1->setVal(par1_value);
       par2->setVal(par2_value);
       par3->setVal(par3_value);
       }

       cout <<  "##################### Down ##########################" << endl;

       //Varie par Down
       while(par = (RooRealVar*) itDown->Next()){ 
	 cout << endl << "par_value Down: " << par->getVal() << endl; 
	 par->setVal(par->getVal() - Sigma*Direction[j] );
	 //cout << "- err: " <<  Sigma*Direction[j] << endl; 
	 cout << endl << "par_value Down - err: " << par->getVal() << endl; 
	 //cout << "sigma: " << Sigma << endl;
	 //cout << "Direction[j]: " << Direction[j] << endl;
       }
		
       //save Down syst. histogram
       std::string histogramName_smoothed_Eigenvec_down = std::string(histogram->GetName()).append("_smoothed_").append("_EigenVec").append(Iindex).append("_Direction").append(Jindex).append("_down");
       TH1* histogram_smoothed_Eigenvec_down = fitFunction->createHistogram(histogramName_smoothed_Eigenvec_down.data(), x, Binning(histogramBinning));
       
       //Retrive parameter value

       if(size < 2){ 
       par0->setVal(par0_value);
       }
       else if(size > 1 && size < 5){ 
       par0->setVal(par0_value);
       par1->setVal(par1_value);
       par2->setVal(par2_value);
       par3->setVal(par3_value);
       }

     }  
   }
   
   cout << endl <<  "XXXXXXXXXXXXXXXXXXXX End syst study XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxXXXXXXX"<<endl;  
   
   
   
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






