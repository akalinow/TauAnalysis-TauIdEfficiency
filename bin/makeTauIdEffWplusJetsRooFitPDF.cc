
/** \executable makeTauIdEffQCDtemplate
 *
 * Obtain QCD template from control region in Data,
 * corrected for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: makeTauIdEffWplusJetsRooFitPDF.cc,v 1.4 2012/06/15 15:16:17 calpas Exp $
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"

using namespace RooFit;

void FitHisto(TH1* htest){
  
  std::string histoNameInput = htest->GetName();
  std::string histoName      = htest->GetName();
  
  double binPDF = htest->GetXaxis()->GetNbins();
  
  //General Variable      
  RooRealVar x("x","x",0,200);

  //////////////Landau convoluted with Gauss function//////////////
  //LG1
  //Landau
  RooRealVar meanl("meanl","mean of Landau",80.,60,90);
  RooRealVar sigmal("sigmal","sigma of Landau",30,20,50);  
  RooLandau landau("landau","landau",x,meanl,sigmal);
  //Gauss
  RooRealVar meang("meang","mean of Gaussian",1); 
  RooRealVar sigmag("sigmag","sigma of Gaussian",1,0.1,10);
  RooGaussian gauss("gauss","gauss",x,meang,sigmag);
  //Landau convoluted Gauss
  RooFFTConvPdf LandauConvGauss("LandauConvGauss","LandauConvGauss",x,landau,gauss);
  

  //////////////Cristal Ball fucntion/////////////////////////////
  
  //CB1
  RooRealVar cbmean1("cbmean1", "cbmean1" , 90.0, 20, 180.0);
  RooRealVar cbsigma1("cbsigma1", "cbsigma1" , 1, 1.0, 40.0); 
  RooRealVar cbsig1("cbsig1", "cbsignal", 800, 20, 80);
  RooRealVar n1("n1","", 0.2); 
  RooRealVar alpha1("alpha1","", 1.3);
  RooCBShape CristalBall1("CristalBall1","CristalBall1",x,cbmean1,cbsigma1,alpha1,n1);

  //////////////////////// 
  
  //CB2 
  RooRealVar cbmean2("cbmean2", "cbmean2" , 70, 20, 180.0);
  RooRealVar cbsigma2("cbsigma2", "cbsigma2" , 10, 1.0, 200.0); 
  RooRealVar cbsig2("cbsig2", "cbsignal2", 1, 20, 80);
  RooRealVar n2("n2","", 1); 
  RooRealVar alpha2("alpha2","", 1); 
  RooCBShape CristalBall2("CristalBall2","CristalBall2",x,cbmean2,cbsigma2,alpha2,n2);
  
  ////////////////////////
  
  //CB3 
  RooRealVar cbmean3("cbmean3", "cbmean3" , 1, 20, 180.0);
  RooRealVar cbsigma3("cbsigma3", "cbsigma3" , 10, 1.0, 200.0); 
  RooRealVar cbsig3("cbsig3", "cbsignal3", 1, 20, 80);
  RooRealVar n3("n3","", 1); 
  RooRealVar alpha3("alpha3","", 1);
  RooCBShape CristalBall3("CristalBall3","CristalBall3",x,cbmean3,cbsigma3,alpha3,n3);
  
  /////////////////////////////////////////////////////////////////
  
  
  //Make hist as a abspdf to be fit
  RooDataHist data("","",x,htest, 1.0); 
  
  RooPlot* frame = x.frame(Title("Fit"));
  
  const char *histoToFitName;
  histoToFitName = histoNameInput.c_str();

  string histoToFitNameSTR = std::string(histoName).data(); 

  const char *histoToSmooth1 = "WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake"; 
  const char *histoToSmooth2 = "WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake";  
  const char *histoToSmooth3 = "WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake";
  const char *histoToSmooth4 = "WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake";
  const char *histoToSmooth5 = "QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all";
  const char *histoToSmooth6 = "QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all";

  // convert CHAR to STR to compare
  string histoToSmooth1STR (histoToSmooth1);
  string histoToSmooth2STR (histoToSmooth2);
  string histoToSmooth3STR (histoToSmooth3);
  string histoToSmooth4STR (histoToSmooth4);
  string histoToSmooth5STR (histoToSmooth5);
  string histoToSmooth6STR (histoToSmooth6);

  //histo for smooth function
  TH1* hist;  

  if( histoToFitNameSTR == histoToSmooth1STR ){
    LandauConvGauss.fitTo(data);
    data.plotOn(frame);
    LandauConvGauss.plotOn(frame,LineColor(kRed)); 
    hist = LandauConvGauss.createHistogram(std::string(histoName).append("_smoothed").data(), x,Binning(binPDF) );  
 }
  
  else if ( histoToFitNameSTR == histoToSmooth2STR){
    CristalBall1.fitTo(data);
    data.plotOn(frame);
    CristalBall1.plotOn(frame,LineColor(kRed)); 
    hist = CristalBall1.createHistogram(std::string(histoName).append("_smoothed").data(),x,Binning(binPDF)); //Bin, min, max
  }
  
  else if (( histoToFitNameSTR == histoToSmooth3STR) ||
	   ( histoToFitNameSTR == histoToSmooth4STR )){
    CristalBall2.fitTo(data); 
    data.plotOn(frame);
    CristalBall2.plotOn(frame,LineColor(kRed)); 
    hist = CristalBall2.createHistogram(std::string(histoName).append("_smoothed").data(),x,Binning(binPDF)); //Bin, min, max
  }
  
  else if (( histoToFitNameSTR == histoToSmooth5STR) ||
	   ( histoToFitNameSTR == histoToSmooth6STR )){
    CristalBall3.fitTo(data);
    data.plotOn(frame);
    CristalBall3.plotOn(frame,LineColor(kRed)); 
    hist = CristalBall3.createHistogram(std::string(histoName).append("_smoothed").data(),x,Binning(binPDF)); //Bin, min, max
  }
  
  
  //Draw frame on canvas
  TCanvas *c = new TCanvas("","Fit",600,600) ;
  gPad->SetLeftMargin(0.15) ; 
  frame->GetYaxis()->SetTitleOffset(1.4) ;
  frame->Draw();
  
  //save canvas in pdf eps and png    
  TString fileNamePDF = std::string(histoName).append("_smoothed.pdf").data();
  TString fileNameEPS = std::string(histoName).append("_smoothed.eps").data();
  TString fileNamePNG = std::string(histoName).append("_smoothed.png").data();
  c->Print(fileNamePDF);
  c->Print(fileNameEPS);
  c->Print(fileNamePNG);
  
  //Get the Pdf of the fit, fill it in a histogram and save
  TFile f(std::string(histoName).append("_smoothed.root").data(),"RECREATE"); 

  //set Bin histoPDF to be the same then the template histogram  
  // x.setBins( htest->Integral());

  hist->Write();
  f.Write();
  f.Close();
}



int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<makeTauIdEffQCDtemplate>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeTauIdEffQCDtemplate");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgMakeTauIdEffQCDtemplate = cfg.getParameter<edm::ParameterSet>("fitTauIdEff"); 
 
  //histogram to smooth
  vstring vHistoToFit = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("HistogramsToFit"); 

  //function to use for the smoothing
  vstring vFitFunctions = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("FitFunctions"); 
 
   fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgMakeTauIdEffQCDtemplate.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory histogramOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;


////////check histo and fit name for fit///////

  map< string, TH1* > htest;

  for ( vstring::const_iterator histogramName = vHistoToFit.begin();
       histogramName != vHistoToFit.end(); ++histogramName ) {

       //Get the Hitsto to be fit 
       htest[histogramName->data()] = (TH1F*)histogramInputFile->Get(histogramName->data());
 
       FitHisto(htest[histogramName->data()] );
   } 


/*
  for ( vstring::const_iterator histogramName = vHistoToFit.begin();
       histogramName != vHistoToFit.end(); ++histogramName ) {

    for ( vstring::const_iterator FitToUse = vFitFunctions.begin();
         FitToUse != vFitFunctions.end(); ++histogramName ) {

      if ( ){

        //Get the Hitsto to be fit 
        htest[histogramName->data()] = (TH1F*)histogramInputFile->Get(histogramName->data());
 
       FitHisto(htest[histogramName->data()] );
       }
     }
   } 
*/

  return 0;
}
