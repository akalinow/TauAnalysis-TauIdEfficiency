
/** \executable makeTauIdEffQCDtemplate
 *
 * Obtain QCD template from control region in Data,
 * corrected for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: makeTauIdEffWplusJetsRooFitPDF.cc,v 1.2 2012/06/15 08:07:06 calpas Exp $
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


typedef std::map<std::string, std::map<std::string, TH1*> > histogramMap;
typedef std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsMap;

struct regionEntryType
{
  regionEntryType(const std::string& regionTakeQCDtemplateFromData, 
		  const std::string& regionWplusJetsSideband, 
		  const std::string& regionStoreQCDtemplate, 
		  const std::string& tauId, const std::string& fitVariable, const std::string& sysUncertainty)
    : regionTakeQCDtemplateFromData_(regionTakeQCDtemplateFromData),
      regionWplusJetsSideband_(regionWplusJetsSideband),
      regionStoreQCDtemplate_(regionStoreQCDtemplate),
      tauId_(tauId),
      fitVariable_(fitVariable),
      sysUncertainty_(sysUncertainty)
  {}
  ~regionEntryType() {}
  
  void makeQCDtemplate(TFileDirectory& dir, 
		       histogramMap3& histograms_data, histogramMap4& histograms_mc, valueMap3& numEvents_mc)
  {
    //std::cout << "<regionEntryType::makeQCDtemplate>" << std::endl;
    //std::cout << " regionWplusJetsSideband = " << regionWplusJetsSideband_ << std::endl;
    //std::cout << " tauId = " << tauId_ << std::endl;
    //std::cout << " fitVariable = " << fitVariable_ << std::endl;
    //std::cout << " sysUncertainty = " << sysUncertainty_ << std::endl;
    
    double numEventsWplusJetsSideband_data = 
      getIntegral(histograms_data[regionWplusJetsSideband_]["EventCounter"][key_central_value], true, true);
    double numEventsWplusJetsSideband_expBgr = 
      (numEvents_mc["ZplusJets"][regionWplusJetsSideband_][sysUncertainty_]
       + numEvents_mc["QCD"][regionWplusJetsSideband_][sysUncertainty_]
       + numEvents_mc["TTplusJets"][regionWplusJetsSideband_][sysUncertainty_]);
    double numEventsWplusJetsSideband_obsWplusJets = numEventsWplusJetsSideband_data - numEventsWplusJetsSideband_expBgr;
    double numEventsWplusJetsSideband_expWplusJets = 
      numEvents_mc["WplusJets"][regionWplusJetsSideband_][sysUncertainty_];    
    //std::cout << "WplusJets sideband: observed = " << numEventsWplusJetsSideband_data << ","
    //	        << " expected background = " << numEventsWplusJetsSideband_expBgr 
    //	        << " --> observed WplusJets contribution = " << numEventsWplusJetsSideband_obsWplusJets << ","
    //	        << " expected = " << numEventsWplusJetsSideband_expWplusJets << std::endl;
    
    double scaleFactorWplusJets = numEventsWplusJetsSideband_obsWplusJets/numEventsWplusJetsSideband_expWplusJets;
    
    TH1* distributionDataQCDsideband = histograms_data[regionTakeQCDtemplateFromData_][fitVariable_][key_central_value];
    TH1* templateZtautauQCDsideband = histograms_mc["ZplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateWplusJetsQCDsideband = histograms_mc["WplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateTTplusJetsQCDsideband = histograms_mc["TTplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    
    TString templateQCDsidebandName = distributionDataQCDsideband->GetName();    
    templateQCDsidebandName.ReplaceAll(regionTakeQCDtemplateFromData_, regionStoreQCDtemplate_);
    if ( sysUncertainty_ != key_central_value ) templateQCDsidebandName.Append("_").Append(sysUncertainty_);
    std::cout << "creating histogram = " << templateQCDsidebandName << std::endl;
    TH1* templateQCDsideband_obsQCD = 0;
    TAxis* xAxis = distributionDataQCDsideband->GetXaxis(); 
    if ( xAxis->GetXbins() && xAxis->GetXbins()->GetSize() >= 2 ) {
      const TArrayD* binning = xAxis->GetXbins();      
      templateQCDsideband_obsQCD = 
	dir.make<TH1D>(templateQCDsidebandName, 
		       distributionDataQCDsideband->GetTitle(), binning->GetSize() - 1, binning->GetArray());
    } else {
      templateQCDsideband_obsQCD = 
	dir.make<TH1D>(templateQCDsidebandName, 
		       distributionDataQCDsideband->GetTitle(), xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax());
    }
    templateQCDsideband_obsQCD->Add(distributionDataQCDsideband, +1.);
    templateQCDsideband_obsQCD->Add(templateZtautauQCDsideband, -1.);
    templateQCDsideband_obsQCD->Add(templateWplusJetsQCDsideband, -1.*scaleFactorWplusJets);
    templateQCDsideband_obsQCD->Add(templateTTplusJetsQCDsideband, -1.);
  }

  std::string regionTakeQCDtemplateFromData_;
  std::string regionWplusJetsSideband_;
  std::string regionStoreQCDtemplate_;

  std::string tauId_;
  std::string fitVariable_;
  
  std::string sysUncertainty_;
};


/////////////////Start smooth function///////////////////////////////////


void FitHisto(TH1* htest){
  
  //the Append() works with string 
  //will keep the genuine name of the histo
  std::string histoNameInput = htest->GetName();
  //will then have append !!
  std::string histoName     = htest->GetName();
  std::string histoNamePDF  = htest->GetName();
  std::string histoNameEPS  = htest->GetName();
  std::string histoNamePNG  = htest->GetName();
  
  cout <<"hitoName" << histoNameInput << endl;
  
  double binPDF = htest->GetXaxis()->GetNbins();
  
  //General Variable      
  RooRealVar x("x","mvis",0,200);
  //RooRealVar x("x","mT",0,200);
  
  

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
  
  //Choose Fit depending on the histo name
  
  RooPlot* frame = x.frame(Title("Fit"));
  
  const char *histoToFitName;
  histoToFitName = histoNameInput.c_str();
  // convert CHAR to STR
  string histoToFitNameSTR (histoToFitName);
  
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
 
    hist = LandauConvGauss.createHistogram("SmoothHisto",x,Binning(binPDF)); //Bin, min, max
  }
  
  else if ( histoToFitNameSTR == histoToSmooth2STR){
    CristalBall1.fitTo(data);
    data.plotOn(frame);
    CristalBall1.plotOn(frame,LineColor(kRed)); 

    hist = CristalBall1.createHistogram("SmoothHisto",x,Binning(binPDF)); //Bin, min, max
  }
  
  else if (( histoToFitNameSTR == histoToSmooth3STR) ||
	   ( histoToFitNameSTR == histoToSmooth4STR )){
    CristalBall2.fitTo(data); 
    data.plotOn(frame);
    CristalBall2.plotOn(frame,LineColor(kRed)); 
 
    hist = CristalBall2.createHistogram("SmoothHisto",x,Binning(binPDF)); //Bin, min, max
  }
  
  else if (( histoToFitNameSTR == histoToSmooth5STR) ||
	   ( histoToFitNameSTR == histoToSmooth6STR )){
    CristalBall3.fitTo(data);
    data.plotOn(frame);
    CristalBall3.plotOn(frame,LineColor(kRed)); 
    
    hist = CristalBall3.createHistogram("SmoothHisto",x,Binning(binPDF)); //Bin, min, max
  }
  
  
  //Draw frame on canvas
  TCanvas *c = new TCanvas("","Fit",600,600) ;
  gPad->SetLeftMargin(0.15) ; 
  frame->GetYaxis()->SetTitleOffset(1.4) ;
  frame->Draw();
  
  //now the histoName variable has the string histo...+smooth!!
  std::string histoNameSTR = histoName.append("_smoothed.root");
  const char *histoNameCHAR;
  histoNameCHAR=histoNameSTR.c_str();
  
  
  //save canvas in pdf eps and png    
  TString fileNamePDF = histoNamePDF.append(".pdf");
  TString fileNameEPS = histoNameEPS.append(".eps"); 
  TString fileNamePNG = histoNamePNG.append(".png");
  c->Print(fileNamePDF);
  c->Print(fileNameEPS);
  c->Print(fileNamePNG);
  
  //Get the Pdf of the fit, fill it in a histogram and save
  TFile f(histoNameCHAR,"RECREATE"); //TFile works with char
  
  //set Bin histoPDF to be the same then the template histogram  
  // x.setBins( htest->Integral());

  hist->Write();
  
  f.Write();
  f.Close();
  cout <<"END Smoothing "<< histoName << endl; 
  
}

/////////////////End smooth function///////////////////////////////////



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

  //name of the output smooth histogram
  vstring vHistoFitted = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("HistoFitted"); 

  //name of the control Fit plot 
  vstring vControlFitPlot = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("ControlPlots"); 

  vstring tauIds = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("tauId"); 
  vstring tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");
  tauIdValues.push_back("D");


  vstring fitVariables = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("fitVariable");
  add_string_uniquely(fitVariables, "EventCounter"); // CV: for normalization purposes, always add 'EventCounter'
  add_string_uniquely(fitVariables, "diTauMt");      // CV: add 'diTauMt' in order to take QCD template for region 'D' from data



  vstring sysUncertainties = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("sysUncertainties");
  vstring sysUncertainties_expanded;
  sysUncertainties_expanded.push_back(key_central_value);
  for ( vstring::const_iterator sysUncertainty = sysUncertainties.begin();
	sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Up"));
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Down"));
  }
  
  vstring sysUncertainties_data;
  sysUncertainties_data.push_back(key_central_value);



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

  vstring processes;
  processes.push_back(std::string("ZplusJets"));
  processes.push_back(std::string("QCD"));
  processes.push_back(std::string("WplusJets"));
  processes.push_back(std::string("TTplusJets"));



  map< string, TH1* > htest;

  for ( vstring::const_iterator histogramName = vHistoToFit.begin();
       histogramName != vHistoToFit.end(); ++histogramName ) {

       //Get the Hitsto to be fit 
       htest[histogramName->data()] = (TH1F*)histogramInputFile->Get(histogramName->data());

   //int pos = vHistoToFit.at(*histogramName);

     //htest[] = (TH1*)histogramInputFile->Get(histogramName->data());



       //Fit the histo with the Good function
       //TH1F* htest = (TH1F*)histogramInputFile->Get(histogramName->data());


      FitHisto(htest[histogramName->data()] );

   } 



  for ( vstring::const_iterator tauId = tauIds.begin(); // for tauDiscrHPScombLooseDBcorr
	tauId != tauIds.end(); ++tauId ) {

    vstring regions;
    std::vector<regionEntryType*> regionEntries;

    for ( vstring::const_iterator tauIdValue = tauIdValues.begin(); // for passed, failed, D 
	  tauIdValue != tauIdValues.end(); ++tauIdValue ) {

      std::string regionTakeQCDtemplateFromData = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionTakeQCDtemplateFromData_").append(*tauIdValue));
      std::string regionWplusJetsSideband = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionWplusJetsSideband_").append(*tauIdValue));      
      std::string regionStoreQCDtemplate = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionStoreQCDtemplate_").append(*tauIdValue));
      
      if ( regionTakeQCDtemplateFromData == "" || regionWplusJetsSideband == "" || regionStoreQCDtemplate == "" ) continue;
      
      add_string_uniquely(regions, regionTakeQCDtemplateFromData);
      //add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("+"));
      //add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("-"));	
      add_string_uniquely(regions, regionWplusJetsSideband);
      //add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("+"));
      //add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("-"));

      for ( vstring::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin();
	      sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  // all tau chargee
	  regionEntryType* regionEntry = 
	    new regionEntryType(regionTakeQCDtemplateFromData, 
				regionWplusJetsSideband, 
				regionStoreQCDtemplate, 
				*tauId, *fitVariable, *sysUncertainty);
	  regionEntries.push_back(regionEntry);

	  // tau+ candidates only
	  //regionEntryType* regionEntry_plus = 
	  //  new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("+"), 
	  //			  std::string(regionWplusJetsSideband).append("+"), 
	  //			  std::string(regionStoreQCDtemplate).append("+"), 
	  //			  *tauId, *fitVariable, *sysUncertainty);
	  //regionEntries.push_back(regionEntry_plus);
	  //
	  // tau- candidates only
	  //regionEntryType* regionEntry_minus = 
	  //  new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("-"), 
	  //			  std::string(regionWplusJetsSideband).append("-"), 
	  //			  std::string(regionStoreQCDtemplate).append("-"), 
	  //			  *tauId, *fitVariable, *sysUncertainty);
	  //regionEntries.push_back(regionEntry_minus);
	}
      }
    }

    histogramMap4 histograms_mc; // key = (process, region, observable, central value/systematic uncertainty)
    valueMap3 numEvents_mc; // key = (process, region, central value/systematic uncertainty)
    for ( vstring::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      loadHistograms(histograms_mc[*process], histogramInputDirectory, 
		     *process, regions, *tauId, fitVariables, sysUncertainties_expanded);
      
      for ( vstring::const_iterator region = regions.begin();
	    region != regions.end(); ++region ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin();
	      sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  numEvents_mc[*process][*region][*sysUncertainty] = 
	    getIntegral(histograms_mc[*process][*region]["EventCounter"][*sysUncertainty], true, true);
	}
      }
    }
    
    histogramMap3 histograms_data; // key = (region, observable, key_central_value)
    loadHistograms(histograms_data, histogramInputDirectory, 
		   "Data", regions, *tauId, fitVariables, sysUncertainties_data);

    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {
      (*regionEntry)->makeQCDtemplate(histogramOutputDirectory, histograms_data, histograms_mc, numEvents_mc);
    }
  }

  delete histogramInputFile;

//--print time that it took macro to run
  std::cout << "finished executing makeTauIdEffQCDtemplate macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("makeTauIdEffQCDtemplate");

  return 0;
}
