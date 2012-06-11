
/** \executable makeTauIdEffQCDtemplate
 *
 * Obtain QCD template from control region in Data,
 * corrected for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.8 $
 *
 * $Id: makeTauIdEffQCDtemplate.cc,v 1.8 2012/05/25 08:17:27 veelken Exp $
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
#include "TCanvas.h"
//RooFit
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooCBShape.h"

#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"

#include "RooFitResult.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooFitResult.h"


using namespace RooFit;

typedef std::map<std::string, std::map<std::string, TH1*> > histogramMap;
typedef std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsMap;

struct regionEntryType
{
//Define the region with parameter pass throug the config 
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

//make template histo for each bkg
  void makeQCDtemplate(TFileDirectory& dir, 
  		       histogramMap3& histograms_data, histogramMap4& histograms_mc, valueMap3& numEvents_mc)

  {
    //std::cout << "<regionEntryType::makeWplusJetstemplate>" << std::endl;
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
    
    //TH1* templateWplusJetsQCDsideband = histograms_mc["WplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateWplusJetsQCDsideband = histograms_mc["WplusJets"][regionWplusJetsSideband_][fitVariable_][tauId_];


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
    //templateQCDsideband_obsQCD->Add(distributionDataQCDsideband, +1.);
    //templateQCDsideband_obsQCD->Add(templateZtautauQCDsideband, -1.);
    //templateQCDsideband_obsQCD->Add(templateWplusJetsQCDsideband, -1.*scaleFactorWplusJets);
    //templateQCDsideband_obsQCD->Add(templateTTplusJetsQCDsideband, -1.);

}

  std::string regionTakeQCDtemplateFromData_;
  std::string regionWplusJetsSideband_;
  std::string regionStoreQCDtemplate_;
 
  std::string tauId_;
  std::string fitVariable_;
  
  std::string sysUncertainty_;
 };


//void WplusJetsRooFitPDF(TH1 &htest){
void WplusJetsRooFitPDF(){

  cout << "start RooFit"<< endl;

  TFile* fh = new TFile("/afs/cern.ch/work/c/calpas/CMSSW/FWliteHisto/analyzeTauIdEffHistograms_all_corrQCDtemplates_2012May12v1_2.root");
  //WplusJets
  //TH1* htest = (TH1F*)fh->Get("WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake"); //OK fit LG1
  //TH1* htest = (TH1F*)fh->Get("WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake"); //OK fit CB1 
  //TH1* htest = (TH1F*)fh->Get("WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake"); //OK fit CB2
  //TH1* htest = (TH1F*)fh->Get("WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake");//OK fit CB2
  //QCD
  //TH1* htest = (TH1F*)fh->Get("QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all"); // OK fit CB3
  TH1* htest = (TH1F*)fh->Get("QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all"); //OK fit CB3
 
  double binPDF = htest->GetXaxis()->GetNbins();
  cout << "binPDF: " << binPDF << endl;   

 
  //General Variable      
  RooRealVar x("x","mvis",0,200);
  //RooRealVar x("x","mT",0,200);
  
  //Landau convoluted with Gauss function
  /*
  //OK LG1
  //Landau
  RooRealVar meanl("meanl","mean of Landau",80.,60,90);
  RooRealVar sigmal("sigmal","sigma of Landau",30,20,50);  //the level of the pick
  RooLandau landau("landau","landau",x,meanl,sigmal);
  //Gauss
  RooRealVar meang("meang","mean of Gaussian",1); //minimum of the tail of the pic?!
  RooRealVar sigmag("sigmag","sigma of Gaussian",1,0.1,10);
  RooGaussian gauss("gauss","gauss",x,meang,sigmag);
  //Landau convoluted Gauss
  RooFFTConvPdf LandauConvGauss("LandauConvGauss","LandauConvGauss",x,landau,gauss);
  */

  //Crystal Ball function
 
  /*
  // OK CB1
  RooRealVar cbmean("cbmean", "cbmean" , 90.0, 20, 180.0);
  RooRealVar cbsigma("cbsigma", "cbsigma" , 1, 1.0, 40.0); //pick wide
  RooRealVar cbsig("cbsig", "cbsignal", 800, 20, 80);
  RooRealVar n("n","", 0.2); // increase the level of the back tail
  RooRealVar alpha("alpha","", 1.3);
  RooCBShape CristalBall("CristalBall","CristalBall",x,cbmean,cbsigma,alpha,n);
  */

  /*
  // OK CB2
  RooRealVar cbmean("cbmean", "cbmean" , 50, 20, 180.0);
  RooRealVar cbsigma("cbsigma", "cbsigma" , 1, 1.0, 200.0); //pick wide
  RooRealVar cbsig("cbsig", "cbsignal", 60, 20, 80);
  RooRealVar n("n","", 1); // increase the level of the back tail
  RooRealVar alpha("alpha","", 1);
  RooCBShape CristalBall("CristalBall","CristalBall",x,cbmean,cbsigma,alpha,n);
  */

  
  //OK CB3
  RooRealVar cbmean("cbmean", "cbmean" , 1, 20, 180.0);
  RooRealVar cbsigma("cbsigma", "cbsigma" , 1, 1.0, 200.0); //pick wide
  RooRealVar cbsig("cbsig", "cbsignal", 60, 20, 80);
  RooRealVar n("n","", 1); // increase the level of the back tail
  RooRealVar alpha("alpha","", 1);
  RooCBShape CristalBall("CristalBall","CristalBall",x,cbmean,cbsigma,alpha,n);
  

  //Make hist as a abspdf to be fit
  RooDataHist data("","",x,htest, 1.0); 

  //Fit histo
  //LandauConvGauss.fitTo(data);
  CristalBall.fitTo(data); 
  
 
  //Draw histo and fit 
  //RooPlot* frame = x.frame(Title("WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake"));
  //RooPlot* frame = x.frame(Title("WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake"));
  //RooPlot* frame = x.frame(Title("WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake"));
  //RooPlot* frame = x.frame(Title("WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake"));
  //RooPlot* frame = x.frame(Title("QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all"));
  RooPlot* frame = x.frame(Title("QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all"));
  
  data.plotOn(frame);
  //LandauConvGauss.plotOn(frame,LineColor(kRed)); 
  CristalBall.plotOn(frame,LineColor(kRed));
  
  //Draw frame on canvas
  TCanvas *cWplusJets = new TCanvas("cWplusJets","WplusJets",600,600) ;
  gPad->SetLeftMargin(0.15) ; 
  frame->GetYaxis()->SetTitleOffset(1.4) ;
  frame->Draw();
  
  
  //save canvas in pdf eps and png    
  // TString fileNamePDF = "WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake.pdf";
  // TString fileNameEPS = "WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake.eps"; 
  // TString fileNamePNG = "WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake.png";
  // cWplusJets->Print(fileNamePDF);
  // cWplusJets->Print(fileNameEPS);
  // cWplusJets->Print(fileNamePNG);
  
  // TString fileNamePDF = "WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake.pdf";
  // TString fileNameEPS = "WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake.eps"; 
  // TString fileNamePNG = "WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake.png";
  // cWplusJets->Print(fileNamePDF);
  // cWplusJets->Print(fileNameEPS);
  // cWplusJets->Print(fileNamePNG);
  
  // TString fileNamePDF = "WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.pdf";
  // TString fileNameEPS = "WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.eps"; 
  // TString fileNamePNG = "WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.png";
  // cWplusJets->Print(fileNamePDF);
  // cWplusJets->Print(fileNameEPS);
  // cWplusJets->Print(fileNamePNG);
  
  // TString fileNamePDF = "WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.pdf";
  // TString fileNameEPS = "WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.eps"; 
  // TString fileNamePNG = "WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.png";
  // cWplusJets->Print(fileNamePDF);
  // cWplusJets->Print(fileNameEPS);
  // cWplusJets->Print(fileNamePNG);

  // TString fileNamePDF = "QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all.pdf";
  // TString fileNameEPS = "QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all.eps"; 
  // TString fileNamePNG = "QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all.png";
  // cWplusJets->Print(fileNamePDF);
  // cWplusJets->Print(fileNameEPS);
  // cWplusJets->Print(fileNamePNG);

  TString fileNamePDF = "QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all.pdf";
  TString fileNameEPS = "QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all.eps"; 
  TString fileNamePNG = "QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all.png";
  cWplusJets->Print(fileNamePDF);
  cWplusJets->Print(fileNameEPS);
  cWplusJets->Print(fileNamePNG);
  
  //Get the Pdf of the fit, fill it in a histogram and save
  //TFile f("WplusJets_C1f_diTauVisMass_tauDiscrHPScombLooseDBcorr_failed_JetToTauFake.root","RECREATE");
  //TFile f("WplusJets_D_diTauMt_tauDiscrHPScombLooseDBcorrAndMuonVeto_all_JetToTauFake.root","RECREATE");
  //TFile f("WplusJets_A_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.root","RECREATE");
  //TFile f("WplusJets_B_diTauMt_tauDiscrHPScombLooseDBcorr_all_JetToTauFake.root","RECREATE");
  //TFile f("QCD_A_diTauMt_tauDiscrHPScombLooseDBcorr_all.root","RECREATE");
  TFile f("QCD_B_diTauMt_tauDiscrHPScombLooseDBcorr_all.root","RECREATE");
  
  //set Bin histoPDF to be the same then the template histogram  
  // x.setBins( htest->Integral());
  //TH1* hist = LandauConvGauss.createHistogram("WplusJets",x,Binning(binPDF)); //Bin, min, max
  TH1* hist = CristalBall.createHistogram("WplusJets",x,Binning(binPDF));

  //LandauConvGauss.Write();
  CristalBall.Write();   
  
  f.Write();
  f.Close();
  cout <<"END RooFit" << endl; 

  }



int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<makeTauIdEffWplusJetsRoFitPDF>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeTauIdEffWplusJetsRooFitPDF");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("makeTauIdEffWplusJetsRooFitPDF") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process"); // Get the parameter of process via PSet

  edm::ParameterSet cfgMakeTauIdEffQCDtemplate = cfg.getParameter<edm::ParameterSet>("fitTauIdEff"); // Get the parameter for fitTauIdEff process

  vstring tauIds = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("tauId"); // TauId is part of fitTauId parameter
  vstring tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");
  tauIdValues.push_back("D");

  vstring fitVariables = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("fitVariable");
  //Betty: add_string_uniquely look for a string "diTauMt" for ex, and if it doesn't exist in the vector fitVariable, it will add it,
  //Betty: At this point is has "diTauMt" passed via the cfg, so now we add "EventConter" to the vector's list.
  //add_string_uniquely(fitVariables, "EventCounter"); // CV: for normalization purposes, always add 'EventCounter'
  //add_string_uniquely(fitVariables, "diTauMt");      // CV: add 'diTauMt' in order to take QCD template for region 'D' from data

  vstring JetToTauFakes = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("JetToTauFakeID");

  vstring sysUncertainties = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("sysUncertainties");
  //Betty: extend the key_central_value name (syUncertainties) with up and down
  vstring sysUncertainties_expanded;
  sysUncertainties_expanded.push_back(key_central_value);
  for ( vstring::const_iterator sysUncertainty = sysUncertainties.begin();
	sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Up")); // Here we make vector of sysUncertainty Up and Down
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Down"));
  }
  
  vstring sysUncertainties_data;
  sysUncertainties_data.push_back(key_central_value); // For data there is no sysUncertainty

   //Betty: Get the hitogram root file via the cfg
   fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauIdEffWplusJetsRooFitPDF") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgMakeTauIdEffQCDtemplate.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("makeTauIdEffWplusJetsRooFitPDF") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory histogramOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;

  vstring processes;
  processes.push_back(std::string("WplusJets"));
  //processes.push_back(std::string("QCD"));
  //processes.push_back(std::string("ZplusJets"));
  //processes.push_back(std::string("TTplusJets"));
  //cout << "process" << processes[0] <<endl;  //WplusJets processes

//loop over tauId = HPS working point
  for ( vstring::const_iterator tauId = tauIds.begin(); //tauId = only tauDiscrHPScombLooseDBcorr (cfg)
	tauId != tauIds.end(); ++tauId ) { 
    vstring regions; //vector of string which will contain the name of the region 
    std::vector<regionEntryType*> regionEntries; //Make a regionEntries vector of struct{regionEntryType}

    for ( vstring::const_iterator tauIdValue = tauIdValues.begin(); //tauIdValues = passed, failed, D
	  tauIdValue != tauIdValues.end(); ++tauIdValue ) {

    //Betty: Create region name for tauId and for each tauIdValue, so for tauDiscrHPScombLooseDBcorr we have passed, failed and D name create  
      //std::string regionTakeQCDtemplateFromData = 
      //cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionTakeQCDtemplateFromData_").append(*tauIdValue));
      std::string regionTakeQCDtemplateFromData = 
      cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("region_").append(*tauIdValue));

      //std::string regionWplusJetsSideband = 
      //cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionWplusJetsSideband_").append(*tauIdValue));
      std::string regionWplusJetsSideband = 
      cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("region_").append(*tauIdValue));

      //std::string regionStoreQCDtemplate = 
      //cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionStoreQCDtemplate_").append(*tauIdValue));
      std::string regionStoreQCDtemplate = 
      cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("region_").append(*tauIdValue));

      
      if ( regionTakeQCDtemplateFromData == "" || regionWplusJetsSideband == "" || regionStoreQCDtemplate == "" ) continue;
      
      //Betty: here we add the name of each region in an empty vector call region. This is done with the function add_string_uniquely: see tauIdEffAuxFunctions.h 
      add_string_uniquely(regions, regionTakeQCDtemplateFromData);
      //add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("+"));
      //add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("-"));	
      add_string_uniquely(regions, regionWplusJetsSideband);
      //add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("+"));
      //add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("-"));

      for ( vstring::const_iterator fitVariable = fitVariables.begin(); // fitVariables = only diTauVisMass or diTauMt see cfg
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin(); // uncertainty with extension Up Down define above
	      sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  // all tau chargee
	  regionEntryType* regionEntry = 
	    new regionEntryType(regionTakeQCDtemplateFromData, 
				regionWplusJetsSideband, 
				regionStoreQCDtemplate, 
				*tauId, *fitVariable, *sysUncertainty);
	  regionEntries.push_back(regionEntry); // fill a vector of struct regionEntries with struct variable defined before  

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
    for ( vstring::const_iterator process = processes.begin(); // loop over WJets, QCD, ZJets, TTJets
	  process != processes.end(); ++process ) {

    //std::cout << "<regionEntryType::makeWplusJetstemplate>" << std::endl;
    //std::cout << " regionWplusJetsSideband = " << regionWplusJetsSideband_ << std::endl;

    //cout << " process = " << *process << endl;
    //std::cout << " sysUncertainty = " << sysUncertainty_ << std::endl;


      //loadHistograms(histograms_mc[*process], histogramInputDirectory, 
		     //*process, regions, *tauId, fitVariables, sysUncertainties_expanded);
      loadHistograms(histograms_mc[*process], histogramInputDirectory, 
		     *process, regions, *tauId, fitVariables, JetToTauFakes);  // Here we load the histogram with the specique name. I create a JetToTauFakes tag in the cfg. It is then pass to the loadHistograms function (in place of the syst variable), and there the histograms_mc is the histogram. 

      for ( vstring::const_iterator region = regions.begin();
	    region != regions.end(); ++region ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin();
	   sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  //numEvents_mc[*process][*region][*sysUncertainty] = 
	    //getIntegral(histograms_mc[*process][*region]["EventCounter"][*sysUncertainty], true, true);
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
  std::cout << "finished executing makeTauIdEffWplusJetsRooFitPDF macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("makeTauIdEffWplusJetsRooFitPDF");

  return 0;
}
