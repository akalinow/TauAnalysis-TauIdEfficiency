//Roofit + FWLite macro to fit the plots produced by HistoPlotter

// STL & system
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

using namespace std;

//Root include
#include "TH1F.h"
#include "TMath.h"
#include "TDirectory.h"
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
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

using namespace TMath;

//Roofit Include
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooCmdArg.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooProduct.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "RooFitResult.h"

//FWLite Include
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

using namespace edm;

typedef vector<string> vstring;
typedef vector<ParameterSet> vpset;

inline double getIntegral(const TH1* histogram, bool inclUnderflowBin=false, bool inclOverflowBin=false)
{
  //--------------------------------------------------------------------------------
  // Compute integral of histogram including/excluding underflow and overflow bins
  //--------------------------------------------------------------------------------
  int firstBin = ( inclUnderflowBin ) ?                            0 :                      1;
  int lastBin  = ( inclOverflowBin  ) ? (histogram->GetNbinsX() + 1) : histogram->GetNbinsX();
  double integral = 0;
  for ( int iBin = firstBin; iBin <= lastBin; ++iBin )
    integral += histogram->GetBinContent(iBin);

  return integral;
}



map<string, TH1*> loadHistograms( TFile* inputFile, const ParameterSet& channel, const vstring& regions,
				 const vstring& tauIds, const vpset& fitVariables, const vstring& sysShifts)
{
//--------------------------------------------------------------------------------
// Load template histograms/distributions observed in data from ROOT file
// Dir Structure: channel/tauid.name/region/plotVariable.name/plotVariable.name+"Distribution"(+"_"+sysShift)
// map = (channel/tauid.name/region/observable_sysShift)
//--------------------------------------------------------------------------------

  map<string, TH1*> retVal;
  const string & channelName = channel.getParameter<string>("name");
  const string & channelLabel = channel.getParameter<string>("label");
  const int & channelColor = channel.getParameter<int>("color");

  TH1* tmpHisto = dynamic_cast<TH1*>(inputFile->Get( (channelName+"/AnalyzedEvents").c_str() ) );
  if ( !tmpHisto ) {
    cout << "Error in <loadHistograms>: failed to load histogram = " << (channelName+"/AnalyzedEvents") 
	 << " from file = " << inputFile->GetName() << " --> aborting !!";
    assert(0);
  }
  retVal[(channelName+"/AnalyzedEvents").c_str()] = tmpHisto;

  for ( vstring::const_iterator tauId = tauIds.begin(); tauId != tauIds.end(); ++tauId ) {
    for ( vstring::const_iterator region = regions.begin(); region != regions.end(); ++region ) {
      for ( vpset::const_iterator fitVariable = fitVariables.begin(); fitVariable != fitVariables.end(); ++fitVariable ) {
	const string& varname = fitVariable->getParameter<string>("name");
	const int& rebinFactor = fitVariable->getParameter<int>("rebinFactor");
	const string& xAxisTitle = fitVariable->getParameter<string>("xAxisTitle");
	for ( vstring::const_iterator sysShift = sysShifts.begin(); sysShift != sysShifts.end(); ++sysShift ) {
	  string histoName = channelName+"/"+*tauId+"/"+*region+"/"+varname+"/"+varname+"Distribution";
	  if ( *sysShift != "CENTRAL_VALUE" ) histoName.append("_").append(*sysShift);
	  tmpHisto = dynamic_cast<TH1*>(inputFile->Get( histoName.c_str() ) );
	  tmpHisto->GetXaxis()->SetTitle(xAxisTitle.c_str());
	  tmpHisto->SetTitle(channelLabel.c_str());
	  tmpHisto->SetLineColor(channelColor);
	  if ( !tmpHisto ) {
	    cout << "Error in <loadHistograms>: failed to load histogram = " << histoName 
		 << " from file = " << inputFile->GetName() << " --> aborting !!";
	    assert(0);
	  }
	  tmpHisto->Rebin(rebinFactor);
	  retVal[(channelName+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()] = tmpHisto;
	} //SysShift loop
      } //tauId loop
    } //region loop
  } //var loop

  return retVal;
}

map<string, TH1*> sumHistograms( vstring channelNames,map<string, TH1*> allDistributions,string newChannelName, const vstring& regions, const vstring& tauIds, const vpset& fitVariables, const vstring& sysShifts)
{  
  /*----------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sums up all the template histograms to create a fake data sample to check if the fitting procedure manages to converge to the simulated values
   ----------------------------------------------------------------------------------------------------------------------------------------------------------------*/
  map<string, TH1*> retVal;
  for ( vstring::const_iterator tauId = tauIds.begin(); tauId != tauIds.end(); ++tauId ) {
    for ( vstring::const_iterator region = regions.begin(); region != regions.end(); ++region ) {
      for ( vpset::const_iterator fitVariable = fitVariables.begin(); fitVariable != fitVariables.end(); ++fitVariable ) {
	const string& varname = fitVariable->getParameter<string>("name");
	for ( vstring::const_iterator sysShift = sysShifts.begin(); sysShift != sysShifts.end(); ++sysShift ) {
	  string histoName = newChannelName+"/"+*tauId+"/"+*region+"/"+varname+"/"+varname+"Distribution";
	  if ( *sysShift != "CENTRAL_VALUE" ) histoName.append("_").append(*sysShift);
	  bool first = true;
	  for(vstring::iterator channel = channelNames.begin(); channel != channelNames.end(); channel++){
	    if(first){
	      retVal[(newChannelName+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()] = (TH1*) allDistributions[(*channel+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()]->Clone(histoName.c_str());
	      first = false;
	      retVal[(newChannelName+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()]->SetLineColor(0);
	    }
	    else{
	      retVal[(newChannelName+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()]->Add(allDistributions[(*channel+"/"+*tauId+"/"+*region+"/"+varname+"_"+*sysShift).c_str()]);
	    }
	  } //channels Loop
	} //SysShift loop
      } //tauId loop
    } //region loop
  } //var loop
  return retVal;
}

void DrawHistograms(TDirectory* dir, map<string,map<string, TH1* > > channelByChannelTemplateDistro, map<string, TH1* > toFitDistro,
		    const string& toFitChannelName, const string& region, const string& tauId, const string& varname, const string& sysShift)
{
  /*------------------------------------------------------------------------------------------------------------------------------------------
    Takes the key and plots the (fake) data points over a stacked histogram with all the components separated, writes it in the directory given
   ------------------------------------------------------------------------------------------------------------------------------------------*/
  dir->cd();
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

  string key = "/"+tauId+"/"+region+"/"+varname+"_"+sysShift;
  TLegend legend(0.64, 0.64, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);

  TH1* datahisto = toFitDistro[toFitChannelName+key];
  datahisto->SetLineColor(1);
  datahisto->SetMarkerColor(1);
  datahisto->SetMarkerStyle(20);
  
  datahisto->SetStats(false);
  datahisto->GetXaxis()->SetTitleOffset(1.15);
  legend.AddEntry(datahisto,"Data","p");
  
  THStack smSum(("smSum_"+sysShift).c_str(), ("smSum_"+sysShift).c_str());
  for( map<string,map<string, TH1* > >::iterator channelMap=channelByChannelTemplateDistro.begin(); channelMap!=channelByChannelTemplateDistro.end(); channelMap++){
    TH1* templateH = channelMap->second[channelMap->first+key];
    templateH->SetStats(false);
    smSum.Add(templateH);
    legend.AddEntry(templateH,templateH->GetTitle(),"f");
  }//template loop
  
  datahisto->GetYaxis()->SetRangeUser(0,1.4*TMath::Max(smSum.GetMaximum(), datahisto->GetMaximum()));
  datahisto->Draw("ep1");
  smSum.Draw("hist same");
  
  TPaveText cmsPreliminaryLabel(0.140, 0.860, 0.46, 0.900, "NDC");
  cmsPreliminaryLabel.AddText("CMS Preliminary 2011");
  cmsPreliminaryLabel.SetTextAlign(13);
  cmsPreliminaryLabel.SetTextSize(0.045);
  cmsPreliminaryLabel.SetFillStyle(0);
  cmsPreliminaryLabel.SetBorderSize(0);
  cmsPreliminaryLabel.Draw();
  
  /*  TPaveText cmsLuminosityLabel(0.145, 0.8075, 0.460, 0.8475, "NDC");
  cmsLuminosityLabel.AddText("#sqrt{s} = 7 TeV, L = 36 pb^{-1}");
  cmsLuminosityLabel.SetTextAlign(13);
  cmsLuminosityLabel.SetTextSize(0.045);
  cmsLuminosityLabel.SetFillStyle(0);
  cmsLuminosityLabel.SetBorderSize(0);
  cmsLuminosityLabel.Draw();*/

  canvas->SetName(sysShift.c_str());
  canvas->Write();  

  delete canvas;
}

void fitUsingRooFit(TDirectory* tauDir,
		    map<string, TH1*> & distributionsData,
		    string & dataChannelName,
		    map<string,map<string, TH1* > > & channelByChannelTemplateDistro,
		    map<string, double> & numEventsAll,
		    map<string, double>& fittedEvents,
		    const string& tauId, const ParameterSet& fitVariable,
		    const vpset& regionsPset,
		    const ParameterSet & fitParameters,
		    double& effValue, double& effError,
		    const string& sysShift,
		    bool isClosureTest = false) 
{
  //-------------------------------------------------------------------------------
  // Make RooRealVar object representing normalization for MC process 
  // passed as function argument in a given region.
  //
  // The regions are defined as follows:
  //
  //   /---------------\ /---------------\
  //   |               | |               |
  //   |               | |               |
  //   |       A       | |       B       | 0.1 * muonPt < muonIso < 0.3 * muonPt
  //   |               | |               |
  //   |               | |               |
  //   \---------------/ \---------------/
  //
  //   /---------------\ /---------------\
  //   |               | |               |
  //   |               | |               |
  //   |       C       | |       D       |                muonIso < 0.1 * muonPt
  //   |               | |               |
  //   |               | |               |
  //   \---------------/ \---------------/
  //       
  //           OS                SS
  //
  // C1p : Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV && tau id. passed
  // C1f : Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV && tau id. failed
  // C2p : Mt > 40 GeV || (Pzeta - 1.5 PzetaVis) < -20 GeV && tau id. passed
  // C2f : Mt > 40 GeV || (Pzeta - 1.5 PzetaVis) < -20 GeV && tau id. failed
  //
  //-------------------------------------------------------------------------------

  const string& varname = fitVariable.getParameter<string>("name");
  TDirectory* varDir = tauDir->mkdir(varname.c_str());
  cout << "<fitUsingRooFit>:" << endl;
  cout << " performing Fit of variable = " << varname << " for Tau id. = " << tauId << endl;
  cout << " isClosureTest = " << isClosureTest << endl;

  double fitMinABC2D = fitVariable.getParameter<double>("minh");
  double fitMaxABC2D = fitVariable.getParameter<double>("maxh");
  RooRealVar* fitVarABC2D = new RooRealVar("fitVarABC2D", "fitVarABC2D", fitMinABC2D, fitMaxABC2D);
  RooRealVar* fitVarC1 = new RooRealVar("fitVarC1", "fitVarC1", fitMinABC2D, fitMaxABC2D);
  
  string datakey = dataChannelName+"/"+tauId+"/ABCD/"+varname+"_"+sysShift;
  double numEventsDataABCD = -1;
  if(distributionsData.find(datakey.c_str()) != distributionsData.end() ) numEventsDataABCD = distributionsData[datakey]->Integral(); 
  else
    for(map<string, TH1*>::iterator entry = distributionsData.begin(); entry != distributionsData.end(); entry++)
      if(entry->first.find("CENTRAL_VALUE") != string::npos)
	numEventsDataABCD += entry->second->Integral(); //sums up all the components

  cout << "numEventsDataABCD = " << numEventsDataABCD << endl;

  map<string, RooRealVar*> pDiTauCharge_OS_SS;          // key = process
  map<string, RooRealVar*> pMuonIso_tight_loose;        // key = process
  map<string, RooRealVar*> pDiTauKine_Sig_Bgr;          // key = process
  map<string, RooRealVar*> pTauId_passed_failed;        // key = process

  map< string,map<string, RooAbsReal*> > normVals;                    // key = region,process
  map<string, map<string, RooHistPdf*> > pdf; // key = (region ,process/"sum" )
  map< string,TObjArray> pdfArray; //key = region
  map< string,TObjArray> fitPars; //key = region

  for( map < string,map<string, TH1* > > ::iterator entry = channelByChannelTemplateDistro.begin(); entry != channelByChannelTemplateDistro.end(); entry++){

    //SETTING VARIABLES AND THEUR STARTING VALUES
    string prekey = entry->first+"/"+tauId+"/";
    string postkey = "/"+varname+"_"+sysShift;
    double pDiTauCharge_OS_SS0 = 0.5;
    double pMuonIso_tight_loose0 = 0.5;
    double pDiTauKine_Sig_Bgr0 = 0.5;
    double pTauId_passed_failed0 = 0.5;
    if(fitParameters.getParameter<bool>("calculateStartingOSSSratio") ){
      if(fittedEvents.find(prekey+"A"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"ABCD"+postkey) != fittedEvents.end())
	pDiTauCharge_OS_SS0 = (fittedEvents[prekey+"A"+postkey] + fittedEvents[prekey+"C"+postkey])/fittedEvents[prekey+"ABCD"+postkey];
      else
	cout<<"Error! calculateStartingOSSSratio set to true implies region A,C,ABCD GIVEN! Skipping...";
    }
    if(fitParameters.getParameter<bool>("calculateStartingMuonIsoRatio") ){
      if(fittedEvents.find(prekey+"D"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"ABCD"+postkey) != fittedEvents.end())
	pMuonIso_tight_loose0        = (fittedEvents[prekey+"C"+postkey] + fittedEvents[prekey+"D"+postkey])/fittedEvents[prekey+"ABCD"+postkey];
      else
	cout<<"Error! calculateStartingOSSSratio set to true implies region D,C,ABCD GIVEN! Skipping...";
    }
    if(fitParameters.getParameter<bool>("calculateStartingTauKineBkSigRatio") ){
      if(fittedEvents.find(prekey+"C1"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() )
	pDiTauKine_Sig_Bgr0          = fittedEvents[prekey+"C1"+postkey]/fittedEvents[prekey+"C"+postkey];
      else
	cout<<"Error! calculateStartingTauKineBkSigRatio set to true implies region C,C1 GIVEN! Skipping...";
    }
    if(fitParameters.getParameter<bool>("calculateStartingTauIdPasFailRatio") ){
      if(fittedEvents.find(prekey+"C1p"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C2p"+postkey) != fittedEvents.end())
	pTauId_passed_failed0        = ( fitParameters.getParameter<bool>("fitTauIdEffC2") ) ?  (fittedEvents[prekey+"C1p"+postkey] + fittedEvents[prekey+"C2p"+postkey])/fittedEvents[prekey+"C"+postkey] : fittedEvents[prekey+"C1p"+postkey]/fittedEvents[prekey+"C1"+postkey];
      else
	cout<<"Error! calculateStartingOSSSratio set to true implies region C1p,C,C2p GIVEN! Skipping...";
    }
    string nameDiTauCharge_OS_SS   = string("pDiTauCharge_OS_SS").append("_").append(entry->first);
    pDiTauCharge_OS_SS[entry->first]        = new RooRealVar(nameDiTauCharge_OS_SS.data(), nameDiTauCharge_OS_SS.data(), pDiTauCharge_OS_SS0, 0., 1.);
    string nameMuonIso_tight_loose = string("pMuonIso_tight_loose").append("_").append(entry->first);
    pMuonIso_tight_loose[entry->first]      = new RooRealVar(nameMuonIso_tight_loose.data(), nameDiTauCharge_OS_SS.data(), pMuonIso_tight_loose0, 0., 1.);
    string nameDiTauKine_Sig_Bgr   = string("pDiTauKine_Sig_Bgr").append("_").append(entry->first);
     pDiTauKine_Sig_Bgr[entry->first]        = new RooRealVar(nameDiTauKine_Sig_Bgr.data(), nameDiTauKine_Sig_Bgr.data(), pDiTauKine_Sig_Bgr0, 0., 1.);
    string nameTauId_passed_failed = string("pTauId_passed_failed").append("_").append(entry->first);
    pTauId_passed_failed[entry->first]      = new RooRealVar(nameTauId_passed_failed.data(), nameTauId_passed_failed.data(), pTauId_passed_failed0, 0., 1.);

    double normABCD0         = (numEventsAll.find(prekey+"ABCD"+postkey) != numEventsAll.end()) ? numEventsAll[prekey+"ABCD"+postkey] : numEventsDataABCD;
    string nameNormABCD = string("normABCD").append("_").append(entry->first);
    normVals["ABCD"][entry->first]       = new RooRealVar(nameNormABCD.data(), nameNormABCD.data(), normABCD0, 0., numEventsDataABCD); //overall 
    
    for ( vpset::const_iterator region = regionsPset.begin(); region != regionsPset.end(); ++region ) {
      //Key : (channel/tauid.name/region/observable_sysShift)
      const string & regionName = region->getParameter<string>("name");
      if( !(region->getParameter<bool>("toFit")) || regionName == "ABCD" ) continue;
      string key = prekey+regionName+postkey;

      RooRealVar* fitVar = (regionName.find("C1") != string::npos) ? fitVarC1 : fitVarABC2D;
      string templateDataHistName = entry->first+"_"+tauId+"_"+regionName+"_"+varname+"_"+sysShift + "_dataHist";
      RooDataHist* templateDataHist = new RooDataHist(templateDataHistName.data(), templateDataHistName.data(), *fitVar, entry->second[key]);
      string templatePdfName = entry->first+"_"+tauId+"_"+regionName+"_"+varname+"_"+sysShift + "_histPdf";
      pdf[regionName][entry->first] = new RooHistPdf(templatePdfName.data(), templatePdfName.data(), *fitVar, *templateDataHist);
      if(pdfArray.find(regionName) == pdfArray.end())
	pdfArray[regionName] = TObjArray();

      pdfArray[regionName].Add(pdf[regionName][entry->first]);

      TObjArray arguments; 
      string formula = "";
      
      if        ( regionName.find("A") != string::npos ) {
	formula.append(pDiTauCharge_OS_SS[entry->first]->GetName());
	formula.append("*(1-").append(pMuonIso_tight_loose[entry->first]->GetName()).append(")");
      } else if ( regionName.find("B") != string::npos ) {
	formula.append("(1-").append(pDiTauCharge_OS_SS[entry->first]->GetName()).append(")");	  
	formula.append("*(1-").append(pMuonIso_tight_loose[entry->first]->GetName()).append(")");
      } else if ( regionName.find("C") != string::npos ) {
	formula.append(pDiTauCharge_OS_SS[entry->first]->GetName());
	formula.append("*").append(pMuonIso_tight_loose[entry->first]->GetName());
      } else if ( regionName.find("D") != string::npos ) {
	formula.append("(1-").append(pDiTauCharge_OS_SS[entry->first]->GetName()).append(")");	  
	formula.append("*").append(pMuonIso_tight_loose[entry->first]->GetName());
      } else {
	cout << "Error in <makeRooFormulaVar>: undefined region = " << regionName << " --> skipping !!" << endl;
	continue;
      }
      arguments.Add(pDiTauCharge_OS_SS[entry->first]);
      arguments.Add(pMuonIso_tight_loose[entry->first]);
      if      ( regionName.find("1") != string::npos ) {
	formula.append("*").append(pDiTauKine_Sig_Bgr[entry->first]->GetName());
	arguments.Add(pDiTauKine_Sig_Bgr[entry->first]);
      }
      else if ( regionName.find("2") != string::npos ) {
	formula.append("*(1-").append(pDiTauKine_Sig_Bgr[entry->first]->GetName()).append(")");
	arguments.Add(pDiTauKine_Sig_Bgr[entry->first]);
      }
	
      if      ( regionName.find("p") != string::npos ) {
	formula.append("*").append(pTauId_passed_failed[entry->first]->GetName());
	arguments.Add(pTauId_passed_failed[entry->first]);
      }
      else if ( regionName.find("f") != string::npos ) {
	formula.append("*(1-").append(pTauId_passed_failed[entry->first]->GetName()).append(")");
	arguments.Add(pTauId_passed_failed[entry->first]);
      }
      
      string fittedFractionName = string( normVals["ABCD"][entry->first]->GetName()).append("_").append(regionName).append("_fittedFraction");
      RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), fittedEvents[key]/numEventsAll[key]);
      formula += "*"+fittedFractionName;
      arguments.Add(fittedFraction);
      
      formula.append("*").append(normVals["ABCD"][entry->first]->GetName());
      arguments.Add(normVals["ABCD"][entry->first]);
      cout << " formula = " << formula << endl;
      
      normVals[regionName][entry->first] = new RooFormulaVar(("norm_"+regionName+"_"+entry->first).c_str(), ("norm_"+regionName+"_"+entry->first).c_str(), formula.data(), RooArgSet(arguments));
    
      if(fitPars.find(regionName) == fitPars.end()) fitPars[regionName] = TObjArray();
      fitPars[regionName].Add(normVals[regionName][entry->first]);
    }//region
  }//channel template loop

  //Loads roo categories and pdfs
  map<string,RooAddPdf*> pdfSum;
  // CV: due to limitation in RooFit
  //    (cf. http://root.cern.ch/phpBB3/viewtopic.php?f=15&t=9518)
  //     need to construct log-likelihood functions separately for regions { A, B, D } and { C1p, C1f }
  
  //--- build data & model objects for fitting regions A, B, C2p, C2f, D
  RooCategory* fitCategoriesABC2D = new RooCategory("categoriesABC2D", "categoriesABC2D");
  RooCategory* fitCategoriesC1 = new RooCategory("categoriesC1", "categoriesC1");
  for(map<string,TObjArray>::const_iterator entry = fitPars.begin(); entry != fitPars.end(); entry++){
    if(pdfArray.find(entry->first) != pdfArray.end()){
      pdfSum[entry->first] = new RooAddPdf( ("pdfSum"+entry->first).c_str(),("pdfSum"+entry->first).c_str(),RooArgList(pdfArray[entry->first]),RooArgList(entry->second));
      if(entry->first.find("C1") != string::npos)
	fitCategoriesC1->defineType( entry->first.c_str());
      else
	fitCategoriesABC2D->defineType( entry->first.c_str());
    }
    else{
      cout<<"ERROR: region " << entry->first <<" has no entry in template pdf (pdfArray)"<<endl;
      exit(0);
    }
  }

  //Loads Simultaneous fit and histo map
  RooSimultaneous* pdfSimultaneousFitABC2D = new RooSimultaneous("pdfSimultaneousFitABC2D", "pdfSimultaneousFitABC2D", *fitCategoriesABC2D);
  RooSimultaneous* pdfSimultaneousFitC1 = new RooSimultaneous("pdfSimultaneousFitC1", "pdfSimultaneousFitC1", *fitCategoriesC1);
  map<string, TH1*> histogramDataMapABC2D;
  map<string, TH1*> histogramDataMapC1;
  for(map<string,TObjArray>::const_iterator entry = fitPars.begin(); entry != fitPars.end(); entry++){
    if(entry->first.find("C1") != string::npos){
      pdfSimultaneousFitC1->addPdf(*pdfSum[entry->first],entry->first.c_str());     
      histogramDataMapC1[entry->first] = distributionsData[dataChannelName+"/"+tauId+"/"+entry->first+"/"+varname+"_"+sysShift];
    }
    else{
      pdfSimultaneousFitABC2D->addPdf(*pdfSum[entry->first],entry->first.c_str());
      histogramDataMapABC2D[entry->first] = distributionsData[dataChannelName+"/"+tauId+"/"+entry->first+"/"+varname+"_"+sysShift];
    }
  }
  
  //creates dataset
  RooDataHist* dataABC2D = new RooDataHist("dataABC2D", "data", *fitVarABC2D, *fitCategoriesABC2D, histogramDataMapABC2D);
  RooDataHist* dataC1 = new RooDataHist("dataC1", "dataC1", *fitVarC1, *fitCategoriesC1, histogramDataMapC1);

  //Fit Variables  pDiTauCharge_OS_SS   pMuonIso_tight_loose   pDiTauKine_Sig_Bgr   pTauId_passed_failed
  //Sets the variable for some channels to fixed
  const vstring & fixpDiTauCharge_OS_SS = fitParameters.getParameter<vstring>("fixPDiTauCharge_OS_SS");
  for(vstring::const_iterator channel = fixpDiTauCharge_OS_SS.begin(); channel != fixpDiTauCharge_OS_SS.end(); channel++){
    if(pDiTauCharge_OS_SS.find(*channel) != pDiTauCharge_OS_SS.end())
      pDiTauCharge_OS_SS[*channel]->setConstant(true);
    else
      cout << "Cannot set pDiTauCharge_OS_SS["<<*channel<<"] to constant!, Channel not found. Skipping.."<<endl;
  }
  const vstring & fixpMuonIso_tight_loose = fitParameters.getParameter<vstring>("fixpMuonIso_tight_loose");
  for(vstring::const_iterator channel = fixpMuonIso_tight_loose.begin(); channel != fixpMuonIso_tight_loose.end(); channel++){
    if(pMuonIso_tight_loose.find(*channel) != pMuonIso_tight_loose.end())
      pMuonIso_tight_loose[*channel]->setConstant(true);
    else
      cout << "Cannot set pMuonIso_tight_loose["<<*channel<<"] to constant!, Channel not found. Skipping.."<<endl;
  }
  const vstring & fixpDiTauKine_Sig_Bgr = fitParameters.getParameter<vstring>("fixpDiTauKine_Sig_Bgr");
  for(vstring::const_iterator channel = fixpDiTauKine_Sig_Bgr.begin(); channel != fixpDiTauKine_Sig_Bgr.end(); channel++){
    if(pDiTauKine_Sig_Bgr.find(*channel) != pDiTauKine_Sig_Bgr.end())
      pDiTauKine_Sig_Bgr[*channel]->setConstant(true);
    else
      cout << "Cannot set pDiTauKine_Sig_Bgr["<<*channel<<"] to constant!, Channel not found. Skipping.."<<endl;
  }
  const vstring & fixpTauId_passed_failed = fitParameters.getParameter<vstring>("fixpTauId_passed_failed");
  for(vstring::const_iterator channel = fixpTauId_passed_failed.begin(); channel != fixpTauId_passed_failed.end(); channel++){
    if(pTauId_passed_failed.find(*channel) != pTauId_passed_failed.end())
      pTauId_passed_failed[*channel]->setConstant(true);
    else
      cout << "Cannot set pTauId_passed_failed["<<*channel<<"] to constant!, Channel not found. Skipping.."<<endl;
  }

  //--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  //Puts some constraints to the normalization values 
  TObjArray fitConstraintsC1;
  const vpset & fitConstraintsC1Set = fitParameters.getParameter<vpset>("fitConstraintsC1"); //{channel name, sigma} ->gaussian centered in the predicted value and sigma = nsigma*predicted value
  for(vpset::const_iterator constraint = fitConstraintsC1Set.begin(); constraint != fitConstraintsC1Set.end(); constraint++){
    const string & channel = constraint->getParameter<string>("channel");
    if(normVals["ABCD"].find(channel) != normVals["ABCD"].end()){
      const double & nsigma = constraint->getParameter<double>("nsigma");
      string pValueName = string(normVals["ABCD"][channel]->GetName()).append("_constValue");
      RooConstVar* pValue = new RooConstVar(pValueName.data(), pValueName.data(), normVals["ABCD"][channel]->getVal());
      string pErrorName = string(normVals["ABCD"][channel]->GetName()).append("_constError");
      RooConstVar* pError = new RooConstVar(pErrorName.data(), pErrorName.data(), nsigma*normVals["ABCD"][channel]->getVal());  
      string constraintName = string(normVals["ABCD"][channel]->GetName()).append("_constraint");
      RooGaussian* constraint = new RooGaussian(constraintName.data(), constraintName.data(), *normVals["ABCD"][channel], *pValue, *pError);
    
      fitConstraintsC1.Add(constraint);
    }
    else
      cout << "Cannot find normABCD["<<channel<<"] . Skipping.."<<endl;
  }

  RooLinkedList fitOptionsC1;
  fitOptionsC1.Add(new RooCmdArg(RooFit::Extended()));
  fitOptionsC1.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsC1))));

  //Not used yet
  TObjArray fitConstraintsABC2D;
  RooLinkedList fitOptionsABC2D;
  fitOptionsABC2D.Add(new RooCmdArg(RooFit::Extended()));
  
  RooAbsReal* nllABC2D = pdfSimultaneousFitABC2D->createNLL(*dataABC2D, fitOptionsABC2D); 
  RooAbsReal* nllC1 = pdfSimultaneousFitC1->createNLL(*dataC1, fitOptionsC1); 
  RooAddition nll("nll", "nll", RooArgSet(*nllABC2D, *nllC1)); 
  RooMinuit minuit(nll); 
  //RooMinuit minuit(*nllC1);
  minuit.setErrorLevel(1);
  minuit.setNoWarn();
  minuit.setPrintEvalErrors(1);
  minuit.setPrintLevel(0);
  //minuit.setWarnLevel(1);
  minuit.migrad(); 
  minuit.hesse(); 

  //--- unpack covariance matrix of fit parameters
  string fitResultName = string("fitResult").append("_").append(tauId);
  RooFitResult*	fitResult = minuit.save(fitResultName.data(), fitResultName.data());
   
  cout << tauId << ":";
  if ( fitResult->status() == 0 ) cout << " fit converged."          << endl; 
  else                            cout << " fit failed to converge." << endl;

  const RooArgList& fitParameter = fitResult->floatParsFinal();

  int numFitParameter = fitParameter.getSize();
  
  TMatrixD cov(numFitParameter, numFitParameter);
  for ( int iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
    const RooAbsArg* paramI_arg = fitParameter.at(iParameter);
    const RooRealVar* paramI = dynamic_cast<const RooRealVar*>(paramI_arg);    
    double sigmaI = paramI->getError();

    cout << " parameter #" << iParameter << ": " << paramI_arg->GetName() 
	      << " = " << paramI->getVal() << " +/- " << paramI->getError() << endl;
    
    for ( int jParameter = 0; jParameter < numFitParameter; ++jParameter ) {
      const RooAbsArg* paramJ_arg = fitParameter.at(jParameter);
      const RooRealVar* paramJ = dynamic_cast<const RooRealVar*>(paramJ_arg);
      double sigmaJ = paramJ->getError();

      double corrIJ = fitResult->correlation(*paramI_arg, *paramJ_arg);

      cov(iParameter, jParameter) = sigmaI*sigmaJ*corrIJ;
    }
  }

  cov.Print();

  cout << endl;

  map<string,map<string, TH1* > > channelByChannelTemplateDistroNormFit; //key = channel key
  map< string, map<string, double> > normFactors;    // key = region process

  for( map < string,map<string, TH1* > > ::iterator entry = channelByChannelTemplateDistro.begin(); entry != channelByChannelTemplateDistro.end(); entry++){
    string prekey = entry->first+"/"+tauId+"/";
    string postkey = "/"+varname+"_"+sysShift;

    cout << " " << (entry->first) << ":" << endl;
    cout << "  normalization = " << normVals["ABCD"][entry->first]->getVal()
	 << " +/- " << dynamic_cast<RooRealVar*>(normVals["ABCD"][entry->first])->getError()
	 << " (MC exp. = " << numEventsAll[prekey+"ABCD"+postkey] << ")" << endl;

    cout << "  pDiTauCharge_OS_SS = " << pDiTauCharge_OS_SS[entry->first]->getVal() 
	      << " +/- " << pDiTauCharge_OS_SS[entry->first]->getError() ;
    if(fittedEvents.find(prekey+"A"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"ABCD"+postkey) != fittedEvents.end())
      cout << " (MC exp. = " << (fittedEvents[prekey+"A"+postkey] + fittedEvents[prekey+"C"+postkey])/fittedEvents[prekey+"ABCD"+postkey] << ")" << endl;

    cout << "  pMuonIso_tight_loose = " << pMuonIso_tight_loose[entry->first]->getVal() 
	      << " +/- " << pMuonIso_tight_loose[entry->first]->getError() ;
    if(fittedEvents.find(prekey+"D"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"ABCD"+postkey) != fittedEvents.end())
      cout << " (MC exp. = " << (fittedEvents[prekey+"C"+postkey] + fittedEvents[prekey+"D"+postkey])/fittedEvents[prekey+"ABCD"+postkey] << ")" << endl;

    cout << "  pDiTauKine_Sig_Bgr = " << pDiTauKine_Sig_Bgr[entry->first]->getVal() 
	      << " +/- " << pDiTauKine_Sig_Bgr[entry->first]->getError() ;
    if(fittedEvents.find(prekey+"C1"+postkey) != fittedEvents.end() && fittedEvents.find(prekey+"C"+postkey) != fittedEvents.end() )
      cout << " (MC exp. = " <<  fittedEvents[prekey+"C1"+postkey]/fittedEvents[prekey+"C"+postkey] << ")" << endl;

    normFactors["ABCD"][entry->first] = normVals["ABCD"][entry->first]->getVal();

    if(entry->second.find(prekey+"ABCD"+postkey) != entry->second.end()){
      channelByChannelTemplateDistroNormFit[entry->first][prekey+"ABCD"+postkey] = (TH1*)entry->second[prekey+"ABCD"+postkey]->Clone();

      if ( !channelByChannelTemplateDistroNormFit[entry->first][prekey+"ABCD"+postkey]->GetSumw2N() ) channelByChannelTemplateDistroNormFit[entry->first][prekey+"ABCD"+postkey]->Sumw2();

      double integral = getIntegral(channelByChannelTemplateDistroNormFit[entry->first][prekey+"ABCD"+postkey], true, true);
      if ( integral != 0. ) channelByChannelTemplateDistroNormFit[entry->first][prekey+"ABCD"+postkey]->Scale(normFactors["ABCD"][entry->first]/integral);
    }

    for ( vpset::const_iterator region = regionsPset.begin(); region != regionsPset.end(); ++region ) {
      const string & regionName = region->getParameter<string>("name");
      if( !(region->getParameter<bool>("toFit")) || regionName == "ABCD") continue;
      string key = prekey+regionName+postkey;
      double normRegion = normVals["ABCD"][entry->first]->getVal();
      
      if        ( regionName.find("A") != string::npos ) {
	normRegion *= pDiTauCharge_OS_SS[entry->first]->getVal();
	normRegion *= (1. - pMuonIso_tight_loose[entry->first]->getVal());
      } else if ( regionName.find("B") != string::npos ) {
	normRegion *= (1. - pDiTauCharge_OS_SS[entry->first]->getVal());
	normRegion *= (1. - pMuonIso_tight_loose[entry->first]->getVal());
      } else if ( regionName.find("C") != string::npos ) {
	normRegion *= pDiTauCharge_OS_SS[entry->first]->getVal();
	normRegion *= pMuonIso_tight_loose[entry->first]->getVal();
      } else if ( regionName.find("D") != string::npos ) {
	normRegion *= (1. - pDiTauCharge_OS_SS[entry->first]->getVal());
	normRegion *= pMuonIso_tight_loose[entry->first]->getVal();
      }
      
      if      ( regionName.find("1") != string::npos ) normRegion *= pDiTauKine_Sig_Bgr[entry->first]->getVal();
      else if ( regionName.find("2") != string::npos ) normRegion *= (1. - pDiTauKine_Sig_Bgr[entry->first]->getVal());
      
      if      ( regionName.find("p") != string::npos ) normRegion *= pTauId_passed_failed[entry->first]->getVal();
      else if ( regionName.find("f") != string::npos ) normRegion *= (1. - pTauId_passed_failed[entry->first]->getVal());
      
      normFactors[regionName][entry->first] = normRegion;

      cout << "--> "<< entry->first <<" = " << normFactors[regionName][entry->first]
		<< " (MC exp. = " << numEventsAll[key] << ")" << endl;

      channelByChannelTemplateDistroNormFit[entry->first][key] = (TH1*)entry->second[key]->Clone();
      if ( !channelByChannelTemplateDistroNormFit[entry->first][key]->GetSumw2N() ) channelByChannelTemplateDistroNormFit[entry->first][key]->Sumw2();
      double integral = getIntegral(channelByChannelTemplateDistroNormFit[entry->first][key], true, true);
      if ( integral != 0. ) channelByChannelTemplateDistroNormFit[entry->first][key]->Scale(normFactors[regionName][entry->first]/integral);

    }//pset loop
  }//template loop

  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();
 
  for ( vpset::const_iterator region = regionsPset.begin(); region != regionsPset.end(); ++region ) {
    if( !(region->getParameter<bool>("toFit")) ) continue;
    const string & regionName = region->getParameter<string>("name");
    TDirectory *regionDir = varDir->mkdir(regionName.c_str());

    DrawHistograms(regionDir,                             //TDirectory* dir
 		   channelByChannelTemplateDistroNormFit, //map<string,map<string, TH1* > > channelByChannelTemplateDistro
		   distributionsData,                     //map<string, TH1* > toFitDistro
		   dataChannelName,                       //string toFitChannelName
		   regionName,                            // const string& region
		   tauId,                                 //const string& tauId
		   varname,                               // string& varname
		   sysShift);                             //const string& sysShift
  }

  for( map<string,map<string, TH1* > >::const_iterator histomap = channelByChannelTemplateDistroNormFit.begin(); histomap != channelByChannelTemplateDistroNormFit.end(); histomap++)
    for(map<string, TH1* >::const_iterator histo = histomap->second.begin(); histo != histomap->second.end(); histo++)
      delete histo->second;
}

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // parse arguments
  if ( argc < 2 ) {
    cout << "Usage : " << argv[0] << " [histoFitter_cfg.py]" << endl;
    return 0;
  }

 if( !readPSetsFrom(argv[1])->existsAs<ParameterSet>("process") ){
   cout << " ERROR: ParametersSet 'process' is missing in your configuration file" << endl; 
   exit(0);
  }
 // get the python configuration
 const ParameterSet& process_ = readPSetsFrom(argv[1])->getParameter<ParameterSet>("process");

 // now get each parameter
 const ParameterSet& parameters_ = process_.getParameter<ParameterSet>("histoFitter");

 const vpset& channels_ = process_.getParameter<vpset>("channels"); //{name, totalEvevts,[Xsec:MC],dataset(?),inputFileName,fixedRatio(?)-->da vedere cosa ne ho bisogno in fitting,label,histoColor (fill)}
 const double& lumi_ = process_.getParameter<double>("lumi");
 const vpset& regionsPset_ = process_.getParameter<vpset>("regions"); //{name,toFit}
 vstring regions_;
 for(vpset::const_iterator pset = regionsPset_.begin(); pset != regionsPset_.end(); pset++)
   regions_.push_back( pset->getParameter<string>("name") );
 const vstring& tauids_ = process_.getParameter<vstring>("tauIds");
 const vstring& sysShifts_ = process_.getParameter<vstring>("sysShifts");
 const vpset& fitVariables_ = process_.getParameter<vpset>("fitVariables"); //{name, xAxisTitle,rebinFactor}
 const bool & runClosureTest_ = process_.getParameter<bool>("runClosureTest");
 const string& outFile_ = process_.getParameter<string>("outFile");
 const ParameterSet& fitSettings_ = process_.getParameter<ParameterSet>("fitSettings");

 map<string, TH1* > allDistributions;
 map<string,map<string, TH1* > > channelByChannelTemplateDistro;
 map<string, TH1* > toFitDistro;
 vstring channelNames;
 vstring templateNames;
 map<string, TFile*> inputFiles; //channelName
 map<string, double> totalEvents; //channelName
 map<string, float> analyzedEvents;

 for(vpset::const_iterator channel = channels_.begin(); channel != channels_.end(); ++channel){   
   const string& channelName = channel->getParameter<string>("name");
   if(channelName != "data" || !runClosureTest_){
     channelNames.push_back(channelName);
     const string& infileName = channel->getParameter<string>("inputFile");
     totalEvents[channelName] = channel->getParameter<double>("totalEvents");
     TFile *tempFile = new TFile(infileName.c_str());
     if(tempFile != 0){
       if(channelName != "data"){
	 channelByChannelTemplateDistro[channelName] = loadHistograms(  tempFile,*channel, regions_,  tauids_, fitVariables_, sysShifts_);
	 analyzedEvents[channelName] = channelByChannelTemplateDistro[channelName][(channelName+"/AnalyzedEvents").c_str()]->GetBinContent(1);
	 allDistributions.insert(channelByChannelTemplateDistro[channelName].begin(), channelByChannelTemplateDistro[channelName].end());
       }
       else{
	 toFitDistro = loadHistograms(  tempFile,*channel, regions_,  tauids_, fitVariables_, sysShifts_);
	 analyzedEvents[channelName] = toFitDistro[(channelName+"/AnalyzedEvents").c_str()]->GetBinContent(1);
	 allDistributions.insert(toFitDistro.begin(), toFitDistro.end());
       }
       inputFiles[channelName] = tempFile;
     }
     else{
       cout << "Problem loadind file: "<< infileName << endl;
       exit(0);
     }
   }//if(channelName != "data" && runClosureTest_)
 }//channel loop
 
 //scales due to correction factors - dead crab jobs
 double lumiOnData = (runClosureTest_) ? lumi_ : lumi_*analyzedEvents["data"]/totalEvents["data"];
 for(vpset::const_iterator channel = channels_.begin(); channel != channels_.end(); ++channel){ 
   const string& channelName = channel->getParameter<string>("name");
   if(channelName != "data" ){
     const double& xsec = channel->getParameter<double>("xsec");
     double scaleFactor = lumiOnData * xsec / analyzedEvents[channelName];
     for(map<string, TH1* >::iterator histo = channelByChannelTemplateDistro[channelName].begin(); histo != channelByChannelTemplateDistro[channelName].end(); histo++)
       histo->second->Scale(scaleFactor);
   }
 }

 //creates a fake data distribution
 if(runClosureTest_){
   string channelName("fakeData");
   toFitDistro = sumHistograms(  channelNames,allDistributions,channelName, regions_,  tauids_, fitVariables_, sysShifts_);
   channelNames.push_back(channelName);
   allDistributions.insert(toFitDistro.begin(), toFitDistro.end());
 }

 TFile* nofitFile= new TFile(("NoFIT_"+outFile_).c_str(),"recreate");
 string toFitChannelName = (runClosureTest_) ? "fakeData" : "data";
 // Dir Structure: tauid.name/sysShift/plotVariable.name/region/plotVariable.name+"Distribution"(+"_"+sysShift)
 for ( vstring::const_iterator tauId = tauids_.begin(); tauId != tauids_.end(); ++tauId ) {
   TDirectory* tauDir = nofitFile->mkdir(tauId->c_str());
   for ( vstring::const_iterator sysShift = sysShifts_.begin(); sysShift != sysShifts_.end(); ++sysShift ) {
     TDirectory* sysDir = tauDir->mkdir(sysShift->c_str());
     for ( vpset::const_iterator fitVariable = fitVariables_.begin(); fitVariable != fitVariables_.end(); ++fitVariable ) {
	const string& varname = fitVariable->getParameter<string>("name");
	TDirectory* varDir = sysDir->mkdir(varname.c_str());
	for ( vstring::const_iterator region = regions_.begin(); region != regions_.end(); ++region ) {
	  TDirectory* regionDir = varDir->mkdir(region->c_str());
	  string shift = *sysShift;
	  DrawHistograms(regionDir,
			 channelByChannelTemplateDistro,
			 toFitDistro,
			 toFitChannelName,
			 *region,
			 *tauId,
			 varname,
			 shift);
	}//region Loop
     }//var loop
   }//sys loop
 }//tauid loop


 map<string, double> numEventsAll;
 map<string, double> fittedEvents;
 for(map<string, TH1* >::iterator entry = allDistributions.begin(); entry != allDistributions.end(); entry++){
   numEventsAll[entry->first] = getIntegral(entry->second, true, true);
   fittedEvents[entry->first] = getIntegral(entry->second);
   cout << "Histo key: "<<entry->first << "    fitted fraction: "<<fittedEvents[entry->first] <<endl;
 }
 cout<<endl;



 for ( vector<string>::const_iterator tauId = tauids_.begin(); tauId != tauids_.end(); ++tauId ) {
   for ( vector<string>::const_iterator channelName = channelNames.begin(); channelName != channelNames.end(); ++channelName ) {
     for(vstring::const_iterator sysShift=sysShifts_.begin(); sysShift!=sysShifts_.end(); sysShift++){
       // map = (channel/tauid.name/region/observable_sysShift)
       string preKey = *channelName+"/"+*tauId+"/";
       string postKey = "diTauMt_"+*sysShift;
       double numEventsA  = numEventsAll[preKey+"A"+postKey];
       double numEventsB  = numEventsAll[preKey+"B"+postKey];
       double numEventsC  = numEventsAll[preKey+"C"+postKey];
       double numEventsC1 = numEventsAll[preKey+"C1"+postKey];
       double numEventsC2 = numEventsAll[preKey+"C2"+postKey];
       double numEventsD  = numEventsAll[preKey+"D"+postKey];
     
       cout << "process = " << (*channelName) << endl;
       for ( size_t i = 0; i < (channelName->length() + 10); ++i ) {
	 cout << "-";
       }
       cout << endl;
       
       cout << "pDiTauCharge_OS_SS:" << endl;
       cout << " A/B = " << numEventsA/numEventsB << endl;
       cout << " C/D = " << numEventsC/numEventsD << endl;
       cout << "pMuonIso_tight_loose:" << endl;
       cout << " C/(A+C) = " << numEventsC/(numEventsA + numEventsC) << endl;
       cout << " D/(B+D) = " << numEventsD/(numEventsB + numEventsD) << endl;
       cout << "pDiTauKine_Sig_Bgr:" << endl;
       cout << " C1/C = " << numEventsC1/numEventsC << endl;
       
       cout << "pTauId_passed_failed:" << endl;
       
       double numEventsC1p = numEventsAll[preKey+"C1p"+postKey];
       double numEventsC2p = numEventsAll[preKey+"C2p"+postKey];
       
       cout << " " << (*tauId) << ":" << endl;
       cout << "  C1p/C1 = " << numEventsC1p/numEventsC1 << endl;
       cout << "  C2p/C2 = " << numEventsC2p/numEventsC2 << endl;
       
       cout << endl;
     }//sysshift
   }//channel
 }//tauID

 map<string, double> effValues; //tauId,var,sys
 map<string, double> effErrors;

 TFile* fitFile= new TFile(outFile_.c_str(),"recreate");
 for(vstring::const_iterator sysShift=sysShifts_.begin(); sysShift!=sysShifts_.end(); sysShift++){
   TDirectory* sysShiftDir = fitFile->mkdir(sysShift->c_str());
   for ( vector<string>::const_iterator tauId = tauids_.begin(); tauId != tauids_.end(); ++tauId ) {
     TDirectory* tauDir = sysShiftDir->mkdir(tauId->c_str());
     for ( vpset::const_iterator fitVariable = fitVariables_.begin(); fitVariable != fitVariables_.end(); ++fitVariable ) {
       double effValue = 0.;
       double effError = 1.;
       const string shift = *sysShift;



       /*fitUsingRooFit(TDirectory* tauDir,
		    map<string, TH1*> & distributionsData,
		    string & dataChannelName,
		    map<string,map<string, TH1* > > & channelByChannelTemplateDistro,
		    map<string, double> & numEventsAll,
		    map<string, double>& fittedEvents,
		    const string& tauId, 
		    const ParameterSet& fitVariable,
		    const vpset& regionsPset,
		    const ParameterSet & fitParameters,
		    bool fitTauIdEffC2,
		    double& effValue,
		    double& effError,
		    const string& sysShift,
		    bool isClosureTest = false) */


       fitUsingRooFit(tauDir,                         // TDirectory* tauDir,
		      toFitDistro,                    //  map<string, TH1*> & distributionsData,
		      toFitChannelName,               // string & dataChannelName,
		      channelByChannelTemplateDistro, // map<string,map<string, TH1* > > & channelByChannelTemplateDistro,
		      numEventsAll,                   // map<string, double> & numEventsAll,
		      fittedEvents,                   // map<string, double>& fittedEvents,
		      *tauId,                         // const string& tauId, 
		      *fitVariable,                   // const ParameterSet& fitVariable,
		      regionsPset_,                   //  const vpset& regionsPset,
		      fitSettings_,                   // const ParameterSet & fitParameters,
		      effValue,                       // bool fitTauIdEffC2,
		      effError,		              //
		      shift,                          //
		      runClosureTest_);               //        bkmrk

     }//var
   }//tauid
 }//sysshift

 for(map<string, TH1*>::iterator histo = allDistributions.begin(); histo != allDistributions.end(); histo++)
   delete histo->second;

 for(map<string, TFile*>::iterator file = inputFiles.begin(); file != inputFiles.end(); file++){
   file->second->Close();
   delete file->second;
 }
 fitFile->Close();
}

//OLD CODE THAT MIGHT TURN USEFUL IF I CHANGE MIND ;)
/*/-------------------------------Zmumu template-------------------------------
 map<string, TH1* > templates_zmumu;
 const ParameterSet& zmumu_pset = channels_.getParameter<ParameterSet>("Zmumu");
 const string& zmumu_infile = data.getParameter<string>("inputFile");
 string& zmumu_name = const_cast<string&>(data.getParameter<string>("name"));
 const double zmumu_EvtTot = data.getParameter<double>("totalEvents");
 TFile *zmumu_PlotFile = new TFile(infile.c_str());
 if(zmumu_PlotFile != 0)
   distributionsData = loadHistograms(  zmumu_PlotFile,zmumu_name, regions_,  tauids_, fitVariables_, sysShifts_);
 else{
   cout << "Problem loadind file: "<< zmumu_infile << endl;
   exit(0);
 }
 float zmumu_EvtAnalyzed = templates_zmumu["General"][(zmumu_name+"/AnalyzedEvents").c_str()]->GetBinContent(1);
 allDistributions.insert(templates_zmumu.begin(), templates_zmumu.end());
 channels.push_back(zmumu_name);

 //-------------------------------Qcd template-------------------------------
 map<string, TH1* > templates_qcd;
 const ParameterSet& qcd_pset = channels_.getParameter<ParameterSet>("QCD");
 const string& qcd_infile = data.getParameter<string>("inputFile");
 string& qcd_name = const_cast<string&>(data.getParameter<string>("name"));
 const double qcd_EvtTot = data.getParameter<double>("totalEvents");
 TFile *qcd_PlotFile = new TFile(infile.c_str());
 if(qcd_PlotFile != 0)
   distributionsData = loadHistograms(  qcd_PlotFile,qcd_name, regions_,  tauids_, fitVariables_, sysShifts_);
 else{
   cout << "Problem loadind file: "<< qcd_infile << endl;
   exit(0);
 }
 float qcd_EvtAnalyzed = templates_qcd["General"][(qcd_name+"/AnalyzedEvents").c_str()]->GetBinContent(1);
 allDistributions.insert(templates_qcd.begin(), templates_qcd.end());
 channels.push_back(qcd_name);

 //-------------------------------Wplusjets template-------------------------------
 map<string, TH1* > templates_wplusjets;
 const ParameterSet& wplusjets_pset = channels_.getParameter<ParameterSet>("WplusJets");
 const string& wplusjets_infile = data.getParameter<string>("inputFile");
 string& wplusjets_name = const_cast<string&>(data.getParameter<string>("name"));
 const double wplusjets_EvtTot = data.getParameter<double>("totalEvents");
 TFile *wplusjets_PlotFile = new TFile(infile.c_str());
 if(wplusjets_PlotFile != 0)
   distributionsData = loadHistograms(  wplusjets_PlotFile,wplusjets_name, regions_,  tauids_, fitVariables_, sysShifts_);
 else{
   cout << "Problem loadind file: "<< wplusjets_infile << endl;
   exit(0);
 }
 float wplusjets_EvtAnalyzed = templates_wplusjets["General"][(wplusjets_name+"/AnalyzedEvents").c_str()]->GetBinContent(1);
 allDistributions.insert(templates_wplusjets.begin(), templates_wplusjets.end());
 channels.push_back(wplusjets_name);

 //-------------------------------Ttplusjets template-------------------------------
 map<string, TH1* > templates_ttplusjets;
 const ParameterSet& ttplusjets_pset = channels_.getParameter<ParameterSet>("TTplusJets");
 const string& ttplusjets_infile = data.getParameter<string>("inputFile");
 string& ttplusjets_name = const_cast<string&>(data.getParameter<string>("name"));
 const double ttplusjets_EvtTot = data.getParameter<double>("totalEvents");
 TFile *ttplusjets_PlotFile = new TFile(infile.c_str());
 if(ttplusjets_PlotFile != 0)
   distributionsData = loadHistograms(  ttplusjets_PlotFile,ttplusjets_name, regions_,  tauids_, fitVariables_, sysShifts_);
 else{
   cout << "Problem loadind file: "<< ttplusjets_infile << endl;
   exit(0);
 }
 float ttplusjets_EvtAnalyzed = templates_ttplusjets["General"][(ttplusjets_name+"/AnalyzedEvents").c_str()]->GetBinContent(1);
 allDistributions.insert(templates_ttplusjets.begin(), templates_ttplusjets.end());
 channels.push_back(ttplusjets_name);

 map<string, TH1* > distributionsData;
 TFile *dataPlotFile = NULL;
 float dataEvtAnalyzed;
 double& dataEvtTot;
 string& nameData;
 if( !runClosureTest_){
   const ParameterSet& data = channels_.getParameter<ParameterSet>("data");
   const string& infile = data.getParameter<string>("inputFile");
   nameData = const_cast<string&>(data.getParameter<string>("name"));
   channels.push_back(nameData);
   dataEvtTot = const_cast<double&>(data.getParameter<double>("totalEvents") ); 
   dataPlotFile = new TFile(infile.c_str());
   if(dataPlotFile != 0)
     distributionsData = loadHistograms(  dataPlotFile,name, regions_,  tauids_, fitVariables_, sysShifts_);
   else{
     cout << "Problem loadind file: "<< infile << endl;
     exit(0);
   }
   dataEvtAnalyzed = distributionsData["General"][(name+"/AnalyzedEvents").c_str()]->GetBinContent(1);
   allDistributions.insert(templates_zmumu.begin(), templates_zmumu.end());
   channels.push_back(nameData);
 }
 else{

 }
 if(dataPlotFile != 0){
   dataPlotFile->Close();
   delete dataPlotFile;
 }
 zmumu_PlotFile->Close();
 delete zmumu_PlotFile;
 qcd_PlotFile->Close();
 delete qcd_PlotFile;
 wplusjets_PlotFile->Close();
 delete wplusjets_PlotFile;
 ttplusjets_PlotFile->Close();
 delete ttplusjets_ttplusjets;
*/
 
