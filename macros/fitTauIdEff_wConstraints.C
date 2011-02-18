
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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

std::string getBranchName(const std::string& branchName_prefix,
                          const std::string& branchName_object, const std::string& branchName_observable,
			  const std::string& branchName_suffix)
{
//--------------------------------------------------------------------------------
// Compose full branchName from name of pat::Muon/pat::Tau/diTau collection
// used when producing (ED)Ntuple, observable and branchName "suffix" ("local"/"lxbatch")
//--------------------------------------------------------------------------------

  //std::cout << "<getBranchName>:" << std::endl;

  const std::string ntupleName = "ntupleProducer_tauIdEffNtuple";

  std::string branchName = branchName_prefix;
  branchName.append("_").append(ntupleName);
  branchName.append("#").append(branchName_object);
  branchName.append("#").append(branchName_observable);
  branchName.append("_").append(branchName_suffix).append(".obj");

  return branchName;
}

std::string replace(std::string& str, const std::string& oldpart, const std::string& newpart) 
{
  //std::cout << "<replace>:" << std::endl;

  std::string retVal = str;

  size_t idx = retVal.find(oldpart);
  while ( idx != std::string::npos ) {
    retVal.replace(idx, oldpart.length(), newpart);
    idx += newpart.length();
    if ( idx < retVal.length() )
      idx = retVal.find(oldpart, idx);
    else
      idx = std::string::npos;
    //std::cout << "retVal = " << retVal << std::endl;
  }

  return retVal;
}

std::map<std::string, std::map<std::string, std::string> > makeBranchNameDict(
  const std::vector<std::string>& tauIds,
  const std::string& sysShift, const std::string& branchName_suffix,
  bool applyZrecoilCorrections)
{
//--------------------------------------------------------------------------------
// Make alias --> branchName mapping 
// to be used for treeSelection and histogram filling
//--------------------------------------------------------------------------------

  //std::cout << "<makeBranchNameDict>:" << std::endl;

  typedef std::map<std::string, std::string> branchNameDictEntry;

  std::map<std::string, branchNameDictEntry> retVal;

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {

    std::string branchNameMuon  = "selectedPatMuonsForTauIdEffTrkIPcumulative";
    std::string branchNameTau   = "";
    std::string branchNameDiTau = "";
    if (        (*tauId) == "tauDiscrTaNCfrOnePercent"     ||
	        (*tauId) == "tauDiscrTaNCfrHalfPercent"    ||
	        (*tauId) == "tauDiscrTaNCfrQuarterPercent" ) { // "old" TaNC algorithm
      branchNameTau   = "selectedPatPFTausShrinkingConeElectronVetoCumulative"; // "old" TaNC algorithm
      branchNameDiTau = "selectedMuPFTauShrinkingConePairsForTauIdEffCumulative";  
    } else if ( (*tauId) == "tauDiscrTaNCloose"            || 
		(*tauId) == "tauDiscrTaNCmedium"           || 
		(*tauId) == "tauDiscrTaNCtight"            ) { // "new" TaNC implemented in HPS+TaNC combined algorithm
      branchNameTau   = "selectedPatPFTausHPSpTaNCCaloMuonVetoCumulative";    
      branchNameDiTau = "selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative"; 
    } else if ( (*tauId) == "tauDiscrIsolationLoose"       ||
		(*tauId) == "tauDiscrIsolationMedium"      ||
		(*tauId) == "tauDiscrIsolationTight"       ) { // "old" HPS algorithm
      branchNameTau   = "selectedPatPFTausHPSElectronVetoCumulative";           // "old" HPS algorithm 
      branchNameDiTau = "selectedMuPFTauHPSpairsForTauIdEffCumulative";
    } else if ( (*tauId) == "tauDiscrHPSloose"             ||
		(*tauId) == "tauDiscrHPSmedium"            || 
		(*tauId) == "tauDiscrHPStight"             ) { // "new" HPS implemented in HPS+TaNC combined algorithm
      branchNameTau   = "selectedPatPFTausHPSpTaNCCaloMuonVetoCumulative"; 
      branchNameDiTau = "selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative";
    } else {
      std::cout << "Error in <makeBranchNameDict>: invalid tauId = " << (*tauId) << " --> aborting !!";
      assert(0);
    }
    
    if ( applyZrecoilCorrections ) branchNameDiTau = replace(branchNameDiTau, "Cumulative", "ZllRecoilCorrectedCumulative");
    
    if (        sysShift == "CENTRAL_VALUE"              ) {
      // nothing to be done yet.
    } else if ( sysShift == "SysTauJetEnUp"              ||
		sysShift == "SysTauJetEnDown"            ) {
      branchNameTau   = replace(branchNameTau,   "Cumulative", std::string(sysShift).append("Cumulative"));
      branchNameDiTau = replace(branchNameDiTau, "Cumulative", std::string(sysShift).append("Cumulative"));
    } else if ( sysShift == "SysJetEnUp"                 ||
		sysShift == "SysJetEnDown"               ) {
      branchNameTau   = replace(branchNameTau,   "Cumulative", std::string(sysShift).append("Cumulative"));
      branchNameDiTau = replace(branchNameDiTau, "Cumulative", std::string(sysShift).append("Cumulative"));
    } else if ( sysShift == "SysZllRecoilCorrectionUp"   ||
		sysShift == "SysZllRecoilCorrectionDown" ) {
      branchNameTau   = replace(branchNameTau,   "Cumulative", std::string(sysShift).append("Cumulative"));
      branchNameDiTau = replace(branchNameDiTau, "Cumulative", std::string(sysShift).append("Cumulative"));
    } else assert(0);
    
    branchNameDictEntry branchNames;

    branchNames["event"] = getBranchName("double", "", "event", branchName_suffix);
    branchNames["ls"] = getBranchName("double", "", "ls", branchName_suffix);
    branchNames["run"] = getBranchName("double", "", "run", branchName_suffix);
    
    branchNames["muonPt"] = getBranchName("double", branchNameMuon, "pt", branchName_suffix);
    branchNames["muonEta"] = getBranchName("double", branchNameMuon, "eta", branchName_suffix);
    branchNames["muonLooseIsoPtSum04"] = getBranchName("double", branchNameMuon, "ptSumLooseIsolation04", branchName_suffix);
    branchNames["muonLooseIsoPtSum06"] = getBranchName("double", branchNameMuon, "ptSumLooseIsolation06", branchName_suffix);
    
    branchNames["tauPt"] = getBranchName("doubles", branchNameTau, "pt", branchName_suffix);
    branchNames["tauEta"] = getBranchName("doubles", branchNameTau, "eta", branchName_suffix);
    branchNames["tauLooseIsoPtSum04"] = getBranchName("doubles", branchNameTau, "ptSumLooseIsolation04", branchName_suffix);
    branchNames["tauLooseIsoPtSum06"] = getBranchName("doubles", branchNameTau, "ptSumLooseIsolation06", branchName_suffix);
    branchNames["tauNumChargedParticles"] = getBranchName("doubles", branchNameTau, "numChargedParticles", branchName_suffix);
    branchNames["tauNumParticles"] = getBranchName("doubles", branchNameTau, "numParticles", branchName_suffix);
    branchNames["tauJetPt"] = getBranchName("doubles", branchNameTau, "jetPt", branchName_suffix);
    branchNames["tauJetEta"] = getBranchName("doubles", branchNameTau, "jetEta", branchName_suffix);
    branchNames["tauJetWidth"] = getBranchName("doubles", branchNameTau, "jetWidth", branchName_suffix);
    branchNames["tauDiscrTaNCfrOnePercent"] = getBranchName("doubles", branchNameTau, "byTaNCfrOnePercent", branchName_suffix);
    branchNames["tauDiscrTaNCfrHalfPercent"] = getBranchName("doubles", branchNameTau, "byTaNCfrHalfPercent", branchName_suffix);
    branchNames["tauDiscrTaNCfrQuarterPercent"] = getBranchName("doubles", branchNameTau, "byTaNCfrQuarterPercent", branchName_suffix);
    branchNames["tauDiscrTaNCloose"] = getBranchName("doubles", branchNameTau, "byTaNCloose", branchName_suffix);
    branchNames["tauDiscrTaNCmedium"] = getBranchName("doubles", branchNameTau, "byTaNCmedium", branchName_suffix);
    branchNames["tauDiscrTaNCtight"] = getBranchName("doubles", branchNameTau, "byTaNCtight", branchName_suffix);
    branchNames["tauDiscrIsolationLoose"] = getBranchName("doubles", branchNameTau, "byIsolationLoose", branchName_suffix);
    branchNames["tauDiscrIsolationMedium"] = getBranchName("doubles", branchNameTau, "byIsolationMedium", branchName_suffix);
    branchNames["tauDiscrIsolationTight"] = getBranchName("doubles", branchNameTau, "byIsolationTight", branchName_suffix);
    branchNames["tauDiscrHPSloose"] = getBranchName("doubles", branchNameTau, "byHPSloose", branchName_suffix);
    branchNames["tauDiscrHPSmedium"] = getBranchName("doubles", branchNameTau, "byHPSmedium", branchName_suffix);
    branchNames["tauDiscrHPStight"] = getBranchName("doubles", branchNameTau, "byHPStight", branchName_suffix);
    
    branchNames["diTauCharge"] = getBranchName("double", branchNameDiTau, "charge", branchName_suffix);
    branchNames["diTauMt"] = getBranchName("double", branchNameDiTau, "Mt", branchName_suffix);
    branchNames["diTauPzeta"] = getBranchName("double", branchNameDiTau, "pZeta", branchName_suffix);
    branchNames["diTauPzetaVis"] = getBranchName("double", branchNameDiTau, "pZetaVis", branchName_suffix);
    branchNames["diTauHt"] = getBranchName("double", branchNameDiTau, "Ht", branchName_suffix);
    branchNames["diTauSVfitMass1"] = getBranchName("double", branchNameDiTau, "SVfitMass1", branchName_suffix);
    branchNames["diTauSVfitMass2"] = getBranchName("double", branchNameDiTau, "SVfitMass2", branchName_suffix);
    branchNames["diTauVisMass"] = getBranchName("double", branchNameDiTau, "visMass", branchName_suffix);
    branchNames["diTauVisMassFromJet"] = getBranchName("double", branchNameDiTau, "visMassFromJet", branchName_suffix);    

    retVal[*tauId] = branchNames;
  }

  return retVal;
}


void writeRunLumiSectionEventNumberFile(const std::string& process, TTree* tree, const std::string& treeSelection,
					const std::string& region, 
					const std::string& tauId, const std::string& tauIdValue,
					std::map<std::string, std::string>& branchNames)
{
  //std::cout << "<writeRunLumiSectionEventNumberFile>:" << std::endl;

  std::string outputFileName = std::string("selEvents_").append(process);
  outputFileName.append("_").append(region);
  outputFileName.append("_").append(tauId).append("_").append(tauIdValue).append(".txt");
  ofstream* outputFile = new ofstream(outputFileName.data());

  std::string drawCommand = branchNames["run"];
  drawCommand.append(":").append(branchNames["ls"]);
  drawCommand.append(":").append(branchNames["event"]);

  tree->Draw(drawCommand.data(), treeSelection.data());

  TPolyMarker3D* tmpPolyMarker = dynamic_cast<TPolyMarker3D*>(gPad->GetPrimitive("TPolyMarker3D"));
  if ( !tmpPolyMarker ) {
    std::cout << "Error in <writeRunLumiSectionEventNumberFile>: failed to create TPolyMarker3D --> skipping !!" << std::endl;
    return;
  }
 
  double run, ls, event;

  int numEvents = tmpPolyMarker->GetN();
  for ( int iEvent = 0 ; iEvent < numEvents; ++iEvent ) {
    tmpPolyMarker->GetPoint(iEvent, event, ls, run); // NOTE: order of run, ls, event triplets is reversed !!
    
    *outputFile << TMath::Nint(run) << ":" << TMath::Nint(ls) << ":" << TMath::Nint(event) << std::endl;
  }
  
  delete outputFile;
}

std::string getKey(const std::string& observable, const std::string& tauId, const std::string tauIdValue = "all")
{
  //std::cout << "<getKey>:" << std::endl;

  std::string key = std::string(observable).append("_").append(tauId).append("_").append(tauIdValue);
  return key;
}

double getIntegral(const TH1* histogram, bool inclUnderflowBin, bool inclOverflowBin)
{
//--------------------------------------------------------------------------------
// Compute integral of histogram including/excluding underflow and overflow bins
//--------------------------------------------------------------------------------

  //std::cout << "<getIntegral>:" << std::endl;
  //std::cout << " histogram        = " << histogram << std::endl;
  //std::cout << " name             = " << histogram->GetName() << std::endl;
  //std::cout << " inclUnderflowBin = " << inclUnderflowBin << std::endl;
  //std::cout << " inclOverflowBin  = " << inclOverflowBin << std::endl;

  int firstBin = ( inclUnderflowBin ) ?                            0 :                      1;
  int lastBin  = ( inclOverflowBin  ) ? (histogram->GetNbinsX() + 1) : histogram->GetNbinsX();
  
  //std::cout << " firstBin = " << firstBin << ", lastBin = " << lastBin << std::endl;

  double integral = 0;
  for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    //std::cout << " binContent(" << iBin << ") = " << binContent << std::endl;
    integral += binContent;
  }

  //std::cout << "--> integral = " << integral << std::endl;

  return integral;
}

std::map<std::string, TH1*> makeHistograms(
  const std::string& process, const std::string& region, double weight,
  TTree* tree, const std::string& treeSelection, 
  const std::string& tauId, const std::vector<std::string>& tauIdValues,
  const std::vector<std::string>& observables,
  std::map<std::string, std::string>& branchNames, const std::string& sysShift = "CENTRAL_VALUE",
  bool saveRunLumiSectionEventNumbers = false)
{
//--------------------------------------------------------------------------------
// Fill histograms with (ED)NTuple entries passing treeSelection
//--------------------------------------------------------------------------------

  //std::cout << "<makeHistograms>:" << std::endl;
  //std::cout << " process = " << process << std::endl;
  //std::cout << " region = " << region << std::endl;

  std::map<std::string, TH1*> retVal;

//--- prepare "dummy" tau id. selection
//    in case no tau id. selection has been passed as function argument
  std::vector<std::string> tauIdValues_local;
  if ( tauIdValues.size() > 0 ) {
    tauIdValues_local = tauIdValues;
  } else {
    tauIdValues_local.push_back("all");
  }

  for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues_local.begin();
	tauIdValue != tauIdValues_local.end(); ++tauIdValue ) {

//--- add kinematic cuts common to all regions
//   (cuts that might as well have been applied during (ED)Ntuple production)
    std::string extTreeSelection = treeSelection;
    if ( extTreeSelection != "" ) extTreeSelection.append(" && ");
    extTreeSelection.append(branchNames["muonPt"]).append(" > 20.");
    extTreeSelection.append(" && ").append(branchNames["tauPt"]).append(" > 20.");
    extTreeSelection.append(" && abs(").append(branchNames["tauEta"]).append(") < 2.3");
    extTreeSelection.append(" && ").append(branchNames["tauLooseIsoPtSum06"]).append(" < 2.5");
    extTreeSelection.append(" && ").append(branchNames["diTauHt"]).append(" > 40.");

//--- add tau id. passed/failed selection
    if      ( (*tauIdValue) == "passed" ) extTreeSelection.append(" && ").append(branchNames[tauId]).append(" > 0.5");
    else if ( (*tauIdValue) == "failed" ) extTreeSelection.append(" && ").append(branchNames[tauId]).append(" < 0.5");
    else if ( (*tauIdValue) == "all"    ) {}
    else assert(0);
    
    std::cout << " treeSelection = " << extTreeSelection << std::endl;

    for ( std::vector<std::string>::const_iterator observable = observables.begin();
	  observable != observables.end(); ++observable ) {
      
      std::string histogramName = std::string(process).append("_").append(region).append("_").append(*observable);
      histogramName.append("_").append(tauId).append("_").append(*tauIdValue);
      if ( sysShift != "CENTRAL_VALUE" ) histogramName.append("_").append(sysShift);
      
      int numBins;
      double min, max;
      
      if ( (*observable) == "diTauMt" ) {
	numBins = 16;
	min = 0.;
	max = 80.;
      } else if ( (*observable) == "diTauVisMass"        ||
		  (*observable) == "diTauVisMassFromJet" ) {
	numBins = 18;
	min = 20.;
	max = 200.;
      } else if ( (*observable) == "diTauHt"         || 
		  (*observable) == "diTauSVfitMass1" || 
		  (*observable) == "diTauSVfitMass2" ) {
	numBins = 18;
	min = 20.;
	max = 200.;
      } else if ( (*observable) == "muonPt"  ||
		  (*observable) == "tauPt"   ||
		  (*observable) == "tauJetPt" ) {
	numBins = 12;
	min =  0.;
	max = 60.;
      } else if ( (*observable) == "muonEta" ) {
	numBins = 21;
	min = -2.1;
	max = +2.1;
      } else if ( (*observable) == "tauEta"   ||
		  (*observable) == "tauJetEta" ) {
	numBins = 23;
	min = -2.3;
	max = +2.3;
      } else if ( (*observable) == "tauNumChargedParticles" ) {
	numBins = 15;
	min =  -0.5;
	max = +14.5;
      } else if ( (*observable) == "tauNumParticles"       ) {
	numBins = 25;
	min =  -0.5;
	max = +24.5;
      } else if ( (*observable) == "tauJetWidth" ) {
	numBins = 20;
	min = 0.;
	max = 0.50;
      } else {
	std::cout << "Error in <makeHistograms>: undefined observable = " << (*observable) << " --> skipping !!" << std::endl;
	return retVal;
      }
      
      TH1* histogram = new TH1F(histogramName.data(), histogramName.data(), numBins, min, max);
      
      std::string drawCommand = std::string(branchNames[*observable]).append(">>+").append(histogramName); 
      tree->Draw(drawCommand.data(), extTreeSelection.data());
      
      if ( !histogram->GetSumw2N() ) histogram->Sumw2();
      histogram->Scale(weight);

      double integral = getIntegral(histogram, true, true);
      double fittedFraction = ( integral > 0. ) ? getIntegral(histogram, false, false)/integral : -1.; 
      std::cout << "histogram = " << histogramName << ":" 
		<< " entries = " << histogram->GetEntries() << ", integral = " << integral 
		<< " (fitted fraction = " << fittedFraction << ")" << std::endl;
      
      std::string key = getKey(*observable, tauId, *tauIdValue);	
      if ( histogram != 0 ) retVal[key] = histogram;

      if ( histogram->GetEntries() > 0 && saveRunLumiSectionEventNumbers && sysShift == "CENTRAL_VALUE" ) 
	writeRunLumiSectionEventNumberFile(process, tree, extTreeSelection, region, tauId, *tauIdValue, branchNames);
    }
  }

  std::cout << std::endl;

  return retVal;
}

TH1* normalize(const TH1* histogram, double norm = 1.)
{
//--------------------------------------------------------------------------------
// Normalize histogram passed as function argument to unit area
// (for use as shape template)
//--------------------------------------------------------------------------------

  //std::cout << "<normalize>:" << std::endl;

  TH1* retVal = (TH1*)histogram->Clone();

  if ( !retVal->GetSumw2N() ) retVal->Sumw2();

  double integral = getIntegral(retVal, true, true);
  if ( integral != 0. ) retVal->Scale(norm/integral);
    
  return retVal;
}

void applyStyleOption(TH1* histogram, const std::string& histogramTitle,
		      const std::string& xAxisTitle, const std::string& yAxisTitle = "Number of Events")
{
  //std::cout << "<applyStyleOption>:" << std::endl;

  histogram->SetStats(false);

  //histogram->SetTitle(histogramTitle.data());
  histogram->SetTitle("");
    
  histogram->GetXaxis()->SetTitle(xAxisTitle.data());
  histogram->GetXaxis()->SetTitleOffset(1.15);
  //histogram->GetXaxis()->SetTitleSize(0.05); 
  //histogram->GetXaxis()->SetLabelSize(0.05);

  histogram->GetYaxis()->SetTitle(yAxisTitle.data());
  histogram->GetYaxis()->SetTitleOffset(1.20);
  //histogram->GetYaxis()->SetTitleSize(0.05); 
  //histogram->GetYaxis()->SetLabelSize(0.05);
}

void drawHistograms(TH1* histogramZtautau, double normZtautau,
		    TH1* histogramZmumu, double normZmumu,
		    TH1* histogramQCD, double normQCD,
		    TH1* histogramWplusJets, double normWplusJets,
		    TH1* histogramTTplusJets, double normTTplusJets,
		    TH1* histogramData,
		    const std::string& histogramTitle, const std::string& xAxisTitle, 
		    const std::string& outputFileName, const std::string& sysShift = "CENTRAL_VALUE")
{
//--------------------------------------------------------------------------------
// Make control plots of sum(MC) versus Data.
// If normalization factors are passed as function argument,
// normalize all MC distributions accordingly;
// else assume MC distributions passed as function arguments
// are already properly normalized (by cross-section)
//--------------------------------------------------------------------------------

  //std::cout << "<drawHistograms>:" << std::endl;
  //std::cout << " Ztautau:    histogram = " << histogramZtautau    << ", norm = " << normZtautau    << std::endl;
  //std::cout << " Zmumu:      histogram = " << histogramZmumu      << ", norm = " << normZmumu      << std::endl;
  //std::cout << " QCD:        histogram = " << histogramQCD        << ", norm = " << normQCD        << std::endl;
  //std::cout << " WplusJets:  histogram = " << histogramWplusJets  << ", norm = " << normWplusJets  << std::endl;
  //std::cout << " TTplusJets: histogram = " << histogramTTplusJets << ", norm = " << normTTplusJets << std::endl;
  //std::cout << " Data:       histogram = " << histogramData << std::endl;
    
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

//--- scale template histograms to given normalization factors
  TH1* templateZtautau    = ( normZtautau    > 0. ) ? normalize(histogramZtautau,    normZtautau)    : histogramZtautau;
  applyStyleOption(templateZtautau, histogramTitle, xAxisTitle);
  templateZtautau->SetFillColor(628);
  
  TH1* templateZmumu      = ( normZmumu      > 0. ) ? normalize(histogramZmumu,      normZmumu)      : histogramZmumu;
  applyStyleOption(templateZmumu, histogramTitle, xAxisTitle);
  templateZmumu->SetFillColor(596);

  TH1* templateQCD        = ( normQCD        > 0. ) ? normalize(histogramQCD,        normQCD)        : histogramQCD;
  applyStyleOption(templateQCD, histogramTitle, xAxisTitle);
  templateQCD->SetFillColor(797);

  TH1* templateWplusJets  = ( normWplusJets  > 0. ) ? normalize(histogramWplusJets,  normWplusJets)  : histogramWplusJets;
  applyStyleOption(templateWplusJets, histogramTitle, xAxisTitle);
  templateWplusJets->SetFillColor(856);

  TH1* templateTTplusJets = ( normTTplusJets > 0. ) ? normalize(histogramTTplusJets, normTTplusJets) : histogramTTplusJets;
  applyStyleOption(templateTTplusJets, histogramTitle, xAxisTitle);
  templateTTplusJets->SetFillColor(618);

  applyStyleOption(histogramData, histogramTitle, xAxisTitle);
  histogramData->SetLineColor(1);
  histogramData->SetMarkerColor(1);
  histogramData->SetMarkerStyle(20);

//--- draw histograms for individual processes
  THStack smSum("smSum", "smSum");
  smSum.Add(templateTTplusJets);
  smSum.Add(templateZmumu);  
  smSum.Add(templateWplusJets);
  smSum.Add(templateQCD);
  smSum.Add(templateZtautau);

  smSum.SetTitle(templateZtautau->GetTitle());
  smSum.SetMaximum(1.4*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));

  smSum.Draw("hist");
  histogramData->SetStats(false);
  histogramData->Draw("ep1same");
  //smSum.Draw("axissame");

  TLegend legend(0.64, 0.64, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  legend.AddEntry(histogramData,      "Data",                            "p");
  legend.AddEntry(templateZtautau,    "Z #rightarrow #tau^{+} #tau^{-}", "f");
  legend.AddEntry(templateQCD,        "QCD",                             "f");
  legend.AddEntry(templateWplusJets,  "W + jets",                        "f");
  legend.AddEntry(templateZmumu,      "Z #rightarrow #mu^{+} #mu^{-}",   "f");
  legend.AddEntry(templateTTplusJets, "t#bar{t} + jets",                 "f");
  legend.Draw();

  TPaveText cmsPreliminaryLabel(0.135, 0.865, 0.46, 0.905, "NDC");
  cmsPreliminaryLabel.AddText("CMS Preliminary");
  cmsPreliminaryLabel.SetTextAlign(13);
  cmsPreliminaryLabel.SetTextSize(0.040);
  cmsPreliminaryLabel.SetFillStyle(0);
  cmsPreliminaryLabel.SetBorderSize(0);
  cmsPreliminaryLabel.Draw();

  TPaveText cmsLuminosityLabel(0.135, 0.8125, 0.46, 0.8525, "NDC");
  cmsLuminosityLabel.AddText("L = 36.1pb^{-1}");
  cmsLuminosityLabel.SetTextAlign(13);
  cmsLuminosityLabel.SetTextSize(0.035);
  cmsLuminosityLabel.SetFillStyle(0);
  cmsLuminosityLabel.SetBorderSize(0);
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  if ( sysShift != "CENTRAL_VALUE" ) outputFilePath.append(sysShift).append("/");
  gSystem->mkdir(outputFilePath.data(), true);
  canvas->Print(outputFilePath.append(outputFileName).data());

//--- draw histograms for all background processes summed plus Ztautau signal 
  TH1* templateSMbgSum = (TH1*)templateZmumu->Clone();
  templateSMbgSum->Add(templateQCD);
  templateSMbgSum->Add(templateWplusJets);
  templateSMbgSum->Add(templateTTplusJets);
  applyStyleOption(templateSMbgSum, histogramTitle, xAxisTitle);
  templateSMbgSum->SetFillColor(42);

  templateZtautau->SetFillColor(46);

  THStack smSum2("smSum", "smSum");
  smSum2.Add(templateSMbgSum);
  smSum2.Add(templateZtautau);

  smSum2.SetTitle(templateZtautau->GetTitle());
  smSum2.SetMaximum(1.4*TMath::Max(smSum2.GetMaximum(), histogramData->GetMaximum()));
	
  smSum2.Draw("hist");
  histogramData->Draw("ep1same");
  //smSum2.Draw("axissame");

  TLegend legend2(0.64, 0.74, 0.89, 0.89, "", "brNDC"); 
  legend2.SetBorderSize(0);
  legend2.SetFillColor(0);
  
  legend2.AddEntry(histogramData,   "Data",                            "p");
  legend2.AddEntry(templateZtautau, "Z #rightarrow #tau^{+} #tau^{-}", "f");
  legend2.AddEntry(templateSMbgSum, "#Sigma Backgrounds",              "f");
  legend2.Draw();

  cmsPreliminaryLabel.Draw();
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath2 = std::string("./plots/");
  if ( sysShift != "CENTRAL_VALUE" ) outputFilePath2.append(sysShift).append("/");
  gSystem->mkdir(outputFilePath2.data(), true);
  TString outputFileName2 = outputFileName;
  outputFileName2.ReplaceAll(".", "_smBgSum.");
  canvas->Print(outputFilePath2.append(outputFileName2).data());

//--- draw histograms for all signal and background processes summed
  TH1* templateSMsum = (TH1*)templateSMbgSum->Clone();
  templateSMsum->Add(templateZtautau);
  applyStyleOption(templateSMsum, histogramTitle, xAxisTitle);
  templateSMsum->SetFillColor(10);
  templateSMsum->SetLineColor(1);
  templateSMsum->SetLineWidth(2);

  templateSMsum->SetTitle(templateZtautau->GetTitle());
  templateSMsum->SetMaximum(1.4*TMath::Max(templateSMsum->GetMaximum(), histogramData->GetMaximum()));
	
  templateSMsum->Draw("hist");
  histogramData->Draw("ep1same");
  //templateSMsum->Draw("axissame");

  TLegend legend3(0.64, 0.79, 0.89, 0.89, "", "brNDC"); 
  legend3.SetBorderSize(0);
  legend3.SetFillColor(0);
  
  legend3.AddEntry(histogramData, "Data", "p");
  legend3.AddEntry(templateSMsum, "MC",   "l");
  legend3.Draw();
  
  cmsPreliminaryLabel.Draw();
  cmsLuminosityLabel.Draw();

  canvas->Update();
  std::string outputFilePath3 = std::string("./plots/");
  if ( sysShift != "CENTRAL_VALUE" ) outputFilePath3.append(sysShift).append("/");
  gSystem->mkdir(outputFilePath3.data(), true);
  TString outputFileName3 = outputFileName;
  outputFileName3.ReplaceAll(".", "_smSum.");
  canvas->Print(outputFilePath3.append(outputFileName3).data());  

  if ( templateZtautau    != histogramZtautau    ) delete templateZtautau;
  if ( templateZmumu      != histogramZmumu      ) delete templateZmumu;
  if ( templateQCD        != histogramQCD        ) delete templateQCD;
  if ( templateWplusJets  != histogramWplusJets  ) delete templateWplusJets;
  if ( templateTTplusJets != histogramTTplusJets ) delete templateTTplusJets;
  delete templateSMbgSum;
  delete templateSMsum;

  delete canvas;
}

void drawHistograms(std::map<std::string, std::map<std::string, TH1*> >& distributionsData, 
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,    
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions,
		    std::map<std::string, RooAbsReal*> normFactors,                                           
		    const std::string& region, const std::string& observable_key,
		    const std::string& histogramTitle, const std::string& xAxisTitle, 
		    const std::string& outputFileName, const std::string& sysShift = "CENTRAL_VALUE")
{
  //std::cout << "<drawHistograms (wrapper)>:" << std::endl;

  drawHistograms(templatesAll["Ztautau"][region][observable_key], 
		 normFactors["Ztautau"]->getVal()*fittedFractions["Ztautau"][region][observable_key],
		 templatesAll["Zmumu"][region][observable_key], 
		 normFactors["Zmumu"]->getVal()*fittedFractions["Zmumu"][region][observable_key],
		 templatesAll["QCD"][region][observable_key], 
		 normFactors["QCD"]->getVal()*fittedFractions["QCD"][region][observable_key],
		 templatesAll["WplusJets"][region][observable_key], 
		 normFactors["WplusJets"]->getVal()*fittedFractions["WplusJets"][region][observable_key],
		 templatesAll["TTplusJets"][region][observable_key], 
		 normFactors["TTplusJets"]->getVal()*fittedFractions["TTplusJets"][region][observable_key],
		 distributionsData[region][observable_key],
		 histogramTitle, xAxisTitle,
		 outputFileName, sysShift);
}

void addToFormula(std::string& formula, const std::string& expression, TObjArray& arguments, RooAbsReal* p)
{
//-------------------------------------------------------------------------------
// Build formula and argument list for RooFormulaVar
//
// NOTE: auxiliary function for makeRooFormulaVar
//
//-------------------------------------------------------------------------------

  //std::cout << "<addToFormula>:" << std::endl;

  if ( expression != "" ) {
    if ( formula != "" ) formula.append(" * ");

    if      ( expression == "regular"  ) formula.append(p->GetName());
    else if ( expression == "inverted" ) formula.append("(1 - ").append(p->GetName()).append(")");
    else assert(0);
    
    arguments.Add(p);
  }
}

RooFormulaVar* makeRooFormulaVar(const std::string& process, const std::string& region,
				 RooAbsReal* norm, double fittedFractionValue,
				 RooAbsReal* pDiTauCharge_OS_SS, RooAbsReal* pDiTauKine_Sig_Bgr, 
				 RooAbsReal* pMuonIso_loose_tight, 
				 RooAbsReal* pTauId_passed_failed)
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

  std::cout << "<makeRooFormulaVar>:" << std::endl;
  std::cout << " building RooFormulaVar expression for process = " << process << ", region = " << region << "." << std::endl;

  RooFormulaVar* retVal = 0;

  std::string exprDiTauCharge = "";
  std::string exprDiTauKine   = "";
  std::string exprMuonIso     = "";
  std::string exprTauId       = "";

  if        ( region == "A"   ) {
    exprDiTauCharge = "regular";
    exprMuonIso     = "regular";
  } else if ( region == "B"   ) {
    exprDiTauCharge = "inverted";
    exprMuonIso     = "regular";
  } else if ( region.find("C") == 0 ) {
    exprDiTauCharge = "regular";
    exprMuonIso     = "inverted"; 
    if      ( region.find("C1") == 0 ) exprDiTauKine = "regular";
    else if ( region.find("C2") == 0 ) exprDiTauKine = "inverted";
    if      ( region == "C1p" || region == "C2p" ) exprTauId = "regular";
    else if ( region == "C1f" || region == "C2f" ) exprTauId = "inverted";
  } else if ( region == "D"   ) {
    exprDiTauCharge = "inverted";
    exprMuonIso     = "inverted";
  } else {
    std::cout << "Error in <makeRooFormulaVar>: undefined region = " << region << " --> skipping !!" << std::endl;
    return retVal;
  }

  std::string fittedFractionName = std::string(norm->GetName()).append("_fittedFraction");
  RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), fittedFractionValue);

  std::string formula = "";
  TObjArray arguments; 
  addToFormula(formula, exprDiTauCharge, arguments, pDiTauCharge_OS_SS);
  addToFormula(formula, exprDiTauKine,   arguments, pDiTauKine_Sig_Bgr);
  addToFormula(formula, exprMuonIso,     arguments, pMuonIso_loose_tight);
  addToFormula(formula, exprTauId,       arguments, pTauId_passed_failed);
  addToFormula(formula, "regular",       arguments, norm);
  addToFormula(formula, "regular",       arguments, fittedFraction);

  std::cout << " formula = " << formula << std::endl;
  //std::cout << " arguments:" << std::endl;
  //for ( int i = 0; i < arguments.GetEntries(); ++i ) {
  //  std::cout << " " << dynamic_cast<TNamed*>(arguments.At(i))->GetName() << std::endl;
  //}

  std::string name = std::string("norm").append(process).append("_").append(region);
  retVal = new RooFormulaVar(name.data(), name.data(), formula.data(), RooArgSet(arguments));

  return retVal;
}

std::vector<std::string> getObservables(const std::string& region, const std::vector<std::string>& fitVariables)
{
  std::vector<std::string> retVal;

  if        ( region == "ABCD" ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
  } else if ( region == "A" ) {
    retVal.push_back(std::string("diTauMt"));
  } else if ( region.find("B")  == 0 ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else if ( region.find("C")  == 0 ) {
    if ( region == "C" || region == "C1" || region == "C2" ) {
      retVal = fitVariables;
      retVal.push_back(std::string("diTauMt"));
      //retVal.push_back(std::string("muonPt"));
      //retVal.push_back(std::string("tauPt"));
    } else if ( region == "C1p" ) {
      retVal = fitVariables;
    } else if ( region == "C1f" ) {
      retVal = fitVariables;
    } else if ( region == "C2p" ) {
      retVal.push_back(std::string("diTauMt"));
    } else if ( region == "C2f" ) {
      retVal.push_back(std::string("diTauMt"));
    }
  } else if ( region == "D" ) {
    retVal = fitVariables;
    retVal.push_back(std::string("diTauMt"));
    //retVal.push_back(std::string("muonPt"));
    //retVal.push_back(std::string("tauPt"));
  } else {
    std::cout << "Error in <getObservables>: undefined region = " << region << " !!" << std::endl;
  }

  return retVal;
}

std::vector<std::string> getTauIdValues(const std::string& region)
{
  std::vector<std::string> retVal;

  if        ( region == "C1p" ) {
    retVal.push_back(std::string("passed"));
  } else if ( region == "C1f" ) {
    retVal.push_back(std::string("failed"));
  } else if ( region == "C2p" ) {
    retVal.push_back(std::string("passed"));
  } else if ( region == "C2f" ) {
    retVal.push_back(std::string("failed"));
  } else {
    retVal.push_back(std::string("all"));
  }

  return retVal;
}

std::map<std::string, TH1*> makeDistributionsInRegion(
  const std::string& process, const std::string& region, double weight,
  TTree* tree, 
  const std::string& tauId, const std::vector<std::string>& fitVariables,
  std::map<std::string, std::string>& branchNames, 
  const std::string& sysShift = "CENTRAL_VALUE",
  bool saveRunLumiSectionEventNumbers = false)
{
//-------------------------------------------------------------------------------
// Make histogram(s) of observables used for determining MC normalization factors
// in different regions
//
// Return value: observable --> histogram mapping
//
// For a definition of the different regions, 
// cf. comments in makeRooFormulaVar function
//-------------------------------------------------------------------------------
  
  std::cout << "<makeDistributionsInRegion>:" << std::endl;
  std::cout << " selecting " << process << " in region " << region << "..." << std::endl;

  std::string treeSelection;
  
  std::string exprMuonIso_loose = std::string(branchNames["muonLooseIsoPtSum04"]).append(" < (0.30*").append(branchNames["muonPt"]).append(")");
  std::string exprMuonIso_tight = std::string(branchNames["muonLooseIsoPtSum04"]).append(" < (0.10*").append(branchNames["muonPt"]).append(")");
  std::string exprMuonIso_loose_not_tight = std::string("(").append(exprMuonIso_loose);
  exprMuonIso_loose_not_tight.append(" && ").append(branchNames["muonLooseIsoPtSum04"]).append(" > (0.10*").append(branchNames["muonPt"]).append(")").append(")");
  std::string exprDiTauCharge_OS = std::string("abs(").append(branchNames["diTauCharge"]).append(") < 0.5");
  std::string exprDiTauCharge_SS = std::string("abs(").append(branchNames["diTauCharge"]).append(") > 1.5");
  std::string exprDiTauCharge_OS_or_SS = std::string("(").append(exprDiTauCharge_OS);
  exprDiTauCharge_OS_or_SS.append(" || ").append(exprDiTauCharge_SS).append(")");
  std::string exprDiTauKine_Sig  = std::string("(").append(branchNames["diTauMt"]).append(" < 40");
  exprDiTauKine_Sig.append(" && ").append("(").append(branchNames["diTauPzeta"]).append(" - 1.5*").append(branchNames["diTauPzetaVis"]).append(") > -20)");
  std::string exprDiTauKine_Bgr  = std::string("(").append(branchNames["diTauMt"]).append(" > 40");
  exprDiTauKine_Bgr.append(" || ").append("(").append(branchNames["diTauPzeta"]).append(" - 1.5*").append(branchNames["diTauPzetaVis"]).append(") < -20)");
    
  if        ( region == "ABCD" ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprDiTauCharge_OS_or_SS);
  } else if ( region == "A" ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprMuonIso_loose_not_tight).append(" && ").append(exprDiTauCharge_OS);
  } else if ( region.find("B")  == 0 ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprMuonIso_loose_not_tight).append(" && ").append(exprDiTauCharge_SS);
    if      ( region.find("B1") == 0 ) treeSelection.append(" && ").append(exprDiTauKine_Sig);
  } else if ( region.find("C")  == 0 ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_OS);
    if      ( region.find("C1") == 0 ) treeSelection.append(" && ").append(exprDiTauKine_Sig);
    else if ( region.find("C2") == 0 ) treeSelection.append(" && ").append(exprDiTauKine_Bgr);
  } else if ( region == "D" ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_SS);
  } else {
    std::cout << "Error in <makeDistribution>: undefined region = " << region << " --> skipping !!" << std::endl;
    return std::map<std::string, TH1*>();
  }

  std::vector<std::string> observables = getObservables(region, fitVariables);
  std::vector<std::string> tauIdValues = getTauIdValues(region);

  return makeHistograms(process, region, weight,
			tree, treeSelection,
			tauId, tauIdValues, observables,
			branchNames, sysShift,
			saveRunLumiSectionEventNumbers);
}

std::map<std::string, std::map<std::string, TH1*> > makeDistributionsAllRegions(
  const std::string& process, double weight,
  TTree* tree, const std::vector<std::string>& regions,
  const std::vector<std::string>& tauIds, const std::vector<std::string>& fitVariables,
  std::map<std::string, std::map<std::string, std::string> >& branchNames, const std::string& sysShift = "CENTRAL_VALUE",
  std::map<std::string, bool>* saveRunLumiSectionEventNumbers = NULL)
{
//-------------------------------------------------------------------------------
// Make histogram(s) of observables used for determining MC normalization factors
// in different regions
//
// Return value: mapping of (region, observable) --> histogram 
//
//   observable = 'Mt'         for regions { A, B, C2p, C2f, D }
//                fitVariables for regions { C1p, C1f }
//
// For a definition of the different regions, 
// cf. comments in makeRooFormulaVar function
//-------------------------------------------------------------------------------
  
  std::cout << "<makeDistributionsAllRegions>:" << std::endl;

  std::map<std::string, std::map<std::string, TH1*> > retVal;
  
  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {	
    for ( std::vector<std::string>::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {

      bool saveRunLumiSectionEventNumbers_region = ( saveRunLumiSectionEventNumbers && tauId == tauIds.begin() ) ?
	(*saveRunLumiSectionEventNumbers)[*region] : false;

      std::map<std::string, TH1*> retVal_region = makeDistributionsInRegion(process, *region, weight,
									    tree, 
									    *tauId, fitVariables,
									    branchNames[*tauId], sysShift,
									    saveRunLumiSectionEventNumbers_region);

      for ( std::map<std::string, TH1*>::iterator distribution = retVal_region.begin();
	    distribution != retVal_region.end(); ++distribution ) {
	retVal[*region][distribution->first] = distribution->second;
	//std::cout << " retVal[" << (*region) << "][" << distribution->first << "] = " << distribution->second << std::endl;
      }
    }
  }

  return retVal;
}

RooHistPdf* makeRooHistPdf(TH1* templateHistogram, RooAbsReal* fitVar)
{
  //std::cout << "<makeRooHistPdf>:" << std::endl;

  std::string templateDataHistName = std::string(templateHistogram->GetName()).append("_dataHist");
  RooDataHist* templateDataHist = new RooDataHist(templateDataHistName.data(), templateDataHistName.data(), *fitVar, templateHistogram);
  std::string templatePdfName = std::string(templateHistogram->GetName()).append("_histPdf");
  RooHistPdf* templatePdf = new RooHistPdf(templatePdfName.data(), templatePdfName.data(), *fitVar, *templateDataHist);
  return templatePdf;
}

RooGaussian* makeFitConstraint(RooAbsReal* p, double value, double error)
{
  //std::cout << "<makeFitConstraint>:" << std::endl;

  std::string pValueName = std::string(p->GetName()).append("_constValue");
  RooConstVar* pValue = new RooConstVar(pValueName.data(), pValueName.data(), value);
  std::string pErrorName = std::string(p->GetName()).append("_constError");
  RooConstVar* pError = new RooConstVar(pErrorName.data(), pErrorName.data(), error);

  std::string constraintName = std::string(p->GetName()).append("_constraint");
  RooGaussian* constraint = new RooGaussian(constraintName.data(), constraintName.data(), *p, *pValue, *pError);
  return constraint;
}

void fitUsingRooFit(std::map<std::string, std::map<std::string, TH1*> >& distributionsData,                         // key = (region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,      // key = (process, region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll,    // key = (process/"sum", region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions, // key = (process, region, observable)
		    const std::vector<std::string>& processes,
		    const std::string& tauId, const std::string& fitVariable, bool fitTauIdEffC2,
		    double& effValue, double& effError,
		    const std::string& sysShift = "CENTRAL_VALUE",
		    std::map<std::string, std::string>* xAxisTitles = NULL)
{
  std::cout << "<fitUsingRooFit>:" << std::endl;
  std::cout << " performing Fit of variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;

  double fitMinABC2D = templatesAll["Ztautau"]["A"][getKey("diTauMt", tauId)]->GetXaxis()->GetXmin();
  double fitMaxABC2D = templatesAll["Ztautau"]["A"][getKey("diTauMt", tauId)]->GetXaxis()->GetXmax();
  RooRealVar* fitVarABC2D = new RooRealVar("fitVarABC2D", "fitVarABC2D", fitMinABC2D, fitMaxABC2D);

  double fitMinC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmin();
  double fitMaxC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmax();
  double fitMinC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmin();
  double fitMaxC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmax();
  assert(fitMinC1p == fitMinC1f && fitMaxC1p == fitMaxC1f);
  RooRealVar* fitVarC1 = new RooRealVar("fitVarC1", "fitVarC1", fitMinC1p, fitMaxC1p);
  
  double numEventsDataABCD = distributionsData["ABCD"][getKey("diTauMt", tauId)]->Integral();
  std::cout << "numEventsDataABCD = " << numEventsDataABCD << std::endl;

  std::map<std::string, RooRealVar*> pDiTauCharge_OS_SS;          // key = process
  std::map<std::string, RooRealVar*> pMuonIso_loose_tight;        // key = process
  std::map<std::string, RooRealVar*> pDiTauKine_Sig_Bgr;          // key = process
  std::map<std::string, RooRealVar*> pTauId_passed_failed;        // key = process

  std::map<std::string, RooAbsReal*> normABCD;                    // key = process
  std::map<std::string, RooAbsReal*> normA;                       // key = process
  std::map<std::string, RooAbsReal*> normB;                       // key = process
  std::map<std::string, RooAbsReal*> normC1;                      // key = process
  std::map<std::string, RooAbsReal*> normC1p;                     // key = process
  std::map<std::string, RooAbsReal*> normC1f;                     // key = process
  std::map<std::string, RooAbsReal*> normC2;                      // key = process
  std::map<std::string, RooAbsReal*> normC2p;                     // key = process
  std::map<std::string, RooAbsReal*> normC2f;                     // key = process
  std::map<std::string, RooAbsReal*> normD;                       // key = process
  
  std::map<std::string, std::map<std::string, RooHistPdf*> > pdf; // key = (process/"sum", region)

  TObjArray pdfsA;
  TObjArray pdfsB;
  TObjArray pdfsC1p;
  TObjArray pdfsC1f;
  TObjArray pdfsC2p;
  TObjArray pdfsC2f;
  TObjArray pdfsC2;
  TObjArray pdfsD;

  TObjArray fitParametersA;
  TObjArray fitParametersB;
  TObjArray fitParametersC1p;
  TObjArray fitParametersC1f;
  TObjArray fitParametersC2p;
  TObjArray fitParametersC2f;
  TObjArray fitParametersC2;
  TObjArray fitParametersD;

  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    std::cout << "process = " << (*process) << ":" << std::endl;

    double numEventsA         = numEventsAll[*process]["A"][getKey("diTauMt", tauId)];
    double fittedFractionA    = fittedFractions[*process]["A"][getKey("diTauMt", tauId)];
    double fittedEventsA      = fittedFractionA*numEventsA;
    double numEventsB         = numEventsAll[*process]["B"][getKey("diTauMt", tauId)];    
    double fittedFractionB    = fittedFractions[*process]["B"][getKey("diTauMt", tauId)];
    double fittedEventsB      = fittedFractionB*numEventsB;
    double numEventsC         = numEventsAll[*process]["C"][getKey("diTauMt", tauId)];
    double fittedFractionC    = fittedFractions[*process]["C"][getKey("diTauMt", tauId)];
    double fittedEventsC      = fittedFractionC*numEventsC;
    double numEventsC1        = numEventsAll[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedFractionC1   = fittedFractions[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedEventsC1     = fittedFractionC1*numEventsC1;
    double numEventsC1p       = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedFractionC1p  = fittedFractions[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedEventsC1p    = fittedFractionC1p*numEventsC1p;
    double numEventsC1f       = numEventsAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedFractionC1f  = fittedFractions[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedEventsC1f    = fittedFractionC1f*numEventsC1f;
    double numEventsC2p       = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedFractionC2p  = fittedFractions[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedEventsC2p    = fittedFractionC2p*numEventsC2p;
    double numEventsC2f       = numEventsAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    double fittedFractionC2f  = fittedFractions[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    double fittedEventsC2f    = fittedFractionC2f*numEventsC2f;
    double numEventsC2        = numEventsAll[*process]["C2"][getKey("diTauMt", tauId)];
    double fittedFractionC2   = fittedFractions[*process]["C2"][getKey("diTauMt", tauId)];
    double fittedEventsC2     = fittedFractionC2*numEventsC2;
    double numEventsD         = numEventsAll[*process]["D"][getKey("diTauMt", tauId)];
    double fittedFractionD    = fittedFractions[*process]["D"][getKey("diTauMt", tauId)];
    double fittedEventsD      = fittedFractionD*numEventsD;
    double numEventsABCD      = numEventsAll[*process]["ABCD"][getKey("diTauMt", tauId)];
    double fittedEventsABCD   = fittedEventsA + fittedEventsB + fittedEventsC + fittedEventsD;
    double fittedFractionABCD = fittedEventsABCD/numEventsABCD;
    std::cout << " numEventsABCD = " << numEventsABCD << ", fittedFractionABCD = " << fittedFractionABCD << std::endl;

    std::string nameDiTauCharge_OS_SS   = std::string("pDiTauCharge_OS_SS").append("_").append(*process);
    double pDiTauCharge_OS_SS0          = (fittedEventsA + fittedEventsC)/fittedEventsABCD;
    pDiTauCharge_OS_SS[*process]        = new RooRealVar(nameDiTauCharge_OS_SS.data(), 							 
							 nameDiTauCharge_OS_SS.data(), pDiTauCharge_OS_SS0, 0., 1.);
    std::string nameMuonIso_loose_tight = std::string("pMuonIso_loose_tight").append("_").append(*process);
    double pMuonIso_loose_tight0        = (fittedEventsA + fittedEventsB)/fittedEventsABCD;
    pMuonIso_loose_tight[*process]      = new RooRealVar(nameMuonIso_loose_tight.data(), 
							 nameDiTauCharge_OS_SS.data(), pMuonIso_loose_tight0, 0., 1.);
    std::string nameDiTauKine_Sig_Bgr   = std::string("pDiTauKine_Sig_Bgr").append("_").append(*process);
    double pDiTauKine_Sig_Bgr0          = fittedEventsC1/fittedEventsC;
    pDiTauKine_Sig_Bgr[*process]        = new RooRealVar(nameDiTauKine_Sig_Bgr.data(), 
							 nameDiTauKine_Sig_Bgr.data(), pDiTauKine_Sig_Bgr0, 0., 1.);
    std::string nameTauId_passed_failed = std::string("pTauId_passed_failed").append("_").append(*process);
    double pTauId_passed_failed0        = ( fitTauIdEffC2 ) ? 
      (fittedEventsC1p + fittedEventsC2p)/fittedEventsC : fittedEventsC1p/fittedEventsC1;
    pTauId_passed_failed[*process]      = new RooRealVar(nameTauId_passed_failed.data(),
							 nameTauId_passed_failed.data(), pTauId_passed_failed0, 0., 1.);

    double numEventsSumABCD = numEventsAll["sum"]["ABCD"][getKey("diTauMt", tauId)];
    std::cout << " numEventsSumABCD = " << numEventsSumABCD << std::endl;

    //double scaleFactorMCtoData = numEventsDataABCD/numEventsSumABCD;
    double scaleFactorMCtoData = 1.;
    std::cout << "--> MC-to-Data scale-factor = " << scaleFactorMCtoData << std::endl;

    std::string nameNormABCD = std::string("normABCD").append("_").append(*process);
    double normABCD0         = scaleFactorMCtoData*numEventsABCD;
    normABCD[*process]       = new RooRealVar(nameNormABCD.data(), nameNormABCD.data(), normABCD0, 0., numEventsDataABCD);

    TH1* templateA     = templatesAll[*process]["A"][getKey("diTauMt", tauId)];
    RooHistPdf* pdfA   = makeRooHistPdf(templateA, fitVarABC2D);
    TH1* templateB     = templatesAll[*process]["B"][getKey("diTauMt", tauId)];
    RooHistPdf* pdfB   = makeRooHistPdf(templateB, fitVarABC2D);
    TH1* templateC1p   = templatesAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    RooHistPdf* pdfC1p = makeRooHistPdf(templateC1p, fitVarC1);
    TH1* templateC1f   = templatesAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    RooHistPdf* pdfC1f = makeRooHistPdf(templateC1f, fitVarC1);    
    TH1* templateC2p   = templatesAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    RooHistPdf* pdfC2p = makeRooHistPdf(templateC2p, fitVarABC2D);
    TH1* templateC2f   = templatesAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    RooHistPdf* pdfC2f = makeRooHistPdf(templateC2f, fitVarABC2D);
    TH1* templateC2    = templatesAll[*process]["C2"][getKey("diTauMt", tauId)];
    RooHistPdf* pdfC2  = makeRooHistPdf(templateC2, fitVarABC2D);
    TH1* templateD     = templatesAll[*process]["D"][getKey("diTauMt", tauId)];
    RooHistPdf* pdfD   = makeRooHistPdf(templateD, fitVarABC2D);

    pdfsA.Add(pdfA);
    pdfsB.Add(pdfB);
    pdfsC1p.Add(pdfC1p);
    pdfsC1f.Add(pdfC1f);
    pdfsC2p.Add(pdfC2p);
    pdfsC2f.Add(pdfC2f);
    pdfsC2.Add(pdfC2);
    pdfsD.Add(pdfD);

    normA[*process]   = makeRooFormulaVar(*process, "A", normABCD[*process], fittedFractionA, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normB[*process]   = makeRooFormulaVar(*process, "B", normABCD[*process], fittedFractionB, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC1[*process]  = makeRooFormulaVar(*process, "C1", normABCD[*process], fittedFractionC1,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC1p[*process] = makeRooFormulaVar(*process, "C1p", normABCD[*process], fittedFractionC1p, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC1f[*process] = makeRooFormulaVar(*process, "C1f", normABCD[*process], fittedFractionC1f, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC2[*process]  = makeRooFormulaVar(*process, "C2", normABCD[*process], fittedFractionC2,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC2p[*process] = makeRooFormulaVar(*process, "C2p", normABCD[*process], fittedFractionC2p,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC2f[*process] = makeRooFormulaVar(*process, "C2f", normABCD[*process], fittedFractionC2f,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);    
    normD[*process]   = makeRooFormulaVar(*process, "D", normABCD[*process], fittedFractionD,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);

    fitParametersA.Add(normA[*process]);
    fitParametersB.Add(normB[*process]);
    fitParametersC1p.Add(normC1p[*process]);
    fitParametersC1f.Add(normC1f[*process]);
    fitParametersC2p.Add(normC2p[*process]);
    fitParametersC2f.Add(normC2f[*process]);
    fitParametersC2.Add(normC2[*process]);
    fitParametersD.Add(normD[*process]);
  }

//--- CV: Monte Carlo simulation underestimates probability pMuonIso_loose_tight0  
//        for muons to be isolated in QCD background events
//       --> adjust start value, to improve convergence of fit
//       (correction factor 1.4 determined "by eye")
  pMuonIso_loose_tight["QCD"]->setVal(1.4*pMuonIso_loose_tight["QCD"]->getVal());

  RooAddPdf* pdfSumA   = new RooAddPdf("pdfSumA",   "pdfSumB",   RooArgList(pdfsA),   RooArgList(fitParametersA));
  RooAddPdf* pdfSumB   = new RooAddPdf("pdfSumB",   "pdfSumB",   RooArgList(pdfsB),   RooArgList(fitParametersB));
  RooAddPdf* pdfSumC1p = new RooAddPdf("pdfSumC1p", "pdfSumC1p", RooArgList(pdfsC1p), RooArgList(fitParametersC1p));
  RooAddPdf* pdfSumC1f = new RooAddPdf("pdfSumC1f", "pdfSumC1f", RooArgList(pdfsC1f), RooArgList(fitParametersC1f));
  RooAddPdf* pdfSumC2p = new RooAddPdf("pdfSumC2p", "pdfSumC2p", RooArgList(pdfsC2p), RooArgList(fitParametersC2p));
  RooAddPdf* pdfSumC2f = new RooAddPdf("pdfSumC2f", "pdfSumC2f", RooArgList(pdfsC2f), RooArgList(fitParametersC2f));
  RooAddPdf* pdfSumC2  = new RooAddPdf("pdfSumC2",  "pdfSumC2",  RooArgList(pdfsC2),  RooArgList(fitParametersC2));
  RooAddPdf* pdfSumD   = new RooAddPdf("pdfSumD",   "pdfSumD",   RooArgList(pdfsD),   RooArgList(fitParametersD));

// CV: due to limitation in RooFit
//    (cf. http://root.cern.ch/phpBB3/viewtopic.php?f=15&t=9518)
//     need to construct log-likelihood functions separately for regions { A, B, D } and { C1p, C1f }

//--- build data & model objects for fitting regions A, B, C2p, C2f, D
  RooCategory* fitCategoriesABC2D = new RooCategory("categoriesABC2D", "categoriesABC2D");
  fitCategoriesABC2D->defineType("A");
  fitCategoriesABC2D->defineType("B");
  if ( fitTauIdEffC2 ) {
    fitCategoriesABC2D->defineType("C2p");
    fitCategoriesABC2D->defineType("C2f");
  } else {
    fitCategoriesABC2D->defineType("C2");
  }
  fitCategoriesABC2D->defineType("D");

  RooSimultaneous* pdfSimultaneousFitABC2D = new RooSimultaneous("pdfSimultaneousFitABC2D", "pdfSimultaneousFitABC2D", *fitCategoriesABC2D);
  pdfSimultaneousFitABC2D->addPdf(*pdfSumA,   "A");
  pdfSimultaneousFitABC2D->addPdf(*pdfSumB,   "B");
  if ( fitTauIdEffC2 ) {
    pdfSimultaneousFitABC2D->addPdf(*pdfSumC2p, "C2p");
    pdfSimultaneousFitABC2D->addPdf(*pdfSumC2f, "C2f");
  } else {
    pdfSimultaneousFitABC2D->addPdf(*pdfSumC2,  "C2");
  }
  pdfSimultaneousFitABC2D->addPdf(*pdfSumD,   "D");

  std::map<std::string, TH1*> histogramDataMapABC2D;
  histogramDataMapABC2D["A"]     = distributionsData["A"][getKey("diTauMt", tauId)];
  histogramDataMapABC2D["B"]     = distributionsData["B"][getKey("diTauMt", tauId)];
  if ( fitTauIdEffC2 ) {
    histogramDataMapABC2D["C2p"] = distributionsData["C2p"][getKey("diTauMt", tauId, "passed")];
    histogramDataMapABC2D["C2f"] = distributionsData["C2f"][getKey("diTauMt", tauId, "failed")];
  } else {
    histogramDataMapABC2D["C2"]  = distributionsData["C2"][getKey("diTauMt", tauId)];
  }
  histogramDataMapABC2D["D"]     = distributionsData["D"][getKey("diTauMt", tauId)];

  RooDataHist* dataABC2D = new RooDataHist("dataABC2D", "data", *fitVarABC2D, *fitCategoriesABC2D, histogramDataMapABC2D);

//--- build data & model objects for fitting regions C1p, C1f
  RooCategory* fitCategoriesC1 = new RooCategory("categoriesC1", "categoriesC1");
  fitCategoriesC1->defineType("C1p");
  fitCategoriesC1->defineType("C1f");

  RooSimultaneous* pdfSimultaneousFitC1 = new RooSimultaneous("pdfSimultaneousFitC1", "pdfSimultaneousFitC1", *fitCategoriesC1);
  pdfSimultaneousFitC1->addPdf(*pdfSumC1p, "C1p");
  pdfSimultaneousFitC1->addPdf(*pdfSumC1f, "C1f");
 
  std::map<std::string, TH1*> histogramDataMapC1;
  histogramDataMapC1["C1p"] = distributionsData["C1p"][getKey(fitVariable, tauId, "passed")];
  histogramDataMapC1["C1f"] = distributionsData["C1f"][getKey(fitVariable, tauId, "failed")];

  RooDataHist* dataC1 = new RooDataHist("dataC1", "dataC1", *fitVarC1, *fitCategoriesC1, histogramDataMapC1);

//--- add "external" constraints
//    on probabilities 
//   o pDiTauCharge_OS_SS
//   o pDiTauKine_Sig_Bgr
//   o pMuonIso_loose_tight
//    separating different regions

  pDiTauCharge_OS_SS["Ztautau"]->setConstant(true);
  pMuonIso_loose_tight["Ztautau"]->setConstant(true);
  pDiTauKine_Sig_Bgr["Ztautau"]->setConstant(true);

  pDiTauCharge_OS_SS["Zmumu"]->setConstant(true);
  pMuonIso_loose_tight["Zmumu"]->setConstant(true);
  pDiTauKine_Sig_Bgr["Zmumu"]->setConstant(true);

  pDiTauCharge_OS_SS["TTplusJets"]->setConstant(true);
  pMuonIso_loose_tight["TTplusJets"]->setConstant(true);
  pDiTauKine_Sig_Bgr["TTplusJets"]->setConstant(true);

//--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  TObjArray fitConstraintsC1;
  fitConstraintsC1.Add(makeFitConstraint(normABCD["Zmumu"],      
					 normABCD["Zmumu"]->getVal(),                 0.5*normABCD["Zmumu"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normABCD["QCD"],        
					 normABCD["QCD"]->getVal(),                   0.5*normABCD["QCD"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normABCD["WplusJets"],  
					 normABCD["WplusJets"]->getVal(),             0.5*normABCD["WplusJets"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normABCD["TTplusJets"], 
					 normABCD["TTplusJets"]->getVal(),            0.5*normABCD["TTplusJets"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["Zmumu"],      
  //					   pTauId_passed_failed["Zmumu"]->getVal(),        1.0*pTauId_passed_failed["Zmumu"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["QCD"],        
  //					   pTauId_passed_failed["QCD"]->getVal(),          1.0*pTauId_passed_failed["QCD"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["WplusJets"],  
  //					   pTauId_passed_failed["WplusJets"]->getVal(),    1.0*pTauId_passed_failed["WplusJets"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["TTplusJets"], 
  //					   pTauId_passed_failed["TTplusJets"]->getVal(),   1.0*pTauId_passed_failed["TTplusJets"]->getVal()));

  RooLinkedList fitOptionsC1;
  fitOptionsC1.Add(new RooCmdArg(RooFit::Extended()));
  fitOptionsC1.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsC1))));

  TObjArray fitConstraintsABC2D;
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["QCD"],         
  //					      pDiTauCharge_OS_SS["QCD"]->getVal(),         0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["QCD"],         
  //					      pDiTauKine_Sig_Bgr["QCD"]->getVal(),         0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_loose_tight["QCD"],       
  //					      pMuonIso_loose_tight["QCD"]->getVal(),       0.25));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["WplusJets"],   
  //					      pDiTauCharge_OS_SS["WplusJets"]->getVal(),   0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["WplusJets"],   
  //					      pDiTauKine_Sig_Bgr["WplusJets"]->getVal(),   0.2));
  //fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_loose_tight["WplusJets"], 
  //					      pMuonIso_loose_tight["WplusJets"]->getVal(), 0.1));

  RooLinkedList fitOptionsABC2D;
  fitOptionsABC2D.Add(new RooCmdArg(RooFit::Extended()));
  //fitOptionsABC2D.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsABC2D))));
  
/*
  fitOptionsC1.Add(new RooCmdArg(RooFit::Save(true)));
  RooFitResult*	fitResult = pdfSimultaneousFitC1->fitTo(*dataC1, fitOptionsC1);
 */
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
  std::string fitResultName = std::string("fitResult").append("_").append(tauId);
  RooFitResult*	fitResult = minuit.save(fitResultName.data(), fitResultName.data());
   
  std::cout << tauId << ":";
  if ( fitResult->status() == 0 ) std::cout << " fit converged."          << std::endl; 
  else                            std::cout << " fit failed to converge." << std::endl;

  const RooArgList& fitParameter = fitResult->floatParsFinal();

  int numFitParameter = fitParameter.getSize();
  
  TMatrixD cov(numFitParameter, numFitParameter);
  for ( int iParameter = 0; iParameter < numFitParameter; ++iParameter ) {
    const RooAbsArg* paramI_arg = fitParameter.at(iParameter);
    const RooRealVar* paramI = dynamic_cast<const RooRealVar*>(paramI_arg);    
    double sigmaI = paramI->getError();

    std::cout << " parameter #" << iParameter << ": " << paramI_arg->GetName() 
	      << " = " << paramI->getVal() << " +/- " << paramI->getError() << std::endl;
    
    for ( int jParameter = 0; jParameter < numFitParameter; ++jParameter ) {
      const RooAbsArg* paramJ_arg = fitParameter.at(jParameter);
      const RooRealVar* paramJ = dynamic_cast<const RooRealVar*>(paramJ_arg);
      double sigmaJ = paramJ->getError();

      double corrIJ = fitResult->correlation(*paramI_arg, *paramJ_arg);

      cov(iParameter, jParameter) = sigmaI*sigmaJ*corrIJ;
    }
  }

  cov.Print();

  std::cout << std::endl;

  std::cout << "Results of fitting variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    double numEventsA         = numEventsAll[*process]["A"][getKey("diTauMt", tauId)];  
    double fittedFractionA    = fittedFractions[*process]["A"][getKey("diTauMt", tauId)];
    double fittedEventsA      = fittedFractionA*numEventsA;
    double numEventsB         = numEventsAll[*process]["B"][getKey("diTauMt", tauId)];    
    double fittedFractionB    = fittedFractions[*process]["B"][getKey("diTauMt", tauId)];
    double fittedEventsB      = fittedFractionB*numEventsB;
    double numEventsC         = numEventsAll[*process]["C"][getKey("diTauMt", tauId)];
    double fittedFractionC    = fittedFractions[*process]["C"][getKey("diTauMt", tauId)];
    double fittedEventsC      = fittedFractionC*numEventsC;
    double numEventsC1        = numEventsAll[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedFractionC1   = fittedFractions[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedEventsC1     = fittedFractionC1*numEventsC1;
    double numEventsC1p       = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedFractionC1p  = fittedFractions[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedEventsC1p    = fittedFractionC1p*numEventsC1p;
    double numEventsC1f       = numEventsAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedFractionC1f  = fittedFractions[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedEventsC1f    = fittedFractionC1f*numEventsC1f;
    double numEventsC2p       = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedFractionC2p  = fittedFractions[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedEventsC2p    = fittedFractionC2p*numEventsC2p;
    double numEventsC2f       = numEventsAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    double fittedFractionC2f  = fittedFractions[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    double fittedEventsC2f    = fittedFractionC2f*numEventsC2f;
    double numEventsC2        = numEventsAll[*process]["C2"][getKey("diTauMt", tauId)];
    double fittedFractionC2   = fittedFractions[*process]["C2"][getKey("diTauMt", tauId)];
    double fittedEventsC2     = fittedFractionC2*numEventsC2;
    double numEventsD         = numEventsAll[*process]["D"][getKey("diTauMt", tauId)];
    double fittedFractionD    = fittedFractions[*process]["D"][getKey("diTauMt", tauId)];
    double fittedEventsD      = fittedFractionD*numEventsD;
    double numEventsABCD      = numEventsAll[*process]["ABCD"][getKey("diTauMt", tauId)];
    double fittedEventsABCD   = fittedEventsA + fittedEventsB + fittedEventsC + fittedEventsD;
    double fittedFractionABCD = fittedEventsABCD/numEventsABCD;
    
    std::cout << " " << (*process) << ":" << std::endl;
    std::cout << "  normalization = " << normABCD[*process]->getVal() 
	      << " +/- " << dynamic_cast<RooRealVar*>(normABCD[*process])->getError() 
	      << " (MC exp. = " << numEventsABCD << ")" << std::endl;
    std::cout << "  pDiTauCharge_OS_SS = " << pDiTauCharge_OS_SS[*process]->getVal() 
	      << " +/- " << pDiTauCharge_OS_SS[*process]->getError() 
	      << " (MC exp. = " << (fittedEventsA + fittedEventsC)/fittedEventsABCD << ")" << std::endl;
    std::cout << "  pMuonIso_loose_tight = " << pMuonIso_loose_tight[*process]->getVal() 
	      << " +/- " << pMuonIso_loose_tight[*process]->getError() 
	      << " (MC exp. = " << (fittedEventsA + fittedEventsB)/fittedEventsABCD << ")" << std::endl;
    std::cout << "  pDiTauKine_Sig_Bgr = " << pDiTauKine_Sig_Bgr[*process]->getVal() 
	      << " +/- " << pDiTauKine_Sig_Bgr[*process]->getError() 
	      << " (MC exp. = " << fittedEventsC1/fittedEventsC << ")" << std::endl;
    double pTauId_passed_failedMCexp = ( fitTauIdEffC2 ) ? 
      (fittedEventsC1p + fittedEventsC2p)/fittedEventsC : fittedEventsC1p/fittedEventsC1;
    std::cout << "  pTauId_passed_failed = " << pTauId_passed_failed[*process]->getVal() 
	      << " +/- " << pTauId_passed_failed[*process]->getError() 
	      << " (MC exp. = " << pTauId_passed_failedMCexp << ")" << std::endl;
    std::cout << "--> A = " << normA[*process]->getVal() 
	      << " (MC exp. = " << numEventsA << ")" << std::endl;
    std::cout << "--> B = " << normB[*process]->getVal()  
	      << " (MC exp. = " << numEventsB << ")" << std::endl;
    std::cout << "--> C = " << normC1p[*process]->getVal() + normC1f[*process]->getVal() + normC2[*process]->getVal()
	      << " (MC exp. = " << numEventsC << ")" << std::endl;
    std::cout << "--> C1 = " << normC1p[*process]->getVal() + normC1f[*process]->getVal()  
	      << " (MC exp. = " << numEventsC1 << ")" << std::endl;
    std::cout << "--> C1p = " << normC1p[*process]->getVal()  
	      << " (MC exp. = " << numEventsC1p << ")" << std::endl;
    std::cout << "--> C1f = " << normC1f[*process]->getVal()
	      << " (MC exp. = " << numEventsC1f << ")" << std::endl;
    std::cout << "--> C2 = " << normC2[*process]->getVal() 
	      << " (MC exp. = " << numEventsC2 << ")" << std::endl;
    if ( fitTauIdEffC2 ) {
      std::cout << "--> C2p = " << normC2p[*process]->getVal() 
    	        << " (MC exp. = " << numEventsC2p << ")" << std::endl;
      std::cout << "--> C2f = " << normC2f[*process]->getVal()  
    	        << " (MC exp. = " << numEventsC2f << ")" << std::endl;
    }
    std::cout << "--> D = " << normD[*process]->getVal()  
	      << " (MC exp. = " << numEventsD << ")" << std::endl;
  }

  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();

//--- make control plots for sum(MC) scaled by normalization determined by fit versus Data 
//    for Mt, fitVariable distributions in different regions
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normABCD, "ABCD", getKey(fitVariable, tauId),
		 std::string("All Events: ").append(fitVariable).append(" (scaled by normalization det. by fit)"),
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff_ABCD_").append(fitVariable).append("_fitted.pdf"), sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normABCD, "ABCD", getKey("diTauMt", tauId),
		 "All Events: M_{T} (scaled by normalization det. by fit)", xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		 "controlPlotsTauIdEff_ABCD_Mt_fitted.pdf", sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normA, "A", getKey("diTauMt", tauId),
		 "Region A: M_{T} (scaled by normalization det. by fit)", xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		 "controlPlotsTauIdEff_A_Mt_fitted.pdf", sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normB, "B", getKey("diTauMt", tauId),
		 "Region B: M_{T} (scaled by normalization det. by fit)", xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		 "controlPlotsTauIdEff_B_Mt_fitted.pdf", sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normC1, "C1", getKey(fitVariable, tauId),		 
		 std::string("Region C1: ").append(fitVariable).append(" (scaled by normalization det. by fit)"), 
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff_C1_").append(fitVariable).append("_fitted.pdf"), sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normC1p, "C1p", getKey(fitVariable, tauId, "passed"),
		 std::string("Region C1p: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff_C1p_").append(fitVariable).append("_").append(tauId).append("_fitted.pdf"), sysShift);
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normC1f, "C1f", getKey(fitVariable, tauId, "failed"),
		 std::string("Region C1f: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff_C1f_").append(fitVariable).append("_").append(tauId).append("_fitted.pdf"), sysShift);
  if ( fitTauIdEffC2 ) {
    drawHistograms(distributionsData, templatesAll, fittedFractions, 
		   normC2p, "C2p", getKey("diTauMt", tauId, "passed"),
		   std::string("Region C2p: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		   xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		   std::string("controlPlotsTauIdEff_C2p_").append(fitVariable).append("_").append(tauId).append("_fitted.pdf"), sysShift);
    drawHistograms(distributionsData, templatesAll, fittedFractions, 
		   normC2f, "C2f", getKey("diTauMt", tauId, "failed"),
		   std::string("Region C2f: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		   xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		   std::string("controlPlotsTauIdEff_C2f_").append(fitVariable).append("_").append(tauId).append("_fitted.pdf"), sysShift);
  } else {
    drawHistograms(distributionsData, templatesAll, fittedFractions, 
		   normC2, "C2", getKey("diTauMt", tauId),
		   "Region C2: M_{T} (scaled by normalization det. by fit)", xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		   "controlPlotsTauIdEff_C2_Mt_fitted.pdf", sysShift);
  }
  drawHistograms(distributionsData, templatesAll, fittedFractions, 
		 normD, "D", getKey("diTauMt", tauId),
		 "Region D: M_{T} (scaled by normalization det. by fit)", xAxisTitles ? (*xAxisTitles)["diTauMt"] : "",
		 "controlPlotsTauIdEff_D_Mt_fitted.pdf", sysShift);
}

void addFileNames(TChain* chain, const std::string& inputFilePath, const std::string& sampleName, const std::string& jobId)
{
//--------------------------------------------------------------------------------
// Compose full fileName from inputFilePath, sampleName, jobId
// and add files to TChain given as function argument
//--------------------------------------------------------------------------------

  //std::cout << "<addFileNames>:" << std::endl;

  std::string fileNames = std::string(inputFilePath);
  fileNames.append("tauIdEffMeasEDNtuple_").append(sampleName).append("_").append(jobId).append("*.root");
  
  chain->Add(fileNames.data());
}

void printFileInfo(TChain* chain, const std::string& chainName)
{
  //std::cout << "<printFileInfo>:" << std::endl;

  TObjArray* files = chain->GetListOfFiles();

  int numFiles = files->GetEntries();

  std::cout << " " << chainName << " has " << numFiles << " files:" << std::endl;

  int numEntriesTotal = 0;

  for ( int iFile = 0; iFile < numFiles; ++iFile ) {
    TNamed* fileName = dynamic_cast<TNamed*>(files->At(iFile));
    TFile* file = new TFile(fileName->GetTitle());

    TTree* tree = dynamic_cast<TTree*>(file->Get(chain->GetName()));

    int numEntries = tree->GetEntries();

    std::cout << " " << fileName->GetTitle() << " (" << numEntries << " entries)" << std::endl;

    numEntriesTotal += numEntries;

    delete file;
  }

  std::cout << "--> " << numEntriesTotal << " entries in total." << std::endl;
  
}

std::map<std::string, std::map<std::string, TH1*> > loadHistograms(
  TFile* inputFile, const std::string& process, const std::vector<std::string>& regions,
  const std::vector<std::string>& tauIds, const std::vector<std::string>& fitVariables, const std::string& sysShift)
{
//--------------------------------------------------------------------------------
// Load template histograms/distributions observed in data from ROOT file
//--------------------------------------------------------------------------------

  std::map<std::string, std::map<std::string, TH1*> > retVal;

  for ( std::vector<std::string>::const_iterator region = regions.begin();
	region != regions.end(); ++region ) {
    
    std::vector<std::string> observables = getObservables(*region, fitVariables);
    std::vector<std::string> tauIdValues = getTauIdValues(*region);
    
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {
	for ( std::vector<std::string>::const_iterator observable = observables.begin();
	      observable != observables.end(); ++observable ) {

//--- obtain template for QCD background from data,
//    from SS && Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV sideband
	  std::string histogramName;
	  //if ( process == "QCD" && ((*region) == "C1" || (*region) == "C1p" || (*region) == "C1f") ) {
	  //  histogramName = std::string("Data").append("_").append("B1").append("_").append(*observable);
	  //  histogramName.append("_").append(*tauId).append("_").append("all");
	  //  if ( sysShift != "CENTRAL_VALUE" ) histogramName.append("_").append(sysShift);
	  //} else {
	    histogramName = std::string(process).append("_").append(*region).append("_").append(*observable);
	    histogramName.append("_").append(*tauId).append("_").append(*tauIdValue);
	    if ( sysShift != "CENTRAL_VALUE" ) histogramName.append("_").append(sysShift);
	  //}
 
	  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
	  if ( !histogram ) {
	    std::cout << "Error in <loadHistograms>: failed to load histogram = " << histogramName 
		      << " from file = " << inputFile->GetName() << " --> aborting !!";
	    assert(0);
	  }

	  std::string key = getKey(*observable, *tauId, *tauIdValue);	
	  if ( histogram != 0 ) retVal[*region][key] = histogram;
	}
      }
    }
  }

  return retVal;
}

void saveHistograms(TFile* outputFile, std::map<std::string, std::map<std::string, TH1*> >& histograms)
{
//--------------------------------------------------------------------------------
// Write template histograms/distributions observed in data into ROOT file
//--------------------------------------------------------------------------------

  for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = histograms.begin();
	region != histograms.end(); ++region ) {
    for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	  key != region->second.end(); ++key ) {
      TH1* histogram = key->second;
      histogram->Write();
    }
  }
}

void fitTauIdEff_wConstraints()
{
  std::cout << "<fitTauIdEff_wConstraints>:" << std::endl;

  gROOT->SetBatch(true);

//--- keep track of time it took the macro to execute
  TBenchmark clock;
  clock.Start("fitTauIdEff_wConstraints");

  std::string inputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas/2011Feb03b/";
  inputFilePath.append("user/v/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas/");

  const std::string jobId = "2011Feb03bV2";

  //const std::string branchName_suffix = "local";
  const std::string branchName_suffix = "lxbatch";

  bool runClosureTest = false;
  //bool runClosureTest = true;

  bool runQuickTest = false;
  //bool runQuickTest = true;

  //bool fitTauIdEffC2 = false;
  bool fitTauIdEffC2 = true;

  //const std::string histogramFileName = "fitTauIdEff_wConstraints_mcClosure_2011Feb10.root";
  //const std::string histogramFileName = "fitTauIdEff_wConstraints_data_2011Feb10fixedQCD.root";
  const std::string histogramFileName = "fitTauIdEff_wConstraints_data_2011Feb18.root";

  bool loadHistogramsFromFile = false;
  //bool loadHistogramsFromFile = true;

  bool saveHistogramsToFile = (!loadHistogramsFromFile);

  std::vector<std::string> sysUncertainties;
  sysUncertainties.push_back(std::string("SysTauJetEnUp"));
  sysUncertainties.push_back(std::string("SysTauJetEnDown"));
  sysUncertainties.push_back(std::string("SysJetEnUp"));
  sysUncertainties.push_back(std::string("SysJetEnDown"));
  sysUncertainties.push_back(std::string("SysZllRecoilCorrectionUp"));
  sysUncertainties.push_back(std::string("SysZllRecoilCorrectionDown"));
  bool runSysUncertainties = false;
  //bool runSysUncertainties = true;

  std::vector<std::string> regions;
  regions.push_back(std::string("ABCD"));
  regions.push_back(std::string("A"));
  regions.push_back(std::string("B"));
  regions.push_back(std::string("B1"));
  regions.push_back(std::string("C"));
  regions.push_back(std::string("C1"));
  regions.push_back(std::string("C1p"));
  regions.push_back(std::string("C1f"));
  regions.push_back(std::string("C2"));
  regions.push_back(std::string("C2p"));
  regions.push_back(std::string("C2f"));
  regions.push_back(std::string("D"));

  std::map<std::string, bool> saveRunLumiSectionEventNumbers;
  saveRunLumiSectionEventNumbers["A"] = false;
  saveRunLumiSectionEventNumbers["B"] = false;
  saveRunLumiSectionEventNumbers["C"] = false;
  saveRunLumiSectionEventNumbers["C1"] = false;
  saveRunLumiSectionEventNumbers["C1p"] = true;
  saveRunLumiSectionEventNumbers["C1f"] = false;
  saveRunLumiSectionEventNumbers["C2"] = false;
  saveRunLumiSectionEventNumbers["C2p"] = false;
  saveRunLumiSectionEventNumbers["C2f"] = false;
  saveRunLumiSectionEventNumbers["D"] = false;

  double corrFactorData = 0.965;
  double dataIntLumi = 36.2*corrFactorData;

  std::string sampleZtautau = "ZtautauPU156bx";
  double weightFactorZtautau = dataIntLumi*1666/2568490;       // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  double corrFactorZtautau = 1.000*1.000;                      // first  number: correction for event looses during skimming
                                                               // second number: correction for event looses during Ntuple production
  std::string sampleZmumu = "Zmumu_pythia";
  double weightFactorZmumu = dataIntLumi*1666/1998931;         // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  double corrFactorZmumu = 1.000*1.000;
  std::string sampleQCD = "PPmuXptGt20Mu15";
  double weightFactorQCD = dataIntLumi*0.2966*1.e+9*2.855e-4/29504866; // xSection (LO) / numEvents (PYTHIA PPmuXptGt20Mu15 sample)
  double corrFactorQCD = 1.101*1.000;
  //std::string sampleQCD = "PPmuXptGt20Mu10";
  //double weightFactorQCD = dataIntLumi*0.2966*1.e+9*1.18e-3/8063288; // xSection (LO) / numEvents (PYTHIA PPmuXptGt20Mu10 sample)
  //double corrFactorQCD = 1.000*1.000; 
  std::string sampleWplusJets = "WplusJets_madgraph";
  double weightFactorWplusJets = dataIntLumi*31314/15168266;   // W --> l nu xSection (FEWZ @ NNLO) / numEvents (MadGraph sample)
  double corrFactorWplusJets = 1.000*1.000;
  std::string sampleTTplusJets = "TTplusJets_madgraph";
  double weightFactorTTplusJets = dataIntLumi*157.5/1164640;   // inclusive TTbar xSection (MCFM @ NLO) / numEvents (MadGraph sample)
  double corrFactorTTplusJets = 1.009*1.000;

  std::vector<std::string> tauIds;
  tauIds.push_back(std::string("tauDiscrTaNCfrOnePercent"));  // "old" TaNC algorithm
  tauIds.push_back(std::string("tauDiscrTaNCfrHalfPercent"));
  tauIds.push_back(std::string("tauDiscrTaNCfrQuarterPercent"));
  //tauIds.push_back(std::string("tauDiscrTaNCloose")); // "new" TaNC implemented in HPS+TaNC combined algorithm
  //tauIds.push_back(std::string("tauDiscrTaNCmedium"));
  //tauIds.push_back(std::string("tauDiscrTaNCtight"));
  //tauIds.push_back(std::string("tauDiscrIsolationLoose"));    // "old" HPS algorithm
  //tauIds.push_back(std::string("tauDiscrIsolationMedium"));   
  //tauIds.push_back(std::string("tauDiscrIsolationTight"));
  tauIds.push_back(std::string("tauDiscrHPSloose"));  // "new" HPS implemented in HPS+TaNC combined algorithm
  tauIds.push_back(std::string("tauDiscrHPSmedium"));
  tauIds.push_back(std::string("tauDiscrHPStight"));

  std::vector<std::string> fitVariables;
  //fitVariables.push_back("diTauHt");
  //fitVariables.push_back("diTauSVfitMass1");
  //fitVariables.push_back("diTauSVfitMass2");
  //fitVariables.push_back("diTauVisMass");
  fitVariables.push_back("diTauVisMassFromJet");
  //fitVariables.push_back("muonPt");
  //fitVariables.push_back("muonEta");
  //fitVariables.push_back("tauPt");
  //fitVariables.push_back("tauEta");
  //fitVariables.push_back("tauNumChargedParticles");
  //fitVariables.push_back("tauNumParticles");
  //fitVariables.push_back("tauJetWidth"); CV: signal/background normalizations --> tau id. efficiencies obtained by using jetWidth variable
  //                                           are **very** different from values obtained by using all other variables
  //                                          --> there seems to be a problem in modeling jetWidth variable
  //                                          --> do not use jetWidth variable for now

  std::vector<std::string> sysShifts;
  sysShifts.push_back("CENTRAL_VALUE");
  if ( runSysUncertainties ) sysShifts.insert(sysShifts.end(), sysUncertainties.begin(), sysUncertainties.end());

  TFile* histogramInputFile = 0;
  if ( loadHistogramsFromFile ) histogramInputFile = new TFile(histogramFileName.data());

  for ( std::vector<std::string>::const_iterator sysShift = sysShifts.begin();
	sysShift != sysShifts.end(); ++sysShift ) {
    std::cout << "running fit for sysShift = " << (*sysShift) << "..." << std::endl;

//--- initialize alias --> branchName mapping
    typedef std::map<std::string, std::string> branchNameDictEntry;
    std::map<std::string, branchNameDictEntry> branchNamesData           
      = makeBranchNameDict(tauIds, "CENTRAL_VALUE", branchName_suffix, false);
    std::map<std::string, branchNameDictEntry> branchNamesMCwZrecoilCorr 
      = makeBranchNameDict(tauIds, *sysShift,       branchName_suffix, true);
    std::map<std::string, branchNameDictEntry> branchNamesMCwoZrecoilCorr 
      = makeBranchNameDict(tauIds, *sysShift,       branchName_suffix, false);
   
//--- define x-axis titles
    std::map<std::string, std::string> xAxisTitles;
    xAxisTitles["diTauCharge"]          = "Charge(#mu + #tau_{had})";
    xAxisTitles["diTauMt"]              = "M_{T}^{#muMET} [GeV]";
    xAxisTitles["diTauHt"]              = "P_{T}^{#mu} + P_{T}^{#tau} + MET [GeV]";
    xAxisTitles["diTauSVfitMass1"]      = "M^{#tau#tau} [GeV]";
    xAxisTitles["diTauSVfitMass2"]      = xAxisTitles["diTauSVfitMass1"];
    xAxisTitles["diTauVisMass"]         = "M_{vis}^{#mu#tau} [GeV]";
    xAxisTitles["diTauVisMassFromJet"]  = xAxisTitles["diTauVisMass"];

    std::map<std::string, std::map<std::string, TH1*> > distributionsData; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      if ( runClosureTest ) {
	std::cout << ">>> NOTE: RUNNING CLOSURE TEST <<<" << std::endl;
	distributionsData = loadHistograms(histogramInputFile, "sum", regions, tauIds, fitVariables, *sysShift);
      } else {
	distributionsData = loadHistograms(histogramInputFile, "Data", regions, tauIds, fitVariables, *sysShift);
      }
    } else {
      TChain* chainData_2011RunA = new TChain("Events");
      TChain* chainData_2011RunB = new TChain("Events");
      if ( !runQuickTest ) { addFileNames(chainData_2011RunA, inputFilePath, "data_Mu_Run2010A_Nov4ReReco", jobId);
  	                     addFileNames(chainData_2011RunB, inputFilePath, "data_Mu_Run2010B_Nov4ReReco", jobId); } 
      else                 { chainData_2011RunA->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010A_Nov4ReReco_2011Feb03b_0_8d77.root").data());
	                     chainData_2011RunB->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010B_Nov4ReReco_2011Feb03b_0_5508.root").data()); }
      printFileInfo(chainData_2011RunA, "chainData_2011RunA");
      printFileInfo(chainData_2011RunB, "chainData_2011RunB");
      TChain* chainData = new TChain("Events");
      chainData->Add(chainData_2011RunA);
      chainData->Add(chainData_2011RunB);
      printFileInfo(chainData, "chainData");
      distributionsData = makeDistributionsAllRegions("Data", 1.0, chainData, regions,
						      tauIds, fitVariables, branchNamesData, *sysShift, &saveRunLumiSectionEventNumbers);
      delete chainData;
      delete chainData_2011RunA;
      delete chainData_2011RunB;
    } 

    std::map<std::string, std::map<std::string, TH1*> > templatesZtautau; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      templatesZtautau = loadHistograms(histogramInputFile, "Ztautau", regions, tauIds, fitVariables, *sysShift);
    } else {
      TChain* chainZtautau = new TChain("Events");
      if ( !runQuickTest ) addFileNames(chainZtautau, inputFilePath, sampleZtautau, jobId);
      else                 chainZtautau->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_ZtautauPU156bx_2011Feb03b_0_53ca.root").data());
      printFileInfo(chainZtautau, "chainZtautau");
      templatesZtautau = makeDistributionsAllRegions("Ztautau", weightFactorZtautau*corrFactorZtautau, chainZtautau, regions, 
						     tauIds, fitVariables, branchNamesMCwZrecoilCorr, *sysShift);
      delete chainZtautau;
    }

    std::map<std::string, std::map<std::string, TH1*> > templatesZmumu; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      templatesZmumu = loadHistograms(histogramInputFile, "Zmumu", regions, tauIds, fitVariables, *sysShift);
    } else {
      TChain* chainZmumu = new TChain("Events");
      if ( !runQuickTest ) addFileNames(chainZmumu, inputFilePath, sampleZmumu, jobId);
      else                 chainZmumu->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Zmumu_pythia_2011Feb03b_0_2f1a.root").data());
      printFileInfo(chainZmumu, "chainZmumu");
      templatesZmumu = makeDistributionsAllRegions("Zmumu", weightFactorZmumu*corrFactorZmumu, chainZmumu, regions, 
						   tauIds, fitVariables, branchNamesMCwZrecoilCorr, *sysShift);
      delete chainZmumu;
    }

    std::map<std::string, std::map<std::string, TH1*> > templatesQCD; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      templatesQCD = loadHistograms(histogramInputFile, "QCD", regions, tauIds, fitVariables, *sysShift);
    } else {
      TChain* chainQCD = new TChain("Events");
      if ( !runQuickTest ) addFileNames(chainQCD, inputFilePath, sampleQCD, jobId);
      else                 chainQCD->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Feb03b_0_35f8.root").data());
      printFileInfo(chainQCD, "chainQCD");
      templatesQCD = makeDistributionsAllRegions("QCD", weightFactorQCD*corrFactorQCD, chainQCD, regions, 
						 tauIds, fitVariables, branchNamesMCwoZrecoilCorr, *sysShift);
      delete chainQCD;
    }

    std::map<std::string, std::map<std::string, TH1*> > templatesWplusJets; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      templatesWplusJets = loadHistograms(histogramInputFile, "WplusJets", regions, tauIds, fitVariables, *sysShift);
    } else {
      TChain* chainWplusJets = new TChain("Events");
      if ( !runQuickTest ) addFileNames(chainWplusJets, inputFilePath, sampleWplusJets, jobId);
      else                 chainWplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Feb03b_0_e872.root").data());
      printFileInfo(chainWplusJets, "chainWplusJets");
      templatesWplusJets = makeDistributionsAllRegions("WplusJets", weightFactorWplusJets*corrFactorWplusJets, chainWplusJets, regions, 
						       tauIds, fitVariables, branchNamesMCwoZrecoilCorr, *sysShift);
      delete chainWplusJets;
    }

    std::map<std::string, std::map<std::string, TH1*> > templatesTTplusJets; // key = (region, observable)
    if ( loadHistogramsFromFile ) { 
      templatesTTplusJets = loadHistograms(histogramInputFile, "TTplusJets", regions, tauIds, fitVariables, *sysShift);
    } else {
      TChain* chainTTplusJets = new TChain("Events");
      if ( !runQuickTest ) addFileNames(chainTTplusJets, inputFilePath, sampleTTplusJets, jobId);
      else                 chainTTplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Feb03b_0_0f74.root").data());
      printFileInfo(chainTTplusJets, "chainTTplusJets");
      templatesTTplusJets = makeDistributionsAllRegions("TTplusJets", weightFactorTTplusJets*corrFactorTTplusJets, chainTTplusJets, regions, 
							tauIds, fitVariables, branchNamesMCwoZrecoilCorr, *sysShift);
      delete chainTTplusJets;
    }

    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > > templatesAll; // key = (process, region, observable)
    templatesAll["Ztautau"]    = templatesZtautau;
    templatesAll["Zmumu"]      = templatesZmumu;
    templatesAll["QCD"]        = templatesQCD;
    templatesAll["WplusJets"]  = templatesWplusJets;
    templatesAll["TTplusJets"] = templatesTTplusJets;
    
    std::vector<std::string> processes;
    processes.push_back(std::string("Ztautau"));
    processes.push_back(std::string("Zmumu"));
    processes.push_back(std::string("QCD"));
    processes.push_back(std::string("WplusJets"));
    processes.push_back(std::string("TTplusJets"));
    
//--- closure test: fit sum(MC) instead of Data
    if ( runClosureTest ) {
      std::cout << ">>> NOTE: RUNNING CLOSURE TEST <<<" << std::endl;
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	      region != distributionsData.end(); ++region ) {
	  for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
		key != region->second.end(); ++key ) {
	    TH1* histogram = templatesAll[*process][region->first][key->first];
	    
	    TH1* histogramSum = templatesAll["sum"][region->first][key->first];
	    if ( !histogramSum ) {
	      std::string histogramName = histogram->GetName();
	      std::string histogramSumName = std::string("sum").append(std::string(histogramName, histogramName.find("_")));	    
	      templatesAll["sum"][region->first][key->first] = (TH1*)histogram->Clone(histogramSumName.data());
	    } else {	    
	      histogramSum->Add(histogram);
	    }
	  }
	}
      }
      
      std::cout << std::endl;
      
      distributionsData = templatesAll["sum"];
    }

//--- obtain template for QCD background from data,
//    from SS && Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV sideband
/*
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	std::string key_all = getKey(*fitVariable, *tauId);
	std::string key_passed = getKey(*fitVariable, *tauId, "passed");
	std::string key_failed = getKey(*fitVariable, *tauId, "failed");

	std::cout << "templatesQCD['C1'][" << key_all << "] = " << templatesQCD["C1"][key_all] << std::endl;
	std::cout << "templatesQCD['C1p'][" << key_passed << "] = " << templatesQCD["C1p"][key_passed] << std::endl;
	std::cout << "templatesQCD['C1f'][" << key_failed << "] = " << templatesQCD["C1f"][key_failed] << std::endl;
	std::cout << "distributionsData['B1'][" << key_all << "] = " << distributionsData["B1"][key_all] << std::endl;

	std::string histogramNameQCD_C1 = templatesQCD["C1"][key_all]->GetName();
	double normQCD_C1 = getIntegral(templatesQCD["C1"][key_all], true, true);
	templatesQCD["C1"][key_all] = normalize(distributionsData["B1"][key_all], normQCD_C1);
	templatesQCD["C1"][key_all]->SetName(histogramNameQCD_C1.data());

	std::string histogramNameQCD_C1p = templatesQCD["C1p"][key_passed]->GetName();
	double normQCD_C1p = getIntegral(templatesQCD["C1p"][key_passed], true, true);	
	templatesQCD["C1p"][key_passed] = normalize(distributionsData["B1"][key_all], normQCD_C1p);
	templatesQCD["C1p"][key_passed]->SetName(histogramNameQCD_C1p.data());

	std::string histogramNameQCD_C1f = templatesQCD["C1f"][key_failed]->GetName();
	double normQCD_C1f = getIntegral(templatesQCD["C1f"][key_failed], true, true);
	templatesQCD["C1f"][key_failed] = normalize(distributionsData["B1"][key_all], normQCD_C1f);
	templatesQCD["C1f"][key_failed]->SetName(histogramNameQCD_C1f.data());
      }
    }

    templatesAll["QCD"] = templatesQCD;
 */
//--- save histograms
    if ( saveHistogramsToFile ) {
      TFile* histogramOutputFile = new TFile(histogramFileName.data(), "RECREATE");
      saveHistograms(histogramOutputFile, templatesZtautau);
      saveHistograms(histogramOutputFile, templatesZmumu);
      saveHistograms(histogramOutputFile, templatesQCD);
      saveHistograms(histogramOutputFile, templatesWplusJets);
      saveHistograms(histogramOutputFile, templatesTTplusJets);
      saveHistograms(histogramOutputFile, distributionsData);
      delete histogramOutputFile;
    }

//--- make control plots for sum(MC) scaled by cross-sections versus Data 
//    for Mt, fitVariable distributions in different regions
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	drawHistograms(templatesZtautau[region->first][key->first], -1.,
		       templatesZmumu[region->first][key->first], -1.,
		       templatesQCD[region->first][key->first], -1.,
		       templatesWplusJets[region->first][key->first], -1.,
		       templatesTTplusJets[region->first][key->first], -1.,
		       distributionsData[region->first][key->first],
		       std::string("Region ").append(region->first).append(": ").append(key->first).append(" (scaled by cross-section)"),
		       xAxisTitles[key->first],
		       std::string("controlPlotsTauIdEff_").append(region->first).append("_").append(key->first).append(".pdf"), *sysShift);
      }
    }

//--- print MC expectations for probabilities 
//   o pDiTauCharge_OS_SS
//   o pDiTauKine_Sig_Bgr
//   o pMuonIso_loose_tight
//   o pTauId_passed_failed
//    separating different regions

    std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll;    // key = (process/"sum", region, observable)
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions; // key = (process, region, observable)

    for ( std::vector<std::string>::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	    region != distributionsData.end(); ++region ) {
	for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	      key != region->second.end(); ++key ) {
	  numEventsAll[*process][region->first][key->first] = getIntegral(templatesAll[*process][region->first][key->first], true, true);
	  numEventsAll["sum"][region->first][key->first] += numEventsAll[*process][region->first][key->first];

	  fittedFractions[*process][region->first][key->first] = 
	    getIntegral(templatesAll[*process][region->first][key->first], false, false)/numEventsAll[*process][region->first][key->first];
	  //std::cout << "fittedFractions[" << (*process) << "][" << region->first << "][" << key->first << "] = "
	  //	      << fittedFractions[*process][region->first][key->first] << std::endl;
	}
      }
    }
    
    std::cout << std::endl;

    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	double numEventsA  = numEventsAll[*process]["A"][getKey("diTauMt", *tauId)];
	double numEventsB  = numEventsAll[*process]["B"][getKey("diTauMt", *tauId)];
	double numEventsC  = numEventsAll[*process]["C"][getKey("diTauMt", *tauId)];
	double numEventsC1 = numEventsAll[*process]["C1"][getKey("diTauMt", *tauId)];
	double numEventsC2 = numEventsAll[*process]["C2"][getKey("diTauMt", *tauId)];
	double numEventsD  = numEventsAll[*process]["D"][getKey("diTauMt", *tauId)];
	
	std::cout << "process = " << (*process) << std::endl;
	for ( size_t i = 0; i < (process->length() + 10); ++i ) {
	  std::cout << "-";
	}
	std::cout << std::endl;
	
	std::cout << "pDiTauCharge_OS_SS:" << std::endl;
	std::cout << " A/B = " << numEventsA/numEventsB << std::endl;
	std::cout << " C/D = " << numEventsC/numEventsD << std::endl;
	std::cout << "pMuonIso_loose_tight:" << std::endl;
	std::cout << " C/(A+C) = " << numEventsC/(numEventsA + numEventsC) << std::endl;
	std::cout << " D/(B+D) = " << numEventsD/(numEventsB + numEventsD) << std::endl;
	std::cout << "pDiTauKine_Sig_Bgr:" << std::endl;
	std::cout << " C1/C = " << numEventsC1/numEventsC << std::endl;
	
	std::cout << "pTauId_passed_failed:" << std::endl;

	double numEventsC1p = numEventsAll[*process]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
	double numEventsC2p = numEventsAll[*process]["C2p"][getKey("diTauMt", *tauId, "passed")];
	  
	std::cout << " " << (*tauId) << ":" << std::endl;
	std::cout << "  C1p/C1 = " << numEventsC1p/numEventsC1 << std::endl;
	std::cout << "  C2p/C2 = " << numEventsC2p/numEventsC2 << std::endl;
	
	std::cout << std::endl;
      }
    }

    std::map<std::string, std::map<std::string, double> > effValues; // key = (tauId, fitVariable)
    std::map<std::string, std::map<std::string, double> > effErrors; // key = (tauId, fitVariable)

    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	double effValue = 0.;
	double effError = 1.;
	fitUsingRooFit(distributionsData, templatesAll, numEventsAll, fittedFractions,
		       processes,
		       *tauId, *fitVariable, fitTauIdEffC2,
		       effValue, effError,		       
		       *sysShift,
		       &xAxisTitles);
	effValues[*tauId][*fitVariable] = effValue;
	effErrors[*tauId][*fitVariable] = effError;
      }
    }
    
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      std::cout << "Efficiency of Tau id. = " << (*tauId) << ":" << std::endl;
      
      for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	std::cout << " fitVariable = " << (*fitVariable) << ":" 
		  << " result = " << effValues[*tauId][*fitVariable]*100. << " +/- " << effErrors[*tauId][*fitVariable]*100. << "%" << std::endl;
	
	double numEventsC    = numEventsAll["Ztautau"]["C"][getKey("diTauMt", *tauId)];
        double numEventsC1   = numEventsAll["Ztautau"]["C1"][getKey("diTauMt", *tauId)];
	double numEventsC1p  = numEventsAll["Ztautau"]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
	double numEventsC2p  = numEventsAll["Ztautau"]["C2p"][getKey("diTauMt", *tauId, "passed")];

	double tauIdEffMCexp = ( fitTauIdEffC2 ) ? (numEventsC1p + numEventsC2p)/numEventsC : numEventsC1p/numEventsC1;
        std::cout << "(Monte Carlo prediction = " << tauIdEffMCexp*100. << "%)" << std::endl;
      }
    }
  }

  if ( loadHistogramsFromFile ) delete histogramInputFile;

//--print time that it took macro to run
  std::cout << "finished executing fitTauIdEff_wConstraints macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  std::cout << " #sysShifts    = " << sysShifts.size() << std::endl;
  clock.Show("fitTauIdEff_wConstraints");
}
