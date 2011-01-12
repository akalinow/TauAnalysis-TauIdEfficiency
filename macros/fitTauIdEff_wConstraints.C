
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

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TFractionFitter.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TBenchmark.h>
#include <TSystem.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

std::string getBranchName(const std::string& branchName_object, const std::string& branchName_observable,
			  const std::string& branchName_suffix)
{
//--------------------------------------------------------------------------------
// Compose full branchName from name of pat::Muon/pat::Tau/diTau collection
// used when producing (ED)Ntuple, observable and branchName "suffix" ("local"/"lxbatch")
//--------------------------------------------------------------------------------

  return std::string(branchName_object).append("#").append(branchName_observable).append("_").append(branchName_suffix).append(".obj");
}

std::map<std::string, std::string> makeBranchNameDict(const std::string& sysShift, const std::string& branchName_suffix)
{
//--------------------------------------------------------------------------------
// Make alias --> branchName mapping 
// to be used for treeSelection and histogram filling
//--------------------------------------------------------------------------------

  std::map<std::string, std::string> branchNames;
  
  const std::string ntupleName = "double_ntupleProducer_tauIdEffNtuple";

  std::string branchNameMuon = std::string(ntupleName).append("#selectedPatMuonsTrkIPcumulative");
  std::string branchNameTau = std::string(ntupleName).append("#selectedPatTausForMuTauEcalCrackVetoCumulative");
  std::string branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoCumulative");

  if (        sysShift == "CENTRAL_VALUE"              ) {
    // nothing to be done yet.
  } else if ( sysShift == "SysTauJetEnUp"              ) {
    branchNameTau   = std::string(ntupleName).append("#selectedPatTausForMuTauEcalCrackVetoSysTauJetEnUpCumulative");
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysTauJetEnUpCumulative");
  } else if ( sysShift == "SysTauJetEnDown"            ) { 
    branchNameTau   = std::string(ntupleName).append("#selectedPatTausForMuTauEcalCrackVetoSysTauJetEnDownCumulative");
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysTauJetEnDownCumulative");
  } else if ( sysShift == "SysJetEnUp"                 ) {
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysJetEnUpCumulative");
  } else if ( sysShift == "SysJetEnDown"               ) {
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysJetEnDownCumulative");
  } else if ( sysShift == "SysZllRecoilCorrectionUp"   ) {
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysZllRecoilCorrectionUpCumulative");
  } else if ( sysShift == "SysZllRecoilCorrectionDown" ) {
    branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsAntiOverlapVetoSysZllRecoilCorrectionDownCumulative");
  } else assert(0);

  branchNames["muonPt"] = getBranchName(branchNameMuon, "pt", branchName_suffix);
  branchNames["muonEta"] = getBranchName(branchNameMuon, "eta", branchName_suffix);
  branchNames["muonLooseIsoPtSum04"] = getBranchName(branchNameMuon, "ptSumLooseIsolation04", branchName_suffix);
  branchNames["muonLooseIsoPtSum06"] = getBranchName(branchNameMuon, "ptSumLooseIsolation06", branchName_suffix);

  branchNames["tauPt"] = getBranchName(branchNameTau, "pt", branchName_suffix);
  branchNames["tauEta"] = getBranchName(branchNameTau, "eta", branchName_suffix);
  branchNames["tauLooseIsoPtSum04"] = getBranchName(branchNameTau, "ptSumLooseIsolation04", branchName_suffix);
  branchNames["tauLooseIsoPtSum06"] = getBranchName(branchNameTau, "ptSumLooseIsolation06", branchName_suffix);
  branchNames["tauNumChargedParticles"] = getBranchName(branchNameTau, "numChargedParticles", branchName_suffix);
  branchNames["tauNumParticles"] = getBranchName(branchNameTau, "numParticles", branchName_suffix);
  branchNames["tauJetPt"] = getBranchName(branchNameTau, "jetPt", branchName_suffix);
  branchNames["tauJetEta"] = getBranchName(branchNameTau, "jetEta", branchName_suffix);
  branchNames["tauJetWidth"] = getBranchName(branchNameTau, "jetWidth", branchName_suffix);
  branchNames["tauDiscrHPSloose"] = getBranchName(branchNameTau, "byHPSloose", branchName_suffix);
  branchNames["tauDiscrHPSmedium"] = getBranchName(branchNameTau, "byHPSmedium", branchName_suffix);
  branchNames["tauDiscrHPStight"] = getBranchName(branchNameTau, "byHPStight", branchName_suffix);
  branchNames["tauDiscrTaNCloose"] = getBranchName(branchNameTau, "byTaNCloose", branchName_suffix);
  branchNames["tauDiscrTaNCmedium"] = getBranchName(branchNameTau, "byTaNCmedium", branchName_suffix);
  branchNames["tauDiscrTaNCtight"] = getBranchName(branchNameTau, "byTaNCtight", branchName_suffix);

  branchNames["diTauCharge"] = getBranchName(branchNameDiTau, "charge", branchName_suffix);
  branchNames["diTauMt"] = getBranchName(branchNameDiTau, "Mt", branchName_suffix);
  branchNames["diTauPzeta"] = getBranchName(branchNameDiTau, "pZeta", branchName_suffix);
  branchNames["diTauPzetaVis"] = getBranchName(branchNameDiTau, "pZetaVis", branchName_suffix);
  branchNames["diTauHt"] = getBranchName(branchNameDiTau, "Ht", branchName_suffix);
  branchNames["diTauSVfitMass1"] = getBranchName(branchNameDiTau, "SVfitMass1", branchName_suffix);
  branchNames["diTauSVfitMass2"] = getBranchName(branchNameDiTau, "SVfitMass2", branchName_suffix);
  branchNames["diTauVisMass"] = getBranchName(branchNameDiTau, "visMass", branchName_suffix);

  return branchNames;
}

std::string getKey(const std::string& observable, const std::string& tauId = "any", const std::string tauIdValue = "all")
{
  std::string key = std::string(observable).append("_").append(tauId).append("_").append(tauIdValue);
  return key;
}

std::map<std::string, TH1*> makeHistograms(const std::string& process, const std::string& region, double weight,
					   TTree* tree, const std::string& treeSelection, 
					   const std::vector<std::string>* tauIds, const std::vector<std::string>* tauIdValues,
					   const std::vector<std::string>& observables,
					   std::map<std::string, std::string>& branchNames, const std::string& sysShift = "CENTRAL_VALUE")
{
//--------------------------------------------------------------------------------
// Fill histograms with (ED)NTuple entries passing treeSelection
//--------------------------------------------------------------------------------

  std::cout << "<makeHistograms>:" << std::endl;
  std::cout << " treeSelection = " << treeSelection << std::endl;

  std::map<std::string, TH1*> retVal;

//--- prepare "dummy" tau id. selection
//    in case no tau id. selection has been passed as function argument
  std::vector<std::string> tauIds_local;
  std::vector<std::string> tauIdValues_local;
  if ( tauIds && tauIdValues ) {
    tauIds_local = (*tauIds);
    tauIdValues_local = (*tauIdValues);
  } else {
    tauIds_local.push_back("any");
    tauIdValues_local.push_back("all");
  }

//--- apply kinematic cuts common to all regions
//   (cuts that might as well have been applied during (ED)Ntuple production)
  std::string extTreeSelection = treeSelection;
  if ( extTreeSelection != "" ) extTreeSelection.append(" && ");
  extTreeSelection.append(branchNames["muonPt"]).append(" > 20.");
  extTreeSelection.append(" && ").append(branchNames["tauPt"]).append(" > 20.");
  extTreeSelection.append(" && abs(").append(branchNames["tauEta"]).append(") < 2.3");
  extTreeSelection.append(" && ").append(branchNames["tauLooseIsoPtSum06"]).append(" < 2.5");

  TTree* selEventsTree = 0;
  if ( extTreeSelection != "" ) {
    selEventsTree = tree->CopyTree(extTreeSelection.data());
  } else {
    selEventsTree = tree;
  }

  //selEventsTree->Print();

  std::cout << "process = " << process << " has " << selEventsTree->GetEntries() << " entries." << std::endl;
  std::cout << " (weighted sum = " << selEventsTree->GetEntries()*weight << ")" << std::endl;

  for ( std::vector<std::string>::const_iterator tauId = tauIds_local.begin();
	tauId != tauIds_local.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues_local.begin();
	  tauIdValue != tauIdValues_local.end(); ++tauIdValue ) {

      std::string treeSelectionForTauIdValue;
      if      ( (*tauIdValue) == "passed" ) treeSelectionForTauIdValue = std::string(branchNames[*tauId]).append(" > 0.5");
      else if ( (*tauIdValue) == "failed" ) treeSelectionForTauIdValue = std::string(branchNames[*tauId]).append(" < 0.5");
      else if ( (*tauIdValue) == "all"    ) treeSelectionForTauIdValue = "";
      else assert(0);
      
      TTree* selEventsTreeForTauIdValue = 0;
      if ( treeSelectionForTauIdValue != "" ) {
	selEventsTreeForTauIdValue = selEventsTree->CopyTree(treeSelectionForTauIdValue.data());
      } else {
	selEventsTreeForTauIdValue = selEventsTree;
      }

      std::cout << "--> " << selEventsTreeForTauIdValue->GetEntries() << " entries selected" 
		<< " for tauId = " << (*tauId) << ", value = " << (*tauIdValue) << std::endl;

      for ( std::vector<std::string>::const_iterator observable = observables.begin();
	    observable != observables.end(); ++observable ) {
	
	std::string histogramName = std::string(process).append("_").append(region).append("_").append(*observable);
	histogramName.append("_").append(*tauId).append("_").append(*tauIdValue);
	if ( sysShift != "CENTRAL_VALUE" ) histogramName.append("_").append(sysShift);

	int numBins;
	double min, max;

	if ( (*observable) == "diTauMt" ) {
	  numBins = 16;
	  min = 0.;
	  max = 80.;
	} else if ( (*observable) == "diTauVisMass" ) {
	  numBins = 20;
	  min = 20.;
	  max = 120.;
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
	  std::cout << "Error in <makeHistograms>: undefined observable = " << (*observable) << " --> skipping !!";
	  return retVal;
	}
	
	TH1* histogram = new TH1F(histogramName.data(), histogramName.data(), numBins, min, max);
	
        std::string drawCommand = std::string(branchNames[*observable]).append(">>+").append(histogramName); 
	selEventsTreeForTauIdValue->Draw(drawCommand.data());
	
	std::cout << "histogram = " << histogramName << ": integral = " << histogram->Integral() << std::endl;

	if ( !histogram->GetSumw2N() ) histogram->Sumw2();
	histogram->Scale(weight);

	std::string key = getKey(*observable, *tauId, *tauIdValue);	
	if ( histogram != 0 ) retVal[key] = histogram;
	
        if ( selEventsTreeForTauIdValue != selEventsTree ) delete selEventsTreeForTauIdValue;
      }
    }
  }

  std::cout << std::endl;

  if ( selEventsTree != tree ) delete selEventsTree;

  return retVal;
}

TH1* normalize(const TH1* histogram, double norm = 1.)
{
//--------------------------------------------------------------------------------
// Normalize histogram passed as function argument to unit area
// (for use as shape template)
//--------------------------------------------------------------------------------

  TH1* retVal = (TH1*)histogram->Clone();

  if ( !retVal->GetSumw2N() ) retVal->Sumw2();

  if ( retVal->Integral() != 0. ) retVal->Scale(norm/retVal->Integral());
    
  return retVal;
}

void drawHistograms(TH1* histogramQCD, double normQCD,
		    TH1* histogramWplusJets, double normWplusJets,
		    TH1* histogramZmumu, double normZmumu,
		    TH1* histogramZtautau, double normZtautau,
		    TH1* histogramData,
		    const std::string& title,
		    const std::string& outputFileName, const std::string& sysShift = "CENTRAL_VALUE")
{
//--------------------------------------------------------------------------------
// Make control plots of sum(MC) versus Data.
// If normalization factors are passed as function argument,
// normalize all MC distributions accordingly;
// else assume MC distributions passed as function arguments
// are already properly normalized (by cross-section)
//--------------------------------------------------------------------------------

  //std::cout << "<drawHistograms (1)>:" << std::endl;
  //std::cout << " normQCD = " << normQCD << std::endl;
  //std::cout << " normWplusJets = " << normWplusJets << std::endl;
  //std::cout << " normZmumu = " << normZmumu << std::endl;
  //std::cout << " normZtautau = " << normZtautau << std::endl;
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  THStack smSum("smSum", "smSum");

  TH1* templateQCD = ( normQCD > 0. ) ? normalize(histogramQCD, normQCD) : histogramQCD;
  templateQCD->SetFillColor(797);
  smSum.Add(templateQCD);

  TH1* templateWplusJets = ( normWplusJets > 0. ) ? normalize(histogramWplusJets, normWplusJets) : histogramWplusJets;
  templateWplusJets->SetFillColor(856);
  smSum.Add(templateWplusJets);

  TH1* templateZmumu = ( normZmumu > 0. ) ? normalize(histogramZmumu, normZmumu) : histogramZmumu;
  templateZmumu->SetFillColor(596);
  smSum.Add(templateZmumu);

  TH1* templateZtautau = ( normZtautau > 0. ) ? normalize(histogramZtautau, normZtautau) : histogramZtautau;
  templateZtautau->SetFillColor(628);
  smSum.Add(templateZtautau);

  histogramData->SetLineColor(1);
  histogramData->SetMarkerColor(1);
  histogramData->SetMarkerStyle(20);

  if ( title != "" )
    smSum.SetTitle(title.data());
  else
    smSum.SetTitle(templateZtautau->GetTitle());
  smSum.SetMaximum(1.4*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));
	
  smSum.Draw("hist");
  histogramData->SetStats(false);
  histogramData->Draw("ep1same");

  TLegend legend(0.64, 0.69, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  legend.AddEntry(templateZtautau, "Z #rightarrow #tau^{+} #tau^{-}", "f");
  legend.AddEntry(templateZmumu, "Z #rightarrow #mu^{+} #mu^{-}", "f");
  legend.AddEntry(templateWplusJets, "W + jets", "f");
  legend.AddEntry(templateQCD, "QCD", "f");
  legend.AddEntry(histogramData, "Data", "p");
  legend.Draw();

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  if ( sysShift != "CENTRAL_VALUE" ) outputFilePath.append(sysShift).append("/");
  gSystem->mkdir(outputFilePath.data(), true);
  canvas->Print(outputFilePath.append(outputFileName).data());

  if ( templateQCD       != histogramQCD       ) delete templateQCD;
  if ( templateWplusJets != histogramWplusJets ) delete templateWplusJets;
  if ( templateZmumu     != histogramZmumu     ) delete templateZmumu;
  if ( templateZtautau   != histogramZtautau   ) delete templateZtautau;

  delete canvas;
}

void drawHistograms(TH1* histogram_passed, TH1* histogram_failed, 
		    const std::string& tauId,
		    const std::string& outputFileName, const std::string& sysShift = "CENTRAL_VALUE")
{
//--------------------------------------------------------------------------------
// Make control plots of MC distribution
// in tau id. passed versus tau id. failed regions
//--------------------------------------------------------------------------------

  //std::cout << "<drawHistograms (2)>:" << std::endl;
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  TH1* template_passed = normalize(histogram_passed, 1.);
  template_passed->SetLineColor(3);
  template_passed->SetMarkerColor(3);
  template_passed->SetMarkerStyle(20);

  TH1* template_failed = normalize(histogram_failed, 1.);
  template_failed->SetLineColor(2);
  template_failed->SetMarkerColor(2);
  template_failed->SetMarkerStyle(24);

  template_passed->SetMaximum(1.4*TMath::Max(template_passed->GetMaximum(), template_failed->GetMaximum()));
  template_passed->SetStats(false);
  template_passed->Draw("ep1");

  template_failed->Draw("ep1same");

  TLegend legend(0.64, 0.69, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  
  legend.AddEntry(template_passed, std::string(tauId).append(" passed").data(), "p");
  legend.AddEntry(template_failed, std::string(tauId).append(" failed").data(), "p");
  legend.Draw();

  canvas->Update();
  std::string outputFilePath = std::string("./plots/");
  if ( sysShift != "CENTRAL_VALUE" ) outputFilePath.append(sysShift).append("/");
  gSystem->mkdir(outputFilePath.data(), true);
  canvas->Print(outputFilePath.append(outputFileName).data());

  delete template_passed;
  delete template_failed;

  delete canvas;
}

void addToFormula(std::string& formula, const std::string& expression, TObjArray& arguments, RooRealVar* p)
{
//-------------------------------------------------------------------------------
// Build formula and argument list for RooFormulaVar
//
// NOTE: auxiliary function for makeRooFormulaVar
//
//-------------------------------------------------------------------------------

  if ( expression != "" ) {
    if ( formula != "" ) formula.append(" * ");

    if      ( expression == "regular"  ) formula.append(p->GetName());
    else if ( expression == "inverted" ) formula.append("(1 - ").append(p->GetName()).append(")");
    else assert(0);
    
    arguments.Add(p);
  }
}

RooFormulaVar* makeRooFormulaVar(const std::string& process, const std::string& region,
				 RooRealVar* norm, 
				 RooRealVar* pDiTauCharge_OS_SS, RooRealVar* pDiTauKine_Sig_Bgr, 
				 RooRealVar* pMuonIso_loose_tight, 
				 RooRealVar* pTauId_passed_failed)
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

  std::cout << "building RooFormulaVar expression for process = " << process << ", region = " << region << "." << std::endl;

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
    std::cout << "Error in <makeRooFormulaVar>: undefined region = " << region << " --> skipping !!";
    return retVal;
  }

  std::string formula = "";
  TObjArray arguments; 
  addToFormula(formula, exprDiTauCharge, arguments, pDiTauCharge_OS_SS);
  addToFormula(formula, exprDiTauKine,   arguments, pDiTauKine_Sig_Bgr);
  addToFormula(formula, exprMuonIso,     arguments, pMuonIso_loose_tight);
  addToFormula(formula, exprTauId,       arguments, pTauId_passed_failed);
  addToFormula(formula, "regular",       arguments, norm);

  std::cout << " formula = " << formula << std::endl;
  //std::cout << " arguments:" << std::endl;
  //for ( int i = 0; i < arguments.GetEntries(); ++i ) {
  //  std::cout << " " << dynamic_cast<TNamed*>(arguments.At(i))->GetName() << std::endl;
  //}

  std::string name = std::string("norm").append(process).append("_").append(region);
  retVal = new RooFormulaVar(name.data(), name.data(), formula.data(), RooArgSet(arguments));

  return retVal;
}

std::map<std::string, TH1*> makeDistributionsInRegion(const std::string& process, const std::string& region, double weight,
						      TTree* tree, 
						      const std::vector<std::string>& tauIds,
						      const std::vector<std::string>& fitVariables,
						      std::map<std::string, std::string>& branchNames, const std::string& sysShift = "CENTRAL_VALUE")
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
  
  std::cout << "selecting " << process << " in region " << region << "..." << std::endl;

  std::string treeSelection;
  
  std::string exprMuonIso_loose  = std::string(branchNames["muonLooseIsoPtSum06"]).append(" > (0.10*").append(branchNames["muonPt"]).append(")");
  exprMuonIso_loose.append(" && ").append(branchNames["muonLooseIsoPtSum06"]).append(" < (0.30*").append(branchNames["muonPt"]).append(")");
  std::string exprMuonIso_tight  = std::string(branchNames["muonLooseIsoPtSum06"]).append(" < (0.10*").append(branchNames["muonPt"]).append(")");
  std::string exprDiTauCharge_OS = std::string("abs(").append(branchNames["diTauCharge"]).append(") < 0.5");
  std::string exprDiTauCharge_SS = std::string("abs(").append(branchNames["diTauCharge"]).append(") > 1.5");
  std::string exprDiTauKine_Sig  = std::string("(").append(branchNames["diTauMt"]).append(" < 40");
  exprDiTauKine_Sig.append(" && ").append("(").append(branchNames["diTauPzeta"]).append(" - 1.5*").append(branchNames["diTauPzetaVis"]).append(") > -20)");
  std::string exprDiTauKine_Bgr  = std::string("(").append(branchNames["diTauMt"]).append(" > 40");
  exprDiTauKine_Bgr.append(" || ").append("(").append(branchNames["diTauPzeta"]).append(" - 1.5*").append(branchNames["diTauPzetaVis"]).append(") < -20)");

  std::vector<std::string> observables;

  const std::vector<std::string>* tauIds_local = 0;
  std::vector<std::string> tauIdValues_local;

  if        ( region == "A" ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprDiTauCharge_OS);
    observables.push_back(std::string("diTauMt"));
  } else if ( region == "B" ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprDiTauCharge_SS);
    observables.push_back(std::string("diTauMt"));
  } else if ( region.find("C") == 0 ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_OS);
    if      ( region.find("C1") == 0 ) treeSelection.append(" && ").append(exprDiTauKine_Sig);
    else if ( region.find("C2") == 0 ) treeSelection.append(" && ").append(exprDiTauKine_Bgr);

    if ( region == "C" || region == "C1" || region == "C2" ) {
      observables.push_back(std::string("diTauMt"));
    } else if ( region == "C1p" ) {
      observables = fitVariables;
      tauIds_local = &tauIds;
      tauIdValues_local.push_back(std::string("passed"));
    } else if ( region == "C1f" ) {
      observables = fitVariables;
      tauIds_local = &tauIds;
      tauIdValues_local.push_back(std::string("failed"));
    } else if ( region == "C2p" ) {
      observables.push_back(std::string("diTauMt"));
      tauIds_local = &tauIds;
      tauIdValues_local.push_back(std::string("passed"));
    } else if ( region == "C2f" ) {
      observables.push_back(std::string("diTauMt"));
      tauIds_local = &tauIds;
      tauIdValues_local.push_back(std::string("failed"));
    }
  } else if ( region == "D" ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_SS);
    observables.push_back(std::string("diTauMt"));    
  } else {
    std::cout << "Error in <makeDistribution>: undefined region = " << region << " --> skipping !!";
    return std::map<std::string, TH1*>();
  }

  return makeHistograms(process, region, weight,
			tree, treeSelection,
			tauIds_local, &tauIdValues_local,
			observables,
			branchNames, sysShift);
}

std::map<std::string, std::map<std::string, TH1*> > makeDistributionsAllRegions(const std::string& process, double weight,
										TTree* tree, 
										const std::vector<std::string>& tauIds,
										const std::vector<std::string>& fitVariables,
										std::map<std::string, std::string>& branchNames, const std::string& sysShift = "CENTRAL_VALUE")
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
  
  std::map<std::string, std::map<std::string, TH1*> > retVal;
  
  std::vector<std::string> regions;
  regions.push_back(std::string("A"));
  regions.push_back(std::string("B"));
  regions.push_back(std::string("C"));
  regions.push_back(std::string("C1"));
  regions.push_back(std::string("C1p"));
  regions.push_back(std::string("C1f"));
  regions.push_back(std::string("C2"));
  regions.push_back(std::string("C2p"));
  regions.push_back(std::string("C2f"));
  regions.push_back(std::string("D"));

  for ( std::vector<std::string>::const_iterator region = regions.begin();
	region != regions.end(); ++region ) {
    retVal[*region] = makeDistributionsInRegion(process, *region, weight,
						tree, 
						tauIds,
						fitVariables,
						branchNames, sysShift);
  }

  return retVal;
}

RooHistPdf* makeRooHistPdf(TH1* templateHistogram, RooRealVar* fitVar)
{
  std::string templateDataHistName = std::string(templateHistogram->GetName()).append("_dataHist");
  RooDataHist* templateDataHist = new RooDataHist(templateDataHistName.data(), templateDataHistName.data(), *fitVar, templateHistogram);
  std::string templatePdfName = std::string(templateHistogram->GetName()).append("_histPdf");
  RooHistPdf* templatePdf = new RooHistPdf(templatePdfName.data(), templatePdfName.data(), *fitVar, *templateDataHist);
  return templatePdf;
}

RooGaussian* makeFitConstraint(RooRealVar* p, double value, double error)
{
  std::string pValueName = std::string(p->GetName()).append("_constValue");
  RooConstVar* pValue = new RooConstVar(pValueName.data(), pValueName.data(), value);
  std::string pErrorName = std::string(p->GetName()).append("_constError");
  RooConstVar* pError = new RooConstVar(pErrorName.data(), pErrorName.data(), error);

  std::string constraintName = std::string(p->GetName()).append("_constraint");
  RooGaussian* constraint = new RooGaussian(constraintName.data(), constraintName.data(), *p, *pValue, *pError);
  return constraint;
}

void fitUsingRooFit(std::map<std::string, std::map<std::string, TH1*> >& distributionsData,                      // key = (region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,   // key = (process, region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll, // key = (process, region, observable)
		    const std::vector<std::string>& processes,
		    const std::string& tauId, const std::string& fitVariable,
		    double& effValue, double& effError,
		    const std::string& sysShift = "CENTRAL_VALUE") 
{
  std::cout << "performing Fit of variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;

  double fitMinABC2D = templatesAll["Ztautau"]["A"][getKey("diTauMt")]->GetXaxis()->GetXmin();
  double fitMaxABC2D = templatesAll["Ztautau"]["A"][getKey("diTauMt")]->GetXaxis()->GetXmax();
  RooRealVar* fitVarABC2D = new RooRealVar("fitVarABC2D", "fitVarABC2D", fitMinABC2D, fitMaxABC2D);

  double fitMinC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmin();
  double fitMaxC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmax();
  double fitMinC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmin();
  double fitMaxC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmax();
  assert(fitMinC1p == fitMinC1f && fitMaxC1p == fitMaxC1f);
  RooRealVar* fitVarC1 = new RooRealVar("fitVarC1", "fitVarC1", fitMinC1p, fitMaxC1p);
  
  double numEventsDataA = distributionsData["A"][getKey("diTauMt")]->Integral();
  double numEventsDataB = distributionsData["B"][getKey("diTauMt")]->Integral();
  double numEventsDataC = distributionsData["C"][getKey("diTauMt")]->Integral();
  double numEventsDataD = distributionsData["D"][getKey("diTauMt")]->Integral();
  double numEventsDataABCD = numEventsDataA + numEventsDataB + numEventsDataC + numEventsDataD;
  std::cout << "numEventsDataABCD = " << numEventsDataABCD << std::endl;

  std::map<std::string, RooRealVar*> pDiTauCharge_OS_SS;          // key = process
  std::map<std::string, RooRealVar*> pMuonIso_loose_tight;        // key = process
  std::map<std::string, RooRealVar*> pDiTauKine_Sig_Bgr;          // key = process
  std::map<std::string, RooRealVar*> pTauId_passed_failed;        // key = process

  std::map<std::string, RooRealVar*> normABCD;                    // key = process
  std::map<std::string, RooFormulaVar*> normA;                    // key = process
  std::map<std::string, RooFormulaVar*> normB;                    // key = process
  std::map<std::string, RooFormulaVar*> normC1p;                  // key = process
  std::map<std::string, RooFormulaVar*> normC1f;                  // key = process
  std::map<std::string, RooFormulaVar*> normC2p;                  // key = process
  std::map<std::string, RooFormulaVar*> normC2f;                  // key = process
  std::map<std::string, RooFormulaVar*> normD;                    // key = process
  
  std::map<std::string, std::map<std::string, RooHistPdf*> > pdf; // key = (process/"sum", region)

  TObjArray pdfsA;
  TObjArray pdfsB;
  TObjArray pdfsC1p;
  TObjArray pdfsC1f;
  TObjArray pdfsC2p;
  TObjArray pdfsC2f;
  TObjArray pdfsD;

  TObjArray fitParametersA;
  TObjArray fitParametersB;
  TObjArray fitParametersC1p;
  TObjArray fitParametersC1f;
  TObjArray fitParametersC2p;
  TObjArray fitParametersC2f;
  TObjArray fitParametersD;

  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    double numEventsA    = numEventsAll[*process]["A"][getKey("diTauMt")];
    double numEventsB    = numEventsAll[*process]["B"][getKey("diTauMt")];
    double numEventsC    = numEventsAll[*process]["C"][getKey("diTauMt")];
    double numEventsC1   = numEventsAll[*process]["C1"][getKey("diTauMt")];
    double numEventsC1p  = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double numEventsC2p  = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double numEventsABCD = numEventsAll[*process]["ABCD"][getKey("diTauMt")];
    std::cout << "numEventsABCD = " << numEventsABCD << std::endl;

    std::string nameDiTauCharge_OS_SS = std::string("pDiTauCharge_OS_SS").append("_").append(*process);
    pDiTauCharge_OS_SS[*process] = new RooRealVar(nameDiTauCharge_OS_SS.data(), nameDiTauCharge_OS_SS.data(), (numEventsA + numEventsC)/numEventsABCD, 0., 1.);
    std::string nameMuonIso_loose_tight = std::string("pMuonIso_loose_tight").append("_").append(*process);
    pMuonIso_loose_tight[*process] = new RooRealVar(nameMuonIso_loose_tight.data(), nameDiTauCharge_OS_SS.data(), (numEventsA + numEventsB)/numEventsABCD, 0., 1.);
    std::string nameDiTauKine_Sig_Bgr = std::string("pDiTauKine_Sig_Bgr").append("_").append(*process);
    pDiTauKine_Sig_Bgr[*process] = new RooRealVar(nameDiTauKine_Sig_Bgr.data(), nameDiTauKine_Sig_Bgr.data(), numEventsC1/numEventsC, 0., 1.);
    std::string nameTauId_passed_failed = std::string("pTauId_passed_failed").append("_").append(*process);
    pTauId_passed_failed[*process] = new RooRealVar(nameTauId_passed_failed.data(), nameTauId_passed_failed.data(), (numEventsC1p + numEventsC2p)/numEventsC, 0., 1.);

    double numEventsSumABCD = numEventsAll["sum"]["ABCD"][getKey("diTauMt")];
    std::cout << "numEventsSumABCD = " << numEventsSumABCD << std::endl;

    std::string nameNormABCD = std::string("normABCD").append("_").append(*process);
    normABCD[*process] = new RooRealVar(nameNormABCD.data(), nameNormABCD.data(), numEventsDataABCD*numEventsABCD/numEventsSumABCD, 0., numEventsDataABCD);

    TH1* templateA   = templatesAll[*process]["A"][getKey("diTauMt")];
    RooHistPdf* pdfA   = makeRooHistPdf(templateA, fitVarABC2D);
    TH1* templateB   = templatesAll[*process]["B"][getKey("diTauMt")];
    RooHistPdf* pdfB   = makeRooHistPdf(templateB, fitVarABC2D);
    TH1* templateC1p = templatesAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    RooHistPdf* pdfC1p = makeRooHistPdf(templateC1p, fitVarC1);
    TH1* templateC1f = templatesAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    RooHistPdf* pdfC1f = makeRooHistPdf(templateC1f, fitVarC1);
    TH1* templateC2p = templatesAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    RooHistPdf* pdfC2p = makeRooHistPdf(templateC2p, fitVarABC2D);
    TH1* templateC2f = templatesAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    RooHistPdf* pdfC2f = makeRooHistPdf(templateC2f, fitVarABC2D);
    TH1* templateD   = templatesAll[*process]["D"][getKey("diTauMt")];
    RooHistPdf* pdfD   = makeRooHistPdf(templateD, fitVarABC2D);

    pdfsA.Add(pdfA);
    pdfsB.Add(pdfB);
    pdfsC1p.Add(pdfC1p);
    pdfsC1f.Add(pdfC1f);
    pdfsC2p.Add(pdfC2p);
    pdfsC2f.Add(pdfC2f);
    pdfsD.Add(pdfD);

    normA[*process]   = makeRooFormulaVar(*process, "A", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normB[*process]   = makeRooFormulaVar(*process, "B", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC1p[*process] = makeRooFormulaVar(*process, "C1p", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC1f[*process] = makeRooFormulaVar(*process, "C1f", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC2p[*process] = makeRooFormulaVar(*process, "C2p", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normC2f[*process] = makeRooFormulaVar(*process, "C2f", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);
    normD[*process]   = makeRooFormulaVar(*process, "D", normABCD[*process], 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_loose_tight[*process], pTauId_passed_failed[*process]);

    fitParametersA.Add(normA[*process]);
    fitParametersB.Add(normB[*process]);
    fitParametersC1p.Add(normC1p[*process]);
    fitParametersC1f.Add(normC1f[*process]);
    fitParametersC2p.Add(normC2p[*process]);
    fitParametersC2f.Add(normC2f[*process]);
    fitParametersD.Add(normD[*process]);
  }

  RooAddPdf* pdfSumA   = new RooAddPdf("pdfSumA", "pdfSumB", RooArgList(pdfsA), RooArgList(fitParametersA));
  RooAddPdf* pdfSumB   = new RooAddPdf("pdfSumB", "pdfSumB", RooArgList(pdfsB), RooArgList(fitParametersB));
  RooAddPdf* pdfSumC1p = new RooAddPdf("pdfSumC1p", "pdfSumC1p", RooArgList(pdfsC1p), RooArgList(fitParametersC1p));
  RooAddPdf* pdfSumC1f = new RooAddPdf("pdfSumC1f", "pdfSumC1f", RooArgList(pdfsC1f), RooArgList(fitParametersC1f));
  RooAddPdf* pdfSumC2p = new RooAddPdf("pdfSumC2p", "pdfSumC2p", RooArgList(pdfsC2p), RooArgList(fitParametersC2p));
  RooAddPdf* pdfSumC2f = new RooAddPdf("pdfSumC2f", "pdfSumC2f", RooArgList(pdfsC2f), RooArgList(fitParametersC2f));
  RooAddPdf* pdfSumD   = new RooAddPdf("pdfSumD", "pdfSumD", RooArgList(pdfsD), RooArgList(fitParametersD));

// CV: due to limitation in RooFit
//    (cf. http://root.cern.ch/phpBB3/viewtopic.php?f=15&t=9518)
//     need to construct log-likelihood functions separately for regions { A, B, D } and { C1p, C1f }

//--- build data & model objects for fitting regions A, B, C2p, C2f, D
  RooCategory* fitCategoriesABC2D = new RooCategory("categoriesABC2D", "categoriesABC2D");
  fitCategoriesABC2D->defineType("A");
  fitCategoriesABC2D->defineType("B");
  fitCategoriesABC2D->defineType("C2p");
  fitCategoriesABC2D->defineType("C2f");
  fitCategoriesABC2D->defineType("D");

  RooSimultaneous* pdfSimultaneousFitABC2D = new RooSimultaneous("pdfSimultaneousFitABC2D", "pdfSimultaneousFitABC2D", *fitCategoriesABC2D);
  pdfSimultaneousFitABC2D->addPdf(*pdfSumA, "A");
  pdfSimultaneousFitABC2D->addPdf(*pdfSumB, "B");
  pdfSimultaneousFitABC2D->addPdf(*pdfSumC2p, "C2p");
  pdfSimultaneousFitABC2D->addPdf(*pdfSumC2f, "C2f");
  pdfSimultaneousFitABC2D->addPdf(*pdfSumD, "D");

  std::map<std::string, TH1*> histogramDataMapABC2D;
  histogramDataMapABC2D["A"]  = distributionsData["A"][getKey("diTauMt")];
  histogramDataMapABC2D["B"]  = distributionsData["B"][getKey("diTauMt")];
  histogramDataMapABC2D["C2p"]  = distributionsData["C2p"][getKey("diTauMt", tauId, "passed")];
  histogramDataMapABC2D["C2f"]  = distributionsData["C2f"][getKey("diTauMt", tauId, "failed")];
  histogramDataMapABC2D["D"]  = distributionsData["D"][getKey("diTauMt")];

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

//--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  TObjArray fitConstraintsABC2D;
  fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["QCD"], pDiTauCharge_OS_SS["QCD"]->getVal(), 0.1));
  fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["QCD"], pDiTauKine_Sig_Bgr["QCD"]->getVal(), 0.1));
  fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_loose_tight["QCD"], pMuonIso_loose_tight["QCD"]->getVal(), 0.1));
  fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["WplusJets"], pDiTauCharge_OS_SS["WplusJets"]->getVal(), 0.1));
  fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["WplusJets"], pDiTauKine_Sig_Bgr["WplusJets"]->getVal(), 0.1));
  fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_loose_tight["WplusJets"], pMuonIso_loose_tight["WplusJets"]->getVal(), 0.1));

  RooLinkedList fitOptionsABC2D;
  fitOptionsABC2D.Add(new RooCmdArg(RooFit::Extended()));
  fitOptionsABC2D.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsABC2D))));
  
  TObjArray fitConstraintsC1;
  fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["QCD"], pTauId_passed_failed["QCD"]->getVal(), 0.05));
  fitConstraintsC1.Add(makeFitConstraint(pTauId_passed_failed["WplusJets"], pTauId_passed_failed["WplusJets"]->getVal(), 0.05));

  RooLinkedList fitOptionsC1;
  fitOptionsC1.Add(new RooCmdArg(RooFit::Extended()));
  fitOptionsC1.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsC1))));

  //pdfSimultaneousFit->fitTo(*data, fitOptions);

  RooAbsReal* nllABC2D = pdfSimultaneousFitABC2D->createNLL(*dataABC2D, fitOptionsABC2D); 
  RooAbsReal* nllC1 = pdfSimultaneousFitC1->createNLL(*dataC1, fitOptionsC1); 
  RooAddition nll("nll", "nll", RooArgSet(*nllABC2D, *nllC1)); 
  RooMinuit minuit(nll); 
  minuit.setErrorLevel(1);
  minuit.setNoWarn();
  minuit.setPrintEvalErrors(1);
  minuit.setPrintLevel(0);
  //minuit.setWarnLevel(1);
  minuit.migrad(); 
  minuit.hesse(); 

  std::cout << "Results of fitting variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    double numEventsA    = numEventsAll[*process]["A"][getKey("diTauMt")];
    double numEventsB    = numEventsAll[*process]["B"][getKey("diTauMt")];
    double numEventsC    = numEventsAll[*process]["C"][getKey("diTauMt")];
    double numEventsC1   = numEventsAll[*process]["C1"][getKey("diTauMt")];
    double numEventsC1p  = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double numEventsC1f  = numEventsAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double numEventsC2   = numEventsAll[*process]["C2"][getKey("diTauMt")];
    double numEventsC2p  = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double numEventsC2f  = numEventsAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
    double numEventsD    = numEventsAll[*process]["D"][getKey("diTauMt")];
    double numEventsABCD = numEventsAll[*process]["ABCD"][getKey("diTauMt")];
    
    std::cout << " " << (*process) << ":" << std::endl;
    std::cout << "  normalization = " << normABCD[*process]->getVal() << " +/- " << normABCD[*process]->getError() 
	      << " (MC exp. = " << numEventsABCD << ")" << std::endl;
    std::cout << "  pDiTauCharge_OS_SS = " << pDiTauCharge_OS_SS[*process]->getVal() << " +/- " << pDiTauCharge_OS_SS[*process]->getError() 
	      << " (MC exp. = " << (numEventsA + numEventsC)/numEventsABCD << ")" << std::endl;
    std::cout << "  pMuonIso_loose_tight = " << pMuonIso_loose_tight[*process]->getVal() << " +/- " << pMuonIso_loose_tight[*process]->getError() 
	      << " (MC exp. = " << (numEventsA + numEventsB)/numEventsABCD << ")" << std::endl;
    std::cout << "  pDiTauKine_Sig_Bgr = " << pDiTauKine_Sig_Bgr[*process]->getVal() << " +/- " << pDiTauKine_Sig_Bgr[*process]->getError() 
	      << " (MC exp. = " << numEventsC1/numEventsC << ")" << std::endl;
    std::cout << "  pTauId_passed_failed = " << pTauId_passed_failed[*process]->getVal() << " +/- " << pTauId_passed_failed[*process]->getError() 
	      << " (MC exp. = " << (numEventsC1p + numEventsC2p)/numEventsC << ")" << std::endl;
    std::cout << "--> A = " << normA[*process]->getVal() 
	      << " (MC exp. = " << numEventsA << ")" << std::endl;
    std::cout << "--> B = " << normB[*process]->getVal()  
	      << " (MC exp. = " << numEventsB << ")" << std::endl;
    std::cout << "--> C = " << normC1p[*process]->getVal() + normC1f[*process]->getVal() + normC2p[*process]->getVal() + normC2f[*process]->getVal() 
	      << " (MC exp. = " << numEventsC << ")" << std::endl;
    std::cout << "--> C1 = " << normC1p[*process]->getVal() + normC1f[*process]->getVal()  
	      << " (MC exp. = " << numEventsC1 << ")" << std::endl;
    std::cout << "--> C1p = " << normC1p[*process]->getVal()  
	      << " (MC exp. = " << numEventsC1p << ")" << std::endl;
    std::cout << "--> C1f = " << normC1f[*process]->getVal()
	      << " (MC exp. = " << numEventsC1f << ")" << std::endl;
    std::cout << "--> C2 = " << normC2p[*process]->getVal() + normC2f[*process]->getVal() 
	      << " (MC exp. = " << numEventsC2 << ")" << std::endl;
    std::cout << "--> C2p = " << normC2p[*process]->getVal() 
	      << " (MC exp. = " << numEventsC2p << ")" << std::endl;
    std::cout << "--> C2f = " << normC2f[*process]->getVal()  
	      << " (MC exp. = " << numEventsC2f << ")" << std::endl;
    std::cout << "--> D = " << normD[*process]->getVal()  
	      << " (MC exp. = " << numEventsD << ")" << std::endl;
  }

  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();

//--- make control plots for sum(MC) scaled by normalization determined by fit versus Data 
//    for Mt, fitVariable distributions in different regions
  drawHistograms(templatesAll["QCD"]["A"][getKey("diTauMt")], normA["QCD"]->getVal(),
		 templatesAll["WplusJets"]["A"][getKey("diTauMt")], normA["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["A"][getKey("diTauMt")], normA["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["A"][getKey("diTauMt")], normA["Ztautau"]->getVal(),
		 distributionsData["A"][getKey("diTauMt")],
		 std::string("Region A: M_{T} (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_A_Mt_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["B"][getKey("diTauMt")], normB["QCD"]->getVal(),
		 templatesAll["WplusJets"]["B"][getKey("diTauMt")], normB["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["B"][getKey("diTauMt")], normB["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["B"][getKey("diTauMt")], normB["Ztautau"]->getVal(),
		 distributionsData["B"][getKey("diTauMt")],
		 std::string("Region B: M_{T} (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_B_Mt_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["C1p"][getKey(fitVariable, tauId, "passed")], normC1p["QCD"]->getVal(),
		 templatesAll["WplusJets"]["C1p"][getKey(fitVariable, tauId, "passed")], normC1p["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["C1p"][getKey(fitVariable, tauId, "passed")], normC1p["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")], normC1p["Ztautau"]->getVal(),
		 distributionsData["C1p"][getKey(fitVariable, tauId, "passed")],
		 std::string("Region C1p: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_C1p_").append(fitVariable).append("_").append(tauId).append("_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["C1f"][getKey(fitVariable, tauId, "failed")], normC1f["QCD"]->getVal(),
		 templatesAll["WplusJets"]["C1f"][getKey(fitVariable, tauId, "failed")], normC1f["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["C1f"][getKey(fitVariable, tauId, "failed")], normC1f["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")], normC1f["Ztautau"]->getVal(),
		 distributionsData["C1f"][getKey(fitVariable, tauId, "failed")],
		 std::string("Region C1f: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_C1f_").append(fitVariable).append("_").append(tauId).append("_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["C2p"][getKey("diTauMt", tauId, "passed")], normA["QCD"]->getVal(),
		 templatesAll["WplusJets"]["C2p"][getKey("diTauMt", tauId, "passed")], normA["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["C2p"][getKey("diTauMt", tauId, "passed")], normA["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["C2p"][getKey("diTauMt", tauId, "passed")], normA["Ztautau"]->getVal(),
		 distributionsData["C2p"][getKey("diTauMt", tauId, "passed")],
		 std::string("Region C2p: M_{T}, ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_C2p_Mt_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["C2f"][getKey("diTauMt", tauId, "failed")], normB["QCD"]->getVal(),
		 templatesAll["WplusJets"]["C2f"][getKey("diTauMt", tauId, "failed")], normB["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["C2f"][getKey("diTauMt", tauId, "failed")], normB["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["C2f"][getKey("diTauMt", tauId, "failed")], normB["Ztautau"]->getVal(),
		 distributionsData["C2f"][getKey("diTauMt", tauId, "failed")],
		 std::string("Region C2f: M_{T}, ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_C2f_Mt_fitted.png"), sysShift);
  drawHistograms(templatesAll["QCD"]["D"][getKey("diTauMt")], normD["QCD"]->getVal(),
		 templatesAll["WplusJets"]["D"][getKey("diTauMt")], normD["WplusJets"]->getVal(),
		 templatesAll["Zmumu"]["D"][getKey("diTauMt")], normD["Zmumu"]->getVal(),
		 templatesAll["Ztautau"]["D"][getKey("diTauMt")], normD["Ztautau"]->getVal(),
		 distributionsData["D"][getKey("diTauMt")],
		 std::string("Region D: M_{T} (scaled by normalization det. by fit)"),
		 std::string("controlPlotsTauIdEff_D_Mt_fitted.png"), sysShift);
}

std::string getFileNames(const std::string& inputFilePath, const std::string& sampleName, const std::string& jobId)
{
//--------------------------------------------------------------------------------
// Compose full fileName from inputFilePath, sampleName, jobId
//--------------------------------------------------------------------------------

  return std::string(inputFilePath).append("tauIdEffMeasEDNtuple_").append(sampleName).append("_").append(jobId).append("*.root");
}

void fitTauIdEff_wConstraints()
{
  gROOT->SetBatch(true);

//--- keep track of time it took the macro to execute
  TBenchmark clock;
  clock.Start("fitTauIdEff_wConstraints");

  //const std::string inputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas/2011Jan06_lxbatch/";
  const std::string inputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/2011Jan08_lxbatch/";

  const std::string jobId = "2011Jan08_lxbatch";

  //const std::string branchName_suffix = "local";
  const std::string branchName_suffix = "lxbatch";

  //bool runClosureTest = false;
  bool runClosureTest = true;

  bool runQuickTest = false;
  //bool runQuickTest = true;

  std::vector<std::string> sysUncertainties;
  sysUncertainties.push_back(std::string("SysTauJetEnUp"));
  sysUncertainties.push_back(std::string("SysTauJetEnDown"));
  sysUncertainties.push_back(std::string("SysJetEnUp"));
  sysUncertainties.push_back(std::string("SysJetEnDown"));
  sysUncertainties.push_back(std::string("SysZllRecoilCorrectionUp"));
  sysUncertainties.push_back(std::string("SysZllRecoilCorrectionDown"));
  //bool runSysUncertainties = false;
  bool runSysUncertainties = true;

  double dataIntLumi = 36.2;

  double weightFactorZtautau = dataIntLumi*1666/1994719;     // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  double corrFactorZtautau = 1.0;
  double weightFactorZmumu = dataIntLumi*1666/1998931;       // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  double corrFactorZmumu = 1.0;
  double weightFactorQCD = dataIntLumi*0.2966*1.e+9*2.855e-4/29504866; // xSection (LO) / numEvents (PYTHIA PPmuXptGt20Mu15 sample)
  double corrFactorQCD = 1.0;
  double weightFactorWplusJets = dataIntLumi*31314/15168266; // W --> l nu xSection (FEWZ @ NNLO) / numEvents (MadGraph sample)
  double corrFactorWplusJets = 1.0;
  //double weightFactorTTplusJets = dataIntLumi*157.5/1164640; // inclusive TTbar xSection (MCFM @ NLO) / numEvents (MadGraph sample)
  //double corrFactorTTplusJets = 1.0;

  std::vector<std::string> tauIds;
  tauIds.push_back(std::string("tauDiscrTaNCloose"));
  //tauIds.push_back(std::string("tauDiscrTaNCmedium"));
  //tauIds.push_back(std::string("tauDiscrTaNCtight"));
  //tauIds.push_back(std::string("tauDiscrHPSloose"));
  //tauIds.push_back(std::string("tauDiscrHPSmedium"));
  //tauIds.push_back(std::string("tauDiscrHPStight"));

  std::vector<std::string> fitVariables;
  //fitVariables.push_back("diTauHt");
  //fitVariables.push_back("diTauSVfitMass1");
  //fitVariables.push_back("diTauSVfitMass2");
  fitVariables.push_back("diTauVisMass");
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

  for ( std::vector<std::string>::const_iterator sysShift = sysShifts.begin();
	sysShift != sysShifts.end(); ++sysShift ) {
    std::cout << "running fit for sysShift = " << (*sysShift) << "..." << std::endl;

//--- initialize alias --> branchName mapping
    std::map<std::string, std::string> branchNamesData = makeBranchNameDict("CENTRAL_VALUE", branchName_suffix);
    std::map<std::string, std::string> branchNamesMC   = makeBranchNameDict(*sysShift,       branchName_suffix);

    TChain* chainData = new TChain("Events");
    if ( !runQuickTest ) { chainData->Add(getFileNames(inputFilePath, "data_Mu_Run2010A_Nov4ReReco", jobId).data());
                           chainData->Add(getFileNames(inputFilePath, "data_Mu_Run2010B_Nov4ReReco", jobId).data()); } 
    else                 { chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010A_Nov4ReReco_2011Jan08_lxbatch_0_3b65.root").data());
                           chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010B_Nov4ReReco_2011Jan08_lxbatch_0_359e.root").data()); }
    std::cout << "chainData has " << chainData->GetListOfFiles()->GetEntries() << " files, containing " << chainData->GetEntries() << " events." << std::endl;
    std::map<std::string, std::map<std::string, TH1*> > distributionsData = // key = (region, observable)
      makeDistributionsAllRegions("Data", 1.0, chainData, 
				  tauIds, fitVariables, branchNamesData, *sysShift);

    TChain* chainZtautau = new TChain("Events");
    if ( !runQuickTest ) chainZtautau->Add(getFileNames(inputFilePath, "Ztautau_powheg", jobId).data());
    else                 chainZtautau->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Ztautau_powheg_2011Jan08_lxbatch_0_97ae.root").data());
    std::cout << "chainZtautau has " << chainZtautau->GetListOfFiles()->GetEntries() << " files, containing " << chainZtautau->GetEntries() << " events." << std::endl;
    std::map<std::string, std::map<std::string, TH1*> > templatesZtautau = // key = (region, observable)
      makeDistributionsAllRegions("Ztautau", weightFactorZtautau*corrFactorZtautau, chainZtautau, 
				  tauIds, fitVariables, branchNamesMC, *sysShift);
    
    TChain* chainZmumu = new TChain("Events");
    if ( !runQuickTest ) chainZmumu->Add(getFileNames(inputFilePath, "Zmumu_powheg", jobId).data());
    else                 chainZmumu->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Zmumu_powheg_2011Jan08_lxbatch_0_ea01.root").data());
    std::cout << "chainZmumu has " << chainZmumu->GetListOfFiles()->GetEntries() << " files, containing " << chainZmumu->GetEntries() << " events." << std::endl;
    std::map<std::string, std::map<std::string, TH1*> > templatesZmumu = // key = (region, observable)
      makeDistributionsAllRegions("Zmumu", weightFactorZmumu*corrFactorZmumu, chainZmumu, 
				  tauIds, fitVariables, branchNamesMC, *sysShift);
    
    TChain* chainQCD = new TChain("Events");
    if ( !runQuickTest ) chainQCD->Add(getFileNames(inputFilePath, "PPmuXptGt20Mu15", jobId).data());
    else                 chainQCD->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jan08_lxbatch_0_6255.root").data());
    std::cout << "chainQCD has " << chainQCD->GetListOfFiles()->GetEntries() << " files, containing " << chainQCD->GetEntries() << " events." << std::endl;
    std::map<std::string, std::map<std::string, TH1*> > templatesQCD = // key = (region, observable)
      makeDistributionsAllRegions("QCD", weightFactorQCD*corrFactorQCD, chainQCD, 
				  tauIds, fitVariables, branchNamesMC, *sysShift);
    
    TChain* chainWplusJets = new TChain("Events");
    if ( !runQuickTest ) chainWplusJets->Add(getFileNames(inputFilePath, "WplusJets_madgraph", jobId).data());
    else                 chainWplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jan08_lxbatch_0_da28.root").data());
    std::cout << "chainWplusJets has " << chainWplusJets->GetListOfFiles()->GetEntries() << " files, containing " << chainWplusJets->GetEntries() << " events." << std::endl;
    std::map<std::string, std::map<std::string, TH1*> > templatesWplusJets = // key = (region, observable)
      makeDistributionsAllRegions("WplusJets", weightFactorWplusJets*corrFactorWplusJets, chainWplusJets, 
				  tauIds, fitVariables, branchNamesMC, *sysShift);
    
    //TChain* chainTTplusJets = new TChain("Events");
    //if ( !runQuickTest ) chainTTplusJets->Add(getFileNames(inputFilePath, "TTplusJets_madgraph", jobId).data());
    //else                 chainTTplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jan08_lxbatch_0_f075.root").data());
    //std::cout << "chainTTplusJets has " << chainTTplusJets->GetListOfFiles()->GetEntries() << " files, containing " << chainTTplusJets->GetEntries() << " events." << std::endl;
    //std::map<std::string, std::map<std::string, TH1*> > templatesTTplusJets = // key = (region, observable)
    //  makeDistributionsAllRegions("TTplusJets", weightFactorTTplusJets*corrFactorTTplusJets, chainTTplusJets, 
    //                              tauIds, fitVariables, branchNamesMC, *sysShift);
    
    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > > templatesAll; // key = (process, region, observable)
    templatesAll["Ztautau"] = templatesZtautau;
    templatesAll["Zmumu"] = templatesZmumu;
    templatesAll["QCD"] = templatesQCD;
    templatesAll["WplusJets"] = templatesWplusJets;
    //templatesAll["TTplusJets"] = templatesTTplusJets;
    
    std::vector<std::string> processes;
    processes.push_back(std::string("Ztautau"));
    processes.push_back(std::string("Zmumu"));
    processes.push_back(std::string("QCD"));
    processes.push_back(std::string("WplusJets"));
    //processes.push_back(std::string("TTplusJets"));
    
//--- closure test: fit sum(MC) instead of Data
    if ( runClosureTest ) {
      std::cout << "NOTE: RUNNING CLOSURE TEST !!" << std::endl;
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	      region != distributionsData.end(); ++region ) {
	  for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
		key != region->second.end(); ++key ) {
	    TH1* histogram = templatesAll[*process][region->first][key->first];
	    std::cout << " histogram = " << histogram->GetName() << ": integral = " << histogram->Integral() << std::endl;
	    
	    TH1* histogramSum = templatesAll["sum"][region->first][key->first];
	    if ( !histogramSum ) {
	      std::string histogramName = histogram->GetName();
	      std::string histogramSumName = std::string("sum").append(std::string(histogramName, histogramName.find("_")));	    
	      templatesAll["sum"][region->first][key->first] = (TH1*)histogram->Clone(histogramSumName.data());
	      std::cout << "--> creating new sum(MC) histogram = " << histogramSumName << ": integral = " << templatesAll["sum"][region->first][key->first]->Integral() << std::endl;
	    } else {	    
	      histogramSum->Add(histogram);
	      std::cout << "--> adding histogram to sum(MC) histogram = " << histogramSum->GetName() << ": integral = " << histogramSum->Integral() << std::endl;
	    }
	  }
	}
      }
      
      std::cout << std::endl;
      
      distributionsData = templatesAll["sum"];
    }

//--- make control plots for sum(MC) scaled by cross-sections versus Data 
//    for Mt, fitVariable distributions in different regions
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	drawHistograms(templatesQCD[region->first][key->first], -1.,
		       templatesWplusJets[region->first][key->first], -1.,
		       templatesZmumu[region->first][key->first], -1.,
		       templatesZtautau[region->first][key->first], -1.,
		       distributionsData[region->first][key->first],
		       std::string("Region ").append(region->first).append(": ").append(key->first).append(" (scaled by cross-section)"),
		       std::string("controlPlotsTauIdEff_").append(region->first).append("_").append(key->first).append(".png"), *sysShift);
      }
    }
    
//--- print MC expectations for probabilities 
//   o pDiTauCharge_OS_SS
//   o pDiTauKine_Sig_Bgr
//   o pMuonIso_loose_tight
//   o pTauId_passed_failed
//    separating different regions

    std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll; // key = (process/"sum", region, observable)
    for ( std::vector<std::string>::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	    region != distributionsData.end(); ++region ) {
	for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	      key != region->second.end(); ++key ) {
	  numEventsAll[*process][region->first][key->first] = templatesAll[*process][region->first][key->first]->Integral();
	  numEventsAll["sum"][region->first][key->first] += numEventsAll[*process][region->first][key->first];
	}
      }
      
      numEventsAll[*process]["ABCD"][getKey("diTauMt")] = numEventsAll[*process]["A"][getKey("diTauMt")] + numEventsAll[*process]["B"][getKey("diTauMt")] 
	                                                 + numEventsAll[*process]["C"][getKey("diTauMt")] + numEventsAll[*process]["D"][getKey("diTauMt")];
      numEventsAll["sum"]["ABCD"][getKey("diTauMt")] += numEventsAll[*process]["ABCD"][getKey("diTauMt")];
    }
    
    for ( std::vector<std::string>::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      double numEventsA  = numEventsAll[*process]["A"][getKey("diTauMt")];
      double numEventsB  = numEventsAll[*process]["B"][getKey("diTauMt")];
      double numEventsC  = numEventsAll[*process]["C"][getKey("diTauMt")];
      double numEventsC1 = numEventsAll[*process]["C1"][getKey("diTauMt")];
      double numEventsC2 = numEventsAll[*process]["C2"][getKey("diTauMt")];
      double numEventsD  = numEventsAll[*process]["D"][getKey("diTauMt")];
      
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
      std::cout << " D/(B+D) = " << numEventsC/(numEventsB + numEventsD) << std::endl;
      std::cout << "pDiTauKine_Sig_Bgr:" << std::endl;
      std::cout << " C1/C = " << numEventsC1/numEventsC << std::endl;
      
      std::cout << "pTauId_passed_failed:" << std::endl;
      for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	    tauId != tauIds.end(); ++tauId ) {
	double numEventsC1p = numEventsAll[*process]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
	double numEventsC2p = numEventsAll[*process]["C2p"][getKey("diTauMt", *tauId, "passed")];
	
	std::cout << " " << (*tauId) << ":" << std::endl;
	std::cout << "  C1p/C1 = " << numEventsC1p/numEventsC1 << std::endl;
	std::cout << "  C2p/C2 = " << numEventsC2p/numEventsC2 << std::endl;
      }
      
      std::cout << std::endl;
    }

    std::map<std::string, std::map<std::string, double> > effValues; // key = (tauId, fitVariable)
    std::map<std::string, std::map<std::string, double> > effErrors; // key = (tauId, fitVariable)

    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	double effValue = 0.;
	double effError = 1.;
	fitUsingRooFit(distributionsData, templatesAll, numEventsAll,
		       processes,
		       *tauId, *fitVariable,
		       effValue, effError,
		       *sysShift);
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
	
	double numEventsC = numEventsAll["Ztautau"]["C"][getKey("diTauMt")];
	double numEventsC1p = numEventsAll["Ztautau"]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
	double numEventsC2p = numEventsAll["Ztautau"]["C2p"][getKey("diTauMt", *tauId, "passed")];
	std::cout << "(Monte Carlo prediction = " << ((numEventsC1p + numEventsC2p)/numEventsC)*100. << "%)" << std::endl;
      }
    }

    delete chainData;
    delete chainZtautau;
    delete chainZmumu;
    delete chainQCD;
    delete chainWplusJets;
    //delete chainTTplusJets;
  }

//--print time that it took macro to run
  std::cout << "finished executing fitTauIdEff_wConstraints macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  std::cout << " #sysShifts    = " << sysShifts.size() << std::endl;
  clock.Show("fitTauIdEff_wConstraints");
}
