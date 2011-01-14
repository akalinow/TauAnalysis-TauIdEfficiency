
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

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

std::map<std::string, TH1*> makeHistograms(const std::string& process, double weight,
					   TTree* tree, const std::string& treeSelection, 
					   const std::vector<std::string>& tauIds, const std::vector<std::string>& tauIdValues,
					   const std::vector<std::string>& variables)
{
  std::map<std::string, TH1*> retVal;
  
  const std::string ntupleName = "doubles_ntupleProducer_tauIdEffNtuple";

  const std::string branchNameDiTau = std::string(ntupleName).append("#selectedMuTauPairsPzetaDiffCumulative");
  const std::string branchNameTau = std::string(ntupleName).append("#selectedPatTausForMuTauEcalCrackVetoCumulative");

  std::map<std::string, std::string> branchNames;
  
  branchNames["Ht"] = std::string(branchNameDiTau).append("#Ht_local.obj");
  branchNames["SVfitMass1"] = std::string(branchNameDiTau).append("#SVfitMass1_local.obj");
  branchNames["SVfitMass2"] = std::string(branchNameDiTau).append("#SVfitMass2_local.obj");
  branchNames["visMass"] = std::string(branchNameDiTau).append("#visMass_local.obj");
  
  branchNames["muonPt"] = std::string(branchNameDiTau).append("#muonPt_local.obj");
  branchNames["muonEta"] = std::string(branchNameDiTau).append("#muonEta_local.obj");
  branchNames["tauPt"] = std::string(branchNameDiTau).append("#tauPt_local.obj");
  branchNames["tauEta"] = std::string(branchNameDiTau).append("#tauEta_local.obj");

  branchNames["tauLooseIsolationPt"] = std::string(branchNameDiTau).append("#tauLooseIsolationPt_local.obj");
  
  branchNames["numChargedParticles"] = std::string(branchNameTau).append("#numChargedParticles_local.obj");
  branchNames["numParticles"] = std::string(branchNameTau).append("#numParticles_local.obj");
  branchNames["jetWidth"] = std::string(branchNameTau).append("#jetWidth_local.obj");
  
  branchNames["discrHPSloose"] = std::string(branchNameTau).append("#byHPSloose_local.obj");
  branchNames["discrHPSmedium"] = std::string(branchNameTau).append("#byHPSmedium_local.obj");
  branchNames["discrHPStight"] = std::string(branchNameTau).append("#byHPStight_local.obj");
  branchNames["discrTaNCloose"] = std::string(branchNameTau).append("#byTaNCloose_local.obj");
  branchNames["discrTaNCmedium"] = std::string(branchNameTau).append("#byTaNCmedium_local.obj");
  branchNames["discrTaNCtight"] = std::string(branchNameTau).append("#byTaNCtight_local.obj");

  std::string tmpTreeSelection = std::string(branchNames["muonPt"]).append(" > 20.");
  tmpTreeSelection.append(" && ").append(branchNames["tauLooseIsolationPt"]).append(" < 2.5");

  TTree* selEventsTree = 0;
  if ( tmpTreeSelection != "" ) {
    selEventsTree = tree->CopyTree(tmpTreeSelection.data());
  } else {
    selEventsTree = tree;
  }

  //selEventsTree->Print();

  std::cout << "process = " << process << " has " << selEventsTree->GetEntries() << " entries." << std::endl;
  std::cout << " (weighted sum = " << selEventsTree->GetEntries()*weight << ")" << std::endl;

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	  tauIdValue != tauIdValues.end(); ++tauIdValue ) {

      std::string treeSelectionForTauIdValue;
      if      ( (*tauIdValue) == "passed" ) treeSelectionForTauIdValue = branchNames[*tauId].append(" > 0.5");
      else if ( (*tauIdValue) == "failed" ) treeSelectionForTauIdValue = branchNames[*tauId].append(" < 0.5");
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

      for ( std::vector<std::string>::const_iterator variable = variables.begin();
	    variable != variables.end(); ++variable ) {
	
	std::string histogramName = std::string(process).append("_").append(*variable);
	histogramName.append("_").append(*tauId).append("_").append(*tauIdValue);

	int numBins;
	double min, max;

	if ( (*variable) == "Ht"         || 
	     (*variable) == "SVfitMass1" || 
	     (*variable) == "SVfitMass2" || 
	     (*variable) == "visMass"   ) {
	  numBins = 18;
	  min = 20.;
	  max = 200.;
	} else if ( (*variable) == "muonPt" ||
		    (*variable) == "tauPt" ) {
	  numBins = 20;
	  min =   0.;
	  max = 100.;
	} else if ( (*variable) == "muonEta" ||
		    (*variable) == "tauEta" ) {
	  numBins = 25;
	  min = -2.5;
	  max = +2.5;
	} else if ( (*variable) == "numChargedParticles" ||
		    (*variable) == "numParticles"       ) {
	  numBins = 25;
	  min =  -0.5;
	  max = +24.5;
	} else if ( (*variable) == "jetWidth" ) {
	  numBins = 20;
	  min = 0.;
	  max = 0.50;
	} else {
	  std::cout << "Error in <makeHistograms>: undefined variable = " << (*variable) << " --> skipping !!";
	  return retVal;
	}

	TH1* histogram = new TH1F(histogramName.data(), histogramName.data(), numBins, min, max);

        std::string drawCommand = std::string(branchNames[*variable]).append(">>+").append(histogramName); 
	selEventsTreeForTauIdValue->Draw(drawCommand.data());

	if ( !histogram->GetSumw2N() ) histogram->Sumw2();
	histogram->Scale(weight);

	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);

	if ( histogram != 0 ) retVal[key] = histogram;

	// CV: deleting selEventsTreeForTauIdValue causes segmentation violation ?!
        //if ( selEventsTreeForTauIdValue != selEventsTree ) delete selEventsTreeForTauIdValue;
      }
    }
  }

  if ( selEventsTree != tree ) delete selEventsTree;

  return retVal;
}

TH1* normalize(const TH1* histogram, double norm = 1.)
{
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
		    const std::string& outputFileName)
{
  //std::cout << "<drawHistograms>:" << std::endl;
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

  //smSum.SetTitle(templateZtautau->GetTitle());
  smSum.SetTitle("");
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
  canvas->Print(outputFilePath.append(outputFileName).data());

  if ( templateQCD       != histogramQCD       ) delete templateQCD;
  if ( templateWplusJets != histogramWplusJets ) delete templateWplusJets;
  if ( templateZmumu     != histogramZmumu     ) delete templateZmumu;
  if ( templateZtautau   != histogramZtautau   ) delete templateZtautau;

  delete canvas;
}

void drawHistograms(TH1* histogram_passed, TH1* histogram_failed, 
		    const std::string& tauId,
		    const std::string& outputFileName)
{
  //std::cout << "<drawHistograms>:" << std::endl;
  
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

  template_passed->SetTitle("");
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
  canvas->Print(outputFilePath.append(outputFileName).data());

  delete template_passed;
  delete template_failed;

  delete canvas;
}

void fitUsingTFractionFitter(std::map<std::string, TH1*>& histogramsData,
			     std::map<std::string, TH1*>& histogramsZtautau,
			     std::map<std::string, TH1*>& histogramsZmumu,
			     std::map<std::string, TH1*>& histogramsWplusJets,
			     std::map<std::string, TH1*>& histogramsQCD,
			     const std::vector<std::string>& tauIds, const std::vector<std::string>& tauIdValues,
			     const std::vector<std::string>& fitVariables,
			     std::map<std::string, double>& effValues,
			     std::map<std::string, double>& effErrors)
{
  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {

      double passedValueZtautau = 0.;
      double passedErrorZtautau = 0.;
      double failedValueZtautau = 0.;
      double failedErrorZtautau = 0.;

      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {

	std::cout << "performing Fit of variable = " << (*variable)
		  << " for Tau id. = " << (*tauId) << ", value = '"<< (*tauIdValue) << "'" << std::endl;
	
	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);
     
	TObjArray templates;
	templates.SetOwner(true);
	templates.Add(normalize(histogramsQCD[key]));
	templates.Add(normalize(histogramsWplusJets[key]));
	templates.Add(normalize(histogramsZmumu[key]));
	templates.Add(normalize(histogramsZtautau[key]));
	
	TH1* histogramData = histogramsData[key];
	
	TFractionFitter* fit = new TFractionFitter(histogramData, &templates);
	for ( int iProcess = 0; iProcess < templates.GetEntries(); ++iProcess ) {
	  fit->Constrain(iProcess, 0., 1.);
	}
	Int_t fitStatus = fit->Fit(); 
	cout << "fit status = " << fitStatus << endl;
	if ( fitStatus == 0 ) { // check that fit converged

	  TVirtualFitter* fitAlgorithm = fit->GetFitter();

	  double fitValueQCD = fitAlgorithm->GetParameter(0)*histogramData->Integral();
	  double fitErrorQCD = fitAlgorithm->GetParError(0)*histogramData->Integral();
	  double fitValueWplusJets = fitAlgorithm->GetParameter(1)*histogramData->Integral();
	  double fitErrorWplusJets = fitAlgorithm->GetParError(1)*histogramData->Integral();
	  double fitValueZmumu = fitAlgorithm->GetParameter(2)*histogramData->Integral();
	  double fitErrorZmumu = fitAlgorithm->GetParError(2)*histogramData->Integral();
	  double fitValueZtautau = fitAlgorithm->GetParameter(3)*histogramData->Integral();
	  double fitErrorZtautau = fitAlgorithm->GetParError(3)*histogramData->Integral();

	  std::cout << "Results of fitting variable = " << (*variable)
		    << " for Tau id. = " << (*tauId) << ", value = '"<< (*tauIdValue) << "':" << std::endl;
	  std::cout << " Ztautau: normalization = " << fitValueZtautau << " +/- " << fitErrorZtautau << std::endl;
	  std::cout << " Zmumu: normalization = " << fitValueZmumu << " +/- " << fitErrorZmumu << std::endl; 
	  std::cout << " WplusJets: normalization = " << fitValueWplusJets << " +/- " << fitErrorWplusJets << std::endl;
	  std::cout << " QCD: normalization = " << fitValueQCD << " +/- " << fitErrorQCD << std::endl; 

	  if ( (*tauIdValue) == "passed" ) {
	    passedValueZtautau = fitValueZtautau;
	    passedErrorZtautau = fitErrorZtautau;
	  } else if ( (*tauIdValue) == "failed" ) {
	    failedValueZtautau = fitValueZtautau;
	    failedErrorZtautau = fitErrorZtautau;
	  }
	}
      }

      double numerator = passedValueZtautau;
      double denominator = passedValueZtautau + failedValueZtautau;

      double effValue = numerator/denominator;
      double effError2 = TMath::Power((failedValueZtautau/denominator)*(passedErrorZtautau/denominator), 2)
	                + TMath::Power((passedValueZtautau/denominator)*(failedErrorZtautau/denominator), 2);

      std::string key = std::string(*variable);
      key.append("_").append(*tauId);

      effValues[key] = effValue;
      effErrors[key] = TMath::Sqrt(effError2);
    }
  }
}

void fitUsingRooFit(std::map<std::string, TH1*>& histogramsData,
		    std::map<std::string, TH1*>& histogramsZtautau,
		    std::map<std::string, TH1*>& histogramsZmumu,
		    std::map<std::string, TH1*>& histogramsWplusJets,
		    std::map<std::string, TH1*>& histogramsQCD,
		    const std::vector<std::string>& tauIds, const std::vector<std::string>& tauIdValues,
		    const std::vector<std::string>& fitVariables,
		    std::map<std::string, double>& effValues,
		    std::map<std::string, double>& effErrors)
{
  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      
      std::cout << "performing Fit of variable = " << (*variable)
		<< " for Tau id. = " << (*tauId) << std::endl;
            
      std::string key = std::string(*variable);
      key.append("_").append(*tauId);

      double fitMin = histogramsZtautau[std::string(key).append("_").append("all")]->GetXaxis()->GetXmin();
      double fitMax = histogramsZtautau[std::string(key).append("_").append("all")]->GetXaxis()->GetXmax();
      
      RooRealVar* fitVar = new RooRealVar("fitVar", "fitVar", fitMin, fitMax);

      double numEventsData = histogramsData[std::string(key).append("_").append("all")]->Integral();

      RooRealVar* fr = new RooRealVar("fr", "fr", 0.05, 0., 1.);

      TH1* histogramQCD_passed = histogramsQCD[std::string(key).append("_").append("passed")];
      RooDataHist* tmpHistQCD_passed = 
	new RooDataHist("tmpHistQCD_passed", "tmpHistQCD_passed", *fitVar, histogramQCD_passed);
      RooHistPdf* templateQCD_passed = 
	new RooHistPdf("templateQCD_passed", "templateQCD_passed", *fitVar, *tmpHistQCD_passed);
      TH1* histogramQCD_failed = histogramsQCD[std::string(key).append("_").append("failed")];
      RooDataHist* tmpHistQCD_failed = 
	new RooDataHist("tmpHistQCD_failed", "tmpHistQCD_failed", *fitVar, histogramQCD_failed);
      RooHistPdf* templateQCD_failed = 
	new RooHistPdf("templateQCD_failed", "templateQCD_failed", *fitVar, *tmpHistQCD_failed);
      RooRealVar* normQCD = new RooRealVar("normQCD", "normQCD", 0.40*numEventsData, 0., numEventsData);
      RooRealVar* frFactorQCD = new RooRealVar("frFactorQCD", "frFactorQCD", 1., 0., 10.);
      RooFormulaVar* normQCD_passed = 
	new RooFormulaVar("normQCD_passed", "normQCD_passed", 
			  "fr*frFactorQCD*normQCD", RooArgSet(*normQCD, *fr, *frFactorQCD));
      RooFormulaVar* normQCD_failed = 
	new RooFormulaVar("normQCD_failed", "normQCD_failed", 
			  "(1.0 - fr*frFactorQCD)*normQCD", RooArgSet(*normQCD, *fr, *frFactorQCD));

      TH1* histogramWplusJets_passed = histogramsWplusJets[std::string(key).append("_").append("passed")];
      RooDataHist* tmpHistWplusJets_passed =
	new RooDataHist("tmpHistWplusJets_passed", "tmpHistWplusJets_passed", *fitVar, histogramWplusJets_passed);
      RooHistPdf* templateWplusJets_passed =
	new RooHistPdf("templateWplusJets_passed", "templateWplusJets_passed", *fitVar, *tmpHistWplusJets_passed);
      TH1* histogramWplusJets_failed = histogramsWplusJets[std::string(key).append("_").append("failed")];
      RooDataHist* tmpHistWplusJets_failed =
	new RooDataHist("tmpHistWplusJets_failed", "tmpHistWplusJets_failed", *fitVar, histogramWplusJets_failed);
      RooHistPdf* templateWplusJets_failed =
	new RooHistPdf("templateWplusJets_failed", "templateWplusJets_failed", *fitVar, *tmpHistWplusJets_failed);
      RooRealVar* normWplusJets = new RooRealVar("normWplusJets", "normWplusJets", 0.25*numEventsData, 0., numEventsData);
      RooRealVar* frFactorWplusJets = new RooRealVar("frFactorWplusJets", "frFactorWplusJets", 1., 0., 10.);
      RooFormulaVar* normWplusJets_passed =
	new RooFormulaVar("normWplusJets_passed", "normWplusJets_passed",
			  "fr*frFactorWplusJets*normWplusJets", RooArgSet(*normWplusJets, *fr, *frFactorWplusJets));
      RooFormulaVar* normWplusJets_failed = 
	new RooFormulaVar("normWplusJets_failed", "normWplusJets_failed",
			  "(1.0 - fr*frFactorWplusJets)*normWplusJets", RooArgSet(*normWplusJets, *fr, *frFactorWplusJets));

      TH1* histogramZmumu_passed = histogramsZmumu[std::string(key).append("_").append("passed")];
      RooDataHist* tmpHistZmumu_passed = 
	new RooDataHist("tmpHistZmumu_passed", "tmpHistZmumu_passed", *fitVar, histogramZmumu_passed);
      RooHistPdf* templateZmumu_passed = 
	new RooHistPdf("templateZmumu_passed", "templateZmumu_passed", *fitVar, *tmpHistZmumu_passed);
      TH1* histogramZmumu_failed = histogramsZmumu[std::string(key).append("_").append("failed")];
      RooDataHist* tmpHistZmumu_failed = 
	new RooDataHist("tmpHistZmumu_failed", "tmpHistZmumu_failed", *fitVar, histogramZmumu_failed);
      RooHistPdf* templateZmumu_failed = 
	new RooHistPdf("templateZmumu_failed", "templateZmumu_failed", *fitVar, *tmpHistZmumu_failed);
      RooRealVar* normZmumu = new RooRealVar("normZmumu", "normZmumu", 0.10*numEventsData, 0., numEventsData);
      RooRealVar* frFactorZmumu = new RooRealVar("frFactorZmumu", "frFactorZmumu", 1., 0., 10.);
      RooFormulaVar* normZmumu_passed = 
	new RooFormulaVar("normZmumu_passed", "normZmumu_passed", 
			  "fr*frFactorZmumu*normZmumu", RooArgSet(*normZmumu, *fr, *frFactorZmumu));
      RooFormulaVar* normZmumu_failed = 
	new RooFormulaVar("normZmumu_failed", "normZmumu_failed", 
			  "(1.0 - fr*frFactorZmumu)*normZmumu", RooArgSet(*normZmumu, *fr, *frFactorZmumu));

      TH1* histogramZtautau_passed = histogramsZtautau[std::string(key).append("_").append("passed")];
      RooDataHist* tmpHistZtautau_passed = 
	new RooDataHist("tmpHistZtautau_passed", "tmpHistZtautau_passed", *fitVar, histogramZtautau_passed);
      RooHistPdf* templateZtautau_passed = 
	new RooHistPdf("templateZtautau_passed", "templateZtautau_passed", *fitVar, *tmpHistZtautau_passed);
      TH1* histogramZtautau_failed = histogramsZtautau[std::string(key).append("_").append("failed")];
      RooDataHist* tmpHistZtautau_failed = 
	new RooDataHist("tmpHistZtautau_failed", "tmpHistZtautau_failed", *fitVar, histogramZtautau_failed);
      RooHistPdf* templateZtautau_failed = 
	new RooHistPdf("templateZtautau_failed", "templateZtautau_failed", *fitVar, *tmpHistZtautau_failed);
      RooRealVar* normZtautau = new RooRealVar("normZtautau", "normZtautau", 0.25*numEventsData, 0., numEventsData);
      RooRealVar* effZtautau = new RooRealVar("effZtautau", "effZtautau", 0.50, 0., 1.);
      RooProduct* normZtautau_passed = 
	new RooProduct("normZtautau_passed", "normZtautau_passed", RooArgSet(*normZtautau, *effZtautau));
      RooFormulaVar* normZtautau_failed = 
	new RooFormulaVar("normZtautau_failed", "normZtautau_failed", 
			  "(1.0 - effZtautau)*normZtautau", RooArgSet(*normZtautau, *effZtautau));

      TObjArray templates_passed;
      templates_passed.Add(templateQCD_passed);
      templates_passed.Add(templateWplusJets_passed);
      templates_passed.Add(templateZmumu_passed);
      templates_passed.Add(templateZtautau_passed);
	
      TObjArray parameters_passed;
      parameters_passed.Add(normQCD_passed);
      parameters_passed.Add(normWplusJets_passed);
      parameters_passed.Add(normZmumu_passed);
      parameters_passed.Add(normZtautau_passed);
      
      RooAddPdf* pdfSMsumMC_passed 
	= new RooAddPdf("pdfSMsumMC_passed", "pdfSMsumMC_passed", RooArgList(templates_passed), RooArgList(parameters_passed));

      TObjArray templates_failed;
      templates_failed.Add(templateQCD_failed);
      templates_failed.Add(templateWplusJets_failed);
      templates_failed.Add(templateZmumu_failed);
      templates_failed.Add(templateZtautau_failed);

      TObjArray parameters_failed;
      parameters_failed.Add(normQCD_failed);
      parameters_failed.Add(normWplusJets_failed);
      parameters_failed.Add(normZmumu_failed);
      parameters_failed.Add(normZtautau_failed);

      RooAddPdf* pdfSMsumMC_failed 
	= new RooAddPdf("pdfSMsumMC_failed", "pdfSMsumMC_failed", RooArgList(templates_failed), RooArgList(parameters_failed));

      TObjArray templates_all;
      templates_all.Add(templateQCD_passed);
      templates_all.Add(templateQCD_failed);
      templates_all.Add(templateWplusJets_passed);
      templates_all.Add(templateWplusJets_failed);
      templates_all.Add(templateZmumu_passed);
      templates_all.Add(templateZmumu_failed);
      templates_all.Add(templateZtautau_passed);
      templates_all.Add(templateZtautau_failed);

      TObjArray parameters_all;
      parameters_all.Add(normQCD_passed);
      parameters_all.Add(normQCD_failed);
      parameters_all.Add(normWplusJets_passed);
      parameters_all.Add(normWplusJets_failed);
      parameters_all.Add(normZmumu_passed);
      parameters_all.Add(normZmumu_failed);
      parameters_all.Add(normZtautau_passed);
      parameters_all.Add(normZtautau_failed);

      RooAddPdf* pdfSMsumMC_all 
	= new RooAddPdf("pdfSMsumMC_all", "pdfSMsumMC_all", RooArgList(templates_all), RooArgList(parameters_all));
      
      RooCategory* fitCategories = new RooCategory("categories", "categories");
      fitCategories->defineType("passed");
      fitCategories->defineType("failed");
      fitCategories->defineType("all");

      RooSimultaneous* pdfSimultaneousFit = new RooSimultaneous("pdfSimultaneousFit", "pdfSimultaneousFit", *fitCategories);
      pdfSimultaneousFit->addPdf(*pdfSMsumMC_passed, "passed");
      pdfSimultaneousFit->addPdf(*pdfSMsumMC_failed, "failed");
      pdfSimultaneousFit->addPdf(*pdfSMsumMC_all, "all");
      
      std::map<std::string, TH1*> histogramDataMap;
      histogramDataMap["passed"] = histogramsData[std::string(key).append("_").append("passed")];
      histogramDataMap["failed"] = histogramsData[std::string(key).append("_").append("failed")];
      histogramDataMap["all"] = histogramsData[std::string(key).append("_").append("all")];
      
      RooDataHist* data = new RooDataHist("data", "data", *fitVar, *fitCategories, histogramDataMap);

      RooLinkedList fitOptions;
      fitOptions.Add(new RooCmdArg(RooFit::Extended()));
      fitOptions.Add(new RooCmdArg(RooFit::PrintLevel(1)));
      fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(true)));
      fitOptions.Add(new RooCmdArg(RooFit::Warnings(true)));

      RooConstVar* frMean = new RooConstVar("frMean", "frMean", 0.05);
      RooConstVar* frSpread = new RooConstVar("frSpread", "frSpread", 0.025);

      RooGaussian* fr_constraint = new RooGaussian("fr_constraint", "fr_constraint", *fr, *frMean, *frSpread);

      RooConstVar* frFactorMean = new RooConstVar("frFactorMean", "frFactorMean", 1.);
      RooConstVar* frFactorSpread = new RooConstVar("frFactorSpread", "frFactorSpread", 0.5);

      RooGaussian* frFactorQCD_constraint = 
	new RooGaussian("frFactorQCD_constraint", "frFactorQCD_constraint", 
			*frFactorQCD, *frFactorMean, *frFactorSpread);
      RooGaussian* frFactorWplusJets_constraint = 
	new RooGaussian("frFactorWplusJets_constraint", "frFactorWplusJets_constraint", 
			*frFactorWplusJets, *frFactorMean, *frFactorSpread);
      RooGaussian* frFactorZmumu_constraint = 
	new RooGaussian("frFactorZmumu_constraint", "frFactorZmumu_constraint", 
			*frFactorZmumu, *frFactorMean, *frFactorSpread);

      TObjArray constraints;
      constraints.Add(fr_constraint);
      constraints.Add(frFactorQCD_constraint);
      constraints.Add(frFactorWplusJets_constraint);
      constraints.Add(frFactorZmumu_constraint);

      fitOptions.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(constraints))));

      pdfSimultaneousFit->fitTo(*data, fitOptions);
      
      std::cout << "Results of fitting variable = " << (*variable)
		<< " for Tau id. = " << (*tauId) << std::endl;
      std::cout << " Ztautau:" << std::endl;
      std::cout << "  normalization = " << normZtautau->getVal() << " +/- " << normZtautau->getError() << std::endl;
      std::cout << "  efficiency = " << effZtautau->getVal() << " +/- " << effZtautau->getError() << std::endl;
      std::cout << " Zmumu:" << std::endl; 
      std::cout << "  normalization = " << normZmumu->getVal() << " +/- " << normZmumu->getError() << std::endl; 
      double frValueZmumu = fr->getVal()*frFactorZmumu->getVal();
      double frErrorZmumu = frValueZmumu*TMath::Sqrt(TMath::Power(fr->getError()/fr->getVal(), 2) 
						    + TMath::Power(frFactorZmumu->getError()/frFactorZmumu->getVal(), 2));
      std::cout << "  fake-rate = " << frValueZmumu << " +/- " << frErrorZmumu << std::endl;
      std::cout << " WplusJets:" << std::endl;
      std::cout << "  normalization = " << normWplusJets->getVal() << " +/- " << normWplusJets->getError() << std::endl; 
      double frValueWplusJets = fr->getVal()*frFactorWplusJets->getVal();
      double frErrorWplusJets = frValueWplusJets*TMath::Sqrt(TMath::Power(fr->getError()/fr->getVal(), 2) 
						            + TMath::Power(frFactorWplusJets->getError()/frFactorWplusJets->getVal(), 2));
      std::cout << "  fake-rate = " << frValueWplusJets << " +/- " << frErrorWplusJets << std::endl;
      std::cout << " QCD:" << std::endl;
      std::cout << "  normalization = " << normQCD->getVal() << " +/- " << normQCD->getError() << std::endl; 
      double frValueQCD = fr->getVal()*frFactorQCD->getVal();
      double frErrorQCD = frValueQCD*TMath::Sqrt(TMath::Power(fr->getError()/fr->getVal(), 2) 
						+ TMath::Power(frFactorQCD->getError()/frFactorQCD->getVal(), 2));
      std::cout << "  fake-rate = " << frValueQCD << " +/- " << frErrorQCD << std::endl;

      effValues[key] = effZtautau->getVal();
      effErrors[key] = effZtautau->getError();
      
      drawHistograms(histogramQCD_passed, normQCD->getVal()*fr->getVal()*frFactorQCD->getVal(),
		     histogramWplusJets_passed, normWplusJets->getVal()*fr->getVal()*frFactorWplusJets->getVal(),
		     histogramZmumu_passed, normZmumu->getVal()*fr->getVal()*frFactorZmumu->getVal(),
		     histogramZtautau_passed, normZtautau->getVal()*effZtautau->getVal(),
		     histogramDataMap["passed"],
		     std::string("fitTauIdEff_").append(key).append("_passed.pdf"));
      drawHistograms(histogramQCD_failed, normQCD->getVal()*(1.0 - fr->getVal()*frFactorQCD->getVal()),
		     histogramWplusJets_failed, normWplusJets->getVal()*(1.0 - fr->getVal()*frFactorWplusJets->getVal()),
		     histogramZmumu_failed, normZmumu->getVal()*(1.0 - fr->getVal()*frFactorZmumu->getVal()),
		     histogramZtautau_failed, normZtautau->getVal()*(1.0 - effZtautau->getVal()),
		     histogramDataMap["failed"],
		     std::string("fitTauIdEff_").append(key).append("_failed.pdf"));

      drawHistograms(histogramQCD_passed, histogramQCD_failed, *tauId,
		     std::string("fitTauIdEff_").append(key).append("_QCD.pdf"));
      drawHistograms(histogramWplusJets_passed, histogramWplusJets_failed, *tauId,
		     std::string("fitTauIdEff_").append(key).append("_WplusJets.pdf"));
      drawHistograms(histogramZmumu_passed, histogramZmumu_failed, *tauId,
		     std::string("fitTauIdEff_").append(key).append("_Zmumu.pdf"));
      drawHistograms(histogramZtautau_passed, histogramZtautau_failed, *tauId,
		     std::string("fitTauIdEff_").append(key).append("_Ztautau.pdf"));

      delete fr;
      
      delete tmpHistQCD_passed;
      delete templateQCD_passed;
      delete tmpHistQCD_failed;
      delete templateQCD_failed;
      delete normQCD;
      delete frFactorQCD;
      delete normQCD_passed;
      delete normQCD_failed;

      delete tmpHistWplusJets_passed;
      delete templateWplusJets_passed;
      delete tmpHistWplusJets_failed;
      delete templateWplusJets_failed;
      delete normWplusJets;
      delete frFactorWplusJets;
      delete normWplusJets_passed;
      delete normWplusJets_failed;

      delete tmpHistZmumu_passed;
      delete templateZmumu_passed;
      delete tmpHistZmumu_failed;
      delete templateZmumu_failed;
      delete normZmumu;
      delete frFactorZmumu;
      delete normZmumu_passed;
      delete normZmumu_failed;

      delete tmpHistZtautau_passed;
      delete templateZtautau_passed;
      delete tmpHistZtautau_failed;
      delete templateZtautau_failed;
      delete normZtautau;
      delete effZtautau;
      delete normZtautau_passed;
      delete normZtautau_failed;

      delete pdfSMsumMC_passed;
      delete pdfSMsumMC_failed;
      delete pdfSMsumMC_all;
      
      delete fitCategories;

      delete pdfSimultaneousFit;
      
      delete data;
      
      delete fitVar;

      delete frMean;
      delete frSpread;

      delete fr_constraint;

      delete frFactorMean;
      delete frFactorSpread;

      delete frFactorQCD_constraint;
      delete frFactorWplusJets_constraint;
      delete frFactorZmumu_constraint;
    }
  }
}

void fitTauIdEff()
{
  gROOT->SetBatch(true);

  const std::string inputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas/local/";

  const std::string treeSelection = "";

  std::vector<std::string> tauIds;
  tauIds.push_back(std::string("discrTaNCloose"));
  tauIds.push_back(std::string("discrTaNCmedium"));
  tauIds.push_back(std::string("discrTaNCtight"));
  tauIds.push_back(std::string("discrHPSloose"));
  tauIds.push_back(std::string("discrHPSmedium"));
  tauIds.push_back(std::string("discrHPStight"));
  
  std::vector<std::string> tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");
  tauIdValues.push_back("all");

  std::vector<std::string> fitVariables;
  //fitVariables.push_back("Ht");
  //fitVariables.push_back("SVfitMass1");
  //fitVariables.push_back("SVfitMass2");
  fitVariables.push_back("visMass");
  //fitVariables.push_back("muonPt");
  //fitVariables.push_back("muonEta");
  //fitVariables.push_back("tauPt");
  //fitVariables.push_back("tauEta");
  //fitVariables.push_back("numChargedParticles");
  //fitVariables.push_back("numParticles");
  //fitVariables.push_back("jetWidth"); CV: signal/background normalizations --> tau id. efficiencies obtained by using jetWidth variable
  //                                        are **very** different from values obtained by using all other variables
  //                                       --> there seems to be a problem in modeling jetWidth variable
  //                                       --> do not use jetWidth variable for now

  TChain* chainData = new TChain("Events");
  chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010A_Sep17ReReco_Run26.root").data());
  chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010B_Prompt_Run26_*.root").data());
  std::map<std::string, TH1*> histogramsData = 
    makeHistograms("Data", 1.0, chainData, treeSelection, tauIds, tauIdValues, fitVariables);

  TChain* chainZtautau = new TChain("Events");
  chainZtautau->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_ZtautauPU156bx_Run26_*.root").data());
  std::map<std::string, TH1*> histogramsZtautau = 
    makeHistograms("Ztautau", 0.030*1.02*2.63, chainZtautau, treeSelection, tauIds, tauIdValues, fitVariables);

  TChain* chainZmumu = new TChain("Events");
  chainZmumu->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Zmumu_Run26_*.root").data());
  std::map<std::string, TH1*> histogramsZmumu = 
    makeHistograms("Zmumu", 0.026, chainZmumu, treeSelection, tauIds, tauIdValues, fitVariables);

  TChain* chainWplusJets = new TChain("Events");
  chainWplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_WplusJets_Run26_*.root").data());
  double fudgeFactorWplusJets = 1.10;
  std::map<std::string, TH1*> histogramsWplusJets = 
    makeHistograms("WplusJets", 0.12*1.04*fudgeFactorWplusJets, chainWplusJets, treeSelection, tauIds, tauIdValues, fitVariables);

  TChain* chainQCD = new TChain("Events");
  chainQCD->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_Run26_*.root").data());
  double fudgeFactorQCD = 1.35;
  std::map<std::string, TH1*> histogramsQCD = 
    makeHistograms("QCD", 0.10*2.04*1.01*fudgeFactorQCD, chainQCD, treeSelection, tauIds, tauIdValues, fitVariables);
  
  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {

	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);

	drawHistograms(histogramsQCD[key], -1.,
		       histogramsWplusJets[key], -1.,
		       histogramsZmumu[key], -1.,
		       histogramsZtautau[key], -1.,
		       histogramsData[key],
		       std::string("plotsTauIdEff_").append(key).append(".pdf"));
      }
    }
  }

  std::map<std::string, double> effValues;
  std::map<std::string, double> effErrors;
  
  //fitUsingTFractionFitter(histogramsData, histogramsZtautau, histogramsZmumu, histogramsWplusJets, histogramsQCD,
  //		  	    tauIds, tauIdValues, fitVariables, effValues, effErrors);
  fitUsingRooFit(histogramsData, histogramsZtautau, histogramsZmumu, histogramsWplusJets, histogramsQCD,
		 tauIds, tauIdValues, fitVariables, effValues, effErrors);

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    
    std::cout << "Efficiency of Tau id. = " << (*tauId) << ":" << std::endl;
    
    for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	  variable != fitVariables.end(); ++variable ) {

      std::string key = std::string(*variable);
      key.append("_").append(*tauId);

      std::cout << " fitVariable = " << (*variable) << ":" 
		<< " result = " << effValues[key]*100. << " +/- " << effErrors[key]*100. << "%" << std::endl;
      double numeratorZtautau = histogramsZtautau[std::string(key).append("_passed")]->Integral();
      double denominatorZtautau = histogramsZtautau[std::string(key).append("_all")]->Integral();
      std::cout << "(Monte Carlo prediction = " << (numeratorZtautau/denominatorZtautau) << ")" << std::endl;
    }
  }
}
