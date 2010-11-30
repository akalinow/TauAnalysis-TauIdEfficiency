
#include "RooAddPdf.h"
#include "RooCmdArg.h"
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"

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
	  numBins = 25;
	  min = 0.;
	  max = 250.;
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

TH1* normalize(const TH1* histogram)
{
  TH1* retVal = (TH1*)histogram->Clone();

  if ( !retVal->GetSumw2N() ) retVal->Sumw2();

  if ( retVal->Integral() != 0. ) retVal->Scale(1./retVal->Integral());
    
  return retVal;
}

void fitUsingTFractionFitter(std::map<std::string, TH1*>& histogramsData,
			     std::map<std::string, TH1*>& histogramsZtautau,
			     std::map<std::string, TH1*>& histogramsZmumu,
			     std::map<std::string, TH1*>& histogramsWplusJets,
			     std::map<std::string, TH1*>& histogramsQCD,
			     const std::vector<std::string>& tauIds, const std::vector<std::string>& tauIdValues,
			     const std::vector<std::string>& fitVariables,
			     std::map<std::string, double>& fitValues,
			     std::map<std::string, double>& fitErrors)
{
  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
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

	  std::cout << "Fit Parameter:" << std::endl;
	  std::cout << " Ztautau: normalization = " << fitValueZtautau << " +/- " << fitErrorZtautau << std::endl;
	  std::cout << " Zmumu: normalization = " << fitValueZmumu << " +/- " << fitErrorZmumu << std::endl; 
	  std::cout << " WplusJets: normalization = " << fitValueWplusJets << " +/- " << fitErrorWplusJets << std::endl;
	  std::cout << " QCD: normalization = " << fitValueQCD << " +/- " << fitErrorQCD << std::endl; 

	  std::string key = std::string(*variable);
	  key.append("_").append(*tauId).append("_").append(*tauIdValue);

	  fitValues[std::string("Ztautau").append("_").append(key)] = fitValueZtautau;
	  fitErrors[std::string("Ztautau").append("_").append(key)] = fitErrorZtautau;
	  fitValues[std::string("Zmumu").append("_").append(key)] = fitValueZmumu;
	  fitErrors[std::string("Zmumu").append("_").append(key)] = fitErrorZmumu;
	  fitValues[std::string("WplusJets").append("_").append(key)] = fitValueWplusJets;
	  fitErrors[std::string("WplusJets").append("_").append(key)] = fitErrorWplusJets;
	  fitValues[std::string("QCD").append("_").append(key)] = fitValueQCD;
	  fitErrors[std::string("QCD").append("_").append(key)] = fitErrorQCD;
	}
      }
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
		    std::map<std::string, double>& fitValues,
		    std::map<std::string, double>& fitErrors)
{
  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {

	std::cout << "performing Fit of variable = " << (*variable)
		  << " for Tau id. = " << (*tauId) << ", value = "<< (*tauIdValue) << std::endl;
	
	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);
     
	double fitMin = histogramsZtautau[key]->GetXaxis()->GetXmin();
	double fitMax = histogramsZtautau[key]->GetXaxis()->GetXmax();

	RooRealVar* fitVar = new RooRealVar("fitVar", "fitVar", fitMin, fitMax);

        double numEventsData = histogramsData[key]->Integral();
	std::cout << " numEventsData = " << numEventsData << std::endl;

	TH1* histogramQCD = histogramsQCD[key];
	RooDataHist* tmpHistQCD = new RooDataHist("tmpHistQCD", "tmpHistQCD", *fitVar, histogramQCD);
        RooHistPdf* templateQCD = new RooHistPdf("templateQCD", "templateQCD", *fitVar, *tmpHistQCD);
        RooRealVar* normQCD = new RooRealVar("normQCD", "normQCD", 0.25*numEventsData, 0., numEventsData);

	TH1* histogramWplusJets = histogramsWplusJets[key];
        RooDataHist* tmpHistWplusJets = new RooDataHist("tmpHistWplusJets", "tmpHistWplusJets", *fitVar, histogramWplusJets);
        RooHistPdf* templateWplusJets = new RooHistPdf("templateWplusJets", "templateWplusJets", *fitVar, *tmpHistWplusJets);
        RooRealVar* normWplusJets = new RooRealVar("normWplusJets", "normWplusJets", 0.25*numEventsData, 0., numEventsData);

	TH1* histogramZmumu = histogramsZmumu[key];
        RooDataHist* tmpHistZmumu = new RooDataHist("tmpHistZmumu", "tmpHistZmumu", *fitVar, histogramZmumu);
        RooHistPdf* templateZmumu = new RooHistPdf("templateZmumu", "templateZmumu", *fitVar, *tmpHistZmumu);
        RooRealVar* normZmumu = new RooRealVar("normZmumu", "normZmumu", 0.25*numEventsData, 0., numEventsData);

	TH1* histogramZtautau = histogramsZtautau[key];
        RooDataHist* tmpHistZtautau = new RooDataHist("tmpHistZtautau", "tmpHistZtautau", *fitVar, histogramZtautau);
        RooHistPdf* templateZtautau = new RooHistPdf("templateZtautau", "templateZtautau", *fitVar, *tmpHistZtautau);
        RooRealVar* normZtautau = new RooRealVar("normZtautau", "normZtautau", 0.25*numEventsData, 0., numEventsData);

        TObjArray templates;
	templates.Add(templateQCD);
	templates.Add(templateWplusJets);
	templates.Add(templateZmumu);
	templates.Add(templateZtautau);
	
	TObjArray normalizations;
	normalizations.Add(normQCD);
	normalizations.Add(normWplusJets);
	normalizations.Add(normZmumu);
	normalizations.Add(normZtautau);

	RooAddPdf* pdfSMsumMC = new RooAddPdf("pdfSMsumMC", "pdfSMsumMC", RooArgList(templates), RooArgList(normalizations));

	TH1* histogramData = histogramsData[key];
	RooDataHist* data = new RooDataHist("data", "data", *fitVar, histogramData);

	RooLinkedList fitOptions;
	fitOptions.Add(new RooCmdArg(RooFit::Extended()));
	fitOptions.Add(new RooCmdArg(RooFit::PrintLevel(1)));
	fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(true)));
	fitOptions.Add(new RooCmdArg(RooFit::Warnings(true)));

	pdfSMsumMC->fitTo(*data, fitOptions);

	double fitValueQCD = normQCD->getVal();
	double fitErrorQCD = normQCD->getError();
        double fitValueWplusJets = normWplusJets->getVal();
	double fitErrorWplusJets = normWplusJets->getError();
	double fitValueZmumu = normZmumu->getVal();
	double fitErrorZmumu = normZmumu->getError();
        double fitValueZtautau = normZtautau->getVal();
	double fitErrorZtautau = normZtautau->getError();

	std::string key2 = std::string(*variable);
	key2.append("_").append(*tauId).append("_").append(*tauIdValue);
	
	fitValues[std::string("Ztautau").append("_").append(key2)] = fitValueZtautau;
	fitErrors[std::string("Ztautau").append("_").append(key2)] = fitErrorZtautau;
	fitValues[std::string("Zmumu").append("_").append(key2)] = fitValueZmumu;
	fitErrors[std::string("Zmumu").append("_").append(key2)] = fitErrorZmumu;
	fitValues[std::string("WplusJets").append("_").append(key2)] = fitValueWplusJets;
	fitErrors[std::string("WplusJets").append("_").append(key2)] = fitErrorWplusJets;
	fitValues[std::string("QCD").append("_").append(key2)] = fitValueQCD;
	fitErrors[std::string("QCD").append("_").append(key2)] = fitErrorQCD;
	
        delete tmpHistQCD;
	delete templateQCD;
        delete tmpHistWplusJets;
	delete templateWplusJets;
        delete tmpHistZmumu;
	delete templateZmumu;
        delete tmpHistZtautau;
	delete templateZtautau; 

        delete pdfSMsumMC;
        delete data;

        delete fitVar;
      }
    }
  }
}

void fitTauIdEff()
{
  gROOT->SetBatch(true);

  const std::string inputFilePath = "/data1/veelken/CMSSW_3_8_x/ntuples/TauIdEffMeas/local/";

  const std::string treeSelection = "";

  std::vector<std::string> tauIds;
  tauIds.push_back(std::string("discrTaNCmedium"));
  
  std::vector<std::string> tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");
  tauIdValues.push_back("all");

  std::vector<std::string> fitVariables;
  fitVariables.push_back("Ht");
  fitVariables.push_back("SVfitMass1");
  fitVariables.push_back("SVfitMass2");
  fitVariables.push_back("visMass");
  //fitVariables.push_back("muonPt");
  //fitVariables.push_back("muonEta");
  //fitVariables.push_back("tauPt");
  //fitVariables.push_back("tauEta");
  fitVariables.push_back("numChargedParticles");
  fitVariables.push_back("numParticles");
  fitVariables.push_back("jetWidth");

  TChain* chainData = new TChain("Events");
  chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010A_Sep17ReReco_Run26.root").data());
  chainData->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_Mu_Run2010B_Prompt_Run26_*.root").data());
  std::map<std::string, TH1*> histogramsData = 
    makeHistograms("Data", 1.0, chainData, treeSelection, tauIds, tauIdValues, fitVariables);

  TChain* chainZtautau = new TChain("Events");
  chainZtautau->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_ZtautauPU156bx_Run26_*.root").data());
  std::map<std::string, TH1*> histogramsZtautau = 
    makeHistograms("Ztautau", 0.030*1.02, chainZtautau, treeSelection, tauIds, tauIdValues, fitVariables);

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
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 640);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {

	canvas->Clear();

	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);

	THStack smSum("smSum", "smSum");

	TH1* histogramQCD = histogramsQCD[key];
	histogramQCD->SetFillColor(797);
	smSum.Add(histogramQCD);

	TH1* histogramWplusJets = histogramsWplusJets[key];
	histogramWplusJets->SetFillColor(856);
	smSum.Add(histogramWplusJets);

	TH1* histogramZmumu = histogramsZmumu[key];
	histogramZmumu->SetFillColor(596);
	smSum.Add(histogramZmumu);

	TH1* histogramZtautau = histogramsZtautau[key];
	histogramZtautau->SetFillColor(628);
	smSum.Add(histogramZtautau);

	TH1* histogramData = histogramsData[key];
	histogramData->SetLineColor(1);
	histogramData->SetMarkerColor(1);
	histogramData->SetMarkerStyle(20);

	smSum.SetTitle(histogramZtautau->GetTitle());
	smSum.SetMaximum(1.4*TMath::Max(smSum.GetMaximum(), histogramData->GetMaximum()));
	
	smSum.Draw("hist");
	histogramData->SetStats(false);
	histogramData->Draw("ep1same");

	TLegend legend(0.64, 0.69, 0.89, 0.89, "", "brNDC"); 
	legend.SetBorderSize(0);
	legend.SetFillColor(0);

	legend.AddEntry(histogramZtautau, "Z #rightarrow #tau^{+} #tau^{-}", "f");
	legend.AddEntry(histogramZmumu, "Z #rightarrow #mu^{+} #mu^{-}", "f");
	legend.AddEntry(histogramWplusJets, "W + jets", "f");
	legend.AddEntry(histogramQCD, "QCD", "f");
	legend.AddEntry(histogramData, "Data", "p");
	legend.Draw();

	canvas->Update();
	std::string outputFileName = std::string("./plots/").append("plotsTauIdEff_").append(key).append(".png");
	canvas->Print(outputFileName.data());
      }
    }
  }

  std::map<std::string, double> fitValues;
  std::map<std::string, double> fitErrors;
  
  //fitUsingTFractionFitter(histogramsData, histogramsZtautau, histogramsZmumu, histogramsWplusJets, histogramsQCD,
  //		  	    tauIds, tauIdValues, fitVariables, fitValues, fitErrors);
  fitUsingRooFit(histogramsData, histogramsZtautau, histogramsZmumu, histogramsWplusJets, histogramsQCD,
		 tauIds, tauIdValues, fitVariables, fitValues, fitErrors);

  for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	variable != fitVariables.end(); ++variable ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator tauIdValue = tauIdValues.begin();
	    tauIdValue != tauIdValues.end(); ++tauIdValue ) {
	
	std::string key = std::string(*variable);
	key.append("_").append(*tauId).append("_").append(*tauIdValue);

	double fitValueZtautau = fitValues[std::string("Ztautau").append("_").append(key)];
	double fitErrorZtautau = fitErrors[std::string("Ztautau").append("_").append(key)];
	double fitValueZmumu = fitValues[std::string("Zmumu").append("_").append(key)];
	double fitErrorZmumu = fitErrors[std::string("Zmumu").append("_").append(key)];
	double fitValueWplusJets = fitValues[std::string("WplusJets").append("_").append(key)];
	double fitErrorWplusJets = fitErrors[std::string("WplusJets").append("_").append(key)];
	double fitValueQCD = fitValues[std::string("QCD").append("_").append(key)];
	double fitErrorQCD = fitErrors[std::string("QCD").append("_").append(key)];

	std::cout << "Results of fitting variable = " << (*variable)
		  << " for Tau id. = " << (*tauId) << ", value = '"<< (*tauIdValue) << "':" << std::endl;
	std::cout << " Ztautau: normalization = " << fitValueZtautau << " +/- " << fitErrorZtautau << std::endl;
	std::cout << " Zmumu: normalization = " << fitValueZmumu << " +/- " << fitErrorZmumu << std::endl; 
	std::cout << " WplusJets: normalization = " << fitValueWplusJets << " +/- " << fitErrorWplusJets << std::endl;
	std::cout << " QCD: normalization = " << fitValueQCD << " +/- " << fitErrorQCD << std::endl; 
      }
    }
  }

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    
    std::cout << "Efficiency of Tau id. = " << (*tauId) << ":" << std::endl;
    
    for ( std::vector<std::string>::const_iterator variable = fitVariables.begin();
	  variable != fitVariables.end(); ++variable ) {

      std::string key = std::string("Ztautau").append("_").append(*variable);
      key.append("_").append(*tauId);
      
      double passedValue = fitValues[std::string(key).append("_").append("passed")];            
      double passedError = fitErrors[std::string(key).append("_").append("passed")]; 
      double failedValue = fitValues[std::string(key).append("_").append("failed")];
      double failedError = fitErrors[std::string(key).append("_").append("failed")];

      //std::cout << "passed = " << passedValue << " +/- " << passedError << std::endl;
      //std::cout << "failed = " << failedValue << " +/- " << failedError << std::endl;

      double numerator = passedValue;
      double denominator = passedValue + failedValue;

      double effValue = numerator/denominator;
      double effError2 = TMath::Power((failedValue/denominator)*(passedError/denominator), 2)
	                + TMath::Power((passedValue/denominator)*(failedError/denominator), 2);

      std::cout << " fitVariable = " << (*variable) << ":" 
		<< " result = " << effValue*100. << " +/- " << TMath::Sqrt(effError2)*100. << "%" << std::endl;
    }
  }
}
