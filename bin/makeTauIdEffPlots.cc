
#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"

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
  branchName.append("_").append(branchName_suffix);
  branchName.append(".obj");

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
  const std::string& sysShift, const std::string& branchName_suffix)
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
      branchNameTau   = "selectedPatPFTausShrinkingConePFRelIsoCumulative"; // "old" TaNC algorithm
      branchNameDiTau = "selectedMuPFTauShrinkingConePairsForTauIdEffCumulative";  
    } else if ( (*tauId) == "tauDiscrTaNCloose"            || 
		(*tauId) == "tauDiscrTaNCmedium"           || 
		(*tauId) == "tauDiscrTaNCtight"            ) { // "new" TaNC implemented in HPS+TaNC combined algorithm
      branchNameTau   = "selectedPatPFTausHPSpTaNCPFRelIsoCumulative";    
      branchNameDiTau = "selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative"; 
    } else if ( (*tauId) == "tauDiscrIsolationLoose"       ||
		(*tauId) == "tauDiscrIsolationMedium"      ||
		(*tauId) == "tauDiscrIsolationTight"       ) { // "old" HPS algorithm
      branchNameTau   = "selectedPatPFTausHPSPFRelIsoCumulative";           
      branchNameDiTau = "selectedMuPFTauHPSpairsForTauIdEffCumulative";
    } else if ( (*tauId) == "tauDiscrHPSloose"             ||
		(*tauId) == "tauDiscrHPSmedium"            || 
		(*tauId) == "tauDiscrHPStight"             ) { // "new" HPS implemented in HPS+TaNC combined algorithm
      branchNameTau   = "selectedPatPFTausHPSpTaNCPFRelIsoCumulative"; 
      branchNameDiTau = "selectedMuPFTauHPSpTaNCpairsForTauIdEffCumulative";
    } else {
      std::cout << "Error in <makeBranchNameDict>: invalid tauId = " << (*tauId) << " --> aborting !!";
      assert(0);
    }
    
    if (        sysShift == "CENTRAL_VALUE"              ) {
      // nothing to be done yet.
    } else if ( sysShift == "SysTauJetEnUp"              ||
		sysShift == "SysTauJetEnDown"            ) {
      branchNameTau   = replace(branchNameTau,   "Cumulative", std::string(sysShift).append("Cumulative"));
      branchNameDiTau = replace(branchNameDiTau, "Cumulative", std::string(sysShift).append("Cumulative"));
    } else if ( sysShift == "SysJetEnUp"                 ||
		sysShift == "SysJetEnDown"               ) {
      //branchNameTau   = replace(branchNameTau,   "Cumulative", std::string(sysShift).append("Cumulative"));
      branchNameDiTau = replace(branchNameDiTau, "Cumulative", std::string(sysShift).append("Cumulative"));
    } else assert(0);
    
    branchNameDictEntry branchNames;

    branchNames["event"] = getBranchName("double", "", "event", branchName_suffix);
    branchNames["ls"] = getBranchName("double", "", "ls", branchName_suffix);
    branchNames["run"] = getBranchName("double", "", "run", branchName_suffix);
    
    branchNames["genPUreweight"] = getBranchName("double", "addPileupInfo", "vtxMultReweight", branchName_suffix);

    branchNames["muonPt"] = getBranchName("double", branchNameMuon, "pt", branchName_suffix);
    branchNames["muonEta"] = getBranchName("double", branchNameMuon, "eta", branchName_suffix);
    branchNames["muonLooseIsoPtSum04"] = getBranchName("double", branchNameMuon, "ptSumLooseIsolation04", branchName_suffix);
    branchNames["muonLooseIsoPtSum06"] = getBranchName("double", branchNameMuon, "ptSumLooseIsolation06", branchName_suffix);
    branchNames["numMuonsGlobal"] = getBranchName("double", "patMuonsGlobal", "multiplicity", branchName_suffix);
    branchNames["numMuonsStandAlone"] = getBranchName("double", "patMuonsStandAlone", "multiplicity", branchName_suffix);

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
    
    branchNames["diTauCharge"] = getBranchName("double", branchNameDiTau, "chargeTauLeadTrack", branchName_suffix);
    branchNames["diTauMt"] = getBranchName("double", branchNameDiTau, "Mt", branchName_suffix);
    branchNames["diTauPzeta"] = getBranchName("double", branchNameDiTau, "pZeta", branchName_suffix);
    branchNames["diTauPzetaVis"] = getBranchName("double", branchNameDiTau, "pZetaVis", branchName_suffix);
    branchNames["diTauHt"] = getBranchName("double", branchNameDiTau, "Ht", branchName_suffix);
    branchNames["diTauVisMass"] = getBranchName("double", branchNameDiTau, "visMass", branchName_suffix);
    branchNames["diTauVisMassFromJet"] = getBranchName("double", branchNameDiTau, "visMassFromJet", branchName_suffix);    
    branchNames["muVertexZ"] = getBranchName("double", branchNameDiTau, "muVertexZ", branchName_suffix);
    branchNames["tauVertexZ"] = getBranchName("double", branchNameDiTau, "tauVertexZ", branchName_suffix);
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

void makeHistograms(
  const std::string& process, 
  std::map<std::string, TH1*>& histograms,
  const std::string& region, double weight,
  TTree* tree, const std::string& treeSelection, 
  const std::string& tauId, const std::vector<std::string>& tauIdValues,
  const std::vector<std::string>& observables,
  std::map<std::string, std::string>& branchNames, 
  bool applyPUreweighting,
  const std::string& sysShift = "CENTRAL_VALUE",
  bool saveRunLumiSectionEventNumbers = false)
{
//--------------------------------------------------------------------------------
// Fill histograms with (ED)NTuple entries passing treeSelection
//--------------------------------------------------------------------------------

  //std::cout << "<makeHistograms>:" << std::endl;
  //std::cout << " process = " << process << std::endl;
  //std::cout << " region = " << region << std::endl;
  //std::cout << " treeSelection = " << treeSelection << std::endl;
  //std::cout << " applyPUreweighting = " << applyPUreweighting << std::endl;

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
    std::string extTreeSelection;

    if ( applyPUreweighting ) extTreeSelection.append("(");

    if ( treeSelection != "" ) {
      extTreeSelection.append(treeSelection);
      extTreeSelection.append(" && ");
    }
    extTreeSelection.append(branchNames["muonPt"]).append(" > 20");
    extTreeSelection.append(" && ").append(branchNames["tauPt"]).append(" > 20");
    extTreeSelection.append(" && ").append(branchNames["tauJetPt"]).append(" > 20");
    extTreeSelection.append(" && abs(").append(branchNames["tauEta"]).append(") < 2.3");
    extTreeSelection.append(" && abs(").append(branchNames["tauJetEta"]).append(") < 2.3");
    extTreeSelection.append(" && ").append(branchNames["tauLooseIsoPtSum06"]).append(" < 2.5");
    extTreeSelection.append(" && ").append(branchNames["numMuonsStandAlone"]).append(" < 1.5");
    extTreeSelection.append(" && ").append(branchNames["diTauVisMassFromJet"]).append(" > 20");
    extTreeSelection.append(" && ").append(branchNames["diTauVisMassFromJet"]).append(" < 200");
    extTreeSelection.append(" && ").append(branchNames["diTauMt"]).append(" < 80");
    extTreeSelection.append(" && ").append("abs(" + branchNames["muVertexZ"] + "-" + branchNames["tauVertexZ"] + ")").append(" < 0.2");

//--- add tau id. passed/failed selection
    if      ( (*tauIdValue) == "passed" ) extTreeSelection.append(" && ").append(branchNames[tauId]).append(" > 0.5");
    else if ( (*tauIdValue) == "failed" ) extTreeSelection.append(" && ").append(branchNames[tauId]).append(" < 0.5");
    else if ( (*tauIdValue) == "all"    ) {}
    else assert(0);

    if ( applyPUreweighting ) extTreeSelection.append(")").append("*").append(branchNames["genPUreweight"]);

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
	numBins = 36; // CV: use 5 GeV binning, in order to rebin plots to 10 and 15 GeV bin-width
	min = 20.;
	max = 200.;
      } else if ( (*observable) == "diTauHt" ) {
	numBins = 36; // CV: use 5 GeV binning, in order to rebin plots to 10 and 15 GeV bin-width
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
      } else if ( (*observable) == "tauNumParticles" ) {
	numBins = 25;
	min =  -0.5;
	max = +24.5;
      } else if ( (*observable) == "tauJetWidth" ) {
	numBins = 20;
	min = 0.;
	max = 0.50;
      } else {
	std::cout << "Error in <makeHistograms>: undefined observable = " << (*observable) << " --> skipping !!" << std::endl;
	return;
      }
      
      TH1* histogram = new TH1F(histogramName.data(), histogramName.data(), numBins, min, max);

      std::string drawCommand = std::string(branchNames[*observable]).append(">>").append(histogramName); 
      tree->Draw(drawCommand.data(), extTreeSelection.data());

      if ( !histogram->GetSumw2N() ) histogram->Sumw2();
      histogram->Scale(weight);

      double integral = getIntegral(histogram, true, true);
      double fittedFraction = ( integral > 0. ) ? getIntegral(histogram, false, false)/integral : -1.; 
      std::cout << "histogram = " << histogramName << ":" 
		<< " entries = " << histogram->GetEntries() << ", integral = " << integral 
		<< " (fitted fraction = " << fittedFraction << ")" << std::endl;
      
      std::string key = getKey(*observable, tauId, *tauIdValue, sysShift);	
      if ( histogram != 0 ) histograms[key] = histogram;
      std::cout << "--> storing histogram: key = " << key << std::endl;

      if ( histogram->GetEntries() > 0 && saveRunLumiSectionEventNumbers && sysShift == "CENTRAL_VALUE" ) 
	writeRunLumiSectionEventNumberFile(process, tree, extTreeSelection, region, tauId, *tauIdValue, branchNames);
    }
  }

  std::cout << std::endl;
}

void makeDistributionsInRegion(
  const std::string& process, 
  std::map<std::string, TH1*>& histograms,
  const std::string& region, double weight,
  TTree* tree, 
  const std::string& tauId, const std::vector<std::string>& fitVariables,
  std::map<std::string, std::string>& branchNames, 
  bool applyPUreweighting,
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
    
  if        ( region           == "ABCD"            ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprDiTauCharge_OS_or_SS);
  } else if ( region.find("A") != std::string::npos ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprMuonIso_loose_not_tight).append(" && ").append(exprDiTauCharge_OS);
  } else if ( region.find("B") != std::string::npos ) {
    treeSelection.append(exprMuonIso_loose).append(" && ").append(exprMuonIso_loose_not_tight).append(" && ").append(exprDiTauCharge_SS);
  } else if ( region.find("C") != std::string::npos ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_OS);
  } else if ( region.find("D") != std::string::npos ) {
    treeSelection.append(exprMuonIso_tight).append(" && ").append(exprDiTauCharge_SS);
  } else {
    std::cout << "Error in <makeDistribution>: undefined region = " << region << " --> skipping !!" << std::endl;
    return;
  }

  if      ( region.find("1") != std::string::npos ) treeSelection.append(" && ").append(exprDiTauKine_Sig);
  else if ( region.find("2") != std::string::npos ) treeSelection.append(" && ").append(exprDiTauKine_Bgr);
  
  std::vector<std::string> observables = getObservables(region, fitVariables);
  std::vector<std::string> tauIdValues = getTauIdValues(region);

  makeHistograms(process, 
		 histograms,
		 region, weight,
		 tree, treeSelection,
		 tauId, tauIdValues, observables,
		 branchNames, applyPUreweighting,
		 sysShift,
		 saveRunLumiSectionEventNumbers);
}

void makeDistributionsAllRegions(
  const std::string& process, 
  std::map<std::string, std::map<std::string, TH1*> >& histograms,
  double weight,
  TTree* tree, const std::vector<std::string>& regions,
  const std::vector<std::string>& tauIds, const std::vector<std::string>& fitVariables,
  std::map<std::string, std::map<std::string, std::string> >& branchNames, 
  bool applyPUreweighting,
  const std::string& sysShift = "CENTRAL_VALUE",
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

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {	
    for ( std::vector<std::string>::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {

      bool saveRunLumiSectionEventNumbers_region = ( saveRunLumiSectionEventNumbers && tauId == tauIds.begin() ) ?
	(*saveRunLumiSectionEventNumbers)[*region] : false;

      makeDistributionsInRegion(process, 
				histograms[*region],
				*region, weight,
				tree, 
				*tauId, fitVariables,
				branchNames[*tauId], applyPUreweighting,
				sysShift,
				saveRunLumiSectionEventNumbers_region);
    }
  }
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

int main(int argc, const char* argv[])
{
  std::cout << "<makeTauIdEffPlots>:" << std::endl;

  // CV: TTree::Draw causes segmentation violation in case no TCanvas exists ?!
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);

  gROOT->SetBatch(true);

//--- keep track of time it took the macro to execute
  TBenchmark clock;
  clock.Start("makeTauIdEffPlots");

  //std::string inputFilePath = "/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06_v2_pureweight/"; // for Mauro
  //const std::string jobId = "2011Jun06V2";  
  std::string inputFilePath = // for Christian
    //"/data2/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/2011Jun10/user/m/mverzett/tagprobe/Jun06Skim/edntuples_v3/"; // for Christian
    "/data2/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/2011Jun18/user/v/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/";
  //const std::string jobId = "2011Jun06V2";
  const std::string jobId = "2011Jun18V1";  

  //const std::string branchName_suffix = "local";
  const std::string branchName_suffix = "lxbatch";

  bool runQuickTest = false;
  //bool runQuickTest = true;

  const std::string histogramFileName = "fitTauIdEff_wConstraints_2011June18V1_PUreweighted.root";

  //bool applyPUreweightingMC = false;
  bool applyPUreweightingMC = true;

  std::vector<std::string> sysUncertainties;
  sysUncertainties.push_back(std::string("SysTauJetEnUp"));   // needed for diTauVisMassFromJet
  sysUncertainties.push_back(std::string("SysTauJetEnDown"));
  sysUncertainties.push_back(std::string("SysJetEnUp"));      // needed for diTauMt
  sysUncertainties.push_back(std::string("SysJetEnDown"));

  bool runSysUncertainties = false;
  //bool runSysUncertainties = true;

  std::vector<std::string> regions;
  regions.push_back(std::string("ABCD"));
  regions.push_back(std::string("A"));
  regions.push_back(std::string("B"));
  regions.push_back(std::string("B1"));  // QCD enriched control region (SS, loose muon isolation, Mt && Pzeta cuts applied)
  regions.push_back(std::string("C"));
  regions.push_back(std::string("C1"));
  regions.push_back(std::string("C1p"));
  regions.push_back(std::string("C1f"));
  regions.push_back(std::string("C2"));
  regions.push_back(std::string("C2p"));
  regions.push_back(std::string("C2f"));
  regions.push_back(std::string("D"));   // generic background control region (SS, tight muon isolation)
  regions.push_back(std::string("D1"));
  regions.push_back(std::string("D1p"));
  regions.push_back(std::string("D1f"));
  regions.push_back(std::string("D2"));
  regions.push_back(std::string("D2p"));
  regions.push_back(std::string("D2f"));

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

  double corrFactorData = 0.971*0.993*0.975;                     // Data/MC correction for muon trigger, id. and isolation efficiencies
                                                                 // (numbers taken from CMS AN-2011/153)
  double dataIntLumi = 191*corrFactorData;

  //std::string sampleZtautau = "ZtautauPU156bx";
  //double weightFactorZtautau = dataIntLumi*1666/2568490;       // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  //double corrFactorZtautau = 1.000*1.000;                      // first  number: correction for event looses during skimming
                                                                 // second number: correction for event looses during Ntuple production

  std::string sampleZtautau = "Ztautau_powheg";
  double weightFactorZtautau = dataIntLumi*1666/1995369;
  double corrFactorZtautau = 1.995369e+06 / 1.99537e+06;
  std::string sampleZmumu = "Zmumu_powheg";
  double weightFactorZmumu = dataIntLumi*1666/1998931;         // Z --> l+ l- xSection (FEWZ @ NNLO) / numEvents (POWHEG sample)
  double corrFactorZmumu = 1.984154e+06 / 1.97413e+06;
  std::string sampleQCD = "PPmuXptGt20Mu15";
  double weightFactorQCD = dataIntLumi*(0.2966*1.e+09*2.855e-4) / 29434562; // xSection (LO) / numEvents (PYTHIA PPmuXptGt20Mu15 sample)
  double corrFactorQCD = 29434562 / 2.71458e+07;
  std::string sampleWplusJets = "WplusJets_madgraph";
  double weightFactorWplusJets = dataIntLumi*31314/15168266;   // W --> l nu xSection (FEWZ @ NNLO) / numEvents (MadGraph sample)
  double corrFactorWplusJets = 1.5110974e+07 / 1.18019e+07;
  std::string sampleTTplusJets = "TTplusJets_madgraph";
  double weightFactorTTplusJets = dataIntLumi*157.5/1164640;   // inclusive TTbar xSection (MCFM @ NLO) / numEvents (MadGraph sample)
  double corrFactorTTplusJets = 1.164208e+06 / 1.16421e+06;

  std::vector<std::string> tauIds;
  //tauIds.push_back(std::string("tauDiscrTaNCfrOnePercent"));  // "old" TaNC algorithm
  //tauIds.push_back(std::string("tauDiscrTaNCfrHalfPercent"));
  //tauIds.push_back(std::string("tauDiscrTaNCfrQuarterPercent"));
  //tauIds.push_back(std::string("tauDiscrTaNCloose")); // "new" TaNC implemented in HPS+TaNC combined algorithm
  //tauIds.push_back(std::string("tauDiscrTaNCmedium"));
  //tauIds.push_back(std::string("tauDiscrTaNCtight"));
  //tauIds.push_back(std::string("tauDiscrIsolationLoose"));    // "old" HPS algorithm
  //tauIds.push_back(std::string("tauDiscrIsolationMedium"));   
  //tauIds.push_back(std::string("tauDiscrIsolationTight"));
  tauIds.push_back(std::string("tauDiscrHPSloose"));  // "new" HPS implemented in HPS+TaNC combined algorithm
  //tauIds.push_back(std::string("tauDiscrHPSmedium"));
  //tauIds.push_back(std::string("tauDiscrHPStight"));

  std::vector<std::string> fitVariables;
  //fitVariables.push_back("diTauHt");
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

  std::map<std::string, std::map<std::string, TH1*> > distributionsData;  // key = (region, observable)
  std::map<std::string, std::map<std::string, TH1*> > templatesZtautau;  
  std::map<std::string, std::map<std::string, TH1*> > templatesZmumu;    
  std::map<std::string, std::map<std::string, TH1*> > templatesQCD; 
  std::map<std::string, std::map<std::string, TH1*> > templatesWplusJets; 
  std::map<std::string, std::map<std::string, TH1*> > templatesTTplusJets; 

  for ( std::vector<std::string>::const_iterator sysShift = sysShifts.begin();
	sysShift != sysShifts.end(); ++sysShift ) {
    std::cout << "filling histograms for sysShift = " << (*sysShift) << "..." << std::endl;

//--- initialize alias --> branchName mapping
    typedef std::map<std::string, std::string> branchNameDictEntry;
    std::map<std::string, branchNameDictEntry> branchNamesData           
      = makeBranchNameDict(tauIds, "CENTRAL_VALUE", branchName_suffix);
    std::map<std::string, branchNameDictEntry> branchNamesMC
      = makeBranchNameDict(tauIds, *sysShift,       branchName_suffix);
   
//--- define x-axis titles
    std::map<std::string, std::string> xAxisTitles;
    xAxisTitles["diTauCharge"]         = "Charge(#mu + #tau_{had})";
    xAxisTitles["diTauMt"]             = "M_{T}^{#muMET} [GeV]";
    xAxisTitles["diTauHt"]             = "P_{T}^{#mu} + P_{T}^{#tau} + MET [GeV]";
    xAxisTitles["diTauVisMass"]        = "M_{vis}^{#mu#tau} [GeV]";
    xAxisTitles["diTauVisMassFromJet"] = xAxisTitles["diTauVisMass"];

    TChain* chainData_2011RunA_v1 = new TChain("Events");
    TChain* chainData_2011RunA_v2 = new TChain("Events");
    if ( !runQuickTest ) { addFileNames(chainData_2011RunA_v1, inputFilePath, "data_SingleMu_Run2011A_PromptReco_v1", jobId);
                           addFileNames(chainData_2011RunA_v2, inputFilePath, "data_SingleMu_Run2011A_PromptReco_v2", jobId); } 
    else                 { chainData_2011RunA_v1->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_SingleMu_Run2011A_PromptReco_v1_2011Jun06V2_0_cab9.root").data());
                           chainData_2011RunA_v2->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_data_SingleMu_Run2011A_PromptReco_v2_2011Jun06V2_0_7893.root").data()); }
    printFileInfo(chainData_2011RunA_v1, "chainData_2011RunA_v1");
    printFileInfo(chainData_2011RunA_v2, "chainData_2011RunA_v2");
    TChain* chainData = new TChain("Events");
    chainData->Add(chainData_2011RunA_v1);
    chainData->Add(chainData_2011RunA_v2);
    printFileInfo(chainData, "chainData");
    makeDistributionsAllRegions("Data", 
                                distributionsData, 1.0, chainData, regions,
                                tauIds, fitVariables, branchNamesData, false, *sysShift, &saveRunLumiSectionEventNumbers);
    delete chainData;
    delete chainData_2011RunA_v1;
    delete chainData_2011RunA_v2;
 
    TChain* chainZtautau = new TChain("Events");
    if ( !runQuickTest ) addFileNames(chainZtautau, inputFilePath, sampleZtautau, jobId);
    else                 chainZtautau->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Ztautau_powheg_2011Jun18V1_0_ba35.root").data());
    printFileInfo(chainZtautau, "chainZtautau");
    makeDistributionsAllRegions("Ztautau", 
                                templatesZtautau, weightFactorZtautau*corrFactorZtautau, chainZtautau, regions, 
				tauIds, fitVariables, branchNamesMC, applyPUreweightingMC, *sysShift);
    delete chainZtautau;

    std::map<std::string, std::map<std::string, TH1*> > templatesZmumu; // key = (region, observable)

    TChain* chainZmumu = new TChain("Events");
    if ( !runQuickTest ) addFileNames(chainZmumu, inputFilePath, sampleZmumu, jobId);
    else                 chainZmumu->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_Zmumu_powheg_2011Jun06V2_0_f76a.root").data());
    printFileInfo(chainZmumu, "chainZmumu");
    makeDistributionsAllRegions("Zmumu", 
                                templatesZmumu, weightFactorZmumu*corrFactorZmumu, chainZmumu, regions, 
                                tauIds, fitVariables, branchNamesMC, applyPUreweightingMC, *sysShift);
    delete chainZmumu;

    TChain* chainQCD = new TChain("Events");
    if ( !runQuickTest ) addFileNames(chainQCD, inputFilePath, sampleQCD, jobId);
    else                 chainQCD->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V2_0_ed96.root").data());
    printFileInfo(chainQCD, "chainQCD");
    makeDistributionsAllRegions("QCD", 
                                templatesQCD, weightFactorQCD*corrFactorQCD, chainQCD, regions, 
		        	tauIds, fitVariables, branchNamesMC, applyPUreweightingMC, *sysShift);
    delete chainQCD;

    TChain* chainWplusJets = new TChain("Events");
    if ( !runQuickTest ) addFileNames(chainWplusJets, inputFilePath, sampleWplusJets, jobId);
    else                 chainWplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V2_0_1127.root").data());
    printFileInfo(chainWplusJets, "chainWplusJets");
    makeDistributionsAllRegions("WplusJets", 
                                templatesWplusJets, weightFactorWplusJets*corrFactorWplusJets, chainWplusJets, regions, 
				tauIds, fitVariables, branchNamesMC, applyPUreweightingMC, *sysShift);
    delete chainWplusJets;

    TChain* chainTTplusJets = new TChain("Events");
    if ( !runQuickTest ) addFileNames(chainTTplusJets, inputFilePath, sampleTTplusJets, jobId);
    else                 chainTTplusJets->Add(std::string(inputFilePath).append("tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V2_0_5a68.root").data());
    printFileInfo(chainTTplusJets, "chainTTplusJets");
    makeDistributionsAllRegions("TTplusJets", 
                                templatesTTplusJets, weightFactorTTplusJets*corrFactorTTplusJets, chainTTplusJets, regions, 
				tauIds, fitVariables, branchNamesMC, applyPUreweightingMC, *sysShift);
    delete chainTTplusJets;
  }
    
//--- save histograms
  TFile* histogramOutputFile = new TFile(histogramFileName.data(), "RECREATE");
  saveHistograms(histogramOutputFile, distributionsData);
  saveHistograms(histogramOutputFile, templatesZtautau);
  saveHistograms(histogramOutputFile, templatesZmumu);
  saveHistograms(histogramOutputFile, templatesQCD);
  saveHistograms(histogramOutputFile, templatesWplusJets);
  saveHistograms(histogramOutputFile, templatesTTplusJets);
  delete histogramOutputFile;

  delete canvas;

//--print time that it took macro to run
  std::cout << "finished executing makeTauIdEffPlots macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  std::cout << " #sysShifts    = " << sysShifts.size() << std::endl;
  clock.Show("makeTauIdEffPlots");
}
