
#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"
#include "TauAnalysis/FittingTools/interface/templateFitAuxFunctions.h"

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

RooFormulaVar* makeRooFormulaVar(const std::string& process, const std::string& region,
				 RooAbsReal* norm, double fittedFractionValue,
				 RooAbsReal* pTauId_passed_failed)
{
  std::cout << "<makeRooFormulaVar>:" << std::endl;
  std::cout << " building RooFormulaVar expression for process = " << process << ", region = " << region << "." << std::endl;

  RooFormulaVar* retVal = 0;

  std::string exprTauId = "";
  
  if      ( region.find("p") != std::string::npos ) exprTauId     = "regular";
  else if ( region.find("f") != std::string::npos ) exprTauId     = "inverted";

  std::string fittedFractionName = std::string(norm->GetName()).append("_").append(region).append("_fittedFraction");
  RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), fittedFractionValue);

  std::string formula = "";
  TObjArray arguments; 
  addToFormula(formula, exprTauId,       arguments, pTauId_passed_failed);
  addToFormula(formula, "regular",       arguments, norm);
  //addToFormula(formula, "regular",       arguments, fittedFraction);

  std::cout << " formula = " << formula << std::endl;
  //std::cout << " arguments:" << std::endl;
  //for ( int i = 0; i < arguments.GetEntries(); ++i ) {
  //  std::cout << "  " << dynamic_cast<TNamed*>(arguments.At(i))->GetName() << ":" 
  //	        << " value = " << dynamic_cast<RooAbsReal*>(arguments.At(i))->getVal() << std::endl;
  //}

  std::string name = std::string("norm").append(process).append("_").append(region);
  retVal = new RooFormulaVar(name.data(), name.data(), formula.data(), RooArgSet(arguments));

  return retVal;
}

double getNormInRegion(std::map<std::string, RooAbsReal*> normC1, 
		       std::map<std::string, RooRealVar*> pTauId_passed_failed,
		       const std::string& process, const std::string& region)
{
  double retVal = normC1[process]->getVal();
  
  if      ( region.find("p") != std::string::npos ) retVal *= pTauId_passed_failed[process]->getVal();
  else if ( region.find("f") != std::string::npos ) retVal *= (1. - pTauId_passed_failed[process]->getVal());

  return retVal;
}

void fitUsingRooFit(std::map<std::string, std::map<std::string, TH1*> >& distributionsData,                         // key = (region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,      // key = (process, region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll,    // key = (process/"sum", region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions, // key = (process, region, observable)
		    const std::vector<std::string>& processes,
		    const std::string& tauId, const std::string& fitVariable, 
		    double& effValue, double& effError, bool& hasFitConverged,
		    bool makeControlPlots, std::map<std::string, std::string>* xAxisTitles = NULL) 
{
  //std::cout << "<fitUsingRooFit>:" << std::endl;
  //std::cout << " performing Fit of variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
      
  double fitMinC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmin();
  double fitMaxC1p = templatesAll["Ztautau"]["C1p"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmax();
  double fitMinC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmin();
  double fitMaxC1f = templatesAll["Ztautau"]["C1f"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmax();
  assert(fitMinC1p == fitMinC1f && fitMaxC1p == fitMaxC1f);
  RooRealVar* fitVarC1 = new RooRealVar("fitVarC1", "fitVarC1", fitMinC1p, fitMaxC1p);

  double numEventsDataC1 = distributionsData["C1"][getKey(fitVariable, tauId)]->Integral();
  if ( makeControlPlots ) std::cout << "numEventsDataC1 = " << numEventsDataC1 << std::endl;

  std::map<std::string, RooRealVar*> pTauId_passed_failed;        // key = process

  std::map<std::string, RooAbsReal*> normC1;                      // key = process
  std::map<std::string, RooAbsReal*> normC1p;                     // key = process
  std::map<std::string, RooAbsReal*> normC1f;                     // key = process

  std::map<std::string, std::map<std::string, RooHistPdf*> > pdf; // key = (process/"sum", region)

  TObjArray pdfsC1p;
  TObjArray pdfsC1f;

  TObjArray fitParametersC1p;
  TObjArray fitParametersC1f;

   for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
     //std::cout << "process = " << (*process) << ":" << std::endl;

    double numEventsC1        = numEventsAll[*process]["C1"][getKey(fitVariable, tauId)];
    double fittedFractionC1   = fittedFractions[*process]["C1"][getKey(fitVariable, tauId)];
    double fittedEventsC1     = fittedFractionC1*numEventsC1;
    double numEventsC1p       = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedFractionC1p  = fittedFractions[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedEventsC1p    = fittedFractionC1p*numEventsC1p;
    double numEventsC1f       = numEventsAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedFractionC1f  = fittedFractions[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedEventsC1f    = fittedFractionC1f*numEventsC1f;
    ///std::cout << " numEventsC1 = " << numEventsC1 << ", fittedFractionC1 = " << fittedFractionC1 << std::endl;

    std::string nameTauId_passed_failed = std::string("pTauId_passed_failed").append("_").append(*process);
    double pTauId_passed_failed0        = fittedEventsC1p/fittedEventsC1;
    pTauId_passed_failed[*process]      = new RooRealVar(nameTauId_passed_failed.data(),
							 nameTauId_passed_failed.data(), pTauId_passed_failed0, 0., 1.);

    double numEventsSumC1 = numEventsAll["sum"]["C1"][getKey("diTauMt", tauId)];
    //std::cout << " numEventsSumC1 = " << numEventsSumC1 << std::endl;
    
    //double scaleFactorMCtoData = numEventsDataC1/numEventsSumC1;
    double scaleFactorMCtoData = 1.;
    //std::cout << "--> MC-to-Data scale-factor = " << scaleFactorMCtoData << std::endl;

    std::string nameNormC1 = std::string("normC1").append("_").append(*process);
    double normC1initial         = scaleFactorMCtoData*numEventsC1;
    normC1[*process]       = new RooRealVar(nameNormC1.data(), nameNormC1.data(), normC1initial, 0., numEventsDataC1);

    TH1* templateC1p   = templatesAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    RooHistPdf* pdfC1p = makeRooHistPdf(templateC1p, fitVarC1);
    TH1* templateC1f   = templatesAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    RooHistPdf* pdfC1f = makeRooHistPdf(templateC1f, fitVarC1);    

    pdfsC1p.Add(pdfC1p);
    pdfsC1f.Add(pdfC1f);

    normC1p[*process] = makeRooFormulaVar(*process, "C1p", normC1[*process], fittedFractionC1p, pTauId_passed_failed[*process]);
    normC1f[*process] = makeRooFormulaVar(*process, "C1f", normC1[*process], fittedFractionC1f, pTauId_passed_failed[*process]);

    fitParametersC1p.Add(normC1p[*process]);
    fitParametersC1f.Add(normC1f[*process]);
  }

  RooAddPdf* pdfSumC1p = new RooAddPdf("pdfSumC1p", "pdfSumC1p", RooArgList(pdfsC1p), RooArgList(fitParametersC1p));
  RooAddPdf* pdfSumC1f = new RooAddPdf("pdfSumC1f", "pdfSumC1f", RooArgList(pdfsC1f), RooArgList(fitParametersC1f));

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

//--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  TObjArray fitConstraintsC1;
  fitConstraintsC1.Add(makeFitConstraint(normC1["Zmumu"],      
					 normC1["Zmumu"]->getVal(),                 2.0*normC1["Zmumu"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normC1["QCD"],        
					 normC1["QCD"]->getVal(),                   1.0*normC1["QCD"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normC1["WplusJets"],  
					 normC1["WplusJets"]->getVal(),             1.0*normC1["WplusJets"]->getVal()));
  fitConstraintsC1.Add(makeFitConstraint(normC1["TTplusJets"], 
					 normC1["TTplusJets"]->getVal(),            1.0*normC1["TTplusJets"]->getVal()));
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
  //fitOptionsC1.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsC1))));
  fitOptionsC1.Add(new RooCmdArg(RooFit::Save(true)));

  //pdfSimultaneousFitC1->printCompactTree();

  RooFitResult*	fitResult = pdfSimultaneousFitC1->fitTo(*dataC1, fitOptionsC1);
   
  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();
  hasFitConverged = (fitResult->status() == 0) ? true : false;

  if ( !makeControlPlots ) return;
  
  std::cout << tauId << ":";
  if ( hasFitConverged ) std::cout << " fit converged."          << std::endl; 
  else                   std::cout << " fit failed to converge." << std::endl;
  
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

  std::map<std::string, double> normFactorC1;   // key = process
  std::map<std::string, double> normFactorC1p;  // key = process
  std::map<std::string, double> normFactorC1f;  // key = process

  std::cout << "Results of fitting variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    double numEventsC1        = numEventsAll[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedFractionC1   = fittedFractions[*process]["C1"][getKey("diTauMt", tauId)];
    double fittedEventsC1     = fittedFractionC1*numEventsC1;
    double numEventsC1p       = numEventsAll[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedFractionC1p  = fittedFractions[*process]["C1p"][getKey(fitVariable, tauId, "passed")];
    double fittedEventsC1p    = fittedFractionC1p*numEventsC1p;
    double numEventsC1f       = numEventsAll[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedFractionC1f  = fittedFractions[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double fittedEventsC1f    = fittedFractionC1f*numEventsC1f;
    
    std::cout << " " << (*process) << ":" << std::endl;
    std::cout << "  normalization = " << normC1[*process]->getVal()
	      << " +/- " << dynamic_cast<RooRealVar*>(normC1[*process])->getError()
	      << " (MC exp. = " << numEventsC1 << ")" << std::endl;
    double pTauId_passed_failedMCexp = fittedEventsC1p/fittedEventsC1;
    std::cout << "  pTauId_passed_failed = " << pTauId_passed_failed[*process]->getVal() 
	      << " +/- " << pTauId_passed_failed[*process]->getError() 
	      << " (MC exp. = " << pTauId_passed_failedMCexp << ")" << std::endl;
    
    normFactorC1[*process]  = getNormInRegion(normC1, pTauId_passed_failed, *process, "C1");
    normFactorC1p[*process] = getNormInRegion(normC1, pTauId_passed_failed, *process, "C1p");
    normFactorC1f[*process] = getNormInRegion(normC1, pTauId_passed_failed, *process, "C1f");
    
    std::cout << "--> C1 = " << normFactorC1[*process]
	      << " (MC exp. = " << numEventsC1 << ")" << std::endl;
    std::cout << "--> C1p = " << normFactorC1p[*process]
	      << " (MC exp. = " << numEventsC1p << ")" << std::endl;
    std::cout << "--> C1f = " << normFactorC1f[*process]
	      << " (MC exp. = " << numEventsC1f << ")" << std::endl;
  }  
 
//--- make control plots for sum(MC) scaled by normalization determined by fit versus Data 
//    for Mt, fitVariable distributions in different regions
  drawHistograms(distributionsData, templatesAll,  
		 normFactorC1, "C1", getKey(fitVariable, tauId),		 
		 std::string("Region C1: ").append(fitVariable).append(" (scaled by normalization det. by fit)"), 
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff2_C1_").append(fitVariable).append("_fitted.png"),
		 true);
  drawHistograms(distributionsData, templatesAll,  
		 normFactorC1p, "C1p", getKey(fitVariable, tauId, "passed"),
		 std::string("Region C1p: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff2_C1p_").append(fitVariable).append("_").append(tauId).append("_fitted.png"),
		 true);
  drawHistograms(distributionsData, templatesAll,  
		 normFactorC1f, "C1f", getKey(fitVariable, tauId, "failed"),
		 std::string("Region C1f: ").append(fitVariable).append(", ").append(tauId).append(" (scaled by normalization det. by fit)"),
		 xAxisTitles ? (*xAxisTitles)[fitVariable] : "",
		 std::string("controlPlotsTauIdEff2_C1f_").append(fitVariable).append("_").append(tauId).append("_fitted.png"),
		 true);
}

int main(int argc, const char* argv[])
{
  std::cout << "<fitTauIdEff>:" << std::endl;

  gROOT->SetBatch(true);

//--- keep track of time it took the macro to execute
  TBenchmark clock;
  clock.Start("fitTauIdEff");

  bool runClosureTest = false;
  //bool runClosureTest = true;

  //const std::string histogramFileName = "fitTauIdEff_wConstraints_2011June10.root";
  const std::string histogramFileName = "fitTauIdEff_wConstraints_2011June06_PUreweighted.root";

  std::vector<sysUncertaintyEntry> sysUncertainties;
  sysUncertainties.push_back(sysUncertaintyEntry("SysTauJetEnUp", "SysTauJetEnDown", "SysTauJetEnDiff"));
  //sysUncertainties.push_back(sysUncertaintyEntry("SysJetEnUp", "SysJetEnDown", "SysJetEnDiff"));

  //bool runSysUncertainties = false;
  bool runSysUncertainties = true;

  unsigned numPseudoExperiments = 10000;  

  std::vector<std::string> regions;
  regions.push_back(std::string("C1"));
  regions.push_back(std::string("C1p"));
  regions.push_back(std::string("C1f"));
  regions.push_back(std::string("B1")); // needed to access QCD template obtained from Data

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
  if ( runSysUncertainties ) {
    for ( std::vector<sysUncertaintyEntry>::const_iterator sysUncertainty = sysUncertainties.begin();
	  sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
      sysShifts.push_back(sysUncertainty->sysNameUp_);
      sysShifts.push_back(sysUncertainty->sysNameDown_);
    }
  }

  typedef std::map<std::string, std::map<std::string, TH1*> > histogramMap; // key = (region, observable + sysShift)
  histogramMap distributionsData;

  histogramMap templatesZtautau;
  histogramMap templatesZmumu;
  histogramMap templatesQCD;
  histogramMap templatesWplusJets;
  histogramMap templatesTTplusJets;

  TFile* histogramInputFile = new TFile(histogramFileName.data());

  for ( std::vector<std::string>::const_iterator sysShift = sysShifts.begin();
	sysShift != sysShifts.end(); ++sysShift ) {
    std::cout << "loading histograms for sysShift = " << (*sysShift) << "..." << std::endl;
  
    loadHistograms(distributionsData,   histogramInputFile, "Data",       regions, tauIds, fitVariables, *sysShift);

    loadHistograms(templatesZtautau,    histogramInputFile, "Ztautau",    regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesZmumu,      histogramInputFile, "Zmumu",      regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesQCD,        histogramInputFile, "QCD",        regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesWplusJets,  histogramInputFile, "WplusJets",  regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesTTplusJets, histogramInputFile, "TTplusJets", regions, tauIds, fitVariables, *sysShift);
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

//--- compute systematic shifts 
  for ( std::vector<sysUncertaintyEntry>::const_iterator sysUncertainty = sysUncertainties.begin();
	sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
    compSysHistograms(templatesAll, *sysUncertainty);
  }

//--- closure test: fit sum(MC) instead of Data
  if ( runClosureTest ) {
    std::cout << ">>> NOTE: RUNNING CLOSURE TEST <<<" << std::endl;

    double mcToDataScaleFactor = 1.; // run closure test for event statistics of current dataset
    //double mcToDataScaleFactor = 10.; // run closure test assuming event statistics to be higher by given factor
    
    sumHistograms(templatesAll, processes, "sum", mcToDataScaleFactor);
    
    std::cout << std::endl;
    
    distributionsData = templatesAll["sum"];
  }

/*
//--- obtain template for QCD background from data,
//    from SS && Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV sideband
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

//--- define x-axis titles
  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["diTauCharge"]         = "Charge(#mu + #tau_{had})";
  xAxisTitles["diTauMt"]             = "M_{T}^{#muMET} [GeV]";
  xAxisTitles["diTauHt"]             = "P_{T}^{#mu} + P_{T}^{#tau} + MET [GeV]";
  xAxisTitles["diTauVisMass"]        = "M_{vis}^{#mu#tau} [GeV]";
  xAxisTitles["diTauVisMassFromJet"] = xAxisTitles["diTauVisMass"];

//--- make control plots for sum(MC) scaled by cross-sections versus Data 
//    for Mt, fitVariable distributions in different regions
  for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	region != distributionsData.end(); ++region ) {
    for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	  key != region->second.end(); ++key ) {
      if ( !isSystematicShift(key->first) )
	drawHistograms(templatesZtautau[region->first][key->first], -1.,
		       templatesZmumu[region->first][key->first], -1.,
		       templatesQCD[region->first][key->first], -1.,
		       templatesWplusJets[region->first][key->first], -1.,
		       templatesTTplusJets[region->first][key->first], -1.,
		       distributionsData[region->first][key->first],
		       std::string("Region ").append(region->first).append(": ").append(key->first).append(" (scaled by cross-section)"),
		       xAxisTitles[key->first],
		       std::string("controlPlotsTauIdEff_").append(region->first).append("_").append(key->first).append(".png"));
    }
  }

  std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll = // key = (process/"sum", region, observable)
    compNumEvents(templatesAll, processes, distributionsData);
  std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions = // key = (process, region, observable)
    compFittedFractions(templatesAll, numEventsAll, processes, distributionsData);

  std::cout << "running fit for central values..." << std::endl;
  
  std::map<std::string, std::map<std::string, double> > effValues;           // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > effErrors;           // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, bool> >   fitConvergenceFlags; // key = (tauId, fitVariable)
  
  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      double effValue = 0.;
      double effError = 1.;
      bool hasFitConverged = false;
      fitUsingRooFit(distributionsData, templatesAll, numEventsAll, fittedFractions,
		     processes,
		     *tauId, *fitVariable, 
		     effValue, effError, hasFitConverged,	       
		     true, &xAxisTitles);
      effValues[*tauId][*fitVariable] = effValue;
      effErrors[*tauId][*fitVariable] = effError;
      fitConvergenceFlags[*tauId][*fitVariable] = hasFitConverged;
    }
  }
      
  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    std::cout << "Efficiency of Tau id. = " << (*tauId) << ":" << std::endl;
    
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      std::cout << " fitVariable = " << (*fitVariable) << ":" 
		<< " result = " << effValues[*tauId][*fitVariable]*100. 
		<< " +/- " << effErrors[*tauId][*fitVariable]*100. << "%";
      if ( fitConvergenceFlags[*tauId][*fitVariable] ) std::cout << " (fit converged successfully)";
      else std::cout << " (fit failed to converge)";
      std::cout << std::endl;      
      
      double numEventsC1   = numEventsAll["Ztautau"]["C1"][getKey("diTauMt", *tauId)];
      double numEventsC1p  = numEventsAll["Ztautau"]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
      
      double tauIdEffMCexp = numEventsC1p/numEventsC1;
      std::cout << "(Monte Carlo prediction = " << tauIdEffMCexp*100. << "%)" << std::endl;
    }
  }

  if ( runSysUncertainties ) {
    
    std::map<std::string, std::map<std::string, TH1*> > effDistributions;            // key = (tauId, fitVariable)
    std::map<std::string, std::map<std::string, TH1*> > fitConvergenceDistributions; // key = (tauId, fitVariable)
    
    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >   templatesAll_fluctuated;
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll_fluctuated;
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions_fluctuated;

    for ( unsigned iPseudoExperiment = 0; iPseudoExperiment < numPseudoExperiments; ++numPseudoExperiments ) {
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	      region != distributionsData.end(); ++region ) {
	  for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
		key != region->second.end(); ++key ) {
	    if ( !isSystematicShift(key->first) ) {
	      TH1* origHistogram = key->second;
	      
	      TH1* fluctHistogram = templatesAll_fluctuated[*process][region->first][key->first];
	      if ( !fluctHistogram ) {
		fluctHistogram = (TH1*)origHistogram->Clone(TString(origHistogram->GetName()).Append("_fluctuated"));
		templatesAll_fluctuated[*process][region->first][key->first] = fluctHistogram;
	      }
	      
	      sampleHistogram_stat(origHistogram, fluctHistogram);
	      
	      for ( std::vector<sysUncertaintyEntry>::const_iterator sysUncertainty = sysUncertainties.begin();
		    sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
		std::string key_sys = std::string(key->first).append("_").append(sysUncertainty->sysNameDiff_);
		TH1* sysHistogram = templatesAll_fluctuated[*process][region->first][key_sys];
		assert(sysHistogram);
		
		sampleHistogram_sys(fluctHistogram, sysHistogram, 1.0, -3.0, +3.0, kCoherent);
	      }
	    }
	  }
	}
      }
      
      numEventsAll_fluctuated = compNumEvents(templatesAll_fluctuated, processes, distributionsData);
      fittedFractions_fluctuated = compFittedFractions(templatesAll_fluctuated, numEventsAll_fluctuated, processes, distributionsData);
      
      for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	    tauId != tauIds.end(); ++tauId ) {
	for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	      fitVariable != fitVariables.end(); ++fitVariable ) {
	  double effValue = 0.;
	  double effError = 1.;
	  bool hasFitConverged = false;
	  fitUsingRooFit(distributionsData, templatesAll_fluctuated, numEventsAll_fluctuated, fittedFractions_fluctuated,
			 processes,
			 *tauId, *fitVariable, 
			 effValue, effError, hasFitConverged,	       
			 false, &xAxisTitles);

	  TH1* effDistribution = effDistributions[*tauId][*fitVariable];
	  if ( !effDistribution ) {
            TString effDistributionName  = Form("effDistribution_%s_%s", tauId->data(), fitVariable->data());
	    TString effDistributionTitle = Form("ToyMC: %s Efficiency for %s", tauId->data(), fitVariable->data());
	    effDistribution = new TH1F(effDistributionName.Data(), effDistributionTitle.Data(), 101, -0.005, +1.005);
	    effDistributions[*tauId][*fitVariable] = effDistribution;
	  }

	  effDistribution->Fill(effValue);
	  
	  TH1* fitConvergenceDistribution = fitConvergenceDistributions[*tauId][*fitVariable];
	  if ( !fitConvergenceDistribution ) {
            TString fitConvergenceDistributionName  = Form("fitConvergenceDistribution_%s_%s", tauId->data(), fitVariable->data());
	    TString fitConvergenceDistributionTitle = Form("ToyMC: %s Efficiency for %s", tauId->data(), fitVariable->data());
	    fitConvergenceDistribution = 
	      new TH1F(fitConvergenceDistributionName.Data(), fitConvergenceDistributionTitle.Data(), 2, -0.005, +1.005);
	    TAxis* xAxis = fitConvergenceDistribution->GetXaxis();
	    xAxis->SetBinLabel(1, "Failure");
	    xAxis->SetBinLabel(2, "Success");
	    fitConvergenceDistributions[*tauId][*fitVariable] = fitConvergenceDistribution;
	  }
	  
	  fitConvergenceDistribution->Fill(hasFitConverged);
	}
      }
    }

    savePseudoExperimentHistograms(effDistributions,            "#varepsilon", ".png");
    savePseudoExperimentHistograms(fitConvergenceDistributions, "Fit status",  ".png");
  }

  delete histogramInputFile;

//--print time that it took macro to run
  std::cout << "finished executing fitTauIdEff macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  std::cout << " #sysShifts    = " << sysShifts.size() 
	    << "(numPseudoExperiments = " << numPseudoExperiments << ")" << std::endl;
  clock.Show("fitTauIdEff");
}
  
 
