
/** \executable fitTauIdEff
 *
 * Determine tau identification efficiency by simultaneous fit of Ztautau signal 
 * plus Zmumu, W + jets, QCD and TTbar background yields in "passed" and "failed" regions 
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.20 $
 *
 * $Id: fitTauIdEff.cc,v 1.20 2011/07/22 17:24:56 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
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
#include "RooIntegralMorph.h"
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
  if      ( region.find("PASSED") != std::string::npos ) exprTauId = "regular";
  else if ( region.find("FAILED") != std::string::npos ) exprTauId = "inverted";

  std::string fittedFractionName = std::string(norm->GetName()).append("_").append(region).append("_fittedFraction");
  //RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), fittedFractionValue);
  // CV: fitted yields for all processes come-out factor 2 too high in closure test
  //    --> compensate by multiplying fittedFraction by fudge-factor 2
  RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), 2.*fittedFractionValue);

  std::string formula = "";
  TObjArray arguments; 
  addToFormula(formula, exprTauId, arguments, pTauId_passed_failed);
  addToFormula(formula, "regular", arguments, norm);
  addToFormula(formula, "regular", arguments, fittedFraction);

  std::cout << " formula = " << formula << std::endl;
  std::cout << " arguments:" << std::endl;
  for ( int i = 0; i < arguments.GetEntries(); ++i ) {
    std::cout << "  " << dynamic_cast<TNamed*>(arguments.At(i))->GetName() << ":" 
	      << " value = " << dynamic_cast<RooAbsReal*>(arguments.At(i))->getVal() << std::endl;
  }

  std::string name = std::string("norm").append(process).append("_").append(region);
  retVal = new RooFormulaVar(name.data(), name.data(), formula.data(), RooArgSet(arguments));

  return retVal;
}

double getNormInRegion(std::map<std::string, RooRealVar*> norm, 
		       std::map<std::string, RooRealVar*> pTauId_passed_failed,
		       const std::string& process, const std::string& region)
{
  double retVal = norm[process]->getVal();
  
  if      ( region.find("PASSED") != std::string::npos ) retVal *= pTauId_passed_failed[process]->getVal();
  else if ( region.find("FAILED") != std::string::npos ) retVal *= (1. - pTauId_passed_failed[process]->getVal());

  return retVal;
}

void fitUsingRooFit(std::map<std::string, std::map<std::string, TH1*> >& distributionsData,                         // key = (region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,      // key = (process, region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll,    // key = (process/"sum", region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions, // key = (process, region, observable)
		    const std::vector<std::string>& processes,
		    const std::string& tauId, const std::string& fitVariable, 		    
		    const std::string& morphSysUncertaintyUp, const std::string& morphSysUncertaintyDown, double sysVariedByNsigma,
		    double& effValue, double& effError, 
		    std::map<std::string, std::map<std::string, double> >& normFactors_fitted, bool& hasFitConverged,
		    int verbosity = 0)
{
  if ( verbosity ) {
    std::cout << "<fitUsingRooFit>:" << std::endl;
    std::cout << " performing Fit of variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
  }
      
  double fitMin_passed = templatesAll["Ztautau"]["PASSED"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmin();
  double fitMax_passed = templatesAll["Ztautau"]["PASSED"][getKey(fitVariable, tauId, "passed")]->GetXaxis()->GetXmax();
  double fitMin_failed = templatesAll["Ztautau"]["FAILED"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmin();
  double fitMax_failed = templatesAll["Ztautau"]["FAILED"][getKey(fitVariable, tauId, "failed")]->GetXaxis()->GetXmax();
  assert(fitMin_passed == fitMin_failed && fitMax_passed == fitMax_failed);
  std::string fitVarName = std::string("fitVar").append("_").append(fitVariable);
  RooRealVar* fitVar = new RooRealVar(fitVarName.data(), fitVarName.data(), fitMin_passed, fitMax_passed);
  if ( verbosity ) std::cout << "range = " << fitMin_passed << ".." << fitMax_passed << std::endl;

  double numEventsData = distributionsData["ALL"][getKey(fitVariable, tauId)]->Integral();
  if ( verbosity ) std::cout << "numEventsData = " << numEventsData << std::endl;

  std::map<std::string, RooRealVar*> pTauId_passed_failed;       // key = process

  std::map<std::string, RooRealVar*> norm;                       // key = process
  std::map<std::string, RooAbsReal*> norm_passed;                // key = process
  std::map<std::string, RooAbsReal*> norm_failed;                // key = process

  std::map<std::string, std::map<std::string, RooAbsPdf*> > pdf; // key = (process/"sum", region)
  std::vector<RooRealVar*> alphas;

  TObjArray pdfs_passed;
  TObjArray pdfs_failed;

  TObjArray fitParameters_passed;
  TObjArray fitParameters_failed;

  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    //std::cout << "process = " << (*process) << ":" << std::endl;
    
    double numEvents             = numEventsAll[*process]["ALL"][getKey(fitVariable, tauId)];
    double fittedFraction        = fittedFractions[*process]["ALL"][getKey(fitVariable, tauId)];
    double fittedEvents          = fittedFraction*numEvents;
    double numEvents_passed      = numEventsAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
    double fittedFraction_passed = fittedFractions[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
    double fittedEvents_passed   = fittedFraction_passed*numEvents_passed;
    double fittedFraction_failed = fittedFractions[*process]["FAILED"][getKey(fitVariable, tauId, "failed")];
    //std::cout << " numEvents = " << numEvents << ", fittedFraction = " << fittedFraction << std::endl;
    
    std::string pTauIdName_passed_failed = std::string("pTauId_passed_failed").append("_").append(*process);
    double pTauId_passed_failed0 = fittedEvents_passed/fittedEvents;
    //std::cout << " pTauId_passed_failed0 = " << pTauId_passed_failed0 << std::endl;
    pTauId_passed_failed[*process] = 
      new RooRealVar(pTauIdName_passed_failed.data(), 
		     pTauIdName_passed_failed.data(), pTauId_passed_failed0, 0., 1.);
    
    //double numEventsSum = numEventsAll["sum"]["ALL"][getKey("diTauMt", tauId)];
    //std::cout << " numEventsSum = " << numEventsSum << std::endl;
    
    //double scaleFactorMCtoData = numEventsData/numEventsSum;
    double scaleFactorMCtoData = 1.;
    //std::cout << "--> MC-to-Data scale-factor = " << scaleFactorMCtoData << std::endl;
    
    std::string nameNorm = std::string("norm").append("_").append(*process);
    double norm_initial = scaleFactorMCtoData*numEvents;
    //std::cout << " norm_initial = " << norm_initial << std::endl;
    norm[*process] = new RooRealVar(nameNorm.data(), nameNorm.data(), norm_initial, 0., 2.*numEventsData);
    
    TH1* template_passed = templatesAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
    //std::cout << "PASSED: numBins = " << template_passed->GetNbinsX() << "," 
    //          << " integral = " << template_passed->Integral() << std::endl;
    RooAbsPdf* pdf_passed = 0;
    if ( morphSysUncertaintyUp != morphSysUncertaintyDown ) {
      TH1* templateUp_passed = templatesAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed", morphSysUncertaintyUp)];
      RooHistPdf* pdfUp_passed = makeRooHistPdf(templateUp_passed, fitVar);
      TH1* templateDown_passed = templatesAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed", morphSysUncertaintyDown)];
      RooHistPdf* pdfDown_passed = makeRooHistPdf(templateDown_passed, fitVar);
      std::string alphaName_passed = std::string(template_passed->GetName()).append("_alpha");
      RooRealVar* alpha_passed = new RooRealVar(alphaName_passed.data(), alphaName_passed.data(), 0.5, 0., 1.);
      alphas.push_back(alpha_passed);
      std::string pdfName_passed = std::string(template_passed->GetName()).append("_pdf");
      pdf_passed = 
	new RooIntegralMorph(pdfName_passed.data(), 
			     pdfName_passed.data(), *pdfUp_passed, *pdfDown_passed, *fitVar, *alpha_passed, true);
    } else {
      pdf_passed = makeRooHistPdf(template_passed, fitVar);
    }
    TH1* template_failed = templatesAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed")];
    //std::cout << "FAILED: numBins = " << template_failed->GetNbinsX() << "," 
    //          << " integral = " << template_failed->Integral() << std::endl;
    RooAbsPdf* pdf_failed = 0;
    if ( morphSysUncertaintyUp != morphSysUncertaintyDown ) {
      TH1* templateUp_failed = templatesAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed", morphSysUncertaintyUp)];
      RooHistPdf* pdfUp_failed = makeRooHistPdf(templateUp_failed, fitVar);
      TH1* templateDown_failed = templatesAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed", morphSysUncertaintyDown)];
      RooHistPdf* pdfDown_failed = makeRooHistPdf(templateDown_failed, fitVar);
      std::string alphaName_failed = std::string(template_failed->GetName()).append("_alpha");
      RooRealVar* alpha_failed = new RooRealVar(alphaName_failed.data(), alphaName_failed.data(), 0.5, 0., 1.);
      alphas.push_back(alpha_failed);
      std::string pdfName_failed = std::string(template_failed->GetName()).append("_pdf");
      pdf_failed = 
	new RooIntegralMorph(pdfName_failed.data(), 
			     pdfName_failed.data(), *pdfUp_failed, *pdfDown_failed, *fitVar, *alpha_failed, true);
    } else {
      pdf_failed = makeRooHistPdf(template_failed, fitVar);
    }
    
    pdf[*process]["PASSED"] = pdf_passed;
    pdfs_passed.Add(pdf_passed);
    pdf[*process]["FAILED"] = pdf_failed;
    pdfs_failed.Add(pdf_failed);

    norm_passed[*process] = makeRooFormulaVar(*process, "PASSED", norm[*process], fittedFraction_passed, pTauId_passed_failed[*process]);
    norm_failed[*process] = makeRooFormulaVar(*process, "FAILED", norm[*process], fittedFraction_failed, pTauId_passed_failed[*process]);
    
    fitParameters_passed.Add(norm_passed[*process]);
    fitParameters_failed.Add(norm_failed[*process]);
  }

  RooAddPdf* pdfSum_passed = new RooAddPdf("pdfSum_passed", "pdfSum_passed", RooArgList(pdfs_passed), RooArgList(fitParameters_passed));
  RooAddPdf* pdfSum_failed = new RooAddPdf("pdfSum_failed", "pdfSum_failed", RooArgList(pdfs_failed), RooArgList(fitParameters_failed));

//--- build data & model objects for fitting "passed" and "failed" regions
  RooCategory* fitCategories = new RooCategory("categories", "categories");
  fitCategories->defineType("PASSED");
  fitCategories->defineType("FAILED");

  RooSimultaneous* pdfSimultaneousFit = new RooSimultaneous("pdfSimultaneousFit", "pdfSimultaneousFit", *fitCategories);
  pdfSimultaneousFit->addPdf(*pdfSum_passed, "PASSED");
  pdfSimultaneousFit->addPdf(*pdfSum_failed, "FAILED");
 
  std::map<std::string, TH1*> histogramDataMap;
  histogramDataMap["PASSED"] = distributionsData["PASSED"][getKey(fitVariable, tauId, "passed")];
  //std::cout << "PASSED: numBins = " << histogramDataMap["PASSED"]->GetNbinsX() << "," 
  //	      << " integral = " << histogramDataMap["PASSED"]->Integral() << std::endl;
  histogramDataMap["FAILED"] = distributionsData["FAILED"][getKey(fitVariable, tauId, "failed")];
  //std::cout << "FAILED: numBins = " << histogramDataMap["FAILED"]->GetNbinsX() << "," 
  //          << " integral = " << histogramDataMap["FAILED"]->Integral() << std::endl;

  RooDataHist* data = new RooDataHist("data", "data", *fitVar, *fitCategories, histogramDataMap);

//--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  TObjArray fitConstraints;
  //fitConstraints.Add(makeFitConstraint(norm["Zmumu"],      
  //					 norm["Zmumu"]->getVal(),                      1.0*norm["Zmumu"]->getVal()));
  //norm["Zmumu"]->setMax(5.0*norm["Zmumu"]->getVal());
  //fitConstraints.Add(makeFitConstraint(norm["QCD"],        
  //					 norm["QCD"]->getVal(),                        1.0*norm["QCD"]->getVal()));
  //fitConstraints.Add(makeFitConstraint(norm["WplusJets"],  
  //					 norm["WplusJets"]->getVal(),                  1.0*norm["WplusJets"]->getVal()));
  //fitConstraints.Add(makeFitConstraint(norm["TTplusJets"], 
  //		        		 norm["TTplusJets"]->getVal(),                 1.0*norm["TTplusJets"]->getVal()));
  //norm["TTplusJets"]->setMax(5.0*norm["TTplusJets"]->getVal());
  //fitConstraints.Add(makeFitConstraint(pTauId_passed_failed["Zmumu"],      
  //					 pTauId_passed_failed["Zmumu"]->getVal(),      1.0*pTauId_passed_failed["Zmumu"]->getVal()));
  //fitConstraints.Add(makeFitConstraint(pTauId_passed_failed["QCD"],        
  //					 pTauId_passed_failed["QCD"]->getVal(),        1.0*pTauId_passed_failed["QCD"]->getVal()));
  //fitConstraints.Add(makeFitConstraint(pTauId_passed_failed["WplusJets"],  
  //					 pTauId_passed_failed["WplusJets"]->getVal(),  1.0*pTauId_passed_failed["WplusJets"]->getVal()));
  //fitConstraints.Add(makeFitConstraint(pTauId_passed_failed["TTplusJets"], 
  //					 pTauId_passed_failed["TTplusJets"]->getVal(), 1.0*pTauId_passed_failed["TTplusJets"]->getVal()));
  for ( std::vector<RooRealVar*>::iterator alpha = alphas.begin();
	alpha != alphas.end(); ++alpha ) {
    fitConstraints.Add(makeFitConstraint(*alpha, 0.5, 0.5/sysVariedByNsigma));
  }

  RooLinkedList fitOptions;
  fitOptions.Add(new RooCmdArg(RooFit::Extended()));
  fitOptions.Add(new RooCmdArg(RooFit::SumW2Error(kTRUE)));
  if ( fitConstraints.GetEntries() > 0 ) 
    fitOptions.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraints))));
  //fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(10)));
  fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(-1)));
  fitOptions.Add(new RooCmdArg(RooFit::Save(true)));

  pdfSimultaneousFit->printCompactTree();

  //RooFitResult* fitResult = pdfSimultaneousFit->fitTo(*data, fitOptions);

  RooAbsReal* nll = pdfSimultaneousFit->createNLL(*data, fitOptions); 
  RooMinuit minuit(*nll);
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
     
  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    normFactors_fitted[*process]["ALL"]    = getNormInRegion(norm, pTauId_passed_failed, *process, "ALL");
    normFactors_fitted[*process]["PASSED"] = getNormInRegion(norm, pTauId_passed_failed, *process, "PASSED");
    normFactors_fitted[*process]["FAILED"] = getNormInRegion(norm, pTauId_passed_failed, *process, "FAILED");
  }
  hasFitConverged = (fitResult->status() == 0) ? true : false;

//--- store fitted/morphed template shapes
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    TH1* template_passed = templatesAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
    TH1* templateFittedShape_passed = ( morphSysUncertaintyUp != morphSysUncertaintyDown ) ?
      compFittedTemplateShape(template_passed, pdf[*process]["PASSED"]) : template_passed;
    templatesAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed").append("_fittedShape")] = templateFittedShape_passed;
    TH1* template_failed = templatesAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed")];
    TH1* templateFittedShape_failed = ( morphSysUncertaintyUp != morphSysUncertaintyDown ) ?
      compFittedTemplateShape(template_failed, pdf[*process]["FAILED"]) : template_failed;
    templatesAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed").append("_fittedShape")] = templateFittedShape_failed;
  }

  if ( verbosity ) {
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

    std::cout << "Results of fitting variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
    for ( std::vector<std::string>::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      double numEvents             = numEventsAll[*process]["ALL"][getKey("diTauMt", tauId)];
      double fittedFraction        = fittedFractions[*process]["ALL"][getKey("diTauMt", tauId)];
      double fittedEvents          = fittedFraction*numEvents;
      double numEvents_passed      = numEventsAll[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
      double fittedFraction_passed = fittedFractions[*process]["PASSED"][getKey(fitVariable, tauId, "passed")];
      double fittedEvents_passed   = fittedFraction_passed*numEvents_passed;
      double numEvents_failed      = numEventsAll[*process]["FAILED"][getKey(fitVariable, tauId, "failed")];
      
      std::cout << " " << (*process) << ":" << std::endl;
      std::cout << "  normalization = " << norm[*process]->getVal()
		<< " +/- " << dynamic_cast<RooRealVar*>(norm[*process])->getError()
		<< " (MC exp. = " << numEvents << ")" << std::endl;
      double pTauId_passed_failedMCexp = fittedEvents_passed/fittedEvents;
      std::cout << "  pTauId_passed_failed = " << pTauId_passed_failed[*process]->getVal() 
		<< " +/- " << pTauId_passed_failed[*process]->getError() 
		<< " (MC exp. = " << pTauId_passed_failedMCexp << ")" << std::endl;
      
      std::cout << "--> ALL = " << normFactors_fitted[*process]["ALL"]
		<< " (MC exp. = " << numEvents << ")" << std::endl;
      std::cout << "--> PASSED = " << normFactors_fitted[*process]["PASSED"]
		<< " (MC exp. = " << numEvents_passed << ")" << std::endl;
      std::cout << "--> FAILED = " << normFactors_fitted[*process]["FAILED"]
		<< " (MC exp. = " << numEvents_failed << ")" << std::endl;
    }  
  }
}

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<fitTauIdEff>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("fitTauIdEff");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("fitTauIdEff") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgFitTauIdEff = cfg.getParameter<edm::ParameterSet>("fitTauIdEff");

  typedef std::vector<std::string> vstring;
  vstring regions = cfgFitTauIdEff.getParameter<vstring>("regions");
  std::string region_passed = cfgFitTauIdEff.getParameter<std::string>("passed_region");
  std::string region_failed = cfgFitTauIdEff.getParameter<std::string>("failed_region");

  vstring tauIds = cfgFitTauIdEff.getParameter<vstring>("tauIds");
  vstring fitVariables = cfgFitTauIdEff.getParameter<vstring>("fitVariables");
  
  bool takeQCDfromData = cfgFitTauIdEff.getParameter<bool>("takeQCDfromData");
  std::string regionQCDtemplateFromData_passed = 
    cfgFitTauIdEff.getParameter<std::string>("regionTakeQCDtemplateFromData_passed");
  std::string regionQCDtemplateFromData_failed = 
    cfgFitTauIdEff.getParameter<std::string>("regionTakeQCDtemplateFromData_failed");
  
  bool allowTemplateMorphing = cfgFitTauIdEff.getParameter<bool>("allowTemplateMorphing");
  double sysVariedByNsigma = cfgFitTauIdEff.getParameter<double>("sysVariedByNsigma");
  std::string morphSysUncertaintyUp = "CENTRAL_VALUE";
  std::string morphSysUncertaintyDown = "CENTRAL_VALUE";
  if ( allowTemplateMorphing ) {
    std::string morphSysUncertainty = cfgFitTauIdEff.getParameter<std::string>("morphSysUncertainty");
    morphSysUncertaintyUp = morphSysUncertainty.append("Up");
    morphSysUncertaintyDown = morphSysUncertainty.append("Down");
  }
  
  std::vector<std::string> loadSysShifts;
  loadSysShifts.push_back("CENTRAL_VALUE");
  vstring loadSysUncertaintyNames = cfgFitTauIdEff.getParameter<vstring>("loadSysUncertainties");
  for ( vstring::const_iterator sysUncertaintyName = loadSysUncertaintyNames.begin();
	sysUncertaintyName != loadSysUncertaintyNames.end(); ++sysUncertaintyName ) {
    loadSysShifts.push_back(std::string(*sysUncertaintyName).append("Up"));
    loadSysShifts.push_back(std::string(*sysUncertaintyName).append("Down"));
  }

  bool runClosureTest = cfgFitTauIdEff.getParameter<bool>("runClosureTest");
  
  bool runPseudoExperiments = cfgFitTauIdEff.getParameter<bool>("runPseudoExperiments");
  unsigned numPseudoExperiments = cfgFitTauIdEff.getParameter<unsigned>("numPseudoExperiments");
  std::vector<sysUncertaintyEntry> varySysUncertainties;
  vstring varySysUncertaintyNames = cfgFitTauIdEff.getParameter<vstring>("varySysUncertainties");
  for ( vstring::const_iterator sysUncertaintyName = varySysUncertaintyNames.begin();
	sysUncertaintyName != varySysUncertaintyNames.end(); ++sysUncertaintyName ) {
    varySysUncertainties.push_back(sysUncertaintyEntry(std::string(*sysUncertaintyName).append("Up"),
						       std::string(*sysUncertaintyName).append("Down"),
						       std::string(*sysUncertaintyName).append("Diff")));
  }

  bool makeControlPlots = cfgFitTauIdEff.getParameter<bool>("makeControlPlots");
  
  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("fitTauIdEff") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  typedef std::map<std::string, std::map<std::string, TH1*> > histogramMap; // key = (region, observable + sysShift)
  histogramMap distributionsData;

  histogramMap templatesZtautau;
  histogramMap templatesZmumu;
  histogramMap templatesQCD;
  histogramMap templatesWplusJets;
  histogramMap templatesTTplusJets;

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgFitTauIdEff.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("fitTauIdEff") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";
  
  for ( std::vector<std::string>::const_iterator sysShift = loadSysShifts.begin();
	sysShift != loadSysShifts.end(); ++sysShift ) {
    std::cout << "loading histograms for sysShift = " << (*sysShift) << "..." << std::endl;
    
    if ( !runClosureTest ) {
      loadHistograms(distributionsData, histogramInputDirectory, "Data",       regions, tauIds, fitVariables, *sysShift);
    }

    loadHistograms(templatesZtautau,    histogramInputDirectory, "Ztautau",    regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesZmumu,      histogramInputDirectory, "Zmumu",      regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesQCD,        histogramInputDirectory, "QCD",        regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesWplusJets,  histogramInputDirectory, "WplusJets",  regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesTTplusJets, histogramInputDirectory, "TTplusJets", regions, tauIds, fitVariables, *sysShift);
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

    double mcToDataScaleFactor = 1.; // run closure test for event statistics of current dataset
    //double mcToDataScaleFactor = 10.; // run closure test assuming event statistics to be higher by given factor

    sumHistograms(templatesAll, processes, "sum", mcToDataScaleFactor);

    std::cout << std::endl;
    
    distributionsData = templatesAll["sum"];
  }

//--- obtain template for QCD background from data,
//    from SS && Mt < 40 GeV && (Pzeta - 1.5 PzetaVis) > -20 GeV sideband
  if ( takeQCDfromData ) {
    for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	  tauId != tauIds.end(); ++tauId ) {
      for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	std::string key_all    = getKey(*fitVariable, *tauId);
	std::string key_passed = getKey(*fitVariable, *tauId, "passed");
	std::string key_failed = getKey(*fitVariable, *tauId, "failed");
	
	std::string keyQCD_passed = key_all;
	if ( region_passed.find("p") != std::string::npos ) keyQCD_passed = key_passed;
	std::string keyQCD_failed = key_all;
	if ( region_failed.find("f") != std::string::npos ) keyQCD_failed = key_failed;
	
	std::string keyData_passed = key_all;
	if ( regionQCDtemplateFromData_passed.find("p") != std::string::npos ) keyData_passed = key_passed;
	std::string keyData_failed = key_all;
	if ( regionQCDtemplateFromData_failed.find("f") != std::string::npos ) keyData_failed = key_failed;
	
	std::cout << "templatesQCD[" << region_passed << "][" << key_passed << "] = " 
		  << templatesQCD[region_passed][keyQCD_passed]  << std::endl;
	std::cout << "templatesQCD[" << region_failed << "][" << key_failed << "] = " 
		  << templatesQCD[region_passed][keyQCD_failed]  << std::endl;
	std::cout << "distributionsData" 
		  << "['" << regionQCDtemplateFromData_passed << "'][" << keyData_passed << "] = " 
		  << distributionsData[regionQCDtemplateFromData_passed][keyData_passed] << std::endl;
	std::cout << "distributionsData" 
		  << "['" << regionQCDtemplateFromData_failed << "'][" << keyData_failed << "] = " 
		  << distributionsData[regionQCDtemplateFromData_failed][keyData_failed] << std::endl;
	
	std::string histogramNameQCD_passed = templatesQCD[region_passed][keyQCD_passed]->GetName();
	double normQCD_passed = getIntegral(templatesQCD[region_passed][keyQCD_passed], true, true);	
	templatesQCD[region_passed][keyQCD_passed] = 
	  normalize(distributionsData[regionQCDtemplateFromData_passed][keyData_passed], normQCD_passed);
	templatesQCD[region_passed][keyQCD_passed]->SetName(histogramNameQCD_passed.data());
	
	std::string histogramNameQCD_failed = templatesQCD[region_failed][keyQCD_failed]->GetName();
	double normQCD_failed = getIntegral(templatesQCD[region_failed][keyQCD_failed], true, true);
	templatesQCD[region_passed][keyQCD_failed] = 
	  normalize(distributionsData[regionQCDtemplateFromData_failed][keyData_failed], normQCD_failed);
	templatesQCD[region_failed][keyQCD_failed]->SetName(histogramNameQCD_failed.data());
      }
    }

    templatesAll["QCD"] = templatesQCD;
  }

//--- compute sum of "passed" and "failed" regions
  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      std::string key_all    = getKey(*fitVariable, *tauId);
      std::string key_passed = key_all;
      if ( region_passed.find("p") != std::string::npos ) key_passed = getKey(*fitVariable, *tauId, "passed");
      std::string key_failed = key_all;
      if ( region_failed.find("f") != std::string::npos ) key_failed = getKey(*fitVariable, *tauId, "failed");
   
      TH1* distributionData_passed = distributionsData[region_passed][key_passed];
      TH1* distributionData_failed = distributionsData[region_failed][key_failed];
      TString distributionDataName_all = TString(distributionData_passed->GetName()).ReplaceAll(region_passed.data(), "Pass&Fail");
      TH1* distributionData_all = (TH1*)distributionData_passed->Clone(distributionDataName_all.Data());
      distributionData_all->Add(distributionData_failed);

      distributionsData["PASSED"][key_passed] = distributionData_passed;
      distributionsData["FAILED"][key_failed] = distributionData_failed;
      distributionsData["ALL"][key_all]       = distributionData_all;

      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	TH1* templateProcess_passed = templatesAll[*process][region_passed][key_passed];
	TH1* templateProcess_failed = templatesAll[*process][region_failed][key_failed];
	TString templateProcessName_all = TString(templateProcess_passed->GetName()).ReplaceAll(region_passed.data(), "Pass&Fail");
	TH1* templateProcess_all = (TH1*)templateProcess_passed->Clone();
	templateProcess_all->Add(templateProcess_failed);

	templatesAll[*process]["PASSED"][key_passed] = templateProcess_passed;
	templatesAll[*process]["FAILED"][key_failed] = templateProcess_failed;
	templatesAll[*process]["ALL"][key_all]       = templateProcess_all;
      }
    }
  }

  vstring regionsToDraw;
  regionsToDraw.push_back("PASSED");
  regionsToDraw.push_back("FAILED");
  regionsToDraw.push_back("ALL");

//--- define x-axis titles
  std::map<std::string, std::string> xAxisTitles;
  xAxisTitles["diTauCharge"]         = "Charge(#mu + #tau_{had})";
  xAxisTitles["diTauMt"]             = "M_{T}^{#muMET} / GeV";
  xAxisTitles["diTauHt"]             = "P_{T}^{#mu} + P_{T}^{#tau} + MET / GeV";
  xAxisTitles["diTauVisMass"]        = "M_{vis}^{#mu#tau} / GeV";
  xAxisTitles["diTauVisMassFromJet"] = xAxisTitles["diTauVisMass"];

//--- make control plots plots of Data compared to sum(MC) scaled by cross-sections
//    for muonPt, tauPt, Mt, visMass,... distributions in different regions
  if ( makeControlPlots ) {
    for ( vstring::const_iterator region = regionsToDraw.begin();
	  region != regionsToDraw.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = distributionsData[*region].begin();
	    key != distributionsData[*region].end(); ++key ) {
	if ( !isSystematicShift(key->first) ) {
	  //std::string histogramTitle = std::string("Region ").append(*region).append(": ").append(key->first);
	  //histogramTitle.append(" (scaled by cross-section)");
	  std::string histogramTitle = "";
	  std::string outputFileName = std::string("controlPlotsTauIdEff_");
	  outputFileName.append(*region).append("_").append(key->first).append(".png");
	  drawHistograms(templatesZtautau[*region][key->first], -1.,
			 templatesZmumu[*region][key->first], -1.,
			 templatesQCD[*region][key->first], -1.,
			 templatesWplusJets[*region][key->first], -1.,
			 templatesTTplusJets[*region][key->first], -1.,
			 distributionsData[*region][key->first],
			 histogramTitle, xAxisTitles[getObservable(key->first)],
			 outputFileName);
	}
      }
    }
  }
  
  for ( vstring::const_iterator region = regionsToDraw.begin();
	region != regionsToDraw.end(); ++region ) {
    for ( std::map<std::string, TH1*>::const_iterator key = distributionsData[*region].begin();
	  key != distributionsData[*region].end(); ++key ) {
      std::cout << "numEvents[" << "Data" << "][" << (*region) << "][" << key->first << "] = "
		<< getIntegral(key->second, true, true) << std::endl;
    }
  }

  std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll = // key = (process/"sum", region, observable)
    compNumEvents(templatesAll, processes, distributionsData);
  std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions = // key = (process, region, observable)
    compFittedFractions(templatesAll, numEventsAll, processes, distributionsData);

  std::cout << "running fit for central values..." << std::endl;
  
  std::map<std::string, std::map<std::string, double> > effValues;           // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > effErrors;           // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, double> > > > normFactorsAll_fitted; // key = (process/"sum", region, tauId, fitVariable)
  std::map<std::string, std::map<std::string, bool> >   fitConvergenceFlags; // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > tauIdEffMCexp;       // key = (tauId, fitVariable)

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      double effValue = 0.;
      double effError = 1.;
      std::map<std::string, std::map<std::string, double> > normFactors_fitted;
      bool hasFitConverged = false;
      fitUsingRooFit(distributionsData, templatesAll, numEventsAll, fittedFractions,
		     processes,
		     *tauId, *fitVariable, 
		     morphSysUncertaintyUp, morphSysUncertaintyDown, sysVariedByNsigma,
		     effValue, effError, normFactors_fitted, hasFitConverged,	       
		     1);
      effValues[*tauId][*fitVariable] = effValue;
      effErrors[*tauId][*fitVariable] = effError;
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( vstring::const_iterator region = regionsToDraw.begin();
	      region != regionsToDraw.end(); ++region ) {
	  normFactorsAll_fitted[*process][*region][*tauId][*fitVariable] = normFactors_fitted[*process][*region];
	}
      }
      fitConvergenceFlags[*tauId][*fitVariable] = hasFitConverged;
    }
  }

//--- make control plots of Data compared to sum(MC) scaled by normalization factors determined by fit
//    for muonPt, tauPt, Mt, visMass,... distributions in different regions
  if ( makeControlPlots ) {
    for ( vstring::const_iterator region = regionsToDraw.begin();
	  region != regionsToDraw.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = distributionsData[*region].begin();
	    key != distributionsData[*region].end(); ++key ) {
	for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	      tauId != tauIds.end(); ++tauId ) {
	  if ( !(key->first.find(*tauId) != std::string::npos) ) continue;
	  for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
		fitVariable != fitVariables.end(); ++fitVariable ) {
	    if ( !isSystematicShift(key->first) ) {
	      //std::string histogramTitle = std::string("Region ").append(*region).append(": ").append(key->first);
	      //histogramTitle.append(" (scaled by normalization det. by fit)");
	      std::string histogramTitle = "";
	      std::string outputFileName = std::string("controlPlotsTauIdEff_").append(*region).append("_").append(key->first);
	      outputFileName.append("_fitted_").append(*fitVariable).append(".png");
	      std::string key_mc = std::string(key->first).append("_fittedShape");
	      drawHistograms(
	        templatesZtautau[*region][key_mc], 
		getTemplateNorm_fitted("Ztautau", *region, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesZmumu[*region][key_mc], 
		getTemplateNorm_fitted("Zmumu", *region, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesQCD[*region][key_mc], 
		getTemplateNorm_fitted("QCD", *region, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesWplusJets[*region][key_mc], 
		getTemplateNorm_fitted("WplusJets", *region, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesTTplusJets[*region][key_mc], 
		getTemplateNorm_fitted("TTplusJets", *region, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		distributionsData[*region][key->first],
		histogramTitle, xAxisTitles[getObservable(key->first)],
		outputFileName, true);
	    }
	  }
	}
      }
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
      
      double numEvents        = numEventsAll["Ztautau"]["ALL"][getKey("diTauMt", *tauId)];
      double numEvents_passed = numEventsAll["Ztautau"]["PASSED"][getKey(fitVariables.front(), *tauId, "passed")];
      
      tauIdEffMCexp[*tauId][*fitVariable] = numEvents_passed/numEvents;
      std::cout << "(Monte Carlo prediction = " << tauIdEffMCexp[*tauId][*fitVariable]*100. << "%)" << std::endl;
    }
  }

  if ( runPseudoExperiments ) {
    
//--- compute systematic shifts 
    for ( std::vector<sysUncertaintyEntry>::const_iterator sysUncertainty = varySysUncertainties.begin();
	  sysUncertainty != varySysUncertainties.end(); ++sysUncertainty ) {
      compSysHistograms(templatesAll, *sysUncertainty);
    }

    std::map<std::string, std::map<std::string, TH1*> > effDistributions;            // key = (tauId, fitVariable)
    std::map<std::string, std::map<std::string, TH1*> > fitConvergenceDistributions; // key = (tauId, fitVariable)
    
    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >   templatesAll_fluctuated;
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll_fluctuated;
    std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions_fluctuated;

    for ( unsigned iPseudoExperiment = 0; iPseudoExperiment < numPseudoExperiments; ++numPseudoExperiments ) {
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( vstring::const_iterator region = regionsToDraw.begin();
	      region != regionsToDraw.end(); ++region ) {
	  for ( std::map<std::string, TH1*>::const_iterator key = distributionsData[*region].begin();
		key != distributionsData[*region].end(); ++key ) {
	    if ( !isSystematicShift(key->first) ) {
	      TH1* origHistogram = key->second;
	      
	      TH1* fluctHistogram = templatesAll_fluctuated[*process][*region][key->first];
	      if ( !fluctHistogram ) {
		fluctHistogram = (TH1*)origHistogram->Clone(TString(origHistogram->GetName()).Append("_fluctuated"));
		templatesAll_fluctuated[*process][*region][key->first] = fluctHistogram;
	      }
	      
	      sampleHistogram_stat(origHistogram, fluctHistogram);
	      
	      for ( std::vector<sysUncertaintyEntry>::const_iterator sysUncertainty = varySysUncertainties.begin();
		    sysUncertainty != varySysUncertainties.end(); ++sysUncertainty ) {
		std::string key_sys = std::string(key->first).append("_").append(sysUncertainty->sysNameDiff_);
		TH1* sysHistogram = templatesAll_fluctuated[*process][*region][key_sys];
		assert(sysHistogram);
		
		sampleHistogram_sys(fluctHistogram, sysHistogram, 1.0/sysVariedByNsigma, -1.0, +1.0, kCoherent);
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
	  std::map<std::string, std::map<std::string, double> > normFactors_fitted;
	  bool hasFitConverged = false;
	  fitUsingRooFit(distributionsData, templatesAll_fluctuated, numEventsAll_fluctuated, fittedFractions_fluctuated,
			 processes,
			 *tauId, *fitVariable, 
			 morphSysUncertaintyUp, morphSysUncertaintyDown, sysVariedByNsigma,
			 effValue, effError, normFactors_fitted, hasFitConverged,	       
			 0);

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

//-- save fit results
  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory fitResultOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      std::string fitResultName  = std::string("fitResult_").append(*fitVariable).append("_").append(*tauId);
      std::string fitResultTitle = std::string(*tauId).append(" Efficiency, obtained by fitting ").append(*fitVariable);
      TH1* histogramFitResult = fitResultOutputDirectory.make<TH1F>(fitResultName.data(), fitResultTitle.data(), 1, -0.5, +0.5);
      int fitResultBin = histogramFitResult->FindBin(0.);
      histogramFitResult->SetBinContent(fitResultBin, effValues[*tauId][*fitVariable]);
      histogramFitResult->SetBinError(fitResultBin, effErrors[*tauId][*fitVariable]);

      std::string fitNormName  = std::string("fitNorm_").append(*fitVariable).append("_").append(*tauId);
      std::string fitNormTitle = std::string("Fitted Number of Z #rightarrow #tau^{+} #tau^{-} Events in 'passed' + 'failed' regions");
      TH1* histogramFitNorm = fitResultOutputDirectory.make<TH1F>(fitNormName.data(), fitNormTitle.data(), 1, -0.5, +0.5);
      int fitNormBin = histogramFitNorm->FindBin(0.);
      histogramFitNorm->SetBinContent(fitNormBin, normFactorsAll_fitted["Ztautau"]["ALL"][*tauId][*fitVariable]);

      std::string expResultName  = std::string("expResult_").append(*fitVariable).append("_").append(*tauId);
      std::string expResultTitle = std::string("Expected ").append(*tauId).append(" Efficiency");
      TH1* histogramExpResult = fitResultOutputDirectory.make<TH1F>(expResultName.data(), expResultTitle.data(), 1, -0.5, +0.5);
      int expResultBin = histogramExpResult->FindBin(0.);
      histogramExpResult->SetBinContent(expResultBin, tauIdEffMCexp[*tauId][*fitVariable]);

      std::string expNormName  = std::string("expNorm_").append(*fitVariable).append("_").append(*tauId);
      std::string expNormTitle = std::string("Expected Number of Z #rightarrow #tau^{+} #tau^{-} Events in 'passed' + 'failed' regions");
      TH1* histogramExpNorm = fitResultOutputDirectory.make<TH1F>(expNormName.data(), expNormTitle.data(), 1, -0.5, +0.5);
      int expNormBin = histogramExpNorm->FindBin(0.);
      histogramExpNorm->SetBinContent(expNormBin, numEventsAll["Ztautau"]["ALL"][getKey("diTauMt", *tauId)]);
    }
  }

//--print time that it took macro to run
  std::cout << "finished executing fitTauIdEff macro:" << std::endl;
  std::cout << " #tauIdDiscriminators = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables        = " << fitVariables.size() << std::endl;
  if ( runPseudoExperiments ) {
    std::cout << " #sysUncertainties    = " << varySysUncertainties.size() 
	      << " (numPseudoExperiments = " << numPseudoExperiments << ")" << std::endl;
  }
  clock.Show("fitTauIdEff");

  return 0;
}
  
 
