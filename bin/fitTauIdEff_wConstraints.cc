
/** \executable fitTauIdEff_wConstraints
 *
 * Perform simultaneous fit of Ztautau signal plus Zmumu, W + jets, QCD and TTbar background yields
 * in regions A, B, C1f, C1p, C2 and D, in order to determine tau identification efficiency
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.16 $
 *
 * $Id: fitTauIdEff_wConstraints.cc,v 1.16 2011/07/13 16:54:31 veelken Exp $
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

RooFormulaVar* makeRooFormulaVar(const std::string& process, const std::string& region,
				 RooAbsReal* norm, double fittedFractionValue, unsigned numCategories,
				 RooAbsReal* pDiTauCharge_OS_SS, RooAbsReal* pDiTauKine_Sig_Bgr, 
				 RooAbsReal* pMuonIso_tight_loose, 
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

  //std::cout << "<makeRooFormulaVar>:" << std::endl;
  //std::cout << " building RooFormulaVar expression for process = " << process << ", region = " << region << "." << std::endl;

  RooFormulaVar* retVal = 0;

  std::string exprDiTauCharge = "";
  std::string exprDiTauKine   = "";
  std::string exprMuonIso     = "";
  std::string exprTauId       = "";

  if        ( region.find("A") != std::string::npos ) {
    exprDiTauCharge = "regular";
    exprMuonIso     = "inverted";
  } else if ( region.find("B") != std::string::npos ) {
    exprDiTauCharge = "inverted";
    exprMuonIso     = "inverted";
  } else if ( region.find("C") != std::string::npos ) {
    exprDiTauCharge = "regular";
    exprMuonIso     = "regular"; 
  } else if ( region.find("D") != std::string::npos ) {
    exprDiTauCharge = "inverted";
    exprMuonIso     = "regular";
  } else {
    throw cms::Exception("makeRooFormulaVar") 
      << "Region = " << region << " not defined !!\n";
  }
  
  if      ( region.find("1") != std::string::npos ) exprDiTauKine = "regular";
  else if ( region.find("2") != std::string::npos ) exprDiTauKine = "inverted";
  
  if      ( region.find("p") != std::string::npos ) exprTauId     = "regular";
  else if ( region.find("f") != std::string::npos ) exprTauId     = "inverted";

  std::string fittedFractionName = std::string(norm->GetName()).append("_").append(region).append("_fittedFraction");
  //RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), fittedFractionValue);
  // CV: fitted yields for all processes come-out factor N too high in closure test
  //    --> compensate by multiplying fittedFraction by fudge-factor N,
  //        N = number of categories used in simultaneous fit
  RooConstVar* fittedFraction = new RooConstVar(fittedFractionName.data(), fittedFractionName.data(), numCategories*fittedFractionValue);

  std::string formula = "";
  TObjArray arguments; 
  addToFormula(formula, exprDiTauCharge, arguments, pDiTauCharge_OS_SS);
  addToFormula(formula, exprDiTauKine,   arguments, pDiTauKine_Sig_Bgr);
  addToFormula(formula, exprMuonIso,     arguments, pMuonIso_tight_loose);
  addToFormula(formula, exprTauId,       arguments, pTauId_passed_failed);
  addToFormula(formula, "regular",       arguments, norm);
  addToFormula(formula, "regular",       arguments, fittedFraction);

  //std::cout << " formula = " << formula << std::endl;
  //std::cout << " arguments:" << std::endl;
  //for ( int i = 0; i < arguments.GetEntries(); ++i ) {
  //  std::cout << "  " << dynamic_cast<TNamed*>(arguments.At(i))->GetName() << ":" 
  //	        << " value = " << dynamic_cast<RooAbsReal*>(arguments.At(i))->getVal() << std::endl;
  //}

  std::string name = std::string("norm").append(process).append("_").append(region);
  retVal = new RooFormulaVar(name.data(), name.data(), formula.data(), RooArgSet(arguments));

  return retVal;
}

double getNormInRegion(std::map<std::string, RooRealVar*> normABCD, 
		       std::map<std::string, RooRealVar*> pDiTauCharge_OS_SS, 
		       std::map<std::string, RooRealVar*> pMuonIso_tight_loose,
		       std::map<std::string, RooRealVar*> pDiTauKine_Sig_Bgr, 
		       std::map<std::string, RooRealVar*> pTauId_passed_failed,
		       const std::string& process, const std::string& region)
{
  double retVal = normABCD[process]->getVal();

  if ( region.find("ABCD") == std::string::npos) {
    if        ( region.find("A") != std::string::npos ) {
      retVal *= pDiTauCharge_OS_SS[process]->getVal();
      retVal *= (1. - pMuonIso_tight_loose[process]->getVal());
    } else if ( region.find("B") != std::string::npos ) {
      retVal *= (1. - pDiTauCharge_OS_SS[process]->getVal());
      retVal *= (1. - pMuonIso_tight_loose[process]->getVal());
    } else if ( region.find("C") != std::string::npos ) {
      retVal *= pDiTauCharge_OS_SS[process]->getVal();
      retVal *= pMuonIso_tight_loose[process]->getVal();
    } else if ( region.find("D") != std::string::npos ) {
      retVal *= (1. - pDiTauCharge_OS_SS[process]->getVal());
      retVal *= pMuonIso_tight_loose[process]->getVal();
    } else {
      throw cms::Exception("getNormInRegion") 
	<< "Region = " << region << " not defined !!\n";
    }
    
    if      ( region.find("1") != std::string::npos ) retVal *= pDiTauKine_Sig_Bgr[process]->getVal();
    else if ( region.find("2") != std::string::npos ) retVal *= (1. - pDiTauKine_Sig_Bgr[process]->getVal());
    
    if      ( region.find("p") != std::string::npos ) retVal *= pTauId_passed_failed[process]->getVal();
    else if ( region.find("f") != std::string::npos ) retVal *= (1. - pTauId_passed_failed[process]->getVal());
  }

  return retVal;
}

void fitUsingRooFit(std::map<std::string, std::map<std::string, TH1*> >& distributionsData,                         // key = (region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, TH1*> > >& templatesAll,      // key = (process, region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& numEventsAll,    // key = (process/"sum", region, observable)
		    std::map<std::string, std::map<std::string, std::map<std::string, double> > >& fittedFractions, // key = (process, region, observable)
		    const std::vector<std::string>& processes,
		    const std::string& tauId, const std::string& fitVariable, bool fitTauIdEffC2,
		    double& effValue, double& effError, 
		    std::map<std::string, std::map<std::string, double> >& normFactors_fitted, bool& hasFitConverged,
		    int verbosity = 0)
{
  if ( verbosity ) {
    std::cout << "<fitUsingRooFit>:" << std::endl;
    std::cout << " performing Fit of variable = " << fitVariable << " for Tau id. = " << tauId << std::endl;
  }

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
  //std::cout << "numEventsDataABCD = " << numEventsDataABCD << std::endl;

  std::map<std::string, RooRealVar*> pDiTauCharge_OS_SS;          // key = process
  std::map<std::string, RooRealVar*> pMuonIso_tight_loose;        // key = process
  std::map<std::string, RooRealVar*> pDiTauKine_Sig_Bgr;          // key = process
  std::map<std::string, RooRealVar*> pTauId_passed_failed;        // key = process

  std::map<std::string, RooRealVar*> normABCD;                    // key = process
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
  
  unsigned numCategoriesABC2D = ( fitTauIdEffC2 ) ? 5 : 4;
  unsigned numCategoriesC1    = 2;

  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    //std::cout << "process = " << (*process) << ":" << std::endl;

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
    double fittedFractionC1f  = fittedFractions[*process]["C1f"][getKey(fitVariable, tauId, "failed")];
    double numEventsC2p       = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedFractionC2p  = fittedFractions[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
    double fittedEventsC2p    = fittedFractionC2p*numEventsC2p;
    double fittedFractionC2f  = fittedFractions[*process]["C2f"][getKey(fitVariable, tauId, "failed")];
    double fittedFractionC2   = fittedFractions[*process]["C2"][getKey("diTauMt", tauId)];
    double numEventsD         = numEventsAll[*process]["D"][getKey("diTauMt", tauId)];
    double fittedFractionD    = fittedFractions[*process]["D"][getKey("diTauMt", tauId)];
    double fittedEventsD      = fittedFractionD*numEventsD;
    double numEventsABCD      = numEventsAll[*process]["ABCD"][getKey("diTauMt", tauId)];
    double fittedEventsABCD   = fittedEventsA + fittedEventsB + fittedEventsC + fittedEventsD;
    //double fittedFractionABCD = fittedEventsABCD/numEventsABCD;
    //std::cout << " numEventsABCD = " << numEventsABCD << ", fittedFractionABCD = " << fittedFractionABCD << std::endl;

    std::string nameDiTauCharge_OS_SS   = std::string("pDiTauCharge_OS_SS").append("_").append(*process);
    double pDiTauCharge_OS_SS0          = (fittedEventsA + fittedEventsC)/fittedEventsABCD;
    pDiTauCharge_OS_SS[*process]        = new RooRealVar(nameDiTauCharge_OS_SS.data(), 							 
							 nameDiTauCharge_OS_SS.data(), pDiTauCharge_OS_SS0, 0., 1.);
    std::string nameMuonIso_tight_loose = std::string("pMuonIso_tight_loose").append("_").append(*process);
    double pMuonIso_tight_loose0        = (fittedEventsC + fittedEventsD)/fittedEventsABCD;
    pMuonIso_tight_loose[*process]      = new RooRealVar(nameMuonIso_tight_loose.data(), 
							 nameDiTauCharge_OS_SS.data(), pMuonIso_tight_loose0, 0., 1.);
    std::string nameDiTauKine_Sig_Bgr   = std::string("pDiTauKine_Sig_Bgr").append("_").append(*process);
    double pDiTauKine_Sig_Bgr0          = fittedEventsC1/fittedEventsC;
    pDiTauKine_Sig_Bgr[*process]        = new RooRealVar(nameDiTauKine_Sig_Bgr.data(), 
							 nameDiTauKine_Sig_Bgr.data(), pDiTauKine_Sig_Bgr0, 0., 1.);
    std::string nameTauId_passed_failed = std::string("pTauId_passed_failed").append("_").append(*process);
    double pTauId_passed_failed0        = ( fitTauIdEffC2 ) ? 
      (fittedEventsC1p + fittedEventsC2p)/fittedEventsC : fittedEventsC1p/fittedEventsC1;
    pTauId_passed_failed[*process]      = new RooRealVar(nameTauId_passed_failed.data(),
							 nameTauId_passed_failed.data(), pTauId_passed_failed0, 0., 1.);

    //double numEventsSumABCD = numEventsAll["sum"]["ABCD"][getKey("diTauMt", tauId)];
    //std::cout << " numEventsSumABCD = " << numEventsSumABCD << std::endl;

    //double scaleFactorMCtoData = numEventsDataABCD/numEventsSumABCD;
    double scaleFactorMCtoData = 1.;
    //std::cout << "--> MC-to-Data scale-factor = " << scaleFactorMCtoData << std::endl;

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

    normA[*process]   = makeRooFormulaVar(*process, "A", normABCD[*process], fittedFractionA, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normB[*process]   = makeRooFormulaVar(*process, "B", normABCD[*process], fittedFractionB, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC1[*process]  = makeRooFormulaVar(*process, "C1", normABCD[*process], fittedFractionC1, numCategoriesC1, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC1p[*process] = makeRooFormulaVar(*process, "C1p", normABCD[*process], fittedFractionC1p, numCategoriesC1,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC1f[*process] = makeRooFormulaVar(*process, "C1f", normABCD[*process], fittedFractionC1f, numCategoriesC1,
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC2[*process]  = makeRooFormulaVar(*process, "C2", normABCD[*process], fittedFractionC2, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC2p[*process] = makeRooFormulaVar(*process, "C2p", normABCD[*process], fittedFractionC2p, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);
    normC2f[*process] = makeRooFormulaVar(*process, "C2f", normABCD[*process], fittedFractionC2f, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);    
    normD[*process]   = makeRooFormulaVar(*process, "D", normABCD[*process], fittedFractionD, numCategoriesABC2D, 
					  pDiTauCharge_OS_SS[*process], pDiTauKine_Sig_Bgr[*process], pMuonIso_tight_loose[*process], 
					  pTauId_passed_failed[*process]);

    fitParametersA.Add(normA[*process]);
    fitParametersB.Add(normB[*process]);
    fitParametersC1p.Add(normC1p[*process]);
    fitParametersC1f.Add(normC1f[*process]);
    fitParametersC2p.Add(normC2p[*process]);
    fitParametersC2f.Add(normC2f[*process]);
    fitParametersC2.Add(normC2[*process]);
    fitParametersD.Add(normD[*process]);
  }

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

  RooSimultaneous* pdfSimultaneousFitABC2D = 
    new RooSimultaneous("pdfSimultaneousFitABC2D", "pdfSimultaneousFitABC2D", *fitCategoriesABC2D);
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

  RooSimultaneous* pdfSimultaneousFitC1 = 
    new RooSimultaneous("pdfSimultaneousFitC1", "pdfSimultaneousFitC1", *fitCategoriesC1);
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
//   o pMuonIso_tight_loose
//    separating different regions

  pDiTauCharge_OS_SS["Ztautau"]->setConstant(true);
  pMuonIso_tight_loose["Ztautau"]->setConstant(true);
  pDiTauKine_Sig_Bgr["Ztautau"]->setConstant(true);

  //pDiTauCharge_OS_SS["Zmumu"]->setConstant(true);
  pMuonIso_tight_loose["Zmumu"]->setConstant(true);
  pDiTauKine_Sig_Bgr["Zmumu"]->setConstant(true);

  pDiTauCharge_OS_SS["TTplusJets"]->setConstant(true);
  pMuonIso_tight_loose["TTplusJets"]->setConstant(true);
  pDiTauKine_Sig_Bgr["TTplusJets"]->setConstant(true);

//--- set tau id. efficiency to "random" value
  pTauId_passed_failed["Ztautau"]->setVal(0.55);

  TObjArray fitConstraintsC1;
  //fitConstraintsC1.Add(makeFitConstraint(normABCD["Zmumu"],      
  //					   normABCD["Zmumu"]->getVal(),                 1.0*normABCD["Zmumu"]->getVal()));
  //normABCD["Zmumu"]->setMax(5.0*normABCD["Zmumu"]->getVal());
  //fitConstraintsC1.Add(makeFitConstraint(normABCD["QCD"],        
  //					   normABCD["QCD"]->getVal(),                   1.0*normABCD["QCD"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(normABCD["WplusJets"],  
  //					   normABCD["WplusJets"]->getVal(),             1.0*normABCD["WplusJets"]->getVal()));
  //fitConstraintsC1.Add(makeFitConstraint(normABCD["TTplusJets"], 
  //					   normABCD["TTplusJets"]->getVal(),            1.0*normABCD["TTplusJets"]->getVal()));
  //normABCD["TTplusJets"]->setMax(5.0*normABCD["TTplusJets"]->getVal());
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

  //pdfSimultaneousFitC1->printCompactTree();

  TObjArray fitConstraintsABC2D;
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["QCD"],         
  //					      pDiTauCharge_OS_SS["QCD"]->getVal(),         0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["QCD"],         
  //					      pDiTauKine_Sig_Bgr["QCD"]->getVal(),         0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_tight_loose["QCD"],       
  //					      pMuonIso_tight_loose["QCD"]->getVal(),       0.25));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauCharge_OS_SS["WplusJets"],   
  //					      pDiTauCharge_OS_SS["WplusJets"]->getVal(),   0.1));
  //fitConstraintsABC2D.Add(makeFitConstraint(pDiTauKine_Sig_Bgr["WplusJets"],   
  //					      pDiTauKine_Sig_Bgr["WplusJets"]->getVal(),   0.2));
  //fitConstraintsABC2D.Add(makeFitConstraint(pMuonIso_tight_loose["WplusJets"], 
  //					      pMuonIso_tight_loose["WplusJets"]->getVal(), 0.1));

  RooLinkedList fitOptionsABC2D;
  fitOptionsABC2D.Add(new RooCmdArg(RooFit::Extended()));
  //fitOptionsABC2D.Add(new RooCmdArg(RooFit::ExternalConstraints(RooArgSet(fitConstraintsABC2D))));

  //pdfSimultaneousFitABC2D->printCompactTree();

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
   
  effValue = pTauId_passed_failed["Ztautau"]->getVal();
  effError = pTauId_passed_failed["Ztautau"]->getError();
  for ( std::vector<std::string>::const_iterator process = processes.begin();
	process != processes.end(); ++process ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      normFactors_fitted[*process][region->first] = 
	getNormInRegion(normABCD, 
			pDiTauCharge_OS_SS, pMuonIso_tight_loose, pDiTauKine_Sig_Bgr, pTauId_passed_failed, 
			*process, region->first);
    }
  }
  hasFitConverged = (fitResult->status() == 0) ? true : false;
  
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
      double numEventsA         = numEventsAll[*process]["A"][getKey("diTauMt", tauId)];  
      double fittedFractionA    = fittedFractions[*process]["A"][getKey("diTauMt", tauId)];
      double fittedEventsA      = fittedFractionA*numEventsA;
      double numEventsB         = numEventsAll[*process]["B"][getKey("diTauMt", tauId)];    
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
      double numEventsC2p       = numEventsAll[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
      double fittedFractionC2p  = fittedFractions[*process]["C2p"][getKey("diTauMt", tauId, "passed")];
      double fittedEventsC2p    = fittedFractionC2p*numEventsC2p;
      double numEventsC2f       = numEventsAll[*process]["C2f"][getKey("diTauMt", tauId, "failed")];
      double numEventsC2        = numEventsAll[*process]["C2"][getKey("diTauMt", tauId)];
      double numEventsD         = numEventsAll[*process]["D"][getKey("diTauMt", tauId)];
      double fittedFractionD    = fittedFractions[*process]["D"][getKey("diTauMt", tauId)];
      double fittedEventsD      = fittedFractionD*numEventsD;
      double numEventsABCD      = numEventsAll[*process]["ABCD"][getKey("diTauMt", tauId)];
      double fittedFractionABCD = fittedFractions[*process]["ABCD"][getKey("diTauMt", tauId)];
      double fittedEventsABCD   = fittedFractionABCD*numEventsABCD;

      std::cout << " " << (*process) << ":" << std::endl;
      std::cout << "  normalization = " << normABCD[*process]->getVal()
		<< " +/- " << dynamic_cast<RooRealVar*>(normABCD[*process])->getError()
		<< " (MC exp. = " << numEventsABCD << ")" << std::endl;
      std::cout << "  pDiTauCharge_OS_SS = " << pDiTauCharge_OS_SS[*process]->getVal() 
		<< " +/- " << pDiTauCharge_OS_SS[*process]->getError() 
		<< " (MC exp. = " << (fittedEventsA + fittedEventsC)/fittedEventsABCD << ")" << std::endl;
      std::cout << "  pMuonIso_tight_loose = " << pMuonIso_tight_loose[*process]->getVal() 
		<< " +/- " << pMuonIso_tight_loose[*process]->getError() 
		<< " (MC exp. = " << (fittedEventsC + fittedEventsD)/fittedEventsABCD << ")" << std::endl;
      std::cout << "  pDiTauKine_Sig_Bgr = " << pDiTauKine_Sig_Bgr[*process]->getVal() 
		<< " +/- " << pDiTauKine_Sig_Bgr[*process]->getError() 
		<< " (MC exp. = " << fittedEventsC1/fittedEventsC << ")" << std::endl;
      double pTauId_passed_failedMCexp = ( fitTauIdEffC2 ) ? 
	(fittedEventsC1p + fittedEventsC2p)/fittedEventsC : fittedEventsC1p/fittedEventsC1;
      std::cout << "  pTauId_passed_failed = " << pTauId_passed_failed[*process]->getVal() 
		<< " +/- " << pTauId_passed_failed[*process]->getError() 
		<< " (MC exp. = " << pTauId_passed_failedMCexp << ")" << std::endl;
      
      std::cout << "--> A = " << normFactors_fitted[*process]["A"]
		<< " (MC exp. = " << numEventsA << ")" << std::endl;
      std::cout << "--> B = " << normFactors_fitted[*process]["B"]
		<< " (MC exp. = " << numEventsB << ")" << std::endl;
      std::cout << "--> C = " << normFactors_fitted[*process]["C"]
		<< " (MC exp. = " << numEventsC << ")" << std::endl;
      std::cout << "--> C1 = " << normFactors_fitted[*process]["C1"]
		<< " (MC exp. = " << numEventsC1 << ")" << std::endl;
      std::cout << "--> C1p = " << normFactors_fitted[*process]["C1p"]
		<< " (MC exp. = " << numEventsC1p << ")" << std::endl;
      std::cout << "--> C1f = " << normFactors_fitted[*process]["C1f"]
		<< " (MC exp. = " << numEventsC1f << ")" << std::endl;
      std::cout << "--> C2 = " << normFactors_fitted[*process]["C2"]
		<< " (MC exp. = " << numEventsC2 << ")" << std::endl;
      if ( fitTauIdEffC2 ) {
	std::cout << "--> C2p = " << normFactors_fitted[*process]["C2p"]
		  << " (MC exp. = " << numEventsC2p << ")" << std::endl;
	std::cout << "--> C2f = " << normFactors_fitted[*process]["C2f"]
		  << " (MC exp. = " << numEventsC2f << ")" << std::endl;
      }
      std::cout << "--> D = " << normFactors_fitted[*process]["D"]
		<< " (MC exp. = " << numEventsD << ")" << std::endl;
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

  std::cout << "<fitTauIdEff_wConstraints>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("fitTauIdEff_wConstraints");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("fitTauIdEff_wConstraints") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgFitTauIdEff_wConstraints = cfg.getParameter<edm::ParameterSet>("fitTauIdEff_wConstraints");

  typedef std::vector<std::string> vstring;
  vstring regions = cfgFitTauIdEff_wConstraints.getParameter<vstring>("regions");
  std::string regionQCDtemplateFromData_passed = 
    cfgFitTauIdEff_wConstraints.getParameter<std::string>("regionQCDtemplateFromData_passed");
  std::string regionQCDtemplateFromData_failed = 
    cfgFitTauIdEff_wConstraints.getParameter<std::string>("regionQCDtemplateFromData_failed");
  vstring tauIds = cfgFitTauIdEff_wConstraints.getParameter<vstring>("tauIds");
  vstring fitVariables = cfgFitTauIdEff_wConstraints.getParameter<vstring>("fitVariables");
  std::vector<sysUncertaintyEntry> sysUncertainties;
  vstring sysUncertaintyNames = cfgFitTauIdEff_wConstraints.getParameter<vstring>("sysUncertainties");
  for ( vstring::const_iterator sysUncertaintyName = sysUncertaintyNames.begin();
	sysUncertaintyName != sysUncertaintyNames.end(); ++sysUncertaintyName ) {
    sysUncertainties.push_back(sysUncertaintyEntry(std::string(*sysUncertaintyName).append("Up"),
						   std::string(*sysUncertaintyName).append("Down"),
						   std::string(*sysUncertaintyName).append("Diff")));
  }
  bool runClosureTest = cfgFitTauIdEff_wConstraints.getParameter<bool>("runClosureTest");
  bool takeQCDfromData = cfgFitTauIdEff_wConstraints.getParameter<bool>("takeQCDfromData");
  bool fitTauIdEffC2 = cfgFitTauIdEff_wConstraints.getParameter<bool>("fitTauIdEffC2");
  bool runSysUncertainties = cfgFitTauIdEff_wConstraints.getParameter<bool>("runSysUncertainties");
  unsigned numPseudoExperiments = cfgFitTauIdEff_wConstraints.getParameter<unsigned>("numPseudoExperiments");
  bool makeControlPlots = cfgFitTauIdEff_wConstraints.getParameter<bool>("makeControlPlots");

  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("fitTauIdEff_wConstraints") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

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
  std::string directory = cfgFitTauIdEff_wConstraints.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("fitTauIdEff_wConstraints") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";
  
  for ( std::vector<std::string>::const_iterator sysShift = sysShifts.begin();
	sysShift != sysShifts.end(); ++sysShift ) {
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
	
	std::string keyQCDtemplateFromData_passed = key_all;
	if ( regionQCDtemplateFromData_passed.find("p") != std::string::npos ) keyQCDtemplateFromData_passed = key_passed;
	std::string keyQCDtemplateFromData_failed = key_all;
	if ( regionQCDtemplateFromData_failed.find("f") != std::string::npos ) keyQCDtemplateFromData_failed = key_failed;

	std::cout << "templatesQCD['C1']["  << key_all    << "] = " << templatesQCD["C1"][key_all]      << std::endl;
	std::cout << "templatesQCD['C1p'][" << key_passed << "] = " << templatesQCD["C1p"][key_passed]  << std::endl;
	std::cout << "templatesQCD['C1f'][" << key_failed << "] = " << templatesQCD["C1f"][key_failed]  << std::endl;
	std::cout << "distributionsData" 
		  << "['" << regionQCDtemplateFromData_passed << "'][" << keyQCDtemplateFromData_passed << "] = " 
		  << distributionsData[regionQCDtemplateFromData_passed][keyQCDtemplateFromData_passed] << std::endl;
	std::cout << "distributionsData" 
		  << "['" << regionQCDtemplateFromData_failed << "'][" << keyQCDtemplateFromData_failed << "] = " 
		  << distributionsData[regionQCDtemplateFromData_failed][keyQCDtemplateFromData_failed] << std::endl;
	
	std::string histogramNameQCD_C1 = templatesQCD["C1"][key_all]->GetName();
	double normQCD_C1 = getIntegral(templatesQCD["C1"][key_all], true, true);
	templatesQCD["C1"][key_all] = 
	  normalize(distributionsData[regionQCDtemplateFromData_failed][keyQCDtemplateFromData_failed], normQCD_C1);
	templatesQCD["C1"][key_all]->SetName(histogramNameQCD_C1.data());
	
	std::string histogramNameQCD_C1p = templatesQCD["C1p"][key_passed]->GetName();
	double normQCD_C1p = getIntegral(templatesQCD["C1p"][key_passed], true, true);	
	templatesQCD["C1p"][key_passed] = 
	  normalize(distributionsData[regionQCDtemplateFromData_passed][keyQCDtemplateFromData_passed], normQCD_C1p);
	templatesQCD["C1p"][key_passed]->SetName(histogramNameQCD_C1p.data());
	
	std::string histogramNameQCD_C1f = templatesQCD["C1f"][key_failed]->GetName();
	double normQCD_C1f = getIntegral(templatesQCD["C1f"][key_failed], true, true);
	templatesQCD["C1f"][key_failed] = 
	  normalize(distributionsData[regionQCDtemplateFromData_failed][keyQCDtemplateFromData_failed], normQCD_C1f);
	templatesQCD["C1f"][key_failed]->SetName(histogramNameQCD_C1f.data());
      }
    }

    templatesAll["QCD"] = templatesQCD;
  }

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
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	if ( !isSystematicShift(key->first) ) {
	  //std::string histogramTitle = std::string("Region ").append(region->first).append(": ").append(key->first);
	  //histogramTitle.append(" (scaled by cross-section)");
	  std::string histogramTitle = "";
	  std::string outputFileName = std::string("controlPlotsTauIdEff_wConstraints_");
	  outputFileName.append(region->first).append("_").append(key->first).append(".pdf");
	  drawHistograms(templatesZtautau[region->first][key->first], -1.,
			 templatesZmumu[region->first][key->first], -1.,
			 templatesQCD[region->first][key->first], -1.,
			 templatesWplusJets[region->first][key->first], -1.,
			 templatesTTplusJets[region->first][key->first], -1.,
			 distributionsData[region->first][key->first],
			 histogramTitle, xAxisTitles[getObservable(key->first)],
			 outputFileName);
	}
      }
    }
  }

  for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	region != distributionsData.end(); ++region ) {
    for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	  key != region->second.end(); ++key ) {
      std::cout << "numEvents[" << "Data" << "][" << region->first << "][" << key->first << "] = "
		<< getIntegral(key->second, true, true) << std::endl;
    }
  }

  std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsAll = // key = (process/"sum", region, observable)
    compNumEvents(templatesAll, processes, distributionsData);
  std::map<std::string, std::map<std::string, std::map<std::string, double> > > fittedFractions = // key = (process, region, observable)
    compFittedFractions(templatesAll, numEventsAll, processes, distributionsData);

//--- print MC expectations for probabilities 
//   o pDiTauCharge_OS_SS
//   o pDiTauKine_Sig_Bgr
//   o pMuonIso_tight_loose
//   o pTauId_passed_failed
//    separating different regions
  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      double numEventsABCD = numEventsAll[*process]["ABCD"][getKey("diTauMt", *tauId)];
      double numEventsA    = numEventsAll[*process]["A"][getKey("diTauMt", *tauId)];
      double numEventsB    = numEventsAll[*process]["B"][getKey("diTauMt", *tauId)];
      double numEventsC    = numEventsAll[*process]["C"][getKey("diTauMt", *tauId)];
      double numEventsC1   = numEventsAll[*process]["C1"][getKey("diTauMt", *tauId)];
      double numEventsC2   = numEventsAll[*process]["C2"][getKey("diTauMt", *tauId)];
      double numEventsD    = numEventsAll[*process]["D"][getKey("diTauMt", *tauId)];
      
      std::cout << "process = " << (*process) << std::endl;
      for ( size_t i = 0; i < (process->length() + 10); ++i ) {
	std::cout << "-";
      }
      std::cout << std::endl;
      
      std::cout << "numEventsABCD = " << numEventsABCD << std::endl;
      std::cout << "pDiTauCharge_OS_SS:" << std::endl;
      std::cout << " A/B = " << numEventsA/numEventsB << std::endl;
      std::cout << " C/D = " << numEventsC/numEventsD << std::endl;
      std::cout << "pMuonIso_tight_loose:" << std::endl;
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
		     *tauId, *fitVariable, fitTauIdEffC2,
		     effValue, effError, normFactors_fitted, hasFitConverged,		       
		     1);
      effValues[*tauId][*fitVariable] = effValue;
      effErrors[*tauId][*fitVariable] = effError;
      for ( std::vector<std::string>::const_iterator process = processes.begin();
	    process != processes.end(); ++process ) {
	for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	      region != distributionsData.end(); ++region ) {
	  normFactorsAll_fitted[*process][region->first][*tauId][*fitVariable] = normFactors_fitted[*process][region->first];
	}
      }
      fitConvergenceFlags[*tauId][*fitVariable] = hasFitConverged;
    }
  }

//--- make control plots of Data compared to sum(MC) scaled by normalization factors determined by fit
//    for muonPt, tauPt, Mt, visMass,... distributions in different regions
  if ( makeControlPlots ) {
    for ( std::map<std::string, std::map<std::string, TH1*> >::const_iterator region = distributionsData.begin();
	  region != distributionsData.end(); ++region ) {
      for ( std::map<std::string, TH1*>::const_iterator key = region->second.begin();
	    key != region->second.end(); ++key ) {
	for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	      tauId != tauIds.end(); ++tauId ) {
	  if ( !(key->first.find(*tauId) != std::string::npos) ) continue;
	  for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
		fitVariable != fitVariables.end(); ++fitVariable ) {
	    if ( !isSystematicShift(key->first) ) {
	      //std::string histogramTitle = std::string("Region ").append(region->first).append(": ").append(key->first);
	      //histogramTitle.append(" (scaled by normalization det. by fit)");
	      std::string histogramTitle = "";
	      std::string outputFileName = std::string("controlPlotsTauIdEff_wConstraints_");
	      outputFileName.append(region->first).append("_").append(key->first);
	      outputFileName.append("_fitted_").append(*fitVariable).append(".pdf");
  	      drawHistograms(
	        templatesZtautau[region->first][key->first], 
		getTemplateNorm_fitted("Ztautau", region->first, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesZmumu[region->first][key->first], 
		getTemplateNorm_fitted("Zmumu", region->first, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesQCD[region->first][key->first], 
		getTemplateNorm_fitted("QCD", region->first, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesWplusJets[region->first][key->first], 
		getTemplateNorm_fitted("WplusJets", region->first, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		templatesTTplusJets[region->first][key->first], 
		getTemplateNorm_fitted("TTplusJets", region->first, key->first, *tauId, *fitVariable, normFactorsAll_fitted, fittedFractions),
		distributionsData[region->first][key->first],
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
      
      double numEventsC    = numEventsAll["Ztautau"]["C"][getKey("diTauMt", *tauId)];
      double numEventsC1   = numEventsAll["Ztautau"]["C1"][getKey("diTauMt", *tauId)];
      double numEventsC1p  = numEventsAll["Ztautau"]["C1p"][getKey(fitVariables.front(), *tauId, "passed")];
      double numEventsC2p  = numEventsAll["Ztautau"]["C2p"][getKey("diTauMt", *tauId, "passed")];
      
      tauIdEffMCexp[*tauId][*fitVariable] = ( fitTauIdEffC2 ) ? (numEventsC1p + numEventsC2p)/numEventsC : numEventsC1p/numEventsC1;
      std::cout << "(Monte Carlo prediction = " << tauIdEffMCexp[*tauId][*fitVariable]*100. << "%)" << std::endl;
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
	  std::map<std::string, std::map<std::string, double> > normFactors_fitted;
	  bool hasFitConverged = false;
	  fitUsingRooFit(distributionsData, templatesAll_fluctuated, numEventsAll_fluctuated, fittedFractions_fluctuated,
			 processes,
			 *tauId, *fitVariable, fitTauIdEffC2,
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

      std::string fitNormC1Name  = std::string("fitNormC1_").append(*fitVariable).append("_").append(*tauId);
      std::string fitNormC1Title = std::string("Fitted Number of Z #rightarrow #tau^{+} #tau^{-} Events in region C1");
      TH1* histogramFitNormC1 = fitResultOutputDirectory.make<TH1F>(fitNormC1Name.data(), fitNormC1Title.data(), 1, -0.5, +0.5);
      int fitNormC1Bin = histogramFitNormC1->FindBin(0.);
      histogramFitNormC1->SetBinContent(fitNormC1Bin, normFactorsAll_fitted["Ztautau"]["C1"][*tauId][*fitVariable]);

      std::string expResultName  = std::string("expResult_").append(*fitVariable).append("_").append(*tauId);
      std::string expResultTitle = std::string("Expected ").append(*tauId).append(" Efficiency");
      TH1* histogramExpResult = fitResultOutputDirectory.make<TH1F>(expResultName.data(), expResultTitle.data(), 1, -0.5, +0.5);
      int expResultBin = histogramExpResult->FindBin(0.);
      histogramExpResult->SetBinContent(expResultBin, tauIdEffMCexp[*tauId][*fitVariable]);

      std::string expNormC1Name  = std::string("expNormC1_").append(*fitVariable).append("_").append(*tauId);
      std::string expNormC1Title = std::string("Expected Number of Z #rightarrow #tau^{+} #tau^{-} Events in region C1");
      TH1* histogramExpNormC1 = fitResultOutputDirectory.make<TH1F>(expNormC1Name.data(), expNormC1Title.data(), 1, -0.5, +0.5);
      int expNormC1Bin = histogramExpNormC1->FindBin(0.);
      histogramExpNormC1->SetBinContent(expNormC1Bin, numEventsAll["Ztautau"]["C1"][getKey("diTauMt", *tauId)]);
    }
  }
 
//--print time that it took macro to run
  std::cout << "finished executing fitTauIdEff_wConstraints macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  std::cout << " #sysShifts    = " << sysShifts.size() 
	    << " (numPseudoExperiments = " << numPseudoExperiments << ")" << std::endl;
  clock.Show("fitTauIdEff_wConstraints");

  return 0;
}
