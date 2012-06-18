#ifndef TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h
#define TauAnalysis_TauIdEfficiency_tauIdEffAuxFunctions_h

#include "TauAnalysis/DQMTools/interface/histogramAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <TDirectory.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TPaveText.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <time.h>

typedef std::map<std::string, TH1*> histogramMap1;
typedef std::map<std::string, histogramMap1> histogramMap2;
typedef std::map<std::string, histogramMap2> histogramMap3;
typedef std::map<std::string, histogramMap3> histogramMap4;

typedef std::map<std::string, double> valueMap1;
typedef std::map<std::string, valueMap1> valueMap2;
typedef std::map<std::string, valueMap2> valueMap3;

typedef std::vector<std::string> vstring;

const std::string key_central_value = "CENTRAL_VALUE";

//
//-------------------------------------------------------------------------------
//

bool contains_string(const vstring& names, const std::string& nameToFind)
{
  bool retVal = false;
  for ( vstring::const_iterator name = names.begin();
	name != names.end(); ++name ) {
    if ( (*name) == nameToFind ) retVal = true;
  }
  return retVal;
}

void add_string_uniquely(vstring& names, const std::string& nameToAdd)
{
  if ( !contains_string(names, nameToAdd) ) names.push_back(nameToAdd);
}

//
//-------------------------------------------------------------------------------
//

double getIntegral(const TH1* histogram, bool inclUnderflowBin, bool inclOverflowBin)
{
//--------------------------------------------------------------------------------
// Compute integral of histogram including/excluding underflow and overflow bins
//--------------------------------------------------------------------------------

  int firstBin = ( inclUnderflowBin ) ?                            0 :                      1;
  int lastBin  = ( inclOverflowBin  ) ? (histogram->GetNbinsX() + 1) : histogram->GetNbinsX();

  double integral = 0;
  for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    integral += binContent;
  }

  return integral;
}

TH1* normalize(const TH1* histogram, double norm = 1.)
{
//--------------------------------------------------------------------------------
// Normalize histogram passed as function argument to unit area
// (for use as shape template)
//--------------------------------------------------------------------------------

  TH1* retVal = (TH1*)histogram->Clone();

  if ( !retVal->GetSumw2N() ) retVal->Sumw2();

  double integral = getIntegral(retVal, true, true);
  if ( integral != 0. ) retVal->Scale(norm/integral);
    
  return retVal;
}

//
//-------------------------------------------------------------------------------
//

template <typename T>
void applyStyleOption(T* histogram, const std::string& histogramTitle,
		      const std::string& xAxisTitle, const std::string& yAxisTitle = "Events")
{
  histogram->SetTitle(histogramTitle.data());

  if ( histogram->GetXaxis() ) {
    histogram->GetXaxis()->SetTitle(xAxisTitle.data());
    histogram->GetXaxis()->SetTitleOffset(1.15);
  }

  if ( histogram->GetYaxis() ) {
    histogram->GetYaxis()->SetTitle(yAxisTitle.data());
    histogram->GetYaxis()->SetTitleOffset(1.65);
  }
}

void drawCMSprelimaryLabels(double intLumiData, double xOffset = 0.160, double yOffset = 0.8075)
{
  static TPaveText* cmsPreliminaryLabel = 0;
  if ( !cmsPreliminaryLabel ) {
    cmsPreliminaryLabel = new TPaveText(xOffset, yOffset + 0.0525, xOffset + 0.32, yOffset + 0.0925, "NDC");
    cmsPreliminaryLabel->AddText("CMS Preliminary 2012");
    cmsPreliminaryLabel->SetTextAlign(13);
    cmsPreliminaryLabel->SetTextSize(0.045);
    cmsPreliminaryLabel->SetFillStyle(0);
    cmsPreliminaryLabel->SetBorderSize(0);
  }
  cmsPreliminaryLabel->Draw();

  static TPaveText* cmsLuminosityLabel = 0;
  if ( !cmsLuminosityLabel ) {
    cmsLuminosityLabel = new TPaveText(xOffset + 0.005, yOffset, xOffset + 0.32, yOffset + 0.0400, "NDC");
    TString cmsLuminosityLabel_text = Form("#sqrt{s} = 8 TeV, L = %1.2f fb^{-1}", intLumiData);
    cmsLuminosityLabel->AddText(cmsLuminosityLabel_text.Data());
    cmsLuminosityLabel->SetTextAlign(13);
    cmsLuminosityLabel->SetTextSize(0.045);
    cmsLuminosityLabel->SetFillStyle(0);
    cmsLuminosityLabel->SetBorderSize(0);
  }
  cmsLuminosityLabel->Draw();
}

//
//-------------------------------------------------------------------------------
//

vstring getObservables(const std::string& region, const std::vector<std::string>& fitVariables)
{
  vstring retVal;
  if ( region == "ABCD" ) {
    retVal = fitVariables;
    add_string_uniquely(retVal, "diTauMt");
  } else if ( region.find("A") != std::string::npos ) {
    retVal = fitVariables;
    add_string_uniquely(retVal, "diTauMt");
  } else if ( region.find("B") != std::string::npos ) {
    retVal = fitVariables;
    add_string_uniquely(retVal, "diTauMt");
  } else if ( region.find("C") != std::string::npos ) {
    retVal = fitVariables;
    add_string_uniquely(retVal, "diTauMt");
  } else if ( region.find("D") != std::string::npos ) {
    retVal = fitVariables;
    add_string_uniquely(retVal, "diTauMt");
  } else {
    std::cout << "Error in <getObservables>: undefined region = " << region << " !!" << std::endl;
  }
  return retVal;
}

std::string getTauIdValue(const std::string& region)
{
  std::string retVal = "";
  if      ( region.find("_qcd") != std::string::npos ) retVal = "all";
  else if ( region.find("p")    != std::string::npos ) retVal = "passed";
  else if ( region.find("f")    != std::string::npos ) retVal = "failed";
  else                                                 retVal = "all";
  return retVal;
}

void loadHistograms(histogramMap3& histogramMap,
		    TDirectory* inputDirectory, const std::string& process, const vstring& regions,
		    const std::string& tauId, 
		    const vstring& fitVariables, 
		    const vstring& sysUncertainties, bool allowRebinning, bool applySmoothing, const std::string& genMatch = "")
{
//--------------------------------------------------------------------------------
// Load template histograms/distributions observed in data from ROOT file
//--------------------------------------------------------------------------------

  std::cout << "<loadHistograms>:" << std::endl;
  std::cout << " process = " << process << std::endl;
  std::cout << " regions = " << format_vstring(regions) << std::endl;
  std::cout << " tauId = " << tauId << std::endl;
  std::cout << " fitVariables = " << format_vstring(fitVariables) << std::endl;
  std::cout << " sysUncertainties = " << format_vstring(sysUncertainties) << std::endl;
  std::cout << " genMatch = " << genMatch << std::endl;
  
  for ( vstring::const_iterator region = regions.begin();
	region != regions.end(); ++region ) {

    vstring observables = getObservables(*region, fitVariables);
    add_string_uniquely(observables, "EventCounter"); // CV: for normalization purposes, always add 'EventCounter'
    std::string tauIdValue = getTauIdValue(*region);

    for ( vstring::const_iterator observable = observables.begin();
	  observable != observables.end(); ++observable ) {
      for ( vstring::const_iterator sysUncertainty = sysUncertainties.begin();
	    sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {

	std::string histogramName = std::string(process).append("_").append(*region).append("_").append(*observable);
	histogramName.append("_").append(tauId).append("_").append(tauIdValue);
	if ( (*sysUncertainty) != key_central_value ) histogramName.append("_").append(*sysUncertainty);
	if (   genMatch        != ""                ) histogramName.append("_").append(genMatch);
	
	// CV: switch to smoothed histograms if requested
	//    (for now the combinations of processes + regions for which template smoothing is applied
	//     is hardcoded here... this is to be changed later)
	if ( (process == "WplusJets" && (genMatch == "" || genMatch == "JetToTauFake") && 
	      ((*region) == "A" || (*region) == "A_mW" || (*region) == "A_mW" || 
	       (*region) == "B" ||
	       (*region) == "C1p" || (*region) == "C1f" ||
	       (*region) == "D")) ||
	     (process == "QCD" && 
	      ((*region) == "A" || (*region) == "A_mW" || (*region) == "A_mW" || 
	       (*region) == "B")) ) {
	  histogramName.append("_smoothed");
	  histogramName.append("__x"); // CV: this suffix is added by RooFit when running smoothTauIdEffTemplates macro
	}

	if ( histogramMap[*region][*observable].find(*sysUncertainty) != histogramMap[*region][*observable].end() ) {
	  std::cout << "Warning in <loadHistograms>:" 
		    << " histogram = " << histogramName << " already exists --> skipping !!" << std::endl;
	  continue;
	}
	
	TH1* histogram = dynamic_cast<TH1*>(inputDirectory->Get(histogramName.data()));
	if ( !histogram ) {
	  std::cout << "available histograms:" << std::endl;
	  inputDirectory->ls();
	  throw cms::Exception("loadHistograms")  
	    << "Failed to load histogram = " << histogramName << " from file/directory = " << inputDirectory->GetName() << " !!\n";
	}
	
	// CV: rebin histograms to avoid problem with too low Monte Carlo event statistics
	//     and large pile-up reweighting factors in Spring'12 MC production
	//if ( allowRebinning ) {
	//  int numBins = histogram->GetNbinsX();
	//  if      ( (numBins % 2) == 0 ) histogram->Rebin(2);
	//  else if ( (numBins % 3) == 0 ) histogram->Rebin(3);
	//}
	
	// CV: check that contents of all bins are positive,
	//     print warning if not
	int numBins = histogram->GetNbinsX();
	for ( int iBin = 0; iBin <= (numBins + 1); ++iBin ) {
	  double binContent = histogram->GetBinContent(iBin);
	  if ( binContent < 0. ) {
	    double x = histogram->GetBinCenter(iBin);
	    std::cout << "Warning in <loadHistograms>:" 
		      << " histogram = " << histogramName << ":" 
		      << " bin(x = " << x << ") = " << binContent << " --> setting it to 0." << std::endl;
	    histogram->SetBinContent(iBin, 0.);
	  }
	}
	
	if ( histogram != 0 ) {
	  histogramMap[*region][*observable][*sysUncertainty] = histogram;
	  std::cout << "histogramMap[region = " << (*region) << "][observable = " << (*observable) << "]" 
		    << "[sysUncertainty = " << (*sysUncertainty) << "] = " << histogram << std::endl;
	  std::cout << " (name = " << histogram->GetName() << ", integral = " << histogram->Integral() << ")" << std::endl;
	}
      }
    }
  }
}

//
//-------------------------------------------------------------------------------
//

std::pair<double, double> getNumber(TDirectory* inputDirectory, const TString& auxHistogramName, 
				    int auxHistogramBin, int assertAuxHistogramNumBins = 1)
{
  //std::cout << "<getNumber>:" << std::endl;
  //std::cout << " auxHistogramName = " << auxHistogramName << std::endl;
  //std::cout << " auxHistogramBin = " << auxHistogramBin << std::endl;
  //std::cout << " assertAuxHistogramNumBins = " << assertAuxHistogramNumBins << std::endl;

  TH1* histogram = dynamic_cast<TH1*>(inputDirectory->Get(auxHistogramName.Data()));
  if ( !histogram ) 
    throw cms::Exception("getNumber")  
      << "Failed to find histogram = " << auxHistogramName << " in input file/directory = " << inputDirectory->GetName() << " !!\n";
  
  int numBins = histogram->GetNbinsX();
  // CV: check that histogram has the expected number of bins,
  //     if it hasn't, FWLiteTauIdEffPreselNumbers/TauIdEffCutFlowTable has probably changed
  //     and the auxHistogramBin parameter needs to be updated !!
  //    (in particular if auxHistogramBin is negative)
  if ( numBins != assertAuxHistogramNumBins ) 
    throw cms::Exception("getNumber")  
      << "Histogram = " << auxHistogramName << " has incompatible binning:" 
      << " found = " << numBins << ", expected = " << assertAuxHistogramNumBins << " !!\n";

  int x;
  if   ( auxHistogramBin >= 0 ) x = auxHistogramBin;
  else                          x = numBins - TMath::Abs(auxHistogramBin);
  if ( !(x >= 0 && x < numBins) ) 
    throw cms::Exception("getNumber") 
      << "Invalid auxHistogramBin = " << auxHistogramBin << ", expected range = " << (-numBins) << ".." << (numBins - 1) << " !!\n";

  int bin = histogram->FindBin(x);

  std::pair<double, double> retVal;
  retVal.first = histogram->GetBinContent(bin);
  retVal.second = histogram->GetBinError(bin);

  //std::cout << "--> returning " << retVal.first << " +/- " << retVal.second << std::endl;

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

void printTimeStamp()
{
  time_t now_raw;
  time(&now_raw);
  struct tm* now_local = localtime(&now_raw);
  std::cout << "it is NOW: " << asctime(now_local) << std::endl;
}

#endif
