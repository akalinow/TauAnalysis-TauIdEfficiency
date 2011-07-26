
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"

#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TString.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

double square(double x)
{
  return x*x;
}

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<compTauIdEffFinalNumbers>:" << std::endl;  

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("compTauChargeMisIdFinalNumbers");
  
//--- define all "global" numbers
//   (numbers not determined by tau id. efficiency code itself)
  double tauJetEnScaleErr = 0.010; // uncertainty on fitResults arising from tau-jet/jet energy scale uncertainty
  double frRelErr         = 0.20;  // relative uncertainty on jet --> tau fake-rate (for Z + jets events)

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("compTauChargeMisIdFinalNumbers") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgCompTauIdEffNumbers = cfg.getParameter<edm::ParameterSet>("compTauChargeMisIdFinalNumbers");

  typedef std::vector<std::string> vstring;
  vstring tauIds = cfgCompTauIdEffNumbers.getParameter<vstring>("tauIds");
  vstring fitVariables = cfgCompTauIdEffNumbers.getParameter<vstring>("fitVariables");

  std::string region_passed = cfgCompTauIdEffNumbers.getParameter<std::string>("passed_region");
  std::string region_failed = cfgCompTauIdEffNumbers.getParameter<std::string>("failed_region");

  std::string directory = cfgCompTauIdEffNumbers.getParameter<std::string>("directory");
  
  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("compTauChargeMisIdFinalNumbers") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string inputFileName = (*inputFiles.files().begin());

//--- open input file
  TFile* inputFile = TFile::Open(inputFileName.data());
  if ( !inputFile ) 
    throw cms::Exception("compTauChargeMisIdFinalNumbers") 
      << "Failed to open inputFile = " << inputFileName << " !!\n";

  std::string directory_fit = directory;
  TDirectory* inputDirectory_fit = ( directory_fit != "" ) ?
    dynamic_cast<TDirectory*>(inputFile->Get(directory_fit.data())) : inputFile;
  if ( !inputDirectory_fit ) 
    throw cms::Exception("compTauChargeMisIdFinalNumbers") 
      << "Directory = " << directory_fit << " does not exists in input file = " << inputFileName << " !!\n";

  std::string directory_presel = std::string("presel");
  if ( directory != "" ) directory_presel.append("/").append(directory);
  TDirectory* inputDirectory_presel = dynamic_cast<TDirectory*>(inputFile->Get(directory_presel.data()));
  if ( !inputDirectory_presel ) 
    throw cms::Exception("compTauChargeMisIdFinalNumbers") 
      << "Directory = " << directory_presel << " does not exists in input file = " << inputFileName << " !!\n";

  std::map<std::string, std::map<std::string, double> > expTauChargeMisIdRate;     // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > measTauChargeMisIdRate;    // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > measTauChargeMisIdRateErr; // key = (tauId, fitVariable)  

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {

      std::cout << "processing tauId = " << (*tauId) << ", fitVariable = " << (*fitVariable) << "..." << std::endl;

      TString key_passed = Form("_%s_", region_passed.data());
      TString key_failed = Form("_%s_", region_failed.data());

//--- get fitResult and uncertainty on fitResult
      TString fitResultName = Form("fitResult_%s_%s", fitVariable->data(), tauId->data());   
      double fitResult = getNumber(inputDirectory_fit, fitResultName, 0).first;
      double fitStatErr = getNumber(inputDirectory_fit, fitResultName, 0).second;
      std::cout << " fitResult = " << fitResult << " +/- " << fitStatErr << std::endl;

//--- get Monte Carlo expected fitResult
      TString mcExpName = Form("expResult_%s_%s", fitVariable->data(), tauId->data());  
      double mcExp = getNumber(inputDirectory_fit, mcExpName, 0).first;
      std::cout << " mcExp = " << mcExp << std::endl;

//--- compute "purity" correction factor to account for contribution of jet --> tau fakes
//    to Z --> tau+ tau- signal yield
      TString fitNormName = Form("fitNorm_%s_%s", fitVariable->data(), tauId->data());   
      double fitNorm = getNumber(inputDirectory_fit, fitNormName, 0).first;
      double fitNorm_passed = fitNorm*fitResult;
      double fitNorm_failed = fitNorm*(1. - fitResult);

      TString expNormName = Form("expNorm_%s_%s", fitVariable->data(), tauId->data()); 
      double expNorm = getNumber(inputDirectory_fit, expNormName, 0).first;
      double expNorm_passed = expNorm*mcExp;
      double expNorm_failed = expNorm*(1. - mcExp);

      TString matchedTauHadName_passed = Form("Ztautau_%s_%s_TauHadMatched", region_passed.data(), tauId->data());  
      double matchedTauHad_passed = getNumber(inputDirectory_presel, matchedTauHadName_passed, 7, 9).first;
      TString matchedFakeTauName_passed = Form("Ztautau_%s_%s_FakeTauMatched", region_passed.data(), tauId->data());
      double matchedFakeTau_passed = getNumber(inputDirectory_presel, matchedFakeTauName_passed, 7, 9).first;
      double expPurity_passed = matchedTauHad_passed/(matchedTauHad_passed + matchedFakeTau_passed);
      std::cout << " expPurity_passed = " << expPurity_passed << std::endl;
      double expFake_passed = expNorm_passed*(1. - expPurity_passed);

      TString matchedTauHadName_failed = TString(matchedTauHadName_passed).ReplaceAll(key_passed, key_failed);
      double matchedTauHad_failed = getNumber(inputDirectory_presel, matchedTauHadName_failed, 5, 9).first;
      TString matchedFakeTauName_failed = TString(matchedFakeTauName_passed).ReplaceAll(key_passed, key_failed);
      double matchedFakeTau_failed = getNumber(inputDirectory_presel, matchedFakeTauName_failed, 5, 9).first;
      double expPurity_failed = matchedTauHad_failed/(matchedTauHad_failed + matchedFakeTau_failed);
      std::cout << " expPurity_failed = " << expPurity_failed << std::endl;
      double expFake_failed = expNorm_failed*(1. - expPurity_failed);

      double numPurityCorrFactor = 
	(fitNorm_passed - expFake_passed)/((fitNorm_passed - expFake_passed) + (fitNorm_failed - expFake_failed));
      double denomPurityCorrFactor = fitNorm_passed/fitNorm;
      double purityCorrFactor = numPurityCorrFactor/denomPurityCorrFactor;
      double purityCorrFactorErr = 
	purityCorrFactor
       *((expFake_passed + expFake_failed)/((fitNorm_passed - expFake_passed) + (fitNorm_failed - expFake_failed)))
       *frRelErr;
      std::cout << " purityCorrFactor = " << purityCorrFactor << " +/- " << purityCorrFactorErr 
		<< " (" << purityCorrFactorErr/purityCorrFactor << "%)" << std::endl;
      
      double numExpPurityCorrFactor = 
	(expNorm_passed - expFake_passed)/((expNorm_passed - expFake_passed) + (expNorm_failed - expFake_failed));
      double denomExpPurityCorrFactor = expNorm_passed/expNorm;
      double expPurityCorrFactor = numExpPurityCorrFactor/denomExpPurityCorrFactor;
      std::cout << " expPurityCorrFactor = " << expPurityCorrFactor << std::endl;

//--- compute tau id. efficiency and uncertainty
      std::cout << "stat. uncertainty = " << fitStatErr/fitResult << std::endl;

      double totErr2 = square(purityCorrFactorErr/purityCorrFactor)
                      + square(tauJetEnScaleErr)
                      + square(fitStatErr/fitResult);
      std::cout << "total uncertainty = " << TMath::Sqrt(totErr2) << std::endl;

      double totRate = 1. - (purityCorrFactor * fitResult);
      std::cout << "--> tau charge misid. rate = " << totRate << " +/- " << totRate*TMath::Sqrt(totErr2) << std::endl;
      double totRateExp = 1. - (expPurityCorrFactor * mcExp);
      std::cout << " MC exp. (1) = " << totRateExp << std::endl;
      std::cout << "--> Data/MC = " << totRate/totRateExp << " +/- " << (totRate/totRateExp)*TMath::Sqrt(totErr2) << std::endl;

      expTauChargeMisIdRate[*tauId][*fitVariable] = totRateExp;

      measTauChargeMisIdRate[*tauId][*fitVariable] = totRate;
      if ( measTauChargeMisIdRate[*tauId][*fitVariable] < 0. ) measTauChargeMisIdRate[*tauId][*fitVariable] = 0.;
      measTauChargeMisIdRateErr[*tauId][*fitVariable] = measTauChargeMisIdRate[*tauId][*fitVariable]*TMath::Sqrt(totErr2);

//--- CV: cross-check Monte Carlo expected tau charge misidentification rate
//        with FWLiteTauIdEffPreselNumbers output
//       ('cc': reconstructed tau-jet matches hadronic tau decay on generator level,
//              charge of tau-jet matches charge of tau lepton on generator level
//       ('wc': reconstructed tau-jet matches hadronic tau decay on generator level,
//              charge of tau-jet does not match charge of tau lepton on generator level)
      TString totRateExpName_control_cc = Form("Ztautau_%s_%s_TauHadMatchedCorrectCharge", region_passed.data(), tauId->data()); 
      double totRateExp_control_cc = getNumber(inputDirectory_presel, totRateExpName_control_cc, 4, 6).first;
      std::cout << " totRateExp_control_cc = " << totRateExp_control_cc << std::endl;
      TString totRateExpName_control_wc = TString(totRateExpName_control_cc).ReplaceAll("CorrectCharge", "WrongCharge");
      double totRateExp_control_wc = getNumber(inputDirectory_presel, totRateExpName_control_wc, 4, 6).first;
      std::cout << " totRateExp_control_wc = " << totRateExp_control_wc << std::endl;
      double totRateExp_control = 1. - totRateExp_control_cc/(totRateExp_control_cc + totRateExp_control_wc);
      std::cout << " MC exp. (2) = " << totRateExp_control << std::endl;

      std::cout << std::endl;
    }
  }

//-- save expected and measured tau id. efficiency values
  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory outputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {
      std::string expRateName = std::string("expRate_").append(*fitVariable).append("_").append(*tauId);
      std::string expRateTitle = std::string("Expected Tau charge misid. Rate after ").append(*tauId);
      TH1* histogramExpRate = outputDirectory.make<TH1F>(expRateName.data(), expRateTitle.data(), 1, -0.5, +0.5);
      int expRateBin = histogramExpRate->FindBin(0.);
      histogramExpRate->SetBinContent(expRateBin, expTauChargeMisIdRate[*tauId][*fitVariable]);

      std::string measRateName = std::string("measRate_").append(*fitVariable).append("_").append(*tauId);
      std::string measRateTitle = std::string("Measured Tau charge misid. Rate after ").append(*tauId);
      measRateTitle.append(", obtained by fitting ").append(*fitVariable);
      TH1* histogramMeasRate = outputDirectory.make<TH1F>(measRateName.data(), measRateTitle.data(), 1, -0.5, +0.5);
      int measRateBin = histogramMeasRate->FindBin(0.);
      histogramMeasRate->SetBinContent(measRateBin, measTauChargeMisIdRate[*tauId][*fitVariable]);
      histogramMeasRate->SetBinError(measRateBin, measTauChargeMisIdRateErr[*tauId][*fitVariable]);
    }
  }

//--print time that it took macro to run
  std::cout << "finished executing compTauChargeMisIdFinalNumbers macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("compTauChargeMisIdFinalNumbers");

  return 0;
}
