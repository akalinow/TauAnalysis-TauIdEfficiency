
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
  clock.Start("compTauIdEffFinalNumbers");
  
//--- define all "global" numbers
//   (numbers not determined by tau id. efficiency code itself)
  double leadTrackFindingEffErr = 0.039; // uncertainty on track reconstruction efficiency for charged hadrons (taken from TRK-10-002)
  double leadTrackPtEffErr      = 0.010; // relative uncertainty on efficiency to pass leadTrackPt > 5.0 GeV requirement
  double pfLooseIsoEffErr       = 0.025; // difference between efficiencies to pass pfLooseIso for e/mu/tau-jets 
                                         // plus Data-MC difference in e/mu efficiencies (measured in Zee/Zmumu events via Tag & Probe)
  double tauJetEnScaleErr       = 0.010; // uncertainty on fitResults arising from tau-jet/jet energy scale uncertainty
  double leadTrackPtCorrRelErr  = 0.20;  // relative uncertainty on leadTrackPt correction factor
  double pfLooseIsoCorrRelErr   = 0.30;  // relative uncertainty on pfLooseIso correction factor
  double frRelErr               = 0.20;  // relative uncertainty on jet --> tau fake-rate (for Z + jets events)

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("compTauIdEffFinalNumbers") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgCompTauIdEffNumbers = cfg.getParameter<edm::ParameterSet>("compTauIdEffFinalNumbers");

  typedef std::vector<std::string> vstring;
  vstring tauIds = cfgCompTauIdEffNumbers.getParameter<vstring>("tauIds");
  vstring fitVariables = cfgCompTauIdEffNumbers.getParameter<vstring>("fitVariables");

  std::string region_passed = cfgCompTauIdEffNumbers.getParameter<std::string>("passed_region");
  std::string region_failed = cfgCompTauIdEffNumbers.getParameter<std::string>("failed_region");

  std::string directory = cfgCompTauIdEffNumbers.getParameter<std::string>("directory");
  
  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("compTauIdEffFinalNumbers") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string inputFileName = (*inputFiles.files().begin());

//--- open input file
  TFile* inputFile = TFile::Open(inputFileName.data());
  if ( !inputFile ) 
    throw cms::Exception("compTauIdEffFinalNumbers") 
      << "Failed to open inputFile = " << inputFileName << " !!\n";

  std::string directory_fit = directory;
  TDirectory* inputDirectory_fit = ( directory_fit != "" ) ?
    dynamic_cast<TDirectory*>(inputFile->Get(directory_fit.data())) : inputFile;
  if ( !inputDirectory_fit ) 
    throw cms::Exception("compTauIdEffFinalNumbers") 
      << "Directory = " << directory_fit << " does not exists in input file = " << inputFileName << " !!\n";

  std::string directory_presel = std::string("presel");
  if ( directory != "" ) directory_presel.append("/").append(directory);
  TDirectory* inputDirectory_presel = dynamic_cast<TDirectory*>(inputFile->Get(directory_presel.data()));
  if ( !inputDirectory_presel ) 
    throw cms::Exception("compTauIdEffFinalNumbers") 
      << "Directory = " << directory_presel << " does not exists in input file = " << inputFileName << " !!\n";

  std::map<std::string, std::map<std::string, double> > expTauIdEfficiency;     // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > measTauIdEfficiency;    // key = (tauId, fitVariable)
  std::map<std::string, std::map<std::string, double> > measTauIdEfficiencyErr; // key = (tauId, fitVariable)  

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( std::vector<std::string>::const_iterator fitVariable = fitVariables.begin();
	  fitVariable != fitVariables.end(); ++fitVariable ) {

      std::cout << "processing tauId = " << (*tauId) << ", fitVariable = " << (*fitVariable) << "..." << std::endl;

//--- compute efficiencies and uncertainties of lead. track finding, track Pt cut 
//    and loose isolation requirement applied in preselection
      TString effPreselectionName_passed = Form("Ztautau_%s_%s_TauHadMatched", region_passed.data(), tauId->data());  
      TString key_passed = Form("_%s_", region_passed.data());
      TString key_failed = Form("_%s_", region_failed.data());
      TString effPreselectionName_failed = TString(effPreselectionName_passed).ReplaceAll(key_passed, key_failed);
      
      double numLeadTrackFindingEff = 
	getNumber(inputDirectory_presel, effPreselectionName_passed, 1, 9).first
       + getNumber(inputDirectory_presel, effPreselectionName_failed, 1, 9).first;
      double denomLeadTrackFindingEff = 
        getNumber(inputDirectory_presel, effPreselectionName_passed, 0, 9).first
       + getNumber(inputDirectory_presel, effPreselectionName_failed, 0, 9).first;
      double leadTrackFindingEff = numLeadTrackFindingEff/denomLeadTrackFindingEff;
      std::cout << " leadTrackFindingEff = " << leadTrackFindingEff  << " +/- " << leadTrackFindingEffErr << std::endl;

      double numLeadTrackPtEff = 
	getNumber(inputDirectory_presel, effPreselectionName_passed, 2, 9).first
       + getNumber(inputDirectory_presel, effPreselectionName_failed, 2, 9).first;
      double denomLeadTrackPtEff = numLeadTrackFindingEff;
      double leadTrackPtEff = numLeadTrackPtEff/denomLeadTrackPtEff;
      std::cout << " leadTrackPtEff = " << leadTrackPtEff << " +/- " << leadTrackPtEffErr << std::endl;

      double numPFLooseIsoEff = 
	getNumber(inputDirectory_presel, effPreselectionName_passed, 3, 9).first;
       + getNumber(inputDirectory_presel, effPreselectionName_failed, 3, 9).first;
      double denomPFLooseIsoEff = numLeadTrackPtEff;
      double pfLooseIsoEff = numPFLooseIsoEff/denomPFLooseIsoEff;
      std::cout << " pfLooseIsoEff = " << pfLooseIsoEff << " +/- " << pfLooseIsoEffErr << std::endl;

      std::cout << "--> preselection efficiency = " << (leadTrackFindingEff*leadTrackPtEff*pfLooseIsoEff) << std::endl;

//--- get fitResult and uncertainty on fitResult
      TString fitResultName = Form("fitResult_%s_%s", fitVariable->data(), tauId->data());   
      double fitResult = getNumber(inputDirectory_fit, fitResultName, 0).first;
      double fitStatErr = getNumber(inputDirectory_fit, fitResultName, 0).second;
      std::cout << " fitResult = " << fitResult << " +/- " << fitStatErr << std::endl;

//--- get Monte Carlo expected fitResult
      TString mcExpName = Form("expResult_%s_%s", fitVariable->data(), tauId->data());  
      double mcExp = getNumber(inputDirectory_fit, mcExpName, 0).first;
      std::cout << " mcExp = " << mcExp << std::endl;

//--- compute correction factors for loose isolation requirement applied in preselection:
//
//     corrFactor = tauIdEff(MC, leadTrackFinding && leadTrackPtCut && TaNC/HPS discr. passed)
//                 ------------------------------------------------------------------------------------------
//                  tauIdEff(MC, leadTrackFinding && leadTrackPtCut && pfLooseIso06 passed) * tauIdEffFit(MC)
//
//    accounting for the fraction of tau-jets passing the TaNC/HPS discriminators,
//    but failing the loose isolation requirement
//
      TString pfLooseIsoCorrFactorName = Form("Ztautau_%s_%s_TauHadMatchedReversed", region_passed.data(), tauId->data());   
      double numPFLooseIsoCorrFactor = getNumber(inputDirectory_presel, pfLooseIsoCorrFactorName, -5, 9).first;
      //std::cout << " numPFLooseIsoCorrFactor = " << numPFLooseIsoCorrFactor << std::endl;
      double denomPFLooseIsoCorrFactor = getNumber(inputDirectory_presel, pfLooseIsoCorrFactorName, -4, 9).first;
      //std::cout << " denomPFLooseIsoCorrFactor = " << denomPFLooseIsoCorrFactor << std::endl;
      double pfLooseIsoCorrFactor = numPFLooseIsoCorrFactor/denomPFLooseIsoCorrFactor;
      double pfLooseIsoCorrFactorErr = (1./pfLooseIsoCorrFactor  - 1.)*pfLooseIsoCorrRelErr;
      std::cout << " pfLooseIsoCorrFactor = " << pfLooseIsoCorrFactor << " +/- " << pfLooseIsoCorrFactorErr 
		<< " (" << pfLooseIsoCorrFactorErr/pfLooseIsoCorrFactor << "%)" << std::endl;

//--- compute correction factor (specific to HPS)
//    for lead. track Pt cut applied in preselection, 
//    but not included in list of HPS discriminators used by Mike/recommended by Tau POG
//
//     corrFactor = tauIdEff(MC, leadTrackFinding && HPS discr. passed)
//                 ------------------------------------------------------------------------------------------
//                  tauIdEff(MC, leadTrackFinding && leadTrackPtCut && HPS discr. passed)
//
      TString leadTrackPtCorrFactorName = pfLooseIsoCorrFactorName;
      double numLeadTrackPtCorrFactor = getNumber(inputDirectory_presel, pfLooseIsoCorrFactorName, -6, 9).first;
      //std::cout << " numLeadTrackPtCorrFactor = " << numLeadTrackPtCorrFactor << std::endl;
      double denomLeadTrackPtCorrFactor = getNumber(inputDirectory_presel, pfLooseIsoCorrFactorName, -5, 9).first;
      //std::cout << " denomLeadTrackPtCorrFactor = " << denomLeadTrackPtCorrFactor << std::endl;
      double leadTrackPtCorrFactor = numLeadTrackPtCorrFactor/denomLeadTrackPtCorrFactor;
      double leadTrackPtCorrFactorErr = (1./leadTrackPtCorrFactor - 1.)*leadTrackPtCorrRelErr;
      std::cout << " leadTrackPtCorrFactor = " << leadTrackPtCorrFactor << " +/- " << leadTrackPtCorrFactorErr 
		<< " (" << leadTrackPtCorrFactorErr/leadTrackPtCorrFactor << "%)" << std::endl;

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

      double totErr2 = square(leadTrackFindingEffErr)
	              + square(pfLooseIsoEffErr) 
  	              + square(leadTrackPtEffErr)
                      + square(pfLooseIsoCorrFactorErr/pfLooseIsoCorrFactor)
 	              + square(purityCorrFactorErr/purityCorrFactor)
                      + square(tauJetEnScaleErr)
                      + square(fitStatErr/fitResult);
      std::cout << "total uncertainty = " << TMath::Sqrt(totErr2) << std::endl;

      double totEff = 
	leadTrackFindingEff
       * leadTrackPtEff
       * pfLooseIsoEff
       * pfLooseIsoCorrFactor
       * leadTrackPtCorrFactor
       * purityCorrFactor
       * fitResult;
      std::cout << "--> efficiency = " << totEff << " +/- " << totEff*TMath::Sqrt(totErr2) << std::endl;
      double totEffExp = 
	leadTrackFindingEff
       * leadTrackPtEff
       * pfLooseIsoEff
       * pfLooseIsoCorrFactor
       * leadTrackPtCorrFactor
       * expPurityCorrFactor
       * mcExp;
      std::cout << " MC exp. (1) = " << totEffExp << std::endl;
      std::cout << "--> Data/MC = " << totEff/totEffExp << " +/- " << totEff/totEffExp*TMath::Sqrt(totErr2) << std::endl;

      expTauIdEfficiency[*tauId][*fitVariable] = totEffExp;

      measTauIdEfficiency[*tauId][*fitVariable] = totEff;
      if ( measTauIdEfficiency[*tauId][*fitVariable] > 1. ) measTauIdEfficiency[*tauId][*fitVariable] = 1.;
      measTauIdEfficiencyErr[*tauId][*fitVariable] = measTauIdEfficiency[*tauId][*fitVariable]*TMath::Sqrt(totErr2);

//--- CV: cross-check Monte Carlo expected tau id. efficiency
//        with FWLiteTauIdEffPreselNumbers output
//
//    NOTE: quoted efficiencies are for "true" hadronic tau decays of Pt > 20 GeV, abs(eta) < 2.3
//          which are reconstructed as PFJets of Pt > 20 GeV, abs(eta) < 2.3
//         (in order to get "full" efficiencies, need to multiply the quoted efficiencies
//          by the probability to reconstruct PFJet and to have Pt > 20 GeV on reconstruction level)
//
      TString totEffExpName_control_passed = Form("Ztautau_%s_%s_TauHadMatchedReversed", region_passed.data(), tauId->data()); 
      TString totEffExpName_control_failed = TString(totEffExpName_control_passed).ReplaceAll(key_passed, key_failed);
      double numTotEffExp_control = 
	getNumber(inputDirectory_presel, totEffExpName_control_passed, -7, 9).first
       + getNumber(inputDirectory_presel, totEffExpName_control_failed, -7, 9).first;
      //std::cout << " numTotEffExp_control = " << numTotEffExp_control << std::endl;
      double denomTotEffExp_control = 
	getNumber(inputDirectory_presel, totEffExpName_control_passed, 0, 9).first
       + getNumber(inputDirectory_presel, totEffExpName_control_failed, 0, 9).first;
      //std::cout << " denomTotEffExp_control = " << denomTotEffExp_control << std::endl;
      double totEffExp_control = numTotEffExp_control/denomTotEffExp_control;
      std::cout << " MC exp. (2) = " << totEffExp_control << std::endl;

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
      std::string expEffName = std::string("expEff_").append(*fitVariable).append("_").append(*tauId);
      std::string expEffTitle = std::string("Expected efficiency of ").append(*tauId);
      TH1* histogramExpEff = outputDirectory.make<TH1F>(expEffName.data(), expEffTitle.data(), 1, -0.5, +0.5);
      int expEffBin = histogramExpEff->FindBin(0.);
      histogramExpEff->SetBinContent(expEffBin, expTauIdEfficiency[*tauId][*fitVariable]);

      std::string measEffName = std::string("measEff_").append(*fitVariable).append("_").append(*tauId);
      std::string measEffTitle = std::string("Measured Efficiency of ").append(*tauId);
      measEffTitle.append(", obtained by fitting ").append(*fitVariable);
      TH1* histogramMeasEff = outputDirectory.make<TH1F>(measEffName.data(), measEffTitle.data(), 1, -0.5, +0.5);
      int measEffBin = histogramMeasEff->FindBin(0.);
      histogramMeasEff->SetBinContent(measEffBin, measTauIdEfficiency[*tauId][*fitVariable]);
      histogramMeasEff->SetBinError(measEffBin, measTauIdEfficiencyErr[*tauId][*fitVariable]);
    }
  }

//--print time that it took macro to run
  std::cout << "finished executing compTauIdEffFinalNumbers macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("compTauIdEffFinalNumbers");

  return 0;
}
