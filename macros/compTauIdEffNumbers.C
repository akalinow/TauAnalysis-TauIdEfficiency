
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>

#include <TMath.h>

double square(double x)
{
  return x*x;
}

void compTauIdEffNumbers()
{
//-------------------------------------------------------------------------------
// define numbers entering computation of tau id. efficiency "central values"
//-------------------------------------------------------------------------------

  double leadTrackFindingEff = 0.941;
  double leadTrackPtEff      = 0.797;
  double pfLooseIso06Eff     = 0.837;

//--- define correction factors for loose isolation requirement applied in preselection:
//
//     corrFactor = tauIdEff(MC, leadTrackFinding && leadTrackPtCut && TaNC/HPS discr. passed)
//                 ------------------------------------------------------------------------------------------
//                  tauIdEff(MC, leadTrackFinding && leadTrackPtCut && pfLooseIso06 passed) * tauIdEffFit(MC)
//
//    accounting for the fraction of tau-jets passing the TaNC/HPS discriminators,
//    but failing the loose isolation requirement
//
  std::map<std::string, double> pfLooseIsoCorrFactor;
  pfLooseIsoCorrFactor["TaNCloose"]  = 1./(1. - 0.085);
  pfLooseIsoCorrFactor["TaNCmedium"] = 1./(1. - 0.067);
  pfLooseIsoCorrFactor["TaNCtight"]  = 1./(1. - 0.051);
  pfLooseIsoCorrFactor["HPSloose"]   = 1./(1. - 0.041);
  pfLooseIsoCorrFactor["HPSmedium"]  = 1./(1. - 0.032);
  pfLooseIsoCorrFactor["HPStight"]   = 1./(1. - 0.027);

//--- define correction factors (specific to HPS)
//    for lead. track Pt cut applied in preselection, 
//    but not included in list of HPS discriminators used by Mike/recommended by Tau POG
//
//     corrFactor = tauIdEff(MC, leadTrackFinding && HPS discr. passed)
//                 ------------------------------------------------------------------------------------------
//                  tauIdEff(MC, leadTrackFinding && leadTrackPtCut && HPS discr. passed)
//
  std::map<std::string, double> leadTrackPtCorrFactor;
  leadTrackPtCorrFactor["TaNCloose"]  = 1.;
  leadTrackPtCorrFactor["TaNCmedium"] = 1.;
  leadTrackPtCorrFactor["TaNCtight"]  = 1.;
  leadTrackPtCorrFactor["HPSloose"]   = 1./(1. - 0.073);
  leadTrackPtCorrFactor["HPSmedium"]  = 1./(1. - 0.074);
  leadTrackPtCorrFactor["HPStight"]   = 1./(1. - 0.073);

  std::map<std::string, double> fitResult;
  fitResult["TaNCloose"]  = 0.760;
  fitResult["TaNCmedium"] = 0.634;
  fitResult["TaNCtight"]  = 0.554;
  fitResult["HPSloose"]   = 0.704;
  fitResult["HPSmedium"]  = 0.532;
  fitResult["HPStight"]   = 0.330;

  std::map<std::string, double> mcExp; 
  mcExp["TaNCloose"]  = 0.729;
  mcExp["TaNCmedium"] = 0.672;
  mcExp["TaNCtight"]  = 0.562;
  mcExp["HPSloose"]   = 0.701;
  mcExp["HPSmedium"]  = 0.528;
  mcExp["HPStight"]   = 0.357;

//-------------------------------------------------------------------------------
// define numbers entering computation of  tau id. efficiency uncertainties
//-------------------------------------------------------------------------------

  std::vector<std::string> tauIds;
  tauIds.push_back(std::string("TaNCloose"));
  tauIds.push_back(std::string("TaNCmedium"));
  tauIds.push_back(std::string("TaNCtight"));
  tauIds.push_back(std::string("HPSloose"));
  tauIds.push_back(std::string("HPSmedium"));
  tauIds.push_back(std::string("HPStight"));

  double leadTrackFindingEffErr = 0.04;

  double pfLooseIso06Eff1Err = 0.031;

  std::map<std::string, double> leadTrackPtEffErr;
  leadTrackPtEffErr["TaNCloose"]  = 0.;
  leadTrackPtEffErr["TaNCmedium"] = 0.;
  leadTrackPtEffErr["TaNCtight"]  = 0.;
  leadTrackPtEffErr["HPSloose"]   = 0.073;
  leadTrackPtEffErr["HPSmedium"]  = 0.074;
  leadTrackPtEffErr["HPStight"]   = 0.073;

  std::map<std::string, double> pfLooseIso06Eff2Err;
  pfLooseIso06Eff2Err["TaNCloose"]  = 0.085;
  pfLooseIso06Eff2Err["TaNCmedium"] = 0.067;
  pfLooseIso06Eff2Err["TaNCtight"]  = 0.051;
  pfLooseIso06Eff2Err["HPSloose"]   = 0.041;
  pfLooseIso06Eff2Err["HPSmedium"]  = 0.032;
  pfLooseIso06Eff2Err["HPStight"]   = 0.027;

  std::map<std::string, double> tauJetEnScaleErr;
  tauJetEnScaleErr["TaNCloose"]  = 0.004;
  tauJetEnScaleErr["TaNCmedium"] = 0.005;
  tauJetEnScaleErr["TaNCtight"]  = 0.008;
  tauJetEnScaleErr["HPSloose"]   = 0.003;
  tauJetEnScaleErr["HPSmedium"]  = 0.001;
  tauJetEnScaleErr["HPStight"]   = 0.005;

  std::map<std::string, double> fitStatErr;
  fitStatErr["TaNCloose"]  = 0.200;
  fitStatErr["TaNCmedium"] = 0.171;
  fitStatErr["TaNCtight"]  = 0.151;
  fitStatErr["HPSloose"]   = 0.154;
  fitStatErr["HPSmedium"]  = 0.126;
  fitStatErr["HPStight"]   = 0.084;

//-------------------------------------------------------------------------------
// compute tau id. efficiency "central values" and uncertainties
// for all TaNC/HPS working-points
//-------------------------------------------------------------------------------

  for ( std::vector<std::string>::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    std::cout << (*tauId) << ": stat. uncertainty = " << fitStatErr[*tauId]/fitResult[*tauId] << std::endl;

    double totErr2 = square(leadTrackFindingEffErr)
                    + square(pfLooseIso06Eff1Err) 
                    + square(leadTrackPtEffErr[*tauId])
                    + square(pfLooseIso06Eff2Err[*tauId])
                    + square(tauJetEnScaleErr[*tauId])
                    + square(fitStatErr[*tauId]/fitResult[*tauId]);
    std::cout << (*tauId) << ": total uncertainty = " << TMath::Sqrt(totErr2) << std::endl;

    std::cout << (*tauId) << ": loose iso. corr. factor = " << pfLooseIsoCorrFactor[*tauId] << std::endl;
    std::cout << (*tauId) << ": lead. track Pt corr. factor = " << leadTrackPtCorrFactor[*tauId] << std::endl;

    double totEff = leadTrackFindingEff
                   * leadTrackPtEff
                   * pfLooseIso06Eff
                   * pfLooseIsoCorrFactor[*tauId]
                   * leadTrackPtCorrFactor[*tauId]
                   * fitResult[*tauId];
    std::cout << (*tauId) << ": efficiency = " << totEff << " +/- " << totEff*TMath::Sqrt(totErr2) << std::endl;
    double totEffExp = leadTrackFindingEff
                   * leadTrackPtEff
                   * pfLooseIso06Eff
                   * pfLooseIsoCorrFactor[*tauId]
                   * leadTrackPtCorrFactor[*tauId]
                   * mcExp[*tauId];
    std::cout << (*tauId) << ": MC exp. = " << totEffExp << std::endl;
    std::cout << (*tauId) << ": Data/MC = " << totEff/totEffExp << " +/- " << totEff/totEffExp*TMath::Sqrt(totErr2) << std::endl;
  }
}
