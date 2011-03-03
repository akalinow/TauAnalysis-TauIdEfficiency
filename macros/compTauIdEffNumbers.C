
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

  double leadTrackFindingEff = 0.924;
  double leadTrackPtEff      = 0.812;
  double pfLooseIso06Eff     = 0.785;

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
  pfLooseIsoCorrFactor["TaNCloose"]  = 1./(1. - 0.128);
  pfLooseIsoCorrFactor["TaNCmedium"] = 1./(1. - 0.097);
  pfLooseIsoCorrFactor["TaNCtight"]  = 1./(1. - 0.072);
  pfLooseIsoCorrFactor["HPSloose"]   = 1./(1. - 0.054);
  pfLooseIsoCorrFactor["HPSmedium"]  = 1./(1. - 0.041);
  pfLooseIsoCorrFactor["HPStight"]   = 1./(1. - 0.033);

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
  leadTrackPtCorrFactor["HPSloose"]   = 1./(1. - 0.069);
  leadTrackPtCorrFactor["HPSmedium"]  = 1./(1. - 0.069);
  leadTrackPtCorrFactor["HPStight"]   = 1./(1. - 0.068);

  std::map<std::string, double> fitResult;
  fitResult["TaNCloose"]  = 0.813;
  fitResult["TaNCmedium"] = 0.649;
  fitResult["TaNCtight"]  = 0.568;
  fitResult["HPSloose"]   = 0.692;
  fitResult["HPSmedium"]  = 0.534;
  fitResult["HPStight"]   = 0.337;

  std::map<std::string, double> mcExp; 
  mcExp["TaNCloose"]  = 0.726;
  mcExp["TaNCmedium"] = 0.672;
  mcExp["TaNCtight"]  = 0.557;
  mcExp["HPSloose"]   = 0.689;
  mcExp["HPSmedium"]  = 0.519;
  mcExp["HPStight"]   = 0.370;

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

  double leadTrackFindingEffErr = 0.06;

  double pfLooseIso06Eff1Err = 0.02;

  std::map<std::string, double> leadTrackPtEffErr;
  leadTrackPtEffErr["TaNCloose"]  = 0.;
  leadTrackPtEffErr["TaNCmedium"] = 0.;
  leadTrackPtEffErr["TaNCtight"]  = 0.;
  leadTrackPtEffErr["HPSloose"]   = 0.069;
  leadTrackPtEffErr["HPSmedium"]  = 0.069;
  leadTrackPtEffErr["HPStight"]   = 0.068;

  std::map<std::string, double> pfLooseIso06Eff2Err;
  pfLooseIso06Eff2Err["TaNCloose"]  = 0.128;
  pfLooseIso06Eff2Err["TaNCmedium"] = 0.097;
  pfLooseIso06Eff2Err["TaNCtight"]  = 0.072;
  pfLooseIso06Eff2Err["HPSloose"]   = 0.054;
  pfLooseIso06Eff2Err["HPSmedium"]  = 0.041;
  pfLooseIso06Eff2Err["HPStight"]   = 0.033;

  std::map<std::string, double> tauJetEnScaleErr;
  tauJetEnScaleErr["TaNCloose"]  = 0.004;
  tauJetEnScaleErr["TaNCmedium"] = 0.005;
  tauJetEnScaleErr["TaNCtight"]  = 0.008;
  tauJetEnScaleErr["HPSloose"]   = 0.003;
  tauJetEnScaleErr["HPSmedium"]  = 0.001;
  tauJetEnScaleErr["HPStight"]   = 0.005;

  std::map<std::string, double> fitStatErr;
  fitStatErr["TaNCloose"]  = 0.210;
  fitStatErr["TaNCmedium"] = 0.140;
  fitStatErr["TaNCtight"]  = 0.157;
  fitStatErr["HPSloose"]   = 0.150;
  fitStatErr["HPSmedium"]  = 0.118;
  fitStatErr["HPStight"]   = 0.082;

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
