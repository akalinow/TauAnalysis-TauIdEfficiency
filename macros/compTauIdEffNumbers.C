
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

  double leadTrackFindingEff = 0.993;
  double leadTrackPtEff      = 0.927;
  double pfLooseIso06Eff     = 0.869;

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
  pfLooseIsoCorrFactor["HPSloose"]   = 1./0.971;
  pfLooseIsoCorrFactor["HPSmedium"]  = 1./0.976;
  pfLooseIsoCorrFactor["HPStight"]   = 1./0.976;

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
  leadTrackPtCorrFactor["HPSloose"]   = 1./0.931;
  leadTrackPtCorrFactor["HPSmedium"]  = 1./0.927;
  leadTrackPtCorrFactor["HPStight"]   = 1./0.926;

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
  mcExp["HPSloose"]   = 0.853;
  mcExp["HPSmedium"]  = 0.701;
  mcExp["HPStight"]   = 0.591;

//-------------------------------------------------------------------------------
// define numbers entering computation of  tau id. efficiency uncertainties
//-------------------------------------------------------------------------------

  std::vector<std::string> tauIds;
  //tauIds.push_back(std::string("TaNCloose"));
  //tauIds.push_back(std::string("TaNCmedium"));
  //tauIds.push_back(std::string("TaNCtight"));
  tauIds.push_back(std::string("HPSloose"));
  tauIds.push_back(std::string("HPSmedium"));
  tauIds.push_back(std::string("HPStight"));

  double leadTrackFindingEffErr = 0.039;

  double pfLooseIso06Eff1Err = 0.031;

  std::map<std::string, double> leadTrackPtEffErr;
  leadTrackPtEffErr["TaNCloose"]  = 0.;
  leadTrackPtEffErr["TaNCmedium"] = 0.;
  leadTrackPtEffErr["TaNCtight"]  = 0.;
  leadTrackPtEffErr["HPSloose"]   = (1./leadTrackPtCorrFactor["HPSloose"]  - 1.)*0.20;
  leadTrackPtEffErr["HPSmedium"]  = (1./leadTrackPtCorrFactor["HPSmedium"] - 1.)*0.20;
  leadTrackPtEffErr["HPStight"]   = (1./leadTrackPtCorrFactor["HPStight"]  - 1.)*0.20;

  std::map<std::string, double> pfLooseIso06Eff2Err;
  pfLooseIso06Eff2Err["TaNCloose"]  = (1./pfLooseIsoCorrFactor["TaNCloose"]  - 1.)*0.30;
  pfLooseIso06Eff2Err["TaNCmedium"] = (1./pfLooseIsoCorrFactor["TaNCmedium"] - 1.)*0.30;
  pfLooseIso06Eff2Err["TaNCtight"]  = (1./pfLooseIsoCorrFactor["TaNCtight"]  - 1.)*0.30;
  pfLooseIso06Eff2Err["HPSloose"]   = (1./pfLooseIsoCorrFactor["HPSloose"]   - 1.)*0.30;
  pfLooseIso06Eff2Err["HPSmedium"]  = (1./pfLooseIsoCorrFactor["HPSmedium"]  - 1.)*0.30;
  pfLooseIso06Eff2Err["HPStight"]   = (1./pfLooseIsoCorrFactor["HPStight"]   - 1.)*0.30;

  std::map<std::string, double> tauJetEnScaleErr;
  tauJetEnScaleErr["TaNCloose"]  = 0.010;
  tauJetEnScaleErr["TaNCmedium"] = 0.010;
  tauJetEnScaleErr["TaNCtight"]  = 0.010;
  tauJetEnScaleErr["HPSloose"]   = 0.010;
  tauJetEnScaleErr["HPSmedium"]  = 0.010;
  tauJetEnScaleErr["HPStight"]   = 0.010;

  std::map<std::string, double> fitStatErr;
  fitStatErr["TaNCloose"]  = 0.050;
  fitStatErr["TaNCmedium"] = 0.050;
  fitStatErr["TaNCtight"]  = 0.050;
  fitStatErr["HPSloose"]   = 0.050;
  fitStatErr["HPSmedium"]  = 0.050;
  fitStatErr["HPStight"]   = 0.050;

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
