#include "TauAnalysis/TauIdEfficiency/interface/TauPtResHistManager.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/TauIdEfficiency/interface/tauPtResAuxFunctions.h"

#include <TMath.h>

#include <iostream>
#include <iomanip>

TauPtResHistManager::tauPtResManCorrHistograms::tauPtResManCorrHistograms(int level)
  : level_(level)
{}

void TauPtResHistManager::tauPtResManCorrHistograms::bookHistograms(TFileDirectory& dir)
{
  std::ostringstream level_string;
  for ( unsigned iBit = 0; iBit < 8; ++iBit ) {
    if ( level_ & (1 << iBit) ) level_string << (iBit + 1);
  }
  //std::cout << "level = " << level_ << ": level_string = " << level_string.str() << std::endl;

  histogramTauPtRes_                  = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%s", level_string.str().data()),                  "tauPtRes",                               40, 0., 2.);
  histogramTauPtResGenOneProng0Pi0_   = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%sGenOneProng0Pi0", level_string.str().data()),   "tauPtRes (gen. one-prong, 0 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenOneProng1Pi0_   = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%sGenOneProng1Pi0", level_string.str().data()),   "tauPtRes (gen. one-prong, 1 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenOneProng2Pi0_   = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%sGenOneProng2Pi0", level_string.str().data()),   "tauPtRes (gen. one-prong, 2 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenThreeProng0Pi0_ = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%sGenThreeProng0Pi0", level_string.str().data()), "tauPtRes (gen. three-prong, 0 #pi^{0})", 40, 0., 2.);
  histogramTauPtResGenThreeProng1Pi0_ = TauPtResHistManager::book1D(
    dir, Form("tauPtResManCorrLev%sGenThreeProng1Pi0", level_string.str().data()), "tauPtRes (gen. three-prong, 1 #pi^{0})", 40, 0., 2.);
}

void TauPtResHistManager::tauPtResManCorrHistograms::fillHistograms(
       const pat::Tau& patTau, const std::string& genTauDecayMode, double genVisPt, const reco::Vertex& vertex, 
       const reco::tau::RecoTauQualityCuts& qualityCuts, double weight)
{
  double tauPtResManCorr = (patTau.pt() + getTauPtManCorr(patTau, vertex, qualityCuts, level_))/genVisPt;
  histogramTauPtRes_->Fill(tauPtResManCorr, weight);
  if ( genTauDecayMode == "oneProng0Pi0"   ) histogramTauPtResGenOneProng0Pi0_->Fill(tauPtResManCorr, weight);
  if ( genTauDecayMode == "oneProng1Pi0"   ) histogramTauPtResGenOneProng1Pi0_->Fill(tauPtResManCorr, weight);
  if ( genTauDecayMode == "oneProng2Pi0"   ) histogramTauPtResGenOneProng2Pi0_->Fill(tauPtResManCorr, weight);
  if ( genTauDecayMode == "threeProng0Pi0" ) histogramTauPtResGenThreeProng0Pi0_->Fill(tauPtResManCorr, weight);
  if ( genTauDecayMode == "threeProng1Pi0" ) histogramTauPtResGenThreeProng1Pi0_->Fill(tauPtResManCorr, weight);
}

//
//-------------------------------------------------------------------------------
//

TauPtResHistManager::TauPtResHistManager(const edm::ParameterSet& cfg)
  : histogramsTauPtResManCorrLev1_(1),
    histogramsTauPtResManCorrLev2_(2),
    histogramsTauPtResManCorrLev12_(3),
    histogramsTauPtResManCorrLev123_(7),
    histogramsTauPtResManCorrLev14_(9)
{}

TauPtResHistManager::~TauPtResHistManager()
{
// nothing to be done yet...
}

void TauPtResHistManager::bookHistograms(TFileDirectory& dir)
{
  histogramTauPtRes_                  = book1D(dir, "tauPtRes",                  "tauPtRes",                               40, 0., 2.);
  histogramTauPtResGenOneProng0Pi0_   = book1D(dir, "tauPtResGenOneProng0Pi0",   "tauPtRes (gen. one-prong, 0 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenOneProng1Pi0_   = book1D(dir, "tauPtResGenOneProng1Pi0",   "tauPtRes (gen. one-prong, 1 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenOneProng2Pi0_   = book1D(dir, "tauPtResGenOneProng2Pi0",   "tauPtRes (gen. one-prong, 2 #pi^{0})",   40, 0., 2.);
  histogramTauPtResGenThreeProng0Pi0_ = book1D(dir, "tauPtResGenThreeProng0Pi0", "tauPtRes (gen. three-prong, 0 #pi^{0})", 40, 0., 2.);
  histogramTauPtResGenThreeProng1Pi0_ = book1D(dir, "tauPtResGenThreeProng1Pi0", "tauPtRes (gen. three-prong, 1 #pi^{0})", 40, 0., 2.);

  histogramEtaPhiForTauPtResLt05GenOneProng0Pi0_ = 
    book2D(dir, "etaPhiForTauPtResLt05GenOneProng0Pi0", 
	   "#eta vs. #phi (tauPtRes < 0.5, gen. one-prong, 0 #pi^{0})", 25, -2.5, +2.5, 36, -TMath::Pi(), +TMath::Pi());
  histogramEtaPhiForTauPtResLt05GenOneProng1Pi0_ = 
    book2D(dir, "etaPhiForTauPtResLt05GenOneProng1Pi0", 
	   "#eta vs. #phi (tauPtRes < 0.5, gen. one-prong, 1 #pi^{0})", 25, -2.5, +2.5, 36, -TMath::Pi(), +TMath::Pi());
  histogramEtaPhiForTauPtResLt05GenOneProng2Pi0_ = 
    book2D(dir, "etaPhiForTauPtResLt05GenOneProng2Pi0", 
	   "#eta vs. #phi (tauPtRes < 0.5, gen. one-prong, 2 #pi^{0})", 25, -2.5, +2.5, 36, -TMath::Pi(), +TMath::Pi());
  histogramEtaPhiForTauPtResLt05GenThreeProng0Pi0_ = 
    book2D(dir, "etaPhiForTauPtResLt05GenThreeProng0Pi0", 
	   "#eta vs. #phi (tauPtRes < 0.5, gen. three-prong, 0 #pi^{0})", 25, -2.5, +2.5, 36, -TMath::Pi(), +TMath::Pi());
  histogramEtaPhiForTauPtResLt05GenThreeProng1Pi0_ = 
    book2D(dir, "etaPhiForTauPtResLt05GenThreeProng1Pi0", 
	   "#eta vs. #phi (tauPtRes < 0.5, gen. three-prong, 1 #pi^{0})", 25, -2.5, +2.5, 36, -TMath::Pi(), +TMath::Pi());

  histogramsTauPtResManCorrLev1_.bookHistograms(dir);
  histogramsTauPtResManCorrLev2_.bookHistograms(dir);
  histogramsTauPtResManCorrLev12_.bookHistograms(dir);
  histogramsTauPtResManCorrLev123_.bookHistograms(dir);
  histogramsTauPtResManCorrLev14_.bookHistograms(dir);
 
  histogramJetPtRes_                  = book1D(dir, "jetPtRes",                  "jetPtRes",                               40, 0., 2.);
  histogramJetPtResGenOneProng0Pi0_   = book1D(dir, "jetPtResGenOneProng0Pi0",   "jetPtRes (gen. one-prong, 0 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenOneProng1Pi0_   = book1D(dir, "jetPtResGenOneProng1Pi0",   "jetPtRes (gen. one-prong, 1 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenOneProng2Pi0_   = book1D(dir, "jetPtResGenOneProng2Pi0",   "jetPtRes (gen. one-prong, 2 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenThreeProng0Pi0_ = book1D(dir, "jetPtResGenThreeProng0Pi0", "jetPtRes (gen. three-prong, 0 #pi^{0})", 40, 0., 2.);
  histogramJetPtResGenThreeProng1Pi0_ = book1D(dir, "jetPtResGenThreeProng1Pi0", "jetPtRes (gen. three-prong, 1 #pi^{0})", 40, 0., 2.);

  histogramRecVsGenTauDecayMode_ = book2D(dir, "recVsGenTauDecayMode", "rec. vs. gen. Tau decay mode", 20, -0.5, 19.5, 20, -0.5, 19.5);
  setAxisLabelsGenTauDecayMode(histogramRecVsGenTauDecayMode_->GetXaxis());
  setAxisLabelsRecTauDecayMode(histogramRecVsGenTauDecayMode_->GetYaxis());
}

void TauPtResHistManager::fillHistograms(
       const pat::Tau& patTau, const reco::GenParticleCollection& genParticles, const reco::Vertex& vertex, 
       const reco::tau::RecoTauQualityCuts& qualityCuts, double weight)
{
  std::string genTauDecayMode = getGenTauDecayMode(patTau, genParticles);
  if ( !(genTauDecayMode == "oneProng0Pi0"   ||
	 genTauDecayMode == "oneProng1Pi0"   ||
	 genTauDecayMode == "oneProng2Pi0"   ||
	 genTauDecayMode == "threeProng0Pi0" ||
	 genTauDecayMode == "threeProng1Pi0") ) return;

  double genVisPt;
  if ( patTau.genJet() ) genVisPt = patTau.genJet()->pt();
  else return;
  if ( !(genVisPt > 0.) ) return;

  double recTauPt = patTau.pt();
  double tauPtRes = recTauPt/genVisPt;
  histogramTauPtRes_->Fill(tauPtRes, weight);
  if ( genTauDecayMode == "oneProng0Pi0"   ) histogramTauPtResGenOneProng0Pi0_->Fill(tauPtRes, weight);
  if ( genTauDecayMode == "oneProng1Pi0"   ) histogramTauPtResGenOneProng1Pi0_->Fill(tauPtRes, weight);
  if ( genTauDecayMode == "oneProng2Pi0"   ) histogramTauPtResGenOneProng2Pi0_->Fill(tauPtRes, weight);
  if ( genTauDecayMode == "threeProng0Pi0" ) histogramTauPtResGenThreeProng0Pi0_->Fill(tauPtRes, weight);
  if ( genTauDecayMode == "threeProng1Pi0" ) histogramTauPtResGenThreeProng1Pi0_->Fill(tauPtRes, weight);

  if ( tauPtRes < 0.5 ) {
    double eta = genVisPt = patTau.genJet()->eta();
    double phi = genVisPt = patTau.genJet()->phi();
    if ( genTauDecayMode == "oneProng0Pi0"   ) histogramEtaPhiForTauPtResLt05GenOneProng0Pi0_->Fill(eta, phi, weight);
    if ( genTauDecayMode == "oneProng1Pi0"   ) histogramEtaPhiForTauPtResLt05GenOneProng1Pi0_->Fill(eta, phi, weight);
    if ( genTauDecayMode == "oneProng2Pi0"   ) histogramEtaPhiForTauPtResLt05GenOneProng2Pi0_->Fill(eta, phi, weight);
    if ( genTauDecayMode == "threeProng0Pi0" ) histogramEtaPhiForTauPtResLt05GenThreeProng0Pi0_->Fill(eta, phi, weight);
    if ( genTauDecayMode == "threeProng1Pi0" ) histogramEtaPhiForTauPtResLt05GenThreeProng1Pi0_->Fill(eta, phi, weight);
  }

  histogramsTauPtResManCorrLev1_.fillHistograms(patTau, genTauDecayMode, genVisPt, vertex, qualityCuts, weight);
  histogramsTauPtResManCorrLev2_.fillHistograms(patTau, genTauDecayMode, genVisPt, vertex, qualityCuts, weight);
  histogramsTauPtResManCorrLev12_.fillHistograms(patTau, genTauDecayMode, genVisPt, vertex, qualityCuts, weight);
  histogramsTauPtResManCorrLev123_.fillHistograms(patTau, genTauDecayMode, genVisPt, vertex, qualityCuts, weight);
  histogramsTauPtResManCorrLev14_.fillHistograms(patTau, genTauDecayMode, genVisPt, vertex, qualityCuts, weight);

  double recJetPt = patTau.p4Jet().pt();
  double jetPtRes = recJetPt/genVisPt;
  histogramJetPtRes_->Fill(jetPtRes, weight);
  if ( genTauDecayMode == "oneProng0Pi0"   ) histogramJetPtResGenOneProng0Pi0_->Fill(jetPtRes, weight);
  if ( genTauDecayMode == "oneProng1Pi0"   ) histogramJetPtResGenOneProng1Pi0_->Fill(jetPtRes, weight);
  if ( genTauDecayMode == "oneProng2Pi0"   ) histogramJetPtResGenOneProng2Pi0_->Fill(jetPtRes, weight);
  if ( genTauDecayMode == "threeProng0Pi0" ) histogramJetPtResGenThreeProng0Pi0_->Fill(jetPtRes, weight);
  if ( genTauDecayMode == "threeProng1Pi0" ) histogramJetPtResGenThreeProng1Pi0_->Fill(jetPtRes, weight);
    
  int recTauDecayMode = patTau.decayMode();
  histogramRecVsGenTauDecayMode_->Fill(genTauDecayMode.data(), recTauDecayMode, weight);
}

TH1* TauPtResHistManager::book1D(TFileDirectory& dir,
				 const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, min, max);
  return retVal;
}
 
TH2* TauPtResHistManager::book2D(TFileDirectory& dir,
				 const std::string& distribution, const std::string& title, 
				 int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)
{
  TH2* retVal = dir.make<TH2D>(getHistogramName(distribution).data(), title.data(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  return retVal;
}
 
std::string TauPtResHistManager::getHistogramName(const std::string& distribution)
{
  std::string retVal = distribution;
  return retVal;
}
