#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"

#include <TMath.h>

TauIdEffHistManager::TauIdEffHistManager(const edm::ParameterSet& cfg)
{
  process_              = cfg.getParameter<std::string>("process");
  region_               = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_   = cfg.getParameter<std::string>("tauIdDiscriminator");
  label_                = cfg.getParameter<std::string>("label");
  svFitMassHypothesis_  = cfg.exists("svFitMassHypothesis") ?
    cfg.getParameter<std::string>("svFitMassHypothesis") : "";
}

TauIdEffHistManager::~TauIdEffHistManager()
{
// nothing to be done yet...
}

void TauIdEffHistManager::bookHistograms(TFileDirectory& dir)
{
  histogramMuonPt_          = book1D(dir, "muonPt",             "P_{T}^{#mu}",                           40,          0. ,         100.);
  histogramMuonEta_         = book1D(dir, "muonEta",            "#eta_{#mu}",                            50,         -2.5,         +2.5);
  histogramMuonPhi_         = book1D(dir, "muonPhi",            "#phi_{#mu}",                            36, -TMath::Pi(), +TMath::Pi());
  
  histogramTauPt_           = book1D(dir, "tauJetPt",           "P_{T}^{#tau}",                          40,          0. ,         100.);
  histogramTauEta_          = book1D(dir, "tauJetEta",          "#eta_{#tau}",                           50,         -2.5,         +2.5);
  histogramTauPhi_          = book1D(dir, "tauJetPhi",          "#phi_{#tau}",                           36, -TMath::Pi(), +TMath::Pi());
  histogramTauNumTracks_    = book1D(dir, "tauJetNumTracks",    "Num. Tracks #tau-Jet",                  25,         -0.5,         24.5);
  histogramTauNumSelTracks_ = book1D(dir, "tauJetNumSelTracks", "Num. selected Tracks #tau-Jet",         25,         -0.5,         24.5);

  histogramVisMass_         = book1D(dir, "diTauVisMass",       "M_{vis}(#mu + #tau_{had})",             36,         20.0,        200.0);
  if ( svFitMassHypothesis_ != "" )
    histogramSVfitMass_     = book1D(dir, "diTauSVfitMass",     "SVfit Mass",                            42,         40.0,        250.0);
  histogramMt_              = book1D(dir, "diTauMt",            "M_{T}(#mu + MET)",                      24,          0.0,        120.0);
  histogramPzetaDiff_       = book1D(dir, "diTauPzetaDiff",     "P_{#zeta} - 1.5 #cdot P_{#zeta}^{vis}", 24,        -80.0,        +40.0);

  histogramMEt_             = book1D(dir, "met",                "E_{T}^{miss}",                          20,          0.0,        100.0);
  histogramSumEt_           = book1D(dir, "sumEt",              "#Sigma E_{T}^{PF}",                     50,          0.,         500.0);
  histogramNumVertices_     = book1D(dir, "numVertices",        "Num. Vertices",                         20,         -0.5,         19.5);
}

void TauIdEffHistManager::fillHistograms(const PATMuTauPair& muTauPair, size_t numVertices, double weight)
{
  histogramMuonPt_->Fill(muTauPair.leg1()->pt(), weight);
  histogramMuonEta_->Fill(muTauPair.leg1()->eta(), weight);
  histogramMuonPhi_->Fill(muTauPair.leg1()->phi(), weight);
  
  histogramTauPt_->Fill(muTauPair.leg2()->pt(), weight);
  histogramTauEta_->Fill(muTauPair.leg2()->eta(), weight);
  histogramTauPhi_->Fill(muTauPair.leg2()->phi(), weight);
  histogramTauNumTracks_->Fill(muTauPair.leg2()->userFloat("numTracks"), weight);
  histogramTauNumSelTracks_->Fill(muTauPair.leg2()->userFloat("numSelTracks"), weight);
  
  histogramVisMass_->Fill((muTauPair.leg1()->p4() + muTauPair.leg2()->p4()).mass(), weight); 
  if ( svFitMassHypothesis_ != "" ) {
    int errorFlag;
    const NSVfitResonanceHypothesisSummary* svFitSolution = muTauPair.nSVfitSolution(svFitMassHypothesis_, &errorFlag);
    if ( svFitSolution ) histogramSVfitMass_->Fill(svFitSolution->mass(), weight); 
  }
  histogramMt_->Fill(muTauPair.mt1MET(), weight);
  histogramPzetaDiff_->Fill(muTauPair.pZeta() - 1.5*muTauPair.pZetaVis(), weight);

  histogramMEt_->Fill(muTauPair.met()->pt(), weight);
  histogramSumEt_->Fill(muTauPair.met()->sumEt(), weight);
  histogramNumVertices_->Fill(numVertices, weight);
}

void TauIdEffHistManager::scaleHistograms(double factor)
{
  for ( std::vector<TH1*>::iterator histogram = histograms_.begin();
	histogram != histograms_.end(); ++histogram ) {
    if ( !(*histogram)->GetSumw2N() ) (*histogram)->Sumw2(); // CV: compute "proper" errors before scaling histogram
    (*histogram)->Scale(factor);
  }
}

TH1* TauIdEffHistManager::book1D(TFileDirectory& dir,
				 const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, min, max);
  histograms_.push_back(retVal);
  return retVal;
}
 
TH1* TauIdEffHistManager::book1D(TFileDirectory& dir,
				 const std::string& distribution, const std::string& title, int numBins, float* binning)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, binning);
  histograms_.push_back(retVal);
  return retVal;
}
 
std::string TauIdEffHistManager::getHistogramName(const std::string& distribution)
{
  std::string retVal = std::string(process_).append("_").append(region_).append("_").append(distribution);
  retVal.append("_").append(tauIdDiscriminator_).append("_").append(label_);
  return retVal;
}
