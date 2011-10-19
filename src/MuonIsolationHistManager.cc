#include "TauAnalysis/TauIdEfficiency/interface/MuonIsolationHistManager.h"

#include <TMath.h>

MuonIsolationHistManager::MuonIsolationHistManager(const edm::ParameterSet& cfg)
{
// nothing to be done yet...
}

MuonIsolationHistManager::~MuonIsolationHistManager()
{
// nothing to be done yet...
}

void MuonIsolationHistManager::bookHistograms(TFileDirectory& dir)
{
  const int numMuonEtaBins = 16;
  float muonEtaBinning[numMuonEtaBins + 1] = { -2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., +0.3, +0.6, +0.9, +1.2, +1.5, +1.8, +2.1, +2.4 };

  const int numMuonPtBins = 14;
  float muonPtBinning[numMuonPtBins + 1] = { 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 45., 50., 60., 80., 120., 200. };

  histogramMuonPt_          = book1D(dir, "muonPt",             "P_{T}^{#mu}",                numMuonPtBins,  muonPtBinning);
  histogramMuonEta_         = book1D(dir, "muonEta",            "#eta_{#mu}",                 numMuonEtaBins, muonEtaBinning);
  histogramMuonPtVsEta_     = book2D(dir, "muonPtVsEta",        "P_{T}^{#mu} vs. #eta_{#mu}", numMuonEtaBins, muonEtaBinning, numMuonPtBins, muonPtBinning);
  histogramMuonPhi_         = book1D(dir, "muonPhi",            "#phi_{#mu}",                            36, -TMath::Pi(), +TMath::Pi());
  histogramMuonCharge_      = book1D(dir, "muonCharge",         "#mu Charge",                            3,          -1.5,         +1.5);
  histogramMuonAbsIso_      = book1D(dir, "muonAbsIso",         "#mu abs. Isolation",                    100,        -0.01,       +10.);
  histogramMuonRelIso_      = book1D(dir, "muonRelIso",         "#mu rel. Isolation",                    251,        -0.005,       +2.505);
  
  histogramTauJetPt_        = book1D(dir, "tauJetPt",           "P_{T}^{#tau}",                          40,          0.,          100.);
  histogramTauJetEta_       = book1D(dir, "tauJetEta",          "#eta_{#tau}",                           25,         -2.5,         +2.5);
  histogramTauJetPhi_       = book1D(dir, "tauJetPhi",          "#phi_{#tau}",                           36, -TMath::Pi(), +TMath::Pi());

  histogramVisMass_         = book1D(dir, "diTauVisMass",       "M_{vis}(#mu + #tau_{had})",             36,         20.0,        200.0);
  histogramMt_              = book1D(dir, "diTauMt",            "M_{T}(#mu + MET)",                      24,          0.0,        120.0);
  histogramPzetaDiff_       = book1D(dir, "diTauPzetaDiff",     "P_{#zeta} - 1.5 #cdot P_{#zeta}^{vis}", 24,        -80.0,        +40.0);

  histogramMEt_             = book1D(dir, "met",                "E_{T}^{miss}",                          20,          0.0,        100.0);
  histogramSumEt_           = book1D(dir, "sumEt",              "#Sigma E_{T}^{PF}",                     50,          0.,         500.0);
  histogramNumVertices_     = book1D(dir, "numVertices",        "Num. Vertices",                         20,         -0.5,         19.5);
}

void MuonIsolationHistManager::fillHistograms(const PATMuTauPair& muTauPair, size_t numVertices, double muonIsoPtSum, double weight)
{
  histogramMuonPt_->Fill(muTauPair.leg1()->pt(), weight);
  histogramMuonEta_->Fill(muTauPair.leg1()->eta(), weight);
  histogramMuonPtVsEta_->Fill(muTauPair.leg1()->eta(), muTauPair.leg1()->pt(), weight);
  histogramMuonPhi_->Fill(muTauPair.leg1()->phi(), weight);
  histogramMuonCharge_->Fill(muTauPair.leg1()->charge(), weight);
  histogramMuonAbsIso_->Fill(muonIsoPtSum, weight);
  if ( muTauPair.leg1()->pt() > 0. ) histogramMuonRelIso_->Fill(muonIsoPtSum/muTauPair.leg1()->pt(), weight);
  
  histogramTauJetPt_->Fill(muTauPair.leg2()->pt(), weight);
  histogramTauJetEta_->Fill(muTauPair.leg2()->eta(), weight);
  histogramTauJetPhi_->Fill(muTauPair.leg2()->phi(), weight);
  
  histogramVisMass_->Fill((muTauPair.leg1()->p4() + muTauPair.leg2()->p4()).mass(), weight); 
  histogramMt_->Fill(muTauPair.mt1MET(), weight);
  histogramPzetaDiff_->Fill(muTauPair.pZeta() - 1.5*muTauPair.pZetaVis(), weight);

  histogramMEt_->Fill(muTauPair.met()->pt(), weight);
  histogramSumEt_->Fill(muTauPair.met()->sumEt(), weight);
  histogramNumVertices_->Fill(numVertices, weight);
}

TH1* MuonIsolationHistManager::book1D(TFileDirectory& dir,
				      const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, min, max);
  histograms_.push_back(retVal);
  return retVal;
}
 
TH1* MuonIsolationHistManager::book1D(TFileDirectory& dir,
				      const std::string& distribution, const std::string& title, int numBins, float* binning)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, binning);
  histograms_.push_back(retVal);
  return retVal;
}

TH2* MuonIsolationHistManager::book2D(TFileDirectory& dir,
				      const std::string& distribution, const std::string& title, int numBinsX, float* xBinning, int numBinsY, float* yBinning)
{
  TH2* retVal = dir.make<TH2D>(getHistogramName(distribution).data(), title.data(), numBinsX, xBinning, numBinsY, yBinning);
  histograms_.push_back(retVal);
  return retVal;
}
 
std::string MuonIsolationHistManager::getHistogramName(const std::string& distribution)
{
  return distribution;
}
