#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

TauIdEffHistManager::TauIdEffHistManager(const edm::ParameterSet& cfg)
{
  process_            = cfg.getParameter<std::string>("process");
  region_             = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_ = cfg.getParameter<std::string>("tauIdDiscriminator");
  label_              = cfg.getParameter<std::string>("label");
}

TauIdEffHistManager::~TauIdEffHistManager()
{
// nothing to be done yet...
}

void TauIdEffHistManager::bookHistograms(fwlite::TFileService& fs)
{
  histogramMuonPt_       = book1D(fs, "muonPt",          "muonPt",                    40,          0. ,         100.);
  histogramMuonEta_      = book1D(fs, "muonEta",         "muonEta",                   50,         -2.5,         +2.5);
  histogramMuonPhi_      = book1D(fs, "muonPhi",         "muonPhi",                   36, -TMath::Pi(), +TMath::Pi());
  
  histogramTauPt_        = book1D(fs, "tauJetPt",        "tauPt",                     40,          0. ,         100.);
  histogramTauEta_       = book1D(fs, "tauJetEta",       "tauEta",                    50,         -2.5,         +2.5);
  histogramTauPhi_       = book1D(fs, "tauJetPhi",       "tauPhi",                    36, -TMath::Pi(), +TMath::Pi());
  histogramTauNumTracks_ = book1D(fs, "tauJetNumTracks", "tauNumTracks",              25,         -0.5,         24.5);
  
  histogramVisMass_      = book1D(fs, "diTauVisMass",    "M_{vis}(#mu + #tau_{had})", 36,         20.0,        200.0);
  histogramMt_           = book1D(fs, "diTauMt",         "M_{T}(#mu + MET)",          16,          0.0,         80.0);
}

TH1* TauIdEffHistManager::book1D(fwlite::TFileService& fs,
				 const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  return fs.make<TH1F>(getHistogramName(distribution).data(), title.data(), numBins, min, max);
}
 
TH1* TauIdEffHistManager::book1D(fwlite::TFileService& fs,
				 const std::string& distribution, const std::string& title, int numBins, float* binning)
{
  return fs.make<TH1F>(getHistogramName(distribution).data(), title.data(), numBins, binning);
}
 
std::string TauIdEffHistManager::getHistogramName(const std::string& distribution)
{
  std::string retVal = std::string(process_).append("_").append(region_).append("_").append(distribution);
  retVal.append("_").append(tauIdDiscriminator_).append("_").append(label_);
  return retVal;
}
