#include "TauAnalysis/TauIdEfficiency/interface/TauFakeRateHistManager.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

TauFakeRateHistManager::TauFakeRateHistManager(const edm::ParameterSet& cfg)
{
  process_            = cfg.getParameter<std::string>("process");
  region_             = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_ = cfg.getParameter<std::string>("tauIdDiscriminator");
}

TauFakeRateHistManager::~TauFakeRateHistManager()
{
// nothing to be done yet...
}

void TauFakeRateHistManager::bookHistograms(TFileDirectory& dir)
{
  histogramJetPt_       = book1D(dir, "jetPt",       "P_{T}^{jet}",       20,          0. ,         100.);
  histogramJetEta_      = book1D(dir, "jetEta",      "#eta_{jet}",        25,         -2.5,         +2.5);
  histogramJetPhi_      = book1D(dir, "jetPhi",      "#phi_{jet}",        36, -TMath::Pi(), +TMath::Pi());
  
  histogramTauPt_       = book1D(dir, "tauJetPt",    "P_{T}^{#tau}",      20,          0. ,         100.);
  histogramTauEta_      = book1D(dir, "tauJetEta",   "#eta_{#tau}",       25,         -2.5,         +2.5);
  histogramTauPhi_      = book1D(dir, "tauJetPhi",   "#phi_{#tau}",       36, -TMath::Pi(), +TMath::Pi());

  histogramSumEt_       = book1D(dir, "sumEt",       "#Sigma E_{T}^{PF}",  8,        100.,          500.0);
  histogramNumVertices_ = book1D(dir, "numVertices", "Num. Vertices",     20,         -0.5,         19.5);
}

void TauFakeRateHistManager::fillHistograms(const pat::Tau& tauJetCand, size_t numVertices, double sumEt, double weight)
{
  std::cout << "<TauFakeRateHistManager::fillHistograms>:" << std::endl;
  std::cout << " filling " << histogramJetPt_->GetName() << ": jetPt = " << tauJetCand.p4Jet().pt() << "," 
	    << " jetEta = " << tauJetCand.p4Jet().eta() << ", jetPhi = " << tauJetCand.p4Jet().phi() << std::endl;
  std::string histogramJetPtName = histogramJetPt_->GetName();
  if ( tauJetCand.p4Jet().pt() > 20. && tauJetCand.p4Jet().pt() < 25. && 
       histogramJetPtName.find("tauDiscrHPSloose")       != std::string::npos &&
       histogramJetPtName.find("tauDiscrHPSlooseDBcorr") == std::string::npos ) {
    static unsigned numJets_passed = 0;
    if ( histogramJetPtName.find("_P_") != std::string::npos ) {
      ++numJets_passed;
      std::cout << "decayModeFinding = " << tauJetCand.tauID("decayModeFinding") << std::endl;
      std::cout << "byLooseIsolation = " << tauJetCand.tauID("byLooseIsolation") << std::endl;
    }
    static unsigned numJets_failed = 0;
    if ( histogramJetPtName.find("_F_") != std::string::npos ) ++numJets_failed;
    static unsigned numJets_all    = 0;
    if ( histogramJetPtName.find("_A_") != std::string::npos ) ++numJets_all;
    std::cout << "numJets_passed = " << numJets_passed << std::endl;
    std::cout << "numJets_failed = " << numJets_failed << std::endl;
    std::cout << "numJets_all = " << numJets_all << std::endl;
  }
  histogramJetPt_->Fill(tauJetCand.p4Jet().pt(), weight);
  histogramJetEta_->Fill(tauJetCand.p4Jet().eta(), weight);
  histogramJetPhi_->Fill(tauJetCand.p4Jet().phi(), weight);
  
  std::cout << " filling " << histogramTauPt_->GetName() << ": tauPt = " << tauJetCand.pt() << "," 
	    << " tauEta = " << tauJetCand.eta() << ", tauPhi = " << tauJetCand.phi() << std::endl;
  histogramTauPt_->Fill(tauJetCand.pt(), weight);
  histogramTauEta_->Fill(tauJetCand.eta(), weight);
  histogramTauPhi_->Fill(tauJetCand.phi(), weight);
  
  histogramSumEt_->Fill(sumEt, weight);
  histogramNumVertices_->Fill(numVertices, weight);
}

void TauFakeRateHistManager::scaleHistograms(double factor)
{
  for ( std::vector<TH1*>::iterator histogram = histograms_.begin();
	histogram != histograms_.end(); ++histogram ) {
    if ( !(*histogram)->GetSumw2N() ) (*histogram)->Sumw2(); // CV: compute "proper" errors before scaling histogram
    (*histogram)->Scale(factor);
  }
}

TH1* TauFakeRateHistManager::book1D(TFileDirectory& dir,
				    const std::string& distribution, const std::string& title, int numBins, double min, double max)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, min, max);;
  histograms_.push_back(retVal);
  return retVal;
}
 
TH1* TauFakeRateHistManager::book1D(TFileDirectory& dir,
				    const std::string& distribution, const std::string& title, int numBins, float* binning)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, binning);
  histograms_.push_back(retVal);
  return retVal;
}
 
std::string TauFakeRateHistManager::getHistogramName(const std::string& distribution)
{
  std::string retVal = std::string(process_).append("_").append(region_).append("_").append(distribution);
  retVal.append("_").append(tauIdDiscriminator_);
  return retVal;
}
