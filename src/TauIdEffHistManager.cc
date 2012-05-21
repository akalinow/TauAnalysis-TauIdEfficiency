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
  triggerBits_          = cfg.getParameter<vstring>("triggerBits");
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

  histogramPFMEt_           = book1D(dir, "pfMEt",              "pf-E_{T}^{miss}",                       20,          0.0,        100.0);
  histogramPFSumEt_         = book1D(dir, "pfSumEt",            "#Sigma E_{T}^{PF}",                     50,          0.,         500.0);
  histogramCaloMEt_         = book1D(dir, "caloMEt",            "calo-E_{T}^{miss}",                     20,          0.0,        100.0);
  histogramCaloSumEt_       = book1D(dir, "caloSumEt",          "#Sigma E_{T}^{calo}",                   50,          0.,         500.0);

  histogramNumVertices_     = book1D(dir, "numVertices",        "Num. Vertices",                         35,         -0.5,         34.5);

  histogramLogEvtWeight_    = book1D(dir, "logEvtWeight",       "log(Event weight)",                    101,        -5.05,        +5.05);

  for ( vstring::const_iterator triggerBit = triggerBits_.begin();
	triggerBit != triggerBits_.end(); ++triggerBit ) {
    // CV: skip creating histogram in case it already exists
    //    (created for different version of HLT path)
    if ( histogramNumCaloMEt_.find(*triggerBit) != histogramNumCaloMEt_.end() ) continue;
    std::string histogramName  = std::string("numCaloMEt").append("_").append(*triggerBit);
    std::string histogramTitle = std::string(*triggerBit).append(" vs. calo-E_{T}^{miss}");
    histogramNumCaloMEt_[*triggerBit] = book1D(dir, histogramName, histogramTitle, 100, 0., 100.);
    //std::cout << "histogramNumCaloMEt[" << (*triggerBit) << "] = " 
    //          << histogramNumCaloMEt_[*triggerBit] << std::endl;
  }
  histogramDenomCaloMEt_ = book1D(dir, "denomCaloMEt", "calo-E_{T}^{miss} denominator", 100, 0., 100.);
  
  histogramEventCounter_    = book1D(dir, "EventCounter",       "Event Counter",                          1,         -0.5,         +0.5);
}

void TauIdEffHistManager::fillHistograms(const PATMuTauPair& muTauPair, const pat::MET& caloMEt, 
					 size_t numVertices, const std::map<std::string, bool>& triggerBits_passed, double weight)
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

  histogramPFMEt_->Fill(muTauPair.met()->pt(), weight);
  histogramPFSumEt_->Fill(muTauPair.met()->sumEt(), weight);
  histogramCaloMEt_->Fill(caloMEt.pt(), weight);
  histogramCaloSumEt_->Fill(caloMEt.sumEt(), weight);

  histogramNumVertices_->Fill(numVertices, weight);

  if ( weight > 0. ) {
    double logWeight = TMath::Log(weight);
    if      ( logWeight < -5.0 ) logWeight = -5.0;
    else if ( logWeight > +5.0 ) logWeight = +5.0;
    histogramLogEvtWeight_->fill(logWeight);
  }
  
  for ( std::map<std::string, bool>::const_iterator triggerBit_passed = triggerBits_passed.begin();
	triggerBit_passed != triggerBits_passed.end(); ++triggerBit_passed ) {    
    //std::cout << "histogramNumCaloMEt[" << triggerBit_passed->first << "] = " 
    //          << histogramNumCaloMEt_[triggerBit_passed->first] << std::endl;
    std::cout << triggerBit_passed->first << ": " << triggerBit_passed->second << std::endl;
    assert(histogramNumCaloMEt_[triggerBit_passed->first]);
    if ( triggerBit_passed->second ) histogramNumCaloMEt_[triggerBit_passed->first]->Fill(caloMEt.pt(), weight);
  }
  histogramDenomCaloMEt_->Fill(caloMEt.pt(), weight);

  histogramEventCounter_->Fill(0., weight);
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
  if ( !retVal->GetSumw2N() ) retVal->Sumw2();
  histograms_.push_back(retVal);
  return retVal;
}
 
TH1* TauIdEffHistManager::book1D(TFileDirectory& dir,
				 const std::string& distribution, const std::string& title, int numBins, float* binning)
{
  TH1* retVal = dir.make<TH1D>(getHistogramName(distribution).data(), title.data(), numBins, binning);
  if ( !retVal->GetSumw2N() ) retVal->Sumw2();
  histograms_.push_back(retVal);
  return retVal;
}
 
std::string TauIdEffHistManager::getHistogramName(const std::string& distribution)
{
  std::string retVal = std::string(process_).append("_").append(region_).append("_").append(distribution);
  retVal.append("_").append(tauIdDiscriminator_).append("_").append(label_);
  return retVal;
}
