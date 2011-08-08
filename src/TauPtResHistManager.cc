#include "TauAnalysis/TauIdEfficiency/interface/TauPtResHistManager.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"

#include <TMath.h>

TauPtResHistManager::TauPtResHistManager(const edm::ParameterSet& cfg)
{
// nothing to be done yet...
}

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

  histogramTauPtResExclK0s_                  
    = book1D(dir, "tauPtResExclK0s",                  "tauPtResExclK0s",                               40, 0., 2.);
  histogramTauPtResExclK0sGenOneProng0Pi0_   
    = book1D(dir, "tauPtResExclK0sGenOneProng0Pi0",   "tauPtResExclK0s (gen. one-prong, 0 #pi^{0})",   40, 0., 2.);
  histogramTauPtResExclK0sGenOneProng1Pi0_   
    = book1D(dir, "tauPtResExclK0sGenOneProng1Pi0",   "tauPtResExclK0s (gen. one-prong, 1 #pi^{0})",   40, 0., 2.);
  histogramTauPtResExclK0sGenOneProng2Pi0_   
    = book1D(dir, "tauPtResExclK0sGenOneProng2Pi0",   "tauPtResExclK0s (gen. one-prong, 2 #pi^{0})",   40, 0., 2.);
  histogramTauPtResExclK0sGenThreeProng0Pi0_ 
    = book1D(dir, "tauPtResExclK0sGenThreeProng0Pi0", "tauPtResExclK0s (gen. three-prong, 0 #pi^{0})", 40, 0., 2.);
  histogramTauPtResExclK0sGenThreeProng1Pi0_ 
    = book1D(dir, "tauPtResExclK0sGenThreeProng1Pi0", "tauPtResExclK0s (gen. three-prong, 1 #pi^{0})", 40, 0., 2.);

  histogramJetPtRes_                  = book1D(dir, "jetPtRes",                  "jetPtRes",                               40, 0., 2.);
  histogramJetPtResGenOneProng0Pi0_   = book1D(dir, "jetPtResGenOneProng0Pi0",   "jetPtRes (gen. one-prong, 0 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenOneProng1Pi0_   = book1D(dir, "jetPtResGenOneProng1Pi0",   "jetPtRes (gen. one-prong, 1 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenOneProng2Pi0_   = book1D(dir, "jetPtResGenOneProng2Pi0",   "jetPtRes (gen. one-prong, 2 #pi^{0})",   40, 0., 2.);
  histogramJetPtResGenThreeProng0Pi0_ = book1D(dir, "jetPtResGenThreeProng0Pi0", "jetPtRes (gen. three-prong, 0 #pi^{0})", 40, 0., 2.);
  histogramJetPtResGenThreeProng1Pi0_ = book1D(dir, "jetPtResGenThreeProng1Pi0", "jetPtRes (gen. three-prong, 1 #pi^{0})", 40, 0., 2.);

  histogramRecVsGenTauDecayMode_ = book2D(dir, "recVsGenTauDecayMode", "rec. vs. gen. Tau decay mode", 20, -0.5, 19.5, 20, -0.5, 19.5);
  setAxisLabelsGenTauDecayMode(histogramRecVsGenTauDecayMode_->GetXaxis());
  setAxisLabelsRecTauDecayMode(histogramRecVsGenTauDecayMode_->GetYaxis());

  histogramSumPFNeutralHadronPtVsTauPtRes_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtRes", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes", 40, 0., 2., 40, 0., 2.);
  histogramSumPFNeutralHadronPtVsTauPtResGenOneProng0Pi0_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtResGenOneProng0Pi0", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes (gen. one-prong, 0 #pi^{0})", 40, 0., 2., 40, 0., 2.);
  histogramSumPFNeutralHadronPtVsTauPtResGenOneProng1Pi0_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtResGenOneProng1Pi0", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes (gen. one-prong, 1 #pi^{0})", 40, 0., 2., 40, 0., 2.);
  histogramSumPFNeutralHadronPtVsTauPtResGenOneProng2Pi0_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtResGenOneProng2Pi0", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes (gen. one-prong, 2 #pi^{0})", 40, 0., 2., 40, 0., 2.);
  histogramSumPFNeutralHadronPtVsTauPtResGenThreeProng0Pi0_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtResGenThreeProng0Pi0", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes (gen. three-prong, 0 #pi^{0})", 40, 0., 2., 40, 0., 2.);
  histogramSumPFNeutralHadronPtVsTauPtResGenThreeProng1Pi0_ = 
    book2D(dir, "sumPFNeutralHadronPtVsTauPtResGenThreeProng1Pi0", 
	   "#Sigma PFNeutralHadron P_{T} vs. tauPtRes (gen. three-prong, 1 #pi^{0})", 40, 0., 2., 40, 0., 2.);
}

void TauPtResHistManager::fillHistograms(const pat::Tau& patTau, const reco::GenParticleCollection& genParticles, double weight)
{
  std::string genTauDecayMode;
  const reco::GenParticle* genTau = findGenParticle(patTau.p4(), genParticles);
  if      ( genTau          ) genTauDecayMode = getGenTauDecayMode(genTau);
  else if ( patTau.genJet() ) genTauDecayMode = JetMCTagUtils::genTauDecayMode(*patTau.genJet());
  else return;

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
  
  bool isGenK0s = false;
  std::vector<const reco::GenParticle*> genTauJetConstituents = patTau.genJet()->getGenConstituents();
  for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
	genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {
    if ( (*genTauJetConstituent)->pdgId() == 310 ) isGenK0s = true;
  }
  if ( !isGenK0s ) {
    histogramTauPtResExclK0s_->Fill(tauPtRes, weight);
    if ( genTauDecayMode == "oneProng0Pi0"   ) histogramTauPtResExclK0sGenOneProng0Pi0_->Fill(tauPtRes, weight);
    if ( genTauDecayMode == "oneProng1Pi0"   ) histogramTauPtResExclK0sGenOneProng1Pi0_->Fill(tauPtRes, weight);
    if ( genTauDecayMode == "oneProng2Pi0"   ) histogramTauPtResExclK0sGenOneProng2Pi0_->Fill(tauPtRes, weight);
    if ( genTauDecayMode == "threeProng0Pi0" ) histogramTauPtResExclK0sGenThreeProng0Pi0_->Fill(tauPtRes, weight);
    if ( genTauDecayMode == "threeProng1Pi0" ) histogramTauPtResExclK0sGenThreeProng1Pi0_->Fill(tauPtRes, weight);
  }

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

  double pfNeutralHadronPtSum = 0.;
  std::vector<reco::PFCandidatePtr> pfJetConstituents = patTau.pfJetRef()->getPFConstituents();
  for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
    if ( (*pfJetConstituent)->particleId() == reco::PFCandidate::h0 ) pfNeutralHadronPtSum += (*pfJetConstituent)->pt();
  }
  pfNeutralHadronPtSum /= genVisPt;
  histogramSumPFNeutralHadronPtVsTauPtRes_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
  if ( genTauDecayMode == "oneProng0Pi0"   ) 
    histogramSumPFNeutralHadronPtVsTauPtResGenOneProng0Pi0_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
  if ( genTauDecayMode == "oneProng1Pi0"   ) 
    histogramSumPFNeutralHadronPtVsTauPtResGenOneProng1Pi0_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
  if ( genTauDecayMode == "oneProng2Pi0"   ) 
    histogramSumPFNeutralHadronPtVsTauPtResGenOneProng2Pi0_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
  if ( genTauDecayMode == "threeProng0Pi0" ) 
    histogramSumPFNeutralHadronPtVsTauPtResGenThreeProng0Pi0_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
  if ( genTauDecayMode == "threeProng1Pi0" ) 
    histogramSumPFNeutralHadronPtVsTauPtResGenThreeProng1Pi0_->Fill(tauPtRes, pfNeutralHadronPtSum, weight);
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
