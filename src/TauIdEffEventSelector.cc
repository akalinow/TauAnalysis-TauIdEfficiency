#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// define flag for Mt && Pzeta cut (not appplied, Mt && Pzeta cut passed, Mt || Pzeta cut failed) 
// and tau id. discriminators      (no tau id. discriminators applied, all discriminators passed, at least one discriminator failed)
enum { kNotApplied, kSignalLike, kBackgroundLike };

TauIdEffEventSelector::TauIdEffEventSelector(const edm::ParameterSet& cfg)
{
  tauIdDiscriminators_ = cfg.getParameter<vstring>("tauIdDiscriminators");

//--- define default cuts for ABCD regions
  muonPtMin_              =  20.0; 
  muonPtMax_              =  +1.e+3; 
  muonEtaMin_             =  -2.1;
  muonEtaMax_             =  +2.1; 
  muonRelIsoMin_          =  -1.e+3;
  muonRelIsoMax_          =   0.30;
  tauPtMin_               =  20.0;   // CV: cut is actually applied on PFJet Pt, not PFTau Pt
  tauPtMax_               =  +1.e+3; 
  tauEtaMin_              =  -2.1;
  tauEtaMax_              =  +2.1;   
  tauAbsIsoMin_           =  -1.e+3;
  tauAbsIsoMax_           =   2.5;
  muTauPairAbsDzMax_      =  +1.e+3;
  muTauPairChargeMin_     =  -1.e+3;
  muTauPairChargeMax_     =  +1.e+3; 
  MtMin_                  =  -1.e+3;
  MtMax_                  =  40.0;
  PzetaMin_               = -20.0;
  PzetaMax_               =  +1.e+3;
  MtAndPzetaCut_          = kNotApplied;
  tauIdDiscriminatorMin_  =   0.5;
  tauIdDiscriminatorMax_  =  +1.e+3;
  tauIdDiscriminatorCut_  = kNotApplied;

//--- additional cuts to make sure there are no events in underflow/overflow bins
//    of visible and transverse mass distributions
//   (difference in event yields may cause problem with simultaneous fit of visMass and Mt distributions)
  visMassCutoffMin_       =  20.0;
  visMassCutoffMax_       = 200.0;
  MtCutoffMin_            =  -1.e+3;
  MtCutoffMax_            =  80.0;

  region_ = cfg.getParameter<std::string>("region");
  if        ( region_           == "ABCD"            ) {
    // nothing to be done yet...
    // (simply use default cuts)
  } else if ( region_.find("A") != std::string::npos ) {
    muonRelIsoMin_        =   0.10;
    muTauPairChargeMax_   =  -0.5;
  } else if ( region_.find("B") != std::string::npos ) {
    muonRelIsoMin_        =   0.10;
    muTauPairChargeMin_   =  +0.5;
  } else if ( region_.find("C") != std::string::npos ) {
    muonRelIsoMax_        =   0.10;
    muTauPairChargeMax_   =  -0.5;
  } else if ( region_.find("D") != std::string::npos ) {
    muonRelIsoMax_        =   0.10;
    muTauPairChargeMin_   =  +0.5;
  } else {
    throw cms::Exception("TauIdEffEventSelector") 
      << "Invalid region = " << region_ << " !!\n";
  }

  if      ( region_.find("1") != std::string::npos ) MtAndPzetaCut_ = kSignalLike;
  else if ( region_.find("2") != std::string::npos ) MtAndPzetaCut_ = kBackgroundLike;

  if      ( region_.find("p") != std::string::npos ) tauIdDiscriminatorCut_ = kSignalLike;
  else if ( region_.find("f") != std::string::npos ) tauIdDiscriminatorCut_ = kBackgroundLike;
}

TauIdEffEventSelector::~TauIdEffEventSelector()
{
// nothing to be done yet...
}

bool TauIdEffEventSelector::operator()(const PATMuTauPair& muTauPair, pat::strbitset& result)
{
  double muonPt          = muTauPair.leg1()->pt();
  double muonEta         = muTauPair.leg1()->eta();
  double muonIso         = muTauPair.leg1()->userFloat("pfLooseIsoPt04");
  double tauPt           = muTauPair.leg2()->pt();
  double tauEta          = muTauPair.leg2()->eta();
  double tauIso          = muTauPair.leg2()->userFloat("preselLoosePFIsoPt");
  double muTauPairAbsDz  = TMath::Abs(muTauPair.leg1()->vertex().z() - muTauPair.leg2()->vertex().z());
  double muTauPairCharge = muTauPair.leg1()->charge() + muTauPair.leg2()->userFloat("leadTrackCharge");
  double visMass         = (muTauPair.leg1()->p4() + muTauPair.leg2()->p4()).mass();
  double Mt              = muTauPair.mt1MET();
  double Pzeta           = muTauPair.pZeta() - 1.5*muTauPair.pZetaVis();
  
  if ( muonPt          >  muonPtMin_             && muonPt          >  muonPtMax_             &&
       muonEta         >  muonEtaMin_            && muonEta         <  muonEtaMax_            &&
       muonIso         > (muonRelIsoMin_*muonPt) && muonIso         < (muonRelIsoMax_*muonPt) && 
       tauPt           >  tauPtMin_              && tauPt           >  tauPtMax_              &&
       tauEta          >  tauEtaMin_             && tauEta          <  tauEtaMax_             &&
       tauIso          >  tauAbsIsoMin_          && tauIso          <  tauAbsIsoMax_          && 
       muTauPairAbsDz  <  muTauPairAbsDzMax_     &&
       muTauPairCharge >  muTauPairChargeMin_    && muTauPairCharge <  muTauPairChargeMin_    && 
       visMass         >  visMassCutoffMin_      && visMass         <  visMassCutoffMax_      &&
       Mt              >  MtCutoffMin_           && Mt              <  MtCutoffMax_           ) {
    bool MtAndPzetaCut_passed = (Mt > MtMin_ && Mt > MtMax_ && Pzeta > PzetaMin_ && Pzeta < PzetaMax_);

    bool tauIdDiscriminators_passed = true;
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      double tauIdDiscriminator_value = muTauPair.leg2()->tauID(*tauIdDiscriminator);
      if ( !(tauIdDiscriminator_value > tauIdDiscriminatorMin_  && 
	     tauIdDiscriminator_value < tauIdDiscriminatorMax_) ) tauIdDiscriminators_passed = false;
    }

    if ( ( MtAndPzetaCut_ == kNotApplied                                             ||
	  (MtAndPzetaCut_ == kSignalLike             &&  MtAndPzetaCut_passed)       ||
	  (MtAndPzetaCut_ == kBackgroundLike         && !MtAndPzetaCut_passed)       ) &&
	 ( tauIdDiscriminatorCut_ == kNotApplied                                     ||
	  (tauIdDiscriminatorCut_ == kSignalLike     &&  tauIdDiscriminators_passed) ||
	  (tauIdDiscriminatorCut_ == kBackgroundLike && !tauIdDiscriminators_passed) ) ) return true;
  }

  return false;
}
