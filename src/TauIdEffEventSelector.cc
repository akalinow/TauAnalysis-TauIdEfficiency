#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

// define flag for Mt && Pzeta cut (not appplied, Mt && Pzeta cut passed, Mt || Pzeta cut failed) 
// and tau id. discriminators      (no tau id. discriminators applied, all discriminators passed, at least one discriminator failed)
enum { kNotApplied, kSignalLike, kBackgroundLike, kWplusJetBackgroundLike };

// define flag indicating whether to take charge of tau-jet candidate 
// from "leading track" or from all "signal" charged hadrons
enum { kLeadTrackCharge, kSignalChargedHadronSum };

TauIdEffEventSelector::TauIdEffEventSelector(const edm::ParameterSet& cfg)
{
  //std::cout << "<TauIdEffEventSelector::TauIdEffEventSelector>:" << std::endl;

  tauIdDiscriminators_ = cfg.getParameter<vstring>("tauIdDiscriminators");

  std::string tauChargeMode_string = cfg.getParameter<std::string>("tauChargeMode");
  if      ( tauChargeMode_string == "tauLeadTrackCharge"        ) tauChargeMode_ = kLeadTrackCharge;
  else if ( tauChargeMode_string == "tauSignalChargedHadronSum" ) tauChargeMode_ = kSignalChargedHadronSum;
  else throw cms::Exception("TauIdEffEventSelector") 
    << "Invalid configuration parameter 'tauChargeMode' = " << tauChargeMode_string << " !!\n";

  disableTauCandPreselCuts_ = cfg.getParameter<bool>("disableTauCandPreselCuts");

//--- define default cuts for ABCD regions
  numJets_bTaggedMin_       =   0;
  numJets_bTaggedMax_       =   0;
  muonPtMin_                =  17.0; 
  muonPtMax_                =  +1.e+3; 
  muonEtaMin_               =  -2.1;
  muonEtaMax_               =  +2.1; 
  muonRelIsoMin_            =  -1.e+3;
  muonRelIsoMax_            =   0.50;
  tauPtMin_                 =  20.0; // CV: cut is actually applied on PFJet Pt, not PFTau Pt
  tauPtMax_                 =  +1.e+3; 
  tauEtaMin_                =  -2.3;
  tauEtaMax_                =  +2.3;  
  tauLeadTrackPtMin_        =   5.0; 
  tauAbsIsoMin_             =  -1.e+3;
  tauAbsIsoMax_             =   2.5;
  tauChargeMin_             =  -1.e+3;
  tauChargeMax_             =  +1.e+3;
  //muTauPairAbsDzMax_        =  +1.e+3;
  muTauPairAbsDzMax_        =   0.2; // 2mm
  muTauPairChargeProdMin_   =  -1.e+3;
  muTauPairChargeProdMax_   =  +1.e+3; 
  caloMEtPtMin_             =  20.0;
  caloMEtPtMax_             =  +1.e+2;
  pfMEtPtMin_               =  20.0;
  pfMEtPtMax_               =  +1.e+2;
  MtMin_                    =  -1.e+3;
  MtMax_                    =  40.0;
  PzetaDiffMin_             = -20.0;
  PzetaDiffMax_             =  +1.e+3;
  MtAndPzetaDiffCut_        = kNotApplied;
  tauIdDiscriminatorMin_    =   0.5;
  tauIdDiscriminatorMax_    =  +1.e+3;
  tauIdDiscriminatorCut_    = kNotApplied;

  if ( cfg.exists("muonPtMin")         ) muonPtMin_         = cfg.getParameter<double>("muonPtMin");
  if ( cfg.exists("tauLeadTrackPtMin") ) tauLeadTrackPtMin_ = cfg.getParameter<double>("tauLeadTrackPtMin");
  if ( cfg.exists("tauAbsIsoMax")      ) tauAbsIsoMax_      = cfg.getParameter<double>("tauAbsIsoMax");
  if ( cfg.exists("caloMEtPtMin")      ) caloMEtPtMin_      = cfg.getParameter<double>("caloMEtPtMin");
  if ( cfg.exists("pfMEtPtMin")        ) pfMEtPtMin_        = cfg.getParameter<double>("pfMEtPtMin");

//--- additional cuts to make sure there are no events in underflow/overflow bins
//    of visible and transverse mass distributions
//   (difference in event yields may cause problem with simultaneous fit of visMass and Mt distributions)
  visMassCutoffMin_         =  20.0;
  visMassCutoffMax_         = 200.0;
  MtCutoffMin_              =  -1.e+3;
  MtCutoffMax_              = 120.0;

  region_ = cfg.getParameter<std::string>("region");
  //std::cout << " region = " << region_ << std::endl;

  if        ( region_           == "ABCD"            ) {
    // nothing to be done yet...
    // (simply use default cuts)
  } else if ( region_.find("A") != std::string::npos ) {
    muonRelIsoMin_          =   0.20;
    muTauPairChargeProdMax_ =  -0.5;
  } else if ( region_.find("B") != std::string::npos ) {
    muonRelIsoMin_          =   0.20;
    muTauPairChargeProdMin_ =  +0.5;
  } else if ( region_.find("C") != std::string::npos ) {
    muonRelIsoMax_          =   0.10;
    muTauPairChargeProdMax_ =  -0.5;
  } else if ( region_.find("D") != std::string::npos ) {
    muonRelIsoMax_          =   0.10;
    muTauPairChargeProdMin_ =  +0.5;
  } else {
    throw cms::Exception("TauIdEffEventSelector") 
      << "Invalid region = " << region_ << " !!\n";
  }

  if      ( region_.find("1")  != std::string::npos ) MtAndPzetaDiffCut_ = kSignalLike;
  else if ( region_.find("2")  != std::string::npos ) MtAndPzetaDiffCut_ = kBackgroundLike;
  else if ( region_.find("Wj") != std::string::npos ) MtAndPzetaDiffCut_ = kWplusJetBackgroundLike;

  if      ( region_.find("p") != std::string::npos ) tauIdDiscriminatorCut_ = kSignalLike;
  else if ( region_.find("f") != std::string::npos ) tauIdDiscriminatorCut_ = kBackgroundLike;

  if      ( region_.find("+") != std::string::npos ) tauChargeMin_ = +0.5;
  else if ( region_.find("-") != std::string::npos ) tauChargeMax_ = -0.5;
}

TauIdEffEventSelector::~TauIdEffEventSelector()
{
// nothing to be done yet...
}

std::string getCutStatus_string(int cut)
{
  if      ( cut == kNotApplied             ) return "not applied";
  else if ( cut == kSignalLike             ) return "signal-like";
  else if ( cut == kBackgroundLike         ) return "background-like";
  else if ( cut == kWplusJetBackgroundLike ) return "W+jet background-like";
  else assert(0);
}

void printCutValue(const std::string& variable, double value, double min, double max)
{
  std::cout << variable << " = " << value;
  std::cout << " (";
  if   ( min <= -1.e+3 ) std::cout << "-infinity";
  else                   std::cout << min;
  std::cout << "..";
  if   ( max >= +1.e+3 ) std::cout << "+infinity";
  else                   std::cout << max;
  std::cout << ")";
  std::cout << std::endl;
}

bool TauIdEffEventSelector::operator()(const PATMuTauPair& muTauPair, const pat::MET& caloMEt, 
				       size_t numJets_bTagged, pat::strbitset& result)
{
  //std::cout << "<TauIdEffEventSelector::operator()>:" << std::endl;

  double muonPt              = muTauPair.leg1()->pt();
  double muonEta             = muTauPair.leg1()->eta();
  // compute deltaBeta corrected isolation Pt sum "by hand":
  //   muonIsoPtSum = pfChargedParticles(noPileUp) + pfNeutralHadrons + pfGammas - deltaBetaCorr, deltaBetaCorr = 0.5*pfChargedParticlesPileUp
  // ( User1Iso = pfAllChargedHadrons(noPileUp), User2Iso = pfAllChargedHadronsPileUp
  //   as defined in TauAnalysis/TauIdEfficiency/test/commissioning/produceMuonIsolationPATtuple_cfg.py )
  double muonIso             = muTauPair.leg1()->userIsolation(pat::User1Iso) 
  	                      + TMath::Max(0., muTauPair.leg1()->userIsolation(pat::PfNeutralHadronIso) 
                                              + muTauPair.leg1()->userIsolation(pat::PfGammaIso) 
                                              - 0.5*muTauPair.leg1()->userIsolation(pat::User2Iso));
  double tauPt               = muTauPair.leg2()->pt();
  double tauEta              = muTauPair.leg2()->eta();
  double tauLeadTrackPt      = muTauPair.leg2()->userFloat("leadTrackPt");
  double tauIso              = muTauPair.leg2()->userFloat("preselLoosePFIsoPt");
  double tauCharge           = 0.;
  if      ( tauChargeMode_ == kLeadTrackCharge        ) tauCharge = muTauPair.leg2()->userFloat("leadTrackCharge");
  else if ( tauChargeMode_ == kSignalChargedHadronSum ) tauCharge = muTauPair.leg2()->charge();
  else assert(0);
  double muTauPairAbsDz      = TMath::Abs(muTauPair.leg1()->vertex().z() - muTauPair.leg2()->vertex().z());
  double muTauPairChargeProd = muTauPair.leg1()->charge()*tauCharge;
  double visMass             = (muTauPair.leg1()->p4() + muTauPair.leg2()->p4()).mass();
  double caloMEtPt           = caloMEt.pt();
  double pfMEtPt             = muTauPair.met()->pt();
  double Mt                  = muTauPair.mt1MET();
  double PzetaDiff           = muTauPair.pZeta() - 1.5*muTauPair.pZetaVis();

  //printCutValue("tauLeadTrackPt", tauLeadTrackPt, tauLeadTrackPtMin_, +1.e+3);
  //printCutValue("tauIso", tauIso, tauAbsIsoMin_, tauAbsIsoMax_); 
  //printCutValue("muTauPairAbsDz", muTauPairAbsDz, -1.e+3, muTauPairAbsDzMax_);
  //printCutValue("muTauPairChargeProd", muTauPairChargeProd, muTauPairChargeProdMin_, muTauPairChargeProdMax_);
  //printCutValue("Mt", Mt, TMath::Min(MtMin_, MtCutoffMin_), TMath::Max(MtMax_, MtCutoffMax_));
  //printCutValue("PzetaDiff", PzetaDiff, PzetaDiffMin_, PzetaDiffMax_);

  //std::cout << "MtAndPzetaDiffCut = " << getCutStatus_string(MtAndPzetaDiffCut_) << std::endl;
  //std::cout << "tauIdDiscriminatorCut = " << getCutStatus_string(tauIdDiscriminatorCut_) << std::endl;

  if ( numJets_bTagged     >= numJets_bTaggedMin_     && numJets_bTagged         <= numJets_bTaggedMax_     &&
       muonPt              >  muonPtMin_              && muonPt                  <  muonPtMax_              &&
       muonEta             >  muonEtaMin_             && muonEta                 <  muonEtaMax_             &&
       muonIso             > (muonRelIsoMin_*muonPt)  && muonIso                 < (muonRelIsoMax_*muonPt)  && 
       tauPt               >  tauPtMin_               && tauPt                   <  tauPtMax_               &&
       tauEta              >  tauEtaMin_              && tauEta                  <  tauEtaMax_              &&
       tauCharge           >  tauChargeMin_           && tauCharge               <  tauChargeMax_           &&
       ((tauLeadTrackPt    >  tauLeadTrackPtMin_      &&
         tauIso            >  tauAbsIsoMin_           && tauIso                  <  tauAbsIsoMax_           ) ||
        disableTauCandPreselCuts_                                                                           ) &&
       muTauPairAbsDz      <  muTauPairAbsDzMax_      &&
       muTauPairChargeProd >  muTauPairChargeProdMin_ && muTauPairChargeProd     <  muTauPairChargeProdMax_ && 
       visMass             >  visMassCutoffMin_       && visMass                 <  visMassCutoffMax_       &&
       caloMEtPt           >  caloMEtPtMin_           && caloMEtPt               <  caloMEtPtMax_           &&
       pfMEtPt             >  pfMEtPtMin_             && pfMEtPt                 <  pfMEtPtMax_             &&
       Mt                  >  MtCutoffMin_            && Mt                      <  MtCutoffMax_            ) {
    bool MtAndPzetaDiffCut_passed = (Mt > MtMin_ && Mt < MtMax_ && PzetaDiff > PzetaDiffMin_ && PzetaDiff < PzetaDiffMax_);
    //bool MtAndPzetaDiffCut_passed = (Mt > MtMin_ && Mt < MtMax_);

    bool tauIdDiscriminators_passed = true;
    for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	  tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
      double tauIdDiscriminator_value = muTauPair.leg2()->tauID(*tauIdDiscriminator);
      //std::cout << " " << (*tauIdDiscriminator) << ": " << tauIdDiscriminator_value << std::endl;
      if ( !(tauIdDiscriminator_value > tauIdDiscriminatorMin_  && 
	     tauIdDiscriminator_value < tauIdDiscriminatorMax_) ) tauIdDiscriminators_passed = false;
    }

    if ( ( MtAndPzetaDiffCut_ == kNotApplied                                           ||
	  (MtAndPzetaDiffCut_ == kSignalLike             &&  MtAndPzetaDiffCut_passed) ||
	  (MtAndPzetaDiffCut_ == kBackgroundLike         && !MtAndPzetaDiffCut_passed) ||
	  (MtAndPzetaDiffCut_ == kWplusJetBackgroundLike &&  Mt > 60.                ) ) &&
	 ( tauIdDiscriminatorCut_ == kNotApplied                                       ||
	  (tauIdDiscriminatorCut_ == kSignalLike     &&  tauIdDiscriminators_passed)   ||
	  (tauIdDiscriminatorCut_ == kBackgroundLike && !tauIdDiscriminators_passed)   ) ) return true;
  }

  return false;
}
