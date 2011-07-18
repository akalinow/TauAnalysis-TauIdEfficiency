#include "TauAnalysis/TauIdEfficiency/interface/TauFakeRateEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

// define flag for tau id. discriminators 
// (no tau id. discriminators applied, all discriminators passed, at least one discriminator failed)
enum { kNotApplied, kSignalLike, kBackgroundLike };

TauFakeRateEventSelector::TauFakeRateEventSelector(const edm::ParameterSet& cfg)
{
//--- get tau id. discriminators for which the fake-rate is to be determined
  tauIdDiscriminators_ = cfg.getParameter<vstring>("tauIdDiscriminators");

//--- define default preselection criteria for tau-jet candidates
  vstring tauJetCandPreselCriteria_string;
  tauJetCandPreselCriteria_string.push_back("userFloat('jetIdLoose') > 0.5");
  tauJetCandPreselCriteria_string.push_back("tauID('againstElectronLoose') > 0.5");
  tauJetCandPreselCriteria_string.push_back("tauID('againstMuonTight') > 0.5");

//--- check if additional preselection criteria are to be applied
//   (e.g. 'tag'/'probe' flags for HLT single jet trigger matching)
  if ( cfg.exists("tauJetCandPreselCriteria") ) {
    vstring tauJetCandPreselCriteria_custom = cfg.getParameter<vstring>("tauJetCandPreselCriteria");
    tauJetCandPreselCriteria_string.insert(
      tauJetCandPreselCriteria_string.end(), tauJetCandPreselCriteria_custom.begin(), tauJetCandPreselCriteria_custom.end());
  }

//--- create CutString objects for all preselection criteria
  for ( vstring::const_iterator tauJetCandPreselCriterion = tauJetCandPreselCriteria_string.begin();
	tauJetCandPreselCriterion != tauJetCandPreselCriteria_string.end(); ++tauJetCandPreselCriterion ) {
    tauJetCandPreselCriteria_.push_back(new StringCutTauSelector(*tauJetCandPreselCriterion)); 
  }
  
//--- define default Pt and eta cuts for tau-jet candidates
  jetPtMin_  = 20.0; 
  jetPtMax_  = +1.e+3; 
  jetEtaMin_ = -2.3;
  jetEtaMax_ = +2.3; 

  tauIdDiscriminatorMin_ =   0.5;
  tauIdDiscriminatorMax_ =  +1.e+3;
  tauIdDiscriminatorCut_ = kNotApplied;

  if        ( region_.find("A") != std::string::npos ) { 
    // nothing to be done yet...
    // (simply use default cuts)
  } else if ( region_.find("P") != std::string::npos ) { // require that tau-jet candidates passed tau id. criteria
    tauIdDiscriminatorCut_ = kSignalLike;
  } else if ( region_.find("F") != std::string::npos ) { // require that tau-jet candidates fails tau id. criteria    
    tauIdDiscriminatorCut_ = kBackgroundLike;
  } else {
    throw cms::Exception("TauFakeRateEventSelector") 
      << "Invalid region = " << region_ << " !!\n";
  }
}

TauFakeRateEventSelector::~TauFakeRateEventSelector()
{ 
  for ( std::vector<StringCutTauSelector*>::iterator it = tauJetCandPreselCriteria_.begin();
	it != tauJetCandPreselCriteria_.end(); ++it ) {
    delete (*it);
  }
}

bool TauFakeRateEventSelector::operator()(const pat::Tau& tauJetCand, pat::strbitset& result)
{
  //std::cout << "<TauFakeRateEventSelector::operator()>:" << std::endl;

  reco::Candidate::LorentzVector p4Jet;
  if      ( tauJetCand.isCaloTau() ) p4Jet = tauJetCand.caloTauTagInfoRef()->jetRef()->p4();
  else if ( tauJetCand.isPFTau()   ) p4Jet = tauJetCand.pfJetRef()->p4();
  else throw cms::Exception("TauFakeRateHistManager::fillHistograms") 
    << "Tau-jet candidate passed as function argument is neither PFTau nor CaloTau !!";

  double jetPt  = p4Jet.pt();
  double jetEta = p4Jet.eta();

//--- check if tau-jet candidates passes Pt and eta cuts
  if ( jetPt  > jetPtMin_  && jetPt  < jetPtMax_  &&
       jetEta > jetEtaMin_ && jetEta < jetEtaMax_ ) {

//--- check if tau-jet candidates passes preselection criteria
    bool preselCriteria_passed = true;
    for ( std::vector<StringCutTauSelector*>::iterator tauJetCandPreselCriterion = tauJetCandPreselCriteria_.begin();
	  tauJetCandPreselCriterion != tauJetCandPreselCriteria_.end(); ++tauJetCandPreselCriterion ) {
      if ( !(**tauJetCandPreselCriterion)(tauJetCand) ) {
	preselCriteria_passed = false;
	break;
      }
    }

    if ( preselCriteria_passed ) {
      bool tauIdDiscriminators_passed = true;
      for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	    tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
	double tauIdDiscriminator_value = tauJetCand.tauID(*tauIdDiscriminator);
	//std::cout << " " << (*tauIdDiscriminator) << ": " << tauIdDiscriminator_value << std::endl;
	if ( !(tauIdDiscriminator_value > tauIdDiscriminatorMin_  && 
	       tauIdDiscriminator_value < tauIdDiscriminatorMax_) ) tauIdDiscriminators_passed = false;
      }
      
      if ( tauIdDiscriminatorCut_ == kNotApplied                                     ||
	  (tauIdDiscriminatorCut_ == kSignalLike     &&  tauIdDiscriminators_passed) ||
	  (tauIdDiscriminatorCut_ == kBackgroundLike && !tauIdDiscriminators_passed) ) return true;
    }
  }

  return false;
}
