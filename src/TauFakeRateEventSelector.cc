#include "TauAnalysis/TauIdEfficiency/interface/TauFakeRateEventSelector.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

// define flag for tau id. discriminators 
// (no tau id. discriminators applied, all discriminators passed, at least one discriminator failed)
enum { kNotApplied, kSignalLike, kBackgroundLike };

TauFakeRateEventSelector::TauFakeRateEventSelector(const edm::ParameterSet& cfg)
{
//--- get tau id. discriminators for which the fake-rate is to be determined
//    and whether to fill histograms for jets passing or the tau id. discriminators
  tauIdDiscriminators_ = cfg.getParameter<vstring>("tauIdDiscriminators");

  region_ = cfg.getParameter<std::string>("region");

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
    tauJetCandPreselCriteria_.push_back(new StringCutTauSelectorType(*tauJetCandPreselCriterion)); 
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
  for ( std::vector<StringCutTauSelectorType*>::iterator it = tauJetCandPreselCriteria_.begin();
	it != tauJetCandPreselCriteria_.end(); ++ it ) {
    delete (*it);
  }
}

bool TauFakeRateEventSelector::operator()(const pat::Tau& tauJetCand, pat::strbitset& result)
{
  //std::cout << "<TauFakeRateEventSelector::operator()>:" << std::endl;

  double jetPt  = tauJetCand.p4Jet().pt();
  double jetEta = tauJetCand.p4Jet().eta();

//--- check if tau-jet candidates passes Pt and eta cuts
  if ( jetPt  > jetPtMin_  && jetPt  < jetPtMax_  &&
       jetEta > jetEtaMin_ && jetEta < jetEtaMax_ ) {
    //std::cout << "passed Pt and eta cuts." << std::endl;

//--- check if tau-jet candidates passes preselection criteria
    bool preselCriteria_passed = true;
    for ( std::vector<StringCutTauSelectorType*>::iterator tauJetCandPreselCriterion = tauJetCandPreselCriteria_.begin();
	  tauJetCandPreselCriterion != tauJetCandPreselCriteria_.end(); ++tauJetCandPreselCriterion ) {
      //std::cout << "checking " << (*tauJetCandPreselCriterion)->cut_ << std::endl;
      if ( !(*(*tauJetCandPreselCriterion)->selector_)(tauJetCand) ) {
	//std::cout << " failed." << std::endl;
	preselCriteria_passed = false;
	break;
      }
    }

    if ( preselCriteria_passed ) {
      bool tauIdDiscriminators_passed = true;
      for ( vstring::const_iterator tauIdDiscriminator = tauIdDiscriminators_.begin();
	    tauIdDiscriminator != tauIdDiscriminators_.end(); ++tauIdDiscriminator ) {
	double tauIdDiscriminator_value = tauJetCand.tauID(*tauIdDiscriminator);
	//std::cout << "checking " << (*tauIdDiscriminator) << ": " << tauIdDiscriminator_value << std::endl;
	if ( !(tauIdDiscriminator_value > tauIdDiscriminatorMin_  && 
	       tauIdDiscriminator_value < tauIdDiscriminatorMax_) ) {
	  tauIdDiscriminators_passed = false;
	  break;
	}
	//std::cout << " passed." << std::endl;
      }

      tauIdDiscriminators_passed &= (tauJetCand.pt() > 15.0); // require tauPt > 15 GeV, in order to compare with "old" HPS/TaNC results

      //std::cout << "tauIdDiscriminatorCut = ";
      //if      ( tauIdDiscriminatorCut_ == kNotApplied     ) std::cout << "not applied.";
      //else if ( tauIdDiscriminatorCut_ == kSignalLike     ) std::cout << "signal-like.";
      //else if ( tauIdDiscriminatorCut_ == kBackgroundLike ) std::cout << "background-like.";
      //else std::cout << "undefined.";
      //std::cout << std::endl;

      if ( tauIdDiscriminatorCut_ == kNotApplied                                     ||
	  (tauIdDiscriminatorCut_ == kSignalLike     &&  tauIdDiscriminators_passed) ||
	  (tauIdDiscriminatorCut_ == kBackgroundLike && !tauIdDiscriminators_passed) ) return true;
    }
  }

  return false;
}
