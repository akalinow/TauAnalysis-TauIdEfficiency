#ifndef TauAnalysis_TauIdEfficiency_TauFakeRateEventSelector_h
#define TauAnalysis_TauIdEfficiency_TauFakeRateEventSelector_h

/** \class TauFakeRateEventSelector
 *
 * Select tau-jet candidates entering passed/failed samples
 * of jet --> tau fake-rate measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: TauFakeRateEventSelector.h,v 1.4 2011/07/05 17:07:22 veelken Exp $
 *
 */

#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

class TauFakeRateEventSelector : public EventSelector 
{

 public:
  /// constructor
  TauFakeRateEventSelector(edm::ParameterSet const&);

  /// destructor
  virtual ~TauFakeRateEventSelector();

  /// here is where the selection occurs
  bool operator()(const edm::EventBase& event, pat::strbitset& result) { return true; }
  bool operator()(const pat::Tau& tauJetCand, pat::strbitset& result);

  friend class regionEntryType; // allow regionEntryType to overwrite cut values

 private:

  /// specify region in which to select events
  std::string region_;
 
  /// list of tau id. discrimators
  /// (e.g. 'decayModeFinding' && 'byLooseCombinedIsolationDeltaBetaCorr')
  typedef std::vector<std::string> vstring;
  vstring tauIdDiscriminators_;

  double tauIdDiscriminatorMin_;
  double tauIdDiscriminatorMax_;
  int    tauIdDiscriminatorCut_;

  /// list of preselection criteria applied on tau-jet candidates
  /// to enter fake-rate measurement
  /// (e.g. jetId, discriminators against electrons/muons)
  typedef StringCutObjectSelector<pat::Tau> StringCutTauSelector;
  std::vector<StringCutTauSelector*> tauJetCandPreselCriteria_;

  /// Pt and eta cuts
  double jetPtMin_;
  double jetPtMax_;  
  double jetEtaMin_;
  double jetEtaMax_;  
};

#endif
