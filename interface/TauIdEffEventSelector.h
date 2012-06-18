#ifndef TauAnalysis_TauIdEfficiency_TauIdEffEventSelector_h
#define TauAnalysis_TauIdEfficiency_TauIdEffEventSelector_h

/** \class TauIdEffEventSelector
 *
 * Select muon + tau-jet pairs entering signal/control regions
 * of tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.8 $
 *
 * $Id: TauIdEffEventSelector.h,v 1.8 2011/11/06 13:25:26 veelken Exp $
 *
 */

#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class TauIdEffEventSelector : public EventSelector 
{

 public:
  /// constructor
  TauIdEffEventSelector(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffEventSelector();

  /// here is where the selection occurs
  bool operator()(const edm::EventBase&, pat::strbitset&) { return true; }
  bool operator()(const PATMuTauPair&, const pat::MET&, size_t, pat::strbitset&);

  friend class regionEntryType; // allow regionEntryType to overwrite cut values

 private:

  /// specify region in which to select events
  std::string region_;
 
  /// list of tau id. discrimators
  /// (e.g. 'decayModeFinding' && 'byLooseCombinedIsolationDeltaBetaCorr')
  typedef std::vector<std::string> vstring;
  vstring tauIdDiscriminators_;

  /// flag indicating whether to take charge of tau-jet candidate 
  /// from "leading track" or from all "signal" charged hadrons
  int tauChargeMode_;

  /// flag to disable "leading" track Pt > 5 GeV && pfIso < 2.5 GeV cuts applied in preselection of tau-jet candidates
  /// CV: do not apply preselection cuts when measuring tau charge misidentification rate
  bool disableTauCandPreselCuts_;

  /// cuts applied in specified region
  size_t numJets_bTaggedMin_;
  size_t numJets_bTaggedMax_;
  double muonPtMin_;
  double muonPtMax_;  
  double muonEtaMin_;
  double muonEtaMax_;  
  double muonRelIsoMin_;
  double muonRelIsoMax_;
  double tauPtMin_;
  double tauPtMax_;  
  double tauEtaMin_;
  double tauEtaMax_;
  double tauLeadTrackPtMin_;
  double tauAbsIsoMin_;
  double tauAbsIsoMax_;
  double tauChargeMin_;
  double tauChargeMax_;
  double muTauPairAbsDzMax_;
  double muTauPairChargeProdMin_;
  double muTauPairChargeProdMax_;
  double caloMEtPtMin_;
  double caloMEtPtMax_;
  double pfMEtPtMin_;
  double pfMEtPtMax_;
  double MtMin_;
  double MtMax_;
  double PzetaDiffMin_;
  double PzetaDiffMax_;
  int    MtAndPzetaDiffCut_;

  double tauIdDiscriminatorMin_;
  double tauIdDiscriminatorMax_;
  int    tauIdDiscriminatorCut_;

  double visMassCutoffMin_;
  double visMassCutoffMax_;
  double MtCutoffMin_;
  double MtCutoffMax_;
};

#endif
