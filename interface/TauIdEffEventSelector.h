#ifndef TauAnalysis_TauIdEfficiency_TauIdEffEventSelector_h
#define TauAnalysis_TauIdEfficiency_TauIdEffEventSelector_h

/** \class TauIdEffEventSelector
 *
 * Select muon + tau-jet pairs entering signal/control regions
 * of tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: TauIdEffEventSelector.h,v 1.2 2011/07/01 18:30:16 veelken Exp $
 *
 */

#include "PhysicsTools/SelectorUtils/interface/EventSelector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

class TauIdEffEventSelector : public EventSelector 
{

 public:
  /// constructor
  TauIdEffEventSelector(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffEventSelector();

  /// here is where the selection occurs
  bool operator()(const edm::EventBase& event, pat::strbitset& result) { return true; }
  bool operator()(const PATMuTauPair& muTauPair, pat::strbitset& result);

 private:

  /// specify region in which to select events
  std::string region_;
 
  /// list of tau id. discrimators
  /// (e.g. 'decayModeFinding' && 'byLooseCombinedIsolationDeltaBetaCorr')
  typedef std::vector<std::string> vstring;
  vstring tauIdDiscriminators_;

  /// cuts applied in specified region
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
  double muTauPairAbsDzMax_;
  double muTauPairChargeProdMin_;
  double muTauPairChargeProdMax_;
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
