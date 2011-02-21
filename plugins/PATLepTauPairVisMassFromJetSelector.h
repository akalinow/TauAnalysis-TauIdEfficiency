#ifndef TauAnalysis_TauIdEfficiency_PATLepTauPairVisMassFromJetSelector_h
#define TauAnalysis_TauIdEfficiency_PATLepTauPairVisMassFromJetSelector_h

/** \class PATLepTauPairVisMassFromJetSelector
 *
 * Select lepton + tau-jet pairs based on visible invariant mass of lepton + tau-jet candidate,
 * with momentum of the tau-jet candidate taken from reco::PFJet instead of from reco::PFTau object
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATLepTauPairVisMassFromJetSelector.h,v 1.1 2011/01/21 10:08:06 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

template<typename T>
class PATLepTauPairVisMassFromJetSelector 
{
 public:
  
  typedef std::vector<CompositePtrCandidateT1T2MEt<T,pat::Tau> > collection;

  explicit PATLepTauPairVisMassFromJetSelector(const edm::ParameterSet&);
  ~PATLepTauPairVisMassFromJetSelector();

  typename std::vector<const CompositePtrCandidateT1T2MEt<T,pat::Tau> *>::const_iterator begin() const { return selected_.begin(); }
  typename std::vector<const CompositePtrCandidateT1T2MEt<T,pat::Tau> *>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<collection>&, edm::Event&, const edm::EventSetup&);
    
  size_t size() const { return selected_.size(); }

 private:

  std::vector<const CompositePtrCandidateT1T2MEt<T,pat::Tau> *> selected_;

//--- configuration parameters
  edm::InputTag src_;

  double minVisMass_;
  double maxVisMass_;
};

#endif  


