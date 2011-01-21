#ifndef TauAnalysis_TauIdEfficiency_PATLepTauPairVisMassFromJetExtractor_h
#define TauAnalysis_TauIdEfficiency_PATLepTauPairVisMassFromJetExtractor_h

/** \class PATLepTauPairVisMassFromJetExtractor
 *
 * Auxiliary class for extracting visible invariant mass of lepton + tau-jet candidate,
 * with momentum of tau-jet candidate taken from reco::PFJet instead of from reco::PFTau object
 * (used for Ntuple filling)
 *
 * NOTE: the values are extracted from the second leg of the PATElecTauPair/PATMuTauPair object 
 *       specified by the "index" configuration parameter (**first** PATElecTauPair/PATMuTauPair object in case "index" is not specified)
 *       contained in the collection specified by the "src" configuration parameter;
 *       in case the collection of PAT objects is empty, 
 *       a substitute value of -1. is returned by operator()
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATLepTauPairVisMassFromJetExtractor.h,v 1.2 2010/09/28 11:23:28 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

template<typename T>
class PATLepTauPairVisMassFromJetExtractor : public ObjValExtractorBase
{
 public:
  
  explicit PATLepTauPairVisMassFromJetExtractor(const edm::ParameterSet&);
  ~PATLepTauPairVisMassFromJetExtractor();
  
  double operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;

  unsigned index_;
};

#endif  


