#ifndef TauAnalysis_TauIdEfficiency_PATPFTauSelectorForTauIdEff_h
#define TauAnalysis_TauIdEfficiency_PATPFTauSelectorForTauIdEff_h

/** \class PATPFTauSelectorForTauIdEff
 *
 * Preselect (PF)Tau-jet candidates for tau id. efficiency measurement
 *
 * NOTE: Class produces collection of pat::Tau objects the four-vectors of which 
 *       are set to (jet-energy corrected) PFJet four-vectors.
 *       The "original" components of the PFTau four-vector are stored in userFloat variables
 *       'tauPt', 'tauEta', 'tauPhi', 'tauMass'
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATPFTauSelectorForTauIdEff.cc,v 1.1 2010/11/06 14:43:54 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/PtComparator.h"

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"

#include "TauAnalysis/RecoTools/interface/ParticlePFIsolationExtractor.h"

#include <string>

class PATPFTauSelectorForTauIdEff : public edm::EDFilter
{  
 public:

  explicit PATPFTauSelectorForTauIdEff(const edm::ParameterSet&);
  virtual ~PATPFTauSelectorForTauIdEff();

  bool filter(edm::Event&, const edm::EventSetup&);
  
private:

  bool filter_;

  edm::InputTag src_;

  std::string jetEnergyCorrection_;

  double minJetPt_;
  double maxJetEta_;

  reco::tau::RecoTauQualityCuts* trackQualityCuts_;
  double minLeadTrackPt_;
  double maxDzLeadTrack_;
  double maxLeadTrackPFElectronMVA_;
  bool applyECALcrackVeto_;

  double minDeltaRtoNearestMuon_;
  StringCutObjectSelector<pat::Muon>* muonSelection_;
  edm::InputTag srcMuon_;

  ParticlePFIsolationExtractor<pat::Tau>* pfIsolationExtractor_;
  double maxPFIsoPt_;
  edm::InputTag srcPFIsoCandidates_;
  edm::InputTag srcBeamSpot_;
  edm::InputTag srcVertex_;
  edm::InputTag srcRhoFastJet_;

  // utility to sort (PF)tau-jet candidates in order of decreasing Pt
  // (taken from PhysicsTools/PatAlgos/plugins/PATTauProducer.h)
  GreaterByPt<pat::Tau> pfTauPtComparator_;
};

#endif
