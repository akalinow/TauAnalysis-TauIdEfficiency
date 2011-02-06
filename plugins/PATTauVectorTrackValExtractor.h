#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorTrackValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorTrackValExtractor_h

/** \class PATTauVectorTrackValExtractor
 *
 * Auxiliary class for extracting track quantities
 * of leading track, tracks in signal cone or tracks in isolation cone
 * of reconstructed PAT tau objects
 * (used for Ntuple filling)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATTauVectorTrackValExtractor.h,v 1.2 2010/09/28 11:23:39 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include <string>

class PATTauVectorTrackValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorTrackValExtractor(const edm::ParameterSet&);
  ~PATTauVectorTrackValExtractor(){};
    
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- auxiliary functions
  const reco::Track* getTrack_pfTau(const pat::Tau&) const;
  const reco::Track* getTrack_caloTau(const pat::Tau&) const;
  
//--- configuration parameters
  edm::InputTag src_;
  
  enum { kLeadTrack, kSignalConeTracks, kIsolationConeTracks };
  int collection_;

  unsigned index_;

  enum { kNumValidHits, kNumMissingHits, kChi2, kNumDoF, kDeltaZ, kDeltaXY, kPt, kPtErr, kQuality };
  int value_;

  edm::InputTag srcVertex_;
};

#endif  


