#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorExtraValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorExtraValExtractor_h

/** \class PATTauVectorExtraValExtractor
 *
 * Auxiliary class for extracting generator level information
 * matching reconstructed PAT tau objects
 * (used for Ntuple filling)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PATTauVectorExtraValExtractor.h,v 1.3 2010/06/14 08:34:03 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

class PATTauVectorExtraValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorExtraValExtractor(const edm::ParameterSet&);
  ~PATTauVectorExtraValExtractor();
    
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;

  edm::InputTag pfCandSrc_;
  edm::InputTag jetSrc_;

  double jetMinPt_;
  double jetMaxAbsEta_;

  enum { kNumTracksOut, kNumChargedHadrOut, kNumPhotonsOut, 
	 kNearestJetDeltaR, kNearestJetPt, kNearestJetEta, kNearestJetPhi, kNearestJetWidth };

  int value_;
};

#endif  


