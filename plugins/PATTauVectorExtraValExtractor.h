#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorExtraValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorExtraValExtractor_h

/** \class PATTauVectorExtraValExtractor
 *
 * Auxiliary class for extracting generator level information
 * matching reconstructed PAT tau objects
 * (used for Ntuple filling)
 *
 * NOTE: the values are extracted from the PAT object
 *       specified by the "index" configuration parameter (**first** PAT object in case "index" is not specified)
 *       contained in the collection specified by the "src" configuration parameter;
 *       in case the collection of PAT objects is empty, 
 *       a substitute value of -1. is returned by operator()
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: PATTauVectorExtraValExtractor.h,v 1.2 2010/06/10 08:14:07 veelken Exp $
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


