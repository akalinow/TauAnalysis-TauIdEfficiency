#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorJetIdValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorJetIdValExtractor_h

/** \class PATTauVectorJetIdValExtractor 
 *
 * Auxiliary class for extracting jetId bits
 * for reco::CaloJet/reco::PFJet object associated 
 * to reconstructed PAT tau objects
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
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include <string>

template <typename T>
class PATTauVectorJetIdValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorJetIdValExtractor(const edm::ParameterSet&);
  ~PATTauVectorJetIdValExtractor();
    
  std::vector<double> operator()(const edm::Event&) const;

 private:
  
//--- configuration parameters
  edm::InputTag src_;

  T* jetId_;

  edm::InputTag srcJet_;
};

#endif  

