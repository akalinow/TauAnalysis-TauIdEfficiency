#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorDecayModeValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorDecayModeValExtractor_h

/** \class PATTauVectorDecayModeValExtractor 
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
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include <string>

class PATTauVectorDecayModeValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorDecayModeValExtractor(const edm::ParameterSet&);
  ~PATTauVectorDecayModeValExtractor(){};
    
  std::vector<double> operator()(const edm::Event&) const;

 private:
  
//--- configuration parameters
  edm::InputTag src_;
  
  int collection_;

  unsigned index_;

  enum { kEnergy, kPt, kEta, kPhi, kMass };
  int value_;

  edm::InputTag srcPFTau_;
  edm::InputTag srcPFTauDecayMode_;
};

#endif  


