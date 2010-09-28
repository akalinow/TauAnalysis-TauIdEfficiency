#ifndef TauAnalysis_TauIdEfficiency_HLTInfoExtractor_h  
#define TauAnalysis_TauIdEfficiency_HLTInfoExtractor_h

/** \class HLTInfoExtractor
 *
 * Auxiliary class for extracting HLT Info for the event
 * (used for Ntuple filling)
 *
 * \author Michail Bachtis, U.Wisconsin
 *
 * \version $Revision: 1.3 $
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

class HLTInfoExtractor : public ObjValExtractorBase
{
 public:
  
  explicit HLTInfoExtractor(const edm::ParameterSet&);
  ~HLTInfoExtractor();
 
  double operator()(const edm::Event&) const;

 private:
  //--- configuration parameters
  edm::InputTag src_;
  std::string value_;

  int maxWarnings_;
  mutable int numWarnings_;
};

#endif  


