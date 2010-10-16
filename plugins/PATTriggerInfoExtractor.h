#ifndef TauAnalysis_TauIdEfficiency_PATTriggerInfoExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTriggerInfoExtractor_h

/** \class HLTInfoExtractor
 *
 * Auxiliary class for extracting trigger bit and prescale information for the event
 * (used for Ntuple filling)
 *
 * \author Christian, UC Davis
 *
 * \version $Revision: 1.4 $
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <string>

class PATTriggerInfoExtractor : public ObjValExtractorBase
{
 public:
  
  explicit PATTriggerInfoExtractor(const edm::ParameterSet&);
  ~PATTriggerInfoExtractor();
 
  double operator()(const edm::Event&) const;

 private:
  //--- configuration parameters
  edm::InputTag src_;

  std::string hltPath_;

  enum { kBit, kPrescale };
  int value_;

  int maxWarnings_;
  mutable int numWarnings_;

  int cfgError_;
};

#endif  


