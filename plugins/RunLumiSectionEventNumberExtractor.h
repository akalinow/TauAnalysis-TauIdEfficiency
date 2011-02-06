#ifndef TauAnalysis_TauIdEfficiency_PATTriggerInfoExtractor_h  
#define TauAnalysis_TauIdEfficiency_RunLumiSectionEventNumberExtractor_h

/** \class RunLumiSectionEventNumberExtractor
 *
 * Auxiliary class for extracting run number, luminosity section and event number information for the event
 * (used for Ntuple filling)
 *
 * \author Christian, UC Davis
 *
 * \version $Revision: 1.3 $
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <string>

class RunLumiSectionEventNumberExtractor : public ObjValExtractorBase
{
 public:
  
  explicit RunLumiSectionEventNumberExtractor(const edm::ParameterSet&);
  ~RunLumiSectionEventNumberExtractor();
 
  double operator()(const edm::Event&) const;

 private:
//--- configuration parameters

  enum { kRunNumber, kLumiSection, kEventNumber };
  int value_;

  int maxWarnings_;
  mutable int numWarnings_;

  int cfgError_;
};

#endif  


