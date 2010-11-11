#ifndef TauAnalysis_TauIdEfficiency_GenPileUpEventInfoExtractor_h  
#define TauAnalysis_TauIdEfficiency_GenPileUpEventInfoExtractor_h

/** \class GenPileUpEventInfoExtractor
 *
 * Auxiliary class for extracting information on Monte Carlo generator level
 * about multiplicity of mixed-in pile-up interactions
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

class GenPileUpEventInfoExtractor : public ObjValExtractorBase
{
 public:
  
  explicit GenPileUpEventInfoExtractor(const edm::ParameterSet&);
  ~GenPileUpEventInfoExtractor();
 
  double operator()(const edm::Event&) const;

 private:
//--- configuration parameters
  edm::InputTag src_;

  enum { kNumPileUpInteractions };

  int value_;
};

#endif  


