#ifndef TauAnalysis_TauIdEfficiency_GenPhaseSpaceEventInfoExtractor_h  
#define TauAnalysis_TauIdEfficiency_GenPhaseSpaceEventInfoExtractor_h

/** \class GenPhaseSpaceEventInfoExtractor
 *
 * Auxiliary class for extracting information from Monte Carlo generator 
 * (used for Ntuple filling; e.g. PtHat value in events produced by PYTHIA)
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

class GenPhaseSpaceEventInfoExtractor : public ObjValExtractorBase
{
 public:
  
  explicit GenPhaseSpaceEventInfoExtractor(const edm::ParameterSet&);
  ~GenPhaseSpaceEventInfoExtractor();
 
  double operator()(const edm::Event&) const;

 private:
//--- configuration parameters
  edm::InputTag src_;

  enum { kPtHat };

  int value_;
};

#endif  


