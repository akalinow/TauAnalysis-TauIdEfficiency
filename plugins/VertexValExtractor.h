#ifndef TauAnalysis_TauIdEfficiency_VertexValExtractor_h  
#define TauAnalysis_TauIdEfficiency_VertexValExtractor_h

/** \class VertexValExtractor
 *
 * Auxiliary class for extracting the multiplicity of reconstructed 
 * primary event vertices with sum(trackPt) exceeding a configurable threshold
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: VertexValExtractor.h,v 1.2 2010/09/28 11:23:39 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

class VertexValExtractor : public ObjValExtractorBase
{
 public:
  
  explicit VertexValExtractor(const edm::ParameterSet&);
  ~VertexValExtractor();
    
  double operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag srcVertex_;

  double trackSumPtThreshold_;

  int error_;
};

#endif  


