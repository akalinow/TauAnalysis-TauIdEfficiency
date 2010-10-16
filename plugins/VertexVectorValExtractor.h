#ifndef TauAnalysis_TauIdEfficiency_VertexVectorValExtractor_h  
#define TauAnalysis_TauIdEfficiency_VertexVectorValExtractor_h

/** \class VertexVectorValExtractor
 *
 * Auxiliary class for extracting x,y,z coordinates of reconstructed primary event vertices
 *
 * NOTE: the x and y coordinates are computed with respect to the (offline) beam-spot
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

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

class VertexVectorValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit VertexVectorValExtractor(const edm::ParameterSet&);
  ~VertexVectorValExtractor();
    
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag srcVertex_;
  edm::InputTag srcBeamSpot_;
  
  enum { kVertexX, kVertexY, kVertexZ };

  int value_;
};

#endif  


