#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorGenJetValExtractor_h  
#define TauAnalysis_TauIdEfficiency_VertexValExtractor_h

/** \class VertexValExtractor
 *
 * Auxiliary class for extracting x,y,z coordinates of reconstructed primary event vertex
 *
 * NOTE: the x and y coordinates are computed with respect to the (offline) beam-spot
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: VertexValExtractor.h,v 1.1 2010/05/20 17:34:34 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
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
  edm::InputTag srcBeamSpot_;
  
  enum { kVertexX, kVertexY, kVertexZ };

  int value_;
};

#endif  


