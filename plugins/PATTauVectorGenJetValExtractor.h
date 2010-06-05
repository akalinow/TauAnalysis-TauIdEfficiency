#ifndef TauAnalysis_TauIdEfficiency_VectorGenJetValExtractor_h  
#define TauAnalysis_TauIdEfficiency_VectorGenJetValExtractor_h

/** \class VectorGenJetValExtractor
 *
 * Auxiliary class for extracting generator level information
 * matching reconstructed PAT tau objects
 * (used for Ntuple filling)
 *
 * NOTE: the values are extracted from the PAT object
 *       specified by the "index" configuration parameter (**first** PAT object in case "index" is not specified)
 *       contained in the collection specified by the "src" configuration parameter;
 *       in case the collection of PAT objects is empty, 
 *       a substitute value of -1. is returned by operator()
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: VectorGenJetValExtractor.h,v 1.2 2010/05/22 16:47:50 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

template<class T>
class VectorGenJetValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit VectorGenJetValExtractor(const edm::ParameterSet&);
  ~VectorGenJetValExtractor(){};
    
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;
  
  enum { kGenMatch, kGenPt, kGenEta, kGenPhi, kGenDecayMode };

  int value_;
};

#endif  


