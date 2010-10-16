#ifndef TauAnalysis_BgEstimationTools_NumSelObjExtractor_h  
#define TauAnalysis_BgEstimationTools_NumSelObjExtractor_h

/** \class NumSelObjExtractor
 *
 * Auxiliary class for extracting number of PAT objects
 * not overlapping with other objects and passing selection
 * (used for Ntuple filling)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: NumSelObjExtractor.h,v 1.2 2010/09/28 11:23:25 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <vector>

template<typename T>
class NumSelObjExtractor : public ObjValExtractorBase
{
 public:
  
  explicit NumSelObjExtractor(const edm::ParameterSet&);
  ~NumSelObjExtractor();
  
  double operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcNotToBeFiltered_;
  double dRmin_;

  StringCutObjectSelector<T>* cut_;
};

#endif  


