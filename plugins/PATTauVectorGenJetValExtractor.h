#ifndef TauAnalysis_TauIdEfficiency_VectorGenJetValExtractor_h  
#define TauAnalysis_TauIdEfficiency_VectorGenJetValExtractor_h

/** \class VectorGenJetValExtractor
 *
 * Auxiliary class for extracting generator level information
 * matching reconstructed PAT tau objects
 * (used for Ntuple filling)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PATTauVectorGenJetValExtractor.h,v 1.3 2010/06/05 00:50:30 friis Exp $
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


