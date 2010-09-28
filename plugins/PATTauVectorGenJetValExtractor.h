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
 * \version $Revision: 1.6 $
 *
 * $Id: PATTauVectorGenJetValExtractor.h,v 1.6 2010/09/27 10:34:40 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include <vector>

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
  
  edm::InputTag srcGenParticles_;

  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch_;

  enum { kGenMatch, kGenPt, kGenEta, kGenPhi, kGenMass, kGenDecayMode, kGenPdgId };

  int value_;
};

#endif  


