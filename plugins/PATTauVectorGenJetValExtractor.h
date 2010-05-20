#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorGenJetValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorGenJetValExtractor_h

/** \class PATTauVectorGenJetValExtractor
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
 * $Id: PATTauVectorGenJetValExtractor.h,v 1.2 2009/10/25 12:38:14 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

class PATTauVectorGenJetValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorGenJetValExtractor(const edm::ParameterSet&);
  ~PATTauVectorGenJetValExtractor();
  
  unsigned int size() const;
  
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;
  
  enum { kGenMatch, kGenPt, kGenEta, kGenPhi, kGenDecayMode };

  int value_;
};

#endif  


