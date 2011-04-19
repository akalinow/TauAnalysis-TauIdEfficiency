#ifndef TauAnalysis_TauIdEfficiency_PATTauVectorJetCorrMomValExtractor_h  
#define TauAnalysis_TauIdEfficiency_PATTauVectorJetCorrMomValExtractor_h

/** \class PATTauVectorJetCorrMomValExtractor 
 *
 * Auxiliary class for extracting jet enery corrected
 * momentum of reco::CaloJet/reco::PFJet object associated 
 * to reconstructed PAT tau objects
 * (used for Ntuple filling)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATTauVectorJetCorrMomValExtractor.h,v 1.1 2011/02/06 10:41:00 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include <string>

class PATTauVectorJetCorrMomValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit PATTauVectorJetCorrMomValExtractor(const edm::ParameterSet&);
  ~PATTauVectorJetCorrMomValExtractor();
    
  std::vector<double> operator()(const edm::Event&) const;

 private:
  
//--- configuration parameters
  edm::InputTag src_;

  std::string value_;

  StringObjectFunction<pat::Jet> stringObjFunction_;  

  edm::InputTag srcJet_;
};

#endif  

