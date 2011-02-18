#ifndef TauAnalysis_TauIdEfficiency_ObjValEDNtupleProducer_h  
#define TauAnalysis_TauIdEfficiency_ObjValEDNtupleProducer_h

/** \class ObjValEDNtupleProducer
 *
 * Produce an Ntuple of various quantities extracted via 
 * the ObjValExtractor and store it in the edm:Event
 *
 * \author Christian Veelken, Evan Friis, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: ObjValEDNtupleProducer.h,v 1.3 2010/09/28 11:23:39 jkolb Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include <string>
#include <vector>
#include <memory>

class ObjValEDNtupleProducer : public edm::EDProducer
{
  struct ntupleEntryType
  {
    ntupleEntryType(const std::string& ntupleName, ObjValExtractorBase* objValExtractor)
      : ntupleName_(ntupleName), 
	objValExtractor_(objValExtractor) 
    {}
    ~ntupleEntryType() { delete objValExtractor_; }
    std::string ntupleName_;
    ObjValExtractorBase* objValExtractor_;
  };
  
  struct ntupleVectorEntryType
  {
    ntupleVectorEntryType(const std::string& ntupleName, ObjValVectorExtractorBase* objValExtractor, 
			  const std::vector<unsigned>& indices)
      : ntupleName_(ntupleName), 
	objValExtractor_(objValExtractor),
	indices_(indices)
    {
      //std::cout << "<ntupleVectorEntryType>:" << std::endl;
      //std::cout << " indices = " << format_vunsigned(indices_) << std::endl;
    }
    ~ntupleVectorEntryType() { delete objValExtractor_; }
    std::string ntupleName_;
    ObjValVectorExtractorBase* objValExtractor_;
    std::vector<unsigned> indices_;
  };

 public:
  
  explicit ObjValEDNtupleProducer(const edm::ParameterSet&);
  ~ObjValEDNtupleProducer();
  
 private:

  void beginJob();
  void produce(edm::Event&, const edm::EventSetup&);
  void endJob() {}

//--- configuration parameters
  std::string ntupleName_;

//--- internal data-members for handling ntuplees
  std::vector<ntupleEntryType*> ntupleEntries_;
  std::vector<ntupleVectorEntryType*> ntupleVectorEntries_;

  long numEvents_processed_;
};

#endif  


