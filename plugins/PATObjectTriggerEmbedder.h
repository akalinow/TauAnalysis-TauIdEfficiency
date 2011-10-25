#ifndef TauAnalysis_TauIdEfficiency_PATObjectTriggerEmbedder_h
#define TauAnalysis_TauIdEfficiency_PATObjectTriggerEmbedder_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <string>

template <typename T>
class PATObjectTriggerEmbedder : public edm::EDProducer 
{  
  typedef std::vector<T> PATObjectCollection;

 public:

  explicit PATObjectTriggerEmbedder(const edm::ParameterSet&);
  virtual ~PATObjectTriggerEmbedder();
  void produce(edm::Event&, const edm::EventSetup&);

 private:
  
  edm::InputTag src_;

  edm::InputTag matched_;
  
  std::map<std::string, std::string> triggerPathsAndLabels_;

  double dRmax_;
};

#endif
