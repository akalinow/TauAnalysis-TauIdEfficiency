#ifndef TauAnalysis_TauIdEfficiency_TauIdTagAndProbeProducer_h
#define TauAnalysis_TauIdEfficiency_TauIdTagAndProbeProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"

#include <vector>
#include <string>

class TauIdTagAndProbeProducer : public edm::EDProducer 
{  
 public:

  explicit TauIdTagAndProbeProducer(const edm::ParameterSet&);
  virtual ~TauIdTagAndProbeProducer();
  void produce(edm::Event&, const edm::EventSetup&);

 private:
  
  edm::InputTag src_;

  typedef StringCutObjectSelector<pat::TriggerObject> StringCutTriggerObjectSelector;
  typedef std::vector<StringCutTriggerObjectSelector*> vStringCutTriggerObjectSelector;
  std::map<std::string, vStringCutTriggerObjectSelector> triggerPaths_;
};

#endif
