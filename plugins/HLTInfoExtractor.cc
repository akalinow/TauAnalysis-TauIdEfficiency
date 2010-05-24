#include "TauAnalysis/TauIdEfficiency/plugins/HLTInfoExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <string>

HLTInfoExtractor::HLTInfoExtractor(const edm::ParameterSet& cfg)
{

  srcTrigger_  = cfg.getParameter<edm::InputTag>("triggerResults"); 
  valueString_ = cfg.getParameter<std::string>("HLTPath") ;
}

HLTInfoExtractor::~HLTInfoExtractor()
{

}

double  
HLTInfoExtractor::operator()(const edm::Event& evt) const
{
  edm::Handle<edm::TriggerResults> trigResults;
  evt.getByLabel(srcTrigger_, trigResults);

  //get the names of the triggers
  edm::TriggerNames const& triggerNames = evt.triggerNames(*trigResults);
  unsigned int triggerId = triggerNames.triggerIndex(valueString_);
  double value=-1;//Return -1 if error or path not in the menu
  if(triggerId!=triggerNames.size()) {
    if(trigResults->accept(triggerId)) {
      value= 1.0;
    }
    else {
      value = 0.0;
    }
  }

  return value;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, HLTInfoExtractor, "HLTInfoExtractor");
