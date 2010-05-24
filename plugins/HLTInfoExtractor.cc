#include "TauAnalysis/TauIdEfficiency/plugins/HLTInfoExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <string>

HLTInfoExtractor::HLTInfoExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src"); 
  value_ = cfg.getParameter<std::string>("value") ;
}

HLTInfoExtractor::~HLTInfoExtractor()
{}

double  
HLTInfoExtractor::operator()(const edm::Event& evt) const
{
  edm::Handle<edm::TriggerResults> hltResults;
  evt.getByLabel(src_, hltResults);

  //get the names of the triggers
  edm::TriggerNames const& triggerNames = evt.triggerNames(*hltResults);
  unsigned int triggerId = triggerNames.triggerIndex(value_);
  double val = -1.; // return -1 if error or path not in the menu
  if ( triggerId != triggerNames.size() ) {
    if ( hltResults->accept(triggerId) ) {
      val = 1.0;
    } else {
      val = 0.0;
    }
  }
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, HLTInfoExtractor, "HLTInfoExtractor");
