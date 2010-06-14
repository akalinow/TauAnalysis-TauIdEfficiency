#include "TauAnalysis/TauIdEfficiency/plugins/HLTInfoExtractor.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

HLTInfoExtractor::HLTInfoExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src"); 
  value_ = cfg.getParameter<std::string>("value");

  maxWarnings_ = cfg.exists("maxWarnings") ? cfg.getParameter<int>("maxWarnings") : 1;
  numWarnings_ = 0;
}

HLTInfoExtractor::~HLTInfoExtractor()
{}

double  
HLTInfoExtractor::operator()(const edm::Event& evt) const
{
  edm::Handle<edm::TriggerResults> hltResults;
  evt.getByLabel(src_, hltResults);

  // get the names of HLT trigger paths
  edm::TriggerNames const& triggerNames = evt.triggerNames(*hltResults);

  // find index of trigger path specified by "value" configuration parameter;
  // print error message and return -1 in case path not found in trigger menu
  unsigned int triggerId = triggerNames.triggerIndex(value_);

  double val = -1.;
  if ( triggerId >= 0 && triggerId < triggerNames.size() ) {
    if ( hltResults->accept(triggerId) ) {
      val = 1.0;
    } else {
      val = 0.0;
    }
  } else {
    if ( numWarnings_ < maxWarnings_ || maxWarnings_ == -1 ) {
      edm::LogError("HLTInfoExtractor::operator()") 
	<< " Trigger path = " << value_ << " not found in Trigger menu --> skipping !!";

      std::cout << "Trigger paths defined in menu: " << std::endl;
      for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
	    triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
	unsigned int triggerId = triggerNames.triggerIndex(*triggerName);
	if ( triggerId >= 0 && triggerId < triggerNames.size() ) {
	  std::string triggerDecision = ( hltResults->accept(triggerId) ) ? "passed" : "failed";
	  
	  std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
	}
      }

      ++numWarnings_;
    }
  }
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, HLTInfoExtractor, "HLTInfoExtractor");
