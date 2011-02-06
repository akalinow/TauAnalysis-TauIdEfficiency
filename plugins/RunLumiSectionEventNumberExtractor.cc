#include "TauAnalysis/TauIdEfficiency/plugins/RunLumiSectionEventNumberExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

RunLumiSectionEventNumberExtractor::RunLumiSectionEventNumberExtractor(const edm::ParameterSet& cfg)
  : cfgError_(0)
{
  //std::cout << "<RunLumiSectionEventNumberExtractor::RunLumiSectionEventNumberExtractor>:" << std::endl;

  std::string value_string = cfg.getParameter<std::string>("value");
  //std::cout << " value_string = " << value_string << std::endl;

  if      ( value_string == "run"   ) value_ = kRunNumber;
  else if ( value_string == "ls"    ) value_ = kLumiSection;
  else if ( value_string == "event" ) value_ = kEventNumber;
  else {
    edm::LogError("RunLumiSectionEventNumberExtractor") 
      << " Configuration parameter value = " << value_string  << " invalid !!";
    cfgError_ = 1;
  }

  maxWarnings_ = cfg.exists("maxWarnings") ? cfg.getParameter<int>("maxWarnings") : 1;
  numWarnings_ = 0;
}

RunLumiSectionEventNumberExtractor::~RunLumiSectionEventNumberExtractor()
{}

double  
RunLumiSectionEventNumberExtractor::operator()(const edm::Event& evt) const
{
  if ( cfgError_ ) {
    if ( numWarnings_ < maxWarnings_ ) 
      edm::LogError("RunLumiSectionEventNumberExtractor::operator()") 
	<< " Error in Configuration ParameterSet --> skipping !!";
    ++numWarnings_;
    return -1;
  }

  double val = -1.;

  if      ( value_ == kRunNumber   ) val = evt.id().run();
  else if ( value_ == kLumiSection ) val = evt.luminosityBlock();
  else if ( value_ == kEventNumber ) val = evt.id().event();
  else assert(0);
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, RunLumiSectionEventNumberExtractor, "RunLumiSectionEventNumberExtractor");
