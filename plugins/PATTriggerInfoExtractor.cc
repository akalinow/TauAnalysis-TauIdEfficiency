#include "TauAnalysis/TauIdEfficiency/plugins/PATTriggerInfoExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

const std::string value_separator = ":";

PATTriggerInfoExtractor::PATTriggerInfoExtractor(const edm::ParameterSet& cfg)
  : cfgError_(0)
{
  //std::cout << "<PATTriggerInfoExtractor::PATTriggerInfoExtractor>:" << std::endl;
  
  src_ = cfg.getParameter<edm::InputTag>("src"); 

  std::string encodedValue_string = cfg.getParameter<std::string>("value");
  //std::cout << " encodedValue_string = " << encodedValue_string << std::endl;

//--- check that value contains **exactly once** the characted ":",
//    which separates the name of the HLT trigger path from the name of the value (bit, prescale,...)
//    to be extracted for that trigger path
  if ( encodedValue_string.find(value_separator) != std::string::npos                          &&
       encodedValue_string.find(value_separator) == encodedValue_string.rfind(value_separator) ) {
    size_t index_separator = encodedValue_string.find(value_separator);

    hltPath_ = std::string(encodedValue_string, 0, index_separator);
    //std::cout << " hltPath = " << hltPath_ << std::endl;

    std::string decodedValue_string = std::string(encodedValue_string, index_separator + value_separator.length());
    //std::cout << " decodedValue_string = " << decodedValue_string << std::endl;
    if      ( decodedValue_string == "bit"      ) value_ = kBit;
    else if ( decodedValue_string == "prescale" ) value_ = kPrescale;
    else {
      edm::LogError("PATTriggerInfoExtractor") 
	<< " Configuration parameter value = " << encodedValue_string  << " invalid !!";
      cfgError_ = 1;
    }
  } else {
    edm::LogError("PATTriggerInfoExtractor") 
      << " Failed to decode configuration parameter value = " << encodedValue_string  << " !!";
    cfgError_ = 1;
  }

  maxWarnings_ = cfg.exists("maxWarnings") ? cfg.getParameter<int>("maxWarnings") : 1;
  numWarnings_ = 0;
}

PATTriggerInfoExtractor::~PATTriggerInfoExtractor()
{}

double  
PATTriggerInfoExtractor::operator()(const edm::Event& evt) const
{
  if ( cfgError_ ) {
    if ( numWarnings_ < maxWarnings_ ) edm::LogError("PATTriggerInfoExtractor::operator()") 
                                         << " Error in Configuration ParameterSet --> skipping !!";
    ++numWarnings_;
    return -1;
  }

  double val = -1.;

  edm::Handle<pat::TriggerEvent> patTriggerEvent;
  evt.getByLabel(src_, patTriggerEvent);

//--- get pat::TriggerAlgorithm object
//    containing information for specified HLT path
  const pat::TriggerAlgorithm* algorithm = patTriggerEvent->algorithm(hltPath_);
  if ( algorithm ) {
    if      ( value_ == kBit      ) val = algorithm->bit();
    else if ( value_ == kPrescale ) val = algorithm->prescale();
    else assert(0);
  } else {
    if ( numWarnings_ < maxWarnings_ ) {
      edm::LogError("PATTriggerInfoExtractor::operator()") 
	<< " Trigger path = " << hltPath_ << " not found in Trigger menu --> skipping !!";
      std::cout << "available algorithms:" << std::endl;
      const pat::TriggerAlgorithmCollection* algorithms = patTriggerEvent->algorithms();
      for ( pat::TriggerAlgorithmCollection::const_iterator algorithm = algorithms->begin();
	    algorithm != algorithms->end(); ++algorithm ) {
	std::cout << " " << algorithm->name() << " (alias = " << algorithm->alias() << ")" << std::endl;
      }
    }
    ++numWarnings_;
  }
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, PATTriggerInfoExtractor, "PATTriggerInfoExtractor");
