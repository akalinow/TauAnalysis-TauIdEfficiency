#include "TauAnalysis/TauIdEfficiency/plugins/GenPhaseSpaceEventInfoExtractor.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

GenPhaseSpaceEventInfoExtractor::GenPhaseSpaceEventInfoExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src"); 

  std::string value_string = cfg.getParameter<std::string>("value");
  if ( value_string == "ptHat" ) value_ = kPtHat;
  else {
    edm::LogError ("GenPhaseSpaceEventInfoExtractor") 
      << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

GenPhaseSpaceEventInfoExtractor::~GenPhaseSpaceEventInfoExtractor()
{}

double  
GenPhaseSpaceEventInfoExtractor::operator()(const edm::Event& evt) const
{
  //std::cout << "<GenPhaseSpaceEventInfoExtractor::operator()>:" << std::endl;

  edm::Handle<GenEventInfoProduct> genEventInfo;
  evt.getByLabel(src_, genEventInfo);

  double val = -1.;
  if ( genEventInfo.isValid() && genEventInfo->hasBinningValues() ) {
    //std::cout << " ptHat = " << genEventInfo->binningValues()[0] << std::endl;
    val = genEventInfo->binningValues()[0];
  }
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, GenPhaseSpaceEventInfoExtractor, "GenPhaseSpaceEventInfoExtractor");
