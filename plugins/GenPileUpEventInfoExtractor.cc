#include "TauAnalysis/TauIdEfficiency/plugins/GenPileUpEventInfoExtractor.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

GenPileUpEventInfoExtractor::GenPileUpEventInfoExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src"); 

  std::string value_string = cfg.getParameter<std::string>("value");
  if ( value_string == "numPileUpInteractions" ) value_ = kNumPileUpInteractions;
  else {
    edm::LogError ("GenPileUpEventInfoExtractor") 
      << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

GenPileUpEventInfoExtractor::~GenPileUpEventInfoExtractor()
{}

double  
GenPileUpEventInfoExtractor::operator()(const edm::Event& evt) const
{
  //std::cout << "<GenPileUpEventInfoExtractor::operator()>:" << std::endl;

  typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
  edm::Handle<PileupSummaryInfoCollection> genPileUpInfo;
  evt.getByLabel(src_, genPileUpInfo);

  double val = -1.;
  if ( genPileUpInfo.isValid() && genPileUpInfo->size() >= 1 ) {
    val = genPileUpInfo->front().getPU_NumInteractions();
    //std::cout << " numPileUpInteractions = " << val << std::endl;
  }
  
  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, GenPileUpEventInfoExtractor, "GenPileUpEventInfoExtractor");
