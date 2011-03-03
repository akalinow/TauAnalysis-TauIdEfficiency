#include "TauAnalysis/TauIdEfficiency/plugins/VertexValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <string>

const std::string trackSumPtThreshold_keyword = "numVerticesPtGt";

VertexValExtractor::VertexValExtractor(const edm::ParameterSet& cfg)
  : error_(0)
{
  srcVertex_ = cfg.getParameter<edm::InputTag>("src");
  
  std::string value_string = cfg.getParameter<std::string>("value");
  if ( value_string.find(trackSumPtThreshold_keyword) != std::string::npos ) {
    std::string trackSumPtThreshold_string = 
      std::string(value_string, trackSumPtThreshold_keyword.length());
    trackSumPtThreshold_string = replace_string(trackSumPtThreshold_string, "_", ".", 0, 1000, error_);
    assert(error_ == 0);
    trackSumPtThreshold_ = atof(trackSumPtThreshold_string.data());
  } else {
    edm::LogError ("VertexValExtractor") << " Invalid configuration parameter value = " << value_string << " !!";
    error_ = 1;
  }
}

VertexValExtractor::~VertexValExtractor()
{
//--- nothing to be done yet...
}

double VertexValExtractor::operator()(const edm::Event& evt) const
{
  double val = 0.;

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(srcVertex_, vertices);

  if ( !error_ ) {
    std::vector<double> trackPtSums = compTrackPtSums(*vertices);
    val = getNumVerticesPtGtThreshold(trackPtSums, trackSumPtThreshold_);
  }

  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, VertexValExtractor, "VertexValExtractor");
