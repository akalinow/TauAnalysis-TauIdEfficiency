#include "TauAnalysis/TauIdEfficiency/plugins/VertexVectorValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <string>

VertexVectorValExtractor::VertexVectorValExtractor(const edm::ParameterSet& cfg)
{
  srcVertex_ = cfg.getParameter<edm::InputTag>("src");
  srcBeamSpot_ = cfg.getParameter<edm::InputTag>("srcBeamSpot");
  
  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "vertexX" ) value_ = kVertexX;
  else if ( value_string == "vertexY" ) value_ = kVertexY;
  else if ( value_string == "vertexZ" ) value_ = kVertexZ;
  else {
    edm::LogError ("VertexValExtractor") << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

VertexVectorValExtractor::~VertexVectorValExtractor()
{
//--- nothing to be done yet...
}

std::vector<double> VertexVectorValExtractor::operator()(const edm::Event& evt) const
{
   std::vector<double> vec;

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(srcVertex_, vertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_, beamSpot);

  for ( reco::VertexCollection::const_iterator vertex = vertices->begin();
	vertex != vertices->end(); ++vertex ) {
    
    double vec_i = -1.;
    
    if      ( value_ == kVertexX ) vec_i = vertex->x() - beamSpot->x0();
    else if ( value_ == kVertexY ) vec_i = vertex->y() - beamSpot->y0();
    else if ( value_ == kVertexZ ) vec_i = vertex->z();
    
    vec.push_back(vec_i);
  }
  
  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, VertexVectorValExtractor, "VertexVectorValExtractor");
