#include "TauAnalysis/TauIdEfficiency/plugins/VertexValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <string>

VertexValExtractor::VertexValExtractor(const edm::ParameterSet& cfg)
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

VertexValExtractor::~VertexValExtractor()
{
//--- nothing to be done yet...
}

double VertexValExtractor::operator()(const edm::Event& evt) const
{
  double val = 0.;

  edm::Handle<reco::VertexCollection> vertices;
  evt.getByLabel(srcVertex_, vertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  evt.getByLabel(srcBeamSpot_, beamSpot);

//--- take as "the" event vertex the vertex with the highest Pt sum of tracks associated to it;
//    this vertex is the first entry in the collection of (offline) reconstructed vertices
  reco::VertexCollection::const_iterator vertex = vertices->begin();
  if ( vertex != vertices->end() ) {
    
    if      ( value_ == kVertexX ) val = vertex->x() - beamSpot->x0();
    else if ( value_ == kVertexY ) val = vertex->y() - beamSpot->y0();
    else if ( value_ == kVertexZ ) val = vertex->z();
  }

  return val;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, VertexValExtractor, "VertexValExtractor");
