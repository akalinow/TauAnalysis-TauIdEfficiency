#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorTrackValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include <string>

PATTauVectorTrackValExtractor::PATTauVectorTrackValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  std::string collection_string = cfg.getParameter<std::string>("collection");
  if      ( collection_string == "leadTrack"           ) collection_ = kLeadTrack;
  else if ( collection_string == "signalConeTracks"    ) collection_ = kSignalConeTracks;
  else if ( collection_string == "isolationConeTracks" ) collection_ = kIsolationConeTracks;
  else {
    edm::LogError ("PATTauVectorTrackValExtractor") << " Invalid configuration parameter collection = " << collection_string << " !!";
    collection_ = -1;
  }

  if ( collection_ == kSignalConeTracks || collection_ == kIsolationConeTracks ) {
    index_ = cfg.getParameter<unsigned>("index");
  }

  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "numValidHits"   ) value_ = kNumValidHits;
  else if ( value_string == "numMissingHits" ) value_ = kNumMissingHits;
  else if ( value_string == "chi2"           ) value_ = kChi2;
  else if ( value_string == "nDoF"           ) value_ = kNumDoF;
  else if ( value_string == "dz"             ) value_ = kDeltaZ;
  else if ( value_string == "dxy"            ) value_ = kDeltaXY;
  else if ( value_string == "pt"             ) value_ = kPt;
  else if ( value_string == "ptErr"          ) value_ = kPtErr;
  else if ( value_string == "quality"        ) value_ = kQuality;
  else {
    edm::LogError ("PATTauVectorTrackValExtractor") 
      << " Invalid Configuration Parameter 'value' = " << value_string << " !!";
    value_ = -1;
  }

  if ( value_ == kDeltaZ || value_ == kDeltaXY ) {
    srcVertex_ = cfg.getParameter<edm::InputTag>("vertexSrc");
  }
}

const reco::Track* PATTauVectorTrackValExtractor::getTrack_pfTau(const pat::Tau& tau) const
{
  const reco::Track* track = 0;

  if ( collection_ == kLeadTrack ) {
    if ( tau.leadPFChargedHadrCand().isNonnull() ) track = tau.leadPFChargedHadrCand()->trackRef().get();
  } else if ( collection_ == kSignalConeTracks ) {
    const std::vector<reco::PFCandidatePtr>& tauSignalTracks = tau.signalPFChargedHadrCands();
    if ( index_ < tauSignalTracks.size() ) track = tauSignalTracks[index_]->trackRef().get();
  } else if ( collection_ == kSignalConeTracks ) {
    const std::vector<reco::PFCandidatePtr>& tauIsolationTracks = tau.isolationPFChargedHadrCands();
    if ( index_ < tauIsolationTracks.size() ) track = tauIsolationTracks[index_]->trackRef().get();
  } 

  return track;
}

const reco::Track* PATTauVectorTrackValExtractor::getTrack_caloTau(const pat::Tau& tau) const
{
  const reco::Track* track = 0;

  if ( collection_ == kLeadTrack ) {
    if ( tau.leadTrack().isNonnull() ) track = tau.leadTrack().get();
  } else if ( collection_ == kSignalConeTracks ) {
    const reco::TrackRefVector& tauSignalTracks = tau.signalTracks();
    if ( index_ < tauSignalTracks.size() ) track = tauSignalTracks[index_].get();
  } else if ( collection_ == kSignalConeTracks ) {
    const reco::TrackRefVector& tauIsolationTracks = tau.isolationTracks();
    if ( index_ < tauIsolationTracks.size() ) track = tauIsolationTracks[index_].get();
  } 

  return track;
}

std::vector<double> PATTauVectorTrackValExtractor::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;

  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> patTaus;
  evt.getByLabel(src_, patTaus);

  reco::Vertex::Point thePrimaryEventVertexPosition;
  if ( value_ == kDeltaZ ||  value_ == kDeltaXY ) {
    edm::Handle<reco::VertexCollection> primaryEventVertices;
    evt.getByLabel(srcVertex_, primaryEventVertices);
    if ( primaryEventVertices->size() > 0 ) thePrimaryEventVertexPosition = primaryEventVertices->begin()->position();
  }

  unsigned numPatTaus = patTaus->size();
  for ( unsigned iTau = 0; iTau < numPatTaus; ++iTau ) {
    edm::Ptr<pat::Tau> patTauPtr = patTaus->ptrAt(iTau);

    double vec_i = -1.;

    const reco::Track* track = 0;
    if      ( patTauPtr->isPFTau()   ) track = getTrack_pfTau(*patTauPtr);
    else if ( patTauPtr->isCaloTau() ) track = getTrack_caloTau(*patTauPtr);
    
    if ( track ) {
      if      ( value_ == kNumValidHits   ) vec_i = track->numberOfValidHits();
      else if ( value_ == kNumMissingHits ) vec_i = track->numberOfLostHits();    
      else if ( value_ == kChi2           ) vec_i = track->chi2();
      else if ( value_ == kNumDoF         ) vec_i = track->ndof();   
      else if ( value_ == kDeltaZ         ) vec_i = track->dz(thePrimaryEventVertexPosition);
      else if ( value_ == kDeltaXY        ) vec_i = track->dxy(thePrimaryEventVertexPosition);
      else if ( value_ == kPt             ) vec_i = track->pt();
      else if ( value_ == kPtErr          ) vec_i = track->ptError();
//--- fill track quality
//
//    NOTE: the track quality is bit-coded;
//          the meaning of the bits is:
//            0: loose
//            1: tight
//            2: highPurity
//            3: confirmed
//            4: goodIterative
//          (cf. definition of "TrackQuality" enum in DataFormats/TrackReco/interface/TrackBase.h)
//
      else if ( value_ == kQuality        ) vec_i = track->qualityMask();
    }

    vec.push_back(vec_i);
  }

  return vec;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorTrackValExtractor, "PATTauVectorTrackValExtractor");
