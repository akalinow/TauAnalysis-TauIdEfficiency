#include "TauAnalysis/TauIdEfficiency/plugins/NumSelObjExtractor.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <string>

template<typename T>
NumSelObjExtractor<T>::NumSelObjExtractor(const edm::ParameterSet& cfg)
  : cut_(0)
{
  src_ = cfg.getParameter<edm::InputTag>("src");

  if ( cfg.exists("srcNotToBeFiltered") ) {
    srcNotToBeFiltered_ = cfg.getParameter<vInputTag>("srcNotToBeFiltered");
    dRmin_ = cfg.getParameter<double>("dRmin");
  }

  if ( cfg.exists("value") ) {
    std::string cut_string = cfg.getParameter<std::string>("value");
    cut_ = new StringCutObjectSelector<T>(cut_string);
  }
}

template<typename T>
NumSelObjExtractor<T>::~NumSelObjExtractor()
{
  delete cut_;
}

template<typename T>
double NumSelObjExtractor<T>::operator()(const edm::Event& evt) const
{
  unsigned numSelObjects = 0;
  
  typedef edm::View<T> patCollectionType;
  edm::Handle<patCollectionType> patObjects;
  evt.getByLabel(src_, patObjects);

  for ( typename patCollectionType::const_iterator patObject = patObjects->begin();
	patObject != patObjects->end(); ++patObject ) {

    if ( cut_ && (*cut_)(*patObject) == false ) continue;

    bool isOverlap = false;    
    for ( vInputTag::const_iterator srcOverlap = srcNotToBeFiltered_.begin();
	  srcOverlap != srcNotToBeFiltered_.end(); ++srcOverlap ) {

      typedef edm::View<reco::Candidate> overlapCollectionType;
      edm::Handle<overlapCollectionType> overlapObjects;
      evt.getByLabel(*srcOverlap, overlapObjects);

      for ( overlapCollectionType::const_iterator overlapObject = overlapObjects->begin();
	    overlapObject != overlapObjects->end(); ++overlapObject ) {
	if ( deltaR(patObject->p4(), overlapObject->p4()) < dRmin_ ) isOverlap = true;
      }
    }

    if ( isOverlap ) continue;

    ++numSelObjects;
  }

  return numSelObjects;
}

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

typedef NumSelObjExtractor<reco::Candidate> NumSelCandidateExtractor;
typedef NumSelObjExtractor<pat::Electron> NumSelPATElectronExtractor;
typedef NumSelObjExtractor<pat::Muon> NumSelPATMuonExtractor;
typedef NumSelObjExtractor<pat::Tau> NumSelPATTauExtractor;
typedef NumSelObjExtractor<pat::Jet> NumSelPATJetExtractor;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, NumSelCandidateExtractor, "NumSelCandidateExtractor");
DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, NumSelPATElectronExtractor, "NumSelPATElectronExtractor");
DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, NumSelPATMuonExtractor, "NumSelPATMuonExtractor");
DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, NumSelPATTauExtractor, "NumSelPATTauExtractor");
DEFINE_EDM_PLUGIN(ObjValExtractorPluginFactory, NumSelPATJetExtractor, "NumSelPATJetExtractor");
