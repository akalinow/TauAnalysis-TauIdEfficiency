#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorExtraValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "Math/GenVector/VectorUtil.h"

#include <string>

PATTauVectorExtraValExtractor::PATTauVectorExtraValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  pfCandSrc_ = cfg.getParameter<edm::InputTag>("pfCandSrc");
  jetSrc_ = cfg.getParameter<edm::InputTag>("jetSrc");
  jetMinPt_ = cfg.getParameter<double>("jetMinPt");
  jetMaxAbsEta_ = cfg.getParameter<double>("jetMaxAbsEta");
  
  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "numTracksOut"      ) value_ = kNumTracksOut;
  else if ( value_string == "numChargedHadrOut" ) value_ = kNumChargedHadrOut;
  else if ( value_string == "numPhotonsOut"     ) value_ = kNumPhotonsOut;
  else if ( value_string == "nearestJetDR"      ) value_ = kNearestJetDeltaR;
  else if ( value_string == "nearestJetPt"      ) value_ = kNearestJetPt;
  else if ( value_string == "nearestJetEta"     ) value_ = kNearestJetEta;
  else if ( value_string == "nearestJetPhi"     ) value_ = kNearestJetPhi;
  else if ( value_string == "nearestJetWidth"   ) value_ = kNearestJetWidth;
  else {
    edm::LogError ("PATTauVectorExtraValExtractor") << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

PATTauVectorExtraValExtractor::~PATTauVectorExtraValExtractor()
{
//--- nothing to be done yet...
}

std::vector<double> PATTauVectorExtraValExtractor::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;
  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> patTaus;
  evt.getByLabel(src_, patTaus);

  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  evt.getByLabel(pfCandSrc_, pfCandidates);

  std::vector<reco::PFCandidate> pfChargedCands;
  std::vector<reco::PFCandidate> pfChargedHadrCands;
  std::vector<reco::PFCandidate> pfLeptonCands; 
  std::vector<reco::PFCandidate> pfGammaCands; 
    
  for ( reco::PFCandidateCollection::const_iterator pfCandidate = pfCandidates->begin(); 
	pfCandidate != pfCandidates->end(); ++pfCandidate ) {
    switch ( pfCandidate->particleId() ) {
    case reco::PFCandidate::h :
      pfChargedCands.push_back(*pfCandidate); 
      pfChargedHadrCands.push_back(*pfCandidate);
      break;
    case reco::PFCandidate::e  :
    case reco::PFCandidate::mu :
      pfChargedCands.push_back(*pfCandidate); 
      pfLeptonCands.push_back(*pfCandidate);
      break;
    case reco::PFCandidate::gamma :
      pfGammaCands.push_back(*pfCandidate);
      break;
    default :
      break;
    }
  }

  edm::Handle<reco::PFJetCollection> jets;
  evt.getByLabel(jetSrc_, jets); 

  unsigned numPatTaus = patTaus->size();
  for ( unsigned iTau = 0; iTau < numPatTaus; ++iTau ) {
    edm::Ptr<pat::Tau> patTauPtr = patTaus->ptrAt(iTau);
    double vec_i = -1.;
    
    if      ( value_ == kNumTracksOut      ) 
      vec_i = pfChargedCands.size() - patTauPtr->signalTracks().size() - patTauPtr->isolationTracks().size();
    else if ( value_ == kNumChargedHadrOut ) 
      vec_i = pfChargedHadrCands.size() - patTauPtr->signalPFChargedHadrCands().size() - patTauPtr->isolationPFChargedHadrCands().size();
    else if ( value_ == kNumPhotonsOut     ) 
      vec_i = pfGammaCands.size() - patTauPtr->signalPFGammaCands().size() - patTauPtr->isolationPFGammaCands().size();

    if ( value_ == kNearestJetDeltaR || 
	 value_ == kNearestJetPt     || 
	 value_ == kNearestJetEta    || 
	 value_ == kNearestJetPhi    || 
	 value_ == kNearestJetWidth ){
      float dRmin = 99.;
      int nearestJet_index = -1;
      for ( size_t iJet = 0; iJet < jets->size(); ++iJet ) {
	float dR = ROOT::Math::VectorUtil::DeltaR((*jets)[iJet].p4(), patTauPtr->p4());
	if ( dR > 0.5 && dR < dRmin && (*jets)[iJet].p4().Pt() >= jetMinPt_ && fabs((*jets)[iJet].p4().Eta()) <= jetMaxAbsEta_ ){
	  nearestJet_index = (int)iJet;
	  dRmin = dR;
	}
      }
      if ( nearestJet_index != -1 ){
	if      ( value_ == kNearestJetDeltaR   ) vec_i = dRmin;
	else if ( value_ == kNearestJetPt       ) vec_i = (*jets)[nearestJet_index].p4().Pt();
	else if ( value_ == kNearestJetEta      ) vec_i = (*jets)[nearestJet_index].p4().Eta();
	else if ( value_ == kNearestJetPhi      ) vec_i = (*jets)[nearestJet_index].p4().Phi();
	else if ( value_ == kNearestJetWidth ) 
	  vec_i = sqrt((*jets)[nearestJet_index].etaetaMoment() + (*jets)[nearestJet_index].phiphiMoment());
      }
    }

    vec.push_back(vec_i);
  }

  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorExtraValExtractor, "PATTauVectorExtraValExtractor");
