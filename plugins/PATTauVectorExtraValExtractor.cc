#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorExtraValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "Math/GenVector/VectorUtil.h"

#include <string>

PATTauVectorExtraValExtractor::PATTauVectorExtraValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  pfCandsSrc_ = cfg.getParameter<edm::InputTag>("pfCandsSrc");
  jetSrc_ = cfg.getParameter<edm::InputTag>("jetSrc");
  jetMinPt_ = cfg.getParameter<double>("jetMinPt");
  jetMaxAbsEta_ = cfg.getParameter<double>("jetMaxAbsEta");
  
  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "nTracksOut"      ) value_ = kTracksOut;
  else if ( value_string == "nChargedHadrOut" ) value_ = kChargedHadrOut;
  else if ( value_string == "nPhotonsOut"     ) value_ = kPhotonsOut;
  else if ( value_string == "ClosestJetDR"    ) value_ = kClosestJetDeltaR;
  else if ( value_string == "ClosestJetPt"    ) value_ = kClosestJetPt;
  else if ( value_string == "ClosestJetEta"   ) value_ = kClosestJetEta;
  else if ( value_string == "ClosestJetPhi"   ) value_ = kClosestJetPhi;
  else if ( value_string == "ClosestJetJWidth") value_ = kClosestJetJetWidth;
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
  evt.getByLabel(pfCandsSrc_, pfCandidates);

  std::vector<reco::PFCandidate> PFChargedCands;
  std::vector<reco::PFCandidate> PFChargedHadrCands;
  std::vector<reco::PFCandidate> PFLeptonCands; 
  std::vector<reco::PFCandidate> PFGammaCands; 
    
  for(reco::PFCandidateCollection::const_iterator iCand = pfCandidates->begin(); iCand != pfCandidates->end(); iCand++){
    if((*iCand).particleId() == 1){PFChargedCands.push_back(*iCand); PFChargedHadrCands.push_back(*iCand);}
    if((*iCand).particleId() == 2 || (*iCand).particleId() == 2){PFChargedCands.push_back(*iCand); PFLeptonCands.push_back(*iCand);}
    if((*iCand).particleId() == 4)PFGammaCands.push_back(*iCand);
  }

  edm::Handle<reco::PFJetCollection> jets;
  evt.getByLabel(jetSrc_, jets); 

  unsigned numPatTaus = patTaus->size();
  for ( unsigned i = 0; i < numPatTaus; ++i ) {
    edm::Ptr<pat::Tau> patTauPtr = patTaus->ptrAt(i);
    double vec_i = -1.;
    
    if(      value_ == kTracksOut      ) vec_i = PFChargedCands.size() - patTauPtr->signalTracks().size() - patTauPtr->isolationTracks().size();
    else if (value_ == kChargedHadrOut ) vec_i = PFChargedCands.size() - patTauPtr->signalPFChargedHadrCands().size() - patTauPtr->isolationPFChargedHadrCands().size();
    else if (value_ == kPhotonsOut     ) vec_i = PFChargedCands.size() - patTauPtr->signalPFGammaCands().size() - patTauPtr->isolationPFGammaCands().size();

    if(value_ == kClosestJetDeltaR || value_ == kClosestJetPt || value_ == kClosestJetEta || value_ == kClosestJetPhi || value_ == kClosestJetJetWidth){
      float MinDR = 99.;
      int jClosest = -1;
      for (size_t j = 0; j < jets->size(); ++j) {
	float DRtemp = ROOT::Math::VectorUtil::DeltaR((*jets)[j].p4(),patTauPtr->p4());
	if(DRtemp < MinDR && DRtemp > 0.5 && (*jets)[j].p4().Pt() >= jetMinPt_ && fabs((*jets)[j].p4().Eta()) <= jetMaxAbsEta_){
	  jClosest = (int)j;
	  MinDR = DRtemp;
	}
      }
      if(jClosest != -1){
	if(value_ == kClosestJetDeltaR       ) vec_i = MinDR;
	else if(value_ == kClosestJetPt      ) vec_i = (*jets)[jClosest].p4().Pt();
	else if(value_ == kClosestJetEta     ) vec_i = (*jets)[jClosest].p4().Eta();
	else if(value_ == kClosestJetPhi     ) vec_i = (*jets)[jClosest].p4().Phi();
	else if(value_ == kClosestJetJetWidth) vec_i = sqrt((*jets)[jClosest].etaetaMoment() + (*jets)[jClosest].phiphiMoment());
      }
    }

    vec.push_back(vec_i);
  }

  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorExtraValExtractor, "PATTauVectorExtraValExtractor");
