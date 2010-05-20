#include "TauAnalysis/TauIdEfficiency/plugins/PATTauVectorGenJetValExtractor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include <string>

PATTauVectorGenJetValExtractor::PATTauVectorGenJetValExtractor(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  std::string value_string = cfg.getParameter<std::string>("value");
  if      ( value_string == "genMatch"     ) value_ = kGenMatch;
  else if ( value_string == "genPt"        ) value_ = kGenPt;
  else if ( value_string == "genEta"       ) value_ = kGenEta;
  else if ( value_string == "genPhi"       ) value_ = kGenPhi;
  else if ( value_string == "genDecayMode" ) value_ = kGenDecayMode;
  else {
    edm::LogError ("PATTauVectorGenJetValExtractor") << " Invalid configuration parameter value = " << value_string << " !!";
    value_ = -1;
  }
}

PATTauVectorGenJetValExtractor::~PATTauVectorGenJetValExtractor()
{
//--- nothing to be done yet...
}

unsigned int 
PATTauVectorGenJetValExtractor::size() const
{
   return 0;
}

std::vector<double> PATTauVectorGenJetValExtractor::operator()(const edm::Event& evt) const
{
  std::vector<double> vec;
  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> patTaus;
  evt.getByLabel(src_, patTaus);

  unsigned numPatTaus = patTaus->size();
  for ( unsigned i = 0; i < numPatTaus; ++i ) {
    edm::Ptr<pat::Tau> patTauPtr = patTaus->ptrAt(i);

    double vec_i = -1.;

    bool isMatched = ( patTauPtr->genJet() != 0 ) ? true : false;
    
    if      ( value_ == kGenMatch ) vec_i = isMatched;
    else if ( isMatched ) {
      if      ( value_ == kGenPt        ) vec_i = patTauPtr->genJet()->pt();
      else if ( value_ == kGenEta       ) vec_i = patTauPtr->genJet()->eta();
      else if ( value_ == kGenPhi       ) vec_i = patTauPtr->genJet()->phi();
      else if ( value_ == kGenDecayMode ) {	
	std::string genDecayMode_string = JetMCTagUtils::genTauDecayMode(*patTauPtr->genJet());
//--- decode generated tau decay mode
//    ( as defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc )
	if      ( genDecayMode_string == "electron"        ) vec_i = 0;
	else if ( genDecayMode_string == "muon"            ) vec_i = 1;
	else if ( genDecayMode_string == "oneProng0Pi0"    ) vec_i = 2;
	else if ( genDecayMode_string == "oneProng1Pi0"    ) vec_i = 3;
	else if ( genDecayMode_string == "oneProng2Pi0"    ) vec_i = 4;
	else if ( genDecayMode_string == "oneProngOther"   ) vec_i = 5;
	else if ( genDecayMode_string == "threeProng0Pi0"  ) vec_i = 6;
	else if ( genDecayMode_string == "threeProng1Pi0"  ) vec_i = 7;
	else if ( genDecayMode_string == "threeProngOther" ) vec_i = 8;
	else if ( genDecayMode_string == "rare"            ) vec_i = 9;
      }
    } 

    vec.push_back(vec_i);
  }

  return vec;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(ObjValVectorExtractorPluginFactory, PATTauVectorGenJetValExtractor, "PATTauVectorGenJetValExtractor");
