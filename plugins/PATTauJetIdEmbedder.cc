#include "TauAnalysis/TauIdEfficiency/plugins/PATTauJetIdEmbedder.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"

PATTauJetIdEmbedder::PATTauJetIdEmbedder(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("src");
  
  srcJet_ = cfg.getParameter<edm::InputTag>("srcJet");

  jetIds_.push_back(new jetIdType("jetIdMinimal"));
  jetIds_.push_back(new jetIdType("jetIdLoose"));
  jetIds_.push_back(new jetIdType("jetIdLoose_AOD"));
  //jetIds_.push_back(new jetIdType("jetIdMedium"));
  jetIds_.push_back(new jetIdType("jetIdTight"));

  produces<pat::TauCollection>();
}

PATTauJetIdEmbedder::~PATTauJetIdEmbedder()
{
  for ( std::vector<jetIdType*>::iterator it = jetIds_.begin();
	it != jetIds_.end(); ++it ) {
    delete (*it);
  }
}

void PATTauJetIdEmbedder::produce(edm::Event& evt, const edm::EventSetup& es)
{
  typedef edm::View<pat::Tau> patTauCollectionType;
  edm::Handle<patTauCollectionType> inputTaus;
  evt.getByLabel(src_, inputTaus);

  edm::Handle<pat::JetCollection> patJets;
  evt.getByLabel(srcJet_, patJets);

  std::auto_ptr<pat::TauCollection> outputTaus(new pat::TauCollection() );
  outputTaus->reserve(inputTaus->size());

  for ( patTauCollectionType::const_iterator inputTau = inputTaus->begin(); 
	inputTau != inputTaus->end(); ++inputTau ) {
    pat::Tau outputTau(*inputTau);

    const pat::Jet* patJet = getJet_Tau(*inputTau, *patJets);
    
    for ( std::vector<jetIdType*>::const_iterator jetId = jetIds_.begin();
	  jetId != jetIds_.end(); ++jetId ) {
      float jetId_passed = -1.;
      if ( patJet ) {
	if ( inputTau->isCaloTau() ) {
	  pat::strbitset bits = (*jetId)->caloJetId_->getBitTemplate();
	  jetId_passed = (*(*jetId)->caloJetId_)(*patJet, bits);
	} else if ( inputTau->isPFTau() ) {
	  pat::strbitset bits = (*jetId)->pfJetId_->getBitTemplate();
	  jetId_passed = (*(*jetId)->pfJetId_)(*patJet, bits);
	}	
      }
      outputTau.addUserFloat((*jetId)->value_, jetId_passed);
    }

    outputTaus->push_back(outputTau);
  }

  evt.put(outputTaus);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATTauJetIdEmbedder);

