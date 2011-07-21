#include "TauAnalysis/TauIdEfficiency/plugins/TauIdTagAndProbeProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

typedef std::vector<std::string> vstring;
typedef std::vector<pat::Tau> vPatTaus;

TauIdTagAndProbeProducer::TauIdTagAndProbeProducer(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("source");

  edm::ParameterSet cfgTriggerPaths = cfg.getParameter<edm::ParameterSet>("triggerPaths");
  vstring triggerPathsNames = cfgTriggerPaths.getParameterNamesForType<vstring>();
  for ( vstring::const_iterator triggerPathsName = triggerPathsNames.begin();
	triggerPathsName != triggerPathsNames.end(); ++triggerPathsName ) {
    vstring triggerPathSelections = cfgTriggerPaths.getParameter<vstring>(*triggerPathsName);
    for ( vstring::const_iterator triggerPathSelection = triggerPathSelections.begin();
	  triggerPathSelection != triggerPathSelections.end(); ++triggerPathSelection ) {
      triggerPaths_[*triggerPathsName].push_back(new StringCutTriggerObjectSelector(*triggerPathSelection));
    }
  }
  
  // register products
  produces<vPatTaus>();
}

TauIdTagAndProbeProducer::~TauIdTagAndProbeProducer()
{
  for ( std::map<std::string, vStringCutTriggerObjectSelector>::iterator it1 = triggerPaths_.begin();
	it1 != triggerPaths_.end(); ++it1) {
    for ( vStringCutTriggerObjectSelector::iterator it2 = it1->second.begin();
	  it2 != it1->second.end(); ++it2 ) {
      delete (*it2);
    }
  }
}

namespace {
  struct TauInfo 
  {
    TauInfo()
      : tau_(0), 
	tauPt_(0.)
    {}
    TauInfo(const pat::Tau& tau)
      : tau_(&tau)
    {
      // get Pt of jet associated to tau
      if ( tau.isPFTau() ) {
	tauPt_ = tau.pfJetRef()->pt();
      } else if ( tau.isCaloTau() ) {
	tauPt_ = tau.caloTauTagInfoRef()->jetRef()->pt();
      } else {
	throw cms::Exception("TauInfo") 
	  << "pat::Tau is of neither PF nor Calo type !!\n";
      }
    }

    const pat::Tau* tau_;
    double tauPt_;

    std::map<std::string, bool> matchesTriggerObject_;
  };

  // Sorting predicate
  bool tauInfoDescendingPt(const TauInfo& a, const TauInfo& b)
  {
    return ( a.tauPt_ > b.tauPt_ );
  }
}

void TauIdTagAndProbeProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  // output products
  std::auto_ptr<vPatTaus> output(new vPatTaus());
  
  edm::Handle<edm::View<pat::Tau> > sourceView;
  evt.getByLabel(src_, sourceView);
  
  size_t inputSize = sourceView->size();
  std::vector<TauInfo> tauInfos(inputSize);

  // count how many taus are matched to a trigger object
  std::map<std::string, unsigned> nTriggers;
  for ( size_t iTau = 0; iTau < inputSize; ++iTau ){
    // get tau and check if it matches trigger
    const pat::Tau& tau = sourceView->at(iTau);

    TauInfo myTauInfo(tau);

    for ( std::map<std::string, vStringCutTriggerObjectSelector>::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {
      bool matchesSelTriggerObject = false;
      const pat::TriggerObjectStandAlone* matchedTriggerObjectStandAlone = 0;
      for ( pat::TriggerObjectStandAloneCollection::const_iterator matchedTrigger = tau.triggerObjectMatches().begin();
	    matchedTrigger != tau.triggerObjectMatches().end(); ++matchedTrigger ) {
	bool passesSelection = true;
	for ( vStringCutTriggerObjectSelector::const_iterator triggerPathSelection = triggerPath->second.begin();
	      triggerPathSelection != triggerPath->second.end(); ++triggerPathSelection ) {
	  if ( !(**triggerPathSelection)(*matchedTrigger) ) passesSelection = false;
	}

	if ( passesSelection ) {
	  matchesSelTriggerObject = true;
	  if ( matchedTriggerObjectStandAlone == 0                         || 
	       matchedTriggerObjectStandAlone->pt() < matchedTrigger->pt() ) matchedTriggerObjectStandAlone = &(*matchedTrigger);	  
	}
      }

      myTauInfo.matchesTriggerObject_[triggerPath->first] = matchesSelTriggerObject;

      //std::cout << "tau: pt = " << tau.pt() << ", eta = " << tau.eta() << ", phi = " << tau.phi();
      //if ( matchesSelTriggerObject ) 
      //  std::cout << " --> matches " << triggerPath->first << " "
      //	    << "(pt = " << matchedTriggerObjectStandAlone->pt() << "," 
      //	    << " eta = " << matchedTriggerObjectStandAlone->eta() << ","
      //	    << " phi = " << matchedTriggerObjectStandAlone->phi() << ")";
      //else 
      //  std::cout << " --> does not match " << triggerPath->first;
      //std::cout << std::endl;

      if ( matchesSelTriggerObject ) ++nTriggers[triggerPath->first];
    }

    tauInfos[iTau] = myTauInfo;
  }
  
  // sort by descending pt
  sort(tauInfos.begin(), tauInfos.end(), tauInfoDescendingPt);
  
  for ( size_t iTau = 0; iTau < inputSize; ++iTau ) {
    // make our own copy of the tau
    pat::Tau newTau = (*tauInfos[iTau].tau_);

    // store pt-ordered index
    newTau.addUserFloat("pt_index", iTau);
    
    for ( std::map<std::string, vStringCutTriggerObjectSelector>::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {

      //std::cout << "triggerPath = " << triggerPath->first << std::endl;

      // set probe flag
      bool isProbeTau = false;
      if ( nTriggers[triggerPath->first] >= 2 ){
	// At least two triggers, so this is always a probe
	isProbeTau = true;
      } else if ( nTriggers[triggerPath->first] == 1 && !tauInfos[iTau].matchesTriggerObject_[triggerPath->first] ) {
	// If there is only one trigger, but it is not this tau, then 
	// this is a probe object
	isProbeTau = true;
      }

      std::string triggerLabel;
      if ( triggerPath->first.rfind("_") == std::string::npos ) triggerLabel = triggerPath->first;
      else triggerLabel = std::string(triggerPath->first, triggerPath->first.rfind("_") + 1);
      
      std::string probeLabel = std::string("probe").append(triggerLabel);
      //std::cout << "probeLabel = " << probeLabel << ": " << isProbeTau << std::endl;
      newTau.addUserFloat(probeLabel.data(), isProbeTau);
      
      // set tag flag
      std::string tagLabel = std::string("tag").append(triggerLabel);
      //std::cout << "tagLabel = " << tagLabel << ": " << tauInfos[iTau].matchesTriggerObject_[triggerPath->first] << std::endl;
      newTau.addUserFloat(tagLabel.data(), tauInfos[iTau].matchesTriggerObject_[triggerPath->first]);      
    }

    // copy to output collection
    output->push_back(newTau);
  }
  
  // store in event
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauIdTagAndProbeProducer);
