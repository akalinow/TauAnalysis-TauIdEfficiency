#include "TauAnalysis/TauIdEfficiency/plugins/TauIdTagAndProbeProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"




#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"



typedef std::vector<pat::Tau> vPatTaus;

TauIdTagAndProbeProducer::TauIdTagAndProbeProducer(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("source");

  edm::ParameterSet cfgTriggerPaths = cfg.getParameter<edm::ParameterSet>("triggerPaths");
  vstring triggerPathsNames = cfgTriggerPaths.getParameterNamesForType<vstring>();
  for ( vstring::const_iterator triggerPathsName = triggerPathsNames.begin();
	triggerPathsName != triggerPathsNames.end(); ++triggerPathsName ) {
    triggerPaths_[*triggerPathsName] = cfgTriggerPaths.getParameter<vstring>(*triggerPathsName);
  }
  
  // register products
  produces<vPatTaus>();
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

  edm::InputTag hltResultsSource("TriggerResults", "", "HLT");
  if ( hltResultsSource.label() != "" ) {
    edm::Handle<edm::TriggerResults> hltResults;
    evt.getByLabel(hltResultsSource, hltResults);
    if ( hltResults.isValid() ) {    
      const edm::TriggerNames& triggerNames = evt.triggerNames(*hltResults);
      for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
            triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
        unsigned int index = triggerNames.triggerIndex(*triggerName);
        if ( index < triggerNames.size() ) {
          std::string triggerDecision = ( hltResults->accept(index) ) ? "passed" : "failed";
     
          std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
        }
      }

      std::cout << "HLT Decisions:" << std::endl;
    
      vstring hltPathsToPrint;
      hltPathsToPrint.push_back("HLT_Jet30_v1");
      hltPathsToPrint.push_back("HLT_Jet30_v2");
      hltPathsToPrint.push_back("HLT_Jet30_v3");
      hltPathsToPrint.push_back("HLT_Jet30_v4");
      hltPathsToPrint.push_back("HLT_Jet30_v5");
      hltPathsToPrint.push_back("HLT_Jet30_v6");
      
      for ( std::vector<std::string>::const_iterator hltPath = hltPathsToPrint.begin();
	    hltPath != hltPathsToPrint.end(); ++hltPath ) {
        unsigned int index = triggerNames.triggerIndex(*hltPath);
        if ( index < triggerNames.size() ) {
  	  std::string hltDecision = ( hltResults->accept(index) ) ? "passed" : "failed";	
	  std::cout << " " << (*hltPath) << " " << hltDecision << std::endl;
        } else {
	  edm::LogError ("printEventTriggerInfo") << " Undefined trigger Path = " << (*hltPath) << " --> skipping !!";
	  continue;
        }
      }
    }
    
    std::cout << std::endl;
  }




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

    std::cout << "matched Trigger paths: " << std::endl;
    std::cout << tau.triggerObjectMatches().size() << std::endl;
    for ( pat::TriggerObjectStandAloneCollection::const_iterator matchedTrigger = tau.triggerObjectMatches().begin();
	  matchedTrigger != tau.triggerObjectMatches().end(); ++matchedTrigger ) {
      std::cout << "break-point 1 reached" << std::endl;
      std::cout << "pt = " << matchedTrigger->pt() << std::endl;
      vstring matchedTriggerPaths = matchedTrigger->pathNames();
      for ( vstring::const_iterator matchedTriggerPath = matchedTriggerPaths.begin();
	    matchedTriggerPath != matchedTriggerPaths.end(); ++matchedTriggerPath ) {
	std::cout << " " << (*matchedTriggerPath) << std::endl;
      }
    }

    TauInfo myTauInfo(tau);

    for ( std::map<std::string, vstring>::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {
      bool matchesTriggerObject = false;
      for ( vstring::const_iterator triggerPath_version = triggerPath->second.begin();
	    triggerPath_version != triggerPath->second.end(); ++triggerPath_version ) {
	if ( tau.triggerObjectMatchesByPath(*triggerPath_version).size() >= 1 ) matchesTriggerObject = true;
      }

      myTauInfo.matchesTriggerObject_[triggerPath->first] = matchesTriggerObject;
      if ( matchesTriggerObject ) ++nTriggers[triggerPath->first];
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
    
    for ( std::map<std::string, vstring>::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {

      std::cout << "triggerPath = " << triggerPath->first << std::endl;

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
      std::cout << "probeLabel = " << probeLabel << ": " << isProbeTau << std::endl;
      newTau.addUserFloat(probeLabel.data(), isProbeTau);
      
      // set tag flag
      std::string tagLabel = std::string("tag").append(triggerLabel);
      std::cout << "tagLabel = " << tagLabel << ": " << tauInfos[iTau].matchesTriggerObject_[triggerPath->first] << std::endl;
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
