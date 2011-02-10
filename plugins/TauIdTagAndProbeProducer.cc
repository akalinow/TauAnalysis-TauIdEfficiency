#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


/*
 * TauIdTagAndProbeProducer
 *
 * Author: Evan K. Friis (UC Davis), Christian Veelken (UC Davis)
 *
 */

typedef std::vector<std::string> vstring;

class TauIdTagAndProbeProducer : public edm::EDProducer 
{  
 public:
  struct TauInfo 
  {
    TauInfo()
      : tau_(0), tauPt_(0.)
    {}
    TauInfo(const pat::Tau& tau)
      : tau_(&tau)
    {
      // get Pt of jet associated to tau
      if ( tau.isPFTau() ) {
	tauPt_ = tau.pfTauTagInfoRef()->pfjetRef()->pt();
      } else if ( tau.isCaloTau() ) {
	tauPt_ = tau.caloTauTagInfoRef()->jetRef()->pt();
      } else {
	edm::LogError("Invalid Tau Type")  << "pat::Tau is neither PF nor Calo !!";
      }
    }
    const pat::Tau* tau_;
    double tauPt_;
    std::map<std::string, bool> matchesTriggerObject_;
  };
  
  typedef std::vector<pat::Tau> vPatTaus;
  explicit TauIdTagAndProbeProducer(const edm::ParameterSet& pset);
  virtual ~TauIdTagAndProbeProducer(){}
  void produce(edm::Event&, const edm::EventSetup&);
  
private:
  edm::InputTag src_;
  vstring triggerPaths_;
};

TauIdTagAndProbeProducer::TauIdTagAndProbeProducer(const edm::ParameterSet& pset)
{
  src_ = pset.getParameter<edm::InputTag>("source");
  triggerPaths_ = pset.getParameter<vstring>("triggerPaths");
  
  // register products
  produces<vPatTaus>();
}

namespace {
  // Sorting predicate
  bool tauInfoDescendingPt(const TauIdTagAndProbeProducer::TauInfo& a, 
			   const TauIdTagAndProbeProducer::TauInfo& b)
  {
    return ( a.tauPt_ > b.tauPt_ );
  }
}

void
TauIdTagAndProbeProducer::produce(edm::Event& evt, const edm::EventSetup& es)
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

    for ( vstring::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {
      bool matchesTriggerObject = tau.triggerObjectMatchesByPath(*triggerPath).size();
      myTauInfo.matchesTriggerObject_[*triggerPath] = matchesTriggerObject;
      if ( matchesTriggerObject ) ++nTriggers[*triggerPath];
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
    
    for ( vstring::const_iterator triggerPath = triggerPaths_.begin();
	  triggerPath != triggerPaths_.end(); ++triggerPath ) {

      std::cout << "triggerPath = " << (*triggerPath) << std::endl;

      // set probe flag
      bool isProbeTau = false;
      if ( nTriggers[*triggerPath] > 1 ){
	// At least two triggers, so this is always a probe
	isProbeTau = true;
      } else if ( nTriggers[*triggerPath] == 1 && !tauInfos[iTau].matchesTriggerObject_[*triggerPath] ) {
	// If there is only one trigger, but it is not this tau, then 
	// this is a probe object
	isProbeTau = true;
      }

      std::string triggerLabel;
      if   ( triggerPath->rfind("_") == std::string::npos ) triggerLabel = (*triggerPath);
      else triggerLabel = std::string(*triggerPath, triggerPath->rfind("_") + 1);
      
      std::string probeLabel = std::string("probe").append(triggerLabel);
      std::cout << "probeLabel = " << probeLabel << std::endl;
      newTau.addUserFloat(probeLabel.data(), isProbeTau);
      
      // set tag flag
      std::string tagLabel = std::string("tag").append(triggerLabel);
      std::cout << "tagLabel = " << tagLabel << std::endl;
      newTau.addUserFloat(tagLabel.data(), tauInfos[iTau].matchesTriggerObject_[*triggerPath]);      
    }

    // copy to output collection
    output->push_back(newTau);
  }
  
  // store in event
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauIdTagAndProbeProducer);
