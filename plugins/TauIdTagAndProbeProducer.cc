#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"


/*
 * TauIdTagAndProbeProducer
 *
 * Author: Evan K. Friis (UC Davis), Christian Veelken (UC Davis)
 *
 */

using namespace edm;
using namespace std;

class TauIdTagAndProbeProducer : public EDProducer {

   struct TauInfo {
      const pat::Tau tau;
      bool matchesTriggerObject;
   };

   public:
      typedef std::vector<pat::Tau> vPatTaus;
      explicit TauIdTagAndProbeProducer(const ParameterSet& pset);
      virtual ~TauIdTagAndProbeProducer(){}
      void produce(Event&, const EventSetup&);

   private:
      InputTag src_;
      string triggerPath_;
};

TauIdTagAndProbeProducer::TauIdTagAndProbeProducer(const ParameterSet& pset)
{
   src_ = pset.getParameter<InputTag>("source");
   triggerPath_ = pset.getParameter<string>("triggerPath");

   // register products
   produces<vPatTaus>();
}

namespace {
   // Sorting predicate
   bool tauInfoDescendingPt(const TauIdTagAndProbeProducer::TauInfo &a, 
         const TauIdTagAndProbeProducer::TauInfo &b)
   {
      return (a.pt() > b.pt());
   }
}

void
TauIdTagAndProbeProducer::produce(Event &evt, const EventSetup &es)
{
   // output products
   auto_ptr<vPatTaus> output(new vPatTaus());

   Handle<View<pat::Tau> > sourceView;
   evt.getByLabel(src_, sourceView);

   size_t inputSize = sourceView.size();
   vector<TauInfo> tauInfos(inputSize);
   // count how many taus are matched to a trigger object
   unsigned int nTriggers = 0;
   for(size_t iTau = 0; iTau < inputSize; ++iTau)
   {
      // get tau and check if it matches trigger
      const pat::Tau* tau = sourceView[iTau].get();
      bool matchesTriggerObject = tau->triggerObjectMatchesByPath(triggerPath_);
      TauInfo myTauInfo;
      myTauInfo.tau = *tau;
      myTauInfo.matchesTriggerObject = matchesTriggerObject;
      tauInfos[iTau] = myTauInfo;
      if(matchesTriggerObject)
         nTriggers++;
   }

   // sort by descending pt
   sort(tauInfos.begin(), tauInfos.end(), tauInfoDescendingPt);

   for(size_t iTau = 0; iTau < inputSize; ++iTau)
   {
      // make our own copy of the tau
      pat::Tau newTau = tauInfos[iTau].tau;
      // store pt-ordered index
      newTau.addUserFloat("pt_index", iTau);

      // set probe flag
      bool isProbeTau = false;
      // At least two triggers, so this is always a probe
      if(nTriggers > 1)
      {
         isProbeTau = true;
      }
      // If there is only one trigger, but it is not this tau, then 
      // this is a probe object
      else if (nTriggers == 1 && !tauInfos[iTau].matchesTriggerObject)
      {
         isProbeTau = true;
      }
      newTau.addUserFloat("probe", isProbeTau);

      // set tag flag
      newTau.addUserFloat("tag", tauInfos[iTau].matchesTriggerObject);

      // copy to output collection
      output->push_back(newTau);
   }

   // store in event
   evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauIdTagAndProbeProducer);
