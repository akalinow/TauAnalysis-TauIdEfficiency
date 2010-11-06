
/** \class MCEmbeddingTagAndProbeProducer
 *
 * Flag tau-jets embedded into Zmumu events as tag and probe, 
 * depending on whether the muon replaced by the tau-jet 
 * was isolated (tag, probe flags both set) or not (only probe flag set)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: MCEmbeddingTagAndProbeProducer.cc,v 1.4 2010/03/29 17:08:41 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class MCEmbeddingTagAndProbeProducer : public edm::EDProducer 
{  
 public:
  explicit MCEmbeddingTagAndProbeProducer(const edm::ParameterSet& pset);
  virtual ~MCEmbeddingTagAndProbeProducer() {}

  void produce(edm::Event&, const edm::EventSetup&);
  
private:
  edm::InputTag src_;
  
  edm::InputTag srcTags_;
  edm::InputTag srcProbes_;

  double dRmatch_;
};

MCEmbeddingTagAndProbeProducer::MCEmbeddingTagAndProbeProducer(const edm::ParameterSet& cfg)
{
  src_ = cfg.getParameter<edm::InputTag>("source");

  srcTags_ = cfg.getParameter<edm::InputTag>("srcTags");
  srcProbes_ = cfg.getParameter<edm::InputTag>("srcTags");

  dRmatch_ = cfg.getParameter<double>("dRmatch");

  // register products
  produces<pat::TauCollection>();
}

bool matches(const pat::Tau& patTau, const std::vector<reco::Candidate::LorentzVector>& flags, double dRmatch)
{
  for ( std::vector<reco::Candidate::LorentzVector>::const_iterator flag = flags.begin();
	flag != flags.end(); ++flag ) {
    float dR = reco::deltaR(patTau.p4(), *flag);
    if ( dR < dRmatch ) return true;
  }
  return false;
}

void MCEmbeddingTagAndProbeProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  // output products
  std::auto_ptr<pat::TauCollection> output(new pat::TauCollection());
  
  edm::Handle<pat::TauCollection> patTaus;
  evt.getByLabel(src_, patTaus);

  typedef std::vector<reco::Candidate::LorentzVector> vLorentzVector;
  edm::Handle<vLorentzVector> tags;
  evt.getByLabel(srcTags_, tags);
  edm::Handle<vLorentzVector> probes;
  evt.getByLabel(srcProbes_, probes);

  for ( pat::TauCollection::const_iterator patTau = patTaus->begin();
	patTau != patTaus->end(); ++patTau ) {
    // create new tau object
    pat::Tau newTau = (*patTau);

    // check tag & probe flags
    bool isTagTau = matches(*patTau, *tags, dRmatch_);
    bool isProbeTau = matches(*patTau, *probes, dRmatch_);

    // add tag & probe flags to tau object
    newTau.addUserFloat("tag", isTagTau);
    newTau.addUserFloat("probe", isProbeTau);
    
    // copy to output collection
    output->push_back(newTau);
  }
  
  // store in event
  evt.put(output);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MCEmbeddingTagAndProbeProducer);
