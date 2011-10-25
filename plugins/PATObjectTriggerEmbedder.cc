#include "TauAnalysis/TauIdEfficiency/plugins/PATObjectTriggerEmbedder.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

template <typename T>
PATObjectTriggerEmbedder<T>::PATObjectTriggerEmbedder(const edm::ParameterSet& cfg)
  : src_(cfg.getParameter<edm::InputTag>("src")),
    matched_(cfg.getParameter<edm::InputTag>("matched")),
    dRmax_(cfg.getParameter<double>("dRmax"))
{
  typedef std::vector<std::string> vstring;
  vstring triggerPaths = cfg.getParameter<vstring>("triggerPaths");
  for ( vstring::const_iterator triggerPath = triggerPaths.begin();
	triggerPath != triggerPaths.end(); ++triggerPath ) {
    std::string label = (*triggerPath);
    if ( label.find("_v") != std::string::npos ) label = std::string(label, 0, label.find("_v"));
    triggerPathsAndLabels_[*triggerPath] = label;
  }
  
typedef std::vector<std::string> vstring;
  std::map<std::string, std::string> triggerPathsAndLabels_;

  produces<PATObjectCollection>();
}

template <typename T>
PATObjectTriggerEmbedder<T>::~PATObjectTriggerEmbedder()
{
// nothing to be done yet...
}

template <typename T>
void PATObjectTriggerEmbedder<T>::produce(edm::Event& evt, const edm::EventSetup& es)
{
  std::auto_ptr<PATObjectCollection> patObjects_output(new PATObjectCollection());
  
  edm::Handle<PATObjectCollection> patObjects_input;
  evt.getByLabel(src_, patObjects_input);

  edm::Handle<pat::TriggerObjectStandAloneCollection > patTriggerObjects;
  evt.getByLabel(matched_, patTriggerObjects);

  for ( typename PATObjectCollection::const_iterator patObject_input = patObjects_input->begin();
	patObject_input != patObjects_input->end(); ++ patObject_input ) {
    T patObject_output(*patObject_input);

    std::map<std::string, float> patObjectUserFloats;
    for ( std::map<std::string, std::string>::const_iterator triggerPathAndLabel = triggerPathsAndLabels_.begin();
	  triggerPathAndLabel != triggerPathsAndLabels_.end(); ++triggerPathAndLabel ) {
      patObjectUserFloats[triggerPathAndLabel->second] = 0.;
    }

    for ( pat::TriggerObjectStandAloneCollection::const_iterator patTriggerObject = patTriggerObjects->begin();
	  patTriggerObject != patTriggerObjects->end(); ++ patTriggerObject ) {
      if ( deltaR(patTriggerObject->p4(), patObject_input->p4()) < dRmax_ ) {
	for ( std::map<std::string, std::string>::const_iterator triggerPathAndLabel = triggerPathsAndLabels_.begin();
	      triggerPathAndLabel != triggerPathsAndLabels_.end(); ++triggerPathAndLabel ) {
          if ( patTriggerObject->hasPathName(triggerPathAndLabel->first, true, false) ) {
	    patObjectUserFloats[triggerPathAndLabel->second] = 1.;
	  }
	}
      }
    }
    
    for ( std::map<std::string, float>::const_iterator patObjectUserFloat = patObjectUserFloats.begin();
	  patObjectUserFloat != patObjectUserFloats.end(); ++patObjectUserFloat ) {
      //std::cout << "setting userFloat(" << patObjectUserFloat->first << ") = " << patObjectUserFloat->second << std::endl;
      patObject_output.addUserFloat(patObjectUserFloat->first, patObjectUserFloat->second);
    }

    patObjects_output->push_back(patObject_output);
  }
  
  evt.put(patObjects_output);
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

typedef PATObjectTriggerEmbedder<pat::Electron> PATElectronTriggerEmbedder;
typedef PATObjectTriggerEmbedder<pat::Muon> PATMuonTriggerEmbedder;
typedef PATObjectTriggerEmbedder<pat::Tau> PATTauTriggerEmbedder;
typedef PATObjectTriggerEmbedder<pat::Jet> PATJetTriggerEmbedder;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATElectronTriggerEmbedder);
DEFINE_FWK_MODULE(PATMuonTriggerEmbedder);
DEFINE_FWK_MODULE(PATTauTriggerEmbedder);
DEFINE_FWK_MODULE(PATJetTriggerEmbedder);
