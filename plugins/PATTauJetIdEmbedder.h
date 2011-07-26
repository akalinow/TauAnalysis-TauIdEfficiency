#ifndef TauAnalysis_TauIdEfficiency_PATTauJetIdEmbedder_h  
#define TauAnalysis_TauIdEfficiency_PATTauJetIdEmbedder_h

/** \class PATTauJetIdEmbedder
 *
 * Auxiliary class for extracting jetId bits
 * for reco::CaloJet/reco::PFJet object associated 
 * to reconstructed PAT tau objects.
 * The extracted jetId bits are added as userFloats 
 * to the pat::Tau object.
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: PATTauJetIdEmbedder.h,v 1.1 2011/02/06 10:41:00 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <string>
#include <vector>

class PATTauJetIdEmbedder : public edm::EDProducer
{
 public:
  
  explicit PATTauJetIdEmbedder(const edm::ParameterSet&);
  ~PATTauJetIdEmbedder();
    
  void produce(edm::Event&, const edm::EventSetup&);

 private:
  
//--- configuration parameters
  edm::InputTag src_;
  
  edm::InputTag srcJet_;
 
  struct jetIdType
  {
    jetIdType(const std::string& value)
      : value_(value),
	caloJetId_(0),
	pfJetId_(0)
    {
      if      ( value_ == "jetIdMinimal"   ) quality_ = "MINIMAL";
      else if ( value_ == "jetIdLoose"     ) quality_ = "LOOSE";
      else if ( value_ == "jetIdLoose_AOD" ) quality_ = "LOOSE_AOD";
      else if ( value_ == "jetIdMedium"    ) quality_ = "MEDIUM";
      else if ( value_ == "jetIdTight"     ) quality_ = "TIGHT";
      else throw cms::Exception("jetIdType") 
	<< "Inalid 'value' = " << value << " passed as function parameter !!\n";

      edm::ParameterSet cfgJetId;
      cfgJetId.addParameter<std::string>("version", "PURE09");
      cfgJetId.addParameter<std::string>("quality", quality_);
  
      caloJetId_ = new JetIDSelectionFunctor(cfgJetId);
      pfJetId_ = new PFJetIDSelectionFunctor(cfgJetId);
    }
    ~jetIdType()
    {
      delete caloJetId_;
      delete pfJetId_;
    }

    std::string value_;

    std::string quality_;

    JetIDSelectionFunctor* caloJetId_;
    PFJetIDSelectionFunctor* pfJetId_;
  };

  std::vector<jetIdType*> jetIds_;
};

#endif  

