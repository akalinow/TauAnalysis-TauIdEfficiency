#ifndef TauAnalysis_TauIdEfficiency_PATCandViewCountAnalyzer_h  
#define TauAnalysis_TauIdEfficiency_PATCandViewCountAnalyzer_h

 /** \class PATCandViewCountAnalyzer
  *
  * Print number of entries in collections of particles
  * (module to be used for debugging purposes)
  * 
  * \author Christian Veelken, UC Davis
  *
  * \version $Revision: 1.2 $
  *
  * $Id: PATCandCountAnalyzer.h,v 1.2 2010/09/28 11:23:39 jkolb Exp $
  *
  */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include <string>
#include <vector>

class PATCandViewCountAnalyzer : public edm::EDAnalyzer 
{
 public: 

  explicit PATCandViewCountAnalyzer(const edm::ParameterSet&);
  ~PATCandViewCountAnalyzer();
  
  void beginJob() {}
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();

 private:

  struct candCollectionEntryType
  {
    candCollectionEntryType(const edm::InputTag& src, int minNumEntries = -1, int maxNumEntries = -1)
      : src_(src),
        minNumEntries_(minNumEntries),
        maxNumEntries_(maxNumEntries),
	numEventsProcessed_(0),
        numEventsPassed_(0)
    {
      numEntries_.resize(12);
    }
    ~candCollectionEntryType() {}
    void analyze(const edm::Event& evt)
    {
      typedef edm::View<reco::Candidate> CandidateView;
      edm::Handle<CandidateView> candidates;
      evt.getByLabel(src_, candidates);
      int numCandidates = candidates->size();
      //std::cout << " " << src_.label() << ":" 
      //	  << " " << numCandidates << std::endl;
      if ( numCandidates >= 0 && numCandidates <= 10 ) ++numEntries_[numCandidates];
      else if ( numCandidates > 10 ) ++numEntries_[11];
      ++numEventsProcessed_;
      if ( numCandidates >= minNumEntries_ && numCandidates <= maxNumEntries_ ) ++numEventsPassed_;      
      std::vector<int> numEntries_;
    }
    void endJob() 
    {
      std::cout << " " << src_.label() << ":" 
		<< " processed = " << numEventsProcessed_ << ", passed = " << numEventsPassed_ << std::endl;
      std::cout << "numEntries = " << format_vint(numEntries_) << std::endl;
    }
    edm::InputTag src_;
    int minNumEntries_;
    int maxNumEntries_;
    std::vector<int> numEntries_;
    long numEventsProcessed_;
    long numEventsPassed_;
  };

  std::vector<candCollectionEntryType> candCollections_;
};

#endif  


