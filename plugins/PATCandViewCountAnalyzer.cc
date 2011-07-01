#include "TauAnalysis/TauIdEfficiency/plugins/PATCandViewCountAnalyzer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

PATCandViewCountAnalyzer::PATCandViewCountAnalyzer(const edm::ParameterSet& cfg)
{
  //std::cout << "<PATCandCountAnalyzer::PATCandCountAnalyzer>:" << std::endl;

  int minNumEntries = ( cfg.exists("minNumEntries") ) ? cfg.getParameter<int>("minNumEntries") : -1;
  int maxNumEntries = ( cfg.exists("maxNumEntries") ) ? cfg.getParameter<int>("maxNumEntries") : -1;

  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag collections = cfg.getParameter<vInputTag>("src");
  for ( vInputTag::const_iterator collection = collections.begin();
	collection != collections.end(); ++collection ) {
    candCollections_.push_back(candCollectionEntryType(*collection, minNumEntries, maxNumEntries));
  }
}

PATCandViewCountAnalyzer::~PATCandViewCountAnalyzer()
{
// nothing to be done yet...
}

void PATCandViewCountAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{  
  //std::cout << "<PATCandViewCountAnalyzer::analyze>:" << std::endl; 

  for ( std::vector<candCollectionEntryType>::iterator candCollection = candCollections_.begin();
	candCollection != candCollections_.end(); ++candCollection ) {
    candCollection->analyze(evt);
  }
}

void PATCandViewCountAnalyzer::endJob()
{  
  std::cout << "<PATCandViewCountAnalyzer::endJob>:" << std::endl; 

  for ( std::vector<candCollectionEntryType>::iterator candCollection = candCollections_.begin();
	candCollection != candCollections_.end(); ++candCollection ) {
    candCollection->endJob();
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(PATCandViewCountAnalyzer);

