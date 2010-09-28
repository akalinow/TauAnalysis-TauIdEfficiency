#ifndef TauAnalysis_TauIdEfficiency_PATTauCollectionDiffAnalyzer_h  
#define TauAnalysis_TauIdEfficiency_PATTauCollectionDiffAnalyzer_h

 /** \class PATTauCollectionDiffAnalyzer
  *
  * Compare content of two (PAT)Tau collections
  * (moduleto be used for debugging purposes)
  * 
  * \author Christian Veelken, UC Davis
  *
  * \version $Revision: 1.1 $
  *
  * $Id: PATTauCollectionDiffAnalyzer.h,v 1.1 2010/06/20 16:00:56 veelken Exp $
  *
  */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include <string>

class PATTauCollectionDiffAnalyzer : public edm::EDAnalyzer 
{
 public: 
  explicit PATTauCollectionDiffAnalyzer(const edm::ParameterSet&);
  ~PATTauCollectionDiffAnalyzer();
  
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob() {}

 private:
  edm::InputTag patTauSource1_;
  edm::InputTag patTauSource2_;

  std::string patTauSelection_;
  StringCutObjectSelector<pat::Tau>* patTauSelector_;

  double dRmatch_;

  std::string dqmDirectory_;

  int dqmError_;

  struct histogramCollectionType
  {
    MonitorElement* jetPt_;
    MonitorElement* jetEta_;
  };

  histogramCollectionType histogramsPatTaus1_;
  histogramCollectionType histogramsPatTaus2_;

  histogramCollectionType histogramsPatTausSel1_;
  histogramCollectionType histogramsPatTausSel2_;
};

#endif  


