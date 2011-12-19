
/** \executable FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.26 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.26 2011/11/14 13:57:00 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffHistManager.h"
#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"
#include "TauAnalysis/RecoTools/interface/PATObjectLUTvalueExtractorFromKNN.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TF1.h>

#include <vector>
#include <string>
#include <fstream>

typedef std::vector<std::string> vstring;

struct histManagerEntryType
{
  histManagerEntryType(const edm::ParameterSet& cfg, bool fillGenMatchHistograms, 
		       const std::string& binVariable = "", double min = -0.5, double max = +0.5)
    : binVariable_(binVariable),
      min_(min),
      max_(max),
      fillGenMatchHistograms_(fillGenMatchHistograms)
  {
    histManager_ = new TauIdEffHistManager(cfg);

    if ( fillGenMatchHistograms_ ) {
      const std::string& label = cfg.getParameter<std::string>("label");

      edm::ParameterSet cfgJetToTauFake(cfg);
      std::string labelJetToTauFake = std::string(label).append("_").append("JetToTauFake");
      cfgJetToTauFake.addParameter<std::string>("label", labelJetToTauFake);
      histManagerJetToTauFake_ = new TauIdEffHistManager(cfgJetToTauFake);

      edm::ParameterSet cfgMuToTauFake(cfg);
      std::string labelMuToTauFake = std::string(label).append("_").append("MuToTauFake");
      cfgMuToTauFake.addParameter<std::string>("label", labelMuToTauFake);
      histManagerMuToTauFake_ = new TauIdEffHistManager(cfgMuToTauFake);

      edm::ParameterSet cfgGenTau(cfg);
      std::string labelGenTau = std::string(label).append("_").append("GenTau");
      cfgGenTau.addParameter<std::string>("label", labelGenTau);
      histManagerGenTau_ = new TauIdEffHistManager(cfgGenTau);
    }
  }
  ~histManagerEntryType() {}
  void bookHistograms(TFileDirectory& dir)
  {
    histManager_->bookHistograms(dir);

    if ( fillGenMatchHistograms_ ) {
      histManagerJetToTauFake_->bookHistograms(dir);
      histManagerMuToTauFake_->bookHistograms(dir);
      histManagerGenTau_->bookHistograms(dir);
    }
  }
  void fillHistograms(double x, const PATMuTauPair& muTauPair, const pat::MET& caloMEt, 
		      size_t numVertices, const std::map<std::string, bool>& plot_triggerBits_passed, int genMatchType, double weight)
  {
    if ( x > min_ && x <= max_ ) {
      histManager_->fillHistograms(muTauPair, caloMEt, numVertices, plot_triggerBits_passed, weight);

      if ( fillGenMatchHistograms_ ) {
	if      ( genMatchType == kJetToTauFakeMatched ) 
	  histManagerJetToTauFake_->fillHistograms(muTauPair, caloMEt, numVertices, plot_triggerBits_passed, weight);
	else if ( genMatchType == kMuToTauFakeMatched  ) 
	  histManagerMuToTauFake_->fillHistograms(muTauPair, caloMEt, numVertices, plot_triggerBits_passed, weight);
	else if ( genMatchType == kGenTauHadMatched    ||
		  genMatchType == kGenTauOtherMatched  ) 
	  histManagerGenTau_->fillHistograms(muTauPair, caloMEt, numVertices, plot_triggerBits_passed, weight);
      }
    }
  }

  std::string binVariable_;

  double min_;
  double max_;

  TauIdEffHistManager* histManager_;

  bool fillGenMatchHistograms_;

  TauIdEffHistManager* histManagerJetToTauFake_;
  TauIdEffHistManager* histManagerMuToTauFake_;
  TauIdEffHistManager* histManagerGenTau_;
};

struct regionEntryType
{
  regionEntryType(fwlite::TFileService& fs,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName, const std::string& sysShift,
		  const edm::ParameterSet& cfgBinning, const std::string& svFitMassHypothesis, 
		  const std::string& tauChargeMode, bool disableTauCandPreselCuts, const edm::ParameterSet& cfgEventSelCuts, 
		  bool fillGenMatchHistograms, const vstring& plot_triggerBits,
		  const std::string& selEventsFileName)
    : process_(process),
      region_(region),
      tauIdDiscriminators_(tauIdDiscriminators),
      tauIdName_(tauIdName),
      sysShift_(sysShift),
      selector_(0),
      histogramsUnbinned_(0),
      numMuTauPairs_selected_(0),
      numMuTauPairsWeighted_selected_(0.),
      selEventsFile_(0)
  {
    edm::ParameterSet cfgSelector = cfgEventSelCuts;
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators_);
    cfgSelector.addParameter<std::string>("region", region_);
    cfgSelector.addParameter<std::string>("tauChargeMode", tauChargeMode);
    cfgSelector.addParameter<bool>("disableTauCandPreselCuts", disableTauCandPreselCuts);

    selector_ = new TauIdEffEventSelector(cfgSelector);

    edm::ParameterSet cfgHistManager;
    cfgHistManager.addParameter<std::string>("process", process_);
    cfgHistManager.addParameter<std::string>("region", region_);
    cfgHistManager.addParameter<std::string>("tauIdDiscriminator", tauIdName_);
    if      ( region.find("p") != std::string::npos ) label_ = "passed";
    else if ( region.find("f") != std::string::npos ) label_ = "failed";
    else                                              label_ = "all";
    if ( sysShift_ != "CENTRAL_VALUE" ) label_.append("_").append(sysShift_);
    cfgHistManager.addParameter<std::string>("label", label_);
    cfgHistManager.addParameter<std::string>("svFitMassHypothesis", svFitMassHypothesis);
    cfgHistManager.addParameter<vstring>("triggerBits", plot_triggerBits);

    histogramsUnbinned_ = new histManagerEntryType(cfgHistManager, fillGenMatchHistograms);
    histogramsUnbinned_->bookHistograms(fs);

    typedef std::vector<edm::ParameterSet> vParameterSet;
    vstring binVariableNames = cfgBinning.getParameterNamesForType<vParameterSet>();
    for ( vstring::const_iterator binVariableName = binVariableNames.begin();
	  binVariableName != binVariableNames.end(); ++binVariableName ) {
      vParameterSet cfgBinVariableBins = cfgBinning.getParameter<vParameterSet>(*binVariableName);
      for ( vParameterSet::const_iterator cfgBinVariableBin = cfgBinVariableBins.begin();
	    cfgBinVariableBin != cfgBinVariableBins.end(); ++cfgBinVariableBin ) {
	double min = cfgBinVariableBin->getParameter<double>("min");
	double max = cfgBinVariableBin->getParameter<double>("max");
	histManagerEntryType* histManagerEntry = 
	  new histManagerEntryType(cfgHistManager, fillGenMatchHistograms, 
				   *binVariableName, min, max);
	std::string dir_string = cfgBinVariableBin->getParameter<std::string>("subdir");
	TFileDirectory dir = fs.mkdir(dir_string);
	histManagerEntry->bookHistograms(dir);
	histogramEntriesBinned_.push_back(histManagerEntry);
      }
    }

    if ( selEventsFileName != "" ) {
      size_t idx = selEventsFileName.rfind(".");
      if ( idx != std::string::npos ) {
	std::string selEventsFileName_region = std::string(selEventsFileName, 0, idx);
	selEventsFileName_region.append("_").append(region_);
	selEventsFileName_region.append(std::string(selEventsFileName, idx));
	//std::cout << "selEventsFileName_region = " << selEventsFileName_region << std::endl;
	selEventsFile_ = new std::ofstream(selEventsFileName_region.data(), std::ios::out);
      } else throw cms::Exception("regionEntryType")
	  << "Invalid selEventsFileName = " << selEventsFileName << " !!\n";
    }
  }
  ~regionEntryType()
  {
    delete selector_;

    delete histogramsUnbinned_;

    for ( std::vector<histManagerEntryType*>::iterator it = histogramEntriesBinned_.begin();
	  it != histogramEntriesBinned_.end(); ++it ) {
      delete (*it);
    }
    
    delete selEventsFile_;
  }
  void analyze(const fwlite::Event& evt, const PATMuTauPair& muTauPair, const pat::MET& caloMEt, 
	       size_t numVertices, const std::map<std::string, bool>& plot_triggerBits_passed, int genMatchType, double evtWeight)
  {
    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, caloMEt, evtSelFlags) ) {
//--- fill histograms for "inclusive" tau id. efficiency measurement
      histogramsUnbinned_->fillHistograms(0., muTauPair, caloMEt, numVertices, plot_triggerBits_passed, genMatchType, evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of 
//   o tau-jet transverse momentum
//   o tau-jet pseudo-rapidity
//   o reconstructed vertex multiplicity
//   o sumEt
//   o ...
      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = histogramEntriesBinned_.begin();
	    histManagerEntry != histogramEntriesBinned_.end(); ++histManagerEntry ) {
	double x = 0.;
	if      ( (*histManagerEntry)->binVariable_ == "tauPt"       ) x = muTauPair.leg2()->pt();
	else if ( (*histManagerEntry)->binVariable_ == "tauAbsEta"   ) x = TMath::Abs(muTauPair.leg2()->eta());
	else if ( (*histManagerEntry)->binVariable_ == "numVertices" ) x = numVertices;
	else if ( (*histManagerEntry)->binVariable_ == "sumEt"       ) x = muTauPair.met()->sumEt();
	else throw cms::Exception("regionEntryType::analyze")
	  << "Invalid binVariable = " << (*histManagerEntry)->binVariable_ << " !!\n";
	(*histManagerEntry)->fillHistograms(x, muTauPair, caloMEt, numVertices, plot_triggerBits_passed, genMatchType, evtWeight);
      }
 
      if ( selEventsFile_ ) 
	(*selEventsFile_) << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() << std::endl;

      ++numMuTauPairs_selected_;
      numMuTauPairsWeighted_selected_ += evtWeight;
    }
  }

  std::string process_;
  std::string region_;
  vstring tauIdDiscriminators_;
  std::string tauIdName_;
  std::string sysShift_;
  std::string label_;

  bool applyMuonIsoWeights;
  bool appyTauFakeRateWeights;

  TauIdEffEventSelector* selector_;

  histManagerEntryType* histogramsUnbinned_;
  std::vector<histManagerEntryType*> histogramEntriesBinned_;

  int numMuTauPairs_selected_;
  double numMuTauPairsWeighted_selected_;

  std::ofstream* selEventsFile_;
};

void checkL1Bits(const fwlite::Event& evt,
		 const vstring& l1Bits, 
		 const edm::InputTag& srcL1GtReadoutRecord, const edm::InputTag& srcL1GtObjectMapRecord,
		 std::map<std::string, bool>* l1Bits_passed, bool* anyL1bit_passed)
{
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  evt.getByLabel(srcL1GtReadoutRecord, l1GtReadoutRecord);
  edm::Handle<L1GlobalTriggerObjectMapRecord> l1GtObjectMapRecord;
  evt.getByLabel(srcL1GtObjectMapRecord, l1GtObjectMapRecord);

  DecisionWord l1GtDecision = l1GtReadoutRecord->decisionWord();
  const std::vector<L1GlobalTriggerObjectMap>& l1GtObjectMaps = l1GtObjectMapRecord->gtObjectMap();

  if ( anyL1bit_passed ) {
    (*anyL1bit_passed) = false;
  }

  for ( vstring::const_iterator l1Bit = l1Bits.begin();
	l1Bit != l1Bits.end(); ++l1Bit ) {	  
    bool isL1bit_passed = false;
    for ( std::vector<L1GlobalTriggerObjectMap>::const_iterator l1GtObjectMap = l1GtObjectMaps.begin();
	  l1GtObjectMap != l1GtObjectMaps.end(); ++l1GtObjectMap ) {
      std::string l1Bit_idx = (*l1GtObjectMap).algoName();
      int idx = (*l1GtObjectMap).algoBitNumber();	    
      if ( l1Bit_idx == (*l1Bit) ) isL1bit_passed = l1GtDecision[idx];
    }

    if ( l1Bits_passed ) {
      (*l1Bits_passed)[*l1Bit] = isL1bit_passed;
    }

    if ( anyL1bit_passed ) {
      if ( isL1bit_passed ) (*anyL1bit_passed) = true;
    }
  }
}

std::string getHLTpath_key(const std::string& hltPath)
{
  std::string key = hltPath;
  size_t idx = hltPath.find_last_of("_v");
  if ( idx != std::string::npos ) {
    // CV: std::string.find_last_of returns position of 'v' character
    //    --> need to decrease 'idx' by one in order to eliminate trailing underscore
    key = std::string(hltPath, 0, idx - 1);
  }
  return key;
}

void checkHLTpaths(const fwlite::Event& evt, 
		   const vstring& hltPaths,
		   const edm::InputTag& srcHLTresults,
		   std::map<std::string, bool>* hltPaths_passed, bool* anyHLTpath_passed)
{
  edm::Handle<edm::TriggerResults> hltResults;
  evt.getByLabel(srcHLTresults, hltResults);

  const edm::TriggerNames& triggerNames = evt.triggerNames(*hltResults);

  if ( hltPaths_passed ) {
    for ( vstring::const_iterator hltPath = hltPaths.begin();
	  hltPath != hltPaths.end(); ++hltPath ) {
      std::string key = getHLTpath_key(*hltPath);
      (*hltPaths_passed)[key] = false;
    }
  }
  
  if ( anyHLTpath_passed ) {
    (*anyHLTpath_passed) = false;
  }

  for ( vstring::const_iterator hltPath = hltPaths.begin();
	hltPath != hltPaths.end(); ++hltPath ) {
    bool isHLTpath_passed = false;
    unsigned int idx = triggerNames.triggerIndex(*hltPath);
    if ( idx < triggerNames.size() ) {
      isHLTpath_passed = hltResults->accept(idx);
    }

    if ( isHLTpath_passed ) {
      if ( hltPaths_passed ) {
	std::string key = getHLTpath_key(*hltPath);
	(*hltPaths_passed)[key] = isHLTpath_passed;
      }
    
      if ( anyHLTpath_passed ) {
	(*anyHLTpath_passed) = true;
      }
    }
  }
}

//-------------------------------------------------------------------------------
//
// Integral of Crystal Ball function for fitting trigger efficiency turn-on curves
// (code from Pascal Paganini)
//
double integralCrystalBall(double m, double m0, double sigma, double alpha, double n, double norm) 
{
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  
  double sig = fabs((double)sigma);
  
  double t = (m - m0)/sig;
  
  if (alpha < 0) t = -t;
  
  double absAlpha = fabs(alpha / sig);
  double a = TMath::Power(n/absAlpha, n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  
  if ( a >= std::numeric_limits<double>::max() ) return -1.;
  
  double approxErf;
  double arg = absAlpha / sqrt2;
  if      ( arg >  5. ) approxErf =  1;
  else if ( arg < -5. ) approxErf = -1;
  else                  approxErf = erf(arg);
  
  double leftArea = (1 + approxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  
  if ( t <= absAlpha ) {
    arg = t / sqrt2;
    if      ( arg >  5.) approxErf =  1;
    else if ( arg < -5.) approxErf = -1;
    else                 approxErf = erf(arg);
    return norm * (1 + approxErf) * sqrtPiOver2 / area;
  } else {
    return norm * (leftArea +  a * (1/TMath::Power(t-b, n - 1) - 1/TMath::Power(absAlpha - b, n - 1)) / (1 - n)) / area;
  }
}

Double_t integralCrystalBall_data_div_mc(Double_t* x, Double_t* par)
{
  if  ( x[0] < 100. ) {
    double data = integralCrystalBall(x[0], par[0], par[1], par[2], par[3], par[4]);
    double mc   = integralCrystalBall(x[0], par[5], par[6], par[7], par[8], par[9]);
    return ( mc > 0. ) ? (data/mc) : 1.;
  } else {
    return 1.;
  }
}
//-------------------------------------------------------------------------------

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauIdEffAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteTauIdEffAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauIdEffAnalyzer = cfg.getParameter<edm::ParameterSet>("tauIdEffAnalyzer");

  edm::InputTag srcMuTauPairs = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcMuTauPairs");
  bool requireUniqueMuTauPair = ( cfgTauIdEffAnalyzer.exists("") ) ? 
    cfgTauIdEffAnalyzer.getParameter<bool>("requireUniqueMuTauPair") : false;
  std::string svFitMassHypothesis = cfgTauIdEffAnalyzer.getParameter<std::string>("svFitMassHypothesis");
  std::string tauChargeMode = cfgTauIdEffAnalyzer.getParameter<std::string>("tauChargeMode");
  bool disableTauCandPreselCuts = cfgTauIdEffAnalyzer.getParameter<bool>("disableTauCandPreselCuts");
  edm::ParameterSet cfgEventSelCuts = cfgTauIdEffAnalyzer.getParameter<edm::ParameterSet>("eventSelCuts");
  edm::InputTag srcHLTresults = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcHLTresults");
  edm::InputTag srcL1GtReadoutRecord = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcL1GtReadoutRecord"); 
  edm::InputTag srcL1GtObjectMapRecord = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcL1GtObjectMapRecord");
  vstring hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("hltPaths");
  vstring l1Bits = cfgTauIdEffAnalyzer.getParameter<vstring>("l1Bits");
  edm::InputTag srcCaloMEt = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcCaloMEt");
  edm::InputTag srcGoodMuons = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGoodMuons");
  edm::InputTag srcVertices = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcVertices");
  edm::InputTag srcGenParticles = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGenParticles");
  bool fillGenMatchHistograms = cfgTauIdEffAnalyzer.getParameter<bool>("fillGenMatchHistograms");
  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch = cfgTauIdEffAnalyzer.getParameter<vint>("skipPdgIdsGenParticleMatch");  
  vstring plot_l1Bits = cfgTauIdEffAnalyzer.getParameter<vstring>("plot_l1Bits");
  vstring plot_hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("plot_hltPaths");
  vstring plot_triggerBits;
  plot_triggerBits.insert(plot_triggerBits.end(), plot_l1Bits.begin(), plot_l1Bits.end());
  for ( vstring::const_iterator plot_hltPath = plot_hltPaths.begin();
	plot_hltPath != plot_hltPaths.end(); ++plot_hltPath ) {
    plot_triggerBits.push_back(getHLTpath_key(*plot_hltPath));
  }
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffAnalyzer.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffAnalyzer.exists("sysShift") ?
    cfgTauIdEffAnalyzer.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";
  edm::InputTag srcEventCounter = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcEventCounter");

  PATMuonLUTvalueExtractorFromKNN* muonIsoProbExtractor = 0;
  if ( cfgTauIdEffAnalyzer.exists("muonIsoProbExtractor") ) {
    edm::ParameterSet cfgMuonIsoProbExtractor = cfgTauIdEffAnalyzer.getParameter<edm::ParameterSet>("muonIsoProbExtractor");
    muonIsoProbExtractor = new PATMuonLUTvalueExtractorFromKNN(cfgMuonIsoProbExtractor);
  }

  TF1* triggerEffCorrection = new TF1("triggerEffCorrection", &integralCrystalBall_data_div_mc, 0., 1.e+6, 10);
  double triggerEffCorr_parameter[] = {
    3.08698e+01, 4.12547e+00, 7.67225e+00, 1.55325e+02, 1.00826e+00,
    3.07315e+01, -4.16879e+00, 1.37289e+01, 1.02238e+00, 1.03602e+00
  };
  for ( int iPar = 0; iPar < triggerEffCorrection->GetNpar(); ++iPar ) {
    triggerEffCorrection->SetParameter(iPar, triggerEffCorr_parameter[iPar]);
  }

  std::string selEventsFileName = ( cfgTauIdEffAnalyzer.exists("selEventsFileName") ) ? 
    cfgTauIdEffAnalyzer.getParameter<std::string>("selEventsFileName") : "";

  fwlite::InputSource inputFiles(cfg); 
  edm::ParameterSet cfgInputSource = cfg.getParameter<edm::ParameterSet>("fwliteInput");
  int firstRun = cfgInputSource.getParameter<int>("firstRun");
  int lastRun = cfgInputSource.getParameter<int>("lastRun");
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauIdEffAnalyzer.getParameter<std::string>("process");
  std::cout << " process = " << process << std::endl;
  std::string processType = cfgTauIdEffAnalyzer.getParameter<std::string>("type");
  std::cout << " type = " << processType << std::endl;
  bool isData = (processType == "Data");
  vstring regions = cfgTauIdEffAnalyzer.getParameter<vstring>("regions");
  edm::ParameterSet cfgBinning = cfgTauIdEffAnalyzer.getParameter<edm::ParameterSet>("binning");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauIdEffAnalyzer.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");

      // all tau charges
      regionEntryType* regionEntry = 
	new regionEntryType(fs, process, *region, tauIdDiscriminators, tauIdName, 
			    sysShift, cfgBinning, svFitMassHypothesis, 
			    tauChargeMode, disableTauCandPreselCuts, cfgEventSelCuts, 
			    fillGenMatchHistograms, plot_triggerBits, selEventsFileName);
      regionEntries.push_back(regionEntry);

      // tau+ candidates only
      //regionEntryType* regionEntry_plus = 
      //  new regionEntryType(fs, process, std::string(*region).append("+"), tauIdDiscriminators, tauIdName, 
      //		      sysShift, cfgBinning, svFitMassHypothesis, 
      //		      tauChargeMode, disableTauCandPreselCuts, cfgEventSelCuts, 
      //		      fillGenMatchHistograms, plot_triggerBits, selEventsFileName);
      //regionEntries.push_back(regionEntry_plus);
      //
      // tau- candidates only
      //regionEntryType* regionEntry_minus = 
      //  new regionEntryType(fs, process, std::string(*region).append("-"), tauIdDiscriminators, tauIdName, 
      //		      sysShift, cfgBinning, svFitMassHypothesis, 
      //		      tauChargeMode, disableTauCandPreselCuts, cfgEventSelCuts, 
      //		      fillGenMatchHistograms, plot_triggerBits, selEventsFileName);
      //regionEntries.push_back(regionEntry_minus);
    }
  }

  edm::ParameterSet cfgSelectorABCD = cfgEventSelCuts;
  cfgSelectorABCD.addParameter<vstring>("tauIdDiscriminators", vstring());
  cfgSelectorABCD.addParameter<std::string>("region", "ABCD");
  cfgSelectorABCD.addParameter<std::string>("tauChargeMode", tauChargeMode);
  cfgSelectorABCD.addParameter<bool>("disableTauCandPreselCuts", disableTauCandPreselCuts);

  TauIdEffEventSelector* selectorABCD = new TauIdEffEventSelector(cfgSelectorABCD);

//--- book "dummy" histogram counting number of processed events
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 3, -0.5, +2.5);
  histogramEventCounter->GetXaxis()->SetBinLabel(1, "all Events (DBS)");      // CV: bin numbers start at 1 (not 0) !!
  histogramEventCounter->GetXaxis()->SetBinLabel(2, "processed by Skimming");
  histogramEventCounter->GetXaxis()->SetBinLabel(3, "analyzed in PAT-tuple");
  
  int allEvents_DBS = cfgTauIdEffAnalyzer.getParameter<int>("allEvents_DBS");
  if ( allEvents_DBS > 0 ) {
    histogramEventCounter->SetBinContent(1, allEvents_DBS);
  } else {
    histogramEventCounter->SetBinContent(1, -1.);
  }
  
  double xSection = cfgTauIdEffAnalyzer.getParameter<double>("xSection");
  double intLumiData = cfgTauIdEffAnalyzer.getParameter<double>("intLumiData");

  int    numEvents_processed                     = 0; 
  double numEventsWeighted_processed             = 0.;
  int    numEvents_passedTrigger                 = 0;
  double numEventsWeighted_passedTrigger         = 0.;
  int    numEvents_passedDiMuonVeto              = 0;
  double numEventsWeighted_passedDiMuonVeto      = 0.;
  int    numEvents_passedDiMuTauPairVeto         = 0;
  double numEventsWeighted_passedDiMuTauPairVeto = 0.;

  edm::RunNumber_t lastLumiBlock_run = -1;
  edm::LuminosityBlockNumber_t lastLumiBlock_ls = -1;

  double intLumiData_analyzed = 0.;
  edm::InputTag srcLumiProducer = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcLumiProducer");

  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauIdEffAnalyzer") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << "opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
      double evtWeight = 1.0;
      for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	    srcWeight != srcWeights.end(); ++srcWeight ) {
	edm::Handle<double> weight;
	evt.getByLabel(*srcWeight, weight);
	evtWeight *= (*weight);
      }

//--- apply trigger efficiency correction (to MC only)
      typedef std::vector<pat::MET> PATMETCollection;
      edm::Handle<PATMETCollection> caloMEt;
      evt.getByLabel(srcCaloMEt, caloMEt);
      if ( caloMEt->size() != 1 )
	throw cms::Exception("FWLiteTauIdEffAnalyzer")
	  << "Failed to find unique CaloMEt object !!\n";
      //std::cout << " " << srcCaloMEt.label() << ": " << caloMEt->front().pt() << std::endl;

      if ( !isData ) {
	double triggerEffCorrection_value = triggerEffCorrection->Eval(caloMEt->front().pt());
	//std::cout << " triggerEffCorrection_value = " << triggerEffCorrection_value << std::endl;
	evtWeight *= triggerEffCorrection_value;
      }

//--- check if current event is within specified run-range
//   (this check is important in case triggers changed or became active/inactive **during** a data-taking period)
      if ( (firstRun != -1 && (int)evt.id().run() < firstRun) || 
	   (lastRun  != -1 && (int)evt.id().run() > lastRun ) ) continue;

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

      //std::cout << "processing run = " << evt.id().run() << ":" 
      //	  << " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << std::endl;

//--- check if new luminosity section has started;
//    if so, retrieve number of events contained in this luminosity section before skimming
      if ( !(evt.id().run() == lastLumiBlock_run && evt.luminosityBlock() == lastLumiBlock_ls) ) {
	const fwlite::LuminosityBlock& ls = evt.getLuminosityBlock();
	edm::Handle<edm::MergeableCounter> numEvents_skimmed;
	ls.getByLabel(srcEventCounter, numEvents_skimmed);
	if ( numEvents_skimmed.isValid() ) histogramEventCounter->Fill(1, numEvents_skimmed->value);
	lastLumiBlock_run = evt.id().run();
	lastLumiBlock_ls = evt.luminosityBlock();

	if ( isData ) {
	  edm::Handle<LumiSummary> lumiSummary;
	  edm::InputTag srcLumiProducer("lumiProducer");
	  ls.getByLabel(srcLumiProducer, lumiSummary);
	  intLumiData_analyzed += lumiSummary->intgRecLumi();
	}
      }

//--- fill "dummy" histogram counting number of processed events
      histogramEventCounter->Fill(2);

//--- check that event has passed triggers
      bool anyHLTpath_passed = false;
      if ( hltPaths.size() == 0 ) {
	anyHLTpath_passed = true;
      } else {
	checkHLTpaths(evt, hltPaths, srcHLTresults, NULL, &anyHLTpath_passed);
      }

      if ( !anyHLTpath_passed ) continue;

      bool anyL1bit_passed = false;
      if ( l1Bits.size() == 0 ) {
	anyL1bit_passed = true;
      } else {
	checkL1Bits(evt, l1Bits, srcL1GtReadoutRecord, srcL1GtObjectMapRecord, NULL, &anyL1bit_passed);
      }

      if ( !anyL1bit_passed ) continue;

      ++numEvents_passedTrigger;
      numEventsWeighted_passedTrigger += evtWeight;

//--- require event to contain only one "good quality" muon
      typedef std::vector<pat::Muon> PATMuonCollection;
      edm::Handle<PATMuonCollection> goodMuons;
      evt.getByLabel(srcGoodMuons, goodMuons);
      size_t numGoodMuons = goodMuons->size();
	
      if ( !(numGoodMuons <= 1) ) continue;
      ++numEvents_passedDiMuonVeto;
      numEventsWeighted_passedDiMuonVeto += evtWeight;

//--- require event to contain exactly one muon + tau-jet pair
//    passing the selection criteria for region "ABCD"
      edm::Handle<PATMuTauPairCollection> muTauPairs;
      evt.getByLabel(srcMuTauPairs, muTauPairs);           

      unsigned numMuTauPairsABCD = 0;
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {
	pat::strbitset evtSelFlags;
	if ( selectorABCD->operator()(*muTauPair, caloMEt->front(), evtSelFlags) ) ++numMuTauPairsABCD;
      }
      
      if ( !(numMuTauPairsABCD <= 1) ) continue;
      ++numEvents_passedDiMuTauPairVeto;
      numEventsWeighted_passedDiMuTauPairVeto += evtWeight;

//--- determine number of vertices reconstructed in the event
//   (needed to parametrize dependency of tau id. efficiency on number of pile-up interactions)
      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      size_t numVertices = vertices->size();
      
//--- check L1 bits for trigger efficiency control plots
      std::map<std::string, bool> plot_triggerBits_passed;
      if ( plot_l1Bits.size() > 0 ) {
	checkL1Bits(evt, plot_l1Bits, srcL1GtReadoutRecord, srcL1GtObjectMapRecord, &plot_triggerBits_passed, NULL);
      }	
      if ( plot_hltPaths.size() > 0 ) {
	checkHLTpaths(evt, plot_hltPaths, srcHLTresults, &plot_triggerBits_passed, NULL);
      }

//--- iterate over collection of muon + tau-jet pairs:
//    check which region muon + tau-jet pair is selected in,
//    fill histograms for that region
      if ( requireUniqueMuTauPair && muTauPairs->size () > 1 ) continue;  
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {

//--- determine type of particle matching reconstructed tau-jet candidate
//    on generator level (used in case of Ztautau or Zmumu Monte Carlo samples only,
//    in order to distinguish between jet --> tau fakes, muon --> tau fakes and genuine taus)
	int genMatchType = kUnmatched;
	if ( fillGenMatchHistograms ) {
	  edm::Handle<reco::GenParticleCollection> genParticles;
	  evt.getByLabel(srcGenParticles, genParticles);	  
	  genMatchType = getGenMatchType(*muTauPair, *genParticles);
	}

	for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	      regionEntry != regionEntries.end(); ++regionEntry ) {	  
	  double evtWeight_region = evtWeight;
	  if ( muonIsoProbExtractor && (*regionEntry)->region_.find("_mW") != std::string::npos ) 
	    evtWeight_region *= (*muonIsoProbExtractor)(*muTauPair->leg1());
	  (*regionEntry)->analyze(evt, *muTauPair, caloMEt->front(), 
				  numVertices, plot_triggerBits_passed, genMatchType, evtWeight_region);

	  //pat::strbitset evtSelFlags;
	  //if ( (*regionEntry)->region_ == "D1p" && (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	  //  std::cout << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() 
	  //	        << " (weight = " << evtWeight << ")" << std::endl;
	  //}
	}
      }
    }

//--- close input file
    delete inputFile;
  }

//--- scale histograms taken from Monte Carlo simulation
//    according to cross-section times luminosity
  if ( !isData ) {
    double mcScaleFactor = (intLumiData*xSection)/(double)allEvents_DBS;
    std::cout << " intLumiData = " << intLumiData << std::endl;
    std::cout << " xSection = " << xSection << std::endl;
    std::cout << " allEvents_DBS = " << allEvents_DBS << std::endl;
    std::cout << "--> scaling histograms by factor = " << mcScaleFactor
	      << " according to cross-section times luminosity." << std::endl;

//--- apply correction to scale-factor in order to account for events lost, 
//    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs
    std::cout << " allEvents_processed = " << histogramEventCounter->GetBinContent(2) << std::endl;
    double lostStatCorrFactor = 1.;
    if ( histogramEventCounter->GetBinContent(1) > histogramEventCounter->GetBinContent(2) && 
	 histogramEventCounter->GetBinContent(2) > 0.                                      ) {
      lostStatCorrFactor = histogramEventCounter->GetBinContent(1)/histogramEventCounter->GetBinContent(2);
      std::cout << "--> scaling histograms by additional factor = " << lostStatCorrFactor
		<< " to account for events lost," << std::endl; 
      std::cout << "    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs." << std::endl;
    }

    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {  
      (*regionEntry)->histogramsUnbinned_->histManager_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = (*regionEntry)->histogramEntriesBinned_.begin();
	    histManagerEntry != (*regionEntry)->histogramEntriesBinned_.end(); ++histManagerEntry ) {
	(*histManagerEntry)->histManager_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
	(*histManagerEntry)->histManagerJetToTauFake_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
	(*histManagerEntry)->histManagerMuToTauFake_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
        (*histManagerEntry)->histManagerGenTau_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      }
    }
  }

  std::cout << "<FWLiteTauIdEffAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;
  std::cout << " numEvents_passedTrigger: " << numEvents_passedTrigger 
	    << " (weighted = " << numEventsWeighted_passedTrigger << ")" << std::endl;
  std::cout << " numEvents_passedDiMuonVeto: " << numEvents_passedDiMuonVeto 
	    << " (weighted = " << numEventsWeighted_passedDiMuonVeto << ")" << std::endl;
  std::cout << " numEvents_passedDiMuTauPairVeto: " << numEvents_passedDiMuTauPairVeto
	    << " (weighted = " << numEventsWeighted_passedDiMuTauPairVeto << ")" << std::endl;
  std::string lastTauIdName = "";
  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    if ( (*regionEntry)->tauIdName_ != lastTauIdName ) 
      std::cout << " numMuTauPairs_selected, " << (*regionEntry)->tauIdName_ << std::endl;
    std::cout << "  region " << (*regionEntry)->region_ << ":" 
	      << " " << (*regionEntry)->numMuTauPairs_selected_ 
	      << " (weighted = " << (*regionEntry)->numMuTauPairsWeighted_selected_ << ")" << std::endl;
    lastTauIdName = (*regionEntry)->tauIdName_;
  }

  delete muonIsoProbExtractor;

  delete triggerEffCorrection;

//--- close ASCII files containing run + event numbers of events selected in different regions
  for ( std::vector<regionEntryType*>::iterator it = regionEntries.begin();
	it != regionEntries.end(); ++it ) {
    delete (*it);
  }
  
  if ( isData ) {
    std::cout << " intLumiData = " << intLumiData << " pb" << std::endl;
    // CV: luminosity is recorded in some 'weird' units,
    //     needs to be multiplied by factor 0.10 in order to be in units of pb^-1
    std::cout << " intLumiData_analyzed = " << intLumiData_analyzed*1.e-6*0.10 << " pb" << std::endl;
  }

  clock.Show("FWLiteTauIdEffAnalyzer");

  return 0;
}
