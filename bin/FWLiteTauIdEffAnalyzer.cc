
/** \executable FWLiteTauIdEffAnalyzer
 *
 * Apply event selections for ABCD regions 
 * and fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.31 $
 *
 * $Id: FWLiteTauIdEffAnalyzer.cc,v 1.31 2012/06/18 19:15:34 veelken Exp $
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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
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
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

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
      fillGenMatchHistograms_(fillGenMatchHistograms),
      histManagerJetToTauFake_(0),
      histManagerMuToTauFake_(0),
      histManagerGenTau_(0)
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
		      size_t numJets, size_t numJets_bTagged,
		      size_t numVertices, const std::map<std::string, bool>& plot_triggerBits_passed, int genMatchType, double weight)
  {
    if ( x > min_ && x <= max_ ) {
      histManager_->fillHistograms(muTauPair, caloMEt, 
				   numJets, numJets_bTagged, 
				   numVertices, plot_triggerBits_passed, weight);

      if ( fillGenMatchHistograms_ ) {
	if      ( genMatchType == kJetToTauFakeMatched ) 
	  histManagerJetToTauFake_->fillHistograms(muTauPair, caloMEt, 
						   numJets, numJets_bTagged, 
						   numVertices, plot_triggerBits_passed, weight);
	else if ( genMatchType == kMuToTauFakeMatched  ) 
	  histManagerMuToTauFake_->fillHistograms(muTauPair, caloMEt, 
						  numJets, numJets_bTagged, 
						  numVertices, plot_triggerBits_passed, weight);
	else if ( genMatchType == kGenTauHadMatched    ||
		  genMatchType == kGenTauOtherMatched  ) 
	  histManagerGenTau_->fillHistograms(muTauPair, caloMEt, 
					     numJets, numJets_bTagged, 
					     numVertices, plot_triggerBits_passed, weight);
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
		  bool fillGenMatchHistograms, bool fillControlPlots, const vstring& plot_triggerBits,
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
    cfgHistManager.addParameter<bool>("fillControlPlots", fillControlPlots);
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
	       size_t numJets, size_t numJets_bTagged,
	       size_t numVertices, const std::map<std::string, bool>& plot_triggerBits_passed, int genMatchType, double evtWeight)
  {
    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, caloMEt, numJets_bTagged, evtSelFlags) ) {
//--- fill histograms for "inclusive" tau id. efficiency measurement
      histogramsUnbinned_->fillHistograms(0., muTauPair, caloMEt, 
					  numJets, numJets_bTagged,
					  numVertices, plot_triggerBits_passed, genMatchType, evtWeight);

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
	(*histManagerEntry)->fillHistograms(x, muTauPair, caloMEt, 
					    numJets, numJets_bTagged,
					    numVertices, plot_triggerBits_passed, genMatchType, evtWeight);
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
  vstring hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("hltPaths");
  edm::InputTag srcCaloMEt = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcCaloMEt");
  edm::InputTag srcGoodMuons = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGoodMuons");
  edm::InputTag srcJets = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcJets");
  edm::ParameterSet cfgJetId;
  cfgJetId.addParameter<std::string>("version", "FIRSTDATA");
  cfgJetId.addParameter<std::string>("quality", "LOOSE");
  PFJetIDSelectionFunctor jetId(cfgJetId);
  edm::InputTag srcVertices = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcVertices");
  edm::InputTag srcGenParticles = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcGenParticles");
  bool fillGenMatchHistograms = cfgTauIdEffAnalyzer.getParameter<bool>("fillGenMatchHistograms");
  bool fillControlPlots = cfgTauIdEffAnalyzer.getParameter<bool>("fillControlPlots");
  typedef std::vector<int> vint;
  vint skipPdgIdsGenParticleMatch = cfgTauIdEffAnalyzer.getParameter<vint>("skipPdgIdsGenParticleMatch");  
  vstring plot_hltPaths = cfgTauIdEffAnalyzer.getParameter<vstring>("plot_hltPaths");
  vstring plot_triggerBits;
  for ( vstring::const_iterator plot_hltPath = plot_hltPaths.begin();
	plot_hltPath != plot_hltPaths.end(); ++plot_hltPath ) {
    plot_triggerBits.push_back(getHLTpath_key(*plot_hltPath));
  }
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffAnalyzer.getParameter<vInputTag>("weights");
  double minWeight = cfgTauIdEffAnalyzer.getParameter<double>("minWeight");
  double maxWeight = cfgTauIdEffAnalyzer.getParameter<double>("maxWeight");
  std::string sysShift = cfgTauIdEffAnalyzer.exists("sysShift") ?
    cfgTauIdEffAnalyzer.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";
  int shiftCaloMEtResponse = 0; //  0: no shift applied
                                // +1: shift reconstructed CaloMEt by +15% and corresponding turn-on curve for L1_ETM20 trigger
                                // -1: shift reconstructed CaloMEt by -15% and corresponding turn-on curve for L1_ETM20 trigger
  if      ( sysShift == "CaloMEtResponseUp"   ) shiftCaloMEtResponse = +1;
  else if ( sysShift == "CaloMEtResponseDown" ) shiftCaloMEtResponse = -1;
  
  edm::InputTag srcEventCounter = cfgTauIdEffAnalyzer.getParameter<edm::InputTag>("srcEventCounter");

  PATMuonLUTvalueExtractorFromKNN* muonIsoProbExtractor = 0;
  bool applyMuonIsoWeights = false;
  if ( cfgTauIdEffAnalyzer.exists("muonIsoProbExtractor") ) {
    edm::ParameterSet cfgMuonIsoProbExtractor = cfgTauIdEffAnalyzer.getParameter<edm::ParameterSet>("muonIsoProbExtractor");
    muonIsoProbExtractor = new PATMuonLUTvalueExtractorFromKNN(cfgMuonIsoProbExtractor);
    applyMuonIsoWeights = cfgTauIdEffAnalyzer.getParameter<bool>("applyMuonIsoWeights");
  }

  TF1* triggerEffCorrection = new TF1("triggerEffCorrection", &integralCrystalBall_data_div_mc, 0., 1.e+6, 10);
  double triggerEffCorr_parameter_data[] = {
    // CV: HLT_IsoMu15_eta2p1_L1ETM20 efficiency correction parameters for 2012 run A
    3.56934e+01, 8.95771e+00, 3.81006e+01, 3.29700e-02, 9.88290e-01, // Zmumu Data
  };
  double triggerEffCorr_parameter_mc[] = {
    3.51723e+01, 9.44549e+00, 2.13250e+01, 1.64246e+02, 1.00001e+00  // Zmumu Spring'12 MC (pile-up reweighted)
  };
  double triggerEffCorr_parameter_mc_shiftUp[] = {
    2.99172e+01, 8.04091e+00, 4.13981e+01, 2.10489e+00, 9.99198e-01 // Zmumu Spring'12 MC (raw CaloMEt shifted by +15% wrt. corrected CaloMEt)
  };
  double triggerEffCorr_parameter_mc_shiftDown[] = {
    4.04481e+01, 1.08623e+01, 2.45109e+01, 1.64435e+02, 1.00001e+00 // Zmumu Spring'12 MC (raw CaloMEt shifted by -15% wrt. corrected CaloMEt)
  };
  int numParameter_data_or_mc = (triggerEffCorrection->GetNpar() / 2);
  assert(2*numParameter_data_or_mc == triggerEffCorrection->GetNpar());
  for ( int iPar = 0; iPar < (2*numParameter_data_or_mc); ++iPar ) {
    if ( iPar <= numParameter_data_or_mc ) { // data
      triggerEffCorrection->SetParameter(iPar, triggerEffCorr_parameter_data[iPar]);
    } else {                                 // mc (either central value or CaloMEt response shifted up/down)
      if      ( shiftCaloMEtResponse ==  0 ) triggerEffCorrection->SetParameter(iPar, triggerEffCorr_parameter_mc[iPar - numParameter_data_or_mc]);
      else if ( shiftCaloMEtResponse == +1 ) triggerEffCorrection->SetParameter(iPar, triggerEffCorr_parameter_mc_shiftUp[iPar - numParameter_data_or_mc]);
      else if ( shiftCaloMEtResponse == -1 ) triggerEffCorrection->SetParameter(iPar, triggerEffCorr_parameter_mc_shiftDown[iPar - numParameter_data_or_mc]);
      else assert(0);      
    }
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
			    fillGenMatchHistograms, fillControlPlots, plot_triggerBits, selEventsFileName);
      regionEntries.push_back(regionEntry);

      // tau+ candidates only
      //regionEntryType* regionEntry_plus = 
      //  new regionEntryType(fs, process, std::string(*region).append("+"), tauIdDiscriminators, tauIdName, 
      //		      sysShift, cfgBinning, svFitMassHypothesis, 
      //		      tauChargeMode, disableTauCandPreselCuts, cfgEventSelCuts, 
      //		      fillGenMatchHistograms, fillControlPlots, plot_triggerBits, selEventsFileName);
      //regionEntries.push_back(regionEntry_plus);
      //
      // tau- candidates only
      //regionEntryType* regionEntry_minus = 
      //  new regionEntryType(fs, process, std::string(*region).append("-"), tauIdDiscriminators, tauIdName, 
      //		      sysShift, cfgBinning, svFitMassHypothesis, 
      //		      tauChargeMode, disableTauCandPreselCuts, cfgEventSelCuts, 
      //		      fillGenMatchHistograms, fillControlPlots, plot_triggerBits, selEventsFileName);
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
      if ( evtWeight < minWeight ) evtWeight = minWeight;
      if ( evtWeight > maxWeight ) evtWeight = maxWeight;
      
//--- apply trigger efficiency correction (to MC only)
      typedef std::vector<pat::MET> PATMETCollection;
      edm::Handle<PATMETCollection> caloMETs;
      evt.getByLabel(srcCaloMEt, caloMETs);
      if ( caloMETs->size() != 1 )
	throw cms::Exception("FWLiteTauIdEffAnalyzer")
	  << "Failed to find unique CaloMEt object !!\n";
      pat::MET caloMEt = caloMETs->front();
      //std::cout << " " << srcCaloMEt.label() << ": " << caloMEt.pt() << std::endl;
      if ( shiftCaloMEtResponse != 0 ) {
        reco::Candidate::LorentzVector caloMEtP4 = caloMEt.p4();
        if      ( shiftCaloMEtResponse == +1 ) caloMEtP4 *= 1.15;
        else if ( shiftCaloMEtResponse == -1 ) caloMEtP4 *= 0.85;
        else assert(0);
        caloMEt.setP4(caloMEtP4);
      }

      if ( !isData ) {
	double triggerEffCorrection_value = triggerEffCorrection->Eval(caloMEt.pt());
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

      unsigned numMuTauPairsABCD = 0; // Note: no b-jet veto applied
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {
	pat::strbitset evtSelFlags;
	if ( selectorABCD->operator()(*muTauPair, caloMEt, 0, evtSelFlags) ) ++numMuTauPairsABCD;
      }
      
      if ( !(numMuTauPairsABCD <= 1) ) continue;
      ++numEvents_passedDiMuTauPairVeto;
      numEventsWeighted_passedDiMuTauPairVeto += evtWeight;

      edm::Handle<pat::JetCollection> jets;
      evt.getByLabel(srcJets, jets);         
      
//--- determine number of vertices reconstructed in the event
//   (needed to parametrize dependency of tau id. efficiency on number of pile-up interactions)
      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      size_t numVertices = vertices->size();
      
//--- check L1 bits for trigger efficiency control plots
      std::map<std::string, bool> plot_triggerBits_passed;
      if ( plot_hltPaths.size() > 0 ) {
	checkHLTpaths(evt, plot_hltPaths, srcHLTresults, &plot_triggerBits_passed, NULL);
      }

//--- iterate over collection of muon + tau-jet pairs:
//    check which region muon + tau-jet pair is selected in,
//    fill histograms for that region
      if ( requireUniqueMuTauPair && muTauPairs->size () > 1 ) continue;  
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair ) {

//--- require event to contain to b-jets
//   (not overlapping with muon or tau-jet candidate)
	size_t numJets         = 0;
	size_t numJets_bTagged = 0;
	for ( pat::JetCollection::const_iterator jet = jets->begin();
	      jet != jets->end(); ++jet ) {
	  if ( deltaR(jet->p4(), muTauPair->leg1()->p4()) > 0.5 &&
	       deltaR(jet->p4(), muTauPair->leg2()->p4()) > 0.5 &&
               jet->pt() > 30. && TMath::Abs(jet->eta()) < 2.4 && jetId(*jet) ) {
	    ++numJets;
	    if ( jet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.679 ) ++numJets_bTagged; // "medium" WP
	  }
	}

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
	  if ( muonIsoProbExtractor && applyMuonIsoWeights && (*regionEntry)->region_.find("_mW") != std::string::npos ) 
	    evtWeight_region *= (*muonIsoProbExtractor)(*muTauPair->leg1());
	  (*regionEntry)->analyze(evt, *muTauPair, caloMEt, 
				  numJets, numJets_bTagged,
				  numVertices, plot_triggerBits_passed, genMatchType, evtWeight_region);
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
    if ( histogramEventCounter->GetBinContent(1) != histogramEventCounter->GetBinContent(2) && 
	 histogramEventCounter->GetBinContent(2) > 0.                                       ) {
      lostStatCorrFactor = histogramEventCounter->GetBinContent(1)/histogramEventCounter->GetBinContent(2);
      std::cout << "--> scaling histograms by additional factor = " << lostStatCorrFactor
		<< " to account for events lost," << std::endl; 
      std::cout << "    due to aborted skimming/crab or PAT-tuple production/lxbatch jobs." << std::endl;
    }

    std::vector<histManagerEntryType*> histManagerEntriesToScale;
    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {  
      histManagerEntriesToScale.push_back((*regionEntry)->histogramsUnbinned_);
      for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = (*regionEntry)->histogramEntriesBinned_.begin();
	    histManagerEntry != (*regionEntry)->histogramEntriesBinned_.end(); ++histManagerEntry ) {
	histManagerEntriesToScale.push_back(*histManagerEntry);
      }
    }
    for ( std::vector<histManagerEntryType*>::iterator histManagerEntry = histManagerEntriesToScale.begin();
	  histManagerEntry != histManagerEntriesToScale.end(); ++histManagerEntry ) {
      (*histManagerEntry)->histManager_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      if ( (*histManagerEntry)->histManagerJetToTauFake_ ) 
	(*histManagerEntry)->histManagerJetToTauFake_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      if ( (*histManagerEntry)->histManagerMuToTauFake_ )
	(*histManagerEntry)->histManagerMuToTauFake_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
      if ( (*histManagerEntry)->histManagerGenTau_ )
	(*histManagerEntry)->histManagerGenTau_->scaleHistograms(mcScaleFactor*lostStatCorrFactor);
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
