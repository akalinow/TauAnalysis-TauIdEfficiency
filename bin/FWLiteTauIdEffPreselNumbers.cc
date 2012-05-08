
/** \executable FWLiteTauIdEffPreselNumbers
 *
 * Determine efficiency of "leading" track finding, leading track Pt and loose (PF)isolation requirements 
 * applied in preselection of tau-jet candidates considered for tau id. efficiency measurement.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.16 $
 *
 * $Id: FWLiteTauIdEffPreselNumbers.cc,v 1.16 2012/02/02 09:03:32 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/Event.h"

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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffCutFlowTable.h"
#include "TauAnalysis/TauIdEfficiency/interface/tauIdEffAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

typedef std::vector<std::string> vstring;
typedef std::vector<bool> vbool;

struct cutFlowEntryType
{
  cutFlowEntryType(const edm::ParameterSet& cfg)
  {
    binVariable_ = cfg.getParameter<std::string>("binVariable");

    std::string label = cfg.getParameter<std::string>("label");

    vstring selectionNamesReversed = cfg.getParameter<vstring>("selectionNamesReversed");

    edm::ParameterSet cfgTauHadMatched = cfg;
    std::string labelTauHadMatched = std::string(label).append("TauHadMatched");
    cfgTauHadMatched.addParameter<std::string>("label", labelTauHadMatched);
    cutFlowTauHadMatched_ = new TauIdEffCutFlowTable(cfgTauHadMatched);
    edm::ParameterSet cfgTauHadMatchedReversed = cfg;
    std::string labelTauHadMatchedReversed = std::string(label).append("TauHadMatchedReversed");
    cfgTauHadMatchedReversed.addParameter<std::string>("label", labelTauHadMatchedReversed);
    cfgTauHadMatchedReversed.addParameter<vstring>("selectionNames", selectionNamesReversed);
    cutFlowTauHadMatchedReversed_ = new TauIdEffCutFlowTable(cfgTauHadMatchedReversed);

    edm::ParameterSet cfgFakeTauMatched = cfg;
    std::string labelFakeTauMatched = std::string(label).append("FakeTauMatched");
    cfgFakeTauMatched.addParameter<std::string>("label", labelFakeTauMatched);
    cutFlowFakeTauMatched_ = new TauIdEffCutFlowTable(cfgFakeTauMatched);
    edm::ParameterSet cfgFakeTauMatchedReversed = cfg;
    std::string labelFakeTauMatchedReversed = std::string(label).append("FakeTauMatchedReversed");
    cfgFakeTauMatchedReversed.addParameter<std::string>("label", labelFakeTauMatchedReversed);
    cfgFakeTauMatchedReversed.addParameter<vstring>("selectionNames", selectionNamesReversed);
    cutFlowFakeTauMatchedReversed_ = new TauIdEffCutFlowTable(cfgFakeTauMatchedReversed);

    edm::ParameterSet cfgNoMatchingApplied = cfg;
    std::string labelNoMatchingApplied = std::string(label).append("NoMatchingApplied");
    cfgNoMatchingApplied.addParameter<std::string>("label", labelNoMatchingApplied);
    cutFlowNoMatchingApplied_ = new TauIdEffCutFlowTable(cfgNoMatchingApplied);
    edm::ParameterSet cfgNoMatchingAppliedReversed = cfg;
    std::string labelNoMatchingAppliedReversed = std::string(label).append("NoMatchingAppliedReversed");
    cfgNoMatchingAppliedReversed.addParameter<std::string>("label", labelNoMatchingAppliedReversed);
    cfgNoMatchingAppliedReversed.addParameter<vstring>("selectionNames", selectionNamesReversed);
    cutFlowNoMatchingAppliedReversed_ = new TauIdEffCutFlowTable(cfgNoMatchingAppliedReversed);
  }
  ~cutFlowEntryType()
  {
    delete cutFlowTauHadMatched_;
    delete cutFlowTauHadMatchedReversed_;
    delete cutFlowFakeTauMatched_;
    delete cutFlowFakeTauMatchedReversed_;
    delete cutFlowNoMatchingApplied_;
    delete cutFlowNoMatchingAppliedReversed_;
  }
  void bookCutFlowTables(TFileDirectory& dir)
  {
    cutFlowTauHadMatched_->bookCutFlowTable(dir);
    cutFlowTauHadMatchedReversed_->bookCutFlowTable(dir);
    cutFlowFakeTauMatched_->bookCutFlowTable(dir);
    cutFlowFakeTauMatchedReversed_->bookCutFlowTable(dir);
    cutFlowNoMatchingApplied_->bookCutFlowTable(dir);
    cutFlowNoMatchingAppliedReversed_->bookCutFlowTable(dir);
  }
  void fillCutFlowTables(double x, 
			 const vbool& selectionFlags, const vbool& selectionFlagsReversed, 
			 int genMatchType, double genTauCharge, double recTauCharge,
			 double evtWeight)
  {
    if        ( genMatchType == kGenTauHadMatched   ) {
      cutFlowTauHadMatched_->fillCutFlowTable(x, selectionFlags, evtWeight);
      cutFlowTauHadMatchedReversed_->fillCutFlowTable(x, selectionFlagsReversed, evtWeight);
    } else if ( genMatchType == kGenTauOtherMatched ) {
      cutFlowFakeTauMatched_->fillCutFlowTable(x, selectionFlags, evtWeight);
      cutFlowFakeTauMatchedReversed_->fillCutFlowTable(x, selectionFlagsReversed, evtWeight);
    } 

    cutFlowNoMatchingApplied_->fillCutFlowTable(x, selectionFlags, evtWeight);
    cutFlowNoMatchingAppliedReversed_->fillCutFlowTable(x, selectionFlagsReversed, evtWeight);
  }

  std::string binVariable_;

  TauIdEffCutFlowTable* cutFlowTauHadMatched_;
  TauIdEffCutFlowTable* cutFlowTauHadMatchedReversed_;

  TauIdEffCutFlowTable* cutFlowFakeTauMatched_;
  TauIdEffCutFlowTable* cutFlowFakeTauMatchedReversed_;
  TauIdEffCutFlowTable* cutFlowNoMatchingApplied_;
  TauIdEffCutFlowTable* cutFlowNoMatchingAppliedReversed_;
};

struct regionEntryType
{
  regionEntryType(fwlite::TFileService& fs,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName, const std::string& sysShift,
                  const edm::ParameterSet& cfgBinning,
		  const edm::ParameterSet& cfgEventSelCuts)
    : process_(process),
      region_(region),
      tauIdDiscriminators_(tauIdDiscriminators),
      tauIdName_(tauIdName),      
      sysShift_(sysShift),
      numPreselCuts_(8),
      numTauIdDiscriminators_(tauIdDiscriminators.size()),      
      selector_(0),
      cutFlowUnbinned_(0),
      numMuTauPairs_selected_(0),
      numMuTauPairsWeighted_selected_(0.)
  {
    edm::ParameterSet cfgSelector = cfgEventSelCuts;
    cfgSelector.addParameter<std::string>("region", region_);
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators_);
    cfgSelector.addParameter<std::string>("tauChargeMode", "tauLeadTrackCharge");
    cfgSelector.addParameter<bool>("disableTauCandPreselCuts", false);

    selector_ = new TauIdEffEventSelector(cfgSelector);

//--- disable preselection cuts applied on tau-jet candidates
    selector_->tauLeadTrackPtMin_      =  -1.e+3; 
    selector_->tauAbsIsoMax_           =  +1.e+3;
    selector_->muTauPairAbsDzMax_      =  +1.e+3;
    selector_->muTauPairChargeProdMin_ =  -1.e+3;
    selector_->muTauPairChargeProdMax_ =  +1.e+3; 

    muTauPairChargeProdMin_ = -1.e+3;
    muTauPairChargeProdMax_ = +1.e+3; 
    if        ( region_           == "ABCD"            ) {
      // nothing to be done yet...
      // (simply use default cuts)
    } else if ( region_.find("A") != std::string::npos ) {
      muTauPairChargeProdMax_ = -0.5;
    } else if ( region_.find("B") != std::string::npos ) {
      muTauPairChargeProdMin_ = +0.5;
    } else if ( region_.find("C") != std::string::npos ) {
      muTauPairChargeProdMax_ = -0.5;
    } else if ( region_.find("D") != std::string::npos ) {
      muTauPairChargeProdMin_ = +0.5;
    } else {
      throw cms::Exception("regionEntryType") 
	<< "Invalid region = " << region_ << " !!\n";
    }
           
    edm::ParameterSet cfgCutFlowTable = cfgBinning;
    cfgCutFlowTable.addParameter<std::string>("process", process_);
    cfgCutFlowTable.addParameter<std::string>("region", region_);
    cfgCutFlowTable.addParameter<std::string>("tauIdDiscriminator", tauIdName_);
    label_ = "";
    if ( sysShift != "CENTRAL_VALUE" ) label_.append("_").append(sysShift_);
    cfgCutFlowTable.addParameter<std::string>("label", label_);

    selectionNames_.resize(numPreselCuts_ + numTauIdDiscriminators_);
    selectionNames_[0] = std::string(region);
    selectionNames_[1] = "leadTrackFinding";
    selectionNames_[2] = "leadTrackPtCut";
    selectionNames_[3] = "loosePFIso";
    selectionNames_[4] = "eVeto";
    selectionNames_[5] = "muVeto";
    selectionNames_[6] = "muTauPairAbsDz";
    selectionNames_[7] = "muTauPairChargeProd";
    for ( int iTauIdDiscriminator = 0; iTauIdDiscriminator < numTauIdDiscriminators_; ++iTauIdDiscriminator ) {
      selectionNames_[numPreselCuts_ + iTauIdDiscriminator] = tauIdDiscriminators_[iTauIdDiscriminator];
    }
    cfgCutFlowTable.addParameter<vstring>("selectionNames", selectionNames_);

    selectionNamesReversed_.resize(selectionNames_.size());
    selectionNamesReversed_[0] = selectionNames_[0];
    for ( int iTauIdDiscriminator = 0; iTauIdDiscriminator < numTauIdDiscriminators_; ++iTauIdDiscriminator ) {
      selectionNamesReversed_[1 + iTauIdDiscriminator] = selectionNames_[numPreselCuts_ + iTauIdDiscriminator];
    }
    for ( int iPreselCut = 1; iPreselCut < numPreselCuts_; ++iPreselCut ) {
      selectionNamesReversed_[numTauIdDiscriminators_ + iPreselCut] = selectionNames_[iPreselCut];
    }
    cfgCutFlowTable.addParameter<vstring>("selectionNamesReversed", selectionNamesReversed_);

    tauIdFlags_.resize(numPreselCuts_ + numTauIdDiscriminators_);
    tauIdFlagsReversed_.resize(numPreselCuts_ + numTauIdDiscriminators_);

    TFileDirectory dir = fs.mkdir("presel");

    edm::ParameterSet cfgUnbinned_bin0;
    cfgUnbinned_bin0.addParameter<std::string>("subdir", "");
    cfgUnbinned_bin0.addParameter<double>("min", -0.5);
    cfgUnbinned_bin0.addParameter<double>("max", +0.5);
    typedef std::vector<edm::ParameterSet> vParameterSet;
    vParameterSet cfgUnbinned_bins;
    cfgUnbinned_bins.push_back(cfgUnbinned_bin0);
    edm::ParameterSet cfgUnbinned = cfgCutFlowTable;
    cfgUnbinned.addParameter<vParameterSet>("binning", cfgUnbinned_bins);
    cfgUnbinned.addParameter<std::string>("binVariable", "");
    cutFlowUnbinned_ = new cutFlowEntryType(cfgUnbinned);
    cutFlowUnbinned_->bookCutFlowTables(dir);

    vstring binVariableNames = cfgBinning.getParameterNamesForType<vParameterSet>();
    for ( vstring::const_iterator binVariableName = binVariableNames.begin();
	  binVariableName != binVariableNames.end(); ++binVariableName ) {
      vParameterSet cfgBinVariableBins = cfgBinning.getParameter<vParameterSet>(*binVariableName);
      edm::ParameterSet cfgBinVariable = cfgCutFlowTable;
      cfgBinVariable.addParameter<vParameterSet>("binning", cfgBinVariableBins);
      cfgBinVariable.addParameter<std::string>("binVariable", *binVariableName);
      cutFlowEntryType* cutFlowEntry = new cutFlowEntryType(cfgBinVariable);
      cutFlowEntry->bookCutFlowTables(dir);
      cutFlowEntriesBinned_.push_back(cutFlowEntry);
    }
  }
  ~regionEntryType()
  {
    delete selector_;

    delete cutFlowUnbinned_;
    for ( std::vector<cutFlowEntryType*>::iterator it = cutFlowEntriesBinned_.begin();
	  it != cutFlowEntriesBinned_.end(); ++it ) {
      delete (*it);
    }
  }
  void analyze(const PATMuTauPair& muTauPair, 
	       int genMatchType, double genTauCharge, double recTauCharge,
	       const pat::MET& caloMEt, 
	       size_t numVertices, 
	       double evtWeight)
  {
    pat::strbitset evtSelFlags;
    if ( selector_->operator()(muTauPair, caloMEt, evtSelFlags) ) {

//--- set flags indicating whether tau-jet candidate passes 
//    "leading" track finding, leading track Pt and loose (PF)isolation requirements 
//    plus tau id. discriminators
      tauIdFlags_[0] = true;
      tauIdFlags_[1] = (muTauPair.leg2()->userFloat("hasLeadTrack")       > 0.5);
      tauIdFlags_[2] = (muTauPair.leg2()->userFloat("leadTrackPt")        > 5.0);      
      tauIdFlags_[3] = (muTauPair.leg2()->userFloat("preselLoosePFIsoPt") < 2.5);
      tauIdFlags_[4] = (muTauPair.leg2()->userFloat("PFElectronMVA")      < 0.6);
      tauIdFlags_[5] = (muTauPair.leg2()->userFloat("dRnearestMuon")      > 0.5);
      tauIdFlags_[6] = (TMath::Abs(muTauPair.leg1()->vertex().z() - muTauPair.leg2()->vertex().z()) < 0.2);
      double muTauPairChargeProd = muTauPair.leg1()->charge()*muTauPair.leg2()->userFloat("leadTrackCharge");
      tauIdFlags_[7] = (muTauPairChargeProd > muTauPairChargeProdMin_ && muTauPairChargeProd < muTauPairChargeProdMax_);
      for ( int iTauIdDiscriminator = 0; iTauIdDiscriminator < numTauIdDiscriminators_; ++iTauIdDiscriminator ) {
	const std::string tauIdDiscriminator = tauIdDiscriminators_[iTauIdDiscriminator];
	//std::cout << " tauIdDiscriminator = " << tauIdDiscriminator << ":" 
	//	    << " " << muTauPair.leg2()->tauID(tauIdDiscriminator.data()) << std::endl;
	tauIdFlags_[numPreselCuts_ + iTauIdDiscriminator] = (muTauPair.leg2()->tauID(tauIdDiscriminator.data()) > 0.5);
      }
      //std::cout << "tauIdFlags = " << format_vbool(tauIdFlags_) << std::endl;
      
      tauIdFlagsReversed_[0] = tauIdFlags_[0];
      for ( int iTauIdDiscriminator = 0; iTauIdDiscriminator < numTauIdDiscriminators_; ++iTauIdDiscriminator ) {
	tauIdFlagsReversed_[1 + iTauIdDiscriminator] = tauIdFlags_[numPreselCuts_ + iTauIdDiscriminator];
      }
      for ( int iPreselCut = 1; iPreselCut < numPreselCuts_; ++iPreselCut ) {
	tauIdFlagsReversed_[numTauIdDiscriminators_ + iPreselCut] = tauIdFlags_[iPreselCut];
      }
      //std::cout << "tauIdFlagsReversed = " << format_vbool(tauIdFlagsReversed_) << std::endl;

//--- fill histograms for "inclusive" tau id. efficiency measurement
      cutFlowUnbinned_->fillCutFlowTables(0., 
					  tauIdFlags_, tauIdFlagsReversed_, 
					  genMatchType, genTauCharge, recTauCharge,
					  evtWeight);

//--- fill histograms for tau id. efficiency measurement as function of:
//   o tau-jet transverse momentum
//   o tau-jet pseudo-rapidity
//   o reconstructed vertex multiplicity
//   o sumEt
//   o ...
      for ( std::vector<cutFlowEntryType*>::iterator cutFlowEntry = cutFlowEntriesBinned_.begin();
	    cutFlowEntry != cutFlowEntriesBinned_.end(); ++cutFlowEntry ) {
	double x = 0.;
	if      ( (*cutFlowEntry)->binVariable_ == "tauPt"       ) x = muTauPair.leg2()->pt();
	else if ( (*cutFlowEntry)->binVariable_ == "tauAbsEta"   ) x = TMath::Abs(muTauPair.leg2()->eta());
	else if ( (*cutFlowEntry)->binVariable_ == "numVertices" ) x = numVertices;
	else if ( (*cutFlowEntry)->binVariable_ == "sumEt"       ) x = muTauPair.met()->sumEt();
	else throw cms::Exception("regionEntryType::analyze")
	  << "Invalid binVariable = " << (*cutFlowEntry)->binVariable_ << " !!\n";
	(*cutFlowEntry)->fillCutFlowTables(x, 
					   tauIdFlags_, tauIdFlagsReversed_, 
					   genMatchType, genTauCharge, recTauCharge,
					   evtWeight);
      }

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
  vstring selectionNames_;
  vstring selectionNamesReversed_;
  
  int numPreselCuts_;
  int numTauIdDiscriminators_;

  TauIdEffEventSelector* selector_;
 
  double muTauPairChargeProdMin_;
  double muTauPairChargeProdMax_;

  vbool tauIdFlags_;
  vbool tauIdFlagsReversed_;
  
  cutFlowEntryType* cutFlowUnbinned_;
  std::vector<cutFlowEntryType*> cutFlowEntriesBinned_;

  int numMuTauPairs_selected_;
  double numMuTauPairsWeighted_selected_;
};

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauIdEffPreselNumbers>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauIdEffPreselNumbers");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteTauIdEffPreselNumbers") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauIdEffPreselNumbers = cfg.getParameter<edm::ParameterSet>("compTauIdEffPreselNumbers");

  edm::InputTag srcMuTauPairs = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcMuTauPairs");
  edm::InputTag srcGenParticles = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcGenParticles");
  edm::ParameterSet cfgEventSelCuts = cfgTauIdEffPreselNumbers.getParameter<edm::ParameterSet>("eventSelCuts");
  edm::InputTag srcHLTresults = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcHLTresults");
  edm::InputTag srcL1GtReadoutRecord = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcL1GtReadoutRecord"); 
  edm::InputTag srcL1GtObjectMapRecord = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcL1GtObjectMapRecord");
  vstring hltPaths = cfgTauIdEffPreselNumbers.getParameter<vstring>("hltPaths");
  edm::InputTag srcCaloMEt = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcCaloMEt");
  edm::InputTag srcGoodMuons = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcGoodMuons");
  edm::InputTag srcVertices = cfgTauIdEffPreselNumbers.getParameter<edm::InputTag>("srcVertices");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauIdEffPreselNumbers.getParameter<vInputTag>("weights");
  std::string sysShift = cfgTauIdEffPreselNumbers.exists("sysShift") ?
    cfgTauIdEffPreselNumbers.getParameter<std::string>("sysShift") : "CENTRAL_VALUE";

  std::string selEventsFileName = ( cfgTauIdEffPreselNumbers.exists("selEventsFileName") ) ?
    cfgTauIdEffPreselNumbers.getParameter<std::string>("selEventsFileName") : "";

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for different ABCD regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauIdEffPreselNumbers.getParameter<std::string>("process");
  vstring regions = cfgTauIdEffPreselNumbers.getParameter<vstring>("regions");
  edm::ParameterSet cfgBinning = cfgTauIdEffPreselNumbers.getParameter<edm::ParameterSet>("binning");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauIdEffPreselNumbers.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");
      regionEntryType* regionEntry = 
	new regionEntryType(fs, process, *region, tauIdDiscriminators, tauIdName, sysShift, cfgBinning, cfgEventSelCuts);
      regionEntries.push_back(regionEntry);
    }
  }

  edm::ParameterSet cfgSelectorABCD = cfgEventSelCuts;
  cfgSelectorABCD.addParameter<vstring>("tauIdDiscriminators", vstring());
  cfgSelectorABCD.addParameter<std::string>("region", "ABCD");
  cfgSelectorABCD.addParameter<std::string>("tauChargeMode", "tauLeadTrackCharge");
  cfgSelectorABCD.addParameter<bool>("disableTauCandPreselCuts", false);
  TauIdEffEventSelector* selectorABCD = new TauIdEffEventSelector(cfgSelectorABCD);

  std::ofstream* selEventsFile = ( selEventsFileName != "" ) ?
    new std::ofstream(selEventsFileName.data(), std::ios::out) : 0;

  int    numEvents_processed                     = 0; 
  double numEventsWeighted_processed             = 0.;
  int    numEvents_passedTrigger                 = 0;
  double numEventsWeighted_passedTrigger         = 0.;
  int    numEvents_passedDiMuonVeto              = 0;
  double numEventsWeighted_passedDiMuonVeto      = 0.;
  int    numEvents_passedDiMuTauPairVeto         = 0;
  double numEventsWeighted_passedDiMuTauPairVeto = 0.;
  
  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauIdEffPreselNumbers") 
	<< "Failed to open inputFile = " << (*inputFileName) << " !!\n";

    std::cout << "opening inputFile = " << (*inputFileName);
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("Events"));
    if ( tree ) std::cout << " (" << tree->GetEntries() << " Events)";
    std::cout << std::endl;

    fwlite::Event evt(inputFile);
    for ( evt.toBegin(); !(evt.atEnd() || maxEvents_processed); ++evt ) {

      //std::cout << "processing run = " << evt.id().run() << ":" 
      //	  << " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << std::endl;

//--- compute event weight
//   (pile-up reweighting, Data/MC correction factors,...)
      double evtWeight = 1.0;
      for ( vInputTag::const_iterator srcWeight = srcWeights.begin();
	    srcWeight != srcWeights.end(); ++srcWeight ) {
	edm::Handle<double> weight;
	evt.getByLabel(*srcWeight, weight);
	evtWeight *= (*weight);
      }

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

//--- check that event has passed triggers
      bool anyHLTpath_passed = false;
      if ( hltPaths.size() == 0 ) {
	anyHLTpath_passed = true;
      } else {
	edm::Handle<edm::TriggerResults> hltResults;
	evt.getByLabel(srcHLTresults, hltResults);
	
	const edm::TriggerNames& triggerNames = evt.triggerNames(*hltResults);
	
	for ( vstring::const_iterator hltPath = hltPaths.begin();
	      hltPath != hltPaths.end(); ++hltPath ) {
	  bool isHLTpath_passed = false;
	  unsigned int idx = triggerNames.triggerIndex(*hltPath);
	  if ( idx < triggerNames.size() ) isHLTpath_passed = hltResults->accept(idx);
	  if ( isHLTpath_passed ) anyHLTpath_passed = true;
	}
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

      typedef std::vector<pat::MET> PATMETCollection;
      edm::Handle<PATMETCollection> caloMEt;
      evt.getByLabel(srcCaloMEt, caloMEt);
      if ( caloMEt->size() != 1 )
	throw cms::Exception("FWLiteTauIdEffPreselNumbers")
	  << "Failed to find unique CaloMEt object !!\n";

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

//--- iterate over collection of muon + tau-jet pairs:
//    check which region muon + tau-jet pair is selected in
//    and whether reconstructed tau-jet matches "true" hadronic tau decay on generator level or is fake,
//    count number of "true" and fake taus selected in all regions
      edm::Handle<reco::GenParticleCollection> genParticles;
      evt.getByLabel(srcGenParticles, genParticles);

      int muTauPairIdx = 0;
      for ( PATMuTauPairCollection::const_iterator muTauPair = muTauPairs->begin();
	    muTauPair != muTauPairs->end(); ++muTauPair, ++muTauPairIdx ) {
	double genTauCharge, recTauCharge;
	int genMatchType = getGenMatchType(*muTauPair, *genParticles, &genTauCharge, &recTauCharge);
	for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	      regionEntry != regionEntries.end(); ++regionEntry ) {   
/*
          if ( (*regionEntry)->region_ == "C1" ) {
            pat::strbitset evtSelFlags;
	    if ( genMatchType == kGenTauHadMatched || (*regionEntry)->selector_->operator()(*muTauPair, evtSelFlags) ) {
	      std::cout << "muTauPair #" << muTauPairIdx << std::endl;
	      std::cout << " leg1: Pt = " << muTauPair->leg1()->pt() << "," 
	    	        << " eta = " << muTauPair->leg1()->eta() << ", phi = " << muTauPair->leg1()->phi() << std::endl;
	      std::cout << " leg2: Pt = " << muTauPair->leg2()->pt() << "," 
	      	        << " eta = " << muTauPair->leg2()->eta() << ", phi = " << muTauPair->leg2()->phi();
	      if      ( genMatchType == kTauHadMatched  ) std::cout << " ('true' hadronic tau decay)";
	      else if ( genMatchType == kFakeTauMatched ) std::cout << " (tau fake)"; 
	      else                                        std::cout << " (no gen. match)"; 
	      std::cout << std::endl;
	    }
          }
 */
	  (*regionEntry)->analyze(*muTauPair, 
				  genMatchType, genTauCharge, recTauCharge, 
				  caloMEt->front(),
				  numVertices, 
				  evtWeight);
	}

	if ( genMatchType == kGenTauHadMatched &&
	     muTauPair->leg2()->tauID("decayModeFinding") > 0.5 &&
	     muTauPair->leg2()->tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 &&
	     muTauPair->leg2()->userFloat("hasLeadTrack") > 0.5 &&
	     muTauPair->leg2()->userFloat("leadTrackPt") > 5.0 &&     
	     muTauPair->leg2()->userFloat("preselLoosePFIsoPt") > 2.5 ) {
	  //std::cout << "run = " << evt.id().run() << "," 
	  //	      << " ls = " << evt.luminosityBlock() << ", event = " << evt.id().event() << ":" << std::endl;
	  
	  if ( selEventsFile ) (*selEventsFile) << evt.id().run() << ":" << evt.luminosityBlock() << ":" << evt.id().event() << std::endl;
	}
      }
    }

//--- close input file
    delete inputFile;
  }

  std::cout << "<FWLiteTauIdEffPreselNumbers>:" << std::endl;
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
    std::cout << "  region " << (*regionEntry)->region_ << ":" << std::endl;
    std::cout << "   entries = " << (*regionEntry)->numMuTauPairs_selected_ 
	      << " (weighted = " << (*regionEntry)->numMuTauPairsWeighted_selected_ << ")" << std::endl;
    TauIdEffCutFlowTable* cutFlowTableTauHadMatched = (*regionEntry)->cutFlowUnbinned_->cutFlowTauHadMatched_;
    double effPreselection = cutFlowTableTauHadMatched->getCutFlowNumber(0, 7)/
                             cutFlowTableTauHadMatched->getCutFlowNumber(0, 0);
    std::cout << "  eff(preselection) = " << effPreselection << std::endl;
    double effTauId = cutFlowTableTauHadMatched->getCutFlowNumber(0, 7 + (*regionEntry)->tauIdDiscriminators_.size())/
                      cutFlowTableTauHadMatched->getCutFlowNumber(0, 7);
    std::cout << "  eff(tauId) = " << effTauId << std::endl;
    TauIdEffCutFlowTable* cutFlowTableFakeTauMatched = (*regionEntry)->cutFlowUnbinned_->cutFlowFakeTauMatched_;
    double purity = cutFlowTableTauHadMatched->getCutFlowNumber(0, 7)/
                   (cutFlowTableTauHadMatched->getCutFlowNumber(0, 7) + cutFlowTableFakeTauMatched->getCutFlowNumber(0, 7));
    std::cout << "  purity = " << purity << std::endl;
  }

  clock.Show("FWLiteTauIdEffPreselNumbers");

  return 0;
}
