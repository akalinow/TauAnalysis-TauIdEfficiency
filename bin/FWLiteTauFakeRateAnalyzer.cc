
/** \executable FWLiteTauFakeRateAnalyzer
 *
 * Apply event selections for 
 *  o QCD multi-jet
 *  o QCD muon enriched
 *  o W + jets enriched
 *  o Zmumu enriched
 * event samples and fill histograms for jet --> tau fake-rate measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.6 $
 *
 * $Id: FWLiteTauFakeRateAnalyzer.cc,v 1.6 2012/02/02 09:03:32 veelken Exp $
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

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/TauIdEfficiency/interface/TauFakeRateEventSelector.h"
#include "TauAnalysis/TauIdEfficiency/interface/TauFakeRateHistManager.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>

typedef std::vector<std::string> vstring;
typedef StringCutObjectSelector<pat::Tau> StringCutPatTauSelector;

struct regionEntryType
{
  regionEntryType(TFileDirectory& dir,
		  const std::string& process, const std::string& region, 
		  const vstring& tauIdDiscriminators, const std::string& tauIdName)
    : process_(process),
      region_(region),
      tauIdDiscriminators_(tauIdDiscriminators),
      tauIdName_(tauIdName),
      selector_(0),
      histograms_(0),
      numTauJetCands_processed_(0),
      numTauJetCandsWeighted_processed_(0.),
      numTauJetCands_selected_(0),
      numTauJetCandsWeighted_selected_(0.)
  {
    //std::cout << "<regionEntryType>:" << std::endl;
    //std::cout << " region = " << region << std::endl;

    edm::ParameterSet cfgSelector;
    cfgSelector.addParameter<vstring>("tauIdDiscriminators", tauIdDiscriminators_);
    cfgSelector.addParameter<std::string>("region", region_);

    selector_ = new TauFakeRateEventSelector(cfgSelector);

    edm::ParameterSet cfgHistManager;
    cfgHistManager.addParameter<std::string>("process", process_);
    cfgHistManager.addParameter<std::string>("region", region_);
    cfgHistManager.addParameter<std::string>("tauIdDiscriminator", tauIdName_);

    histograms_ = new TauFakeRateHistManager(cfgHistManager);
    histograms_->bookHistograms(dir);
  }
  ~regionEntryType()
  {
    delete selector_;

    delete histograms_;
  }
  void analyze(const pat::Tau& tauJetCand, size_t numVertices, double sumEt, double evtWeight)
  {
    ++numTauJetCands_processed_;
    numTauJetCandsWeighted_processed_ += evtWeight;

    pat::strbitset evtSelFlags;
    if ( selector_->operator()(tauJetCand, evtSelFlags) ) {
      histograms_->fillHistograms(tauJetCand, numVertices, sumEt, evtWeight);

      ++numTauJetCands_selected_;
      numTauJetCandsWeighted_selected_ += evtWeight;
    }
  }

  std::string process_;
  std::string region_;
  vstring tauIdDiscriminators_;
  std::string tauIdName_;

  TauFakeRateEventSelector* selector_;

  TauFakeRateHistManager* histograms_;

  double numTauJetCands_processed_;
  double numTauJetCandsWeighted_processed_;
  double numTauJetCands_selected_;
  double numTauJetCandsWeighted_selected_;
};

int main(int argc, char* argv[]) 
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<FWLiteTauFakeRateAnalyzer>:" << std::endl;

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("FWLiteTauFakeRateAnalyzer");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("FWLiteTauFakeRateAnalyzer") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgTauFakeRateAnalyzer = cfg.getParameter<edm::ParameterSet>("tauFakeRateAnalyzer");

  edm::InputTag srcTauJetCandidates = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcTauJetCandidates");
  edm::InputTag srcVertices = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcVertices");
  edm::InputTag srcMET = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcMET");
  edm::InputTag srcTrigger = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcTrigger");
  vstring hltPaths = cfgTauFakeRateAnalyzer.getParameter<vstring>("hltPaths");
  typedef std::vector<edm::InputTag> vInputTag;
  vInputTag srcWeights = cfgTauFakeRateAnalyzer.getParameter<vInputTag>("weights");
  edm::InputTag srcEventCounter = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcEventCounter");

  fwlite::InputSource inputFiles(cfg); 
  int maxEvents = inputFiles.maxEvents();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

//--- initialize selections and histograms
//    for P(assed)/F(ailed) and A(ll) regions
  std::vector<regionEntryType*> regionEntries;

  std::string process = cfgTauFakeRateAnalyzer.getParameter<std::string>("process");
  std::cout << " process = " << process << std::endl;
  std::string processType = cfgTauFakeRateAnalyzer.getParameter<std::string>("type");
  std::cout << " type = " << processType << std::endl;
  bool isData = (processType == "Data");
  std::string evtSel = cfgTauFakeRateAnalyzer.getParameter<std::string>("evtSel");
  TFileDirectory dir = fs.mkdir(evtSel);
  vstring regions = cfgTauFakeRateAnalyzer.getParameter<vstring>("regions");
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIdDiscriminators = cfgTauFakeRateAnalyzer.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauIdDiscriminator = cfgTauIdDiscriminators.begin();
	cfgTauIdDiscriminator != cfgTauIdDiscriminators.end(); ++cfgTauIdDiscriminator ) {
    for ( vstring::const_iterator region = regions.begin();
	  region != regions.end(); ++region ) {
      vstring tauIdDiscriminators = cfgTauIdDiscriminator->getParameter<vstring>("discriminators");
      std::string tauIdName = cfgTauIdDiscriminator->getParameter<std::string>("name");
      regionEntryType* regionEntry = new regionEntryType(dir, process, *region, tauIdDiscriminators, tauIdName);
      regionEntries.push_back(regionEntry);
    }
  }
  std::vector<StringCutPatTauSelector*> tauJetCandSelection;
  vstring tauJetCandSelection_string = cfgTauFakeRateAnalyzer.getParameter<vstring>("tauJetCandSelection");
  for ( vstring::const_iterator tauJetCandSelCriterion = tauJetCandSelection_string.begin();
	tauJetCandSelCriterion != tauJetCandSelection_string.end(); ++tauJetCandSelCriterion ) {
    tauJetCandSelection.push_back(new StringCutPatTauSelector(*tauJetCandSelCriterion));
  }

//--- book "dummy" histogram counting number of processed events
  TH1* histogramEventCounter = fs.make<TH1F>("numEventsProcessed", "Number of processed Events", 3, -0.5, +2.5);
  histogramEventCounter->GetXaxis()->SetBinLabel(1, "all Events (DBS)");      // CV: bin numbers start at 1 (not 0) !!
  histogramEventCounter->GetXaxis()->SetBinLabel(2, "processed by Skimming");
  histogramEventCounter->GetXaxis()->SetBinLabel(3, "analyzed in PAT-tuple");
  
  int allEvents_DBS = cfgTauFakeRateAnalyzer.getParameter<int>("allEvents_DBS");
  if ( allEvents_DBS > 0 ) {
    histogramEventCounter->SetBinContent(1, cfgTauFakeRateAnalyzer.getParameter<int>("allEvents_DBS"));
  } else {
    histogramEventCounter->SetBinContent(1, -1.);
  }
  
  double xSection = cfgTauFakeRateAnalyzer.getParameter<double>("xSection");
  double intLumiData = cfgTauFakeRateAnalyzer.getParameter<double>("intLumiData");

  int    numEvents_processed             = 0; 
  double numEventsWeighted_processed     = 0.;
  int    numEvents_passedTrigger         = 0;
  double numEventsWeighted_passedTrigger = 0.;

  edm::RunNumber_t lastLumiBlock_run = -1;
  edm::LuminosityBlockNumber_t lastLumiBlock_ls = -1;

  double intLumiData_analyzed = 0.;
  edm::InputTag srcLumiProducer = cfgTauFakeRateAnalyzer.getParameter<edm::InputTag>("srcLumiProducer");

  bool maxEvents_processed = false;
  for ( vstring::const_iterator inputFileName = inputFiles.files().begin();
	inputFileName != inputFiles.files().end() && !maxEvents_processed; ++inputFileName ) {

//--- open input file
    TFile* inputFile = TFile::Open(inputFileName->data());
    if ( !inputFile ) 
      throw cms::Exception("FWLiteTauFakeRateAnalyzer") 
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

//--- quit event loop if maximal number of events to be processed is reached 
      ++numEvents_processed;
      numEventsWeighted_processed += evtWeight;
      if ( maxEvents > 0 && numEvents_processed >= maxEvents ) maxEvents_processed = true;

//--- check that event has passed triggers
//
//    CV: In order to match pile-up distribution in Monte Carlo with Data,
//        weight Data events by probability to pass trigger prescales.
//        Note that this assumes that the prescales of all HLT paths are uncorrelated
//       (assumption is not valid in case HLT paths share L1 conditions and those L1 conditions are prescaled)
//     
      edm::Handle<pat::TriggerEvent> hltEvent;
      evt.getByLabel(srcTrigger, hltEvent);
  
      bool isTriggered = false;
      double probFailedPrescale = 1.;
      for ( vstring::const_iterator hltPathName = hltPaths.begin();
	    hltPathName != hltPaths.end() && !isTriggered; ++hltPathName ) {
	//std::cout << "hltPathName = " << (*hltPathName) << std::endl;
	if ( (*hltPathName) == "*" ) { // check for wildcard character "*" that accepts all events
	  isTriggered = true;
	  probFailedPrescale = 0.;
	  break;
	} else {
	  const pat::TriggerPath* hltPath = hltEvent->path(*hltPathName);
	  if ( hltPath && hltPath->wasAccept() ) {
	    isTriggered = true;
	    unsigned hltPrescale = hltPath->prescale();
	    if ( hltPrescale < 1 ) hltPrescale = 1;
	    //std::cout << "HLT path = " << hltPath->name() << ": prescale = " << hltPrescale << std::endl;
	    double probFailedL1Prescale = 1.;
	    const pat::L1SeedCollection& l1Seeds = hltPath->l1Seeds();
	    for ( pat::L1SeedCollection::const_iterator l1Seed_status = l1Seeds.begin();
		  l1Seed_status != l1Seeds.end(); ++l1Seed_status ) {
	      const std::string& l1SeedName = l1Seed_status->second;
	      //std::cout << "l1SeedName = " << l1SeedName << std::endl;	      
	      const pat::TriggerAlgorithm* l1Seed = hltEvent->algorithm(l1SeedName);
	      if ( !l1Seed ) {
		//std::cout << "Failed to access L1 seed = " << l1SeedName << "," 
		//	    << " needed for HLT path = " << hltPath->name() << " !!" << std::endl;
		vstring l1SeedNames;
		const pat::TriggerAlgorithmCollection* l1Seeds = hltEvent->algorithms();
		if ( l1Seeds ) {
		  for ( pat::TriggerAlgorithmCollection::const_iterator l1Seed = l1Seeds->begin();
			l1Seed != l1Seeds->end(); ++l1Seed ) {
		    l1SeedNames.push_back(l1Seed->name());
		  }
		  //std::cout << "Available L1 seeds = " << format_vstring(l1SeedNames) << std::endl;
		} else {
		  //std::cout << "No L1 seeds available in pat::TriggerEvent !!" << std::endl;
		}
		continue;
	      }
	      bool l1Passed = l1Seed->gtlResult();
	      if ( l1Passed ) {
		unsigned l1Prescale = l1Seed->prescale();
		if ( l1Prescale < 1 ) l1Prescale = 1;
		//std::cout << " L1 seed = " << l1Seed->name() << ": prescale = " << l1Prescale << std::endl;
		if ( l1Prescale <= 1 ) probFailedL1Prescale = 0.;
		else probFailedL1Prescale *= (1. - 1./l1Prescale);
	      }
	    }
	    probFailedPrescale *= (1. - (1./hltPrescale)*(1. - probFailedL1Prescale));
	  }
	}
      }
      
      if ( !isTriggered ) continue;
      ++numEvents_passedTrigger;
      numEventsWeighted_passedTrigger += evtWeight;

      //std::cout << "probFailedPrescale = " << probFailedPrescale << std::endl;
      if ( isData && probFailedPrescale > (1. - 1.e-9) ) continue;

      double prescaleCorrFactor = ( isData ) ? 1./(1. - probFailedPrescale) : 1.;
      //std::cout << "prescaleCorrFactor = " << prescaleCorrFactor << std::endl;
      evtWeight *= prescaleCorrFactor;

//--- determine number of vertices reconstructed in the event
//   (needed to parametrize dependency of jet --> tau fake-rate on number of pile-up interactions)
      edm::Handle<reco::VertexCollection> vertices;
      evt.getByLabel(srcVertices, vertices);
      size_t numVertices = vertices->size();

//--- determine sumEt of all particles in the event
//   (needed to parametrize dependency of jet --> tau fake-rate on level of hadronic activity)
      edm::Handle<pat::METCollection> METs;
      evt.getByLabel(srcMET, METs);
      if ( !(METs->size() == 1) )
	throw cms::Exception("FWLiteTauFakeRateAnalyzer") 
	  << "Failed to find unique MET object in the event !!\n";
      double sumEt = METs->begin()->sumEt();

//--- iterate over collection of tau-jet candidates:
//    check if tau-jet candidate passed/fails tau id. criteria,
//    fill corresponding histograms
      edm::Handle<pat::TauCollection> tauJetCandidates;
      evt.getByLabel(srcTauJetCandidates, tauJetCandidates);
      for ( pat::TauCollection::const_iterator tauJetCand = tauJetCandidates->begin();
	    tauJetCand != tauJetCandidates->end(); ++tauJetCand ) {
	bool passesTauJetCandSelection = true;
	for ( std::vector<StringCutPatTauSelector*>::const_iterator tauJetCandSelCriterion = tauJetCandSelection.begin();
	      tauJetCandSelCriterion != tauJetCandSelection.end(); ++tauJetCandSelCriterion ) {
	  if ( !(**tauJetCandSelCriterion)(*tauJetCand) ) {
	    passesTauJetCandSelection = false;
	    break;
	  }
	}

	if ( !passesTauJetCandSelection ) continue;

	for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	      regionEntry != regionEntries.end(); ++regionEntry ) {
	  (*regionEntry)->analyze(*tauJetCand, numVertices, sumEt, evtWeight);
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

    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {  
      (*regionEntry)->histograms_->scaleHistograms(mcScaleFactor);
    }
  }

  std::cout << "<FWLiteTauFakeRateAnalyzer>:" << std::endl;
  std::cout << " numEvents_processed: " << numEvents_processed 
	    << " (weighted = " << numEventsWeighted_processed << ")" << std::endl;
  std::cout << " numEvents_passedTrigger: " << numEvents_passedTrigger 
	    << " (weighted = " << numEventsWeighted_passedTrigger << ")" << std::endl;
  std::string lastTauIdName = "";
  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    if ( (*regionEntry)->tauIdName_ != lastTauIdName ) 
      std::cout << " numTauJetCands_processed/numTauJetCands_selected," 
		<< " " << (*regionEntry)->tauIdName_ << std::endl;
    std::cout << "  region " << (*regionEntry)->region_ << ":" 
	      << " " << (*regionEntry)->numTauJetCands_processed_ 
	      << "/" << (*regionEntry)->numTauJetCands_selected_ 
	      << " (weighted = " << (*regionEntry)->numTauJetCandsWeighted_processed_ 
	      << "/" << (*regionEntry)->numTauJetCandsWeighted_selected_ << ")" << std::endl;
    lastTauIdName = (*regionEntry)->tauIdName_;
  }
  
  if ( isData ) {
    std::cout << " intLumiData (recorded, Trigger prescale corr.) = " << intLumiData << " pb" << std::endl;
    // CV: luminosity is recorded in some 'weird' units,
    //     needs to be multiplied by factor 0.10 in order to be in units of pb^-1
    std::cout << " intLumiData_analyzed (recorded) = " << intLumiData_analyzed*1.e-6*0.10 << " pb" << std::endl;
  }

  clock.Show("FWLiteTauFakeRateAnalyzer");

  return 0;
}
