
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/FWLite/interface/InputSource.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/TauIdEfficiency/bin/tauFakeRateAuxFunctions.h"

#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TString.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<makeTauFakeRatePlots>:" << std::endl;  

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeTauFakeRatePlots");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("makeTauFakeRatePlots") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgMakeTauFakeRatePlots = cfg.getParameter<edm::ParameterSet>("makeTauFakeRatePlots");

  typedef std::vector<std::string> vstring;
  vstring processesToPlot       = cfgMakeTauFakeRatePlots.getParameter<vstring>("processesToPlot");
  vstring tauIdsToPlot          = cfgMakeTauFakeRatePlots.getParameter<vstring>("tauIdsToPlot");
  vstring eventSelectionsToPlot = cfgMakeTauFakeRatePlots.getParameter<vstring>("eventSelectionsToPlot");

//--- read configuration parameters needed to make control plots 
//    comparing jetPt, jetEta, jetPhi distributions observed in analyzed dataset to Monte Carlo predictions
  vstring processNamesSim;
  std::string processNameData;
  std::map<std::string, histogramDrawOptionType> drawOptionsProcesses;
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgProcesses = cfgMakeTauFakeRatePlots.getParameter<vParameterSet>("processes");
  for ( vParameterSet::const_iterator cfgProcess = cfgProcesses.begin();
	cfgProcess != cfgProcesses.end(); ++cfgProcess ) {
    std::string processName = cfgProcess->getParameter<std::string>("name");
    std::string type = cfgProcess->getParameter<std::string>("type");
    if ( type == "Data" || type == "embeddedData" ) processNameData = processName;
    else processNamesSim.push_back(processName);
    drawOptionsProcesses[processName] = histogramDrawOptionType(*cfgProcess);
  }

//--- read configuration parameters needed to make plots
//    comparing fake-rates of different tau id. discriminators
//   (for one event selection at a time)
  std::map<std::string, graphDrawOptionType> drawOptionsTauIds; // key = tauId
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgTauIds = cfgMakeTauFakeRatePlots.getParameter<vParameterSet>("tauIds");
  for ( vParameterSet::const_iterator cfgTauId = cfgTauIds.begin();
	cfgTauId != cfgTauIds.end(); ++cfgTauId ) {
    std::string tauIdName = cfgTauId->getParameter<std::string>("name");
    drawOptionsTauIds[tauIdName] = graphDrawOptionType(*cfgTauId);
  }

//--- read configuration parameters needed to make plots
//    comparing fake-rates of different jets selected in QCD, W + jets and Zmumu events
//   (for one tau id. discriminator at a time)
  std::map<std::string, graphDrawOptionType> drawOptionsEventSelections; // key = eventSelection
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgEventSelections = cfgMakeTauFakeRatePlots.getParameter<vParameterSet>("eventSelections");
  for ( vParameterSet::const_iterator cfgEventSelection = cfgEventSelections.begin();
	cfgEventSelection != cfgEventSelections.end(); ++cfgEventSelection ) {
    std::string eventSelectionName = cfgEventSelection->getParameter<std::string>("name");
    drawOptionsEventSelections[eventSelectionName] = graphDrawOptionType(*cfgEventSelection);
  }

  vstring labels = cfgMakeTauFakeRatePlots.getParameter<vstring>("labels");

  vstring regions;
  regions.push_back("P"); // passed
  regions.push_back("F"); // failed

  vstring observables;
  observables.push_back("jetPt"); 
  observables.push_back("jetEta");
  observables.push_back("jetPhi");
  observables.push_back("sumEt");
  observables.push_back("numVertices");
      
  std::string outputFileName = cfgMakeTauFakeRatePlots.getParameter<std::string>("outputFileName");

  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauFakeRatePlots") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string inputFileName = (*inputFiles.files().begin());
  
//--- open input file
  TFile* inputFile = TFile::Open(inputFileName.data());
  if ( !inputFile ) 
    throw cms::Exception("makeTauFakeRatePlots") 
      << "Failed to open inputFile = " << inputFileName << " !!\n";

//--- load histograms
  histogramMapType5 histogramMap; // key = (eventSelection, process, tauId, observable, region)
  loadHistograms(histogramMap, inputFile, processesToPlot, eventSelectionsToPlot, tauIdsToPlot, observables, regions);

//--- compute histogram for sum of all Standard Model processes
//   (= Monte Carlo prediction to be compared with data)
  sumHistograms(histogramMap, processNamesSim, "sum");

//--- compute fake-rates
  frMapType4 frMap; // key = (eventSelection, process, tauId, observable)
  compFakeRates(histogramMap, frMap, "P", "F");

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.12);

//--- make control plots 
//    comparing jetPt, jetEta, jetPhi distributions observed in analyzed dataset to Monte Carlo predictions
  drawHistograms("jetPt", "P_{T}^{jet} / GeV", 
		 canvas, histogramMap, processNamesSim, processNameData, drawOptionsProcesses,
		 eventSelectionsToPlot, tauIdsToPlot, regions, 
		 labels, outputFileName);
  drawHistograms("jetEta", "#eta_{jet}", 
		 canvas, histogramMap, processNamesSim, processNameData, drawOptionsProcesses,
		 eventSelectionsToPlot, tauIdsToPlot, regions, 
		 labels, outputFileName);
  drawHistograms("jetPhi", "phi_{jet}", 
		 canvas, histogramMap, processNamesSim, processNameData, drawOptionsProcesses,
		 eventSelectionsToPlot, tauIdsToPlot, regions, 
		 labels, outputFileName);
  
//--- plot fake-rates of different tau id. discriminators
//   (for one event selection at a time)
  for ( vstring::const_iterator eventSelectionToPlot = eventSelectionsToPlot.begin();
	eventSelectionToPlot != eventSelectionsToPlot.end(); ++eventSelectionToPlot ) {
    frMapType3 frMapToPlot; // key = (process, tauId, observable)
    for ( vstring::const_iterator tauIdToPlot = tauIdsToPlot.begin();
	  tauIdToPlot != tauIdsToPlot.end(); ++tauIdToPlot ) {
      frMapToPlot["sum"][*tauIdToPlot] = frMap[*eventSelectionToPlot]["sum"][*tauIdToPlot];
      frMapToPlot[processNameData][*tauIdToPlot] = frMap[*eventSelectionToPlot][processNameData][*tauIdToPlot];
    }

    std::string label_eventSelection = drawOptionsEventSelections[*eventSelectionToPlot].legendEntry_;
    
    drawGraphs("jetPt", "P_{T}^{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, tauIdsToPlot, drawOptionsTauIds, 
	       *eventSelectionToPlot, labels, label_eventSelection, outputFileName);
    drawGraphs("jetEta", "#eta_{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, tauIdsToPlot, drawOptionsTauIds, 
	       *eventSelectionToPlot, labels, label_eventSelection, outputFileName);
    drawGraphs("jetPhi", "#phi_{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, tauIdsToPlot, drawOptionsTauIds, 
	       *eventSelectionToPlot, labels, label_eventSelection, outputFileName);
    
    drawGraphs("sumEt", "#Sigma E_{T} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, tauIdsToPlot, drawOptionsTauIds, 
	       *eventSelectionToPlot, labels, label_eventSelection, outputFileName);
    
    drawGraphs("numVertices", "Num. Vertices", 
	       canvas, frMapToPlot, "sum", processNameData, tauIdsToPlot, drawOptionsTauIds, 
	       *eventSelectionToPlot, labels, label_eventSelection, outputFileName);
  }

//--- plot fake-rates of QCD, W + jets and Zmumu events
//   (for one tau id. discriminator at a time)
  for ( vstring::const_iterator tauIdToPlot = tauIdsToPlot.begin();
	tauIdToPlot != tauIdsToPlot.end(); ++tauIdToPlot ) {
    frMapType3  frMapToPlot; // key = (process, tauId, observable)
    for ( vstring::const_iterator eventSelectionToPlot = eventSelectionsToPlot.begin();
	  eventSelectionToPlot != eventSelectionsToPlot.end(); ++eventSelectionToPlot ) {
      frMapToPlot["sum"][*tauIdToPlot] = frMap[*eventSelectionToPlot]["sum"][*tauIdToPlot];
      frMapToPlot[processNameData][*tauIdToPlot] = frMap[*eventSelectionToPlot][processNameData][*tauIdToPlot];
    }

    std::string label_tauId = drawOptionsTauIds[*tauIdToPlot].legendEntry_;
    
    drawGraphs("jetPt", "P_{T}^{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, eventSelectionsToPlot, drawOptionsEventSelections, 
	       *tauIdToPlot, labels, label_tauId, outputFileName);
    drawGraphs("jetEta", "#eta_{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, eventSelectionsToPlot, drawOptionsEventSelections, 
	       *tauIdToPlot, labels, label_tauId, outputFileName);
    drawGraphs("jetPhi", "#phi_{jet} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, eventSelectionsToPlot, drawOptionsEventSelections, 
	       *tauIdToPlot, labels, label_tauId, outputFileName);

    drawGraphs("sumEt", "#Sigma E_{T} / GeV", 
	       canvas, frMapToPlot, "sum", processNameData, eventSelectionsToPlot, drawOptionsEventSelections, 
	       *tauIdToPlot, labels, label_tauId, outputFileName);
    
    drawGraphs("numVertices", "Num. Vertices", 
	       canvas, frMapToPlot, "sum", processNameData, eventSelectionsToPlot, drawOptionsEventSelections, 
	       *tauIdToPlot, labels, label_tauId, outputFileName);
  }
  
  delete canvas;

//--print time that it took macro to run
  std::cout << "finished executing makeTauFakeRatePlots macro:" << std::endl;
  std::cout << " #tauIdDiscr.   = " << tauIdsToPlot.size() << std::endl;
  std::cout << " #evtSelections = " << eventSelectionsToPlot.size() << std::endl;
  clock.Show("makeTauFakeRatePlots");

  return 0;
}
