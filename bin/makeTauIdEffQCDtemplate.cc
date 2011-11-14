
/** \executable makeTauIdEffQCDtemplate
 *
 * Obtain QCD template from control region in Data,
 * corrected for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: makeTauIdEffQCDtemplate.cc,v 1.2 2011/11/06 13:25:23 veelken Exp $
 *
 */

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "TauAnalysis/TauIdEfficiency/bin/tauIdEffAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"

#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TBenchmark.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

typedef std::map<std::string, std::map<std::string, TH1*> > histogramMap;
typedef std::map<std::string, std::map<std::string, std::map<std::string, double> > > numEventsMap;

struct regionEntryType
{
  regionEntryType(const std::string& regionTakeQCDtemplateFromData, 
		  const std::string& regionWplusJetsSideband, 
		  const std::string& regionStoreQCDtemplate, 
		  const std::string& tauId, const std::string& fitVariable, const std::string& tauIdValue, 
		  const std::string& sysShift)
    : regionTakeQCDtemplateFromData_(regionTakeQCDtemplateFromData),
      regionWplusJetsSideband_(regionWplusJetsSideband),
      regionStoreQCDtemplate_(regionStoreQCDtemplate),
      tauId_(tauId),
      fitVariable_(fitVariable),
      tauIdValue_(tauIdValue),
      sysShift_(sysShift)
  {}
  ~regionEntryType() {}

  void makeQCDtemplate(TFileDirectory& dir, 
		       histogramMap& distributionsData, std::map<std::string, histogramMap>& templates, numEventsMap& numEvents)
  {
    //std::cout << "<regionEntryType::makeQCDtemplate>:" << std::endl;
    //std::cout << " fitVariable = " << fitVariable_ << std::endl;
    //std::cout << " tauId = " << tauId_ << std::endl;
    //std::cout << " tauIdValue = " << tauIdValue_ << std::endl;
    //std::cout << " regionWplusJetsSideband = " << regionWplusJetsSideband_ << std::endl;
    //std::cout << " regionTakeQCDtemplateFromData = " << regionTakeQCDtemplateFromData_ << std::endl;

    std::string keyWplusJetsSideband_data = getKey(fitVariable_, tauId_, "all");
    std::string keyWplusJetsSideband_mc = getKey(fitVariable_, tauId_, "all", sysShift_);
    TH1* distributionDataWplusJetsSideband = distributionsData[regionWplusJetsSideband_][keyWplusJetsSideband_data];
    //std::cout << "distributionDataWplusJetsSideband = " << distributionDataWplusJetsSideband << std::endl;
    double numEventsWplusJetsSideband_obsWplusJets = 
      getIntegral(distributionDataWplusJetsSideband, true, true)
     - (numEvents["Ztautau"][regionWplusJetsSideband_][keyWplusJetsSideband_mc]
      + numEvents["Zmumu"][regionWplusJetsSideband_][keyWplusJetsSideband_mc]
      + numEvents["QCD"][regionWplusJetsSideband_][keyWplusJetsSideband_mc]
      + numEvents["TTplusJets"][regionWplusJetsSideband_][keyWplusJetsSideband_mc]);
    double numEventsWplusJetsSideband_expWplusJets = 
      numEvents["WplusJets"][regionWplusJetsSideband_][keyWplusJetsSideband_mc];
    std::cout << "numEventsWplusJets_sideband: observed = " << numEventsWplusJetsSideband_obsWplusJets << "," 
	      << " expected = " << numEventsWplusJetsSideband_expWplusJets << std::endl;
    
    double scaleFactorWplusJets = numEventsWplusJetsSideband_obsWplusJets/numEventsWplusJetsSideband_expWplusJets;

    std::string keyQCDsideband_data = getKey(fitVariable_, tauId_, "all");
    std::string keyQCDsideband_mc = getKey(fitVariable_, tauId_, "all", sysShift_);    
    TH1* distributionDataQCDsideband = distributionsData[regionTakeQCDtemplateFromData_][keyQCDsideband_data];
    //std::cout << "distributionDataQCDsideband = " << distributionDataQCDsideband << std::endl;
    TH1* templateZtautauQCDsideband = templates["Ztautau"][regionTakeQCDtemplateFromData_][keyQCDsideband_mc];
    TH1* templateZmumuQCDsideband = templates["Zmumu"][regionTakeQCDtemplateFromData_][keyQCDsideband_mc];
    TH1* templateWplusJetsQCDsideband = templates["WplusJets"][regionTakeQCDtemplateFromData_][keyQCDsideband_mc];
    TH1* templateTTplusJetsQCDsideband = templates["TTplusJets"][regionTakeQCDtemplateFromData_][keyQCDsideband_mc];

    TString templateQCDsidebandName = distributionDataQCDsideband->GetName();    
    templateQCDsidebandName.ReplaceAll(regionTakeQCDtemplateFromData_, regionStoreQCDtemplate_);
    if ( sysShift_ != "CENTRAL_VALUE" ) templateQCDsidebandName.Append("_").Append(sysShift_);
    TH1* templateQCDsideband_obsQCD = 0;
    TAxis* xAxis = distributionDataQCDsideband->GetXaxis(); 
    if ( xAxis->GetXbins() && xAxis->GetXbins()->GetSize() >= 2 ) {
      const TArrayD* binning = xAxis->GetXbins();      
      templateQCDsideband_obsQCD = 
	dir.make<TH1D>(templateQCDsidebandName, 
		       distributionDataQCDsideband->GetTitle(), binning->GetSize() - 1, binning->GetArray());
    } else {
      templateQCDsideband_obsQCD = 
	dir.make<TH1D>(templateQCDsidebandName, 
		       distributionDataQCDsideband->GetTitle(), xAxis->GetNbins(), xAxis->GetXmin(), xAxis->GetXmax());
    }
    templateQCDsideband_obsQCD->Add(distributionDataQCDsideband, +1.);
    templateQCDsideband_obsQCD->Add(templateZtautauQCDsideband, -1.);
    templateQCDsideband_obsQCD->Add(templateZmumuQCDsideband, -1.);
    templateQCDsideband_obsQCD->Add(templateWplusJetsQCDsideband, -1.*scaleFactorWplusJets);
    templateQCDsideband_obsQCD->Add(templateTTplusJetsQCDsideband, -1.);
  }

  std::string regionTakeQCDtemplateFromData_;
  std::string regionWplusJetsSideband_;
  std::string regionStoreQCDtemplate_;

  std::string tauId_;
  std::string fitVariable_;
  std::string tauIdValue_;
  
  std::string sysShift_;
};

int main(int argc, const char* argv[])
{
//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  std::cout << "<makeTauIdEffQCDtemplate>:" << std::endl;

//--- disable pop-up windows showing graphics output
  gROOT->SetBatch(true);

//--- load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeTauIdEffQCDtemplate");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfgMakeTauIdEffQCDtemplate = cfg.getParameter<edm::ParameterSet>("makeTauIdEffQCDtemplate");

  typedef std::vector<std::string> vstring;
  vstring tauIds = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("tauIds");
  vstring tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");

  vstring fitVariables = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("fitVariables");
  
  std::vector<std::string> loadSysShifts;
  loadSysShifts.push_back("CENTRAL_VALUE");
  vstring loadSysUncertaintyNames = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("loadSysUncertainties");
  for ( vstring::const_iterator sysUncertaintyName = loadSysUncertaintyNames.begin();
	sysUncertaintyName != loadSysUncertaintyNames.end(); ++sysUncertaintyName ) {
    if ( (*sysUncertaintyName) == "sysAddPUsmearing" ) {
      loadSysShifts.push_back(*sysUncertaintyName);
    } else {
      loadSysShifts.push_back(std::string(*sysUncertaintyName).append("Up"));
      loadSysShifts.push_back(std::string(*sysUncertaintyName).append("Down"));
    } 
  }

  vstring regions;
  std::vector<regionEntryType*> regionEntries;
  for ( vstring::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {
    for ( vstring::const_iterator tauIdValue = tauIdValues.begin();
	  tauIdValue != tauIdValues.end(); ++tauIdValue ) {
      for ( vstring::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	for ( std::vector<std::string>::const_iterator sysShift = loadSysShifts.begin();
	      sysShift != loadSysShifts.end(); ++sysShift ) {
	  std::string regionTakeQCDtemplateFromData = 
	    cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionTakeQCDtemplateFromData_").append(*tauIdValue));
	  regions.push_back(regionTakeQCDtemplateFromData);
	  regions.push_back(std::string(regionTakeQCDtemplateFromData).append("+"));
	  regions.push_back(std::string(regionTakeQCDtemplateFromData).append("-"));
	  std::string regionWplusJetsSideband = 
	    cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionWplusJetsSideband_").append(*tauIdValue));
	  regions.push_back(regionWplusJetsSideband);
	  regions.push_back(std::string(regionWplusJetsSideband).append("+"));
	  regions.push_back(std::string(regionWplusJetsSideband).append("-"));
	  std::string regionStoreQCDtemplate = 
	    cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionStoreQCDtemplate_").append(*tauIdValue));

	  // all tau charges
	  regionEntryType* regionEntry = 
	    new regionEntryType(regionTakeQCDtemplateFromData, 
				regionWplusJetsSideband, 
				regionStoreQCDtemplate, 
				*tauId, *fitVariable, *tauIdValue, *sysShift);
	  regionEntries.push_back(regionEntry);
	  
	  // tau+ candidates only
	  regionEntryType* regionEntry_plus = 
	    new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("+"), 
				std::string(regionWplusJetsSideband).append("+"), 
				std::string(regionStoreQCDtemplate).append("+"), 
				*tauId, *fitVariable, *tauIdValue, *sysShift);
	  regionEntries.push_back(regionEntry_plus);
	  
	  // tau- candidates only
	  regionEntryType* regionEntry_minus = 
	    new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("-"), 
				std::string(regionWplusJetsSideband).append("-"), 
				std::string(regionStoreQCDtemplate).append("-"), 
				*tauId, *fitVariable, *tauIdValue, *sysShift);
	  regionEntries.push_back(regionEntry_minus);
	}
      }
    }
  }

  fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  histogramMap distributionsData; // key = (region, observable + sysShift)

  histogramMap templatesZtautau;
  histogramMap templatesZmumu;
  histogramMap templatesQCD;
  histogramMap templatesWplusJets;
  histogramMap templatesTTplusJets;

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgMakeTauIdEffQCDtemplate.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";
  
   for ( std::vector<std::string>::const_iterator sysShift = loadSysShifts.begin();
	sysShift != loadSysShifts.end(); ++sysShift ) {
    std::cout << "loading histograms for sysShift = " << (*sysShift) << "..." << std::endl;
    
    if ( (*sysShift) == "CENTRAL_VALUE" ) {
      loadHistograms(distributionsData, histogramInputDirectory, "Data",       regions, tauIds, fitVariables, *sysShift);
    }

    loadHistograms(templatesZtautau,    histogramInputDirectory, "Ztautau",    regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesZmumu,      histogramInputDirectory, "Zmumu",      regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesQCD,        histogramInputDirectory, "QCD",        regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesWplusJets,  histogramInputDirectory, "WplusJets",  regions, tauIds, fitVariables, *sysShift);
    loadHistograms(templatesTTplusJets, histogramInputDirectory, "TTplusJets", regions, tauIds, fitVariables, *sysShift);
  }

  std::map<std::string, histogramMap> templatesAll; // key = (process, region, observable)
  templatesAll["Ztautau"]    = templatesZtautau;
  templatesAll["Zmumu"]      = templatesZmumu;
  templatesAll["QCD"]        = templatesQCD;
  templatesAll["WplusJets"]  = templatesWplusJets;
  templatesAll["TTplusJets"] = templatesTTplusJets;

  std::vector<std::string> processes;
  processes.push_back(std::string("Ztautau"));
  processes.push_back(std::string("Zmumu"));
  processes.push_back(std::string("QCD"));
  processes.push_back(std::string("WplusJets"));
  processes.push_back(std::string("TTplusJets"));

  numEventsMap numEventsAll = // key = (process/"sum", region, observable)
    compNumEvents(templatesAll, processes, regions, distributionsData);
  
  TFileDirectory histogramOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;

  for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	regionEntry != regionEntries.end(); ++regionEntry ) {
    (*regionEntry)->makeQCDtemplate(histogramOutputDirectory, distributionsData, templatesAll, numEventsAll);
  }

  delete histogramInputFile;

//--print time that it took macro to run
  std::cout << "finished executing makeTauIdEffQCDtemplate macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("makeTauIdEffQCDtemplate");

  return 0;
}
