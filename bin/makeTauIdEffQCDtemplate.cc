
/** \executable makeTauIdEffQCDtemplate
 *
 * Obtain QCD template from control region in Data,
 * corrected for contributions of Ztautau signal plus Zmumu, W + jets and TTbar backgrounds
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: makeTauIdEffQCDtemplate.cc,v 1.4 2011/11/26 15:37:59 veelken Exp $
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

#include <TROOT.h>
#include <TSystem.h>
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
		  const std::string& tauId, const std::string& fitVariable, const std::string& sysUncertainty)
    : regionTakeQCDtemplateFromData_(regionTakeQCDtemplateFromData),
      regionWplusJetsSideband_(regionWplusJetsSideband),
      regionStoreQCDtemplate_(regionStoreQCDtemplate),
      tauId_(tauId),
      fitVariable_(fitVariable),
      sysUncertainty_(sysUncertainty)
  {}
  ~regionEntryType() {}

  void makeQCDtemplate(TFileDirectory& dir, 
		       histogramMap3& histograms_data, histogramMap4& histograms_mc, valueMap3& numEvents_mc)
  {
    //std::cout << "<regionEntryType::makeQCDtemplate>" << std::endl;
    //std::cout << " regionWplusJetsSideband = " << regionWplusJetsSideband_ << std::endl;
    //std::cout << " tauId = " << tauId_ << std::endl;
    //std::cout << " fitVariable = " << fitVariable_ << std::endl;
    //std::cout << " sysUncertainty = " << sysUncertainty_ << std::endl;

    double numEventsWplusJetsSideband_data = 
      getIntegral(histograms_data[regionWplusJetsSideband_]["EventCounter"][key_central_value], true, true);
    double numEventsWplusJetsSideband_expBgr = 
      (numEvents_mc["Ztautau"][regionWplusJetsSideband_][sysUncertainty_]
     + numEvents_mc["Zmumu"][regionWplusJetsSideband_][sysUncertainty_]
     + numEvents_mc["QCD"][regionWplusJetsSideband_][sysUncertainty_]
     + numEvents_mc["TTplusJets"][regionWplusJetsSideband_][sysUncertainty_]);
    double numEventsWplusJetsSideband_obsWplusJets = numEventsWplusJetsSideband_data - numEventsWplusJetsSideband_expBgr;
    double numEventsWplusJetsSideband_expWplusJets = 
      numEvents_mc["WplusJets"][regionWplusJetsSideband_][sysUncertainty_];    
    //std::cout << "WplusJets sideband: observed = " << numEventsWplusJetsSideband_data << ","
    //	        << " expected background = " << numEventsWplusJetsSideband_expBgr 
    //	        << " --> observed WplusJets contribution = " << numEventsWplusJetsSideband_obsWplusJets << ","
    //	        << " expected = " << numEventsWplusJetsSideband_expWplusJets << std::endl;
    
    double scaleFactorWplusJets = numEventsWplusJetsSideband_obsWplusJets/numEventsWplusJetsSideband_expWplusJets;

    TH1* distributionDataQCDsideband = histograms_data[regionTakeQCDtemplateFromData_][fitVariable_][key_central_value];
    TH1* templateZtautauQCDsideband = histograms_mc["Ztautau"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateZmumuQCDsideband = histograms_mc["Zmumu"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateWplusJetsQCDsideband = histograms_mc["WplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];
    TH1* templateTTplusJetsQCDsideband = histograms_mc["TTplusJets"][regionTakeQCDtemplateFromData_][fitVariable_][sysUncertainty_];

    TString templateQCDsidebandName = distributionDataQCDsideband->GetName();    
    templateQCDsidebandName.ReplaceAll(regionTakeQCDtemplateFromData_, regionStoreQCDtemplate_);
    if ( sysUncertainty_ != key_central_value ) templateQCDsidebandName.Append("_").Append(sysUncertainty_);
    std::cout << "creating histogram = " << templateQCDsidebandName << std::endl;
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
  
  std::string sysUncertainty_;
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

  vstring tauIds = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("tauIds");
  vstring tauIdValues;
  tauIdValues.push_back("passed");
  tauIdValues.push_back("failed");

  vstring fitVariables = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("fitVariables");
  add_string_uniquely(fitVariables, "EventCounter"); // CV: for normalization purposes, always add 'EventCounter'

  vstring sysUncertainties = cfgMakeTauIdEffQCDtemplate.getParameter<vstring>("sysUncertainties");
  vstring sysUncertainties_expanded;
  sysUncertainties_expanded.push_back(key_central_value);
  for ( vstring::const_iterator sysUncertainty = sysUncertainties.begin();
	sysUncertainty != sysUncertainties.end(); ++sysUncertainty ) {
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Up"));
    sysUncertainties_expanded.push_back(std::string(*sysUncertainty).append("Down"));
  }
  
  vstring sysUncertainties_data;
  sysUncertainties_data.push_back(key_central_value);

   fwlite::InputSource inputFiles(cfg); 
  if ( inputFiles.files().size() != 1 ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Input file must be unique, got = " << format_vstring(inputFiles.files()) << " !!\n";
  std::string histogramFileName = (*inputFiles.files().begin());

  TFile* histogramInputFile = new TFile(histogramFileName.data());
  std::string directory = cfgMakeTauIdEffQCDtemplate.getParameter<std::string>("directory");
  TDirectory* histogramInputDirectory = ( directory != "" ) ?
    dynamic_cast<TDirectory*>(histogramInputFile->Get(directory.data())) : histogramInputFile;
  if ( !histogramInputDirectory ) 
    throw cms::Exception("makeTauIdEffQCDtemplate") 
      << "Directory = " << directory << " does not exists in input file = " << histogramFileName << " !!\n";

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TFileDirectory histogramOutputDirectory = ( directory != "" ) ?
    fs.mkdir(directory.data()) : fs;

  vstring processes;
  processes.push_back(std::string("Ztautau"));
  processes.push_back(std::string("Zmumu"));
  processes.push_back(std::string("QCD"));
  processes.push_back(std::string("WplusJets"));
  processes.push_back(std::string("TTplusJets"));

  for ( vstring::const_iterator tauId = tauIds.begin();
	tauId != tauIds.end(); ++tauId ) {

    vstring regions;
    std::vector<regionEntryType*> regionEntries;

    for ( vstring::const_iterator tauIdValue = tauIdValues.begin();
	  tauIdValue != tauIdValues.end(); ++tauIdValue ) {

      std::string regionTakeQCDtemplateFromData = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionTakeQCDtemplateFromData_").append(*tauIdValue));
      add_string_uniquely(regions, regionTakeQCDtemplateFromData);
      add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("+"));
      add_string_uniquely(regions, std::string(regionTakeQCDtemplateFromData).append("-"));	  
      std::string regionWplusJetsSideband = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionWplusJetsSideband_").append(*tauIdValue));
      add_string_uniquely(regions, regionWplusJetsSideband);
      add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("+"));
      add_string_uniquely(regions, std::string(regionWplusJetsSideband).append("-"));
      std::string regionStoreQCDtemplate = 
	cfgMakeTauIdEffQCDtemplate.getParameter<std::string>(std::string("regionStoreQCDtemplate_").append(*tauIdValue));
      
      for ( vstring::const_iterator fitVariable = fitVariables.begin();
	    fitVariable != fitVariables.end(); ++fitVariable ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin();
	      sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  // all tau charges
	  regionEntryType* regionEntry = 
	    new regionEntryType(regionTakeQCDtemplateFromData, 
				regionWplusJetsSideband, 
				regionStoreQCDtemplate, 
				*tauId, *fitVariable, *sysUncertainty);
	  regionEntries.push_back(regionEntry);
	  
	  // tau+ candidates only
	  regionEntryType* regionEntry_plus = 
	    new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("+"), 
				std::string(regionWplusJetsSideband).append("+"), 
				std::string(regionStoreQCDtemplate).append("+"), 
				*tauId, *fitVariable, *sysUncertainty);
	  regionEntries.push_back(regionEntry_plus);
	  
	  // tau- candidates only
	  regionEntryType* regionEntry_minus = 
	    new regionEntryType(std::string(regionTakeQCDtemplateFromData).append("-"), 
				std::string(regionWplusJetsSideband).append("-"), 
				std::string(regionStoreQCDtemplate).append("-"), 
				*tauId, *fitVariable, *sysUncertainty);
	  regionEntries.push_back(regionEntry_minus);
	}
      }
    }

    histogramMap4 histograms_mc; // key = (process, region, observable, central value/systematic uncertainty)
    valueMap3 numEvents_mc; // key = (process, region, central value/systematic uncertainty)
    for ( vstring::const_iterator process = processes.begin();
	  process != processes.end(); ++process ) {
      loadHistograms(histograms_mc[*process], histogramInputDirectory, 
		     *process, regions, *tauId, fitVariables, sysUncertainties_expanded);
      
      for ( vstring::const_iterator region = regions.begin();
	    region != regions.end(); ++region ) {
	for ( vstring::const_iterator sysUncertainty = sysUncertainties_expanded.begin();
	      sysUncertainty != sysUncertainties_expanded.end(); ++sysUncertainty ) {
	  numEvents_mc[*process][*region][*sysUncertainty] = 
	    getIntegral(histograms_mc[*process][*region]["EventCounter"][*sysUncertainty], true, true);
	}
      }
    }
    
    histogramMap3 histograms_data; // key = (region, observable, key_central_value)
    loadHistograms(histograms_data, histogramInputDirectory, 
		   "Data", regions, *tauId, fitVariables, sysUncertainties_data);

    for ( std::vector<regionEntryType*>::iterator regionEntry = regionEntries.begin();
	  regionEntry != regionEntries.end(); ++regionEntry ) {
      (*regionEntry)->makeQCDtemplate(histogramOutputDirectory, histograms_data, histograms_mc, numEvents_mc);
    }
  }

  delete histogramInputFile;

//--print time that it took macro to run
  std::cout << "finished executing makeTauIdEffQCDtemplate macro:" << std::endl;
  std::cout << " #tauIdDiscr.  = " << tauIds.size() << std::endl;
  std::cout << " #fitVariables = " << fitVariables.size() << std::endl;
  clock.Show("makeTauIdEffQCDtemplate");

  return 0;
}
