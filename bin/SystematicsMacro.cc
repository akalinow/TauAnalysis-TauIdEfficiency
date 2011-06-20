// STL & system
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iomanip>

using namespace std;

//Root include
#include "TH1F.h"
#include "TMath.h"
#include "TDirectory.h"
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TPolyMarker3D.h>
#include <TPaveText.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TMatrixD.h>

//FWLite Include
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

using namespace edm;

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

struct ratio{
  double Num;
  double Den;
  inline ratio() {Num=0; Den=0;}
};

inline double Ratio(ratio pippo) {return (pippo.Den == 0) ?  -1. : pippo.Num / pippo.Den;}

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // Declare the supported options.
  vector<string> fileNames;
  int report,maxevents;
  bool openfile;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("input-files,f", po::value< vector<string> >(), "Sets the PATTuple files to use")
    //    ("notify-open-file,n",po::value<int>()->default_value(0), "(1/0) Switches the notification of opening file (default 1)")
    ("notify-open-file,n", "Switches the notification of opening file (default true)")
    ("report-every,r",po::value<int>()->default_value(100), " Number of events to process before reporting running again (default 100)")
    ("max-events,m",po::value<int>()->default_value(-1), " Maximum number of events to analyze, -1 for all (default -1)")
    ;

  po::positional_options_description p;
  p.add("input-files", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);    

  if (vm.count("help") || vm.count("h") || argc < 2) {
    cout << desc << "\n";
    return 0;
  }

  openfile = !( vm.count("notify-open-file") || vm.count("n") );
  fileNames = vm["input-files"].as< vector<string> >();
  report =vm["report-every"].as< int >();
  maxevents =vm["max-events"].as< int >();
  
  cout << "Examining " << fileNames.size() << " files"<<endl;
  cout << "Reporting every " << report << " events"<<endl;
  cout << "for " << maxevents << " events"<<endl;
  //  for(int i=1; i<argc; i++)
  //  fileNames.push_back(argv[i]);
  // loop the events
  int ievt=0;

  ratio matchedDecayMode,matchedDecayModeNoVTX, matchedDecayModeLeadTk, matchedDecayModeLeadTkPfIso;
  string names[3];
  names[0] = "HPSLoose";
  names[1] = "HPSMedium";
  names[2] = "HPSTight";
  ratio matchedDecayModeHPSIsoLeadTk[3];
  ratio matchedDecayModeHPSIsoLeadTkPfIso[3];
  ratio c1pRegHPSIsoMatched[3];
  ratio c1fRegHPSIsoMatched[3];
  for(unsigned int iFile=0; iFile<fileNames.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(fileNames[iFile].c_str());
    if( !inFile ){
      cout << "No file matching name:" <<  fileNames[iFile] << " skipping...";
      continue;
    }
    else
      if( openfile)
	cout << "Successfully opened file: "<< fileNames[iFile] << endl;

    fwlite::Event ev(inFile);
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
      // break loop if maximal number of events is reached 
      edm::EventBase const & event = ev;

      // simple event counter
      if(ievt%report==0) 
      	std::cout << "  processing event: " << ievt << std::endl;
     
      if(maxevents != -1 && ievt > maxevents) break;

      edm::Handle<std::vector<pat::Tau> > taus;
      InputTag tauTag("selectedPatPFTausHPSPFRelIsoCumulative");
      event.getByLabel(tauTag, taus);

      edm::Handle<double> weightCollection;
      InputTag weightTag("ntupleProducer","tauIdEffNtuple#addPileupInfo#vtxMultReweight");
      event.getByLabel(weightTag, weightCollection);
      double weight = *weightCollection	;
      //     double weight = 1;

      InputTag ditaucollectionTag("selectedMuPFTauHPSpairsForTauIdEffCumulative");
      Handle< vector<CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> > > muTaus;
      event.getByLabel(ditaucollectionTag,muTaus);
      if( muTaus->size() != taus->size() )
	cout << "pairs: "<< muTaus->size() << " Taus: "<< taus->size()<<endl;
      for(vector<CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> >::const_iterator mutau = muTaus->begin(); mutau != muTaus->end(); ++mutau){
	pat::Muon muon = *(mutau->leg1());
	pat::Tau tau = *(mutau->leg2());
	bool isGenTauLeptonMatched = (tau.genLepton() && TMath::Abs(tau.genLepton()->pdgId()) == 15);
	bool vertexCompatibility = TMath::Abs(tau.vertex().z()-muon.vertex().z()) < 0.2;
	//cout<<"vertex: {"<<tau.vertex().x()<<","<< tau.vertex().y()<<","<< tau.vertex().z()<<"}"<<endl;
	if(tau.vertex().z() ==0 && tau.vertex().x() == 0) {
	  cout << "Warning! IS Tau Vertex set??"<<endl;
	  Handle<reco::VertexCollection> myVertices;
	  InputTag vertexTag("offlinePrimaryVerticesWithBS");
	  event.getByLabel(vertexTag, myVertices);
	  const reco::Vertex& pv = (*myVertices)[0];
	  vertexCompatibility = pv.position() == muon.vertex() && TMath::Abs(tau.leadPFChargedHadrCandsignedSipt()) < 0.2;
	}

	std::string genTauDecayMode = "";
	if ( tau.genJet() != 0 ) {
	  genTauDecayMode = JetMCTagUtils::genTauDecayMode(*tau.genJet());
	}

	bool isGenHadTauDecay = (genTauDecayMode == "oneProng0Pi0"    ||
				 genTauDecayMode == "oneProng1Pi0"    ||
				 genTauDecayMode == "oneProng2Pi0"    ||
				 genTauDecayMode == "oneProngOther"   ||
				 genTauDecayMode == "threeProng0Pi0"  ||
				 genTauDecayMode == "threeProng1Pi0"  ||
				 genTauDecayMode == "threeProngOther" ||
				 genTauDecayMode == "rare");

	bool regionCommonCuts =  muon.pt() > 20. &&
	  tau.pt() > 20. &&
	  tau.pfJetRef()->pt() > 20. &&
	  abs(tau.eta()) < 2.3 &&
	  abs(  tau.pfJetRef()->eta() ) < 2.3 &&
	  (muon.p4() + tau.pfJetRef()->p4()).mass() > 20. &&
	  (muon.p4() + tau.pfJetRef()->p4()).mass() < 200. &&
	  muon.userFloat("pfLooseIsoPt04") < (0.30*muon.pt()) && //In Chris macro was dependent by region but present on all
	  tau.leadPFChargedHadrCand().isAvailable() && //otherwise pair charge is undefined
	  mutau->mt1MET() < 80.;
      //	bool MuonIso_loose = muon.userFloat("pfLooseIsoPt04") < (0.30*muon.pt()); // repeated here to have only this variable
	bool MuonIso_tight = muon.userFloat("pfLooseIsoPt04") < (0.10*muon.pt());
	bool DiTauCharge_OS =  (tau.leadPFChargedHadrCand().isAvailable()) ? (abs(muon.charge() + tau.leadPFChargedHadrCand()->charge()) < 0.5) : false;
	//bool DiTauCharge_SS =  (tau.leadTrack().isNonnull()) ? (abs(muon.charge() + tau.leadTrack()->charge()) > 1.5) : false;
	bool DiTauKine_Sig = mutau->mt1MET() < 40. && (mutau->pZeta() - 1.5*mutau->pZetaVis()) > -20;
	bool c1Reg = regionCommonCuts && MuonIso_tight && DiTauCharge_OS && DiTauKine_Sig;

	bool hasGenPtGt20 =(tau.genJet() != 0 && tau.genJet()->pt() > 20. );
	bool hasGenAbsEtaLs23 = tau.genJet() != 0 && TMath::Abs(tau.genJet()->eta()) < 2.3;
	bool isMatched = isGenTauLeptonMatched && vertexCompatibility && isGenHadTauDecay && hasGenPtGt20 && hasGenAbsEtaLs23;
	bool isMatchedNoVtx = isGenTauLeptonMatched && isGenHadTauDecay && hasGenPtGt20 && hasGenAbsEtaLs23;
	bool hasPFLeadTrackGt5 = (tau.leadPFChargedHadrCand().isAvailable() && tau.leadPFChargedHadrCand()->pt() > 5);
	cout << "lead track dump: isAvailable "<< tau.leadPFChargedHadrCand().isAvailable();
	if(tau.leadPFChargedHadrCand().isAvailable())
	  cout << "  pt" << tau.leadPFChargedHadrCand()->pt();
	cout << endl;
	bool passPfLooseIso = tau.userFloat("pfLooseIsoPt06") < 2.5 ;
	bool passDecayMode = tau.tauID("decayModeFinding") > 0.5;
	bool passLooseIso =  tau.tauID("byLooseIsolation") > 0.5;
	bool passMediumIso =  tau.tauID("byMediumIsolation") > 0.5;
	bool passTightIso =  tau.tauID("byTightIsolation") > 0.5;

	//cout << "byLooseIsolation " << tau.tauID("byLooseIsolation") << "  byMediumIsolation " << tau.tauID("byMediumIsolation") << "  byTightIsolation " << tau.tauID("byTightIsolation") << endl;

	bool passHPSIso[3] = {passLooseIso,passMediumIso,passTightIso}; 

	matchedDecayMode.Num += (isMatched && passDecayMode)*weight;
	matchedDecayMode.Den += (isMatched)*weight;
	matchedDecayModeNoVTX.Num += (isMatchedNoVtx && passDecayMode)*weight;
	matchedDecayModeNoVTX.Den += (isMatchedNoVtx)*weight;
	matchedDecayModeLeadTk.Num += (isMatched && passDecayMode && hasPFLeadTrackGt5)*weight;
	matchedDecayModeLeadTk.Den += (isMatched && passDecayMode)*weight;
	matchedDecayModeLeadTkPfIso.Den += (isMatched && passDecayMode && hasPFLeadTrackGt5 && passPfLooseIso)*weight;
	matchedDecayModeLeadTkPfIso.Num += (isMatched && passDecayMode && hasPFLeadTrackGt5)*weight;

	for(int isoid=0; isoid < 3; isoid++){
	  matchedDecayModeHPSIsoLeadTk[isoid].Num += (isMatched && passDecayMode && passHPSIso[isoid] && hasPFLeadTrackGt5)*weight;
	  matchedDecayModeHPSIsoLeadTk[isoid].Den += (isMatched && passDecayMode && passHPSIso[isoid])*weight;

	  matchedDecayModeHPSIsoLeadTkPfIso[isoid].Num += (isMatched && passDecayMode && passHPSIso[isoid] && hasPFLeadTrackGt5 && passPfLooseIso)*weight;
	  matchedDecayModeHPSIsoLeadTkPfIso[isoid].Den += (isMatched && passDecayMode && passHPSIso[isoid] && hasPFLeadTrackGt5)*weight;

	  c1pRegHPSIsoMatched[isoid].Num += (c1Reg && passHPSIso[isoid] && isMatched)*weight;
	  c1pRegHPSIsoMatched[isoid].Den += (c1Reg && passHPSIso[isoid])*weight;

	  c1fRegHPSIsoMatched[isoid].Num += (c1Reg && !passHPSIso[isoid] && isMatched)*weight;
	  c1fRegHPSIsoMatched[isoid].Den += (c1Reg && !passHPSIso[isoid])*weight;
	}
      }//tau loop
    }//event loop
  if(maxevents != -1 && ievt > maxevents) break;
  }//loop on files

  //ComputeRatio
  cout << "Summary results for SystematicsMacro: "<<endl;
  cout << "   Analyzed "<< ievt <<" events"<<endl;
  cout << "   matchedDecayMode: " << Ratio(matchedDecayMode) << endl;
  cout << "   matchedDecayModeNoVTX: " << Ratio(matchedDecayModeNoVTX) << endl;
  cout << "   matchedDecayModeLeadTk: " << Ratio(matchedDecayModeLeadTk) << endl;
  cout << "   matchedDecayModeLeadTkPfIso: " << Ratio(matchedDecayModeLeadTkPfIso) << endl;

  cout << "   matchedDecayModeHPSIsoLeadTk: ";
  for(int i=0;i<3;i++) cout << Ratio(matchedDecayModeHPSIsoLeadTk[i]) << " ("<< names[i] <<");    ";
  cout<<endl;
  cout << "   matchedDecayModeHPSIsoLeadTkPfIso: ";
  for(int i=0;i<3;i++) cout << Ratio(matchedDecayModeHPSIsoLeadTkPfIso[i]) << " ("<< names[i] <<");    ";
  cout<<endl;
  cout << "   c1pRegHPSIsoMatched: ";
  for(int i=0;i<3;i++) cout << Ratio(c1pRegHPSIsoMatched[i]) << " ("<< names[i] <<");    ";
  cout<<endl;
  cout << "   c1fRegHPSIsoMatched: ";
  for(int i=0;i<3;i++) cout << Ratio(c1fRegHPSIsoMatched[i]) << " ("<< names[i] <<");    ";
  cout<<endl;
}
