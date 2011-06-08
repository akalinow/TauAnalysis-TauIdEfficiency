// -*- C++ -*-
//
// Package:    HistoPlotter
// Class:      HistoPlotter
// 
/**\class HistoPlotter HistoPlotter.cc TauAnalysis/HistoPlotter/src/HistoPlotter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauro Verzetti,15 R-035,+41227678587,
//         Created:  Tue May 17 10:06:53 CEST 2011
// $Id: HistoPlotter.cc,v 1.1 2011/06/06 09:00:23 mverzett Exp $
//
//


// STL & system
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>

using namespace std;

//Root include
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TCanvas.h"

using namespace TMath;

// FWK include
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//FWK classes
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>

using namespace pat;
using namespace edm;

//
// class declaration
//
typedef vector<string> vstring;
typedef vector<ParameterSet> vpset;

class HistoPlotter : public edm::EDAnalyzer {
public:
  explicit HistoPlotter(const edm::ParameterSet&);
  ~HistoPlotter();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  inline string DetachPath(string );
  inline string DetachFileName(string );
  inline string GetProducerName(string );
  bool CommonCuts(const CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> &); //cut common to all the regions
  bool VertexMatch(Ptr<pat::Muon>,Ptr<pat::Tau>);//checks that both mu and tau come from the same vertex
  string AssignRegion(const CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> &,ParameterSet &); //Assigna region in ABCD scheme
  string replace(std::string& , const std::string& , const std::string& );

  // ----------member data ---------------------------

  string outFileName_;
  InputTag evtCounter_;
  TFile *outFile_;
  map<string,TDirectory *> dir_;
  vstring sysUncertainties_;
  vstring regions_;
  bool plotSys_;
  string channel_;
  vpset tauids_; //{name, taucollection, ditaucollection, discName (pat discriminator name)} get them @ line 102. discName is what gets written in the tautuple, redundant? should change name to be = discName
  vpset plotVariables_; //{name, numBins,minh,maxh,expr} GetThem @ line 333 expr is how to get the variable from the pair
  ParameterSet muonCollections_;
  map<string,TH1F*> histoMap_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HistoPlotter::HistoPlotter(const edm::ParameterSet& iConfig):
  outFileName_(iConfig.getParameter<string>("outputFile").c_str()),
  evtCounter_(iConfig.getParameter<InputTag>("evtCounter")),
  channel_(iConfig.getParameter<string>("channel").c_str()),
  sysUncertainties_(iConfig.getParameter<vstring>("sysUncertainties")),
  regions_(iConfig.getParameter<vstring>("regions")),
  tauids_(iConfig.getParameter<vpset>("tauids")),
  plotVariables_(iConfig.getParameter<vpset>("plotVariables")),
  muonCollections_(iConfig.getParameter<ParameterSet>("muonCollections")),
  plotSys_(iConfig.getParameter<bool>("plotSys")),
  histoMap_()
{
  cout<<"<HistoPlotter::HistoPlotter>"<<endl;
  //Creates output file
  outFile_ = new TFile(outFileName_.c_str(),"recreate");
  cout << "created output file: "<< outFileName_.c_str()<<endl;
  //Creates Dir Structure: channel/tauid.name/region/plotVariable.name/plotVariable.name+"Distribution"+sysShift
  dir_[channel_.c_str()] = outFile_->mkdir(channel_.c_str(),channel_.c_str());
  dir_[channel_.c_str()]->cd();
  cout << "crated directory "<< channel_ << endl;
  cout << "tauids_.size() = " << tauids_.size() << endl;
  histoMap_[(channel_+"/AnalyzedEvents").c_str()] = new TH1F((channel_+"/AnalyzedEvents").c_str(),(channel_+"/AnalyzedEvents").c_str(),1,0,1); 
  for(vpset::iterator tauid = tauids_.begin(); tauid != tauids_.end(); tauid++){
    string idname = tauid->getParameter<string>("name");
    string chTau = channel_+"/"+idname;
    dir_[chTau.c_str()] = dir_[channel_.c_str()]->mkdir(idname.c_str(),idname.c_str());
    cout << "crated directory "<< idname << endl;
    for(vstring::iterator region = regions_.begin(); region != regions_.end(); region++){
      string chTauReg = chTau + "/" + *region;
      dir_[chTauReg.c_str()] =  dir_[chTau.c_str()]->mkdir(region->c_str(),region->c_str());
      cout << "crated directory "<< *region << endl;
      for(vpset::iterator plotVariable = plotVariables_.begin(); plotVariable != plotVariables_.end(); plotVariable++){
	string varname = plotVariable->getParameter<string>("name");
	int numBins = plotVariable->getParameter<int>("numBins");
	double min = plotVariable->getParameter<double>("minh");
	double max = plotVariable->getParameter<double>("maxh");
	string chTauRegVar = chTauReg + "/" + varname;
	dir_[chTauRegVar.c_str()] = dir_[chTauReg.c_str()]->mkdir(varname.c_str(),varname.c_str());
	cout << "crated directory "<< varname << endl;
	dir_[chTauRegVar.c_str()]->cd();

	//books Histograms
	string histoname = varname+"Distribution";
	if(sysUncertainties_.size() == 0){
	  histoMap_[(chTauRegVar+"/"+histoname).c_str()] = new TH1F(histoname.c_str(),histoname.c_str(),numBins,min,max);
	  cout<<" Created Histogram: "<<histoname.c_str();
	  sysUncertainties_.push_back("CENTRAL_VALUE");
	}
	else
	  for(vstring::iterator sysShift = sysUncertainties_.begin(); sysShift != sysUncertainties_.end(); sysShift++){
	    if(*sysShift == "CENTRAL_VALUE"){
	      histoMap_[(chTauRegVar+"/"+histoname).c_str()] = new TH1F(histoname.c_str(),histoname.c_str(),numBins,min,max);
	      cout<<" Created Histogram: "<<histoname.c_str() << endl;
	    }
	    else{
	      histoMap_[(chTauRegVar+"/"+histoname+"_"+*sysShift).c_str()] = new TH1F((histoname+"_"+*sysShift).c_str(),(histoname+"_"+*sysShift).c_str(),numBins,min,max);
	      cout<<" Created Histogram: "<<(histoname+"_"+*sysShift) << endl;
	    }
	  }
      }//variable loop
    }//region loop
  }//tauid loop
}


HistoPlotter::~HistoPlotter()
{
  cout<<"<HistoPlotter::~HistoPlotter>"<<endl;
  //Writes the file

  //Delete histograms
  for(std::map<string,TH1F*>::const_iterator histo = histoMap_.begin(); histo != histoMap_.end(); histo++){
    dir_[histo->first.substr(0,histo->first.rfind("/"))]->cd();
    histo->second->Write();
    delete histo->second;
  }
  //cout << "file Write: " <<  outFile_->Write() <<endl;
  outFile_->Close();
}


//
// member functions
//

string HistoPlotter::replace(std::string& str, const std::string& oldpart, const std::string& newpart) 
{
  //std::cout << "<replace>:" << std::endl;

  std::string retVal = str;

  size_t idx = retVal.find(oldpart);
  while ( idx != std::string::npos ) {
    retVal.replace(idx, oldpart.length(), newpart);
    idx += newpart.length();
    if ( idx < retVal.length() )
      idx = retVal.find(oldpart, idx);
    else
      idx = std::string::npos;
    //std::cout << "retVal = " << retVal << std::endl;
  }

  return retVal;
}

void HistoPlotter::endLuminosityBlock(const edm::LuminosityBlock &lumi, const edm::EventSetup &setup)
{
  Handle<MergeableCounter> numEventsCounter;
  lumi.getByLabel(evtCounter_, numEventsCounter);
  
  if (numEventsCounter.isValid()) {
    float oldCont = histoMap_[(channel_+"/AnalyzedEvents").c_str()]->GetBinContent(1);
    oldCont += numEventsCounter->value;
    histoMap_[(channel_+"/AnalyzedEvents").c_str()]->GetBinContent(1,oldCont);
  }
}


bool HistoPlotter::VertexMatch(Ptr<pat::Muon> muon,Ptr<pat::Tau> tau) //returns if muon and tau have the same vertex.
{
  reco::Candidate::Point mounVtx = muon->vertex();
  reco::Candidate::Point tauVtx = muon->vertex();
  if(tauVtx.x() == 0 && tauVtx.y() == 0 && tauVtx.z() == 0){ //if vertex matching is not defined in this code tag
    Handle<reco::VertexCollection> myVertices;
    const reco::Vertex& pv = (*myVertices)[0];
    return pv.position() == mounVtx && abs(tau->leadPFChargedHadrCandsignedSipt()) < 0.2;
  }
  else //if vertex matching between tau and vertex is defined
    return mounVtx == tauVtx;
}

bool HistoPlotter::CommonCuts(const CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> & muTauPair)
{
  Ptr<pat::Muon> muon = muTauPair.leg1();
  Ptr<pat::Tau> tau = muTauPair.leg2();
  return muon->pt() > 20. &&
    tau->pt() > 20. &&
    tau->pfJetRef()->pt() > 20. &&
    abs(tau->eta()) < 2.3 &&
    abs(  tau->pfJetRef()->eta() ) < 2.3 &&
    tau->userFloat("pfLooseIsoPt06") < 2.5 &&
    (muon->p4() + tau->pfJetRef()->p4()).mass() > 20. &&
    (muon->p4() + tau->pfJetRef()->p4()).mass() < 200. &&
    muon->userFloat("pfLooseIsoPt04") < (0.30*muon->pt()) && //In Chris macro was dependent by region but present on all
    tau->leadTrack().isNonnull() && //otherwise pair charge is undefined
    muTauPair.mt1MET() < 80.;
}

string HistoPlotter::AssignRegion(const CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> &muTauPair, ParameterSet & tauid)
{
  Ptr<pat::Muon> muon = muTauPair.leg1();
  Ptr<pat::Tau> tau = muTauPair.leg2();
  bool MuonIso_loose = muon->userFloat("pfLooseIsoPt04") < (0.30*muon->pt()); // repeated here to have only this variable
  bool MuonIso_tight = muon->userFloat("pfLooseIsoPt04") < (0.10*muon->pt());
  bool DiTauCharge_OS =  (tau->leadTrack().isNonnull()) ? (abs(muon->charge() + tau->leadTrack()->charge()) < 0.5) : false;
  bool DiTauCharge_SS =  (tau->leadTrack().isNonnull()) ? (abs(muon->charge() + tau->leadTrack()->charge()) > 1.5) : false;
  if(!DiTauCharge_OS && !DiTauCharge_SS) return "Undefined"; //may happen that no charge is defined for tau
  bool DiTauKine_Sig = muTauPair.mt1MET() < 40. && (muTauPair.pZeta() - 1.5*muTauPair.pZetaVis()) > -20;
  //bool DiTauKine_Bgr  = !DiTauKine_Sig;

  string retstring("");
  //Region definition
  if(MuonIso_loose && !MuonIso_tight && DiTauCharge_OS)
    retstring.append("A");
  else if(MuonIso_loose && !MuonIso_tight && DiTauCharge_SS)
    retstring.append("B");
  else if(MuonIso_tight && DiTauCharge_OS)
    retstring.append("C");
  else if(MuonIso_tight && DiTauCharge_OS)
    retstring.append("D");
  //signal / background-like signal
  if(DiTauKine_Sig)
    retstring.append("1");
  else
    retstring.append("2");

  string discName = tauid.getParameter<string>("discName");
  if(tau->isTauIDAvailable(discName)){
    if(tau->tauID(discName) > 0.5)
      retstring.append("p");
    else if(tau->tauID(discName) < 0.5)
      retstring.append("f");
  }
  //  else
  //  return "Undefined"; //to keep the number of events per region equal
  return retstring;
}

// ------------ method called to for each event  ------------
void HistoPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  InputTag standardMuons = muonCollections_.getParameter<InputTag>("standardMuons"); //selectedPatMuonsForTauIdEffTrkIPcumulative
  //InputTag globalMuons = muonCollections_.getParameter<InputTag>("globalMuons"); //patMuonsGlobal
  InputTag standAloneMuons = muonCollections_.getParameter<InputTag>("standAloneMuons"); //patMuonsStandAlone

  Handle<vector<Muon> > muonHandle;
  iEvent.getByLabel(standardMuons,muonHandle);

  Handle<vector<Muon> > muonStandAloneHandle;
  iEvent.getByLabel(standardMuons,muonStandAloneHandle);
  if(muonStandAloneHandle->size() != 1 || muonHandle->size() != 1)
    return;

  for(vpset::iterator tauid = tauids_.begin(); tauid != tauids_.end(); tauid++){
    string tauidname = tauid->getParameter<string>("name");
    for(vstring::iterator sysShift = sysUncertainties_.begin(); sysShift !=sysUncertainties_.end(); sysShift++){
      string taucollectionbare = tauid->getParameter<string>("taucollection");
      string taucollection = (*sysShift == "CENTRAL_VALUE") ? taucollectionbare : replace(taucollectionbare,   "Cumulative", std::string(*sysShift).append("Cumulative"));
      InputTag taucollectionTag(taucollection.c_str());
      Handle< vector<pat::Tau> > tauHandle;
      iEvent.getByLabel(taucollectionTag,tauHandle);

      string ditaucollectionbare = tauid->getParameter<string>("ditaucollection");
      string ditaucollection = (*sysShift == "CENTRAL_VALUE") ? ditaucollectionbare : replace(ditaucollectionbare,   "Cumulative", std::string(*sysShift).append("Cumulative"));
      InputTag ditaucollectionTag(ditaucollection.c_str());
      Handle< vector<CompositePtrCandidateT1T2MEt<pat::Muon,pat::Tau> > > ditauHandle;
      iEvent.getByLabel(ditaucollectionTag,ditauHandle);
      
      for( vector<CompositePtrCandidateT1T2MEt<Muon,Tau> >::const_iterator muTauPair= ditauHandle->begin() ; muTauPair != ditauHandle->end(); muTauPair++ ){
	Ptr<pat::Muon> muon = muTauPair->leg1();
	Ptr<pat::Tau> tau = (*muTauPair).leg2();
	if(!VertexMatch( muon, tau)) continue;
	if( CommonCuts(*muTauPair) ) {
	  string assignedRegion = AssignRegion(*muTauPair,*tauid);
	  if(assignedRegion != "Undefined"){
	    for(vstring::iterator region = regions_.begin(); region != regions_.end(); region++){
	      if(*region == "ABCD" || assignedRegion.find(*region) != string::npos){
		for(vpset::iterator plotVariable = plotVariables_.begin(); plotVariable != plotVariables_.end(); plotVariable++){
		  // DIR STRUCT: channel/tauid.name/region/plotVariable.name/plotVariable.name+"Distribution"+sysShift
		  string varname = plotVariable->getParameter<string>("name");
		  string histoname = channel_+"/"+tauidname+"/"+*region+"/"+varname+"/"+varname+"Distribution";
		  if(*sysShift != "CENTRAL_VALUE")  histoname.append("_"+*sysShift);
		  StringObjectFunction< CompositePtrCandidateT1T2MEt<Muon,Tau> > vtxFunction_(plotVariable->getParameter<string>("expr"));
		  double varValue = vtxFunction_(*muTauPair);
		  histoMap_[histoname]->Fill(varValue);
		}//variable loop
	      }//if(*region == "ABCD" || assignedRegion.find(*region) != string::npos){
	    }//region loop
	  }//if(assignedRegion != "Undefined")
	}//if(CommonCuts(*muTauPair))
      } //loop over pairs    
    }//loop over sysShifts
  }//loop over tau ID
}


// ------------ method called once each job just before starting event loop  ------------
void HistoPlotter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HistoPlotter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HistoPlotter);
