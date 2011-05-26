#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <TH2F.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
/////////////////////Added By Raman///////////////////////////////////
#include "/cmshome/khurana/ForFWLite/CMSSW_4_1_4/src/PhysicsTools/FWLite/bin/Book_Fill_histograms.cc"
#include <TH1I.h>
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/Math/interface/deltaPhi.h"
////////////////////Added upto this comment//////////////////////////


#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"


int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  //  using reco::Muon;

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // initialize command line parser
  
  ///////////////////////////////////////////////////////////////////////////
  /////////////Function to Book histograms//////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  Book_Fill_histograms *a = new Book_Fill_histograms();      
  a->bookHistograms();
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");
  // set defaults
  parser.integerValue ("maxEvents"  ) = 100000;
  parser.integerValue ("outputEvery") =   1;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";
  
  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");
  
  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  
  //  TH1F* tauPhi_ = dir.make<TH1F>("tauPhi" , "phi" ,   100,  -5.,   5.);  
  TH1F* tauMass_= dir.make<TH1F>("tauMass", "mass",    10,  0.0,  10.);

  TH1F *tauPt_;
  tauPt_ = dir.make<TH1F>("tauPt" , "pt"  ,   900,   0., 300.);
  
  TH1F* tauEta_;
  tauEta_ = dir.make<TH1F>("tauEta"  , "eta"  ,   100,   -3.0, 3.0);
  
  TH1F* tauPhi_;
  tauPhi_  = dir.make<TH1F>("tauPhi"  , "phi"  ,   100,   -5.0, 5.0);
    
  TH1F* PionSumPt_;
  PionSumPt_  = dir.make<TH1F>("PionSumPt"  , "PionSumPt"  ,   100,   0.0, 150.0);
  
  TH1F* PhotonSumPt_;
  PhotonSumPt_  = dir.make<TH1F>("PhotonSumPt"  , "PhotonSumPt"  ,   100,   -5.0, 150.0);
  
  TH1F* genTauPt_;
  genTauPt_  = dir.make<TH1F>("genTauPt" , "genTauPt"  ,   900,   0.0, 300.0);
  
  TH1F* genTauEta_;
  genTauEta_  = dir.make<TH1F>("genTauEta"  , "genTauEta"  ,   100,   -5.0, 5.0);
  
  TH1F* genTauPhi_;
  genTauPhi_  = dir.make<TH1F>("genTauPhi"  , "genTauPhi"  ,   100,   -5.0, 5.0);

  TH1F* gentaudecaymode_;
  gentaudecaymode_  = dir.make<TH1F>("gentaudecaymode"  , "gentaudecaymode"  ,   15,   0 , 15.0);

  TH1F* recotaudecaymode_;
  recotaudecaymode_  = dir.make<TH1F>("recotaudecaymode"  , "recotaudecaymode"  ,   15,   0.0 , 15.0);
  
  TH2F* reco_vs_gen_decay_mode_;
  reco_vs_gen_decay_mode_ =  dir.make<TH2F>("reco_vs_gen_taudecaymode"  , "reco_vs_gen_taudecaymode;decay mode_{reconstructed};decay mode_{generated}"  ,   15,   0.0 , 15.0 ,   15,   0.0 , 15.0);
  
  TH1I* isoparticleID_;
  isoparticleID_  = dir.make<TH1I>("isoparticleID"  , "isoparticleID"  ,   10,   -0.5, 9.5);
  //for(int ii=0; ii<15;ii++)
  //{
  //  char tmphistname[100];
  //  sprintf(tmphistname,"_%d",ii+1);
  //  std::string histname(tmphistname);
  //isoparticleID_[ii]  = dir.make<TH1I>(("isoparticleID"+histname).c_str()  , "isoparticleID"  ,   10,   -0.5, 9.5);
  //}
  
  TH1F* PFJetPt_;
  PFJetPt_  = dir.make<TH1F>("PFJetPt"  , "PFJetPt"  ,   900,   0.0, 300.0);
  
  TH1F* PFJetPx_;
  PFJetPx_  = dir.make<TH1F>("PFJetPx" , "PFJetPx"  ,   300,   -150.0, 150.0);
    
  TH1F*PFJetPy_;
  PFJetPy_ = dir.make<TH1F>("PFJetPy"  , "PFJetPy"  ,   300,   -150.0, 150.0);
  
  TH1F* PFJetPz_;
  PFJetPz_  = dir.make<TH1F>("PFJetPz" , "PFJetPz"  ,   300,   -150.0, 150.0);
  
  TH1F* PFJetPhi_;
  PFJetPhi_ = dir.make<TH1F>("PFJetPhi"  , "PFJetPhi"  ,   100,   -5.0, 5.0);
  
  TH1F* PFJetEta_;
  PFJetEta_  = dir.make<TH1F>("PFJetEta"  , "PFJetEta"  ,   100,   -3.0, 3.0);
  
  TH1F* PFJetEnergy_;
  PFJetEnergy_  = dir.make<TH1F>("PFJetEnergy"  , "PFJetEnergy"  ,   500,   0.0, 500.0);
  
  TH1F* PtRes_;
  PtRes_  = dir.make<TH1F>("PtRes" , "PtRes"  ,   100,   0.0, 5.0);
  
  TH1F* PTResOneProng0Pi0_;
  PTResOneProng0Pi0_  = dir.make<TH1F>("PTResOneProng0Pi0" , "PTResOneProng0Pi0"  ,   100,   0.0, 5.0);
  
  TH1F* PTResOneProng1Pi0_;
  PTResOneProng1Pi0_  = dir.make<TH1F>("PTResOneProng1Pi0" , "PTResOneProng1Pi0"  ,   100,   0.0, 5.0);
  
  TH1F* PTResOneProng2Pi0_;
  PTResOneProng2Pi0_  = dir.make<TH1F>("PTResOneProng2Pi0" , "PTResOneProng2Pi0"  ,   100,   0.0, 5.0);
  
  TH1F* PTResThreeProng0Pi0_;
  PTResThreeProng0Pi0_  = dir.make<TH1F>("PTResThreeProng0Pi0" , "PTResThreeProng0Pi0"  ,   100,   0.0, 5.0);
  
  TH1F* PTResThreeProng1Pi0_;
  PTResThreeProng1Pi0_  = dir.make<TH1F>("PTResThreeProng1Pi0" , "PTResThreeProng1Pi0"  ,   100,   0.0, 5.0);
  
  TH1F* matched_PtRes_;
  matched_PtRes_ = dir.make<TH1F>("Matched_PtRes" , "Matched_PtRes; Pt_{Res};;"  ,   100,   0.0, 5.0);

  TH1F* PFCPtSum_;
  PFCPtSum_  = dir.make<TH1F>("PFCPtSum"  , "PFCPtSum"  ,   100,   0.0, 300.0);
  
  TH1F* PFGPtSum_;
  PFGPtSum_  = dir.make<TH1F>("PFGPtSum"  , "PFGPtSum"  ,   100,   0.0, 300.0);
  
  TH1F* PFNPtSum_;
  PFNPtSum_  = dir.make<TH1F>("PFNPtSum"  , "PFNPtSum"  ,   100,   0.0, 300.0);
    
  TH1F* PFCandidateSumPt_;
  PFCandidateSumPt_ = dir.make<TH1F>("PFCandidateSumPt"  , "PFCandidateSumPt"  ,   100,   0.0, 300.0);
  
  TH1F* PFCandidateSumPt_div_genJetPt_;
  PFCandidateSumPt_div_genJetPt_  = dir.make<TH1F>("PFCandidateSumPt_div_genJetPt"  , "PFCandidateSumPt_div_genJetPt"  ,   100,   0.0, 5.0);
    
  TH1F* PFCisoPtSum_;
  PFCisoPtSum_ = dir.make<TH1F>("PFCisoPtSum"  , "PFCisoPtSum"  ,   100,   0.0, 300.0);
  
  TH1F* PFGisoPtSum_;
  PFGisoPtSum_ = dir.make<TH1F>("PFGisoPtSum"  , "PFGisoPtSum"  ,   100,   0.0, 300.0);
    
  TH1F* PFNisoPtSum_;
  PFNisoPtSum_  = dir.make<TH1F>("PFNisoPtSum"  , "PFNisoPtSum"  ,   100,   0.0, 300.0);
  
  TH1F* PFCandidateisoSumPt_;
  PFCandidateisoSumPt_  = dir.make<TH1F>("PFCandidateisoSumPt" , "PFCandidateisoSumPt"  ,   100,   0.0, 300.0);
  
  TH1F* OneMinusPtRes_;
  OneMinusPtRes_ = dir.make<TH1F>("OneMinusPtRes"  , "OneMinusPtRes"  ,   100,   -5.0, 5.0);
  
  TH1F* PFCandidateisoSumPtDivGenTauPt_;
  PFCandidateisoSumPtDivGenTauPt_  = dir.make<TH1F>("PFCandidateisoSumPtDivGenTauPt" , "PFCandidateisoSumPtDivGenTauPt_" , 100, -5.0, 5.0);
  
  TH1F* GenMinusPatPt_;
  GenMinusPatPt_  = dir.make<TH1F>("GenMinusPatPt"  , "GenMinusPatPt" , 100, 0.0, 300.0);
  
  TH1F* PTResPFJet_[15];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      PTResPFJet_[ibin]  = dir.make<TH1F>(("PTResPFJet"+histname).c_str()  , "PTResPFJet" , 300, 0.0, 10.0);
    }

  TH2F* reco_vs_gen_decay_mode_in_bin_[15];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      reco_vs_gen_decay_mode_in_bin_[ibin]  = dir.make<TH2F>(("reco_vs_gen_decay_mode_in_bin_"+histname).c_str()  , "reco_vs_gen_decay_mode;decaymode_{reco};decaymode_{gen}" , 15, 0.0, 15.0,15,0,15);
    }
  
  TH2F* SumPFN_vs_SumPFG_[15];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      SumPFN_vs_SumPFG_[ibin] = dir.make<TH2F>(("SumPFN_vs_SumPFG"+histname).c_str()  , "SumPFN_vs_SumPFG" ,1000,0.0,10, 1000, 0.0, 10.0);
    }
  
  TH1F* SigPlusIso_div_genJetPt_;
  SigPlusIso_div_genJetPt_ = dir.make<TH1F>("SigPlusIso_div_genJetPt" , "SigPlusIso_div_genJetPt_;Sig Plus Iso Pt Sum; # of events;" , 100,   0.0, 5.0);
					    
  TH1F* PTRESPFJET_;
  PTRESPFJET_  = dir.make<TH1F>("PTRESPFJET" , "PTRESPFJET"  ,   100,   0.0, 5.0);

  TH1F* checkbugptres_;
  checkbugptres_  = dir.make<TH1F>("checkbugptres", "checkbugptres"  ,   900,   0.0, 300.0);

  TH1F* checkgentaupt_;
  checkgentaupt_  = dir.make<TH1F>("checkgentaupt", "checkgentaupt"  ,   900,   0.0, 300.0);
				   
  TH1F* checkpattaupt_;
  checkpattaupt_  = dir.make<TH1F>("checkpattaupt", "checkpattaupt"  ,   900,   0.0, 300.0);

  TH1F* nPi0_;
  nPi0_  = dir.make<TH1F>("nPi0", "nPi0"  ,   100,   0.0, 5.0);

  TH1F* dPhipi0_;
  dPhipi0_  = dir.make<TH1F>("dPhipi0", "dPhipi0"  ,   100,   -3.5 , 3.5);
				   
  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------      
      fwlite::Event ev(inFile);
      
      Int_t itau;
      float PionSumPt;
      float PhotonSumPt;
      float tmpPFCPtSum;
      float tmpPFGPtSum;
      float tmpPFNPtSum;
      
      float tmpPFCisoPtSum;
      float tmpPFGisoPtSum;
      float tmpPFNisoPtSum;
      
      
      float genTauPt ;
      float tauPt;
      float PFCisoPtSum;
      float PFGisoPtSum;
      float PFNisoPtSum;
      int gendecaymodeInt =-999;
      float PtRes;
      int recotaudecaymode;
      int nPi0;
      float pi0Phi[10];
      float PFJetPt ;
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
	  std::cout << "  processing event: " << ievt << std::endl;

	using namespace edm;
	using namespace std;
	using namespace reco;
	using namespace pat;
	//////////////////////////////////////////////////////////////
	//////////////////////////Clean Pat Taus//////////////////////
	//////////////////////////////////////////////////////////////
	edm::Handle<std::vector<pat::Tau> >  tauHandle;
	event.getByLabel(std::string("selectedPatTaus"),tauHandle);
	itau =0;
	/*
	  tauPt =-999;
	genTauPt =-999;
	gendecaymodeInt =-999;
	PFJetPt = -999;
	PtRes = -999;
	recotaudecaymode =-999;
	*/
	//std::string gendecaymodestring;
	for(std::vector<pat::Tau>::const_iterator tauItr=tauHandle->begin(); tauItr!=tauHandle->end(); ++tauItr)
	  {
	    tauPt =-999;
	    genTauPt =-999;
	    gendecaymodeInt =-999;
	    PFJetPt = -999;
	    PtRes = -999;
	    recotaudecaymode =-999;
	    if(tauItr->genJet() )
	      {
		std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(*tauItr->genJet());
		if (  (genTauDecayMode == "oneProng0Pi0"   ||
		       genTauDecayMode == "oneProng1Pi0"   ||
		       genTauDecayMode == "oneProng2Pi0"   ||
		       genTauDecayMode == "threeProng0Pi0" ||
		       genTauDecayMode == "threeProng1Pi0" ))
		  {
		    //cout<<" gen tau jet pt in fwlite ="<<tauItr->genJet()->pt()<<endl;
		    genTauPt = tauItr->genJet()->pt();
		    genTauPt_->Fill(tauItr->genJet()->pt()) ;
		    genTauEta_->Fill(tauItr->genJet()->eta());
		    genTauPhi_->Fill(tauItr->genJet()->phi());
		    if(genTauDecayMode=="oneProng0Pi0")
		      gendecaymodeInt=0;
		    if(genTauDecayMode=="oneProng1Pi0")
		      gendecaymodeInt=1;
		    if(genTauDecayMode=="oneProng2Pi0")
		      gendecaymodeInt=2;
		    if(genTauDecayMode=="threeProng0Pi0")
		      gendecaymodeInt=10;
		    if(genTauDecayMode=="threeProng1Pi0")
		      gendecaymodeInt=11;
		    gentaudecaymode_->Fill(gendecaymodeInt);
		    if( (tauItr->genJet()->getGenConstituents()).size() !=0 )
		      {
			std::vector <const GenParticle*>  genParticleList = tauItr->genJet()->getGenConstituents();
			int ipion =0;
			int iphoton = 0;
			int ipi0 = 0;
			PionSumPt=0.0;
			PhotonSumPt=0.0;
			for (int iList = 0; iList <int(genParticleList.size()); iList++)
			  {
			    if(abs(genParticleList[iList]->pdgId()) == 211 || abs(genParticleList[iList]->pdgId()) == 321)
			      {
				PionSumPt+= genParticleList[iList]->pt();
			      }
			    
			    if(abs(genParticleList[iList]->pdgId()) == 22)
			      {
				PhotonSumPt+=genParticleList[iList]->pt();
			      }
			    
			    if(abs(genParticleList[iList]->pdgId()) == 111  )
			      {
				pi0Phi[ipi0]=genParticleList[iList]->phi();
				ipi0++;
			      }
			    nPi0 = ipi0;
			    nPi0_->Fill(nPi0);
			    if(nPi0 >=2)
			      {
				dPhipi0_->Fill(deltaPhi(pi0Phi[0],pi0Phi[1]));
			      }
			  }//for (int iList = 0; iList <int(genParticleList.size()); iList++)...{
			PionSumPt_->Fill(PionSumPt);  
			PhotonSumPt_->Fill(PhotonSumPt);
		      }//if( (tauItr->genJet()->getGenConstituents()).size() !=0 )...{
		  }//genTauDecayMode == "threeProng1Pi0" ))...{
		//////////////////////////////////////////////////////
		///////////////gen level info finished///////////////
		////////////////////////////////////////////////////
		const reco::PFCandidateRefVector &chargedref = tauItr->signalPFChargedHadrCands() ;
		tmpPFCPtSum = 0.0;
		for(int i=0; i<(int)chargedref.size(); i++)
		  {
		    tmpPFCPtSum = tmpPFCPtSum + chargedref[i]->pt();
		  }
		PFCPtSum_->Fill(tmpPFCPtSum);
		tmpPFGPtSum = 0.0;
		const reco::PFCandidateRefVector &gammaref = tauItr->signalPFGammaCands();
		for(int i=0; i<(int)gammaref.size(); i++)
		  {
		    tmpPFGPtSum = tmpPFGPtSum + gammaref[i]->pt();
		  }
		tmpPFNPtSum  = 0.0;
		const reco::PFCandidateRefVector &neutralref = tauItr->signalPFNeutrHadrCands();
		for(int i=0; i<(int)neutralref.size(); i++)
		  {
		    tmpPFNPtSum = tmpPFNPtSum + neutralref[i]->pt();
		  }
		PFNPtSum_->Fill(tmpPFNPtSum);
		//////////////////////////////////////////////////////////////////////
		//////////Sum of PF Candidates in Sognal Region///////////////////////
		//////////////////////////////////////////////////////////////////////
		PFCandidateSumPt_->Fill(tmpPFGPtSum+tmpPFCPtSum);
		PFCandidateSumPt_div_genJetPt_->Fill((tmpPFGPtSum+tmpPFCPtSum)/genTauPt);
		
		/////////////////////////////////////////////////////////////////////////////////
		/////////////Isolation Candidates////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////
		const reco::PFCandidateRefVector &chargedisoref = tauItr->isolationPFChargedHadrCands() ;
		tmpPFCisoPtSum = 0.0;
		for(int i=0; i<(int)chargedisoref.size(); i++)
		  {
		    tmpPFCisoPtSum = tmpPFCisoPtSum + chargedisoref[i]->pt() ;
		  }
		PFCisoPtSum = tmpPFCisoPtSum;
		PFCisoPtSum_->Fill(tmpPFCisoPtSum);
		const reco::PFCandidateRefVector &gammaisoref = tauItr->isolationPFGammaCands();
		tmpPFGisoPtSum = 0.0;
		for(int i=0; i<(int)gammaisoref.size(); i++)
		  {
		    tmpPFGisoPtSum = tmpPFGisoPtSum + gammaisoref[i]->pt() ;
		  }
		PFGisoPtSum = tmpPFGisoPtSum;
		PFGisoPtSum_->Fill(tmpPFGisoPtSum);
		const reco::PFCandidateRefVector &neutralisoref = tauItr->isolationPFNeutrHadrCands();
		tmpPFNisoPtSum = 0.0;
		for(int i=0; i<(int)neutralisoref.size(); i++)
		  {
		    tmpPFNisoPtSum = tmpPFNisoPtSum + neutralisoref[i]->pt() ; 
		  }
		PFNisoPtSum = tmpPFNisoPtSum;
		PFNisoPtSum_->Fill(tmpPFNisoPtSum);
		////////////////////////////////////////////////////////////////////
		//////////////////////Sum of iso candidates////////////////////////
		///////////////////////////////////////////////////////////////////
		PFCandidateisoSumPt_->Fill(tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum);
		PFCandidateisoSumPtDivGenTauPt_->Fill((tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum)/genTauPt);
		const PFCandidateRefVector &pfcandidateisoref = tauItr->isolationPFCands();
		int size = pfcandidateisoref.size();
		for (int i = 0; i<size ; i++)
		  {
		    isoparticleID_->Fill(pfcandidateisoref[i]->particleId());
		  }
				
		/////////////////////////////////////////////////////////////////// 
		///////////////Particle Flow Jets////////////////////////////////// 
		/////////////////////////////////////////////////////////////////// 
		
		if ( tauItr->pfJetRef().isNonnull() )
		  {
		    PFJetPt=tauItr->pfJetRef()->pt();
		    
		    //cout<<"pt of ref jet in fwlite ="<<tauItr->pfJetRef()->pt()<<endl;
		    PFJetPt_->Fill(tauItr->pfJetRef()->pt());
		    PFJetPx_->Fill(tauItr->pfJetRef()->px());
		    PFJetPy_->Fill(tauItr->pfJetRef()->py());
		    PFJetPz_->Fill(tauItr->pfJetRef()->pz());
		    PFJetPhi_->Fill(tauItr->pfJetRef()->phi());
		    PFJetEta_->Fill(tauItr->pfJetRef()->eta());
		    PFJetEnergy_->Fill(tauItr->pfJetRef()->energy());
		  }
		
		
		if( tauItr->tauID("decayModeFinding")            &&
		    tauItr->tauID("byLooseIsolation")            &&
		    tauItr->tauID("againstElectronTight")        &&
		    tauItr->tauID("againstMuonTight")            &&  
		    abs(tauItr->eta()) < 2.3                     
		    )
		  {
		    if((tauItr->pt() > 20.0 && tauItr->pt() < 40.0)     ) 
		      {    
			tauPt = tauItr->pt();
			tauPt_->Fill(tauItr->pt());
			tauEta_->Fill(tauItr->eta());
			tauPhi_->Fill(tauItr->phi());
			recotaudecaymode = tauItr->decayMode();
			recotaudecaymode_->Fill(tauItr->decayMode());
			reco_vs_gen_decay_mode_->Fill(recotaudecaymode,gendecaymodeInt);
			if(recotaudecaymode == gendecaymodeInt)
			  {
			    matched_PtRes_->Fill(tauPt/genTauPt);
			  }
			PtRes = tauPt/genTauPt;
			PtRes_->Fill(tauPt/genTauPt);
			
			if(gendecaymodeInt==0)
			  PTResOneProng0Pi0_->Fill(PtRes);
			if(gendecaymodeInt==1)
			  PTResOneProng1Pi0_->Fill(PtRes);
			if(gendecaymodeInt==2)
			  PTResOneProng2Pi0_->Fill(PtRes);
			if(gendecaymodeInt==10)
			  PTResThreeProng0Pi0_->Fill(PtRes);
			if(gendecaymodeInt==11)
			  PTResThreeProng1Pi0_->Fill(PtRes);
			
			OneMinusPtRes_->Fill(1.0-(tauPt/genTauPt));
			GenMinusPatPt_->Fill(abs(genTauPt - tauPt));
			SigPlusIso_div_genJetPt_->Fill((tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum+tmpPFGPtSum+tmpPFCPtSum)/genTauPt);
			
		      }//tauItr->tauID("againstMuonTight")  )...{			    
		  }//if((tauItr->pt() > 20.0 && tauItr->pt() < 40.0))...{
		
		/*		    if ( tauItr->pfJetRef().isNonnull() )
		  {
		  PFJetPt=tauItr->pfJetRef()->pt();
		  //cout<<"pt of ref jet in fwlite ="<<tauItr->pfJetRef()->pt()<<endl;
		  PFJetPt_->Fill(tauItr->pfJetRef()->pt());
		  PFJetPx_->Fill(tauItr->pfJetRef()->px());
		  PFJetPy_->Fill(tauItr->pfJetRef()->py());
		  PFJetPz_->Fill(tauItr->pfJetRef()->pz());
		  PFJetPhi_->Fill(tauItr->pfJetRef()->phi());
		  PFJetEta_->Fill(tauItr->pfJetRef()->eta());
		  PFJetEnergy_->Fill(tauItr->pfJetRef()->energy());
		  }
		*/
		if (PFJetPt !=-999 && genTauPt !=-999)
		  {
		    checkbugptres_->Fill(PFJetPt);
		    checkgentaupt_->Fill(genTauPt);
		    checkpattaupt_->Fill(tauPt);
		    PTRESPFJET_->Fill(PFJetPt/genTauPt);
		    float binsize = 0.1;
		    for(int i=0; i<10;i++)
		      {
			if( PtRes >= i*binsize  &&  PtRes < (i+1)*binsize )
			  {
			    PTResPFJet_[i]->Fill((PFJetPt/genTauPt));
			    float delta = abs(tauPt-genTauPt);
			    SumPFN_vs_SumPFG_[i]->Fill((PFNisoPtSum)/delta,(PFGisoPtSum)/delta);
			    reco_vs_gen_decay_mode_in_bin_[i]->Fill(recotaudecaymode,gendecaymodeInt);
			  }
		      }
		  }
	      }//if(tauItr->genJet() )...{
	  }//for(std::vector<pat::Tau>::const_iterator tauItr=tauHandle->begin(); tauItr!=tauHandle->end(); ++tauItr)...{
	
	//	Book_Fill_histograms object1;
	//object1.bookHistograms();
	
      }//for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
      // close input file
      inFile->Close();
    }// if( inFile ){
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }// for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
  return 0;
}// int main(int argc, char* argv[])...{
