#include "TProfile.h"
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
#include <TH1I.h>
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
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
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");
  // set defaults
  parser.integerValue ("maxEvents"  ) = -1;
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
  reco_vs_gen_decay_mode_ =  dir.make<TH2F>("reco_vs_gen_taudecaymode"  , "reco_vs_gen_taudecaymode;decay mode_{generated};decay mode_{reconstructed}"  ,   15,   0.0 , 15.0 ,   15,   0.0 , 15.0);
  
  TH1I* isoparticleID_;
  isoparticleID_  = dir.make<TH1I>("isoparticleID"  , "isoparticleID"  ,   10,   -0.5, 9.5);

  TH1F* PFJetPt_;
  PFJetPt_  = dir.make<TH1F>("PFJetPt"  , "PFJetPt"  ,   900,   0.0, 300.0);
  
  TH1F* PFJetPhi_;
  PFJetPhi_ = dir.make<TH1F>("PFJetPhi"  , "PFJetPhi"  ,   100,   -5.0, 5.0);
  
  TH1F* PFJetEta_;
  PFJetEta_  = dir.make<TH1F>("PFJetEta"  , "PFJetEta"  ,   100,   -3.0, 3.0);
  
  TH1F* PtRes_;
  PtRes_  = dir.make<TH1F>("PtRes" , "PtRes"  ,   100,   0.0, 5.0);

  TH1F* PtResNew_;
  PtResNew_  = dir.make<TH1F>("PtResNew" , "PtResNew"  ,   100,   0.0, 5.0);
  
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
      reco_vs_gen_decay_mode_in_bin_[ibin]  = dir.make<TH2F>(("reco_vs_gen_decay_mode_in_bin_"+histname).c_str()  , "reco_vs_gen_decay_mode;decaymode_{gen};decaymode_{reco}" , 15, 0.0, 15.0,15,0,15);
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
  
  TProfile*  PtRatio_vs_ParticleID_ ;
  PtRatio_vs_ParticleID_ = dir.make<TProfile>("SumPtRatio_vs_ParticleID_inbins_" , "ParticleId vs Pt_{Particle}/Pt_{GenTauJet};ParticleId;SumPt_{Particle in iso}/Pt_{GenTauJet}" ,10,0.0,10, 0.0, 10.0);
    
    TProfile* PtRatio_vs_ParticleID_inbins_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      PtRatio_vs_ParticleID_inbins_[ibin] = dir.make<TProfile>(("SumPtRatio_vs_ParticleID_inbins_"+histname).c_str()  , "ParticleId vs Pt_{Particle}/Pt_{GenTauJet};ParticleId;SumPt_{Particle in iso}/Pt_{GenTauJet}" ,10,0.0,10, 0.0, 10.0);
    }
  
  TH2F* RecoGammaSum_vs_GenGammaSum_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      RecoGammaSum_vs_GenGammaSum_[ibin]  = dir.make<TH2F>(("RecoGammaSum_vs_GenGammaSum_"+histname).c_str()  , "SumPt_{gamma reco} vs SumPt_{gamma gen};SumPt_{gamma reco};SumPt_{gamma gen}" , 2000, 0.0, 75.0,2000, 0, 75);
    }  

  TH2F* RecoGammaSum_tot_vs_GenGammaSum_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      RecoGammaSum_tot_vs_GenGammaSum_[ibin]  = dir.make<TH2F>(("RecoGammaSum_tot_vs_GenGammaSum_"+histname).c_str()  , "SumPt_{gamma reco} vs SumPt_{gamma gen};SumPt_{gamma reco};SumPt_{gamma gen}" , 2000, 0.0, 75.0,2000, 0, 75);
    }  
  
  TH1F* Ratio_RecSumPtCgdHad_GenJetCgdHad_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      Ratio_RecSumPtCgdHad_GenJetCgdHad_[ibin]  = dir.make<TH1F>(("Ratio_RecSumPtCgdHad_GenJetCgdHad_"+histname).c_str()  , "SumPt_{chgd had} / SumPt_{gen Jet chgd had};SumPt_{reco chgd had}/SumPt_{genJet chgd had};# of events" , 500, 0, 5.0);
    }  
  

  TH1F* Ratio_RecSumPtGam_GenJetSumPtGam_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      Ratio_RecSumPtGam_GenJetSumPtGam_[ibin]  = dir.make<TH1F>(("Ratio_RecSumPtGam_GenJetSumPtGam_"+histname).c_str()  , "SumPt_{#gamma reco} / SumPt_{#gamma gen};SumPt_{#gamma reco}/SumPt_{#gamma gen};# of events" , 500, 0, 5.0);
    }  

  TH1F* Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_[ibin]  = dir.make<TH1F>(("Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_"+histname).c_str()  , "(SumPt_{#gamma gen} - SumPt_{#gamma reco})/Pt_{Gen Tau Jet} ; (SumPt_{#gamma gen}-SumPt_{#gamma reco})/Pt_{Gen Tau Jet};# of events" , 500, -1.0, 1.0);
    }  

  TH1F* Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_[10];
  for(int ibin=0; ibin <10; ibin++)
    {
      char tmphistname[100];
      sprintf(tmphistname,"_%d_to_%d",ibin,ibin+1);
      std::string histname(tmphistname);
      Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_[ibin]  = dir.make<TH1F>(("Ratio_GenJetPtSumChsHad_Minus_GenJetSumPtChdHad_"+histname).c_str()  , "(SumPt_{Chd Had gen} - SumPt_{Chd Had reco})/Pt_{Gen Tau Jet} ; (SumPt_{Chd Had gen}-SumPt_{Chd Had reco})/Pt_{Gen Tau Jet};# of events" , 500, -0.5, 0.5);
    }  


  TH1F* Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[5];
  std::string decaymode[5]={"OneP0Pi0", "OneP1Pi0", "OneP2Pi0", "ThreeP0Pi0", "ThreeP1Pi0"};
  for(int ibin=0; ibin <5; ibin++)
    {
      std::string histname(decaymode[ibin]);
      Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[ibin]  = dir.make<TH1F>(("Ratio_GenJetPtSumChsHad_Minus_GenJetSumPtChdHad_"+histname).c_str()  , "(SumPt_{Chd Had gen} - SumPt_{Chd Had reco})/Pt_{Gen Tau Jet} ; (SumPt_{Chd Had gen}-SumPt_{Chd Had reco})/Pt_{Gen Tau Jet};# of events" , 500, -0.5, 0.5);
    }  
  
  
  TH1F* Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[5];
  //std::string decaymode[5]={"OneP0Pi0", "OneP1Pi0", "OneP2Pi0", "ThreeP0Pi0", "ThreeP1Pi0"};
  for(int ibin=0; ibin <5; ibin++)
    {
      std::string histname(decaymode[ibin]);
      Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[ibin]  = dir.make<TH1F>(("Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_"+histname).c_str()  , "(SumPt_{#gamma gen} - SumPt_{#gamma reco})/Pt_{Gen Tau Jet} ; (SumPt_{#gamma gen}-SumPt_{#gamma reco})/Pt_{Gen Tau Jet};# of events" , 500, -1.0, 1.0);
    }  

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
      Int_t itau                                                        ;
      float PionSumPt          , PhotonSumPt                            ;
      float tmpPFCPtSum        , tmpPFGPtSum         , tmpPFNPtSum      ;
      float tmpPFCisoPtSum     , tmpPFGisoPtSum       , tmpPFNisoPtSum  ;
      float genTauPt           , tauPt                                  ;
      float PFCisoPtSum        , PFGisoPtSum         , PFNisoPtSum      ;
      float PtRes              , PFJetPt                                ;
      int   recotaudecaymode   , gendecaymodeInt =-999                  ;
      
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
	/*edm::Handle<std::vector<pat::Tau> > patTausAK5GenJetMatchedHandle;
	  event.getByLabel(std::string("patTausAK5GenJetMatched"), patTausAK5GenJetMatchedHandle);
	  for(std::vector<pat::Tau>::const_iterator tauAK5GenJetMatchedItr=patTausAK5GenJetMatchedHandle->begin(); tauAK5GenJetMatchedItr!=patTausAK5GenJetMatchedHandle->end(); ++tauAK5GenJetMatchedItr)
	  {
	  if( tauAK5GenJetMatchedItr->pt() > 20.0 &&  tauAK5GenJetMatchedItr->pt() < 40.0 )      {
	  if(tauAK5GenJetMatchedItr->genJet() ) {
	  if ( tauAK5GenJetMatchedItr->pfJetRef().isNonnull() ) {
	  std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(*tauAK5GenJetMatchedItr->genJet());		  
	  cout<<"Matched Jet Pt = "<<tauAK5GenJetMatchedItr->pfJetRef()->pt()<<endl;
	  cout<<"new decay mode is =  "<<genTauDecayMode<<endl;
	  
	  //}
	  }
	  }
	  }
	  }
	*/

	//////////////////////////////////////////////////////////////
	//////////////////////////Clean Pat Taus//////////////////////
	//////////////////////////////////////////////////////////////
	edm::Handle<std::vector<pat::Tau> >  tauHandle;
	event.getByLabel(std::string("cleanPatTaus"),tauHandle);
	itau =0;
	
	for(std::vector<pat::Tau>::const_iterator tauItr=tauHandle->begin(); tauItr!=tauHandle->end(); ++tauItr)
	  {
	    tauPt =-999;
	    genTauPt =-999;
	    gendecaymodeInt =-999;
	    PFJetPt = -999;
	    PtRes = -999;
	    recotaudecaymode =-999;
	    //Studying only those taus which have 20.0GeV  < Pt < 40GeV , there should be matched gen Jet and PFJet corresponding to the HPS Tau under study.
	    if((tauItr->pt() > 20.0 && tauItr->pt() < 40.0)) {    
	      if(tauItr->genJet() )     {
		if ( tauItr->pfJetRef().isNonnull() ) {
		  std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(*tauItr->genJet());
		  //check that tau decay mode should be hadronic //
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
		    // store gen tau decay mode in terms of integer, enum can be used (may be enum is easier to use)
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
		    //if gen jet constituents are available then loop over them
		    if( (tauItr->genJet()->getGenConstituents()).size() !=0 )
		      {
			std::vector <const GenParticle*>  genParticleList = tauItr->genJet()->getGenConstituents();
			PionSumPt=0.0;
			PhotonSumPt=0.0;
			//loop over gen jet constituents
			for (int iList = 0; iList <int(genParticleList.size()); iList++)
			  {
			    //look for charged hadrons, check if there is sme mistake in pdgId, summing Pt of all charged hadrons.
			    // if(abs(genParticleList[iList]->pdgId()) == 211 || abs(genParticleList[iList]->pdgId()) == 321)
			    if(abs(genParticleList[iList]->charge()) > 0.5 )
			      PionSumPt+= genParticleList[iList]->pt();
			    
			    //look for photons coming from pi0, summing Pt of all Photons
			    if(abs(genParticleList[iList]->pdgId()) == 22)
			      PhotonSumPt+=genParticleList[iList]->pt();
			    
			  }//for (int iList = 0; iList <int(genParticleList.size()); iList++)...{
			
			//fill Pt Sum of charged hadron &  Pt Sum of Photons coming from pi0
			PionSumPt_->Fill(PionSumPt);  
			PhotonSumPt_->Fill(PhotonSumPt);
		      }//if( (tauItr->genJet()->getGenConstituents()).size() !=0 )...{
		    
		    //////////////////////////////////////////////////////
		    ///////////////PFCandidates in signal region//////////
		    /////////////////////////////////////////////////////
		    tmpPFCPtSum = 0.0;
		    tmpPFGPtSum = 0.0;
		    tmpPFNPtSum  = 0.0;
		    
		    // charged candidates in signal region
		    const reco::PFCandidateRefVector &chargedref = tauItr->signalPFChargedHadrCands() ;
		    if( chargedref.size() > 0 ) {
		      for(int i=0; i<(int)chargedref.size(); i++)
			tmpPFCPtSum = tmpPFCPtSum + chargedref[i]->pt();
		      PFCPtSum_->Fill(tmpPFCPtSum);
		    }
		    //Gamma chandidates in Signal region
		    const reco::PFCandidateRefVector &gammaref = tauItr->signalPFGammaCands();
		    if ( gammaref.size() > 0) {
		      for(int i=0; i<(int)gammaref.size(); i++)
			tmpPFGPtSum = tmpPFGPtSum + gammaref[i]->pt();
		      PFGPtSum_->Fill(tmpPFGPtSum);
		    }
		    //Neutral hadrons in signal region, ideally it should be ZERO, and has been checked that it is zero.
		    const reco::PFCandidateRefVector &neutralref = tauItr->signalPFNeutrHadrCands();
		    if ( neutralref.size() > 0 ) {
		      for(int i=0; i<(int)neutralref.size(); i++)
			tmpPFNPtSum = tmpPFNPtSum + neutralref[i]->pt();
		      PFNPtSum_->Fill(tmpPFNPtSum);
		    }
		    //Sum of Gamma and Charged hadrons in Signal Region
		    PFCandidateSumPt_->Fill(tmpPFGPtSum+tmpPFCPtSum);
		    //Ratio of (Sum of Gamma and Charged hadrons in Signal Region) and GenTauJetPt , Ideally it should be
		    //equal to PtRes of Reco Jet, and it has been checked. 
		    PFCandidateSumPt_div_genJetPt_->Fill((tmpPFGPtSum+tmpPFCPtSum)/genTauPt);
		    
		    ///////////////////////////////////////////////////////////////
		    /////////////PFCandidates in isolation region//////////////////
		    ///////////////////////////////////////////////////////////////
		    tmpPFCisoPtSum = 0.0;
		    tmpPFGisoPtSum = 0.0;
		    tmpPFNisoPtSum = 0.0;
		    //Charged hadrons in Isolation Region
		    const reco::PFCandidateRefVector &chargedisoref = tauItr->isolationPFChargedHadrCands() ;
		    for(int i=0; i<(int)chargedisoref.size(); i++)
		      tmpPFCisoPtSum = tmpPFCisoPtSum + chargedisoref[i]->pt() ;
		    PFCisoPtSum = tmpPFCisoPtSum;
		    PFCisoPtSum_->Fill(tmpPFCisoPtSum);
		    //Charged hadrons in Isolation Region
		    const reco::PFCandidateRefVector &gammaisoref = tauItr->isolationPFGammaCands();
		    for(int i=0; i<(int)gammaisoref.size(); i++)
		      tmpPFGisoPtSum = tmpPFGisoPtSum + gammaisoref[i]->pt() ;
		    PFGisoPtSum = tmpPFGisoPtSum;
		    PFGisoPtSum_->Fill(tmpPFGisoPtSum);
		    //Charged hadrons in Isolation Region
		    const reco::PFCandidateRefVector &neutralisoref = tauItr->isolationPFNeutrHadrCands();
		    for(int i=0; i<(int)neutralisoref.size(); i++)
		      tmpPFNisoPtSum = tmpPFNisoPtSum + neutralisoref[i]->pt() ; 
		    PFNisoPtSum = tmpPFNisoPtSum;
		    PFNisoPtSum_->Fill(tmpPFNisoPtSum);
		    //Sum of (Charged hadrons, Neutral hadrons, Photons) in Isolation Region
		    PFCandidateisoSumPt_->Fill(tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum);
		    //Ratio of ( Sum of (Charged hadrons, Neutral hadrons, Photons) in Isolation Region ) and GenTauJetPt
		    PFCandidateisoSumPtDivGenTauPt_->Fill((tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum)/genTauPt);
		    /////////////////////////////////////////////////////////////////// 
		    ///////////////Particle Flow Jets////////////////////////////////// 
		    /////////////////////////////////////////////////////////////////// 
		    PFJetPt=tauItr->pfJetRef()->pt();
		    PFJetPt_->Fill(tauItr->pfJetRef()->pt());
		    PFJetPhi_->Fill(tauItr->pfJetRef()->phi());
		    PFJetEta_->Fill(tauItr->pfJetRef()->eta());
		    
		    std::vector <reco::PFCandidatePtr> pfJetConstituents = tauItr->pfJetRef()->getPFConstituents();
		    float PFJetConst_gammaSumPt = 0;
		    for ( int iParticle = 0; iParticle < int(pfJetConstituents.size()); iParticle++)
		      {
			if(pfJetConstituents[iParticle]->particleId()==reco::PFCandidate::gamma)  //particle id 4 is reco::PFCandidate::gamma
			  {
			    PFJetConst_gammaSumPt += pfJetConstituents[iParticle]->pt() ;
			  }
		      }
		    PtResNew_->Fill( (PFJetPt - tmpPFGPtSum + PFJetConst_gammaSumPt)/genTauPt );
		    
		    //for debuging
		    //cout<<"photon is found with Pt "<<pfJetConstituents[iParticle]->pt()<<endl;
		    //cout<<"particle number = "<<iParticle<<"  with particle id = "<<pfJetConstituents[iParticle]->particleId()<<" and   particle Pt  = "<<pfJetConstituents[iParticle]->pt()<<endl;
		    //cout<<" Photon Sum Pt in PFJet = "<<PFJetConst_gammaSumPt<<endl;
		    //cout<<" particleId() = "<<pfJetConstituents.size()<<endl;
		    /*cout<<" PfJet Charged hadron Pt = "<<tauItr->pfJetRef()->chargedHadronEnergy()<<endl;
		      cout<<" PfJet Neutral hadron Pt = "<<tauItr->pfJetRef()->neutralHadronEnergy()<<endl;
		    cout<<" PfJet Gamma          Pt = "<<tauItr->pfJetRef()->photonEnergy()<<endl;
		    cout<<" PfJet Total          Pt = "<<( tauItr->pfJetRef()->photonEnergy() + tauItr->pfJetRef()->neutralHadronEnergy()  +  tauItr->pfJetRef()->chargedHadronEnergy() )<<endl;
		    */

		    
		    /////////////////////////////////////////////////////////////////// 
		    ///////////////////////Pat Tau Jet///////////////////////////////// 
		    /////////////////////////////////////////////////////////////////// 
		    tauPt = tauItr->pt();
		    tauPt_->Fill(tauItr->pt());
		    tauEta_->Fill(tauItr->eta());
		    tauPhi_->Fill(tauItr->phi());
		    recotaudecaymode = tauItr->decayMode();
		    recotaudecaymode_->Fill(tauItr->decayMode());
		    reco_vs_gen_decay_mode_->Fill(gendecaymodeInt , recotaudecaymode);
		    // fill PtRes when gen decay mode and reco decay mode are same
		    if(recotaudecaymode == gendecaymodeInt)
		      matched_PtRes_->Fill(tauPt/genTauPt);
		    // fill PtRes for all selected taus which are matched to gen tau Jet in eta, phi
		    PtRes = tauPt/genTauPt;
		    PtRes_->Fill(tauPt/genTauPt);
		    //fill PtRes according to gen tau decay mode
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
		    // These three histo are for checking bugs, keeping them alive in code for some time. 
		    OneMinusPtRes_->Fill(1.0-(tauPt/genTauPt));
		    GenMinusPatPt_->Fill(abs(genTauPt - tauPt));
		    SigPlusIso_div_genJetPt_->Fill((tmpPFNisoPtSum + tmpPFGisoPtSum + tmpPFCisoPtSum+tmpPFGPtSum+tmpPFCPtSum)/genTauPt);
		    //requiring this, if condition is not satisfied then default value is -999 will be filled in histo. 
		    //check if this is required or not .. ?
		    if (PFJetPt !=-999 && genTauPt !=-999 && tauPt != -999)
		      {
			PTRESPFJET_->Fill(PFJetPt/genTauPt);
			float binsize = 0.1;
			for(int i=0; i<10;i++)
			  {
			    if( PtRes >= i*binsize  &&  PtRes < (i+1)*binsize )
			      {
				PTResPFJet_[i]->Fill((PFJetPt/genTauPt));
				float delta = abs(tauPt-genTauPt);
				//sum of PF Gamma vs Sum of PF Neutral hadron sum with proper normalisation
				//this variable may be not needed later on.
				SumPFN_vs_SumPFG_[i]->Fill((PFNisoPtSum)/delta,(PFGisoPtSum)/delta);
				//fill (reco tau decay mode) vs (gen tau decay mode) in bins of PtRes
				reco_vs_gen_decay_mode_in_bin_[i]->Fill(gendecaymodeInt , recotaudecaymode);
				//Reco Gamma Sum Pt in signal vs Gamma Sum Pt in GenJet  matched to tau.
				RecoGammaSum_vs_GenGammaSum_[i]->Fill(PhotonSumPt,tmpPFGPtSum);
				//Reco Gamma Sum Pt in (signal + isolation list)  vs Gamma Sum Pt in GenJet  matched to tau.
				RecoGammaSum_tot_vs_GenGammaSum_[i]->Fill(PhotonSumPt, ( tmpPFGPtSum + tmpPFGisoPtSum ) );
				//(sumPt of patTau signalPFChargedHadrCands()) /  (sum Pt of patTau->genJet() constituents with abs(charge) > 0.5 )
				Ratio_RecSumPtCgdHad_GenJetCgdHad_[i]->Fill(tmpPFCPtSum/PionSumPt);
				//(sumPt of patTau signalPFGammaCands()) /  (sum Pt of patTau->genJet() constituents with pdgId == 22 )
				Ratio_RecSumPtGam_GenJetSumPtGam_[i]->Fill(tmpPFGPtSum/PhotonSumPt);
				//( (sum Pt of patTau->genJet() constituents with pdgId == 22 ) - (sumPt of patTau signalPFGammaCands()) ) / genTauPt
				Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_[i]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
				// ( (sum Pt of patTau->genJet() constituents with abs(charge) > 0.5 ) - (sumPt of patTau signalPFChargedHadrCands()) ) / genTauPt
				Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_[i]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			      }
			  }
			if(gendecaymodeInt==0){
			  Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[0]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			  Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[0]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
			}
			if(gendecaymodeInt==1) {
			  Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[1]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			  Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[1]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
			}
			if(gendecaymodeInt==2) {
			  Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[2]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			  Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[2]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
			}
			if(gendecaymodeInt==10) {
			  Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[3]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			  Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[3]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
			}
			if(gendecaymodeInt==11) {
			  Ratio_GenJetPtSumChdHad_Minus_GenJetSumPtChdHad_decaymode_[4]->Fill( (PionSumPt - tmpPFCPtSum)/genTauPt ) ;
			  Ratio_GenJetPtSumGam_Minus_GenJetSumPtGam_decaymode_[4]->Fill( (PhotonSumPt - tmpPFGPtSum)/genTauPt ) ;
			}
		      }
		    
		    
		    //requiring this, if condition is not satisfied then default value is -999 will be filled in histo. 
		    //check if this is required or not .. ?
		    //this if statement is for filling PtRatio_vs_ParticleID in bins of PtRes.
		    if (PFJetPt !=-999 && genTauPt !=-999 && tauPt != -999)
		      {
			float binsize = 0.1;
			const PFCandidateRefVector &pfcandidateisoref = tauItr->isolationPFCands();
			std::vector<double> sumPFCandidatePt(8);
			for (int i = 0; i<int(pfcandidateisoref.size()) ; i++)
			  sumPFCandidatePt[pfcandidateisoref[i]->particleId()] += pfcandidateisoref[i]->pt();
			for ( unsigned i = 0; i < 8; ++i )
			  {
			    PtRatio_vs_ParticleID_->Fill(i , sumPFCandidatePt[i]/genTauPt);
			    for(int j=0; j<10;j++)
			      {
				if( PtRes >= j*binsize  &&  PtRes < (j+1)*binsize )
				  {
				    PtRatio_vs_ParticleID_inbins_[j]->Fill(i , sumPFCandidatePt[i]/ genTauPt);
				  }
			      }
			  }
		      }//if (PFJetPt !=-999 && genTauPt !=-999 && tauPt != -999)
		    }//genTauDecayMode == "threeProng1Pi0" ))...{
		}//if ( tauItr->pfJetRef().isNonnull() ) {
	      }//if(tauItr->genJet() )...{	
	    }//if((tauItr->pt() > 20.0 && tauItr->pt() < 40.0))...{
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
