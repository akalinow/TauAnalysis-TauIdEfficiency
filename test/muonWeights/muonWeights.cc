#include "muonWeights.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TTree.h"

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::ChannelSpecifics(){

  initializeLeptonCorrections();
  initializePileUpCorrections();
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::~ChannelSpecifics(){


}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializePileUpCorrections(){

  TFile::SetCacheFileDir("/tmp/");
  std::string dataPUFileName = "http://akalinow.web.cern.ch/akalinow/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_PileUp.root";

  TFile *puDataFile_ = TFile::Open(dataPUFileName.c_str(),"CACHEREAD");

   std::string hName = "pileup";
   hPUWeight = (TH1F*)puDataFile_->Get(hName.c_str());

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializeLeptonCorrections(){

  TFile::SetCacheFileDir("/tmp/");
  std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v16_5.root";
  TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");

  RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");
  muon_id_iso_scalefactor = scaleWorkspace->function("m_idiso0p15_desy_ratio");
  muon_trg_scalefactor = scaleWorkspace->function("m_trgMu22OR_eta2p1_desy_ratio");
  muon_trk_scalefactor = scaleWorkspace->function("m_trk_ratio");

  m_eta = (RooRealVar*)scaleWorkspace->var("m_eta");
  m_pt = (RooRealVar*)scaleWorkspace->var("m_pt");
  m_iso = (RooRealVar*)scaleWorkspace->var("m_iso");

  std::cout<<"m_eta->getVal(): "<<m_eta->getVal()<<std::endl;
  std::cout<<"muon_id_iso_scalefactor->getVal(): "<<muon_id_iso_scalefactor->getVal()<<std::endl;
  m_eta->setVal(1.4);
  std::cout<<"m_eta->getVal(): "<<m_eta->getVal()<<std::endl;
  std::cout<<"muon_id_iso_scalefactor->getVal(): "<<muon_id_iso_scalefactor->getVal()<<std::endl;
    

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getLeptonCorrection(float pt, float eta, float iso){

  m_eta->setVal(eta);
  m_pt->setVal(pt);
  m_iso->setVal(iso);
  
  float id_iso = muon_id_iso_scalefactor->getVal();
  float trg = muon_trg_scalefactor->getVal();
  float trk = muon_trk_scalefactor->getVal();

  //std::cout<<"pt: "<<pt<<" iso SF: "<<id_iso<<" "<<trg<<" "<<trk<<std::endl;
  
  return id_iso*trg*trk;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getPUWeight(float nPU){

  int iBinPU = hPUWeight->FindBin(nPU);
  return hPUWeight->GetBinContent(iBinPU);
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void  ChannelSpecifics::addDataWeightBranch(std::string path) {

  std::string fileName = path+"/tnpZ_Data.root";
  TFile *file = TFile::Open(fileName.c_str());
  TTree *tree = (TTree*)file->Get("tpTree/fitter_tree");

  Float_t weight = 1.0;
  
  TFile *fOut = new TFile("tnpZ_DatawithWeights.root", "RECREATE");
  fOut->mkdir("tpTree")->cd();
  TTree *tOut = tree->CloneTree(0);
  tOut->Branch("weight", &weight, "weight/F");
  
  Long64_t nentries = tree->GetEntries();
  
  for (Long64_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    tOut->Fill(); 
  }
  
  tOut->AutoSave();
  fOut->Close();
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void  ChannelSpecifics::addMCWeightBranch(std::string path) {

  std::string fileName = path+"/tnpZ_MC.root";
  
  TFile *file = TFile::Open(fileName.c_str());
  TTree *tree = (TTree*)file->Get("tpTree/fitter_tree");

  ////////Calculate the PU weights
  tree->Draw("pair_truePileUp>>hMC(1000,0,100)","pair_truePileUp>0 && pair_truePileUp<100");
  TH1F *hMC = (TH1F*)gDirectory->Get("hMC");

  hPUWeight->Scale(1.0/hPUWeight->Integral());
  hMC->Scale(1.0/hMC->Integral());
  hPUWeight->Divide(hMC);
  //////////////////////////////////

  Float_t weight, tagPt, tagEta, tagIsolation;
  Float_t nPU;

  TFile *fOut = new TFile("tnpZ_MCwithWeights.root", "RECREATE");
  fOut->mkdir("tpTree")->cd();
  TTree *tOut = tree->CloneTree(0);
  tOut->Branch("weight", &weight, "weight/F");

  tree->SetBranchAddress("tag_pt", &tagPt);
  tree->SetBranchAddress("tag_eta", &tagEta);
  tree->SetBranchAddress("tag_isolation", &tagIsolation);
  tree->SetBranchAddress("pair_truePileUp", &nPU);

  Long64_t nentries = tree->GetEntries();
  
  for (Long64_t i=0;i<nentries;i++) {
    tree->GetEntry(i);
    weight = getPUWeight(nPU);
    //weight = getPUWeight(nPU)*getLeptonCorrection(tagPt, tagEta, tagIsolation);
    //std::cout<<"weight: "<<weight<<std::endl;
    tOut->Fill(); 
  }

tOut->AutoSave();
fOut->Close();
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
