
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

void train_kNNmuonIsolation()
{
  gROOT->SetBatch(true);

  TString inputFilePath  = "/data1/veelken/tmp/muonIsoStudy/v4/";
  TString inputFileName  = "analyzeMuonIsoHistograms_PPmuXptGt20Mu15_2011Oct10v3.root";

  TString inputFileName_full = TString(inputFilePath).Append(inputFileName);
  TFile* inputFile = TFile::Open(inputFileName_full.Data());

  TString allTreeName = "HLT_IsoMu12_loose05_tight01/all/ntuple";
  TTree* allTree = dynamic_cast<TTree*>(inputFile->Get(allTreeName.Data()));
  assert(allTree);

  Float_t muonPt, muonEta, muonIso, sumEt, weight;

  allTree->SetBranchAddress("muonPt", &muonPt);
  allTree->SetBranchAddress("muonEta", &muonEta);
  allTree->SetBranchAddress("muonIso", &muonIso);
  //allTree->SetBranchAddress("sumEt", &sumEt);
  allTree->SetBranchAddress("weight", &weight);

  Float_t logMuonPt, absMuonEta, logSumEt;

  TTree* passTree = new TTree("passTree", "passTree");
  passTree->Branch("logMuonPt",  &logMuonPt,  "logMuonPt/F");
  passTree->Branch("absMuonEta", &absMuonEta, "absMuonEta/F");
  //passTree->Branch("logSumEt",   &logSumEt,   "logSumEt/F");
  passTree->Branch("weight",     &weight,     "weight/F");

  TTree* failTree = new TTree("failTree", "failTree");
  failTree->Branch("logMuonPt",  &logMuonPt,  "logMuonPt/F");
  failTree->Branch("absMuonEta", &absMuonEta, "absMuonEta/F");
  //failTree->Branch("logSumEt",   &logSumEt,   "logSumEt/F");
  failTree->Branch("weight",     &weight,     "weight/F");

  int numEntries = allTree->GetEntries();
  for ( int iEntry = 0; iEntry < numEntries; ++iEntry ) {
    allTree->GetEntry(iEntry);
    
    logMuonPt = ( muonPt > 1. ) ? TMath::Log(muonPt) : 0.;
    absMuonEta = TMath::Abs(muonEta);
    //logSumEt = ( sumEt > 1. ) ? TMath::Log(sumEt) : 0.;

    if ( muonIso < (0.10*muonPt) ) passTree->Fill();
    else if ( muonIso > (0.20*muonPt) ) failTree->Fill();
  }

  TString outputFileName = "train_kNNmuonIsolation_plots.root"; // for TMVA control plots

  TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
  
  TMVA::Tools::Instance();
  TMVA::Factory* factory = new TMVA::Factory("train_kNNmuonIsolation", outputFile, "!V:!Silent:Color:DrawProgressBar");
  factory->AddVariable("logMuonPt",  "P_{T}^{#mu}",  "GeV/c", 'F');
  factory->AddVariable("absMuonEta", "#eta_{#mu}",   "",      'F');
  //factory->AddVariable("logSumEt",   "#Sigma E_{T}", "GeV/c", 'F');
  
  factory->SetWeightExpression("weight");

  factory->AddSignalTree(passTree);
  factory->AddBackgroundTree(failTree);

  TCut cutS = "";
  TCut cutB = "";

  factory->PrepareTrainingAndTestTree(cutS, cutB, 
    "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=1:nTest_Background=1:SplitMode=Random:NormMode=NumEvents:!V");
  factory->BookMethod(TMVA::Types::kKNN, "kNN", "nkNN=20:UseWeight=False");
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();  

  outputFile->Close();
  
  delete factory;
}

