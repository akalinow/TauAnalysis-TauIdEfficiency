//#include <Wrapper.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include <vector>
#include <string>
#include <iostream>
//#include "DataFormats/Common/interface/MergeableCounter.h"

using namespace std;

int extractNevents(const char *filename)
{
  cout<< "<extractNevents>" << endl;
  TFile *_file0 = TFile::Open(filename);
  TTree *lumiblock = (TTree*) _file0->Get("LuminosityBlocks");
  //edm::Wrapper<edm::MergeableCounter> counter;
  //edm::MergeableCounter counter;
  //lumiblock->SetBranchAddress("edmMergeableCounter_totalEventsProcessed__skimTauIdEffSample.",&counter);

  int numbins = 10000000;
  TH1I *pippo = new TH1I("pippo","pippo",numbins,0,numbins);
  lumiblock->Draw("edmMergeableCounter_totalEventsProcessed__skimTauIdEffSample.obj.value>>pippo");
  while(pippo->GetBinContent(numbins+1) != 0){
    cout << "Overflow Bin Filled! Increasing size and repeating" << endl;
    numbins = numbins*2;
    delete pippo;
    pippo = new TH1I("pippo","pippo",numbins,0,numbins);
    lumiblock->Draw("edmMergeableCounter_totalEventsProcessed__skimTauIdEffSample.obj.value>>pippo");
  }
  double evtsProcessed =0;
  for(int i=1; i<=numbins; i++){
    //    if(i%1000 == 0) cout << "Analyzing Bin " << i << endl;
    evtsProcessed += pippo->GetBinContent(i)*(i-1);
  }
  /*for(int i=0; i<lumiblock->GetEntries(); i++){
    lumiblock->GetEntry(i);
    evtsProcessed += counter.value;
    }*/

  //cout << "Total events processed: " << evtsProcessed << endl;
  delete pippo;
  _file0->Close();
  return evtsProcessed;
}

void GetNumEvts()
{
  vector<string> sampleZtautau;
  vector<string> sampleZmumu;
  vector<string> sampleWplusJets;
  vector<string> sampleTTplusJets;
  vector<string> sampleQCD;

  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_0_ed96.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_10_f966.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_11_b5f9.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_12_0638.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_13_2bd2.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_14_1b17.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_15_1fa3.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_16_cc2b.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_17_66bf.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_18_e49c.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_19_331a.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_1_acad.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_20_c912.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_21_0883.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_22_ac76.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_23_4165.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_24_6c69.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_25_75cd.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_26_a935.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_27_d74e.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_28_5443.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_29_764f.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_2_d8ae.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_30_31f5.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_31_895a.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_32_6a7b.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_33_27cb.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_34_3e7f.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_35_454d.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_36_6fe4.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_37_e79d.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_38_b4f3.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_39_e9b3.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_3_b9e4.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_40_aa46.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_4_3d8b.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_5_64a8.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_6_3263.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_7_6fd2.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_8_3fab.root");
  sampleQCD.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_PPmuXptGt20Mu15_2011Jun06V1_9_2c94.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_0_5a68.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_10_c4b4.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_1_bbf6.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_2_9d14.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_3_4731.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_4_3cf9.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_5_3cb3.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_6_9c60.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_7_f284.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_8_2aea.root");
  sampleTTplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_TTplusJets_madgraph_2011Jun06V1_9_5a51.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_0_1127.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_10_ca4f.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_11_0a10.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_12_8690.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_13_bdd6.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_14_fc06.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_15_1380.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_16_0d7c.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_17_893e.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_18_93c3.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_1_da7e.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_2_2b63.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_3_ef53.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_4_7007.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_5_f82f.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_6_7bf8.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_7_4e02.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_8_b256.root");
  sampleWplusJets.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_WplusJets_madgraph_2011Jun06V1_9_3204.root");
  sampleZmumu.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Zmumu_powheg_2011Jun06V1_0_f76a.root");
  sampleZmumu.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Zmumu_powheg_2011Jun06V1_1_2dd5.root");
  sampleZmumu.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Zmumu_powheg_2011Jun06V1_2_e27d.root");
  sampleZtautau.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Ztautau_powheg_2011Jun06V1_0_ba35.root");
  sampleZtautau.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Ztautau_powheg_2011Jun06V1_1_3576.root");
  sampleZtautau.push_back("/nfs/data4/verzetti/tagprobe/NTuples/NTuplesJun06/tauIdEffMeasEDNtuple_Ztautau_powheg_2011Jun06V1_2_8965.root");

  double totevts = 0;
  for(vector<string>::const_iterator pippo = sampleZtautau.begin(); pippo != sampleZtautau.end(); pippo++){
    double partEvents = extractNevents(pippo->c_str());
    totevts += partEvents;
    cout << pippo->c_str() << " Events: " << partEvents <<endl;
  }
  cout<< "sampleZtautau: "<< totevts <<endl;


  totevts = 0;
  for(vector<string>::const_iterator pippo = sampleZmumu.begin(); pippo != sampleZmumu.end(); pippo++) totevts += extractNevents(pippo->c_str());
  cout<< "sampleZmumu: "<< totevts <<endl;


  totevts = 0;
  for(vector<string>::const_iterator pippo = sampleWplusJets.begin(); pippo != sampleWplusJets.end(); pippo++) totevts += extractNevents(pippo->c_str());
  cout<< "sampleWplusJets: "<< totevts <<endl;


  totevts = 0;
  for(vector<string>::const_iterator pippo = sampleTTplusJets.begin(); pippo != sampleTTplusJets.end(); pippo++) totevts += extractNevents(pippo->c_str());
  cout<< "sampleTTplusJets: "<< totevts <<endl;

  totevts = 0;
  for(vector<string>::const_iterator pippo = sampleQCD.begin(); pippo != sampleQCD.end(); pippo++) totevts += extractNevents(pippo->c_str());
  cout<< "sampleQCD: "<< totevts <<endl;
}
