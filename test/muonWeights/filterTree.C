void filterTree(const std::string & fileName){

  vector<std::string> aSelections = {
    "mass>60",
    "mass<120",
    "abseta<2.3",
    "alternatLorentzVectPt>20",
    "alternatLorentzVectEta>-2.3",
    "alternatLorentzVectEta<2.3",
    "tag_pt>25",
    "tag_triggerMatch>0.5",
    "tag_dB<0.004",
    "pair_dz>-0.01",
    "pair_dz<0.01",
    "pair_deltaR>0.5",
    "pair_deltaR<5",
    "pair_probeMultiplicity==1",
    "pair_BestZ>0.5",
    "pair_MET<2500",
    "pair_MTtag<40",
    "pair_MTprobe<4000",
    "decayModeFindingNewDMs>0.5",
    "byTightIsolationMVArun2v1DBnewDMwLT2017v2>0.5",    
  };

  TCut filterCut = "";
  for(auto aSelection : aSelections){
      std::cout<<aSelection<<std::endl;
      filterCut += aSelection.c_str();
    }

  filterCut.Print();

  TFile *file = TFile::Open(fileName.c_str());
  TTree *tree = (TTree*)file->Get("tpTree/fitter_tree");
  std::cout<<"Original tree size:"<<tree->GetEntries()<<std::endl;

  std::string newFileName = fileName;
  size_t pos = fileName.find(".root");
  size_t len = std::string::npos;
  newFileName.replace(pos,len, "_filtered.root");
    
  TFile newFile(newFileName.c_str(), "recreate");
  TDirectory *aDir = newFile.mkdir("tpTree");
  aDir->cd();
  TTree *filteredTree = tree->CopyTree(filterCut);
  std::cout<<"Filtered tree size:"<<filteredTree->GetEntries()<<std::endl;

  filteredTree->SetBranchStatus("*",0);
  for (auto activeBranchName : {"mass", "abseta", "againstMuonLoose3", "againstMuonTight3", "weight", "mcTrue"}){
    filteredTree->SetBranchStatus(activeBranchName,1);
  }
  
  TTree *filtered_slimmedTree = filteredTree->CloneTree(-1,"");
  std::cout<<"Slimmed tree size:"<<filtered_slimmedTree->GetEntries()<<std::endl;

  filteredTree->SetDirectory(0);
  newFile.Close();
  
}
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
void filterAll(){

  std::string path = "/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP_2017/Mu2Tau_2017_v9/";
  std::string fileName = path+"/tnpZ_MCwithWeights.root";

  fileName = path+"/tnpZ_MCwithWeights.root";
  filterTree(fileName);

  fileName = path+"/tnpZ_Data.root";
  filterTree(fileName);
  
}
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
