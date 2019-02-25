void addWeightBranch() {
TFile *_file0 = TFile::Open("/home/akalinow/scratch/CMS/TauID/Crab/Data/TauID_TnP_2017/v5_Mu2Tau_2017/tnpZ_MC.root"); 
TTree *T = (TTree*)_file0->Get("tpTree/fitter_tree"); 

Float_t weight, mass;
TFile *fOut = new TFile("tnpZ_MCwithWeights.root", "RECREATE");
fOut->mkdir("tpTree")->cd();
TTree *tOut = T->CloneTree(0);
tOut->Branch("weight", &weight, "weight/F");
T->SetBranchAddress("mass", &mass);

Long64_t nentries = T->GetEntries(); 
for (Long64_t i=0;i<nentries;i++) {
	T->GetEntry(i);
	if (mass > 100) weight = 1;
	else weight = 0; 
	tOut->Fill(); 
	}

tOut->AutoSave();

fOut->Close();

}
