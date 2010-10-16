void makeTauIsolationPlots()
{
  TFile* inputFile = TFile::Open("../test/patTauIsolationAnalyzer.root");

  TObjArray ptThresholds;
  ptThresholds.Add(new TObjString("0_50GeV"));
  ptThresholds.Add(new TObjString("1_00GeV"));
  ptThresholds.Add(new TObjString("1_50GeV"));
  ptThresholds.Add(new TObjString("2_00GeV"));
  ptThresholds.Add(new TObjString("2_50GeV"));
  ptThresholds.Add(new TObjString("3_00GeV"));
  unsigned numPtThresholds = ptThresholds.GetEntries();

  TObjArray sigConeSizes;
  sigConeSizes.Add(new TObjString("0_05dRsig"));
  sigConeSizes.Add(new TObjString("0_10dRsig"));
  sigConeSizes.Add(new TObjString("0_15dRsig"));
  sigConeSizes.Add(new TObjString("0_20dRsig"));
  sigConeSizes.Add(new TObjString("0_25dRsig"));
  sigConeSizes.Add(new TObjString("0_30dRsig"));
  unsigned numSigConeSizes = sigConeSizes.GetEntries();

  TString dqmDirectory = "";

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);

//-------------------------------------------------------------------------------
// show distributions of isolation Pt sums **before** tight tau id. is applied
//-------------------------------------------------------------------------------
  TPostScript* psBeforeTauId = new TPostScript("patTauIsolationPlots_beforeTauId.ps", 112);
  for ( unsigned iSigConeSize = 0; iSigConeSize < numSigConeSizes; ++iSigConeSize ) {
    TObjString* sigConeSize = (TObjString*)sigConeSizes.At(iSigConeSize);
//--- show plots for PFCandidate Pt > XX GeV
//   (same threshold for all types of PFCandidates)
    for ( unsigned iPtThreshold = 0; iPtThreshold < numPtThresholds; ++iPtThreshold ) {
      TObjString* ptThreshold = (TObjString*)ptThresholds.At(iPtThreshold);

      TString meName = TString("PFCandIso").Append(sigConeSize->GetString()).Append(ptThreshold->GetString()).Append("_matched");

      showTauIsolation(inputFile, dqmDirectory, meName, 
		       "TauPFIsolationQuantities/beforeTauId", "MuonPFIsolationQuantities", "ElectronPFIsolationQuantities",
		       canvas, psBeforeTauId, "beforeTauId", true);
    }
//--- show plots for PFChargedHadron Pt > 1.0 GeV, PFGamma Pt > 1.5 GeV
    TString meName = TString("PFCandIso").Append(sigConeSize->GetString()).Append("1_00_1_50GeV").Append("_matched");
    showTauIsolation(inputFile, dqmDirectory, meName, 
		     "TauPFIsolationQuantities/beforeTauId", "MuonPFIsolationQuantities", "ElectronPFIsolationQuantities",
		     canvas, psBeforeTauId, "beforeTauId", true);
  }
  delete psBeforeTauId;

//-------------------------------------------------------------------------------
// show distributions of isolation Pt sums **after** tight tau id. is applied
//-------------------------------------------------------------------------------
  TPostScript* psAfterTauId = new TPostScript("patTauIsolationPlots_afterTauId.ps", 112);
  for ( unsigned iSigConeSize = 0; iSigConeSize < numSigConeSizes; ++iSigConeSize ) {
    TObjString* sigConeSize = (TObjString*)sigConeSizes.At(iSigConeSize);
//--- show plots for PFCandidate Pt > XX GeV
//   (same threshold for all types of PFCandidates)
    for ( unsigned iPtThreshold = 0; iPtThreshold < numPtThresholds; ++iPtThreshold ) {
      TObjString* ptThreshold = (TObjString*)ptThresholds.At(iPtThreshold);

      TString meName = TString("PFCandIso").Append(sigConeSize->GetString()).Append(ptThreshold->GetString()).Append("_matched");

      showTauIsolation(inputFile, dqmDirectory, meName, 
		       "TauPFIsolationQuantities/afterTauId", "", "",
		       canvas, psAfterTauId, "afterTauId", true);
    }
//--- show plots for PFChargedHadron Pt > 1.0 GeV, PFGamma Pt > 1.5 GeV
    TString meName = TString("PFCandIso").Append(sigConeSize->GetString()).Append("1_00_1_50GeV").Append("_matched");
    showTauIsolation(inputFile, dqmDirectory, meName, 
		     "TauPFIsolationQuantities/afterTauId", "", "",
		     canvas, psAfterTauId, "afterTauId", true);
  }
  delete psAfterTauId;

  delete canvas;

  delete inputFile;
}

TH1* getMonitorElement(TFile* inputFile, const TString& dqmDirectory, const char* dqmSubDirectory, const TString& meName)
{
  TString meName_full = TString("DQMData").Append("/");
  if ( dqmDirectory != "") meName_full.Append(dqmDirectory).Append("/");
  meName_full.Append(dqmSubDirectory).Append("/").Append(meName);
  std::cout << "meName_full = " << meName_full << std::endl;
  
  TH1* me = (TH1*)inputFile->Get(meName_full);
  std::cout << "me = " << me <<  std::endl;
  
  //if ( !me->GetSumw2() ) me->Sumw2();
  me->Sumw2();

  me->Rebin(2);

  me->Scale(1./me->Integral());

  me->SetMaximum(1.);
  me->SetStats(false);

  return me;
}

void showTauIsolation(TFile* inputFile, const TString& dqmDirectory, const TString& meName,
		      const char* dqmSubDirectoryTauJet, const char* dqmSubDirectoryMuon, const char* dqmSubDirectoryElectron,
		      TCanvas* canvas, TPostScript* ps, const char* outputFileLabel, bool useLogScale)
{
  canvas->SetLogy(useLogScale);

  TLegend legend(0.74, 0.71, 0.89, 0.89, "", "brNDC"); 
  legend.SetBorderSize(0);
  legend.SetFillColor(0);

  TH1* meTauJet = getMonitorElement(inputFile, dqmDirectory, dqmSubDirectoryTauJet, meName);
  meTauJet->SetMarkerStyle(20);
  meTauJet->SetMarkerColor(kRed);
  meTauJet->SetLineColor(kRed);
  meTauJet->Draw("e1p");
  legend.AddEntry(meTauJet, "#tau-Jet", "p");

  if ( dqmSubDirectoryMuon != "" ) {
    TH1* meMuon = getMonitorElement(inputFile, dqmDirectory, dqmSubDirectoryMuon, meName);
    meMuon->SetMarkerStyle(21);
    meMuon->SetMarkerColor(kBlue);
    meMuon->SetLineColor(kBlue);
    meMuon->Draw("e1psame");
    legend.AddEntry(meMuon, "#mu", "p");
  }

  if ( dqmSubDirectoryElectron != "" ) {
    TH1* meElectron = getMonitorElement(inputFile, dqmDirectory, dqmSubDirectoryElectron, meName);
    meElectron->SetMarkerStyle(23);
    meElectron->SetMarkerColor(kGreen);
    meElectron->SetLineColor(kGreen);
    meElectron->Draw("e1psame");
    legend.AddEntry(meElectron, "e", "p");
  }

  legend.Draw();

  canvas->Update();
  TString outputFileName = TString("plot").Append(meTauJet->GetName()).Append("_").Append(outputFileLabel).Append(".png");
  canvas->Print(outputFileName.Data());
  //ps->NewPage();
}
		     

