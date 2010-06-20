ttree_draw()
{
  TFile* file = TFile::Open("tauIdEffEDNtuple_qcdDiJet.root");
  
  TTree* tree = (TTree*)file->Get("Events");
  //tree->Print();

  //TString algorithm = "patPFTausDijetTagAndProbeShrinkingCone";
  TString algorithm = "patPFTausDijetTagAndProbeHPS";

  TString selection;
  selection.Append("abs(tauIdEffNtuple#").Append(algorithm).Append("#jetEta) < 2.5");
  selection.Append(" && tauIdEffNtuple#").Append(algorithm).Append("#jetPt > 10.");
  selection.Append(" && tauIdEffNtuple#").Append(algorithm).Append("#probe > 0.5");

  TString expression;
  expression.Append("tauIdEffNtuple#").Append(algorithm).Append("#jetPt");

  tree->Draw(expression.Data(), selection.Data());

  delete file;
}
