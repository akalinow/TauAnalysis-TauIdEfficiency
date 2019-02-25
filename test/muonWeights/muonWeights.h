#ifndef muonWeights_H
#define muonWeights_H

#include <string>
#include "TH1F.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"

class ChannelSpecifics {

 public:
  
  ChannelSpecifics();
  
  ~ChannelSpecifics();

  void initializeLeptonCorrections();
  void initializePileUpCorrections();

  float getLeptonCorrection(float eta, float pt, float iso);
  float getPUWeight(float nPU);

  void addDataWeightBranch(std::string path);
  void addMCWeightBranch(std::string path);

 private:

  TH1F *hPUWeight;
  
  RooAbsReal *muon_id_iso_scalefactor, *muon_trg_scalefactor, *muon_trk_scalefactor;
  RooRealVar *m_eta, *m_pt, *m_iso;

};

#endif
