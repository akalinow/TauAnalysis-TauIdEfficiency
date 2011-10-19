#ifndef TauAnalysis_TauIdEfficiency_MuonIsolationHistManager_h
#define TauAnalysis_TauIdEfficiency_MuonIsolationHistManager_h

/** \class MuonIsolationHistManager
 *
 * Fill histograms for measurement of probability of muons selected in QCD background events
 * to pass tight isolation, given that muon satisfies the loose isolation criteria applied on trigger level
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.7 $
 *
 * $Id: MuonIsolationHistManager.h,v 1.7 2011/08/10 16:23:07 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <TH1.h>
#include <TH2.h>

class MuonIsolationHistManager
{

 public:
  /// constructor
  MuonIsolationHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~MuonIsolationHistManager();

  /// book and fill histograms
  void bookHistograms(TFileDirectory&);
  void fillHistograms(const PATMuTauPair&, size_t, double, double);
  
 protected:

  TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, double, double);
  TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, float*);
  TH2* book2D(TFileDirectory&, const std::string&, const std::string&, int, float*, int, float*);

  std::string getHistogramName(const std::string&);

 private:

  TH1* histogramMuonPt_;
  TH1* histogramMuonEta_;
  TH2* histogramMuonPtVsEta_;
  TH1* histogramMuonPhi_;
  TH1* histogramMuonCharge_;
  TH1* histogramMuonAbsIso_;
  TH1* histogramMuonRelIso_;

  TH1* histogramTauJetPt_;
  TH1* histogramTauJetEta_;
  TH1* histogramTauJetPhi_;

  TH1* histogramVisMass_;
  TH1* histogramMt_;
  TH1* histogramPzetaDiff_;

  TH1* histogramMEt_;
  TH1* histogramSumEt_;
  TH1* histogramNumVertices_;
  
  std::vector<TH1*> histograms_;
};

#endif
