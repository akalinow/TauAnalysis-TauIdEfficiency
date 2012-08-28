#ifndef TauAnalysis_TauIdEfficiency_TauIdEffHistManager_h
#define TauAnalysis_TauIdEfficiency_TauIdEffHistManager_h

/** \class TauIdEffHistManager
 *
 * Fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.13 $
 *
 * $Id: TauIdEffHistManager.h,v 1.13 2012/06/18 19:15:35 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include <TH1.h>

class TauIdEffHistManager
{

 public:
  /// constructor
  TauIdEffHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffHistManager();

  /// book and fill histograms
  void bookHistograms(TFileDirectory&);
  void fillHistograms(const PATMuTauPair&, const pat::MET&, size_t, size_t, size_t, const std::map<std::string, bool>&, double);
  
  /// scale all bin-contents/bin-errors by factor given as function argument
  /// (to account for events lost, due to aborted skimming/crab or PAT-tuple production/lxbatch jobs)
  void scaleHistograms(double);

 protected:

  TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, double, double);
  TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, float*);

  std::string getHistogramName(const std::string&);

 private:

  /// specify process, region, tauIdDiscriminator and label
  /// used to uniquely identify histograms
  std::string process_; 
  std::string region_;
  std::string tauIdDiscriminator_; 
  std::string label_;

  /// specify name of SVfit mass hypothesis to be used
  /// as template variable
  std::string svFitMassHypothesis_;
  
  /// CV: only book histograms needed to run tau id. efficiency fit;
  ///     drop all control plots
  bool fillControlPlots_;

  typedef std::vector<std::string> vstring;
  vstring triggerBits_;
  
  TH1* histogramMuonPt_;
  TH1* histogramMuonEta_;
  TH1* histogramMuonPhi_;

  TH1* histogramTauPt_;
  TH1* histogramTauEta_;
  TH1* histogramTauPhi_;
  TH1* histogramTauNumTracks_;
  TH1* histogramTauNumSelTracks_;

  TH1* histogramVisMass_;
  TH1* histogramSVfitMass_;
  TH1* histogramMt_;
  TH1* histogramPzetaDiff_;
  TH1* histogramDPhi_;

  TH1* histogramNumJets_;
  TH1* histogramNumJetsBtagged_;

  TH1* histogramPFMEt_;
  TH1* histogramPFSumEt_;
  TH1* histogramCaloMEt_;
  TH1* histogramCaloSumEt_;

  TH1* histogramNumVertices_;

  TH1* histogramLogEvtWeight_;

  std::map<std::string, TH1*> histogramNumCaloMEt_; // key = L1 bit; numerator for trigger efficiency control plots
  TH1* histogramDenomCaloMEt_;                      // denominator for trigger efficiency control plots

  TH1* histogramEventCounter_;
  
  std::vector<TH1*> histograms_;
};

#endif
