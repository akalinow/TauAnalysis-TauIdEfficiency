#ifndef TauAnalysis_TauIdEfficiency_TauIdEffHistManager_h
#define TauAnalysis_TauIdEfficiency_TauIdEffHistManager_h

/** \class TauIdEffHistManager
 *
 * Fill histograms for tau id. efficiency measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TauIdEffHistManager.h,v 1.1 2011/07/01 10:41:48 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

class TauIdEffHistManager
{

 public:
  /// constructor
  TauIdEffHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffHistManager();

  /// book and fill histograms
  void bookHistograms(fwlite::TFileService&);
  void fillHistograms(const PATMuTauPair&, double);

 protected:

  TH1* book1D(fwlite::TFileService& fs, const std::string&, const std::string&, int, double, double);
  TH1* book1D(fwlite::TFileService& fs, const std::string&, const std::string&, int, float*);

  std::string getHistogramName(const std::string&);

 private:

  /// specify process, region, tauIdDiscriminator and label
  /// used to uniquely identify histograms
  std::string process_; 
  std::string region_;
  std::string tauIdDiscriminator_; 
  std::string label_;
 
  TH1* histogramMuonPt_;
  TH1* histogramMuonEta_;
  TH1* histogramMuonPhi_;

  TH1* histogramTauPt_;
  TH1* histogramTauEta_;
  TH1* histogramTauPhi_;
  TH1* histogramTauNumTracks_;

  TH1* histogramVisMass_;
  TH1* histogramMt_;
};

#endif
