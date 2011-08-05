#ifndef TauAnalysis_TauIdEfficiency_TauFakeRateHistManager_h
#define TauAnalysis_TauIdEfficiency_TauFakeRateHistManager_h

/** \class TauFakeRateHistManager
 *
 * Fill histograms for jet --> tau fake-rate measurement
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TauFakeRateHistManager.h,v 1.1 2011/07/18 16:40:45 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include <TH1.h>

class TauFakeRateHistManager
{

 public:
  /// constructor
  TauFakeRateHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~TauFakeRateHistManager();

  /// book and fill histograms
  void bookHistograms(TFileDirectory&);
  void fillHistograms(const pat::Tau&, size_t, double, double);
  
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
 
  TH1* histogramJetPt_;
  TH1* histogramJetEta_;
  TH1* histogramJetPhi_;

  TH1* histogramTauPt_;
  TH1* histogramTauEta_;
  TH1* histogramTauPhi_;

  TH1* histogramSumEt_;
  TH1* histogramNumVertices_;
  
  std::vector<TH1*> histograms_;
};

#endif
