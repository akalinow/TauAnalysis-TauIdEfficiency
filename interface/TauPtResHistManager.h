#ifndef TauAnalysis_TauIdEfficiency_TauPtResHistManager_h
#define TauAnalysis_TauIdEfficiency_TauPtResHistManager_h

/** \class TauPtResHistManager
 *
 * Fill histograms for tau energy reconstruction study
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.5 $
 *
 * $Id: TauPtResHistManager.h,v 1.5 2011/07/21 16:37:13 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"

#include <TH1.h>
#include <TH2.h>

class TauPtResHistManager
{

 public:
  /// constructor
  TauPtResHistManager(edm::ParameterSet const&);

  /// destructor
  virtual ~TauPtResHistManager();

  /// book and fill histograms
  void bookHistograms(TFileDirectory&);
  void fillHistograms(const pat::Tau&, const reco::GenParticleCollection&, const reco::Vertex&, 
		      const reco::tau::RecoTauQualityCuts&, double);
  
 protected:

  static TH1* book1D(TFileDirectory&, const std::string&, const std::string&, int, double, double);
  static TH2* book2D(TFileDirectory&, const std::string&, const std::string&, int, double, double, int, double, double);

  static std::string getHistogramName(const std::string&);

 private:

  TH1* histogramTauPtRes_;
  TH1* histogramTauPtResGenOneProng0Pi0_;
  TH1* histogramTauPtResGenOneProng1Pi0_;
  TH1* histogramTauPtResGenOneProng2Pi0_;
  TH1* histogramTauPtResGenThreeProng0Pi0_;
  TH1* histogramTauPtResGenThreeProng1Pi0_;

  TH2* histogramEtaPhiForTauPtResLt05GenOneProng0Pi0_;
  TH2* histogramEtaPhiForTauPtResLt05GenOneProng1Pi0_;
  TH2* histogramEtaPhiForTauPtResLt05GenOneProng2Pi0_;
  TH2* histogramEtaPhiForTauPtResLt05GenThreeProng0Pi0_;
  TH2* histogramEtaPhiForTauPtResLt05GenThreeProng1Pi0_;

  struct tauPtResManCorrHistograms
  {
    tauPtResManCorrHistograms(int);
    ~tauPtResManCorrHistograms() {}

    void bookHistograms(TFileDirectory&);
    void fillHistograms(const pat::Tau&, const std::string&, double, const reco::Vertex&, 
			const reco::tau::RecoTauQualityCuts&, double);

    int level_;

    TH1* histogramTauPtRes_;
    TH1* histogramTauPtResGenOneProng0Pi0_;
    TH1* histogramTauPtResGenOneProng1Pi0_;
    TH1* histogramTauPtResGenOneProng2Pi0_;
    TH1* histogramTauPtResGenThreeProng0Pi0_;
    TH1* histogramTauPtResGenThreeProng1Pi0_;
  };

  tauPtResManCorrHistograms histogramsTauPtResManCorrLev1_;
  tauPtResManCorrHistograms histogramsTauPtResManCorrLev2_;
  tauPtResManCorrHistograms histogramsTauPtResManCorrLev12_;
  tauPtResManCorrHistograms histogramsTauPtResManCorrLev123_;
  tauPtResManCorrHistograms histogramsTauPtResManCorrLev14_;
  
  TH1* histogramJetPtRes_;
  TH1* histogramJetPtResGenOneProng0Pi0_;
  TH1* histogramJetPtResGenOneProng1Pi0_;
  TH1* histogramJetPtResGenOneProng2Pi0_;
  TH1* histogramJetPtResGenThreeProng0Pi0_;
  TH1* histogramJetPtResGenThreeProng1Pi0_;

  TH2* histogramRecVsGenTauDecayMode_;
};

#endif
