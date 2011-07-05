#ifndef TauAnalysis_TauIdEfficiency_TauIdEffCutFlowManager_h
#define TauAnalysis_TauIdEfficiency_TauIdEffCutFlowManager_h

/** \class TauIdEffCutFlowManager
 *
 * Keep track of number of muon + tau-jet pairs passing different preselection criteria 
 * applied on tau-jets considered in for tau id. efficiency measurement.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: TauIdEffBinner.h,v 1.4 2011/07/04 09:51:23 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <TMath.h>

#include <string>
#include <sstream>
#include <vector>

class TauIdEffCutFlowManager
{

 public:
  /// constructor
  TauIdEffCutFlowManager(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffCutFlowManager();

  /// book and fill auxiliary histograms
  /// used to keep track of numbers
  void bookCutFlow(TFileDirectory&, int, double, double);
  void bookCutFlow(TFileDirectory&, int, float*);
  void fillCutFlow(double, const std::vector<bool>&, double);
  
 protected:

  struct bin
  {
    bin(TFileDirectory& dir, const std::string& name, double min, double max, const std::vector<std::string>& selectionNames)
      : min_(min),
	max_(max)
    {
      std::stringstream auxHistogramName;
      auxHistogramName << name;
      int min_int = TMath::Nint(min);
      int max_int = TMath::Nint(max);
      const double epsilon = 1.e-1;
      if ( TMath::Abs(min - min_int) > epsilon ||
	   TMath::Abs(max - max_int) > epsilon ) {
	min_int *= 10;
	max_int *= 10;
      }
      auxHistogramName << min_int << "to" << max_int;
      size_t numSelections = selectionNames.size();
      auxHistogram_ = dir.make<TH1D>(auxHistogramName.str().data(), 
				     auxHistogramName.str().data(), numSelections + 1, -0.5, numSelections + 0.5);
    }
    ~bin() {}

    double min_;
    double max_;

    TH1* auxHistogram_; // auxiliary histogram to store cut-flow numbers
  };

 private:

  /// specify binVariable, process, region, tauIdDiscriminator and label
  /// used to uniquely identify bins/auxiliary histograms
  std::string binVariable_;

  std::string process_; 
  std::string region_;
  std::string tauIdDiscriminator_; 
  std::string label_;

  typedef std::vector<std::string> vstring;
  vstring selectionNames_;
  size_t numSelections_;

  std::string binningName_;
 
  std::vector<bin*> bins_;
};

#endif
