#ifndef TauAnalysis_TauIdEfficiency_TauIdEffCutFlowTable_h
#define TauAnalysis_TauIdEfficiency_TauIdEffCutFlowTable_h

/** \class TauIdEffCutFlowTable
 *
 * Keep track of number of muon + tau-jet pairs passing different preselection criteria 
 * applied on tau-jets considered in for tau id. efficiency measurement.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TauIdEffCutFlowManager.h,v 1.1 2011/07/05 08:32:02 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <TMath.h>
#include <TH1.h>
#include <TAxis.h>

#include <string>
#include <sstream>
#include <vector>

class TauIdEffCutFlowTable
{

 public:
  /// constructor
  TauIdEffCutFlowTable(edm::ParameterSet const&);

  /// destructor
  virtual ~TauIdEffCutFlowTable();

  /// book and fill auxiliary histograms
  /// used to keep track of numbers
  void bookCutFlowTable(TFileDirectory&, int, double, double);
  void bookCutFlowTable(TFileDirectory&, int, float*);
  void fillCutFlowTable(double, const std::vector<bool>&, double);
  
  double getCutFlowNumber(int, int);

 protected:

  struct binType
  {
    binType(TFileDirectory& dir, const std::string& name, double min, double max, const std::vector<std::string>& selectionNames)
      : min_(min),
	max_(max)
    {
      std::stringstream auxHistogramName;
      auxHistogramName << name;
      std::stringstream min_string;
      if ( min < 0. ) min_string << "Minus";      
      int min_int = TMath::Abs(TMath::Nint(min));
      min_string << min_int;
      std::stringstream max_string;
      if      ( max < 0. ) max_string << "Minus";      
      else if ( min < 0. ) max_string << "Plus";  
      int max_int = TMath::Abs(TMath::Nint(max));
      max_string << max_int;
      const double epsilon = 5.e-2;
      if ( TMath::Abs(min - min_int) > epsilon ||
	   TMath::Abs(max - max_int) > epsilon ) {
	min_string << TMath::Abs(TMath::Nint(10*(min - min_int)));
	max_string << TMath::Abs(TMath::Nint(10*(max - max_int)));
      }
      auxHistogramName << min_string.str() << "to" << max_string.str();
      size_t numSelections = selectionNames.size();
      //std::cout << "--> booking auxHistogramName = " << auxHistogramName.str() << std::endl;
      auxHistogram_ = dir.make<TH1D>(auxHistogramName.str().data(), 
				     auxHistogramName.str().data(), numSelections + 1, -0.5, numSelections + 0.5);
      TAxis* auxHistogramAxis = auxHistogram_->GetXaxis();
      for ( size_t iSelection = 0; iSelection < numSelections; ++iSelection ) {
	int auxHistogramBin = auxHistogram_->FindBin(iSelection);
	auxHistogramAxis->SetBinLabel(auxHistogramBin, selectionNames[iSelection].data());
      }
    }
    ~binType() {}

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
 
  std::vector<binType*> bins_;
};

#endif
