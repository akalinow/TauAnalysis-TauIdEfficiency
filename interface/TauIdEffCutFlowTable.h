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
 * $Id: TauIdEffCutFlowTable.h,v 1.1 2011/07/10 15:47:27 veelken Exp $
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
  void bookCutFlowTable(TFileDirectory&);
  void fillCutFlowTable(double, const std::vector<bool>&, double);
  
  double getCutFlowNumber(int, int);

 protected:

  struct binType
  {
    binType(const std::string& name, const std::string& subdir, double min, double max, const std::vector<std::string>& selectionNames)
      : name_(name),
        subdir_(subdir),
        min_(min),
	max_(max),
	selectionNames_(selectionNames)
    {}
    ~binType() {}
    void bookCutFlowTable(TFileDirectory& dir)
    {
      TFileDirectory subdir = dir.mkdir(subdir_.data());

      std::string auxHistogramName = name_;
      size_t numSelections = selectionNames_.size();
      auxHistogram_ = dir.make<TH1D>(auxHistogramName.data(), auxHistogramName.data(), numSelections + 1, -0.5, numSelections + 0.5);

      TAxis* auxHistogramAxis = auxHistogram_->GetXaxis();
      for ( size_t iSelection = 0; iSelection < numSelections; ++iSelection ) {
	int auxHistogramBin = auxHistogram_->FindBin(iSelection);
	auxHistogramAxis->SetBinLabel(auxHistogramBin, selectionNames_[iSelection].data());
      }
    }

    std::string name_;
    std::string subdir_;
    double min_;
    double max_;
    std::vector<std::string> selectionNames_;

    TH1* auxHistogram_; // auxiliary histogram to store cut-flow numbers
  };

 private:

  /// specify process, region, tauIdDiscriminator and label
  /// used to uniquely identify bins/auxiliary histograms
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
