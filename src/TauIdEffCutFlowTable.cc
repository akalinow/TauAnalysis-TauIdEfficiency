#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffCutFlowTable.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

TauIdEffCutFlowTable::TauIdEffCutFlowTable(const edm::ParameterSet& cfg)
{
  binVariable_        = cfg.getParameter<std::string>("binVariable");

  process_            = cfg.getParameter<std::string>("process");
  region_             = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_ = cfg.getParameter<std::string>("tauIdDiscriminator");
  label_              = cfg.getParameter<std::string>("label");

  selectionNames_     = cfg.getParameter<vstring>("selectionNames");
  numSelections_      = selectionNames_.size();
  if ( !(numSelections_ >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowTable") 
      << "Collection of selection names must not be empty !!\n";

  binningName_ = std::string(process_).append("_").append(region_);
  if ( binVariable_ != "" ) binningName_.append("_").append(binVariable_);
  binningName_.append("_").append(tauIdDiscriminator_).append("_").append(label_);
}

TauIdEffCutFlowTable::~TauIdEffCutFlowTable()
{
  for ( std::vector<binType*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    delete (*bin);
  }
}

void TauIdEffCutFlowTable::bookCutFlowTable(TFileDirectory& dir, int numBins, double min, double max)
{
  if ( !(numBins >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowTable") 
      << "Parameter numBins must be postive !!\n";
  
  double binWidth = (max - min)*numBins;
  
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double binMin = min + iBin*binWidth;
    double binMax = binMin + binWidth;

    bins_.push_back(new binType(dir, binningName_, binMin, binMax, selectionNames_));
  }
}

void TauIdEffCutFlowTable::bookCutFlowTable(TFileDirectory& dir, int numBins, float* binning)
{
  if ( !(numBins >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowTable") 
      << "Parameter numBins must be postive !!\n";
  
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double binMin = binning[iBin];
    double binMax = binning[iBin + 1];

    bins_.push_back(new binType(dir, binningName_, binMin, binMax, selectionNames_));
  }
}

void TauIdEffCutFlowTable::fillCutFlowTable(double x, const std::vector<bool>& selectionFlags, double weight)
{
  if ( selectionFlags.size() != numSelections_ )
    throw cms::Exception("TauIdEffCutFlowTable::fillCutFlowTable") 
      << "Parameter selectionFlags must be of length " << numSelections_ << " !!\n";
  
  for ( std::vector<binType*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    if ( x > (*bin)->min_ && x <= (*bin)->max_ ) {
      bool isPassed = true;
      for ( size_t iSelection = 0; iSelection <= numSelections_; ++iSelection ) {
	if ( !selectionFlags[iSelection] ) isPassed = false;
	if ( isPassed ) (*bin)->auxHistogram_->Fill(iSelection, weight);
	else break;
      }
    }
  }
}

double TauIdEffCutFlowTable::getCutFlowNumber(int binIdx, int row)
{
  if ( !(binIdx >= 0 && binIdx < (int)bins_.size()) )
    throw cms::Exception("TauIdEffCutFlowTable::getCutFlowNumber") 
      << "Parameter binIdx outside range 0.." << (bins_.size() - 1) << " !!\n";

  if ( !(row >= 0 && row <= (int)numSelections_) )
    throw cms::Exception("TauIdEffCutFlowTable::getCutFlowNumber") 
      << "Parameter row outside range 0.." << numSelections_ << " !!\n";

  binType* bin = bins_[binIdx];
  int auxHistogramBin = bin->auxHistogram_->FindBin(row);
  return bin->auxHistogram_->GetBinContent(auxHistogramBin);
}
