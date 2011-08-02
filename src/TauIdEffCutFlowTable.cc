#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffCutFlowTable.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

TauIdEffCutFlowTable::TauIdEffCutFlowTable(const edm::ParameterSet& cfg)
{
  //std::cout << "<TauIdEffCutFlowTable::TauIdEffCutFlowTable>:" << std::endl;

  process_            = cfg.getParameter<std::string>("process");
  region_             = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_ = cfg.getParameter<std::string>("tauIdDiscriminator");
  label_              = cfg.getParameter<std::string>("label");
  
  //std::cout << " process = " << process_ << std::endl;
  //std::cout << " region = " << region_ << std::endl;
  //std::cout << " tauIdDiscriminator = " << tauIdDiscriminator_ << std::endl;
  //std::cout << " label = " << label_ << std::endl;

  selectionNames_     = cfg.getParameter<vstring>("selectionNames");
  numSelections_      = selectionNames_.size();
  if ( !(numSelections_ >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowTable") 
      << "Collection of selection names must not be empty !!\n";
  
  binningName_ = std::string(process_).append("_").append(region_);
  binningName_.append("_").append(tauIdDiscriminator_).append("_").append(label_);

  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgBins = cfg.getParameter<vParameterSet>("binning");
  for ( vParameterSet::const_iterator cfgBin = cfgBins.begin();
	cfgBin != cfgBins.end(); ++cfgBin ) {
    std::string subdir = cfgBin->getParameter<std::string>("subdir");
    double min = cfgBin->getParameter<double>("min");
    double max = cfgBin->getParameter<double>("max");
    bins_.push_back(new binType(binningName_, subdir, min, max, selectionNames_));
  }

  if ( !(bins_.size() >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowTable") 
      << "Collection of bins must not be empty !!\n";
}

TauIdEffCutFlowTable::~TauIdEffCutFlowTable()
{
  for ( std::vector<binType*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    delete (*bin);
  }
}

void TauIdEffCutFlowTable::bookCutFlowTable(TFileDirectory& dir)
{
  for ( std::vector<binType*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    if ( (*bin)->subdir_ != "" ) {
      TFileDirectory subdir = dir.mkdir((*bin)->subdir_);
      (*bin)->bookCutFlowTable(subdir);
    } else {
      (*bin)->bookCutFlowTable(dir);
    }
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
