#include "TauAnalysis/TauIdEfficiency/interface/TauIdEffCutFlowManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

TauIdEffCutFlowManager::TauIdEffCutFlowManager(const edm::ParameterSet& cfg)
{
  binVariable_        = cfg.getParameter<std::string>("binVariable");

  process_            = cfg.getParameter<std::string>("process");
  region_             = cfg.getParameter<std::string>("region");
  tauIdDiscriminator_ = cfg.getParameter<std::string>("tauIdDiscriminator");
  label_              = cfg.getParameter<std::string>("label");
  //label_              = "presel";

  selectionNames_     = cfg.getParameter<vstring>("selectionNames");
  numSelections_      = selectionNames_.size();
  if ( !(numSelections_ >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowManager") 
      << "Collection of selection names must not be empty !!\n";

  binningName_ = std::string(process_).append("_").append(region_).append("_").append(binVariable_);
  binningName_.append("_").append(tauIdDiscriminator_).append("_").append(label_);
}

TauIdEffCutFlowManager::~TauIdEffCutFlowManager()
{
  for ( std::vector<bin*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    delete (*bin);
  }
}

void TauIdEffCutFlowManager::bookCutFlow(TFileDirectory& dir, int numBins, double min, double max)
{
  if ( !(numBins >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowManager") 
      << "Parameter numBins must be postive !!\n";
  
  double binWidth = (max - min)*numBins;
  
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double binMin = min + iBin*binWidth;
    double binMax = binMin + binWidth;

    bins_.push_back(new bin(dir, binningName_, binMin, binMax, selectionNames_));
  }
}

void TauIdEffCutFlowManager::bookCutFlow(TFileDirectory& dir, int numBins, float* binning)
{
  if ( !(numBins >= 1) ) 
    throw cms::Exception("TauIdEffCutFlowManager") 
      << "Parameter numBins must be postive !!\n";
  
  for ( int iBin = 0; iBin < numBins; ++iBin ) {
    double binMin = binning[iBin];
    double binMax = binning[iBin + 1];

    bins_.push_back(new bin(dir, binningName_, binMin, binMax, selectionNames_));
  }
}

void TauIdEffCutFlowManager::fillCutFlow(double x, const std::vector<bool>& selectionFlags, double weight)
{
  if ( selectionFlags.size() != numSelections_ )
    throw cms::Exception("TauIdEffCutFlowManager::fillCutFlow") 
      << "Parameter selectionFlags must be of length " << numSelections_ << " !!\n";
  
  for ( std::vector<bin*>::iterator bin = bins_.begin();
	bin != bins_.end(); ++bin ) {
    if ( x > (*bin)->min_ && x <= (*bin)->max_ ) {
      (*bin)->auxHistogram_->Fill(0., weight);
      bool isPassed = true;
      for ( size_t iSelection = 0; iSelection <= numSelections_; ++iSelection ) {
	if ( !selectionFlags[iSelection] ) isPassed = false;
	if ( isPassed ) (*bin)->auxHistogram_->Fill(iSelection + 1., weight);
	else break;
      }
    }
  }
}
  
 
