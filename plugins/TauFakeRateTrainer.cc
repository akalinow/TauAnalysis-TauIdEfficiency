
#include <FWCore/Utilities/interface/Exception.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/MVAComputer/interface/MVAModuleHelper.h"
#include "RecoTauTag/TauTagTools/interface/TauMVADBConfiguration.h"

#include "TFile.h"
#include "TArrayD.h"
#include "TPolyMarker3D.h"

namespace {

// This module doesn't know the difference between tau and jet pt.
const PhysicsTools::AtomicId widthId("JetWidth");
const PhysicsTools::AtomicId tauPtId("Pt");
const PhysicsTools::AtomicId tauEtaId("AbsEta");
const PhysicsTools::AtomicId jetPtId("JetPt");
const PhysicsTools::AtomicId jetEtaId("AbsJetEta");

// Check if a value is Nan and print a warning if so
double checkNan(const PhysicsTools::AtomicId &name, double value) {
  if (std::isnan(value)) {
    std::cerr << "Found a nan in variable: " << name << ", returning zero!"
      << std::endl;
    value = 0.;
  }
  return value;
}

// The Phys tools version of this doesn't compile due to const issues.
struct Extractor {
  const TPolyMarker3D* data_;
  size_t index_;
  double pt_;
  double eta_;
  double width_;
  double weight_;
  void setIndex(size_t index) {
    index_ = index;
    data_->GetPoint(index_, pt_, eta_, width_);
  }
  size_t size() const {
    return data_->Size();
  }
  double compute(const PhysicsTools::AtomicId &name) const {
    if (name == tauPtId || name == jetPtId)
      return checkNan(tauPtId, pt_);
    if (name == tauEtaId || name == jetEtaId)
      return checkNan(tauEtaId, eta_);
    if (name == widthId)
      return checkNan(widthId, width_);
    return -1000;
  }
};

struct ExtractorFiller {
	ExtractorFiller(const PhysicsTools::AtomicId &name) {}
	double operator()(const Extractor &object,
	                  const PhysicsTools::AtomicId &name) const
	{ return object.compute(name); }
};

// Print list of indices to the screen
void printIndices(const TPolyMarker3D* ntuple, const std::string& prefix) {
  double fId, fEvt, fRun;
  std::string sep = ":";
  for (int i = 0; i < ntuple->Size(); ++i) {
    ntuple->GetPoint(i, fId, fEvt, fRun);
    std::cout << prefix
      << sep
      << lrint(fId)
      << sep
      << lrint(fEvt)
      << sep
      << lrint(fRun)
      << std::endl;
  }
}

class TauFakeRateTrainer : public edm::EDAnalyzer {
  public:
    explicit TauFakeRateTrainer(const edm::ParameterSet&);
    ~TauFakeRateTrainer() {}
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    TFile* numeratorFile;
    TFile* denominatorFile;
    std::vector<Extractor> numerators_;
    std::vector<Extractor> denominators_;
    PhysicsTools::MVAModuleHelper<TauTagMVAComputerRcd, Extractor,
        ExtractorFiller> helper_;
};


TauFakeRateTrainer::TauFakeRateTrainer(const edm::ParameterSet &pset)
  :helper_("train") {
  // Load the stuff.
  std::cout << "Opening numerator file" << std::endl;
  numeratorFile = TFile::Open(
      pset.getParameter<std::string>("passing").c_str(), "READ");
  if (!numeratorFile) {
    throw cms::Exception("NumeratorFile") << "Can't open passing file!";
  }
  // Get the weights - they should be identical
  std::cout << "Opening numerator weights object" << std::endl;
  const TArrayD* num_weights;
  numeratorFile->GetObject("weights", num_weights);
  if (!num_weights) {
    throw cms::Exception("Missing ntuple") << "Can't Get() passing ntuples weights";
  }

  std::cout << "Opening denominator file" << std::endl;
  denominatorFile = TFile::Open(
      pset.getParameter<std::string>("failing").c_str(), "READ");
  if (!denominatorFile) {
    throw cms::Exception("denominatorFile") << "Can't open failing file!";
  }
  std::cout << "Opening denominator weights object" << std::endl;
  const TArrayD* den_weights;
  denominatorFile->GetObject("weights", den_weights);
  if (!den_weights) {
    throw cms::Exception("Missing ntuple") << "Can't Get() failing ntuples weights";
  }
  std::cout << "Checking weights" << std::endl;
  // Make sure they are the same
  bool weights_ok = (num_weights->GetSize() == den_weights->GetSize());
  for (int i = 0; i < num_weights->GetSize(); i++) {
    if (num_weights->At(i) != den_weights->At(i))
      weights_ok = false;
  }
  std::cout << "Weights are okay!" << std::endl;
  if (!weights_ok)
    throw cms::Exception("Missing ntuple") << "Weights in numerator and denominator files don't match!";

  // Build each extractor source
  for (int i = 0; i < num_weights->GetSize(); i++) {
    std::cout << "Building sample " << i << std::endl;
    std::stringstream passing_name;
    passing_name << "passing_";
    passing_name << i;
    std::stringstream failing_name;
    failing_name << "failing_";
    failing_name << i;
    // The names of the poly markers holding the indices of the jet events
    std::stringstream passing_index_name;
    passing_index_name << "passing_indices_";
    passing_index_name << i;
    std::stringstream failing_index_name;
    failing_index_name << "failing_indices_";
    failing_index_name << i;

    // Build extractors
    Extractor numerator;
    numerator.weight_ = num_weights->At(i);
    numerator.data_ = dynamic_cast<const TPolyMarker3D*>(
        numeratorFile->Get(passing_name.str().c_str()));
    if (!numerator.data_) {
      throw cms::Exception("Missing ntuple") << "Can't Get() passing ntuple";
    }
    numerator.setIndex(0);

    const TPolyMarker3D* passingIndices = dynamic_cast<const TPolyMarker3D*>(
        numeratorFile->Get(passing_index_name.str().c_str()));
    if (passingIndices) {
      std::cout << "Found passing indices:" << std::endl;
      printIndices(passingIndices, "index_pass");
    }

    Extractor denominator;
    denominator.weight_ = den_weights->At(i);
    denominator.data_ = dynamic_cast<const TPolyMarker3D*>(
        denominatorFile->Get(failing_name.str().c_str()));
    if (!denominator.data_) {
      throw cms::Exception("Missing ntuple") << "Can't Get() failing ntuple";
    }

    const TPolyMarker3D* failingIndices = dynamic_cast<const TPolyMarker3D*>(
        denominatorFile->Get(failing_index_name.str().c_str()));
    if (failingIndices) {
      std::cout << "Found failing indices:" << std::endl;
      printIndices(failingIndices, "index_fail");
    }

    denominator.setIndex(0);

    numerators_.push_back(numerator);
    denominators_.push_back(denominator);
  }
}

void TauFakeRateTrainer::analyze(const edm::Event& evt,
                                 const edm::EventSetup& es) {
  helper_.setEventSetup(es, "trainer");
  // Build our extractor objects

  // Pass the training events to the MVA
  for (size_t source = 0; source < numerators_.size(); ++source) {
    Extractor& numerator = numerators_[source];
    for(size_t i = 0; i < numerator.size(); ++i) {
      numerator.setIndex(i);
      helper_.train(numerator, true, numerator.weight_);
    }
  }

  for (size_t source = 0; source < denominators_.size(); ++source) {
    Extractor& denominator = denominators_[source];
    for(size_t i = 0; i < denominator.size(); ++i) {
      denominator.setIndex(i);
      helper_.train(denominator, false, denominator.weight_);
    }
  }
}

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauFakeRateTrainer);
