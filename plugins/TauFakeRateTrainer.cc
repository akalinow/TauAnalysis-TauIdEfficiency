
#include <FWCore/Utilities/interface/Exception.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/MVAComputer/interface/MVAModuleHelper.h"
#include "RecoTauTag/TauTagTools/interface/TauMVADBConfiguration.h"

#include "TFile.h"
#include "TPolyMarker3D.h"

namespace {

const PhysicsTools::AtomicId ptId("JetPt");
const PhysicsTools::AtomicId etaId("JetEta");
const PhysicsTools::AtomicId widthId("JetWidth");

// The Phys tools version of this doesn't compile due to const issues.
struct Extractor {
  const TPolyMarker3D* data_;
  size_t index_;
  double pt_;
  double eta_;
  double width_;
  void setIndex(size_t index) {
    index_ = index;
    data_->GetPoint(index_, pt_, eta_, width_);
  }
  size_t size() {
    return data_->Size();
  }
  double compute(const PhysicsTools::AtomicId &name) const {
    if (name == ptId)
      return pt_;
    if (name == etaId)
      return eta_;
    if (name == widthId)
      return width_;
    return -1000;
  }
};

struct ExtractorFiller {
	ExtractorFiller(const PhysicsTools::AtomicId &name) {}
	double operator()(const Extractor &object,
	                  const PhysicsTools::AtomicId &name) const
	{ return object.compute(name); }
};


}

class TauFakeRateTrainer : public edm::EDAnalyzer {
  public:
    explicit TauFakeRateTrainer(const edm::ParameterSet&);
    ~TauFakeRateTrainer() {}
    void analyze(const edm::Event& evt, const edm::EventSetup& es);
  private:
    TFile* numeratorFile;
    TFile* denominatorFile;
    const TPolyMarker3D* numerator_tuple;
    const TPolyMarker3D* denominator_tuple;
    PhysicsTools::MVAModuleHelper<TauTagMVAComputerRcd, Extractor,
        ExtractorFiller> helper_;
};


TauFakeRateTrainer::TauFakeRateTrainer(const edm::ParameterSet &pset)
  :helper_("train") {
  // Load the stuff.
  numeratorFile = TFile::Open(
      pset.getParameter<std::string>("passing").c_str(), "READ");
  if (!numeratorFile) {
    throw cms::Exception("NumeratorFile") << "Can't open passing file!";
  }
  numerator_tuple = dynamic_cast<const TPolyMarker3D*>(
      numeratorFile->Get("passing"));
  if (!numerator_tuple) {
    throw cms::Exception("Missing ntuple") << "Can't Get() passing ntuple";
  }
  denominatorFile = TFile::Open(
      pset.getParameter<std::string>("failing").c_str(), "READ");
  if (!denominatorFile) {
    throw cms::Exception("denominatorFile") << "Can't open failing file!";
  }
  denominator_tuple = dynamic_cast<const TPolyMarker3D*>(
      denominatorFile->Get("failing"));
  if (!denominator_tuple) {
    throw cms::Exception("Missing ntuple") << "Can't Get() failing ntuple";
  }
}

void TauFakeRateTrainer::analyze(const edm::Event& evt,
                                 const edm::EventSetup& es) {
  helper_.setEventSetup(es, "trainer");
  // Build our extractor objects
  Extractor numerator;
  numerator.data_ = numerator_tuple;
  numerator.setIndex(0);

  Extractor denominator;
  denominator.data_ = denominator_tuple;
  denominator.setIndex(0);
  // Pass the training events to the MVA
  for(size_t i = 0; i < numerator.size(); ++i) {
    numerator.setIndex(i);
    helper_.train(numerator, true, 1.0);
  }
  for(size_t i = 0; i < denominator.size(); ++i) {
    denominator.setIndex(i);
    helper_.train(denominator, false, 1.0);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauFakeRateTrainer);
