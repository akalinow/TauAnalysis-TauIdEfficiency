#include "TauAnalysis/TauIdEfficiency/interface/tauPtResAuxFunctions.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

std::string getGenTauDecayMode(const pat::Tau& patTau, const reco::GenParticleCollection& genParticles)
{
  const reco::GenParticle* genTau = findGenParticle(patTau.p4(), genParticles);
  if      ( genTau          ) return getGenTauDecayMode(genTau);
  else if ( patTau.genJet() ) return JetMCTagUtils::genTauDecayMode(*patTau.genJet());
  else return "";
}

std::string getPFCandidateType(const reco::PFCandidate& pfCandidate)
{
  reco::PFCandidate::ParticleType pfCandidateType = pfCandidate.particleId();
  if      ( pfCandidateType == reco::PFCandidate::X         ) return "undefined";
  else if ( pfCandidateType == reco::PFCandidate::h         ) return "PFChargedHadron";
  else if ( pfCandidateType == reco::PFCandidate::e         ) return "PFElectron";
  else if ( pfCandidateType == reco::PFCandidate::mu        ) return "PFMuon";
  else if ( pfCandidateType == reco::PFCandidate::gamma     ) return "PFGamma";
  else if ( pfCandidateType == reco::PFCandidate::h0        ) return "PFNeutralHadron";
  else if ( pfCandidateType == reco::PFCandidate::h_HF      ) return "HF_had";
  else if ( pfCandidateType == reco::PFCandidate::egamma_HF ) return "HF_em";
  else assert(0);
}

const reco::TrackBaseRef getTrack(const reco::PFCandidate& cand) 
{
  if      ( cand.trackRef().isNonnull()    ) return reco::TrackBaseRef(cand.trackRef());
  else if ( cand.gsfTrackRef().isNonnull() ) return reco::TrackBaseRef(cand.gsfTrackRef());
  else return reco::TrackBaseRef();
}

//
//-------------------------------------------------------------------------------
//

bool isTauSignalPFCandidate(const pat::Tau& patTau, const reco::PFCandidatePtr& pfJetConstituent)
{
  bool retVal = false;
  
  const reco::PFCandidateRefVector& signalPFCandidates = patTau.signalPFCands();
  for ( reco::PFCandidateRefVector::const_iterator signalPFCandidate = signalPFCandidates.begin();
	signalPFCandidate != signalPFCandidates.end(); ++signalPFCandidate ) {
    if ( pfJetConstituent.key() == signalPFCandidate->key() ) retVal = true;
  }
  
  return retVal;
}

double square(double x)
{
  return x*x;
}

double getTauPtManCorr(const pat::Tau& patTau, const reco::Vertex& vertex, 
		       const reco::tau::RecoTauQualityCuts& qualityCuts, unsigned corrLevel)
{
  double retVal = 0.;

  bool needsManCorr = false;
  std::vector<reco::PFCandidatePtr> pfJetConstituents = patTau.pfJetRef()->getPFConstituents();
  for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
    const reco::TrackBaseRef track = getTrack(**pfJetConstituent);
    if ( track.isNonnull() ) {
      double dIP = TMath::Abs(track->dxy(vertex.position()));
      double dZ = TMath::Abs(track->dz(vertex.position()));
      if ( corrLevel & 1 && track->pt() > 2.0 && TMath::Abs(patTau.eta() - (*pfJetConstituent)->eta()) < 0.10 &&
	   !isTauSignalPFCandidate(patTau, *pfJetConstituent) && (dZ < 0.2 || dIP > 0.10) ) needsManCorr = true;

      double trackPt = track->pt();
      double trackPtErr = track->ptError();
      double caloEn = (*pfJetConstituent)->ecalEnergy() + (*pfJetConstituent)->hcalEnergy();
      const double caloEnRes = 1.00;
      double caloEnErr = caloEnRes*TMath::Sqrt(caloEn);
      double caloEt = caloEn*TMath::Sin((*pfJetConstituent)->theta());
      double caloEtErr = caloEnErr*TMath::Sin((*pfJetConstituent)->theta());
      if ( corrLevel & 2 && trackPt > 2.0 && TMath::Abs(patTau.eta() - (*pfJetConstituent)->eta()) < 0.10 &&
	   trackPt > (caloEt + 2.0*TMath::Sqrt(trackPtErr*trackPtErr + caloEtErr*caloEtErr)) ) needsManCorr = true;
    }
  }

  double unaccountedPFCandPtSum = 0.;
  for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
    if ( TMath::Abs(patTau.eta() - (*pfJetConstituent)->eta()) < 0.10 &&
	 !isTauSignalPFCandidate(patTau, *pfJetConstituent) ) unaccountedPFCandPtSum += (*pfJetConstituent)->pt();
  }
  if ( corrLevel & 4 && unaccountedPFCandPtSum > (0.50*patTau.pt()) ) needsManCorr = true;

  if ( needsManCorr ) {
    for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	  pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
      if ( TMath::Abs(patTau.eta() - (*pfJetConstituent)->eta()) < 0.10 &&
	   TMath::Abs(patTau.phi() - (*pfJetConstituent)->phi()) < 0.50 ) {
	if ( (*pfJetConstituent)->particleId() == reco::PFCandidate::h0 ) {
	  retVal += (*pfJetConstituent)->pt();
	} else {
	  if ( !isTauSignalPFCandidate(patTau, *pfJetConstituent) ) {
	    double caloEn = (*pfJetConstituent)->ecalEnergy() + (*pfJetConstituent)->hcalEnergy();
	    double caloEt = caloEn*TMath::Sin((*pfJetConstituent)->theta());
	    retVal += caloEt;
	  }
	}
      }
    }
  }

  double leadTrackMom = 0.;
  double leadTrackMomErr = 0.;
  double jetCaloEn = 0.;

  for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
    const reco::TrackBaseRef track = getTrack(**pfJetConstituent);
    if ( track.isNonnull() ) {
      double trackPt = track->pt();
      double trackPtErr = track->ptError();
      if ( qualityCuts.filter(**pfJetConstituent) && 
	   trackPtErr < (0.20*trackPt) && track->normalizedChi2() < 5.0 && track->hitPattern().numberOfValidPixelHits() >= 1 &&
	   (trackPt - 3.*trackPtErr) > (*pfJetConstituent)->pt() && trackPt < (3.*patTau.pfJetRef()->pt()) ) {
	if ( track->p() > leadTrackMom ) {
	  leadTrackMom = track->p();
	  leadTrackMomErr = leadTrackMom*(trackPtErr/trackPt);
	}
      }
    }

    double caloEn = (*pfJetConstituent)->ecalEnergy() + (*pfJetConstituent)->hcalEnergy();
    jetCaloEn += caloEn;
  }

  if ( corrLevel & 8 && leadTrackMom > patTau.p() ) {
    const double chargedPionMass = 0.13957; // GeV
    double leadTrackEn = TMath::Sqrt(square(leadTrackMom) + square(chargedPionMass));
    double jetCaloEnErr = 1.00*TMath::Sqrt(TMath::Max(jetCaloEn, leadTrackEn));
    double combEn = ((1./square(jetCaloEnErr))*jetCaloEn + (1./square(leadTrackMomErr))*leadTrackEn)/
                    ((1./square(jetCaloEnErr)) + (1./square(leadTrackMomErr)));
    //retVal = TMath::Max(retVal, leadTrackEn*TMath::Sin(patTau.theta()) - patTau.pt());
    //retVal = TMath::Max(retVal, jetCaloEn*TMath::Sin(patTau.theta()) - patTau.pt());
    retVal = TMath::Max(retVal, combEn*TMath::Sin(patTau.theta()) - patTau.pt());
  }

  return retVal;
}

//
//-------------------------------------------------------------------------------
//

void printPatTau(const pat::Tau& patTau)
{
  std::cout << "<printPatTau>:" << std::endl;
  std::cout << " pt = " << patTau.pt() << "," 
	    << " eta = " << patTau.eta() << ", phi = " << patTau.phi() << std::endl;

  std::cout << "'signal' PFCandidates:" << std::endl;
  const reco::PFCandidateRefVector& signalPFCandidates = patTau.signalPFCands();
  for ( reco::PFCandidateRefVector::const_iterator signalPFCandidate = signalPFCandidates.begin();
	signalPFCandidate != signalPFCandidates.end(); ++signalPFCandidate ) {
    std::cout << getPFCandidateType(**signalPFCandidate) << " #" << signalPFCandidate->key() << ":"  
	      << " pt = " << (*signalPFCandidate)->pt() << "," 
	      << " eta = " << (*signalPFCandidate)->eta() << ", phi = " << (*signalPFCandidate)->phi() << std::endl;
  }
  
  std::cout << "'isolation' PFCandidates:" << std::endl;
  const reco::PFCandidateRefVector& isolationPFCandidates = patTau.isolationPFCands();
  for ( reco::PFCandidateRefVector::const_iterator isolationPFCandidate = isolationPFCandidates.begin();
	isolationPFCandidate != isolationPFCandidates.end(); ++isolationPFCandidate ) {
    std::cout << getPFCandidateType(**isolationPFCandidate) << " #" << isolationPFCandidate->key() << ":" 
	      << " pt = " << (*isolationPFCandidate)->pt() << "," 
	      << " eta = " << (*isolationPFCandidate)->eta() << ", phi = " << (*isolationPFCandidate)->phi() << std::endl;
  }

  if ( patTau.genJet() != 0 ) {
    std::cout << "genTauJet:" << std::endl;
    const reco::GenJet* genTauJet = patTau.genJet();
    std::cout << " pt = " << genTauJet->pt() << "," 
	      << " eta = " << genTauJet->eta() << ", phi = " << genTauJet->phi() << std::endl;
    std::vector<const reco::GenParticle*> genTauJetConstituents = genTauJet->getGenConstituents();
    for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
	  genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {
      std::cout << "PDG id. = " << (*genTauJetConstituent)->pdgId() << ": pt = " << (*genTauJetConstituent)->pt() << "," 
		<< " eta = " << (*genTauJetConstituent)->eta() << ", phi = " << (*genTauJetConstituent)->phi() << std::endl;
    }
  }

  std::cout << std::endl;
}

void printRecoPFJet(const reco::PFJet& recoPFJet, const reco::Vertex& vertex)
{
  std::cout << "<printRecoPFJet>:" << std::endl;
  std::cout << " pt = " << recoPFJet.pt() << "," 
	    << " eta = " << recoPFJet.eta() << ", phi = " << recoPFJet.phi() << std::endl;

  std::vector<reco::PFCandidatePtr> pfJetConstituents = recoPFJet.getPFConstituents();
  for ( std::vector<reco::PFCandidatePtr>::const_iterator pfJetConstituent = pfJetConstituents.begin();
	pfJetConstituent != pfJetConstituents.end(); ++pfJetConstituent ) {
    std::cout << getPFCandidateType(**pfJetConstituent) << " #" << pfJetConstituent->key() << ":" 
	      << " pt = " << (*pfJetConstituent)->pt() << "," 
	      << " eta = " << (*pfJetConstituent)->eta() << ", phi = " << (*pfJetConstituent)->phi();
    if ( (*pfJetConstituent)->charge() != 0. ) std::cout << ", charge = " << (*pfJetConstituent)->charge();
    std::cout << std::endl;
    const reco::TrackBaseRef track = getTrack(**pfJetConstituent);
    if ( track.isNonnull() ) {
      std::cout << "has Track: pt = " << track->pt() << " +/- " << track->ptError() << "," 
		<< " eta = " << track->eta() << ", phi = " << track->phi() << ","
		<< " charge = " << track->charge() << std::endl;
      std::cout << " chi2 = " << track->normalizedChi2() << std::endl;
      std::cout << " dIP = " << TMath::Abs(track->dxy(vertex.position())) << std::endl;
      std::cout << " dZ = " << TMath::Abs(track->dz(vertex.position())) << std::endl;
      std::cout << " vtxAssocWeight = " << vertex.trackWeight(track) << std::endl;
      std::cout << " numPxlHits = " << track->hitPattern().numberOfValidPixelHits() << std::endl;
      std::cout << " numTrkHits = " << track->hitPattern().numberOfValidHits() << std::endl;
    }
    std::cout << " ECAL Et: calibrated = " << (*pfJetConstituent)->ecalEnergy()*TMath::Sin((*pfJetConstituent)->theta()) << "," 
	      << " raw = " << (*pfJetConstituent)->rawEcalEnergy()*TMath::Sin((*pfJetConstituent)->theta()) << std::endl;
    std::cout << " HCAL Et: calibrated = " << (*pfJetConstituent)->hcalEnergy()*TMath::Sin((*pfJetConstituent)->theta()) << "," 
	      << " raw = " << (*pfJetConstituent)->rawHcalEnergy()*TMath::Sin((*pfJetConstituent)->theta()) << std::endl;
  }

  std::cout << std::endl;
}
