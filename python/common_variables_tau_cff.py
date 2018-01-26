import FWCore.ParameterSet.Config as cms

MuonIDFlags2016 = cms.PSet(

    Tight2016   = cms.string("isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10 && " +  
                             "numberOfMatchedStations > 1 && " +
                             "abs(dB) < 0.2 && "+
                             "track.hitPattern.numberOfValidPixelHits > 0 && "+
                             "track.hitPattern.trackerLayersWithMeasurement > 5"),
    #fabs(recoMu.muonBestTrack()->dz(vertex->position())) < 0.5 missing
    #"track.hitPattern.numberOfValidMuonHits > 0" always false in miniAOD

    IsolationValue2016 =  cms.string("(pfIsolationR03.sumChargedHadronPt + max(pfIsolationR03.sumNeutralHadronEt +"
                                     + "pfIsolationR03.sumPhotonEt - 0.5 * pfIsolationR03.sumPUPt, 0.0))/pt()"),

    Isolation2016 =   cms.string("(pfIsolationR03.sumChargedHadronPt + max(pfIsolationR03.sumNeutralHadronEt +"
                                     + "pfIsolationR03.sumPhotonEt - 0.5 * pfIsolationR03.sumPUPt, 0.0))/pt() < 0.1"),
                    

)


#
