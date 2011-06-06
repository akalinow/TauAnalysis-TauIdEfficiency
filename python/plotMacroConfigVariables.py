import FWCore.ParameterSet.Config as cms

#//{name, numBins,min,max,expr} GetThem @ line 333 expr is how to get the variable from the pair


diTauMt = cms.PSet(
    xAxisTitle = cms.string("M_{T}^{#muMET} [GeV]"),
    rebinFactor = cms.int32(1),
    name = cms.string("diTauMt"),
    numBins = cms.int32(16),
    minh = cms.double(0.),
    maxh = cms.double(80.),
    expr = cms.string('mt1MET()')
    )

diTauVisMass = cms.PSet(
    xAxisTitle = cms.string("M_{vis}^{#mu#tau} [GeV]"),
    rebinFactor = cms.int32(1),
    name = cms.string("diTauVisMass"),
    numBins = cms.int32(36),
    minh = cms.double(20.),
    maxh = cms.double(200.),
    expr = cms.string('p4Vis.mass()')
    )
diTauVisMassFromJet = diTauVisMass.clone( expr = cms.string('(leg1().p4() + leg2().jetRef().p4()).mass()'), name = cms.string("diTauVisMassFromJet") )

diTauHt =cms.PSet(
    xAxisTitle = cms.string("P_{T}^{#mu} + P_{T}^{#tau} + MET [GeV]"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg1.pt() + leg2.pt() + met.pt()'),
    name = cms.string("diTauHt"),
    numBins = cms.int32(36),
    minh = cms.double(20.),
    maxh = cms.double(200.)
    )
diTauSVfitMass1 = diTauHt.clone( xAxisTitle = cms.string("M^{#tau#tau} [GeV]"), expr = cms.string("svFitSolution(\'psKine_MEt\').mass()"), name = cms.string("diTauSVfitMass1")  )
diTauSVfitMass2 = diTauSVfitMass1.clone( expr = cms.string("svFitSolution(\'psKine_MEt_ptBalance\').mass()"), name = cms.string("diTauSVfitMass2") )

muonPt = cms.PSet(
    xAxisTitle = cms.string("P_{T}^{#mu}"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg1().pt()'),
    name = cms.string("muonPt"),
    numBins = cms.int32(12),
    minh = cms.double( 0.),
    maxh = cms.double(60.),
    )
tauPt = muonPt.clone( xAxisTitle = cms.string("P_{T}^{#tau}"), expr = cms.string('leg2().pt()'), name = cms.string("tauPt") )   
tauJetPt = tauPt.clone( expr = cms.string('leg2().pfJetRef().pt()'), name = cms.string("tauJetPt") )

muonEta = cms.PSet(
    xAxisTitle = cms.string("#eta^{#mu}"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg1().eta()'),
    name = cms.string("muonEta"), 
    numBins = cms.int32(21),
    minh = cms.double(-2.1),
    maxh = cms.double(+2.1)
    )

tauEta = cms.PSet(
    xAxisTitle = cms.string("#eta^{#tau}"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg2().eta()'),
    name = cms.string("tauEta"),   
    numBins = cms.int32(23),
    minh = cms.double(-2.3),
    maxh = cms.double(+2.3)
    )
tauJetEta = tauEta.clone( expr = cms.string("leg2().pfJetRef().eta()"), name = cms.string("tauJetEta") )

tauNumChargedParticles = cms.PSet(
    xAxisTitle = cms.string("# #tau charged particles"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg2().signalPFChargedHadrCands().size() + leg2().isolationPFChargedHadrCands().size()'),
    name = cms.string("tauNumChargedParticles"), 
    numBins = cms.int32(15),
    minh = cms.double( -0.5),
    maxh = cms.double(+14.5)
    )

tauNumParticles = cms.PSet(
    xAxisTitle = cms.string("# #tau PF Candidates"),
    rebinFactor = cms.int32(1),
    expr = cms.string('leg2().signalPFCands().size() + leg2().isolationPFCands().size()'),
    name = cms.string("tauNumParticles"),       
    numBins = cms.int32(25),
    minh = cms.double( -0.5),
    maxh = cms.double(+24.5)
    )

tauJetWidth = cms.PSet(
    xAxisTitle = cms.string("#tau jet width"),
    rebinFactor = cms.int32(1),
    expr = cms.string('? (leg2().pfJetRef().etaetaMoment() + leg2().pfJetRef().phiphiMoment()) > 0. ? sqrt(leg2().pfJetRef().etaetaMoment() + leg2().pfJetRef().phiphiMoment()) : 0.'),
    name = cms.string("tauJetWidth"), 
    numBins = cms.int32(20),
    minh = cms.double(0.),
    maxh = cms.double(0.50)
    )

