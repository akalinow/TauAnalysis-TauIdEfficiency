import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# select pp collision events in "good" runs
# ( following recommendations given at https://twiki.cern.ch/twiki/bin/view/CMS/Collisions2010Recipes )
#--------------------------------------------------------------------------------

# veto events not recorded during periods of "physics declared"
from HLTrigger.special.hltPhysicsDeclared_cfi import *
hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

#veto events in which not all subdetectors were "on"
#
# NOTE: the DCS partitions for which "on"/"off" status is recorded separately are:
#
#  "EBp", "EBm", "EEp", "EEm", "HBHEa", "HBHEb", "HBHEc", "HF", "HO", "RPC"
#  "DT0", "DTp", "DTm", "CSCp", "CSCm", "CASTOR", "TIBTID", "TOB", "TECp", "TECm"
#  "BPIX", "FPIX", "ESp", "ESm"
#
from DPGAnalysis.Skims.DetStatus_cfi import *
dcsstatus.DetectorType = cms.vstring(
    'EBp', 'EBm', 'EEp', 'EEm',
    ##'ESp', 'ESm',
    'HBHEa', 'HBHEb', 'HBHEc',
    ##'HF', 'HO',
    'DT0', 'DTp', 'DTm', 'CSCp', 'CSCm',
    'TIBTID', 'TOB', 'TECp', 'TECm',
    'BPIX', 'FPIX'
)
dcsstatus.ApplyFilter = cms.bool(True)
dcsstatus.DebugOn = cms.untracked.bool(False)
dcsstatus.AndOr = cms.bool(True)

# veto "scraping" beam events
# (identified by small percentage "good" tracks fitted to primary event vertex)
scrapingBeamsFilter = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

# veto events without a "good" primary vertex
# ( see http://indico.cern.ch/getFile.py/access?subContId=0&contribId=2&resId=0&materialId=slides&confId=83613 )
primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(15),	
    maxd0 = cms.double(2)
)

#--------------------------------------------------------------------------------
# add MET cleaning filters
# (cf. https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters)
#--------------------------------------------------------------------------------

# veto events with halo muons
##from RecoMET.METAnalyzers.CSCHaloFilter_cfi import CSCTightHaloFilter

# veto events with significant RBX/HPD noise activity
# ( see https://twiki.cern.ch/twiki/bin/view/CMS/HcalDPGAnomalousSignals )
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import HBHENoiseFilter

# veto events in which HCAL laser calibration fired
##from RecoMET.METFilters.hcalLaserEventFilter_cfi import hcalLaserEventFilter
##hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
##hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

# veto events with significant energy deposits close to dead/masked ECAL cells
##from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import EcalDeadCellTriggerPrimitiveFilter
##EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

# veto events with (suspected) tracking failures
##goodVerticesForTrackingFailureFilter = cms.EDFilter("VertexSelector",
##  filter = cms.bool(False),
##  src = cms.InputTag("offlinePrimaryVertices"),
##  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
##)
##from RecoMET.METFilters.trackingFailureFilter_cfi import trackingFailureFilter
##trackingFailureFilter.VertexSource = cms.InputTag('goodVerticesForTrackingFailureFilter')

# veto events with anomalous supercrystals in ECAL endcap
##from RecoMET.METFilters.eeBadScFilter_cfi import eeBadScFilter

dataQualityFilters = cms.Sequence(
    hltPhysicsDeclared
   * dcsstatus
   * scrapingBeamsFilter
   * primaryVertexFilter
   ##* CSCTightHaloFilter
   * HBHENoiseFilter
   ##* hcalLaserEventFilter
   ##* goodVerticesForTrackingFailureFilter * trackingFailureFilter
   ##* eeBadScFilter
)    
