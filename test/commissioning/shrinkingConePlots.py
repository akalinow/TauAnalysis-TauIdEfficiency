#!/usr/bin/env python

import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples_cache as samples
import os
import copy
import sys
from optparse import OptionParser

# Build the plot manager.  The plot manager keeps track of all the samples
# and ensures they are correctly normalized w.r.t. luminosity.  See 
# samples.py for available samples.
plotter = PlotManager()

# Add each sample we want to plot/compare
# Uncomment to add QCD
plotter.add_sample(samples.qcd_mc, "QCD MC", **style.QCD_MC_STYLE_HIST)

#plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)

plotter.add_sample(samples.data, "Data (7 TeV)", **style.DATA_STYLE)

# Normalize everything to the data luminosity
plotter.set_integrated_lumi(samples.data.effective_luminosity())

# Build the ntuple maanger
ntuple_manager = samples.data.build_ntuple_manager("tauIdEffNtuple")

shrinking_ntuple = ntuple_manager.get_ntuple(
    "patPFTausDijetTagAndProbeShrinkingCone")
fixed_ntuple = ntuple_manager.get_ntuple(
    "patPFTausDijetTagAndProbeFixedCone")

hlt = ntuple_manager.get_ntuple("TriggerResults")

# Make some plots
canvas = ROOT.TCanvas("blah", "blah", 500, 500)

# Plot # different triggers
#    trigger_results_expr =  hlt.expr('$hltJet15U')
#
#    trigger_results = plotter.distribution(
#        expression=hlt.expr('$hltJet15U'),
#        selection=hlt.expr('1'), # no selection
#        binning = (5, -2.5, 2.5),
#        x_axis_title = "HLT_Jet15U Result",
#        y_min = 1, logy=True
#    )
#    canvas.SaveAs("plots/hltJet15U_result.png")
#    canvas.SaveAs("plots/hltJet15U_result.pdf")


# Basic requirement HLT + Probe object
# N.B. currently disabled, no HLT info in ntuples!
base_selection = hlt.expr('$hltJet15U > 0.5') & shrinking_ntuple.expr('$probe > 0.5') 
#base_selection = shrinking_ntuple.expr('1')

# this actually, slows things down.  Lame!
#samples.data.set_selection( hlt.expr('$hltJet15U > 0.5'), "/afs/cern.ch/user/f/friis/entrylists.root")
#samples.qcd_mc.set_selection( hlt.expr('$hltJet15U > 0.5'), "/afs/cern.ch/user/f/friis/entrylists.root")
#samples.minbias_mc.set_selection( hlt.expr('$hltJet15U > 0.5'), "/afs/cern.ch/user/f/friis/entrylists.root")

# Define some common expressions
expr_jetPt = shrinking_ntuple.expr('$jetPt')
expr_jetEta = shrinking_ntuple.expr('$jetEta')
expr_abs_jetEta = shrinking_ntuple.expr('abs($jetEta)')


# Define some common binnings
#pt_binning_fine = (50, 0, 100)
pt_binning_fine = (0, 5, 10, 15, 20, 25, 35, 45, 60, 80, 100)
eta_binning_fine = (50, -2.5, 2.5)
phi_binning_fine = (50, -3.14, 3.14)
decay_mode_binning = (25, -0.5, 24.5)
discriminator_binning = (200, -0.5, 1.5)

eta_acceptance_cut = (expr_abs_jetEta < 2.5)
min_pt_cut = (expr_jetPt > 10)
basic_kinematic_cut = min_pt_cut & eta_acceptance_cut 
lead_pion_selection = shrinking_ntuple.expr('$byLeadPionPtCut')
lead_track_finding = shrinking_ntuple.expr('$byLeadTrackFinding')

distributions = {
    # Trigger result
    'hltJet15U' : {
        'expression': hlt.expr('$hltJet15U'),
        'selection': hlt.expr("$hltJet15U > -1000"),  # ?
        'binning': (5, -2.5, 2.5),
        'x_axis_title': "HLT_Jet15U Result",
        'y_min': 1, 'logy': True
    },
    'testLeadTrackCorner': {
        'expression': shrinking_ntuple.expr('$leadChargedParticlePt'),
        'selection' : lead_track_finding & \
        shrinking_ntuple.expr('$numChargedParticlesIsoCone == 1') & \
        shrinking_ntuple.expr('$numChargedParticlesSignalCone == 0') & \
        shrinking_ntuple.expr('$leadChargedParticlePt > 0') &\
        (lead_pion_selection < 1),
        'binning': (100, 0, 25),
        'x_axis_title': "Difference between lead track pt and iso sum charged",
    },
    'testLeadTrackCornerFixedDiff': {
        'expression': fixed_ntuple.expr('$leadChargedParticlePt') -\
        fixed_ntuple.expr('$ptSumChargedParticlesIsoCone'),
        'selection' : lead_track_finding & \
        fixed_ntuple.expr('$numChargedParticlesIsoCone == 1') & \
        fixed_ntuple.expr('$numChargedParticlesSignalCone == 0') & \
        fixed_ntuple.expr('$leadChargedParticlePt > 0') &\
        lead_pion_selection,
        'binning': (150, -15, 15),
        'x_axis_title': "Difference between lead track pt and iso sum charged",
    },
    'testLeadTrackCornerFixed': {
        'expression': fixed_ntuple.expr('$leadChargedParticlePt'),
        'selection' : lead_track_finding & \
        fixed_ntuple.expr('$numChargedParticlesIsoCone == 1') & \
        fixed_ntuple.expr('$numChargedParticlesSignalCone == 0') & \
        fixed_ntuple.expr('$leadChargedParticlePt > 0') &\
        lead_pion_selection,
        'binning': (150, 0, 15),
        'x_axis_title': "Difference between lead track pt and iso sum charged",
    },
    'testLeadTrackCornerFixedNoLeadPion': {
        'expression': fixed_ntuple.expr('$leadChargedParticlePt'),
        'selection' : lead_track_finding & \
        fixed_ntuple.expr('$numChargedParticlesIsoCone == 1') & \
        fixed_ntuple.expr('$numChargedParticlesSignalCone == 0') & \
        (lead_pion_selection < 1) &\
        fixed_ntuple.expr('$leadChargedParticlePt > 0'),
        'binning': (150, 0, 15),
        'x_axis_title': "Difference between lead track pt and iso sum charged",
    },
    'jetPt' : {
        'expression': expr_jetPt,
        'selection': base_selection & eta_acceptance_cut,
        'binning': pt_binning_fine, 'y_min': 1e2, 'logy': True,
        'extra_labels': [style.ETA_CUT_LABEL_UPPER_LEFT],
        'x_axis_title': "Jet P_{T} [GeV/c]",
    },
    'jetPt_dm0' : {
        'expression': expr_jetPt,
        'selection': base_selection & eta_acceptance_cut & shrinking_ntuple.expr('$decayMode == 0'),
        'binning': pt_binning_fine, 'y_min': 1, 'logy': True,
        'extra_labels': [style.ETA_CUT_LABEL_UPPER_LEFT],
        'x_axis_title': "Jet P_{T} for Decay Mode 0 [GeV/c]",
    },
    'jetPt_dm0_leadingPion' : {
        'expression': expr_jetPt,
        'selection': base_selection & eta_acceptance_cut & shrinking_ntuple.expr('$decayMode == 0') &\
        shrinking_ntuple.expr('$byLeadPionPtCut'),
        'binning': pt_binning_fine, 'y_min': 1, 'logy': True,
        'extra_labels': [style.ETA_CUT_LABEL_UPPER_LEFT],
        'x_axis_title': "Jet P_{T} for Decay Mode 0 [GeV/c]",
    },
    'jetEta' : {
        'expression': expr_jetEta,
        'selection': base_selection & min_pt_cut,
        'binning': eta_binning_fine, 'logy': False,
        'extra_labels': [style.PT_CUT_LABEL_UPPER_LEFT],
        #'y_min': 1e3, 
        'x_axis_title': "Jet #eta",
    },
    'jetWidth' : {
        'expression': shrinking_ntuple.expr('$jetWidth'),
        'selection' : base_selection & basic_kinematic_cut,
        'binning' : (50, 0, 0.5), 'logy':True,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'x_axis_title': "Jet Width"
    },
    'jetPhi' : {
        'expression': shrinking_ntuple.expr('$jetPhi'),
        'selection': base_selection & basic_kinematic_cut,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': phi_binning_fine, 'logy': False,
        'x_axis_title': "Jet #phi",
    },
    'decayMode' : {
        'expression': shrinking_ntuple.expr('$decayMode'),
        #'selection':base_selection & basic_kinematic_cut & lead_pion_selection,
        'selection':base_selection & basic_kinematic_cut & shrinking_ntuple.expr('$byLeadTrackPtCut'),
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': decay_mode_binning, 'logy': False,
        'x_axis_title': "Decay Mode",
    },

    'decayMode_pt20' : {
        'expression': shrinking_ntuple.expr('$decayMode'),
        'selection':base_selection & basic_kinematic_cut & lead_pion_selection &\
        shrinking_ntuple.expr('$jetPt > 20'),
        'extra_labels': [style.PT_20_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': decay_mode_binning, 'logy': False,
        'x_axis_title': "Decay Mode",
    },

    'decayMode_pt40' : {
        'expression': shrinking_ntuple.expr('$decayMode'),
        'selection':base_selection & basic_kinematic_cut & lead_pion_selection &\
        shrinking_ntuple.expr('$jetPt > 40'),
        'extra_labels': [style.PT_40_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': decay_mode_binning, 'logy': False,
        'x_axis_title': "Decay Mode",
    },

    'byTaNC' : {
        'expression': shrinking_ntuple.expr('$byTaNC'),
        'selection': base_selection & basic_kinematic_cut & lead_pion_selection,
        'binning': discriminator_binning, 'logy': True, 'y_min': 1,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'x_axis_title': 'TaNC output',
    },

    'byTaNC_pt40' : {
        'expression': shrinking_ntuple.expr('$byTaNC'),
        'selection': base_selection & basic_kinematic_cut & lead_pion_selection &\
        shrinking_ntuple.expr('$jetPt > 40'),
        'extra_labels': [style.PT_40_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': discriminator_binning, 'logy': True, 'y_min': 1,
        'x_axis_title': 'TaNC output',
    },

    'leadChargedParticlePt' : {
        'expression': shrinking_ntuple.expr('$leadChargedParticlePt'),
        'selection' : base_selection & basic_kinematic_cut & \
        shrinking_ntuple.expr('$byLeadTrackFinding'),
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (50, 0, 50), 'logy': True, 'y_min': 1,
        'x_axis_title' : "P_{T} of leading charged particle"
    },

    'ptSumChargedParticlesIsoCone' : { 
        'expression' : shrinking_ntuple.expr('$ptSumChargedParticlesIsoCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & lead_track_finding, 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (50, 0, 50), 'logy': True, 'y_min': 1e-1,
        'x_axis_title' : "Isolation charged sum P_{T} [GeV/c]"
    },

    'ptSumPhotonsIsoCone' : { 
        'expression' : shrinking_ntuple.expr('$ptSumPhotonsIsoCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & lead_track_finding, 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (50, 0, 50), 'logy': True, 'y_min': 1e-1,
        'x_axis_title' : "Isolation photons sum P_{T} [GeV/c]"
    },

    'numChargedParticlesIsoCone' : { 
        'expression' : shrinking_ntuple.expr('$numChargedParticlesIsoCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & lead_track_finding, 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (21, -0.5, 20.5), 'logy': True, 'y_min': 1,
        'x_axis_title' : "Number of charged hadrons in isolation annulus"
    },

    'numPhotonsIsoCone' : { 
        'expression' : shrinking_ntuple.expr('$numPhotonsIsoCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & lead_track_finding, 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (21, -0.5, 20.5), 'logy': True, 'y_min': 1,
        'x_axis_title' : "Number of photons in isolation annulus"
    },

    'numChargedParticlesSignalCone' : { 
        'expression' : shrinking_ntuple.expr('$numChargedParticlesSignalCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & shrinking_ntuple.expr('$byLeadTrackPtCut'), 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (21, -0.5, 20.5), 'logy': True, 'y_min': 1,
        'x_axis_title' : "Number of charged hadrons in signal cone"
    },

    'numChargedParticlesSignalConeFixed' : { 
        'expression' : fixed_ntuple.expr('$numChargedParticlesSignalCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & fixed_ntuple.expr('$byLeadTrackPtCut'), 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (21, -0.5, 20.5), 'logy': True, 'y_min': 1,
        'x_axis_title' : "Number of charged hadrons in signal cone"
    },

    'numPhotonsSignalCone' : { 
        'expression' : shrinking_ntuple.expr('$numPhotonsSignalCone'),
        # Require lead pion selection to ensure iso cone is built
        'selection' : base_selection & basic_kinematic_cut & lead_track_finding, 
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
        'binning': (21, -0.5, 20.5), 'logy': True, 'y_min': 1,
        'x_axis_title' : "Number of photons in signal cone"
    },
}

# Add in TaNC distributions for each neural net
for decay_mode in [0, 1, 2, 10, 11]:
    dist_name = 'byTaNC_dm%i' % decay_mode
    distributions[dist_name] = copy.deepcopy(distributions['byTaNC'])
    new_selection = distributions[dist_name]['selection'] & \
            shrinking_ntuple.expr('$decayMode == %i' % decay_mode)
    distributions[dist_name]['selection'] = new_selection

######################################################
####      Plot efficiencies                       ####
######################################################

denominator = base_selection & basic_kinematic_cut

efficiencies = {
    'LeadTrackFinding': {
        'numerator': shrinking_ntuple.expr('$byLeadTrackFinding'),
        'denominator': denominator, 'logy': True
    },
    'LeadTrackPtCut': {
        'numerator': shrinking_ntuple.expr('$byLeadTrackPtCut'),
        'denominator': denominator, 'logy': True
    },
    'LeadPionPtCut': {
        'numerator': shrinking_ntuple.expr('$byLeadPionPtCut'),
        'denominator': denominator, 'logy': True
    },
    'EcalIso': {
        'numerator': shrinking_ntuple.expr('$byEcalIsolationUsingLeadingPion') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'TrackIso': {
        'numerator': shrinking_ntuple.expr('$byTrackIsolationUsingLeadingPion') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'CombinedIso': {
        'numerator': shrinking_ntuple.expr('$byIsolationUsingLeadingPion') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'TaNCTenth': {
        'numerator': shrinking_ntuple.expr('$byTaNCfrTenthPercent') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'TaNCQuarter': {
        'numerator': shrinking_ntuple.expr('$byTaNCfrQuarterPercent') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'TaNCHalf': {
        'numerator': shrinking_ntuple.expr('$byTaNCfrHalfPercent') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
    'TaNCOne': {
        'numerator': shrinking_ntuple.expr('$byTaNCfrOnePercent') & lead_pion_selection,
        'denominator': denominator, 'logy': True
    },
}

# Update the numerator for each efficiency to ensure it is a subset of the
# denominator
for eff_name, eff_info in efficiencies.iteritems(): 
    eff_info['numerator'] = eff_info['numerator'] & eff_info['denominator']

# Define quantities to parameterize efficiency
efficiency_versus = {
    'Pt' : {
        'expression' : expr_jetPt,
        'x_axis_title': "Jet P_{T} [GeV/c]",
        'binning': pt_binning_fine,
        'y_min' : 1e-4, 'y_max' : 5,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
    },
    'Eta' : {
        'expression' : expr_jetEta,
        'x_axis_title': "Jet #eta",
        'binning': eta_binning_fine,
        'y_min' : 1e-4, 'y_max' : 5,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
    },
    'Width' : {
        'expression' : shrinking_ntuple.expr('$jetWidth'),
        'x_axis_title': "Jet width",
        'binning': (50, 0, 0.3),
        'y_min' : 1e-4, 'y_max' : 5,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
    },
    'Phi' : {
        'expression' : shrinking_ntuple.expr('$jetPhi'),
        'x_axis_title': "Jet #phi",
        'binning': phi_binning_fine,
        'y_min' : 1e-4, 'y_max' : 5,
        'extra_labels': [style.PT_ETA_CUT_LABEL_UPPER_LEFT],
    },
}

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)

    if not os.path.isdir('plots'):
        os.mkdir('plots')
        
    parser = OptionParser()

    parser.add_option("-p", "--plot", dest="to_plot", 
                      action="append", default=[], 
                      help="Appendable list of distributions to plot")
    parser.add_option("-v", "--versus", dest="versus", 
                      action="append", default=[], 
                      help="Appendable list of dependent variables for the efficiencies.  Default is all")
    parser.add_option("-a", "--all", default=False, action="store_true")
    
    (opts, args) = parser.parse_args()

    if not opts.versus:
        opts.versus = ["Pt", "Eta", "Width", "Phi"]

    print "Building distributions:", " ".join(opts.to_plot)
    print "For efficiencies, dependent variables are:", " ".join(opts.versus)

    # Print available plots
    if not opts.to_plot and not opts.all:
        print " Available distributions:"
        dists = distributions.keys()
        dists.sort()
        for distribution in dists:
            print " *", distribution
        print " Available efficiencies:"
        dists = efficiencies.keys()
        dists.sort()
        for distribution in dists:
            print " *", distribution
        sys.exit(0)


    # Plot all distributions
    for dist_name, dist_info in distributions.iteritems():
        if dist_name not in opts.to_plot and not opts.all:
            continue
        result = plotter.distribution(**dist_info)
        result['legend'].make_legend().Draw()
        canvas.SaveAs("plots/shrinkingCone_dist_%s.png" % dist_name)
        canvas.SaveAs("plots/shrinkingCone_dist_%s.pdf" % dist_name)

    # Change the style of the QCD from filled histogram to dots
    # name mc_qcd is defined in samples.py
    plotter.update_style("mc_qcd", **style.QCD_MC_STYLE_DOTS)

    # Plot all efficiencies
    for eff_vs, eff_vs_info in efficiency_versus.iteritems():
        if eff_vs not in opts.versus and not opts.all:
            continue
        for eff_name, eff_info in efficiencies.iteritems():
            if eff_name not in opts.to_plot and not opts.all:
                continue
            # Combine the information about the x_axis and the
            # numerator and denominator selections
            plot_info_dict = {}
            plot_info_dict.update(eff_vs_info)
            plot_info_dict.update(eff_info)
            # Make the plot
            eff_result = plotter.efficiency(**plot_info_dict)
            # Add a legend
            eff_result['legend'].make_legend().Draw()
            # Save plots
            canvas.SaveAs("plots/shrinkingCone_%s_eff_vs_%s.png" 
                          % (eff_name, eff_vs))
            canvas.SaveAs("plots/shrinkingCone_%s_eff_vs_%s.pdf" 
                          % (eff_name, eff_vs))

