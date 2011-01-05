import FWCore.ParameterSet.Config as cms

# If this is non-empty only specified numerators will be taken
#numerators_to_make = ['HPSloose', 'BgEst', 'TaNCloose']
numerators_to_make = ['HPSloose']

fakeRateDef = {
    'QCDdiJet1st' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 ' \
                       + ' && $ptIndex == 0 && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNC > 0.9 & abs($charge) == 1',
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },
    'QCDdiJet2nd' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 ' \
                       + ' && $ptIndex == 1 && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNC > 0.9 & abs($charge) == 1',
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },

    'PPmuX' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 && $byLeadTrackFinding > 0.5  && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNC > 0.9 & abs($charge) == 1',
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'samples' : {
            'data'       : { 'name' : 'data_ppmux' },
            'mc'         : { 'name' : 'ppmux_mc' },
            'mcPU156bx'  : { 'name' : 'ppmuxPU156bx_mc' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'ppmuxPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'PPmuXafterLoosePFIso' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 && $byLeadTrackFinding > 0.5 ' \
                       + ' && $ptSumLooseIsolation < 2.5 && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators' : {
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'samples' : {
            'data'       : { 'name' : 'data_ppmux' },
            'mc'         : { 'name' : 'ppmux_mc' },
            'mcPU156bx'  : { 'name' : 'ppmuxPU156bx_mc' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'ppmuxPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'WplusJets' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 && $byLeadTrackFinding > 0.5  && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators'  : {
            # NB tanc medium cut is reversed to remove ZTT contribution
            'BgEst' : '$byTaNC > 0.9 & $byTaNCmedium < 0.5',
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNC > 0.9 & abs($charge) == 1',
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'samples' : {
            'data'       : { 'name' : 'data_wjets' },
            'mc'         : { 'name' : 'wmunu_mc' },
            'mcPU156bx'  : { 'name' : 'wmunuPU156bx_mc' }
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'wjetsPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'WplusJetsAfterLoosePFIso' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 &&  $byLeadTrackFinding > 0.5  && $againstMuon > 0.5 && $againstElectron > 0.5' \
                       + ' && $ptSumLooseIsolation < 2.5',
        'numerators' : {
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'samples' : {
            'data'       : { 'name' : 'data_wjets' },
            'mc'         : { 'name' : 'wmunu_mc' },
            'mcPU156bx'  : { 'name' : 'wmunuPU156bx_mc' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'wjetsPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'Ztautau' : {
        'denominator' : '$pt > 20. && abs($eta) < 2.3 &&  $byLeadTrackFinding > 0.5  && $againstMuon > 0.5 && $againstElectron > 0.5',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCloose'  : '$byTaNCloose > 0.5 & abs($charge) == 1',
            'TaNCmedium' : '$byTaNCmedium > 0.5 & abs($charge) == 1',
            'TaNCtight'  : '$byTaNCtight > 0.5 & abs($charge) == 1',
            'HPSloose'   : '$byHPSloose > 0.5 & abs($charge) == 1',
            'HPSmedium'  : '$byHPSmedium > 0.5 & abs($charge) == 1',
            'HPStight'   : '$byHPStight > 0.5 & abs($charge) == 1',
        },
        'samples' : {
            'mc'         : { 'name' : 'ztautau_mc' },
            'mcPU156bx'  : { 'name' : 'zttPU156bx_mc' },
            'test'  : { 'name' : 'ztautau_test' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'zttPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    }
}
