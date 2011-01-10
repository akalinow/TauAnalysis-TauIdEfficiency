# If this is non-empty only specified numerators will be taken
#numerators_to_make = ['HPSloose', 'BgEst', 'TaNCloose']
#fake_rates_to_make = ['ewkTauIdHPSloose', 'ewkTauIdTaNCloose', 'bgEstTemplateTaNCinverted']
#fake_rates_to_make = ['bgEstTemplateTaNCinverted']
fake_rates_to_make = ['ewkTauIdHPSloose']

def and_of(*cuts):
    " Produce the and of a set of cuts "
    return " && ".join(cuts)

_DENOM_COMMON = and_of(
    "$pt > 20.",
    "abs($eta) < 2.3",
    "$byLeadTrackFinding > 0.5",
    "$againstMuon > 0.5",
    "$againstElectron > 0.5",
)

_DENOM_TANC = and_of(_DENOM_COMMON, "$byLeadTrackPtCut > 0.5")

fakeRateDef = {
    'QCDdiJet1st' : {
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : and_of(_DENOM_COMMON, "$probe > 0.5", "$ptIndex == 0"),
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : and_of(_DENOM_TANC, "$probe > 0.5", "$ptIndex == 0"),
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : and_of(_DENOM_TANC, "$probe > 0.5", "$ptIndex == 0"),
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },
    'QCDdiJet2nd' : {
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : and_of(_DENOM_COMMON, "$probe > 0.5", "$ptIndex == 1"),
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : and_of(_DENOM_TANC, "$probe > 0.5", "$ptIndex == 1"),
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : and_of(_DENOM_TANC, "$probe > 0.5", "$ptIndex == 1"),
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },
    'PPmuX' : {
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : _DENOM_COMMON,
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
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
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : and_of(_DENOM_COMMON, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : and_of(_DENOM_TANC, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : and_of(_DENOM_TANC, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
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
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : _DENOM_COMMON,
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
            'bgEstTemplateTaNCinverted' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCmedium < 0.5 & $byTaNCvloose > 0.5 & abs($charge) == 1'
            },
        },
        'samples' : {
            'data'       : { 'name' : 'data_wjets' },
            'mc'         : { 'name' : 'wmunu_mc' },
            'mcPU156bx'  : { 'name' : 'wmunuPU156bx_mc' },
            'closure'  : { 'name' : 'wmunu_closure' }
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'wjetsPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'WplusJetsAfterLoosePFIso' : {
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : and_of(_DENOM_COMMON, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : and_of(_DENOM_TANC, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : and_of(_DENOM_TANC, "$ptSumLooseIsolation04 < 2.5"),
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
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
        'fake_rates' : {
            'ewkTauIdHPSloose' : {
                'denominator' : and_of(_DENOM_COMMON, "$genDecayMode > 1.5"),
                'numerator' : '$byHPSloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCloose' : {
                'denominator' : and_of(_DENOM_TANC, "$genDecayMode > 1.5"),
                'numerator' : '$byTaNCloose > 0.5 & abs($charge) == 1'
            },
            'ewkTauIdTaNCmedium' : {
                'denominator' : and_of(_DENOM_TANC, "$genDecayMode > 1.5"),
                'numerator' : '$byTaNCmedium > 0.5 & abs($charge) == 1'
            },
            'bgEstTemplateTaNCinverted' : {
                'denominator' : _DENOM_TANC,
                'numerator' : '$byTaNCmedium < 0.5 & $byTaNCvloose > 0.5 & abs($charge) == 1'
            },
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
