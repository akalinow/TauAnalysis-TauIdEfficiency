import FWCore.ParameterSet.Config as cms

fakeRateDef = {
    'QCDdiJet1st' : {
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                       + ' && $ptIndex == 0',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNCvloose > 0.5',
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },
    'QCDdiJet2nd' : {
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                       + ' && $ptIndex == 1',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNCvloose > 0.5',
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
        },
        'hlt' : '$hltJet15Ubit > 0.5',
        'samples' : {
            'data'       : { 'name' : 'data_dijet' },
            'mc'         : { 'name' : 'qcddijet_mc' }
        }
    },
    
    'PPmuX' : {
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5',
        'numerators' : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNCvloose > 0.5',
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
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
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                       + ' && $ptSumLooseIsolation < 2.5',
        'numerators' : {
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
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
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5',
        'numerators'  : {
            'ewkTauId'   : '$byTaNCmedium > 0.5 && abs($charge) == 1',
            'TaNCvloose' : '$byTaNCvloose > 0.5',
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
        },
        'samples' : {
            'data'       : { 'name' : 'data_wjets' },
            'mc'         : { 'name' : 'wjets_mc' },
            'mcPU156bx'  : { 'name' : 'wjetsPU156bx_mc' }
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'wjetsPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'WplusJetsAfterLoosePFIso' : {
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5' \
                       + ' && $ptSumLooseIsolation < 2.5',
        'numerators' : {
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
        },
        'samples' : {
            'data'       : { 'name' : 'data_wjets' },
            'mc'         : { 'name' : 'wjets_mc' },
            'mcPU156bx'  : { 'name' : 'wjetsPU156bx_mc' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'wjetsPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    },
    'Ztautau' : {
        'denominator' : '$jetPt > 20. && abs($jetEta) < 2.3 && $probe > 0.5 && $byLeadTrackFinding > 0.5 && $byLeadTrackPtCut > 0.5',
        'numerators' : {
            'TaNCloose'  : '$byTaNCloose > 0.5',
            'TaNCmedium' : '$byTaNCmedium > 0.5',
            'TaNCtight'  : '$byTaNCtight > 0.5',
            'HPSloose'   : '$byHPSloose > 0.5',
            'HPSmedium'  : '$byHPSmedium > 0.5',
            'HPStight'   : '$byHPStight > 0.5'
        },
        'samples' : {
            'mc'         : { 'name' : 'ztautau_mc' },
            'mcPU156bx'  : { 'name' : 'zttPU156bx_mc' },
            ##'mcPU156bxReweighted' : {
            ##    'name'    : 'zttPU156bx_mc',
            ##    'options' : [ 'applyVertexMultiplicityReweighting' ]
            ##}
        }
    }
}