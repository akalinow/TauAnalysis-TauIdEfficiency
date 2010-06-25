#!/usr/bin/env python
from __future__ import with_statement
import ROOT

from TauAnalysis.TauIdEfficiency.ntauples.PlotManager import PlotManager
import TauAnalysis.TauIdEfficiency.ntauples.styles as style

# Defintion of input files.
import samples_cache as samples
import os, sys

class PlotSession:
    def __init__( self, samplesToPlot=["data","mc"], name = None, baseDir = "plots",
                  formats = ["png", "pdf"] ):
        self.formats = formats     

        self.__initSamples(samplesToPlot)
        
        sessions = self.__possibleSessions(baseDir)
        while not (name in sessions):
            print "(Unkown) name. Please Choose:"
            name = self.__makeMenu(sessions)

        self.selection = sessions[name]["selection"]
        self.selectionDiscription = sessions[name]["description"]
        self.subdir = sessions[name]["subdir"]
        self.description = sessions[name]["description"]
        self.labels=[style.CMS_PRELIMINARY_UPPER_LEFT, style.LUMI_LABEL_UPPER_LEFT,
                     style.PT_ETA_CUT_LABEL_UPPER_LEFT, sessions[name]["label"]]
        
        if not os.path.isdir(self.subdir):
            os.makedirs(self.subdir)
            
        self.plots = {}
        self.canvas = ROOT.TCanvas("blah", "blah", 500, 500)

    def drawAll(self):
        self.drawDecayMode()
        self.drawKineticPlots()
        self.drawMainTrackPlots()
        self.drawMassPlots()
        self.drawOutlierSumKineticPlots()
        self.drawOutlierCountPlots()
        self.drawOutlierMassPlots()
        self.drawOutlierSingleKineticPlots()
        self.drawPi0KineticPlots()
        self.drawTrackKineticPlots()

    def drawEssential(self):
        self.drawDecayMode()
        self.drawKineticPlots()
        self.drawMainTrackPlots()

    def drawComplementaryPlots(self):
        self.drawMassPlots()
        self.drawPi0KineticPlots()
        self.drawTrackKineticPlots()

    def drawOutlierPlots(self):
        self.drawOutlierCountPlots()
        self.drawOutlierSumKineticPlots()
        self.drawOutlierMassPlots()
        self.drawOutlierSingleKineticPlots()

    def drawDecayMode(self):
        plotName = "DecayMode"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCDecayMode'),
            selection= self.selection,
            binning = (21, -0.5, 20.5),
            x_axis_title = "#tau decay mode",
            #y_min = 1e-2,
            logy=False,
            labels = self.labels
            )
        latexTable = self.__makeDecayModeTable(plotName)
        with open("%s/%s.tex"%(self.subdir, plotName), "w") as latexFile:
            latexFile.write(latexTable)
       
        self.__printPlots(plotName)

    def drawTaNCOutputPlots(self):
        plotName = "TaNCOutput"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$byTaNC'),
            selection= self.selection,
            binning = (100, 0.0, 1.0),
            x_axis_title = "TaNC NN output",
            #y_min = 1e-2,
            logy=False,
            labels = self.labels
            )
        self.__printPlots(plotName)

    def drawKineticPlots(self):
        plotName = "TauPt"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCPt'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (100, 0, 100),
            x_axis_title = "#tau p_{T} [GeV/c]",
            #y_min = 1e-2,
            logy=True,
            labels = self.labels
            )
        self.__printPlots(plotName)
        
        plotName = "TauEta"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCEta'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (25, 0, 2.6),
            x_axis_title = "|#eta|",
            #y_min = 1e-2,
            logy=False,
            labels = self.labels
            )            
        self.__printPlots(plotName)

    def drawMainTrackPlots(self):
        plotName = "MainTrackPt"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCMainTrackPt'),
            selection= self.selection & self.__getDecayModeSelection(plotName),            
            binning = (100, 0, 70),
            x_axis_title = "p_{T} of main Track",
            #       y_min = 1e-2,
            logy=True,
            labels = self.labels
            )
        self.__printPlots(plotName)
        
        plotName = "MainTrackAngle"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCMainTrackAngle'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (100, 0, 0.3),
            x_axis_title = "#DeltaR of main Track",
            #        y_min = 1e-2,
            logy=False,
            labels = self.labels
            )                
        self.__printPlots(plotName)       

    def drawMassPlots(self):
        plotName = "InvariantMassOfSignal"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCInvariantMassOfSignal'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (100, 0, 4.0),
            x_axis_title = "inv. mass of signal particles",
            #y_min = 1e-2,
            logy=False,
            labels = self.labels
            )
        self.__printPlots(plotName)

        for dalitz in range(1,3):
            plotName = "Dalitz%s"%dalitz
            self.plots[plotName] = self.plotter.distribution(
                expression=self.main.expr("$TaNCDalitz%s"%(dalitz-1)),
                selection= self.selection & self.__getDecayModeSelection(plotName),
                binning = (100, 0, 4),
                x_axis_title = "Dalitz%s"%dalitz,
                #y_min = 1e-2,
                logy=False,
                labels = self.labels
                )
            self.__printPlots(plotName)
            
    def drawOutlierSumKineticPlots(self):
        for charge in ["","Charged","Neutral"]:
            plotName = "%sOutlierSumPt"%charge
            self.plots[plotName] = self.plotter.distribution(
                expression=self.main.expr('$TaNC%sOutlierSumPt'%(charge)),
                selection= self.selection & self.__getDecayModeSelection(plotName),
                binning = (100, 0, 70),
                x_axis_title = "Sum %s Outlier p_{T} [GeV/c]"%(charge),
                #y_min = 1e-2,
                logy=True,
                labels = self.labels
                )
            self.__printPlots(plotName)
                
    def drawOutlierCountPlots(self):
        for charge in ["","Charged"]:
            plotName = "OutlierN%s"%charge
            self.plots[plotName] = self.plotter.distribution(
                expression=self.main.expr('$TaNCOutlierN%s'%(charge)),
                selection= self.selection & self.__getDecayModeSelection(plotName),
                binning = (41, -0.5, 40.5),
                x_axis_title = "%s Outlier Count"%(charge),
                #y_min = 1e-2,
                logy=True,
                labels = self.labels
                )
            self.__printPlots(plotName)
            
    def drawOutlierMassPlots(self):
        plotName = "OutlierMass"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCOutlierMass'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (100, 0, 70),
            x_axis_title = "Inv. Mass of Outliers [GeV/c]",
            #y_min = 1e-2,
            logy=False,
            labels = self.labels
            )
        self.__printPlots(plotName)
                
    def drawOutlierSingleKineticPlots(self, charges = ["","Charged","Neutral"], numOutlier=range(1,5)):
        for charge in charges:
            for i in numOutlier:
                self.__drawPt(charge, "Outlier", i)
                self.__drawAngle(charge, "Outlier", i)
                
    def drawPi0KineticPlots(self):
        for i in range(1,3):
            self.__drawPt("","PiZero", i)
            self.__drawAngle("","PiZero", i)
            
    def drawTrackKineticPlots(self):
        for i in range(1,3):
            self.__drawPt("","Track", i)
            self.__drawAngle("","Track", i)
            
    def __drawPt(self, name, objects, number):
        plotName = "%s%sPt%s"%(name,objects,number)
        try:
            self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr("$TaNC%s%sPt%s"%(name,objects,(number-1))),
                selection= self.selection & self.__getDecayModeSelection(plotName),
                binning = (100, 0, 25),
                x_axis_title = "p_{T} of %s %s No. %s [GeV/c]"%(name,objects,number),
                #y_min = 1e-2,
                logy=True,
                labels = self.labels
                )            
            self.__printPlots(plotName)
        except StandardError, err:
            print "TaNC%s%sPt%s not found: %s"%(name,objects,number,err)
            
    def __drawAngle(self, name, objects, number):
        plotName = "%s%sAngle%s"%(name,objects,number)
        try:
            self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr("$TaNC%s%sAngle%s"%(name,objects,(number-1))),
                selection= self.selection & self.__getDecayModeSelection(plotName),
                binning = (100, 0, 1.5),
                x_axis_title = "#DeltaR %s %s No. %s"%(name,objects,number),
                #y_min = 1e-2,
                logy=False,
                labels = self.labels
                )            
            self.__printPlots(plotName)
        except StandardError, err:
            print "TaNC%s%sAngle%s not found: %s"%(name,objects,number,err)

    def __getFileNames(self, name):
        for f in self.formats:
            yield "%s/%s.%s"%(self.subdir, name, f)

    def __printPlots(self, name):
        self.plots[name]['legend'].make_legend().Draw() 
        for fileName in self.__getFileNames(name):
            self.canvas.SaveAs(fileName)

    def __getDecayModeSelection(self, plotName):
        modes =[]
        allPlots = ["TauPt","TauEta","DecayMode", "MainTrackPt","OutlierMass"]
        for charge in ["","Charged","Neutral"]:
            for i in range(1,5):
                allPlots.append("%sOutlierPt%s"%(charge, i))
                allPlots.append("%sOutlierAngle%s"%(charge, i))
                allPlots.append("%sOutlierSumPt"%(charge))
                allPlots.append("OutlierN%s"%(charge))
        selectedPlots = {
            "MainTrackAngle":  [1,2,10,11],
            "InvariantMassOfSignal":  [1,2,10,11],
            "Dalitz1": [2,10,11],
            "Dalitz2": [2,10,11],
            "PiZeroPt1":[1,2,11],
            "PiZeroPt2":[2],
            "PiZeroAngle1":[1,2,11],
            "PiZeroAngle2":[2],
            "TrackPt1":[10,11],
            "TrackPt2":[10,11],
            "TrackAngle1":[10,11],
            "TrackAngle2":[10,11],
            }
        if plotName in allPlots:
            modes = [0,1,2,10,11]
        elif plotName in selectedPlots:
            modes = selectedPlots[plotName]

        if modes == []: return self.main.expr("0")
        
        return self.main.expr(" | ".join(["$TaNCDecayMode == %s"%mode for mode in modes]))

    def __makeDecayModeTable(self, plotName):
        from math import sqrt
        latexModes={ 0 :"$\\nu_{\\tau}h^{\\pm}(\\pi^{\\pm})$",
                     1 :"$\\nu_{\\tau}h^{\\pm}(\\pi^{\\pm}\\pi^{0})$",
                     2 :"$\\nu_{\\tau}h^{\\pm}(\\pi^{\\pm}\\pi^{0}\\pi^{0})$",
                     10:"$\\nu_{\\tau}h^{\\pm}h^{\\mp}h^{\\pm}(\\pi^{\\pm}\\pi^{\\mp}\\pi^{\\pm})$",
                     11:"$\\nu_{\\tau}h^{\\pm}h^{\\mp}h^{\\pm}(\\pi^{\\pm}\\pi^{\\mp}\\pi^{\\pm}\\pi^{0})$",
            }
        sampleNames = [ "data", "mc_qcd"]
        
        result = "%%Table for selection '%s': "%self.selectionDiscription
        rows = {}
        for mode in latexModes:
            rows[mode]=[ latexModes[mode]]
        for sampleName in sampleNames:
            result += " %s "%sampleName
            N = 0.
            k = {}
            for mode in latexModes:
                hist = self.plots[plotName]['samples'][sampleName]['plot']
                binNr = hist.GetXaxis().FindBin(mode)
                k[mode] = float(hist.GetBinContent(binNr))
                N += k[mode]
                
            for mode in latexModes:
                #binomial errors for now...
                print sampleName, mode, k[mode], N
                rows[mode].append("%.1f \pm %.2f"%(k[mode]/N*100, 1/N*sqrt(k[mode]*(1-k[mode]/N) )*100))
            print rows
        result +="\n"
        for mode in rows:
            result += " & ".join(rows[mode])+" \\\\\n"
                
        return result

    def __initSamples(self, samplesToPlot):
        from copy import deepcopy
        self.plotter = PlotManager()
        if "data" in samplesToPlot:
            self.plotter.add_sample(samples.data, "Data (7 TeV)", **style.DATA_STYLE)    
        if "qcd" in samplesToPlot or "mc" in samplesToPlot:
            self.plotter.add_sample(samples.qcd_mc, "QCD MC", **style.QCD_MC_STYLE_HIST)
        if "minBias" in samplesToPlot:
            self.plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)
        if "Ztautau" in samplesToPlot or "mc" in samplesToPlot:
            ZTAUTAU_MC_STYLE_HIST = deepcopy(style.QCD_MC_STYLE_HIST)
            ZTAUTAU_MC_STYLE_HIST["line_color"] = ROOT.EColor.kRed - 2
            ZTAUTAU_MC_STYLE_HIST["fill_color"] = ROOT.EColor.kRed - 4
            self.plotter.add_sample(samples.ztautau_mc, "Z#rightarrow#tau#tau MC x 10^{4}", **ZTAUTAU_MC_STYLE_HIST)
            for sample in samples.ztautau_mc.subsamples:
                sample.get_events()
                sample.scaleFactor *= 1e4

        self.plotter.set_integrated_lumi(samples.data.effective_luminosity())
        self.manager = samples.data.build_ntuple_manager("tauIdEffNtuple")
        self.main = self.manager.get_ntuple("patPFTausDijetTagAndProbeShrinkingCone")
        self.hlt = self.manager.get_ntuple("TriggerResults")

    def __makeMenu(self, sessions):
        menuMap = {}
        i = 0
        for sessionName in sessions:
            print "%s: (%s) %s"%(i,sessionName ,sessions[sessionName]["description"])
            menuMap[i] = sessionName
            i+=1
        item = raw_input("choose session name (0-%s) "%(i-1))
        if not int(item) in menuMap:
            name = None
        else:
            name = menuMap[int(item)]
        return name

    def __makeLabel(self, text):
        result = ROOT.TPaveText(0.12, 0.68, 0.45, 0.73, "NDC")
        result.AddText(text)
        result.SetTextAlign(13)
        result.SetTextSize(0.04)
        result.SetFillStyle(0)
        result.SetBorderSize(0)
        return result

    def __possibleSessions(self, baseDir):
        #possible Session definitions
        base_selection = self.hlt.expr('$hltJet15U > 0.5')\
                         & self.main.expr("abs($jetEta) < 2.5 & $jetPt > 10 & $probe > 0.5")
        sessions = {
            "all":{"selection": base_selection,
                   "subdir":"%s/combined/"%baseDir,
                   "description":"All Combinded",
                   "label": self.__makeLabel("")
                   },
            "tanc":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent"),
                    "subdir":"%s/TaNCHalfPercent/combined/"%baseDir,
                    "description":"TaNC HalfPercent Combinded",
                    "label": self.__makeLabel("TaNC 0.5% tune"),
                    },
            "tanc0":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 0"),
                     "subdir":"%s/TaNCHalfPercent/oneProngNoPi0/"%baseDir,
                     "description":"TaNC HalfPercent OneProngNoPi0",
                     "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-} #nu_{#tau}"),
                     },
            "tanc1":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 1"),
                     "subdir":"%s/TaNCHalfPercent/oneProngOnePi0/"%baseDir,
                     "description":"TaNC HalfPercent One Prong One Pi0",
                     "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-} #pi^{0} #nu_{#tau}"),
                     },
            "tanc2":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 2"),
                     "subdir":"%s/TaNCHalfPercent/oneProngTwoPi0/"%baseDir,
                     "description":"TaNC HalfPercent One Prong Two Pi0",
                     "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-} #pi^{0} #pi^{0} #nu_{#tau}"),
                     },
            "tanc10":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 10"),
                      "subdir":"%s/TaNCHalfPercent/ThreeProngNoPi0/"%baseDir,
                      "description":"TaNC HalfPercent OneProngNoPi0",
                      "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-} #pi^{+} #pi^{-} #nu_{#tau}"),
                      },
            "tanc11":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 11"),
                      "subdir":"%s/TaNCHalfPercent/ThreeProngOnePi0/"%baseDir,
                      "description":"TaNC HalfPercent One Prong One Pi0",
                      "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-} #pi^{+} #pi^{-} #pi^{0} #nu_{#tau}"),
                      },
            "pre":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut"),
                    "subdir":"%s/Prediscriminant/combined/"%baseDir,
                    "description":"After Predescriminant Combinded",
                    "label": self.__makeLabel("Prediscriminant"),
                    },
            "pre0":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 0"),
                     "subdir":"%s/Prediscriminant/oneProngNoPi0/"%baseDir,
                     "description":"After Predescriminant OneProngNoPi0",
                     "label": self.__makeLabel("Prediscriminant, #pi^{-} #nu_{#tau}"),
                     },
            "pre1":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 1"),
                     "subdir":"%s/Prediscriminant/oneProngOnePi0/"%baseDir,
                     "description":"After Predescriminant One Prong One Pi0",
                     "label": self.__makeLabel("Prediscriminant, #pi^{-} #pi^{0} #nu_{#tau}"),
                     },
            "pre2":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 2"),
                     "subdir":"%s/Prediscriminant/oneProngTwoPi0/"%baseDir,
                     "description":"After Predescriminant One Prong Two Pi0",
                     "label": self.__makeLabel("Prediscriminant, #pi^{-} #pi^{0} #pi^{0} #nu_{#tau}"),
                     },
            "pre10":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 10"),
                      "subdir":"%s/Prediscriminant/ThreeProngNoPi0/"%baseDir,
                      "description":"After Predescriminant OneProngNoPi0",
                      "label": self.__makeLabel("Prediscriminant, #pi^{-} #pi^{+} #pi^{-} #nu_{#tau}"),
                      },
            "pre11":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 11"),
                      "subdir":"%s/Prediscriminant/ThreeProngOnePi0/"%baseDir,
                      "description":"After Predescriminant One Prong One Pi0",
                      "label": self.__makeLabel("Prediscriminant, #pi^{-} #pi^{+} #pi^{-} #pi^{0} #nu_{#tau}"),
                      },
            "preNot0":{"selection": base_selection & self.main.expr("$byLeadTrackPtCut & $TaNCDecayMode == 11"),
                      "subdir":"%s/Prediscriminant/notOneProngNoPi0/"%baseDir,
                      "description":"After Predescriminant Not mode 0",
                      "label": self.__makeLabel("Prediscriminant, not #pi^{-} #nu_{#tau}"),
                      },
            
            "tancNot0":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode != 0"),
                        "subdir":"%s/TaNCHalfPercent/NotOneProngNoPi0/"%baseDir,
                        "description":"TaNC HalfPercent Not OneProngNoPi0",
                        "label": self.__makeLabel("TaNC 0.5% tune, not #pi^{-}#nu_{#tau}"),
                        },
            "tan0noLL":{"selection": base_selection & self.main.expr("$byTaNCfrHalfPercent & $TaNCDecayMode == 0 & $againstElectron > 0.5& $againstMuon > 0.5"),
                        "subdir":"%s/TaNCHalfPercent/oneProngNoPi0NoLightLeptons/"%baseDir,
                        "description":"TaNC HalfPercent OneProngNoPi0 no light leptons",
                        "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-}#nu_{#tau}, no e, #mu"),
                        },
            }
        return sessions

def main():
    ROOT.gROOT.SetBatch(True)
    #    ["tanc","tanc0","tanc1","tanc2","tanc10","tanc11", "all"]
    #["pre0", "pre1", "pre2", "pre10",  "pre11", "pre" ]
    sessionNames = sys.argv[1:]
    if sessionNames == []: sessionNames=[None]
    for name in sessionNames:
        print "+++++++++++++++++++++++++++++++++++++"
        print "++++++++ starting %s +++++++++++++"%name
        print "+++++++++++++++++++++++++++++++++++++"
        session = PlotSession(name = name, samplesToPlot = ["data","qcd"], baseDir="/afs/cern.ch/user/e/edelhoff/scratch0/pFlow/analysis/CMSSW_3_6_1/src/TauAnalysis/TauIdEfficiency/test/commissioning/plots.v2_1/")
        session.drawDecayMode()
        #session.drawKineticPlots() #...done
        #session.drawMainTrackPlots() #... done
        #session.drawMassPlots() #pre11 pre
        #session.drawOutlierSumKineticPlots() # ... done
        #session.drawOutlierCountPlots() # ... done
        #session.drawOutlierMassPlots() # ... done
        #session.drawOutlierSingleKineticPlots(charges=["Charged"]) 
        #session.drawPi0KineticPlots() # ... done
        #session.drawTrackKineticPlots() #... done
        #session.drawTaNCOutputPlots()
        
        #session.drawEssential() 
        #session.drawComplementaryPlots()
        #session.drawOutlierPlots()
        
        
if __name__ == "__main__":
    main()
    

    
