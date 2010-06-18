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
        self.subdir = sessions[name]["subdir"]
        self.description = sessions[name]["description"]
        self.labels=[style.CMS_PRELIMINARY_UPPER_LEFT, style.LUMI_LABEL_UPPER_LEFT,
                     style.PT_ETA_CUT_LABEL_UPPER_LEFT, sessions[name]["label"]]
        
        if not os.path.isdir(self.subdir):
            os.makedirs(self.subdir)
            
        self.plots = {}
        self.canvas = ROOT.TCanvas("blah", "blah", 500, 500)

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
        self.__printPlots(plotName)

    def drawKineticPlots(self):
        plotName = "TauPt"
        self.plots[plotName] = self.plotter.distribution(
            expression=self.main.expr('$TaNCPt'),
            selection= self.selection & self.__getDecayModeSelection(plotName),
            binning = (100, 0, 75),
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
            binning = (100, 0, 2.6),
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
                x_axis_title = "Sum %s Outlier P_{T} [GeV/c]"%(charge),
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
                binning = (101, -0.5, 100.5),
                x_axis_title = "%s Outlier Count"%(charge),
                #y_min = 1e-2,
                logy=False,
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
                
    def drawOutlierSingleKineticPlots(self):
        for charge in ["","Charged","Neutral"]:
            for i in range(1,5):
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
                binning = (100, 0, 20),
                x_axis_title = "p_t of %s %s No. %s [GeV/c]"%(name,objects,number),
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
                #logy=True,
                labels = self.labels
                )            
            self.__printPlots(plotName)
        except StandardError, err:
            print "TaNC%s%sAngle%s not found: %s"%(name,objects,number,err)
        
    def __printPlots(self, name):
        self.plots[name]['legend'].make_legend().Draw() 
        for format in self.formats:
            self.canvas.SaveAs("%s/%s.%s"%(self.subdir, name, format))

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

    def __initSamples(self, samplesToPlot):
        self.plotter = PlotManager()
        if "data" in samplesToPlot:
            self.plotter.add_sample(samples.data, "Data (7 TeV)", **style.DATA_STYLE)    
        if "qcd" in samplesToPlot or "mc" in samplesToPlot:
            self.plotter.add_sample(samples.qcd_mc, "QCD MC", **style.QCD_MC_STYLE_HIST)
        if "minBias" in samplesToPlot:
            self.plotter.add_sample(samples.minbias_mc, "Minbias MC", **style.MINBIAS_MC_STYLE)
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
                 "label": self.__makeLabel("TaNC 0.5% tune, #pi^{-}#nu_{#tau}"),
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
    mySession = PlotSession(name = None, baseDir="testPlots")

    mySession.drawDecayMode()
    mySession.drawKineticPlots()
    mySession.drawMainTrackPlots()
    mySession.drawMassPlots()
    mySession.drawOutlierSumKineticPlots()
    mySession.drawOutlierCountPlots()
    mySession.drawOutlierMassPlots()
    mySession.drawOutlierSingleKineticPlots()
    mySession.drawPi0KineticPlots()
    mySession.drawTrackKineticPlots()
    
    
if __name__ == "__main__":
    main()
    

    
