import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetStyle("Plain")
import math
import sys

steering = {
    'mc' : {
        'files' : [
            ROOT.TFile("/data1/friis/ZmumuEff/efficiency_mcpu.root", "READ"),
            ROOT.TFile("/data1/friis/ZmumuEff/efficiencySta_mcpu.root", "READ"),
            ROOT.TFile("/data1/friis/ZmumuEff/efficiencyInner_mcpu.root", "READ"),
        ],
        'plots' : {},
    },
    'data' : {
        'files' : [
            ROOT.TFile("/data1/friis/ZmumuEff/efficiency_data.root", "READ"),
            ROOT.TFile("/data1/friis/ZmumuEff/efficiencySta_data.root", "READ"),
            ROOT.TFile("/data1/friis/ZmumuEff/efficiencyInner_data.root", "READ"),
        ],
        'plots' : {},
    },
}

def get_by_type(folder, type):
    for key in folder.GetListOfKeys():
        object = key.ReadObj()
        if object.InheritsFrom(type):
            yield object

def maximum(graph):
    y_bins = graph.GetY()
    y_err_bins = graph.GetEYhigh()
    return max(y_bins[i]+y_err_bins[i] for i in range(graph.GetN()))

def minimum(graph):
    y_bins = graph.GetY()
    y_err_bins = graph.GetEYlow()
    return min(y_bins[i]-y_err_bins[i] for i in range(graph.GetN()))

canvas = ROOT.TCanvas("blah", "blah", 800, 600)

for sample, info in steering.iteritems():
    # Loop over each file
    for file in info['files']:
        # Get the first ttree in the file (the only one)
        folder = list(get_by_type(file, "TDirectory"))[0]
        # Now get all the plot subfolders
        for subfolder in get_by_type(folder, "TDirectory"):
            name = subfolder.GetName()
            # Remove the "tag" variables
            name_key = tuple(datum for datum in name.split("_") if datum != "tag")
            #print name_key
            fit_plot_folder = subfolder.Get("fit_eff_plots")
            if not fit_plot_folder:
                print "ERROR: can't get", name
                continue
            fit_canvas = fit_plot_folder.GetListOfKeys()[0].ReadObj()
            plot = fit_canvas.GetPrimitive("hxy_fit_eff")
            info['plots'][name_key] = plot

print "Making plots"

x_axes = {
    'abseta' : "|#eta|",
    'pt' : "p_{T} [GeV/c]",
    'nVertices' : "N_{vtx}",
}

titles = {
    "id": "VBTF ID Efficiency w.r.t. Global Muons",
    "iso": "Abs. Iso. Efficiency w.r.t. Global + VBTF",
    "reliso": "Rel. Iso. Efficiency w.r.t. Global + VBTF",
    "looseiso":  "Loose Abs. Iso. Efficiency w.r.t. Global + VBTF",
    "loosereliso": "Loose Rel. Iso. Efficiency w.r.t. Global + VBTF",
    "trigA" : "HLT Mu9 w.r.t offline (Runs 132440-147116)",
    "trigB" : "HLT IsoMu9 | Mu11 w.r.t offline (Runs 147196-148058)",
    "trigC" : "HLT IsoMu13 | Mu15 w.r.t offline (Runs 148059-149442)",
    "standAlone" : "STA muon given inner track",
    "linking" : "Global muon given inner and outer track",
    "innerTrack" : "Inner track given STA muon",
    "hltMu9" : "HLT Mu9 w.r.t. offline",
    "trigComp" : "HLT Trigger efficiency (data) compared to MC HLT_Mu9",
}

# Make the special trigger plots
# Find the MC hltMu9 efficiencies
mc_hlt_mu9_plots = [key for key in steering['mc']['plots'].keys()
                    if key[1] == 'hltMu9' and key[2] != 'nVertices']


#print [key[1] for key in steering['data']['plots'].keys()]

lumi_info = {
    'trigA' : 8237412.121,
    'trigB' : 9470933.970,
    'trigC' : 18436654.902,
}

def merge_efficiencies(*inputs):
    # Keep track of the weighted sums in each bin
    nbins = inputs[0][1].GetN()
    y_sum = [0]*nbins
    y_err_high_sum2 = [0]*nbins
    y_err_low_sum2 = [0]*nbins
    for weight, input in inputs:
        for bin in range(nbins):
            y_sum[bin] += weight*input.GetY()[bin]
            y_err_high_sum2[bin] += weight*(input.GetEYhigh()[bin]**2)
            y_err_low_sum2[bin] += weight*(input.GetEYlow()[bin]**2)
    norm = sum(weight for weight, input in inputs)
    normalize = lambda x: x/norm
    y_sum = map(normalize, y_sum)
    y_err_low_sum2 = map(math.sqrt, map(normalize, y_err_low_sum2))
    y_err_high_sum2 = map(math.sqrt, map(normalize, y_err_high_sum2))

    # Just copy our input so we don't have worry about the x_errors
    #output = ROOT.TGraphAsymmErrors(inputs[0][1])
    output = ROOT.TGraphAsymmErrors(nbins)

    for bin, (y, e_high, e_low) in enumerate(
        zip(y_sum, y_err_high_sum2, y_err_low_sum2)):
        output.SetPoint(bin, inputs[0][1].GetX()[bin], y)
        output.SetPointEXlow(bin, inputs[0][1].GetEXlow()[bin])
        output.SetPointEXhigh(bin, inputs[0][1].GetEXhigh()[bin])
        output.SetPointEYlow(bin, e_low)
        output.SetPointEYhigh(bin, e_high)
    return output

for plot in mc_hlt_mu9_plots:
    print plot
    # We are making our special trigger comparison.  For MC it is the same, so
    # just copy it over.
    steering['mc']['plots'][(plot[0], 'trigComp', plot[2])] = \
            steering['mc']['plots'][plot]

    # We want to compare against each one individually as well
    for period in lumi_info.keys():
        steering['mc']['plots'][(plot[0], period, plot[2])] = \
                steering['mc']['plots'][plot]

    # Find the corresponding data plots
    data_plots = steering['data']['plots']
    trigger_plots = [
        (lumi_info[trig], data_plots[(plot[0], trig, plot[2])])
        for trig in lumi_info.keys()
    ]
    new_plot = merge_efficiencies(*trigger_plots)
    steering['data']['plots'][(plot[0], 'trigComp', plot[2])] = new_plot

#sys.exit()

for plot in steering['mc']['plots'].keys():
    #print plot
    mc_plot = steering['mc']['plots'][plot]
    data_plot = steering['data']['plots'][plot]

    mc_plot.SetMarkerStyle(24)
    mc_plot.SetMarkerColor(ROOT.EColor.kRed)
    mc_plot.SetMarkerSize(1)
    mc_plot.Draw("ape")
    mc_plot.SetFillColor(ROOT.EColor.kAzure +7)

    data_plot.SetMarkerStyle(20)
    data_plot.SetMarkerColor(ROOT.EColor.kBlack)
    data_plot.SetMarkerSize(1)
    data_plot.Draw("pe")

    mc_max = maximum(mc_plot)
    data_max = maximum(data_plot)
    the_max = mc_max > data_max and mc_max or data_max

    mc_min = minimum(mc_plot)
    data_min = minimum(data_plot)
    the_min = mc_min < data_min and mc_min or data_min

    #print mc_max, data_max, the_max
    #print mc_min, data_min, the_min

    mc_plot.GetHistogram().SetMaximum(the_max + (the_max-the_min)*0.25)
    mc_plot.GetHistogram().SetMinimum(the_min - (the_max-the_min)*0.25)

    legend = ROOT.TLegend(0.6, 0.15, 0.85, 0.3, "", "brNDC")
    legend.AddEntry(data_plot, "Data", "p")
    legend.AddEntry(mc_plot, "Z#mu#mu MC (156BX)", "p")
    legend.SetFillStyle(0)
    legend.Draw()

    mc_plot.GetHistogram().SetTitle(titles[plot[1]])
    mc_plot.GetHistogram().GetXaxis().SetTitle(x_axes[plot[2]])

    file_name = "_".join(plot)

    canvas.Update()
    canvas.SaveAs(file_name+".png")
    canvas.SaveAs(file_name+".pdf")

    mc_scaled = ROOT.TGraphAsymmErrors(mc_plot.GetN())
    mc_scaled.SetMarkerStyle(24)
    mc_scaled.SetMarkerColor(ROOT.EColor.kRed)
    mc_scaled.SetMarkerSize(1)
    mc_scaled.SetFillColor(ROOT.EColor.kAzure +7)
    mc_scaled.SetFillStyle(2001)


    data_scaled = ROOT.TGraphAsymmErrors(mc_plot.GetN())
    data_scaled.SetMarkerStyle(20)
    data_scaled.SetMarkerColor(ROOT.EColor.kBlack)
    data_scaled.SetMarkerSize(1)

    # Make scaled plot
    for bin in range(mc_plot.GetN()):
        mc_y = mc_plot.GetY()[bin]
        mc_x = mc_plot.GetX()[bin]
        data_y = data_plot.GetY()[bin]
        data_x = data_plot.GetX()[bin]
        mc_err_y_high = mc_plot.GetEYhigh()[bin]
        mc_err_y_low = mc_plot.GetEYlow()[bin]
        data_err_y_high = data_plot.GetEYhigh()[bin]
        data_err_y_low = data_plot.GetEYlow()[bin]
        mc_err_x_high = mc_plot.GetEXhigh()[bin]
        mc_err_x_low = mc_plot.GetEXlow()[bin]
        data_err_x_high = data_plot.GetEXhigh()[bin]
        data_err_x_low = data_plot.GetEXlow()[bin]
        if mc_y > 0:
            data_y /= mc_y
            #print data_y
            mc_err_y_high /= mc_y
            mc_err_y_low /= mc_y
            data_err_y_high /= mc_y
            data_err_y_low /= mc_y
        mc_scaled.SetPoint(bin, mc_x, 1.)
        mc_scaled.SetPointEYhigh(bin,mc_err_y_high)
        mc_scaled.SetPointEYlow(bin,mc_err_y_low)

        mc_scaled.SetPointEXhigh(bin,mc_err_x_high)
        mc_scaled.SetPointEXlow(bin,mc_err_x_low)

        data_scaled.SetPoint(bin, data_x, data_y)
        data_scaled.SetPointEYhigh(bin,data_err_y_high)
        data_scaled.SetPointEYlow(bin,data_err_y_low)
        data_scaled.SetPointEXhigh(bin,data_err_x_high)
        data_scaled.SetPointEXlow(bin,data_err_x_low)

    mc_scaled.Draw("a2")
    data_scaled.Draw("p")

    mc_max = maximum(mc_scaled)
    data_max = maximum(data_scaled)
    the_max = mc_max > data_max and mc_max or data_max

    mc_min = minimum(mc_scaled)
    data_min = minimum(data_scaled)
    the_min = mc_min < data_min and mc_min or data_min

    mc_scaled.GetHistogram().SetMaximum(the_max + (the_max-the_min)*0.25)
    mc_scaled.GetHistogram().SetMinimum(the_min - (the_max-the_min)*0.25)

    mc_scaled.GetHistogram().SetTitle(titles[plot[1]])
    mc_scaled.GetHistogram().GetXaxis().SetTitle(x_axes[plot[2]])

    canvas.Update()
    canvas.SaveAs(file_name + "_scaled" + ".png")
    canvas.SaveAs(file_name + "_scaled" + ".pdf")

for plot in [key for key in steering['data']['plots'].keys() if
             key not in steering['mc']['plots']]:
    data_plot = steering['data']['plots'][plot]
    data_plot.SetMarkerStyle(20)
    data_plot.SetMarkerColor(ROOT.EColor.kBlack)
    data_plot.SetMarkerSize(1)
    data_plot.Draw("ape")

    data_plot.GetHistogram().SetTitle(titles[plot[1]])
    data_plot.GetHistogram().GetXaxis().SetTitle(x_axes[plot[2]])

    file_name = "_".join(plot)

    canvas.Update()
    canvas.SaveAs(file_name+".png")
    canvas.SaveAs(file_name+".pdf")

