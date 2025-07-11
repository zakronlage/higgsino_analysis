import numpy as np
# import matplotlib.pyplot as plt
import ROOT
import ctypes

def cosweight(time):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    return np.cos(omega * time)

def sinweight(time):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    return np.sin(omega * time) 

if __name__ == "__main__":
    c = ROOT.TCanvas()

    #vertically combined dataset
    chain = ROOT.TChain("qredtree"); #all called qredtree, with cycle # after
    chain.Add("SuperReduced_Background_ds3808.root"); #Data from Jul-Aug 2019
    chain.Add("SuperReduced_Background_ds3814.root"); #Jul-Aug 2020
    chain.Add("SuperReduced_Background_ds3824.root"); #Jun-Aug 2022
    df = ROOT.RDataFrame(chain); #special rdf constructor for TChain

    collectFiltered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1")

    withCosWeights = collectFiltered.Define("CosWeighting", cosweight, ["EventTime"]) #defines new column CosWeighting using EventTime, called CosWeight
    cosHisto = withCosWeights.Histo1D(["e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.], "Energy", "CosWeighting")

    cosHisto.SetDefaultSumw2()
    cosHisto.SetLineColor(ROOT.kBlue)
    cosHisto.Draw()

    withSinWeights = collectFiltered.Define("SinWeighting", sinweight, ["EventTime"])
    sinHisto = withSinWeights.Histo1D(["e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.], "Energy", "SinWeighting")
    sinHisto.SetLineColor(ROOT.kMagenta)
    sinHisto.Draw("Same")

    legend = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend.SetHeader("Legend")
    legend.AddEntry(cosHisto, "Energy Weighted (Cos)", "l")
    legend.AddEntry(sinHisto, "Energy Weighted (Sin)", "l")
    legend.Draw()

    c.SaveAs("dsCombined.pdf")


    #k40 peak combined cos only
    c1 = ROOT.TCanvas()
    histo = withCosWeights.Histo1D(["e", "K-40 Peak Cosine Weighting (Combined)", 75, 1400., 1550.], "Energy", "CosWeighting")
    histo.SetLineColor(ROOT.kBlue)
    histo.Draw()

    legend2 = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend2.SetHeader("Legend")
    legend2.AddEntry(histo, "Weighted Energy", "l")
    legend2.Draw()

    c1.SaveAs("k40combinedCosOnly.pdf")


    #k40 peak combined sin only
    c2 = ROOT.TCanvas()
    histo2 = withSinWeights.Histo1D(["e", "K-40 Peak Sine Weighted (Combined)", 75, 1400., 1550.], "Energy", "SinWeighting")
    histo2.SetLineColor(ROOT.kMagenta)
    histo2.DrawClone()

    legend3 = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend3.SetHeader("Legend")
    legend3.AddEntry(histo2, "Weighted Energy", "l")
    legend3.Draw()

    c2.SaveAs("k40combinedSinOnly.pdf")


    #3801 dataset only
    ex = ROOT.RDataFrame("qredtree;5", "SuperReduced_Background_ds3801.root")
    c3 = ROOT.TCanvas()

    appliedCut4 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1")
    withCosWeights2 = appliedCut4.Define("CosWeighting2", cosweight, ["EventTime"])
    histo3 = withCosWeights2.Histo1D(["e", "Singular Weighted Energy Spectrum", 5000, 0., 10000.], "Energy", "CosWeighting2")
    histo3.SetLineColor(ROOT.kBlue)
    histo3.DrawClone()

    withSineWeights2 = appliedCut4.Define("SinWeighting2", sinweight, ["EventTime"])
    histo4 = withSineWeights2.Histo1D(["e", "Singular Weighted Energy Spectrum", 5000, 0., 10000.], "Energy", "SinWeighting2")
    histo4.SetLineColor(ROOT.kMagenta)
    histo4.DrawClone("Same")

    legend4 = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend4.SetHeader("Legend")
    legend4.AddEntry(histo3, "Energy Weighted (Cosine)", "l")
    legend4.AddEntry(histo4, "Energy Weighted (Sine)", "l")
    legend4.Draw()

    c3.SaveAs("cosSinWeighted.pdf")

    #3801 k-40 peak cosine only
    c4 = ROOT.TCanvas()
    cos_k40Peak = withCosWeights2.Histo1D(["eB", "K-40 Peak Cosine Weighted (Singular)", 75, 1400., 1550.], "Energy", "CosWeighting2")
    cos_k40Peak.SetLineColor(ROOT.kBlue)
    cos_k40Peak.DrawClone()

    legend5 = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend5.SetHeader("Legend")
    legend5.AddEntry(cos_k40Peak, "Weighted Energy", "l")
    legend5.Draw()

    c4.SaveAs("k40cosonly.pdf")

    #3801 k-40 peak sine only
    c5 = ROOT.TCanvas()
    sin_k40Peak = withSineWeights2.Histo1D(["eB", "K-40 Peak Sine Weighted (Singular)", 75, 1400., 1550.], "Energy", "SinWeighting2")
    sin_k40Peak.SetLineColor(ROOT.kMagenta)
    sin_k40Peak.DrawClone()

    legend6 = ROOT.TLegend(0.65, 0.75, 1.0, 0.94)
    legend6.SetHeader("Legend")
    legend6.AddEntry(sin_k40Peak, "Weighted Energy", "l")
    legend6.Draw()

    c5.SaveAs("k40sinonly.pdf")

    c6 = ROOT.TCanvas()
    ratio = sinHisto.Clone("ratio")
    ratio.Divide(cosHisto.GetPtr())
    ratio.SetTitle("Sine/Cosine Ratio")
    ratio.SetLineColor(ROOT.kRed)
    ratio.Draw("pE")

    sin_events = sinHisto.Integral()
    cos_events = cosHisto.Integral()
    ratio_value = sin_events / cos_events
    print("ratio = %f / %f = %.2f\n", sin_events, cos_events, ratio_value)
    print("phase = %.2f degrees\n", ROOT.TMath.ATan(ratio_value) * 180.0 / ROOT.TMath.Pi())

    c7 = ROOT.TCanvas()
    ratio2 = sin_k40Peak.Clone("ratio2")
    ratio2.Divide(cos_k40Peak.GetPtr())
    ratio2.SetTitle("Sine/Cosine Ratio (Singular)")
    ratio2.SetLineColor(ROOT.kRed)
    ratio2.Draw("E")