import numpy as np
import matplotlib.pyplot as plt
import ROOT
import ctypes

def sin2weight(time):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    phase = 0
    return (1 - np.cos(2 * omega * (time- phase))) / 2 #equivalent to sin squared weighting

def cos2weight(time):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    phase = 0
    return np.cos(2 * omega * (time - phase))**2 #equivalent to cos squared weighting


#this function is used for integration to demodulate and extract phase & amplitude
def cosDemod(time, params): #requires params argument even if unused since integration
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    time1 = time[0]
    phase = 0
    return np.cos(2 * omega * (time1 - phase))

def sinDemod(time, params):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day
    time1 = time[0]
    phase = 0
    return np.sin(2 * omega * (time1 - phase))

if __name__ == "__main__":
    # Demodulation()
    chain = ROOT.TChain("qredtree")
    chain.Add("SuperReduced_Background_ds3808.root") #Jul-Aug 2019
    chain.Add("SuperReduced_Background_ds3814.root") #Jul-Aug 2020
    chain.Add("SuperReduced_Background_ds3824.root") #Jun-Aug 2022
    df = ROOT.RDataFrame(chain)

    # filtered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1")
    # withSin2Weights = filtered.Define("Sin2Weighting", sin2weight, {"EventTime"}) #defines new column Sin2Weighting using EventTime
    # withCos2Weights = filtered.Define("CosWeighting", cos2weight, {"EventTime"}) #cosine weighting

    ex = ROOT.RDataFrame("qredtree", "SuperReduced_Background_ds3814.root")
    filtered2 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1")
    withSin2Weights2 = filtered2.Define("Sin2Weighting", sin2weight, ["EventTime"])
    withCos2Weights2 = filtered2.Define("CosWeighting", cos2weight, ["EventTime"])

    cols = withSin2Weights2.AsNumpy(["EventTime","Energy", "Sin2Weighting"])  # Extracting the necessary columns
    timeArr = cols["EventTime"]
    energyArr = cols["Energy"]
    weightArr = cols["Sin2Weighting"]

    print("Number of events:", len(timeArr))

    # Plotting the energy distribution with sine squared weighting
    # energyScatter = plt.scatter(timeArr, energyArr, c=weightArr, cmap='viridis', s=1)
    # plt.colorbar(energyScatter, label='Sine Squared Weighting')
    # plt.xlabel('Event Time (s)')
    # plt.ylabel('Energy (keV)')
    # plt.ylim(1280, 1530)

    E_upper = 1530
    E_lower = 1280
    T_lower = np.min(timeArr)
    T_upper = np.max(timeArr)
    eRes = 1 # 1 keV resolution
    tRes = 1000*3600*3 # hr resolution in milliseconds


    eBins = int((E_upper - E_lower)/eRes)
    timeBins = int((T_upper - T_lower)/tRes)
    global EvT_Hist
    EvT_Hist = ROOT.TH2D("EvT_Hist", "Energy vs Time;Event Time/3h (ms);Energy/KeV (KeV)", timeBins, T_lower, T_upper, eBins, E_lower, E_upper)
    for i in range(len(timeArr)):
        EvT_Hist.Fill(timeArr[i], energyArr[i], weightArr[i])
    # c1 = ROOT.TCanvas("c1", "Energy vs Time", 800, 600)
    # c1.Draw()
    EvT_Hist.Draw("Lego")
    # c1.SaveAs("Energy_vs_Time.png")

    etScatter = ROOT.TGraph2D(len(timeArr), np.float64(timeArr), np.float64(energyArr), np.float64(weightArr))

    etScatter.GetYaxis().SetRangeUser(1280, 1530)
    etScatter.GetXaxis().SetTitleOffset(2.5)
    etScatter.GetYaxis().SetTitleOffset(2)
    etScatter.SetTitle("Energy Weighted vs Time (Summer 2020);Unix Time (ms);Energy (keV);Sin Squared Weighting")
    etScatter.SetMarkerStyle(8) #changes from default square to large dot
    etScatter.SetMarkerSize(0.4)
    etScatter.Draw("PCOL0") # plot markers with color, no stats box

    etScatter.GetYaxis().SetRangeUser(1459, 1461) #redefine to narrow gap instead
    total = timeArr.size
    validEnergies = np.array(etScatter.GetY())
    validEnergies = validEnergies[(1459. <= validEnergies) & (validEnergies <= 1461.)] # filter energies in range [1459, 1461]
    nEvents = validEnergies.size
    
    inPhase = ROOT.TF1("f1", cosDemod, timeArr[0], timeArr[-1]); #create a TF1 object for cosine demodulation
    quadrature = ROOT.TF1("f2", sinDemod, timeArr[0], timeArr[-1]); #create a TF1 object for sine demodulation

    err = ctypes.c_double(0.0)
    integ = inPhase.IntegralOneDim(timeArr[0], timeArr[-1],1e-3,1e-3,err) #create a TF1 object for cosine demodulation, 1.E-3, 1.E-3, *err); //1.E-3 is relaxed tolerance
    integ2 = quadrature.IntegralOneDim(timeArr[0], timeArr[-1],1e-3,1e-3,err) # create a TF1 object for cosine demodulation, 1.E-3, 1.E-3, *err);
    iqPhase = ROOT.TMath.ATan2(nEvents * integ2, nEvents * integ)
    #arctan will cancel out constants leaving only phase
    amplitude = ROOT.TMath.Hypot(nEvents * integ, nEvents * integ2)

    print("\n\nThe total number of events in Jul-Aug 2020 dataset is: ", total)
    print("The number of events in the energy range [1459, 1461] is: ", nEvents)
    print("\nThe phase is: ", iqPhase)
    print("The amplitude is: ", amplitude, "\n")