import numpy as np
import matplotlib.pyplot as plt
import ROOT
import ctypes

def sin2weight(time):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day in Hz
    phase = 0
    return (1 - np.cos(2 * omega * (time- phase))) / 2 #equivalent to sin squared weighting

def cos2weight(time):
    omega = (2 * np.pi) / ((3600*23 + 56*60 + 4)) #angular frequency of sidereal day in Hz
    phase = 0
    return np.cos(2 * omega * (time - phase))**2 #equivalent to cos squared weighting


#this function is used for integration to demodulate and extract phase & amplitude
def cosDemod(time, params): #requires params argument even if unused since integration
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day in Hz
    time1 = time[0]
    phase = 0
    return np.cos(2 * omega * (time1 - phase))

def sinDemod(time, params):
    omega = (2 * np.pi) / (3600*23 + 56*60 + 4) #angular frequency of sidereal day in Hz
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
    timeArr = cols["EventTime"]/1000 # Convert time from milliseconds to seconds
    energyArr = cols["Energy"]
    weightArr = cols["Sin2Weighting"]

    print("Number of events:", len(timeArr))

    E_upper = 1530
    E_lower = 1280
    T_lower = timeArr[0] # start time

    timeArr = timeArr - T_lower # normalize time to start at 0
    T_upper = np.max(timeArr) # now T_upper is relative to T_lower
    T_lower = 0 # now T_lower is 0
    sidereal_day = (3600*23 + 56*60 + 4) #seconds in a sidereal day
    

    eRes = 1 # 1 keV resolution
    tRes = 3600*3 # 3hr resolution in seconds

    eBins = int((E_upper - E_lower)/eRes)
    timeBins = int((T_upper - T_lower)/tRes)
    EvT_Hist = ROOT.TH2D("EvT_Hist", "Energy vs Time;Elapsed Time/3h (s);Energy/KeV (KeV)", timeBins, T_lower, T_upper, eBins, E_lower, E_upper)
    for i in range(len(timeArr)):
        EvT_Hist.Fill(timeArr[i], energyArr[i], weightArr[i])
    c1 = ROOT.TCanvas("c1", "Energy vs Time", 800, 600)
    c1.Draw()
    EvT_Hist.Draw("Lego")
    EvT_Hist.GetYaxis().SetRangeUser(1450,1470) #set y-axis range

    c2 = ROOT.TCanvas("TC", "Energy vs Time", 800, 600)
    etScatter = ROOT.TGraph2D(len(timeArr), np.float64(timeArr), np.float64(energyArr), np.float64(weightArr))
    etScatter.GetYaxis().SetRangeUser(1280, 1530)
    etScatter.GetXaxis().SetTitleOffset(2.5)
    etScatter.GetYaxis().SetTitleOffset(2)
    etScatter.SetTitle("Energy Weighted vs Time (Summer 2020);Elapsed Time (s);Energy (keV);Sin Squared Weighting")
    etScatter.SetMarkerStyle(8) #changes from default square to large dot
    etScatter.SetMarkerSize(0.4)
    c2.Draw()
    etScatter.Draw("PCOL0") # plot markers with color, no stats box

    periods = timeArr % sidereal_day # get the periods of the sidereal day

    etScatter.GetYaxis().SetRangeUser(1459, 1461) #redefine to narrow gap instead
    total = timeArr.size
    validEnergies = np.array(etScatter.GetY())
    validEnergies = validEnergies[periods.argmin():periods.argmax()] #filter events in the range of periods
    validEnergies = validEnergies[(1459. <= validEnergies) & (validEnergies <= 1461.)] # filter energies in range [1459, 1461]
    nEvents = validEnergies.size
    
    inPhase = ROOT.TF1("f1", cosDemod, timeArr[0], timeArr[-1]) #create a TF1 object for cosine demodulation
    quadrature = ROOT.TF1("f2", sinDemod, timeArr[0], timeArr[-1]) #create a TF1 object for sine demodulation

    err = ctypes.c_double(5)

    ##should come up with more robust way to find a period, but this works for now
    inPhaseIntegral = inPhase.IntegralOneDim(timeArr[periods.argmin()], timeArr[periods.argmax()],1e-6,1e-6,err)
    quadIntegral = quadrature.IntegralOneDim(timeArr[periods.argmax()], timeArr[periods.argmax()],1e-6,1e-6,err)
    iqPhase = ROOT.TMath.ATan2(nEvents * quadIntegral, nEvents * inPhaseIntegral)
    #arctan will cancel out constants leaving only phase
    amplitude = ROOT.TMath.Hypot(nEvents * inPhaseIntegral, nEvents * quadIntegral)

    print("\n\nThe total number of events in Jul-Aug 2020 dataset is: ", total)
    print("The number of events in the energy range [1459, 1461] is: ", nEvents)
    print("\nThe phase is: ", iqPhase)
    print("The amplitude is: ", amplitude, "\n")