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
    # withSin2Weights = filtered.Define("Sin2Weighting", sin2weight, ["EventTime"]) #defines new column Sin2Weighting using EventTime
    # withCos2Weights = filtered.Define("CosWeighting", cos2weight, ["EventTime"]) #cosine weighting

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
    T_begin = timeArr[0] # start time

    timeArr = timeArr - T_begin # normalize time to start at 0
    T_end = np.max(timeArr) # now T_upper is relative to T_lower
    T_begin = 0 # now T_lower is 0
    sidereal_day = (3600*23 + 56*60 + 4) #seconds in a sidereal day

    energyArr = energyArr[timeArr > 0] # filter out negative times, idk why tf this is happening when processing multiple datasets
    weightArr = weightArr[timeArr > 0] 
    timeArr = timeArr[timeArr > 0]
    

    eRes = 1 # 1 keV resolution
    tRes = 3600*3 # 3hr resolution in seconds

    eBins = int((E_upper - E_lower)/eRes)
    timeBins = int((T_end - T_begin)/tRes)
    EvT_Hist = ROOT.TH2D("EvT_Hist", "Energy vs Time;Elapsed Time/3h (s);Energy/KeV (KeV)", timeBins, -100, T_end, eBins, E_lower, E_upper)
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

    nperiods = timeArr % sidereal_day # get the periods of the sidereal day
    begin_period_index = np.asarray(np.isclose(nperiods,0,atol=10)).nonzero()[0].astype(int) # find the indices where the period starts
    end_period_index = np.asarray(np.isclose(nperiods,sidereal_day,atol=10)).nonzero()[0].astype(int) # find the indices where the period ends

    etScatter.GetYaxis().SetRangeUser(1459, 1461) #redefine to narrow gap instead
    total = timeArr.size
    validEnergies = np.array(etScatter.GetY())
    
    inPhase = ROOT.TF1("f1", cosDemod, timeArr[0], timeArr[-1]) #create a TF1 object for cosine demodulation
    quadrature = ROOT.TF1("f2", sinDemod, timeArr[0], timeArr[-1]) #create a TF1 object for sine demodulation

    err = ctypes.c_double(5)

    ##should come up with more robust way to find a period, but this works for now
    inPhaseIntegral = []
    quadIntegral = []
    nEvents = []
    durationArr = []
    #BAD HACK: checked, and end_period_index is shorter than begin_period_index, so this works
    for i in range(len(end_period_index)):
        begin_time = timeArr[begin_period_index[i]]
        end_time = timeArr[end_period_index[i]]
        duration = abs(end_time - begin_time) #calculate the duration of the period
        filteredEnergies = validEnergies[begin_period_index[i]:end_period_index[i]] # filter energies for the current period
        filteredEnergies = filteredEnergies[(1459. <= filteredEnergies) & (filteredEnergies <= 1461.)] # filter energies in range [1459, 1461]
        if filteredEnergies.size == 0 or not np.isclose(duration, sidereal_day, atol=100): #if no events in the energy range, or period too long, or somehow negative, skip
            continue
         #append the duration of the period
        nEvents.append(filteredEnergies.size) #count number of events in the energy range [1459, 1461] in the current period
        durationArr.append(duration) #append the duration of the period
        inPhaseIntegral.append(inPhase.IntegralOneDim(begin_time, end_time, 1e-6, 1e-6, err))
        quadIntegral.append(quadrature.IntegralOneDim(begin_time, end_time, 1e-6, 1e-6, err))
        
    inPhaseIntegral = np.array(inPhaseIntegral)
    quadIntegral = np.array(quadIntegral)   
    durationArr = np.array(durationArr).astype(float) #convert to float for division


    iqPhase = np.arctan2(inPhaseIntegral, quadIntegral) #calculate phase using arctan2
    #arctan will cancel out constants leaving only phase
    amplitude = np.hypot(nEvents * quadIntegral, nEvents * inPhaseIntegral)

    print("\n\nThe total number of events in Jul-Aug 2020 dataset is: ", total)
    print("The number of events in the energy range [1459, 1461] is: ", nEvents)
    print("\nThe average phase is: ", np.average(iqPhase))
    print("The variation in phase is: ", np.std(iqPhase), "radians\n")
    print("The average amplitude is: ", np.average(amplitude), "keV\n")
    print("The variation in amplitude is: ", np.std(amplitude), "keV\n")
    print("duration of each period is: ", durationArr, "seconds")