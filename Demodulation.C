#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>
#include <ROOT/RVec.hxx>
#include <TGraph2D.h>
#include <numeric>

//this function is used for weighting based on time
double sin2weight(ULong64_t time) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time_double = (double)time;
    double phase = 0;
    return (1 - TMath::Cos(2 * omega * (time_double - phase))) / 2; //equivalent to sin squared weighting
}

//this function is used for integration to demodulate and extract phase & amplitude
double cosDemod(double *time, double *params) { //requires params argument even if unused since integration
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time1 = time[0];
    double phase = 0;
    return TMath::Cos(2 * omega * (time1 - phase));
}

double sinDemod(double *time, double *params) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time1 = time[0];
    double phase = 0;
    return TMath::Sin(2 * omega * (time1 - phase));
}

void Demodulation() {
    auto c = new TCanvas();

    //vertically combined dataset
    auto chain = new TChain("qredtree");
    chain->Add("SuperReduced_Background_ds3808.root"); //Jul-Aug 2019
    chain->Add("SuperReduced_Background_ds3814.root"); //Jul-Aug 2020
    chain->Add("SuperReduced_Background_ds3824.root"); //Jun-Aug 2022
    ROOT::RDataFrame df(*chain); //special rdf constructor for TChain

    auto filtered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    auto withSin2Weights = filtered.Define("Sin2Weighting", sin2weight, {"EventTime"}); //defines new column Sin2Weighting using EventTime

    //try single dataset for now
    ROOT::RDataFrame ex("qredtree", "SuperReduced_Background_ds3814.root");
    auto filtered2 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    auto withSin2Weights2 = filtered2.Define("Sin2Weighting", sin2weight, {"EventTime"});

    auto timeArr = withSin2Weights2.Take<ULong64_t>("EventTime");
    auto energyArr = withSin2Weights2.Take<double>("Energy");
    auto weightArr = withSin2Weights2.Take<double>("Sin2Weighting");
    double startTime = (*timeArr)[0];
    double energyLow = 1300;
    double energyHigh = 1500;
    int total = timeArr->size();

    //weighting by sin squared function
    auto g2d = new TGraph2D(total);
    for (int i = 0; i < total; i++) {
        if (energyLow < (*energyArr)[i] && (*energyArr)[i] < energyHigh) {
            //Take method returns pointer, so dereference first
            double fixTime = ((*timeArr)[i] - startTime) / 1000;
            g2d->SetPoint(i, fixTime, (*energyArr)[i], (*weightArr)[i]);
        }
    }
    g2d->Draw("PCOL0"); // plot markers with color, no stats box
    g2d->GetYaxis()->SetRangeUser(1320, 1480); 
    //slightly less than energyLow and High to avoid fXmax error
    g2d->GetXaxis()->SetTitleOffset(2.5);
    g2d->GetYaxis()->SetTitleOffset(2);
    g2d->SetTitle("Energy Weighted vs Time (Summer 2020);Unix Time (s);Energy (keV);Sin Squared Weighting");
    g2d->SetMarkerStyle(8); //changes from default square to large dot
    g2d->SetMarkerSize(0.4);

    c->Update();
    c->SaveAs("energyScatter.pdf");

    //extract phase and amp
    std::vector<int> eventArr; //dynamically sized std::vector arrays
    std::vector<int> nrgArr;
    std::vector<double> inPhaseArr;
    std::vector<double> quadArr;
    int siderealDay = 3600*23 + 56*60 + 4;
    int currTime = 0;
    int endTime = 0;

    while (endTime < timeArr->size()) {
        endTime += siderealDay;
        int nEvents = 0;
        std::vector<double> nrgAvg;
        for (int i = currTime; i < endTime; i++) { //figures stuff out for each sidereal period
            if (1459 <= (*energyArr)[i] && (*energyArr)[i] <= 1461) {
                nEvents++;
                nrgAvg.push_back((*energyArr)[i]); //adds it to dynamic sized array
            }
        }
        eventArr.push_back(nEvents);

        int avgNrg = std::accumulate(nrgAvg.begin(), nrgAvg.end(), 0);
        nrgArr.push_back(avgNrg / nrgAvg.size());

        Double_t allowedErr = 0.;
        Double_t *err = &allowedErr;
        //creates TF1 object from [0] to [total]
        auto inPhase = new TF1("f1", cosDemod, (*timeArr)[0], (*timeArr)[total]);
        double integVal = inPhase->IntegralOneDim(currTime, endTime, 1.E-6, 1.E-3, *err); //1.E-3 is relaxed tolerance;
        inPhaseArr.push_back(integVal);

        auto quadrature = new TF1("f2", sinDemod, (*timeArr)[0], (*timeArr)[total]);
        double integVal2 = quadrature->IntegralOneDim(currTime, endTime, 1.E-6, 1.E-3, *err);
        quadArr.push_back(integVal2);

        currTime = endTime; //move other pointer forward
    }
    int totalROI = std::accumulate(eventArr.begin(), eventArr.end(), 0);
    int totalROIper = totalROI / eventArr.size(); //total per sidereal period

    double integRes = std::accumulate(inPhaseArr.begin(), inPhaseArr.end(), 0);
    double integRes2 = std::accumulate(quadArr.begin(), quadArr.end(), 0);

    double iqPhase = TMath::ATan2(integRes2, integRes);
    //arctan will cancel out constants leaving only phase
    double amplitude = TMath::Hypot(totalROIper * integRes, totalROIper * integRes2);

    std::cout << "\n\nThe total number of events in Jul-Aug 2020 dataset is: " << total << std::endl;
    std::cout << "The total number of events in the energy range [1459, 1461] is: " << totalROI << std::endl;
    std::cout << "The average number of events in range per sidereal day is: " << totalROIper << std::endl;
    std::cout << "\nThe average phase is: " << iqPhase << std::endl;
    std::cout << "The average amplitude is: " << amplitude << "\n\n" << std::endl;
}