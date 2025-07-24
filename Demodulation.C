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
#include <math.h>

//this function is used for weighting based on time
double sin2weight(ULong64_t time) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time_double = (double)((time - 1594676438368)/ 1000); //time shift array here
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
    // auto chain = new TChain("qredtree");
    // chain->Add("SuperReduced_Background_ds3808.root"); //Jul-Aug 2019
    // chain->Add("SuperReduced_Background_ds3814.root"); //Jul-Aug 2020
    // chain->Add("SuperReduced_Background_ds3824.root"); //Jun-Aug 2022
    // ROOT::RDataFrame df(*chain); //special rdf constructor for TChain

    // auto filtered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    // auto withSin2Weights = filtered.Define("Sin2Weighting", sin2weight, {"EventTime"}); //defines new column Sin2Weighting using EventTime

    //try single dataset for now
    ROOT::RDataFrame ex("qredtree", "SuperReduced_Background_ds3814.root");
    auto filtered2 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    auto withSin2Weights2 = filtered2.Define("Sin2Weighting", sin2weight, {"EventTime"});

    auto timeArr = withSin2Weights2.Take<ULong64_t>("EventTime"); //not time adjusted yet
    auto energyArr = withSin2Weights2.Take<double>("Energy");
    auto weightArr = withSin2Weights2.Take<double>("Sin2Weighting");

    double startTime = (*timeArr)[0];
    double last = ((*timeArr)[timeArr->size() - 1] - startTime) / 1000;
    double energyLow = 1300;
    double energyHigh = 1500;

    //weighting by sin squared function
    int totalArrLen = timeArr->size();
    auto g2d = new TGraph2D(totalArrLen);
    for (int i = 0; i < totalArrLen; i++) {
        if (energyLow < (*energyArr)[i] && (*energyArr)[i] < energyHigh) {
            //Take method returns pointer, so dereference first
            double fixTime = ((*timeArr)[i] - startTime) / 1000;
            if (fixTime < 0) {
                std::cout << fixTime << " " << std::flush;
            }// REMOVE SINCE ONLY TEST
            g2d->SetPoint(i, fixTime, (*energyArr)[i], (*weightArr)[i]); //reminder! og time array still not adjusted
        }
    }
    g2d->Draw("PCOL0"); //plot markers with color, no stats box
    g2d->GetYaxis()->SetRangeUser(1320, 1480); 
    //slightly less than energyLow and High to avoid fXmax error
    g2d->GetXaxis()->SetTitleOffset(2.5);
    g2d->GetYaxis()->SetTitleOffset(2);
    g2d->SetTitle("Energy Weighted vs Time (Summer 2020); Adjusted Time (s);Energy (keV);Sin Squared Weighting");
    g2d->SetMarkerStyle(8); //changes from default square to large dot
    g2d->SetMarkerSize(0.4);

    c->Update();
    c->SaveAs("energyScatter.pdf");

    //extract phase and amp
    std::vector<long double> eventArr; //dynamically sized std::vector arrays
    std::vector<long double> nrgArr;
    std::vector<long double> inPhaseArr;
    std::vector<long double> quadArr;

    int siderealDay = 3600*23 + 56*60 + 4;
    int currTime = 0;
    int endTime = 0;
    //count #events & avg energy for small sidereal duration
    while (endTime < last) {
        endTime += siderealDay;
        std::vector<double> nrgAvg;
        int nEvents = 0;
        for (int i = 0; i < totalArrLen; i++) {
            //checks if valid for current sidereal period
            if ((1459 <= (*energyArr)[i] && (*energyArr)[i] <= 1461)  && (currTime <= (((*timeArr)[i] - startTime) / 1000) && (((*timeArr)[i] - startTime) / 1000) < endTime )) {
                //shift time vals in if clause
                nEvents++;
                nrgAvg.push_back((*energyArr)[i]); //adds it to dynamic sized array
            }
        }
        eventArr.push_back(nEvents);
        double avgNrg = std::accumulate(nrgAvg.begin(), nrgAvg.end(), 0.0) / nrgAvg.size(); //should be 1460 ish

        if (!TMath::IsNaN(avgNrg)) {//only if there are events in the period
            nrgArr.push_back(avgNrg);

            Double_t allowedErr = 5.;
            Double_t *err = &allowedErr;
            //creates TF1 object from 0 to last
            auto inPhase = new TF1("f1", cosDemod, 0, last);
            double integVal = inPhase->Integral(currTime, endTime, 1.E-6); //relative error
            inPhaseArr.push_back(integVal);

            auto quadrature = new TF1("f2", sinDemod, 0, last);
            double integVal2 = quadrature->Integral(currTime, endTime, 1.E-6);
            quadArr.push_back(integVal2);
        }
        currTime = endTime; //move other pointer forward
    }
    int totalROI = std::accumulate(eventArr.begin(), eventArr.end(), 0.0);
    double totalROIper = totalROI / (double) eventArr.size(); //total per sidereal period

    std::vector<long double> phaseArr;
    std::vector<long double> ampArr;
    for (int i = 0; i < eventArr.size(); i++) {
        long double iqPhase = TMath::ATan2(quadArr[i], inPhaseArr[i]);
            //arctan will cancel out constants
        phaseArr.push_back(iqPhase);
        long double amplitude = sqrt((eventArr[i] * inPhaseArr[i]) * (eventArr[i] * inPhaseArr[i]) + (eventArr[i] * quadArr[i]) * (eventArr[i] * quadArr[i]));
        //this lets individual #events be calculated with each integral result, more accurate numbers
        ampArr.push_back(amplitude);
    }
    double finalPhase = (std::accumulate(phaseArr.begin(), phaseArr.end(), 0.0)) / (double) phaseArr.size();
    double finalAmp = (std::accumulate(ampArr.begin(), ampArr.end(), 0.0)) / (double) ampArr.size();
    double avgEnerg = (std::accumulate(nrgArr.begin(), nrgArr.end(), 0.0)) / nrgArr.size();

    std::cout << "\n\nThe total number of events in Jul-Aug 2020 dataset is: " << totalArrLen << std::endl;
    std::cout << "The total number of events in the energy range [1459, 1461] is: " << totalROI << std::endl;
    std::cout << "The average number of events in range per sidereal day is: " << totalROIper << std::endl;

    std::cout << "\nThe average phase is: " << std::fixed << std::setprecision(25) << finalPhase << std::endl;
    std::cout << "The standard deviation of the phase is: " << std::fixed << std::setprecision(25) << TMath::StdDev(phaseArr.begin(), phaseArr.end()) << "\n" << std::endl; 

    std::cout << "The average amplitude is: " << std::fixed << std::setprecision(25) << std::scientific << finalAmp << std::endl;
    std::cout << "The standard deviation of the amplitude is: " << std::fixed << std::setprecision(25) << std::scientific << TMath::StdDev(ampArr.begin(), ampArr.end()) << "\n" << std::endl; 

    std::cout << "The average energy in the range is: " << std::fixed << std::setprecision(25) << avgEnerg << std::endl;
    std::cout << "The standard deviation of the energy is: " << std::fixed << std::setprecision(25) << TMath::StdDev(nrgArr.begin(), nrgArr.end()) << std::endl;

    double siderealdaycount = last / siderealDay;
    std::cout << "\nThere are " << siderealdaycount << " sidereal days in this dataset" << std::endl;
    double solardaycount = (last / (3600*24));
    std::cout << "There are " << solardaycount << " solar days in this dataset\n\n" << std::endl;
}