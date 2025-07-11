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

    //weighting by sin squared function
    auto g2d = new TGraph2D(timeArr->size());
    for (int i = 0; i < timeArr->size(); i++) {
        g2d->SetPoint(i, (*timeArr)[i], (*energyArr)[i], (*weightArr)[i]);
    }
    g2d->GetYaxis()->SetRangeUser(1280, 1530);
    g2d->GetXaxis()->SetTitleOffset(2.5);
    g2d->GetYaxis()->SetTitleOffset(2);
    g2d->SetTitle("Energy Weighted vs Time (Summer 2020);Unix Time (ms);Energy (keV);Sin Squared Weighting");
    g2d->SetMarkerStyle(8); //changes from default square to large dot
    g2d->SetMarkerSize(0.4);
    g2d->Draw("PCOL0"); // plot markers with color, no stats box

    c->SaveAs("energyScatter.pdf");


    //use integral to extract phase and amp
    g2d->GetYaxis()->SetRangeUser(1459, 1461); //redefine to narrow gap instead
    int total = timeArr->size();
    int nEvents = 0;
    double *validEnergies = g2d->GetY();
    //counts number of events
    for (int i = 0; i < g2d->GetN(); ++i) {
        if (validEnergies[i] >= 1459 && validEnergies[i] <= 1461) {
            nEvents++;
        }
    }
    auto inPhase = new TF1("f1", cosDemod, (*timeArr)[0], (*timeArr)[timeArr->size() - 1]);
    auto quadrature = new TF1("f2", sinDemod, (*timeArr)[0], (*timeArr)[timeArr->size() - 1]);

    Double_t allowedErr = 0.;
    Double_t *err = &allowedErr;
    auto integ = inPhase->IntegralOneDim((*timeArr)[0], (*timeArr)[timeArr->size() - 1], 1.E-3, 1.E-3, *err); //1.E-3 is relaxed tolerance
    auto integ2 = quadrature->IntegralOneDim((*timeArr)[0], (*timeArr)[timeArr->size() - 1], 1.E-3, 1.E-3, *err);
    double iqPhase = TMath::ATan2(nEvents * integ2, nEvents * integ); 
    //arctan will cancel out constants leaving only phase
    double amplitude = TMath::Hypot(nEvents * integ, nEvents * integ2);

    std::cout << "\n\nThe total number of events in Jul-Aug 2020 dataset is: " << total << std::endl;
    std::cout << "The number of events in the energy range [1459, 1461] is: " << nEvents << std::endl;
    std::cout << "\nThe phase is: " << iqPhase << std::endl;
    std::cout << "The amplitude is: " << amplitude << "\n" << std::endl;

    
    //try snapshot to a temporary ROOT file with custom columns
    auto snapshot = withSin2Weights2.Snapshot("tempTree", "temp_snapshot.root", {"Energy", "EventTime", "Sin2Weighting"});
    TFile snap("temp_snapshot.root");
    auto tree = (TTree*)snap.Get("tempTree");
    //workaround to draw scatterplot with Draw() by not using RDataFrame

    auto d = new TCanvas();

    tree->Draw("Energy:EventTime", "Sin2Weighting"); 

    d->SaveAs("2dscatter.pdf");
}