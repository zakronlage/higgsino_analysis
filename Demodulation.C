#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>
#include "ROOT/RVec.hxx"
#include "TGraph2D.h"

double sin2weight(ULong64_t time) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time_double = (double)time;
    return (1 - TMath::Cos(2 * omega * time_double)) / 2; //equivalent to sin squared weighting
}

void Demodulation() {
    auto c = new TCanvas();

    //vertically combined dataset
    auto chain = new TChain("qredtree"); //all called qredtree, with cycle # after
    chain->Add("SuperReduced_Background_ds3808.root"); //Data from Jul-Aug 2019
    chain->Add("SuperReduced_Background_ds3814.root"); //Jul-Aug 2020
    chain->Add("SuperReduced_Background_ds3824.root"); //Jun-Aug 2022
    ROOT::RDataFrame df(*chain); //special rdf constructor for TChain

    auto collectFiltered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");

    auto withCosWeights = collectFiltered.Define("CosWeighting", sin2weight, {"EventTime"}); //defines new column CosWeighting using EventTime, called CosWeight
    // auto cosHisto = withCosWeights.Histo1D({"e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "CosWeighting");
    // cosHisto->SetLineColor(kBlue);
    // cosHisto->DrawClone();
    
    // withCosWeights.Draw("Energy : CosWeighting", "CosWeighting");

    ROOT::RDataFrame ex("qredtree", "SuperReduced_Background_ds3814.root");
    auto collectFiltered2 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    // auto filterRange = collectFiltered2->GetXaxis()->SetRangeUser(1400, 1600);
    auto withCosWeights2 = collectFiltered2.Define("CosWeighting2", sin2weight, {"EventTime"});

    auto time = withCosWeights2.Take<ULong64_t>("EventTime"); 
    auto energy = withCosWeights2.Take<double>("Energy");  
    auto weighting = withCosWeights2.Take<double>("CosWeighting2"); 

    auto g2 = new TGraph2D(time->size()); //of time size

    for (size_t i = 0; i < time->size(); ++i) {
        g2->SetPoint(i, (*time)[i], (*energy)[i], (*weighting)[i]); 
    }
    g2->GetYaxis()->SetRangeUser(1500, 1600);

    g2->SetTitle("Energy vs Time (Color = Weight);Time;Energy;CosWeighting");
    g2->SetMarkerStyle(20); 
    g2->SetMarkerSize(0.8);

    // auto c = new TCanvas("c", "Energy vs Time Colored by Weight", 800, 600);
    g2->Draw("PCOL"); //colored scatter plot


    // auto withCosWeights2 = collectFiltered2.Define("CosWeighting2", cosweight, {"EventTime"});

    //try snapshot to a temporary ROOT file
    auto snapshot = withCosWeights2.Snapshot("tempTree", "temp_snapshot.root", {"Energy", "EventTime", "CosWeighting2"});

    TFile snap("temp_snapshot.root");
    auto tree = (TTree*)snap.Get("tempTree");

    // tree->Draw("Energy:EventTime", "CosWeighting2"); 
    //should make sure to do all things using same dataset at end since lazy, better run time


    auto withSinWeights = collectFiltered.Define("SinWeighting", sin2weight, {"EventTime"});
    // auto sinHisto = withSinWeights.Histo1D({"e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "SinWeighting");
    // sinHisto->SetLineColor(kMagenta);
    // sinHisto->DrawClone("Same");

    auto legend = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend->SetHeader("Legend");
    // legend->AddEntry(&*cosHisto, "Energy Weighted (Cos)", "l");
    // legend->AddEntry(&*sinHisto, "Energy Weighted (Sin)", "l");
    legend->Draw();

    // auto string = withSinWeights.Display();
    // string->Print();

    c->SaveAs("energyScatter.pdf");
    //use Draw("X:Y") for scatterplot //energy as y, EventTime as x since energy already modulated

    /*
    //k40 peak combined cos only
    auto c1 = new TCanvas();
    auto histo = withCosWeights.Histo1D({"e", "K-40 Peak Cosine Weighting (Combined)", 75, 1400., 1550.}, "Energy", "CosWeighting");
    histo->SetLineColor(kBlue);
    histo->DrawClone();

    auto legend2 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend2->SetHeader("Legend");
    legend2->AddEntry(&*histo, "Weighted Energy", "l");
    legend2->Draw();

    c1->SaveAs("k40combinedCosOnly.pdf");


    //k40 peak combined sin only
    auto c2 = new TCanvas();
    auto histo2 = withSinWeights.Histo1D({"e", "K-40 Peak Sine Weighted (Combined)", 75, 1400., 1550.}, "Energy", "SinWeighting");
    histo2->SetLineColor(kMagenta);
    histo2->DrawClone();

    auto legend3 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend3->SetHeader("Legend");
    legend3->AddEntry(&*histo2, "Weighted Energy", "l");
    legend3->Draw();

    c2->SaveAs("k40combinedSinOnly.pdf");


    //3801 dataset only
    ROOT::RDataFrame ex("qredtree;5", "SuperReduced_Background_ds3801.root");
    auto c3 = new TCanvas();

    auto appliedCut5 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1 && NumberOfPulses == 1");
    auto withCosWeights2 = appliedCut5.Define("CosWeighting2", cosweight, {"EventTime"});
    auto histo3 = withCosWeights2.Histo1D({"e", "Singular Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "CosWeighting2");
    histo3->SetLineColor(kBlue);
    histo3->DrawClone();

    auto withSineWeights2 = appliedCut5.Define("SinWeighting2", sinweight, {"EventTime"});
    auto histo4 = withSineWeights2.Histo1D({"e", "Singular Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "SinWeighting2");
    histo4->SetLineColor(kMagenta);
    histo4->DrawClone("Same");

    auto legend4 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend4->SetHeader("Legend");
    legend4->AddEntry(&*histo3, "Energy Weighted (Cosine)", "l");
    legend4->AddEntry(&*histo4, "Energy Weighted (Sine)", "l");
    legend4->Draw();

    c3->SaveAs("cosSinWeighted.pdf");

    //toy mc in progress
    auto histo4 = new TH1D("test", "", 20);
    for(int i=0; i < 1650; i++) {
        //histo4->Fill(histo3->GetRandom());
    }
    use cos weight first to make weighted histogram for this random sample
    then sin weight one

    //3801 k-40 peak cosine only
    auto c4 = new TCanvas();
    auto histo5 = withCosWeights2.Histo1D({"eB", "K-40 Peak Cosine Weighted (Singular)", 75, 1400., 1550.}, "Energy", "CosWeighting2");
    histo5->SetLineColor(kBlue);
    histo5->DrawClone();

    auto legend5 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend5->SetHeader("Legend");
    legend5->AddEntry(&*histo5, "Weighted Energy", "l");
    legend5->Draw();

    c4->SaveAs("k40cosonly.pdf");


    //3801 k-40 peak sine only
    auto c5 = new TCanvas();
    auto histo6 = withSineWeights2.Histo1D({"eB", "K-40 Peak Sine Weighted (Singular)", 75, 1400., 1550.}, "Energy", "SinWeighting2");
    histo6->SetLineColor(kMagenta);
    histo6->DrawClone();

    auto legend6 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend6->SetHeader("Legend");
    legend6->AddEntry(&*histo6, "Weighted Energy", "l");
    legend6->Draw();

    c5->SaveAs("k40sinonly.pdf");
    */
}