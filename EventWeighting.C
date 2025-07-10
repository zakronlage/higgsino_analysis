#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>

double cosweight(ULong64_t time) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4); //angular frequency of sidereal day
    double time_double = (double)time;
    return TMath::Cos(omega * time_double);
}

double sinweight(ULong64_t time) {
    double omega = (2 * TMath::Pi()) / (3600*23 + 56*60 + 4);
    double time_double = (double)time;
    return TMath::Sin(omega * time_double);
}

void EventWeighting() {
    auto c = new TCanvas();

    //vertically combined dataset
    auto chain = new TChain("qredtree"); //all called qredtree, with cycle # after
    chain->Add("SuperReduced_Background_ds3808.root"); //Data from Jul-Aug 2019
    chain->Add("SuperReduced_Background_ds3814.root"); //Jul-Aug 2020
    chain->Add("SuperReduced_Background_ds3824.root"); //Jun-Aug 2022
    ROOT::RDataFrame df(*chain); //special rdf constructor for TChain

    auto collectFiltered = df.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1");

    auto withCosWeights = collectFiltered.Define("CosWeighting", cosweight, {"EventTime"}); //defines new column CosWeighting using EventTime, called CosWeight
    auto cosHisto = withCosWeights.Histo1D({"e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "CosWeighting");
    
    cosHisto->SetDefaultSumw2();
    cosHisto->SetLineColor(kBlue);
    cosHisto->DrawClone();

    auto withSinWeights = collectFiltered.Define("SinWeighting", sinweight, {"EventTime"});
    auto sinHisto = withSinWeights.Histo1D({"e", "Combined Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "SinWeighting");
    sinHisto->SetLineColor(kMagenta);
    sinHisto->DrawClone("Same");

    auto legend = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend->SetHeader("Legend");
    legend->AddEntry(&*cosHisto, "Energy Weighted (Cos)", "l");
    legend->AddEntry(&*sinHisto, "Energy Weighted (Sin)", "l");
    legend->Draw();

    c->SaveAs("dsCombined.pdf");


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

    auto appliedCut4 = ex.Filter("AnalysisBaseCut && GoodForAnalysis && PCACut && Multiplicity == 1");
    auto withCosWeights2 = appliedCut4.Define("CosWeighting2", cosweight, {"EventTime"});
    auto histo3 = withCosWeights2.Histo1D({"e", "Singular Weighted Energy Spectrum", 5000, 0., 10000.}, "Energy", "CosWeighting2");
    histo3->SetLineColor(kBlue);
    histo3->DrawClone();

    auto withSineWeights2 = appliedCut4.Define("SinWeighting2", sinweight, {"EventTime"});
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
    /*
    auto histo4 = new TH1D("test", "", 20);
    for(int i=0; i < 1650; i++) {
        //histo4->Fill(histo3->GetRandom());
    }
    use cos weight first to make weighted histogram for this random sample
    then sin weight one
    */

    //3801 k-40 peak cosine only
    auto c4 = new TCanvas();
    auto cos_k40Peak = withCosWeights2.Histo1D({"eB", "K-40 Peak Cosine Weighted (Singular)", 75, 1400., 1550.}, "Energy", "CosWeighting2");
    cos_k40Peak->SetLineColor(kBlue);
    cos_k40Peak->DrawClone();

    auto legend5 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend5->SetHeader("Legend");
    legend5->AddEntry(&*cos_k40Peak, "Weighted Energy", "l");
    legend5->Draw();

    c4->SaveAs("k40cosonly.pdf");


    //3801 k-40 peak sine only
    auto c5 = new TCanvas();
    auto sin_k40Peak = withSineWeights2.Histo1D({"eB", "K-40 Peak Sine Weighted (Singular)", 75, 1400., 1550.}, "Energy", "SinWeighting2");
    sin_k40Peak->SetLineColor(kMagenta);
    sin_k40Peak->DrawClone();

    auto legend6 = new TLegend(0.65, 0.75, 1.0, 0.94);
    legend6->SetHeader("Legend");
    legend6->AddEntry(&*sin_k40Peak, "Weighted Energy", "l");
    legend6->Draw();

    c5->SaveAs("k40sinonly.pdf");

    TCanvas* c6 = new TCanvas();
    TH1D* ratio = ((TH1D*)sinHisto->Clone("ratio"));
    ratio->Divide(cosHisto.GetPtr());
    ratio->SetTitle("Sine/Cosine Ratio");
    ratio->SetLineColor(kRed);
    ratio->Draw("pE");

    double sin_events = sinHisto->Integral();
    double cos_events = cosHisto->Integral();
    double ratio_value = sin_events / cos_events;
    printf("ratio = %f / %f = %.2f\n", sin_events, cos_events, ratio_value);
    printf("phase = %.2f degrees\n", TMath::ATan(ratio_value) * 180.0 / TMath::Pi());

    TCanvas* c7 = new TCanvas();
    TH1D* ratio2 = ((TH1D*)sin_k40Peak->Clone("ratio2"));
    ratio2->Divide(cos_k40Peak.GetPtr());
    ratio2->SetTitle("Sine/Cosine Ratio (Singular)");
    ratio2->SetLineColor(kRed);
    ratio2->Draw("E");
}