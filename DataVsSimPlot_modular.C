// plotSimVsData.C
// ROOT macro to automate plotting of simulation vs data comparison

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

// Helper function to split string by delimiter
std::vector<std::string> SplitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Normalize both histograms to area = 1
void NormalizeHistograms(TH1D* h1, TH1D* h2) {
    if (h1->Integral() != 0) h1->Scale(1.0 / h1->Integral());
    if (h2->Integral() != 0) h2->Scale(1.0 / h2->Integral());
}

// Create ratio histogram (Data / Sim)
TH1D* MakeRatioHistogram(TH1D* hData, TH1D* hSim) {
    TH1D* hRatio = (TH1D*)hData->Clone("hRatio");
    hRatio->Divide(hSim);
    return hRatio;
}

// Plotting function
void PlotComparison(TH1D* hData, TH1D* hSim, std::string varName) {
    TCanvas* c1 = new TCanvas(Form("c_%s", varName.c_str()), Form("Simulation vs Data: %s", varName.c_str()), 900, 700);
    c1->Divide(1, 2);

    // Upper pad: main plot
    TPad* pad_u = (TPad*)c1->cd(1);
    pad_u->SetPad(0.0, 0.3, 1.0, 1.0);
    pad_u->SetTopMargin(0.1);
    pad_u->SetBottomMargin(0.0);
    pad_u->SetLeftMargin(0.12);
    pad_u->SetRightMargin(0.05);
    pad_u->Draw();
    pad_u->cd();

    hData->SetLineColor(kRed);
    hSim->SetLineColor(kBlue);
    hData->SetLineWidth(2);
    hSim->SetLineWidth(2);
    hData->SetTitle(Form("%s; %s; Normalized Yield", varName.c_str(), varName.c_str()));
    hData->Draw("HIST");
    hSim->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(hData, "Data", "l");
    leg->AddEntry(hSim,  "Simulation", "l");
    leg->Draw();

    // Lower pad: ratio plot
    c1->cd(2);
    TPad* pad_d = (TPad*)c1->cd(2);
    pad_d->SetPad(0.0, 0.0, 1.0, 0.3);
    pad_d->SetTopMargin(0.0);
    pad_d->SetBottomMargin(0.35);
    pad_d->SetLeftMargin(0.12);
    pad_d->SetRightMargin(0.05);
    pad_d->Draw();
    pad_d->cd();

    TH1D* hRatio = MakeRatioHistogram(hData, hSim);
    hRatio->SetLineColor(kBlack);
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Data / Sim");
    hRatio->GetXaxis()->SetTitle(varName.c_str());
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->Draw("E1");

    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw();

    c1->SaveAs(Form("compare_%s.pdf", varName.c_str()));
}

// MAIN FUNCTION
void DataVsSimPlot_modular() {
    int run_data, run_sim;
    std::cout << "Enter data run number: "; std::cin >> run_data;
    std::cout << "Enter simulation run number: "; std::cin >> run_sim;

    std::vector<std::string> variableList = {"H.gtr.dp", "H.kin.x_bj", "H.kin.Q2"};

    TFile* fData = TFile::Open(Form("ROOTfiles/data_%d.root", run_data));
    TTree* tData = (TTree*)fData->Get("T");

    TFile* fSim = TFile::Open(Form("ROOTfiles/sim_%d.root", run_sim));
    TTree* tSim = (TTree*)fSim->Get("T");

    for (auto var : variableList) {
        TH1D* hData = new TH1D("hData", "", 100, -1.0, 1.0);  // adjust range later
        TH1D* hSim  = new TH1D("hSim",  "", 100, -1.0, 1.0);

        std::string cutData = Form("%s*(1.0/charge_data*hms_eff_data)", var.c_str());
        std::string cutSim  = Form("%s*(normfac*sweight)", var.c_str());

        tData->Project("hData", var.c_str(), cutData.c_str());
        tSim->Project("hSim",  var.c_str(), cutSim.c_str());

        NormalizeHistograms(hData, hSim);
        PlotComparison(hData, hSim, var);

        delete hData;
        delete hSim;
    }
}
