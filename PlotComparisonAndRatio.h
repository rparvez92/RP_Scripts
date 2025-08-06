#ifndef PLOT_COMPARISON_AND_RATIO_H
#define PLOT_COMPARISON_AND_RATIO_H

#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <string>


// Plotting function
// Use h1 and h2 as hSim and hDataSubDummy respectively
void PlotComparisonAndRatio(TH1D* h1, TH1D* h2, std::string varName) {
    // Find the max values from both histograms
    double max1 = h1 -> GetMaximum();
    double max2 = h2 -> GetMaximum();

    // Find the higher of the two and scale up slightly for margin
    double ymax = std::max(max1, max2) * 1.1;

    //Draw the canvas
    TCanvas* c1 = new TCanvas(Form("c_%s", varName.c_str()), Form("Simulation vs Data: %s", varName.c_str()), 900, 700);
    c1->Divide(1, 2);

    // Upper pad: comparison plot
    TPad* pad_u = (TPad*)c1->cd(1);
    pad_u->SetPad(0.0, 0.3, 1.0, 1.0);
    pad_u->SetTopMargin(0.1);
    pad_u->SetBottomMargin(0.0);
    pad_u->SetLeftMargin(0.1);
    pad_u->SetRightMargin(0.05);
    pad_u->Draw();
    pad_u->cd();

    h1->SetMaximum(ymax); //Set the y-axis range of the first drawn histogram
    h1->SetLineColor(kBlue);
    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetTitleSize(0);
    h1->GetYaxis()->SetTitle("Counts");
    h1->Draw("HIST");

    h2->SetLineColor(kRed);
    h2->GetXaxis()->SetLabelSize(0);
    h2->GetXaxis()->SetTitleSize(0);
    //h2->SetTitle(Form("%s; %s; Normalized Yield", varName.c_str(), varName.c_str()));
    h2->Draw("HIST SAME");

    // Draw Legend
    auto legend = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend->AddEntry(h1, "Simulation", "l");
    legend->AddEntry(h2, "Data-Dummy", "l");
    legend->Draw();

    // Lower pad: ratio plot
    c1->cd(2);
    TPad* pad_d = (TPad*)c1->cd(2);
    pad_d->SetPad(0.0, 0.0, 1.0, 0.3);
    pad_d->SetTopMargin(0.0);
    pad_d->SetBottomMargin(0.35);
    pad_d->SetLeftMargin(0.1);
    pad_d->SetRightMargin(0.05);
    pad_d->Draw();
    pad_d->cd();

    // Make ratio histogram
    TH1D* hRatio = (TH1D*) h1 -> Clone("hRatio"); // clone to avoid modifying original
    hRatio->Divide(h2); // Sim / (Data - Dummy)

    hRatio->SetLineColor(kBlack);
    hRatio->SetStats(0);
    hRatio->SetTitle("");

    // Set the range of y-axis
    hRatio->SetMinimum(0.5);
    hRatio->SetMaximum(1.5);

    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetYaxis()->SetTitleOffset(0.5);
    hRatio->GetYaxis()->SetLabelSize(0.1);
    hRatio->GetYaxis()->SetNdivisions(500);

    hRatio->GetXaxis()->SetTitle(varName.c_str());
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetNdivisions(505);


    // Draw with error bar
    hRatio->Draw("E1");

    // Fix for invisible y-axis tick labels
    gPad->Update(); //Forces ROOT to layout the frame
    hRatio->GetYaxis()->SetLabelOffset(0.01); //Small offset
    hRatio->GetYaxis()->SetNdivisions(505); //5 major, 5 minor ticks
    gPad->RedrawAxis();

    // Draw a dashed line at y=1 
    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->Draw();

    //c1->SaveAs(Form("compare_%s.pdf", varName.c_str()));
    // Show the plots interactively
    c1->Update();
    gSystem->ProcessEvents();
    std::cout << "Close the canvas to continue...\n";
    while (gROOT->GetListOfCanvases()->FindObject(c1)) {
        gSystem->ProcessEvents();
        gSystem->Sleep(100);
    }

}

#endif // PLOT_COMPARISON_AND_RATIO_H
