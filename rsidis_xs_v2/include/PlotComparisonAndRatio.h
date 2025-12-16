#ifndef PLOT_COMPARISON_AND_RATIO_H
#define PLOT_COMPARISON_AND_RATIO_H

#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <algorithm>

#include "Mapping.h"

// PlotComparisonAndRatio
//
// Run command example:
//   root -l -b -q 'YourMacro.C'
//
// Inputs:
//   hData : data histogram after your subtractions (data-pos) - (dummy-posdummy)/wall
//   hSim  : simulation histogram (already normalized as you want)
//   varName : simulation branch name (for x-axis label via BranchToPhysicsMap)
//
// Output:
//   Creates directory "PNGs" in the current working directory if it does not exist,
//   then writes PNG images under outDir/PNGs (or ./PNGs if outDir is empty).

inline void PlotComparisonAndRatio(TH1D* hData, TH1D* hSim, const std::string& varName, const TString& outDir)
{
  if (!hData || !hSim) {
    std::cout << "PlotComparisonAndRatio: null histogram pointer\n";
    return;
  }

  TString pngDir = outDir.Length()>0 ? Form("%s/PNGs", outDir.Data()) : TString("PNGs");
  gSystem->mkdir(pngDir, kTRUE);
  gStyle->SetOptStat(0);

  double maxData = hData->GetMaximum();
  double maxSim  = hSim->GetMaximum();
  double ymax = 1.15 * std::max(maxData, maxSim);
  if (ymax <= 0.0) ymax = 1.0;

  TCanvas* c1 = new TCanvas(Form("c_%s", varName.c_str()),
                            Form("Data vs Simulation %s", varName.c_str()),
                            900, 750);

  TPad* pad_u = new TPad(Form("pad_u_%s", varName.c_str()), "", 0.0, 0.30, 1.0, 1.0);
  pad_u->SetTopMargin(0.07);
  pad_u->SetBottomMargin(0.02);
  pad_u->SetLeftMargin(0.12);
  pad_u->SetRightMargin(0.04);
  pad_u->SetTicks(1, 1);
  pad_u->Draw();
  pad_u->cd();

  hSim->SetLineWidth(2);
  hSim->SetLineColor(kBlue);
  hSim->SetFillStyle(0);
  hSim->SetTitle("");
  hSim->SetMaximum(ymax);
  hSim->SetMinimum(0.0);

  hSim->GetYaxis()->SetTitle("Counts");
  hSim->GetYaxis()->SetTitleSize(0.055);
  hSim->GetYaxis()->SetTitleOffset(0.95);
  hSim->GetYaxis()->SetLabelSize(0.045);
  hSim->GetYaxis()->SetMaxDigits(3);

  hSim->GetXaxis()->SetTitleSize(0.0);
  hSim->GetXaxis()->SetLabelSize(0.0);

  hSim->Draw("HIST");

  hData->SetLineWidth(2);
  hData->SetLineColor(kRed);
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(kRed);
  hData->SetTitle("");
  hData->Draw("E1 SAME");

  int maxBin = hData->GetMaximumBin();
  double peakX = hData->GetXaxis()->GetBinCenter(maxBin);

  TLegend* leg = new TLegend();
  double midX = 0.5 * (hData->GetXaxis()->GetXmin() + hData->GetXaxis()->GetXmax());
  if (peakX < midX) {
    leg->SetX1NDC(0.62); leg->SetX2NDC(0.94);
    leg->SetY1NDC(0.74); leg->SetY2NDC(0.90);
  } else {
    leg->SetX1NDC(0.16); leg->SetX2NDC(0.48);
    leg->SetY1NDC(0.74); leg->SetY2NDC(0.90);
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(hSim,  "Simulation", "l");
  leg->AddEntry(hData, "Data",       "lep");
  leg->Draw();

  c1->cd();

  TPad* pad_d = new TPad(Form("pad_d_%s", varName.c_str()), "", 0.0, 0.0, 1.0, 0.30);
  pad_d->SetTopMargin(0.03);
  pad_d->SetBottomMargin(0.38);
  pad_d->SetLeftMargin(0.12);
  pad_d->SetRightMargin(0.04);
  pad_d->SetTicks(1, 1);
  pad_d->Draw();
  pad_d->cd();

  TH1D* hRatio = (TH1D*)hData->Clone(Form("hRatio_%s", varName.c_str()));
  hRatio->SetDirectory(nullptr);
  hRatio->Divide(hSim);

  hRatio->SetTitle("");
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.7);
  hRatio->SetStats(0);

  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);

  hRatio->GetYaxis()->SetTitle("Data/Sim");
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.45);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetNdivisions(505);

  hRatio->GetXaxis()->SetTitle(BranchToPhysicsMap(varName).c_str());
  hRatio->GetXaxis()->SetTitleSize(0.13);
  hRatio->GetXaxis()->SetTitleOffset(1.15);
  hRatio->GetXaxis()->SetLabelSize(0.11);
  hRatio->GetXaxis()->SetNdivisions(506);

  hRatio->Draw("E1");

  TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                          hRatio->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw();

  c1->SaveAs(Form("%s/%s_comparison.png", pngDir.Data(), varName.c_str()));

  // Optional cross section hook for Q2
  if (varName == "Q2") {
    TCanvas* c2 = new TCanvas(Form("c_xs_%s", varName.c_str()),
                              Form("Cross section %s", varName.c_str()),
                              900, 700);

    c2->SetLeftMargin(0.16);
    c2->SetRightMargin(0.05);
    c2->SetBottomMargin(0.16);
    c2->SetTopMargin(0.06);
    c2->SetTicks(1, 1);

    double sigma_model = 1.012785e-03;

    TH1D* hSigma = (TH1D*)hData->Clone(Form("hSigma_%s", varName.c_str()));
    hSigma->SetDirectory(nullptr);
    hSigma->Divide(hSim);
    hSigma->Scale(sigma_model);

    hSigma->SetLineColor(kBlack);
    hSigma->SetLineWidth(2);
    hSigma->SetMarkerStyle(20);
    hSigma->SetMarkerSize(0.85);
    hSigma->SetStats(0);

    hSigma->GetXaxis()->SetTitle("Q2 [GeV2]");
    hSigma->GetXaxis()->SetTitleSize(0.055);
    hSigma->GetXaxis()->SetTitleOffset(1.25);
    hSigma->GetXaxis()->SetLabelSize(0.045);
    hSigma->GetXaxis()->SetNdivisions(506);

    hSigma->GetYaxis()->SetTitle("dSigma (model scaled)");
    hSigma->GetYaxis()->SetTitleSize(0.055);
    hSigma->GetYaxis()->SetTitleOffset(1.35);
    hSigma->GetYaxis()->SetLabelSize(0.045);
    hSigma->GetYaxis()->SetNdivisions(506);
    hSigma->GetYaxis()->SetMaxDigits(3);

    double y2max = 1.2 * hSigma->GetMaximum();
    if (y2max <= 0.0) y2max = 1.0;
    hSigma->GetYaxis()->SetRangeUser(0.0, y2max);

    hSigma->Draw("E1");
    c2->SaveAs(Form("%s/%s_cross_section.png", pngDir.Data(), varName.c_str()));
  }
}

#endif
