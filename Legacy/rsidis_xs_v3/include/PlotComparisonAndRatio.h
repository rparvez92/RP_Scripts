#ifndef PLOT_COMPARISON_AND_RATIO_H
#define PLOT_COMPARISON_AND_RATIO_H

// PlotComparisonAndRatio.h
// Helpers to produce comparison/ratio and xsec plots.
// Output: PNG files saved in <outDir>/PNGs
//
// No uncommon symbols in comments.

#include <cmath>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TStyle.h"

inline void EnsurePngDir(const TString& outDir)
{
  TString pngDir = outDir + "/PNGs";
  gSystem->mkdir(pngDir, true);
}

inline void SetCommonAxisStyle(TH1* h)
{
  if (!h) return;
  h->SetStats(0);

  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.95);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetNdivisions(506);

  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.95);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetNdivisions(506);
  h->GetYaxis()->SetMaxDigits(3);
}

inline void SetXAxisTitleForVar(TH1* h, const std::string& varName)
{
  if (!h) return;
  if (varName == "Q2")      h->GetXaxis()->SetTitle("Q2 [GeV2]");
  else if (varName == "W")  h->GetXaxis()->SetTitle("W [GeV]");
  else if (varName == "xbj") h->GetXaxis()->SetTitle("Bjorken x");
  else                      h->GetXaxis()->SetTitle(varName.c_str());
}

inline double ComputeYMaxWithErr(const TH1* h)
{
  if (!h) return 1.0;
  double ymax = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double y = h->GetBinContent(i);
    double e = h->GetBinError(i);
    if (!std::isfinite(y) || !std::isfinite(e)) continue;
    ymax = std::max(ymax, y + e);
  }
  if (!(ymax > 0.0)) ymax = 1.0;
  return ymax;
}

inline void PlotComparisonAndRatio(TH1D* hData, TH1D* hSim, const std::string& varName, const TString& outDir)
{
  if (!hData || !hSim) return;
  EnsurePngDir(outDir);
  TString pngDir = outDir + "/PNGs";

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(Form("c_%s_comp", varName.c_str()),
                           Form("%s comparison", varName.c_str()),
                           900, 800);
  c->SetTopMargin(0.06);

  TPad* p1 = new TPad("p1","p1",0.0,0.3,1.0,1.0);
  TPad* p2 = new TPad("p2","p2",0.0,0.0,1.0,0.3);
  p1->SetTopMargin(0.1);
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(0.1);
  p1->SetRightMargin(0.05);

  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.35);
  p2->SetLeftMargin(0.1);
  p2->SetRightMargin(0.05);
  //p2->SetBottomMargin(0.35);
  p1->Draw();
  p2->Draw();

  p1->cd();
  hSim->SetLineWidth(2);
  hSim->SetLineColor(kBlue);

  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.85);
  hData->SetLineColor(kRed);
  hData->SetMarkerColor(kRed);

  SetCommonAxisStyle(hSim);
  hSim->GetXaxis()->SetLabelSize(0.0);
  SetXAxisTitleForVar(hSim, varName);
  hSim->GetYaxis()->SetTitle("Counts");

  double ymax1 = std::max(ComputeYMaxWithErr(hSim), ComputeYMaxWithErr(hData));
  hSim->SetMaximum(1.20 * ymax1);
  hSim->SetMinimum(0.0);

  hSim->Draw("HIST");
  hData->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.15,0.70,0.40,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hSim, "Simulation", "l");
  leg->AddEntry(hData, "Data", "pe");
  leg->Draw();

  p2->cd();
  TH1D* hRatio = (TH1D*)hData->Clone(Form("hRatio_%s", varName.c_str()));
  hRatio->Divide(hSim);
  SetCommonAxisStyle(hRatio);
  SetXAxisTitleForVar(hRatio, varName);
  hRatio->GetYaxis()->SetTitle("Data/Sim");
  hRatio->GetYaxis()->SetTitleOffset(0.4);
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetLabelSize(0.1);

  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetXaxis()->SetTitleOffset(1.0);

  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(0.75);
  hRatio->SetLineColor(kBlack);

  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);

  hRatio->Draw("E1");

  TLine* l1 = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                        hRatio->GetXaxis()->GetXmax(), 1.0);
  l1->SetLineStyle(2);
  l1->Draw();

  c->SaveAs(Form("%s/%s_comparison.png", pngDir.Data(), varName.c_str()));
}

inline void PlotXsec(TH1D* hDataFinal, TH1D* hSim, double sigma_model_ub,
                     const std::string& varName, const TString& outDir)
{
  if (!hDataFinal || !hSim) return;
  EnsurePngDir(outDir);
  TString pngDir = outDir + "/PNGs";

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(Form("c_%s_xsec", varName.c_str()),
                           Form("%s xsec", varName.c_str()),
                           900, 800);
  c->SetTopMargin(0.07);
  c->SetBottomMargin(0.1);
  c->SetRightMargin(0.07);
  c->SetLeftMargin(0.1);

  TH1D* hRatio = (TH1D*)hDataFinal->Clone(Form("hRatio_xsec_%s", varName.c_str()));
  hRatio->Divide(hSim);

  TH1D* hXsec = (TH1D*)hRatio->Clone(Form("hXsec_%s", varName.c_str()));
  hXsec->Scale(sigma_model_ub);

  hXsec->SetMarkerStyle(20);
  hXsec->SetMarkerSize(0.85);
  hXsec->SetLineColor(kBlack);

  SetCommonAxisStyle(hXsec);
  SetXAxisTitleForVar(hXsec, varName);
  hXsec->GetYaxis()->SetTitle("dSigma");
  hXsec->GetYaxis()->SetTitleOffset(1.1);
  hXsec->GetYaxis()->SetTitleSize(0.035);
  hXsec->GetYaxis()->SetLabelSize(0.030);

  hXsec->GetXaxis()->SetTitleOffset(0.95);
  hXsec->GetXaxis()->SetTitleSize(0.035);
  hXsec->GetXaxis()->SetLabelSize(0.030);
  double ymax = ComputeYMaxWithErr(hXsec);
  hXsec->SetMinimum(0.0);
  hXsec->SetMaximum(1.25 * ymax);

  hXsec->Draw("E1");
  c->SaveAs(Form("%s/%s_xsec.png", pngDir.Data(), varName.c_str()));
}

#endif
