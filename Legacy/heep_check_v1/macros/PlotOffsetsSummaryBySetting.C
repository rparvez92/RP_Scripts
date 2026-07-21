// PlotOffsetsSummaryBySetting.C
//
// Purpose:
//   Read OffsetsSummaryBySetting.csv and plot weighted-mean measured offsets
//   vs dp bin center, separately for Setting A and Setting B.
//
// Input:
//   results/tables/OffsetsSummaryBySetting.csv
//
// Output PNGs:
//   results/PNGs/PlotOffsetsSummaryBySetting/offset_summary_settingA.png
//   results/PNGs/PlotOffsetsSummaryBySetting/offset_summary_settingB.png
//
// Layout:
//   One PNG per setting, with 4 pads:
//     upper left:  W
//     upper right: Em
//     lower left:  Pmz
//     lower right: Pmy
//
// Plot rules:
//   - Narrow bins only are plotted as points at dp bin centers
//   - Full bin is not plotted as a point
//   - Full bin weighted mean is shown as a blue dashed horizontal line
//   - Full-bin uncertainty is shown as a transparent blue band
//
// Example:
//   root -l -b -q 'macros/PlotOffsetsSummaryBySetting.C+("results/tables/OffsetsSummaryBySetting.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TBox.h"

struct Row {
  TString setting = "";
  int dp_idx = -999;
  TString dp_label = "";
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  TString var = "";

  int n_valid = 0;
  double weighted_mean_dmu_MeV = 0.0;
  double weighted_mean_err_dmu_MeV = 0.0;
};

static std::vector<TString> SplitCSV(const std::string& line) {
  std::vector<TString> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(TString(item));
  return out;
}

static bool IsFinite(double x) {
  return std::isfinite(x);
}

static TString PrettyVar(const TString& v) {
  if (v == "W")   return "W";
  if (v == "Em")  return "Em";
  if (v == "Pmz") return "Pmz";
  if (v == "Pmy") return "Pmy";
  return v;
}

static void ApplyPadStyle(TPad* p) {
  if (!p) return;
  p->SetGridx();
  p->SetGridy();
  p->SetLeftMargin(0.14);
  p->SetRightMargin(0.05);
  p->SetTopMargin(0.10);
  p->SetBottomMargin(0.14);
}

static bool GetFullBinRef(const std::vector<Row>& rows,
                          const TString& setting,
                          const TString& var,
                          double& yRef,
                          double& eyRef) {
  yRef = std::numeric_limits<double>::quiet_NaN();
  eyRef = std::numeric_limits<double>::quiet_NaN();

  for (const auto& r : rows) {
    if (r.setting != setting) continue;
    if (r.var != var) continue;
    if (r.dp_idx != 0) continue;
    if (r.n_valid <= 0) continue;
    if (!IsFinite(r.weighted_mean_dmu_MeV)) continue;

    yRef = r.weighted_mean_dmu_MeV;
    eyRef = (IsFinite(r.weighted_mean_err_dmu_MeV) ? r.weighted_mean_err_dmu_MeV : 0.0);
    return true;
  }
  return false;
}

static void CollectNarrowBinPoints(const std::vector<Row>& rows,
                                   const TString& setting,
                                   const TString& var,
                                   std::vector<double>& xs,
                                   std::vector<double>& ys,
                                   std::vector<double>& exs,
                                   std::vector<double>& eys) {
  xs.clear();
  ys.clear();
  exs.clear();
  eys.clear();

  for (const auto& r : rows) {
    if (r.setting != setting) continue;
    if (r.var != var) continue;
    if (r.dp_idx <= 0) continue;
    if (r.n_valid <= 0) continue;
    if (!IsFinite(r.weighted_mean_dmu_MeV)) continue;

    double x = 0.5 * (r.dp_lo + r.dp_hi);
    xs.push_back(x);
    ys.push_back(r.weighted_mean_dmu_MeV);
    exs.push_back(0.0);
    eys.push_back(IsFinite(r.weighted_mean_err_dmu_MeV) ? r.weighted_mean_err_dmu_MeV : 0.0);
  }

  std::vector<size_t> order(xs.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](size_t a, size_t b) { return xs[a] < xs[b]; });

  std::vector<double> xs2, ys2, exs2, eys2;
  for (size_t idx : order) {
    xs2.push_back(xs[idx]);
    ys2.push_back(ys[idx]);
    exs2.push_back(exs[idx]);
    eys2.push_back(eys[idx]);
  }
  xs.swap(xs2);
  ys.swap(ys2);
  exs.swap(exs2);
  eys.swap(eys2);
}

static bool ComputeRanges(const std::vector<double>& xs,
                          const std::vector<double>& ys,
                          const std::vector<double>& eys,
                          double refY, double refEy,
                          double& xMin, double& xMax,
                          double& yMin, double& yMax) {
  bool found = false;
  xMin =  1e99;
  xMax = -1e99;
  yMin =  1e99;
  yMax = -1e99;

  for (size_t i = 0; i < xs.size(); ++i) {
    if (!IsFinite(xs[i]) || !IsFinite(ys[i])) continue;
    double ey = (i < eys.size() && IsFinite(eys[i])) ? eys[i] : 0.0;
    xMin = std::min(xMin, xs[i]);
    xMax = std::max(xMax, xs[i]);
    yMin = std::min(yMin, ys[i] - ey);
    yMax = std::max(yMax, ys[i] + ey);
    found = true;
  }

  if (IsFinite(refY)) {
    yMin = std::min(yMin, refY - (IsFinite(refEy) ? refEy : 0.0));
    yMax = std::max(yMax, refY + (IsFinite(refEy) ? refEy : 0.0));
    found = true;
  }

  if (!found) return false;

  if (!(xMax > xMin)) {
    xMin -= 0.5;
    xMax += 0.5;
  }

  double xPad = 0.15 * std::max(1.0, xMax - xMin);
  double yPad = 0.20 * std::max(1.0, yMax - yMin);

  xMin -= xPad;
  xMax += xPad;
  yMin -= yPad;
  yMax += yPad;
  return true;
}

static void DrawReferenceBandAndLine(double xLo, double xHi,
                                     double yRef, double eyRef) {
  if (!IsFinite(yRef)) return;

  if (IsFinite(eyRef) && eyRef > 0.0) {
    TBox* band = new TBox(xLo, yRef - eyRef, xHi, yRef + eyRef);
    band->SetFillColorAlpha(kBlue, 0.18);
    band->SetLineColorAlpha(kBlue, 0.0);
    band->Draw("same");
  }

  TLine* ln = new TLine(xLo, yRef, xHi, yRef);
  ln->SetLineColor(kBlue);
  ln->SetLineStyle(2);
  ln->SetLineWidth(2);
  ln->Draw("same");
}

static void DrawOneVariablePad(TPad* pad,
                               const std::vector<Row>& rows,
                               const TString& setting,
                               const TString& var,
                               bool drawLegend = false) {
  if (!pad) return;
  pad->cd();

  std::vector<double> xs, ys, exs, eys;
  CollectNarrowBinPoints(rows, setting, var, xs, ys, exs, eys);

  double refY = std::numeric_limits<double>::quiet_NaN();
  double refEy = std::numeric_limits<double>::quiet_NaN();
  GetFullBinRef(rows, setting, var, refY, refEy);

  double xMin = 0.0, xMax = 1.0, yMin = -1.0, yMax = 1.0;
  if (!ComputeRanges(xs, ys, eys, refY, refEy, xMin, xMax, yMin, yMax)) {
    xMin = -2.5;
    xMax =  3.5;
    yMin = -1.0;
    yMax =  1.0;
  }

  TH1D* frame = new TH1D(Form("frame_%s_%s", setting.Data(), var.Data()), "", 100, xMin, xMax);
  frame->SetMinimum(yMin);
  frame->SetMaximum(yMax);
  frame->GetXaxis()->SetTitle("#delta bin center (%)");
  frame->GetYaxis()->SetTitle("Weighted mean #Delta#mu (MeV)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->Draw();

  DrawReferenceBandAndLine(xMin, xMax, refY, refEy);

  if (!xs.empty()) {
    TGraphErrors* g = new TGraphErrors((int)xs.size());
    for (int i = 0; i < (int)xs.size(); ++i) {
      g->SetPoint(i, xs[i], ys[i]);
      g->SetPointError(i, exs[i], eys[i]);
    }
    g->SetLineColor(kBlack);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetLineWidth(2);
    g->Draw("P SAME");
  }

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.060);
  lat.DrawLatex(0.16, 0.90, PrettyVar(var));

  if (drawLegend) {
    TLegend* leg = new TLegend(0.58, 0.76, 0.91, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);

    TGraphErrors* dummyPoint = new TGraphErrors(1);
    dummyPoint->SetPoint(0, 0.0, 0.0);
    dummyPoint->SetPointError(0, 0.0, 0.0);
    dummyPoint->SetLineColor(kBlack);
    dummyPoint->SetMarkerColor(kBlack);
    dummyPoint->SetMarkerStyle(20);
    dummyPoint->SetMarkerSize(1.2);
    dummyPoint->SetLineWidth(2);
    leg->AddEntry(dummyPoint, "Narrow bins", "p");

    TLine* dummyLine = new TLine(0.0, 0.0, 1.0, 0.0);
    dummyLine->SetLineColor(kBlue);
    dummyLine->SetLineStyle(2);
    dummyLine->SetLineWidth(2);
    leg->AddEntry(dummyLine, "Full bin", "l");

    leg->Draw();
  }
}

static void MakeOneSettingCanvas(const std::vector<Row>& rows,
                                 const TString& setting,
                                 const TString& outPng) {
  TCanvas* c = new TCanvas(Form("c_%s", setting.Data()),
                           Form("Offset summary Setting %s", setting.Data()),
                           1200, 900);
  c->Divide(2, 2);

  const TString vars[4] = {"W", "Em", "Pmz", "Pmy"};
  for (int i = 0; i < 4; ++i) {
    TPad* p = (TPad*) c->cd(i + 1);
    ApplyPadStyle(p);
    DrawOneVariablePad(p, rows, setting, vars[i], i == 0);
  }

  c->SaveAs(outPng);
  delete c;
}

void PlotOffsetsSummaryBySetting(
    const char* inCsvC = "results/tables/OffsetsSummaryBySetting.csv") {

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);

  TString inCsv(inCsvC);

  std::ifstream fin(inCsv.Data());
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open input CSV: " << inCsv << "\n";
    return;
  }

  TString outDir = "results/PNGs/PlotOffsetsSummaryBySetting";
  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty CSV: " << inCsv << "\n";
    return;
  }

  std::vector<TString> header = SplitCSV(line);
  std::map<std::string,int> col;
  for (int i = 0; i < (int)header.size(); ++i) col[header[i].Data()] = i;

  auto need = [&](const char* name) {
    if (!col.count(name)) {
      std::cerr << "[ERROR] Missing required column: " << name << "\n";
      return false;
    }
    return true;
  };

  if (!need("setting") ||
      !need("dp_idx") ||
      !need("dp_label") ||
      !need("dp_lo") ||
      !need("dp_hi") ||
      !need("var") ||
      !need("n_valid") ||
      !need("weighted_mean_dmu_MeV") ||
      !need("weighted_mean_err_dmu_MeV")) {
    return;
  }

  std::vector<Row> rows;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    std::vector<TString> tok = SplitCSV(line);
    if ((int)tok.size() < (int)header.size()) continue;

    Row r;
    r.setting = tok[col["setting"]];
    r.dp_idx = tok[col["dp_idx"]].Atoi();
    r.dp_label = tok[col["dp_label"]];
    r.dp_lo = tok[col["dp_lo"]].Atof();
    r.dp_hi = tok[col["dp_hi"]].Atof();
    r.var = tok[col["var"]];
    r.n_valid = tok[col["n_valid"]].Atoi();
    r.weighted_mean_dmu_MeV = tok[col["weighted_mean_dmu_MeV"]].Atof();
    r.weighted_mean_err_dmu_MeV = tok[col["weighted_mean_err_dmu_MeV"]].Atof();
    rows.push_back(r);
  }
  fin.close();

  MakeOneSettingCanvas(rows, "A", outDir + "/offset_summary_settingA.png");
  MakeOneSettingCanvas(rows, "B", outDir + "/offset_summary_settingB.png");

  std::cout << "[INFO] Wrote PNGs to: " << outDir << "\n";
}
