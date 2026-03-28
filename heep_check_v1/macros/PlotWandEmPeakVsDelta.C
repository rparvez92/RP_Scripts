// PlotWandEmPeakVsDelta.C
//
// Purpose:
//   Read WandEmPeakVsDelta.csv and produce:
//
//   1) Per-run plots:
//      heep_check_v1/results/PNGs/PlotWandEmPeakVsDelta/settingA_run23839.png
//      ...
//      Each canvas has two pads:
//        upper: W_peak vs delta
//        lower: Em_peak vs delta
//
//   2) Setting-mean plots:
//      heep_check_v1/results/PNGs/PlotWandEmPeakVsDelta/settingA_mean.png
//      heep_check_v1/results/PNGs/PlotWandEmPeakVsDelta/settingB_mean.png
//      Each canvas has two pads:
//        upper: weighted-mean W_peak vs delta
//        lower: weighted-mean Em_peak vs delta
//
// Plot rules:
//   - Narrow bins only are plotted as points at dp centers
//   - Full bin (dp_idx=0) is NOT plotted as a point
//   - Full bin is shown as a blue dashed horizontal line
//   - Full-bin uncertainty is shown as a transparent band
//   - Units are MeV on y-axis
//
// Example:
//   root -l -b -q 'macros/PlotWandEmPeakVsDelta.C+("results/tables/WandEmPeakVsDelta.csv","results/PNGs/PlotWandEmPeakVsDelta",true,true)'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TColor.h"

// -----------------------------------------------------------------------------
// CSV row structure
// -----------------------------------------------------------------------------
struct Row {
  int run = -1;
  TString setting = "";
  double ebeam_MeV = std::numeric_limits<double>::quiet_NaN();

  int dp_idx = -1;
  TString dp_label = "";
  double dp_lo = std::numeric_limits<double>::quiet_NaN();
  double dp_hi = std::numeric_limits<double>::quiet_NaN();
  double dp_center = std::numeric_limits<double>::quiet_NaN();

  TString var = ""; // W or Em
  double entries = std::numeric_limits<double>::quiet_NaN();

  double fit_lo_GeV = std::numeric_limits<double>::quiet_NaN();
  double fit_hi_GeV = std::numeric_limits<double>::quiet_NaN();
  double mu_seed_GeV = std::numeric_limits<double>::quiet_NaN();
  double sig_seed_GeV = std::numeric_limits<double>::quiet_NaN();

  double peak_GeV = std::numeric_limits<double>::quiet_NaN();
  double peak_err_GeV = std::numeric_limits<double>::quiet_NaN();
  double sigma_GeV = std::numeric_limits<double>::quiet_NaN();
  double sigma_err_GeV = std::numeric_limits<double>::quiet_NaN();

  double peak_MeV = std::numeric_limits<double>::quiet_NaN();
  double peak_err_MeV = std::numeric_limits<double>::quiet_NaN();
  double sigma_MeV = std::numeric_limits<double>::quiet_NaN();
  double sigma_err_MeV = std::numeric_limits<double>::quiet_NaN();

  int fit_status = -999;
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = -999;
  double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
  int fit_valid = 0;
};

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static std::vector<TString> SplitCsvLine(const TString& s) {
  std::vector<TString> out;
  std::string line = s.Data();
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(TString(item));
  return out;
}

static bool IsFinite(double x) {
  return std::isfinite(x);
}

static double ParseDouble(const TString& s) {
  TString t = s;
  t = t.Strip(TString::kBoth);
  if (t.Length() == 0) return std::numeric_limits<double>::quiet_NaN();

  TString low = t; low.ToLower();
  if (low == "nan") return std::numeric_limits<double>::quiet_NaN();
  if (low == "-nan") return std::numeric_limits<double>::quiet_NaN();

  return t.Atof();
}

static int ParseInt(const TString& s) {
  TString t = s;
  t = t.Strip(TString::kBoth);
  if (t.Length() == 0) return -999;
  return t.Atoi();
}

static bool ReadCsv(const TString& csvPath, std::vector<Row>& rows) {
  std::ifstream fin(csvPath.Data());
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open CSV: " << csvPath << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty CSV: " << csvPath << "\n";
    return false;
  } // skip header

  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    std::vector<TString> f = SplitCsvLine(line.c_str());
    if (f.size() < 27) {
      std::cerr << "[WARN] Skipping malformed line with " << f.size() << " columns\n";
      continue;
    }

    Row r;
    int i = 0;
    r.run            = ParseInt   (f[i++]);
    r.setting        =            f[i++];
    r.ebeam_MeV      = ParseDouble(f[i++]);

    r.dp_idx         = ParseInt   (f[i++]);
    r.dp_label       =            f[i++];
    r.dp_lo          = ParseDouble(f[i++]);
    r.dp_hi          = ParseDouble(f[i++]);
    r.dp_center      = ParseDouble(f[i++]);

    r.var            =            f[i++];
    r.entries        = ParseDouble(f[i++]);

    r.fit_lo_GeV     = ParseDouble(f[i++]);
    r.fit_hi_GeV     = ParseDouble(f[i++]);
    r.mu_seed_GeV    = ParseDouble(f[i++]);
    r.sig_seed_GeV   = ParseDouble(f[i++]);

    r.peak_GeV       = ParseDouble(f[i++]);
    r.peak_err_GeV   = ParseDouble(f[i++]);
    r.sigma_GeV      = ParseDouble(f[i++]);
    r.sigma_err_GeV  = ParseDouble(f[i++]);

    r.peak_MeV       = ParseDouble(f[i++]);
    r.peak_err_MeV   = ParseDouble(f[i++]);
    r.sigma_MeV      = ParseDouble(f[i++]);
    r.sigma_err_MeV  = ParseDouble(f[i++]);

    r.fit_status     = ParseInt   (f[i++]);
    r.chi2           = ParseDouble(f[i++]);
    r.ndf            = ParseInt   (f[i++]);
    r.chi2_ndf       = ParseDouble(f[i++]);
    r.fit_valid      = ParseInt   (f[i++]);

    rows.push_back(r);
  }

  fin.close();
  return true;
}

static std::vector<int> GetAllRunsSorted(const std::vector<Row>& rows) {
  std::vector<int> runs;
  for (const auto& r : rows) {
    if (std::find(runs.begin(), runs.end(), r.run) == runs.end())
      runs.push_back(r.run);
  }
  std::sort(runs.begin(), runs.end());
  return runs;
}

static TString SettingFromRun(int run) {
  if (run >= 23839 && run <= 23848 && run != 23843) return "A";
  if (run >= 23849 && run <= 23851) return "B";
  return "UNKNOWN";
}

static bool FindRow(const std::vector<Row>& rows,
                    int run, const TString& var, int dp_idx,
                    Row& out) {
  for (const auto& r : rows) {
    if (r.run == run && r.var == var && r.dp_idx == dp_idx) {
      out = r;
      return true;
    }
  }
  return false;
}

static std::vector<Row> SelectRowsBySettingVarDp(const std::vector<Row>& rows,
                                                 const TString& setting,
                                                 const TString& var,
                                                 int dp_idx) {
  std::vector<Row> out;
  for (const auto& r : rows) {
    if (r.setting == setting && r.var == var && r.dp_idx == dp_idx)
      out.push_back(r);
  }
  return out;
}

struct MeanPoint {
  double x = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
  double ey = std::numeric_limits<double>::quiet_NaN();
  int nUsed = 0;
};

static MeanPoint WeightedMeanPeak(const std::vector<Row>& rowsForOneBin) {
  MeanPoint m;
  if (rowsForOneBin.empty()) return m;

  double xref = rowsForOneBin.front().dp_center;
  m.x = xref;

  double sw = 0.0;
  double swy = 0.0;
  int n = 0;

  for (const auto& r : rowsForOneBin) {
    if (r.fit_valid != 1) continue;
    if (!IsFinite(r.peak_MeV) || !IsFinite(r.peak_err_MeV)) continue;
    if (!(r.peak_err_MeV > 0.0)) continue;

    double w = 1.0 / (r.peak_err_MeV * r.peak_err_MeV);
    sw  += w;
    swy += w * r.peak_MeV;
    n++;
  }

  if (sw > 0.0) {
    m.y = swy / sw;
    m.ey = std::sqrt(1.0 / sw);
    m.nUsed = n;
  }

  return m;
}

static void EnsureOutDir(const TString& outDir) {
  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);
}

static void ApplyPadStyle(TPad* p, bool topPad) {
  if (!p) return;
  p->SetGridx();
  p->SetGridy();
  p->SetLeftMargin(0.12);
  p->SetRightMargin(0.05);
  p->SetTopMargin(topPad ? 0.08 : 0.03);
  p->SetBottomMargin(topPad ? 0.08 : 0.16);
}

static void ComputeYRange(const std::vector<double>& ys,
                          const std::vector<double>& eys,
                          double refY, double refEy,
                          double& yMin, double& yMax) {
  double ymin =  1e99;
  double ymax = -1e99;
  bool found = false;

  for (size_t i = 0; i < ys.size(); i++) {
    if (!IsFinite(ys[i])) continue;
    double ey = (i < eys.size() && IsFinite(eys[i])) ? eys[i] : 0.0;
    ymin = std::min(ymin, ys[i] - ey);
    ymax = std::max(ymax, ys[i] + ey);
    found = true;
  }

  if (IsFinite(refY)) {
    double ey = (IsFinite(refEy) ? refEy : 0.0);
    ymin = std::min(ymin, refY - ey);
    ymax = std::max(ymax, refY + ey);
    found = true;
  }

  if (!found) {
    yMin = -1.0;
    yMax =  1.0;
    return;
  }

  double span = ymax - ymin;
  if (!(span > 0.0)) {
    span = std::max(1.0, std::fabs(ymax));
  }
  double pad = 0.20 * span;
  yMin = ymin - pad;
  yMax = ymax + pad;
}

static void DrawReferenceBandAndLine(double xLo, double xHi,
                                     double yRef, double eyRef,
                                     int bandColor = kBlue, int lineColor = kBlue) {
  if (!IsFinite(yRef)) return;

  if (IsFinite(eyRef) && eyRef > 0.0) {
    TBox* band = new TBox(xLo, yRef - eyRef, xHi, yRef + eyRef);
    band->SetFillColorAlpha(bandColor, 0.18);
    band->SetLineColorAlpha(bandColor, 0.0);
    band->Draw("same");
  }

  TLine* ln = new TLine(xLo, yRef, xHi, yRef);
  ln->SetLineColor(lineColor);
  ln->SetLineStyle(2);
  ln->SetLineWidth(2);
  ln->Draw("same");
}

static void StyleGraph(TGraphErrors* g) {
  if (!g) return;
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
}

static void DrawOnePad(TPad* pad,
                       const std::vector<double>& xs,
                       const std::vector<double>& ys,
                       const std::vector<double>& exs,
                       const std::vector<double>& eys,
                       const TString& yTitle,
                       const TString& topText,
                       double refY, double refEy,
                       bool drawLegend = true,
                       bool showXTitle = true) {
  if (!pad) return;
  pad->cd();

  const double xMin = -2.1;
  const double xMax =  3.1;

  double yMin = -1.0;
  double yMax =  1.0;
  ComputeYRange(ys, eys, refY, refEy, yMin, yMax);
  if (!(yMax > yMin)) {
    yMin = -1.0;
    yMax =  1.0;
  }

  TH1* frame = pad->DrawFrame(xMin, yMin, xMax, yMax);
  if (showXTitle) {
    frame->GetXaxis()->SetTitle("#delta bin center (%)");
    frame->GetXaxis()->SetTitleSize(0.055);
  } else {
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetTitleSize(0.0);
  }
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.05);
  frame->GetXaxis()->SetLimits(xMin, xMax);

  DrawReferenceBandAndLine(xMin, xMax, refY, refEy, kBlue, kBlue);

  int n = (int)xs.size();
  TGraphErrors* g = new TGraphErrors(n);
  int ip = 0;
  for (int i = 0; i < n; i++) {
    if (!IsFinite(xs[i]) || !IsFinite(ys[i])) continue;
    double ex = (i < (int)exs.size() && IsFinite(exs[i])) ? exs[i] : 0.0;
    double ey = (i < (int)eys.size() && IsFinite(eys[i])) ? eys[i] : 0.0;
    g->SetPoint(ip, xs[i], ys[i]);
    g->SetPointError(ip, ex, ey);
    ip++;
  }
  g->Set(ip);
  StyleGraph(g);
  g->Draw("P SAME");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.14, 0.93, topText);

  if (drawLegend) {
    TLegend* leg = new TLegend(0.64, 0.78, 0.93, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g, "Narrow bins", "pe");

    TLine* dummyLine = new TLine(0, 0, 1, 0);
    dummyLine->SetLineColor(kBlue);
    dummyLine->SetLineStyle(2);
    dummyLine->SetLineWidth(2);
    leg->AddEntry(dummyLine, "Full-bin mean", "l");

    leg->Draw();  
  }

  pad->Modified();
}

static void MakePerRunPlot(const std::vector<Row>& rows,
                           int run,
                           const TString& outDir,
                           bool verbose = true) {
  TString setting = SettingFromRun(run);

  std::vector<double> xW, yW, exW, eyW;
  std::vector<double> xE, yE, exE, eyE;

  double fullW = std::numeric_limits<double>::quiet_NaN();
  double fullWErr = std::numeric_limits<double>::quiet_NaN();
  double fullE = std::numeric_limits<double>::quiet_NaN();
  double fullEErr = std::numeric_limits<double>::quiet_NaN();

  // Full-bin refs
  {
    Row r;
    if (FindRow(rows, run, "W", 0, r) && r.fit_valid == 1) {
      fullW = r.peak_MeV;
      fullWErr = r.peak_err_MeV;
    }
    if (FindRow(rows, run, "Em", 0, r) && r.fit_valid == 1) {
      fullE = r.peak_MeV;
      fullEErr = r.peak_err_MeV;
    }
  }

  // Narrow bins 1..5
  for (int dp = 1; dp <= 5; dp++) {
    Row rW, rE;
    if (FindRow(rows, run, "W", dp, rW) && rW.fit_valid == 1 &&
        IsFinite(rW.dp_center) && IsFinite(rW.peak_MeV)) {
      xW.push_back(rW.dp_center);
      yW.push_back(rW.peak_MeV);
      exW.push_back(0.0);
      eyW.push_back(IsFinite(rW.peak_err_MeV) ? rW.peak_err_MeV : 0.0);
    }
    if (FindRow(rows, run, "Em", dp, rE) && rE.fit_valid == 1 &&
        IsFinite(rE.dp_center) && IsFinite(rE.peak_MeV)) {
      xE.push_back(rE.dp_center);
      yE.push_back(rE.peak_MeV);
      exE.push_back(0.0);
      eyE.push_back(IsFinite(rE.peak_err_MeV) ? rE.peak_err_MeV : 0.0);
    }
  }

  TCanvas* c = new TCanvas(Form("c_run_%d", run),
                           Form("setting%s_run%d", setting.Data(), run),
                           900, 900);
  c->Divide(1, 2);

  TPad* p1 = (TPad*) c->cd(1);
  ApplyPadStyle(p1, true);
  DrawOnePad(p1, xW, yW, exW, eyW,
             "W peak (MeV)",
             Form("Setting %s, Run %d", setting.Data(), run),
             fullW, fullWErr, true, false);

  TPad* p2 = (TPad*) c->cd(2);
  ApplyPadStyle(p2, false);
  DrawOnePad(p2, xE, yE, exE, eyE,
             "Em peak (MeV)",
             Form("Setting %s, Run %d", setting.Data(), run),
             fullE, fullEErr, false, true);

  TString outPng = outDir + Form("/setting%s_run%d.png", setting.Data(), run);
  c->SaveAs(outPng);

  if (verbose) std::cout << "[INFO] Wrote " << outPng << "\n";
  delete c;
}

static void MakeSettingMeanPlot(const std::vector<Row>& rows,
                                const TString& setting,
                                const TString& outDir,
                                bool verbose = true) {
  std::vector<double> xW, yW, exW, eyW;
  std::vector<double> xE, yE, exE, eyE;

  // Narrow-bin weighted means
  for (int dp = 1; dp <= 5; dp++) {
    auto rowsW = SelectRowsBySettingVarDp(rows, setting, "W", dp);
    auto rowsE = SelectRowsBySettingVarDp(rows, setting, "Em", dp);

    MeanPoint mW = WeightedMeanPeak(rowsW);
    MeanPoint mE = WeightedMeanPeak(rowsE);

    if (IsFinite(mW.x) && IsFinite(mW.y)) {
      xW.push_back(mW.x);
      yW.push_back(mW.y);
      exW.push_back(0.0);
      eyW.push_back(IsFinite(mW.ey) ? mW.ey : 0.0);
    }
    if (IsFinite(mE.x) && IsFinite(mE.y)) {
      xE.push_back(mE.x);
      yE.push_back(mE.y);
      exE.push_back(0.0);
      eyE.push_back(IsFinite(mE.ey) ? mE.ey : 0.0);
    }
  }

  // Full-bin weighted means
  MeanPoint fullW = WeightedMeanPeak(SelectRowsBySettingVarDp(rows, setting, "W", 0));
  MeanPoint fullE = WeightedMeanPeak(SelectRowsBySettingVarDp(rows, setting, "Em", 0));

  TCanvas* c = new TCanvas(Form("c_setting_%s_mean", setting.Data()),
                           Form("setting%s_mean", setting.Data()),
                           900, 900);
  c->Divide(1, 2);

  TPad* p1 = (TPad*) c->cd(1);
  ApplyPadStyle(p1, true);
  DrawOnePad(p1, xW, yW, exW, eyW,
             "Weighted mean W peak (MeV)",
             Form("Setting %s mean over runs", setting.Data()),
             fullW.y, fullW.ey, true, false);

  TPad* p2 = (TPad*) c->cd(2);
  ApplyPadStyle(p2, false);
  DrawOnePad(p2, xE, yE, exE, eyE,
             "Weighted mean Em peak (MeV)",
             Form("Setting %s mean over runs", setting.Data()),
             fullE.y, fullE.ey, false, true);

  // Add n_used text
  p1->cd();
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.037);
    tx.DrawLatex(0.14, 0.86,
      Form("Full-bin weighted mean built from %d run(s)", fullW.nUsed));
  }
  p2->cd();
  {
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.037);
    tx.DrawLatex(0.14, 0.86,
      Form("Full-bin weighted mean built from %d run(s)", fullE.nUsed));
  }

  TString outPng = outDir + Form("/setting%s_mean.png", setting.Data());
  c->SaveAs(outPng);

  if (verbose) std::cout << "[INFO] Wrote " << outPng << "\n";
  delete c;
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------
void PlotWandEmPeakVsDelta(const char* csvPathC =
                             "heep_check_v1/results/tables/WandEmPeakVsDelta.csv",
                           const char* outDirC =
                             "heep_check_v1/results/PNGs/PlotWandEmPeakVsDelta",
                           bool makePerRun = true,
                           bool makeSettingMean = true) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  TString csvPath(csvPathC);
  TString outDir(outDirC);

  EnsureOutDir(outDir);

  std::vector<Row> rows;
  if (!ReadCsv(csvPath, rows)) return;

  if (rows.empty()) {
    std::cerr << "[ERROR] No rows read from CSV.\n";
    return;
  }

  std::vector<int> runs = GetAllRunsSorted(rows);

  if (makePerRun) {
    for (int run : runs) {
      MakePerRunPlot(rows, run, outDir, true);
    }
  }

  if (makeSettingMean) {
    MakeSettingMeanPlot(rows, "A", outDir, true);
    MakeSettingMeanPlot(rows, "B", outDir, true);
  }

  std::cout << "[INFO] Done.\n";
}
