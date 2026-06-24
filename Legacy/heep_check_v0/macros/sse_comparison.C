// File: heep_check_v0/macros/sse_comparison.C
//
// Purpose:
//   Compare SSE from:
//     results/tables/heepcheck_inversion_dpeOnly.csv
//     results/tables/heepcheck_inversion_dpeDpp.csv
//
// Output:
//   results/PNGs/sse_comparison/dpeOnly_vs_dpeDpp.png
//
// Canvas layout (4 pads):
//   1) SSE vs run for dp_idx=0 (full bin), log-y
//   2) SSE(dpeDpp)/SSE(dpeOnly) vs run for dp_idx=0
//   3) SSE ratio vs dp_idx for all matched points
//   4) Scatter: SSE_dpeOnly vs SSE_dpeDpp with y=x line, log-log
//
// Suggested usage from heep_check_v0/:
//   root -l -b -q 'macros/sse_comparison.C()'

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TROOT.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct Row {
  int    run = -1;
  int    dp_idx = -999;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  double sse = 0.0;
};

std::string Trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

std::vector<std::string> SplitCSVSimple(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(Trim(item));
  return out;
}

bool ToInt(const std::string& s, int& val) {
  if (s.empty()) return false;
  char* end = nullptr;
  errno = 0;
  long v = std::strtol(s.c_str(), &end, 10);
  if (errno != 0 || !end || *end != '\0') return false;
  val = static_cast<int>(v);
  return true;
}

bool ToDouble(const std::string& s, double& val) {
  if (s.empty()) return false;
  char* end = nullptr;
  errno = 0;
  val = std::strtod(s.c_str(), &end);
  if (errno != 0 || !end || *end != '\0') return false;
  return true;
}

std::map<std::string, int> HeaderIndex(const std::vector<std::string>& hdr) {
  std::map<std::string, int> idx;
  for (int i = 0; i < (int)hdr.size(); ++i) idx[Trim(hdr[i])] = i;
  return idx;
}

bool EnsureDir(const std::string& dir) {
  try {
    std::filesystem::path p(dir);
    if (!std::filesystem::exists(p)) std::filesystem::create_directories(p);
  } catch (...) {
    return false;
  }
  return true;
}

bool LoadSSECSV(const std::string& path, std::map<std::pair<int,int>, Row>& out) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty CSV: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);

  for (const auto& need : {"run","dp_idx","dp_lo","dp_hi","sse"}) {
    if (!idx.count(need)) {
      std::cerr << "[ERROR] Missing required column '" << need << "' in " << path << "\n";
      return false;
    }
  }

  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto tok = SplitCSVSimple(line);
    if ((int)tok.size() < (int)hdr.size()) continue;

    Row r;
    if (!ToInt(tok[idx.at("run")], r.run)) continue;
    if (!ToInt(tok[idx.at("dp_idx")], r.dp_idx)) continue;
    if (!ToDouble(tok[idx.at("dp_lo")], r.dp_lo)) continue;
    if (!ToDouble(tok[idx.at("dp_hi")], r.dp_hi)) continue;
    if (!ToDouble(tok[idx.at("sse")], r.sse)) continue;

    out[{r.run, r.dp_idx}] = r;
  }

  return true;
}

struct MatchedPoint {
  int    run = -1;
  int    dp_idx = -999;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  double sse_dpeOnly = 0.0;
  double sse_dpeDpp  = 0.0;
};

std::vector<MatchedPoint> BuildMatched(
    const std::map<std::pair<int,int>, Row>& dpeOnly,
    const std::map<std::pair<int,int>, Row>& dpeDpp)
{
  std::vector<MatchedPoint> out;

  for (const auto& kv : dpeOnly) {
    auto it = dpeDpp.find(kv.first);
    if (it == dpeDpp.end()) {
      std::cerr << "[WARN] Missing match in dpeDpp for run=" << kv.first.first
                << ", dp_idx=" << kv.first.second << "\n";
      continue;
    }

    MatchedPoint p;
    p.run = kv.second.run;
    p.dp_idx = kv.second.dp_idx;
    p.dp_lo = kv.second.dp_lo;
    p.dp_hi = kv.second.dp_hi;
    p.sse_dpeOnly = kv.second.sse;
    p.sse_dpeDpp  = it->second.sse;
    out.push_back(p);
  }

  std::sort(out.begin(), out.end(), [](const MatchedPoint& a, const MatchedPoint& b){
    if (a.run != b.run) return a.run < b.run;
    return a.dp_idx < b.dp_idx;
  });

  return out;
}

void SetPadStyle(bool logx=false, bool logy=false) {
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.08);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.14);
  if (logx) gPad->SetLogx(1);
  if (logy) gPad->SetLogy(1);
}

void StyleGraph(TGraph* g, int markerStyle, double markerSize) {
  if (!g) return;
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(markerSize);
  g->SetLineWidth(0);
}

void SetAxisTitles(TGraph* g, const std::string& xtitle, const std::string& ytitle) {
  if (!g) return;
  g->GetXaxis()->SetTitle(xtitle.c_str());
  g->GetYaxis()->SetTitle(ytitle.c_str());
  g->GetXaxis()->SetTitleSize(0.050);
  g->GetYaxis()->SetTitleSize(0.050);
  g->GetXaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetTitleOffset(1.20);
}

void DrawPadTitle(const std::string& text) {
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.052);
  lat.DrawLatex(0.14, 0.93, text.c_str());
}

double MinPositive(const std::vector<double>& v) {
  double m = 0.0;
  for (double x : v) {
    if (x > 0.0 && (m == 0.0 || x < m)) m = x;
  }
  return m;
}

} // namespace

void sse_comparison(
    const char* dpeOnly_csv = "results/tables/heepcheck_inversion_dpeOnly.csv",
    const char* dpeDpp_csv  = "results/tables/heepcheck_inversion_dpeDpp.csv",
    const char* out_png     = "results/PNGs/sse_comparison/dpeOnly_vs_dpeDpp.png")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::string f1 = dpeOnly_csv ? dpeOnly_csv : "results/tables/heepcheck_inversion_dpeOnly.csv";
  const std::string f2 = dpeDpp_csv  ? dpeDpp_csv  : "results/tables/heepcheck_inversion_dpeDpp.csv";
  const std::string outpng = out_png ? out_png : "results/PNGs/sse_comparison/dpeOnly_vs_dpeDpp.png";

  std::filesystem::path outpath(outpng);
  if (!EnsureDir(outpath.parent_path().string())) {
    std::cerr << "[ERROR] Could not create output directory for: " << outpng << "\n";
    return;
  }

  std::map<std::pair<int,int>, Row> mapDpeOnly, mapDpeDpp;
  if (!LoadSSECSV(f1, mapDpeOnly)) return;
  if (!LoadSSECSV(f2, mapDpeDpp)) return;

  auto matched = BuildMatched(mapDpeOnly, mapDpeDpp);
  if (matched.empty()) {
    std::cerr << "[ERROR] No matched (run, dp_idx) rows found.\n";
    return;
  }

  std::cout << "[INFO] Matched " << matched.size() << " (run, dp_idx) points.\n";

  // ---- Pad 1 and 2 data: dp_idx = 0 only
  std::vector<double> x_run, y_only_full, y_dpp_full, y_ratio_full;
  for (const auto& p : matched) {
    if (p.dp_idx != 0) continue;
    x_run.push_back(p.run);
    y_only_full.push_back(p.sse_dpeOnly);
    y_dpp_full.push_back(p.sse_dpeDpp);
    y_ratio_full.push_back((p.sse_dpeOnly != 0.0) ? p.sse_dpeDpp / p.sse_dpeOnly : 0.0);
  }

  // ---- Pad 3 data: all matched points ratio vs dp_idx
  std::vector<double> x_dpidx_all, y_ratio_all;
  for (const auto& p : matched) {
    x_dpidx_all.push_back(p.dp_idx);
    y_ratio_all.push_back((p.sse_dpeOnly != 0.0) ? p.sse_dpeDpp / p.sse_dpeOnly : 0.0);
  }
  const int nMatchedAll = (int)x_dpidx_all.size();

  // ---- Pad 4 data: scatter of SSE_dpeOnly vs SSE_dpeDpp
  std::vector<double> x_scatter, y_scatter;
  for (const auto& p : matched) {
    if (p.sse_dpeOnly > 0.0 && p.sse_dpeDpp > 0.0) {
      x_scatter.push_back(p.sse_dpeOnly);
      y_scatter.push_back(p.sse_dpeDpp);
    }
  }

  TCanvas* c = new TCanvas("c_sse_comparison", "SSE comparison", 1500, 1100);
  c->Divide(2,2);

  // --------------------------------------------------------------------------
  // Pad 1: full-bin raw SSE vs run
  // --------------------------------------------------------------------------
  c->cd(1);
  SetPadStyle(false, true);

  if (!x_run.empty()) {
    TGraph* gOnly = new TGraph((int)x_run.size(), x_run.data(), y_only_full.data());
    TGraph* gDpp  = new TGraph((int)x_run.size(), x_run.data(), y_dpp_full.data());

    StyleGraph(gOnly, 20, 1.1);
    StyleGraph(gDpp,  21, 1.1);

    // Build combined range from both series
    std::vector<double> allY;
    allY.reserve(y_only_full.size() + y_dpp_full.size());
    for (double v : y_only_full) if (v > 0.0) allY.push_back(v);
    for (double v : y_dpp_full)  if (v > 0.0) allY.push_back(v);

    double ymin = *std::min_element(allY.begin(), allY.end());
    double ymax = *std::max_element(allY.begin(), allY.end());

    ymin *= 0.8;
    ymax *= 1.2;

    gOnly->Draw("AP");
    SetAxisTitles(gOnly, "Run", "SSE");
    DrawPadTitle("SSE vs run (dp_idx = 0)");
    gOnly->SetMinimum(ymin);
    gOnly->SetMaximum(ymax);

    gDpp->Draw("P same");

    TLegend* leg = new TLegend(0.16, 0.18, 0.42, 0.32);
    leg->AddEntry(gOnly, "dpeOnly", "p");
    leg->AddEntry(gDpp,  "dpeDpp",  "p");
    leg->Draw();
  }

  // --------------------------------------------------------------------------
  // Pad 2: full-bin ratio vs run
  // --------------------------------------------------------------------------
  c->cd(2);
  SetPadStyle(false, false);

  if (!x_run.empty()) {
    TGraph* gRatio = new TGraph((int)x_run.size(), x_run.data(), y_ratio_full.data());
    StyleGraph(gRatio, 20, 1.1);
    gRatio->Draw("AP");
    SetAxisTitles(gRatio, "Run", "SSE(dpeDpp) / SSE(dpeOnly)");
    DrawPadTitle("SSE ratio vs run (dp_idx = 0)");

    double xmin = *std::min_element(x_run.begin(), x_run.end());
    double xmax = *std::max_element(x_run.begin(), x_run.end());
    TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
  }

  // --------------------------------------------------------------------------
  // Pad 3: ratio vs dp_idx for all matched points
  // --------------------------------------------------------------------------
  c->cd(3);
  SetPadStyle(false, false);

  if (!x_dpidx_all.empty()) {
    TGraph* gRatioAll = new TGraph((int)x_dpidx_all.size(), x_dpidx_all.data(), y_ratio_all.data());
    StyleGraph(gRatioAll, 20, 1.0);
    gRatioAll->Draw("AP");
    SetAxisTitles(gRatioAll, "dp_idx", "SSE(dpeDpp) / SSE(dpeOnly)");
    DrawPadTitle(Form("SSE ratio vs dp_idx (N=%d)", nMatchedAll));

    TLine* line = new TLine(0.0, 1.0, 5.0, 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
  }

  // --------------------------------------------------------------------------
  // Pad 4: scatter raw SSE comparison, log-log
  // --------------------------------------------------------------------------
  c->cd(4);
  SetPadStyle(true, true);

  if (!x_scatter.empty()) {
    TGraph* gSc = new TGraph((int)x_scatter.size(), x_scatter.data(), y_scatter.data());
    StyleGraph(gSc, 20, 1.0);
    gSc->Draw("AP");
    SetAxisTitles(gSc, "SSE(dpeOnly)", "SSE(dpeDpp)");
    DrawPadTitle("SSE scatter: dpeOnly vs dpeDpp");

    double xmin = *std::min_element(x_scatter.begin(), x_scatter.end());
    double xmax = *std::max_element(x_scatter.begin(), x_scatter.end());
    double ymin = *std::min_element(y_scatter.begin(), y_scatter.end());
    double ymax = *std::max_element(y_scatter.begin(), y_scatter.end());

    // log-safe margins
    double xlo = xmin / 1.5;
    double xhi = xmax * 1.5;
    double ylo = ymin / 1.5;
    double yhi = ymax * 1.5;

    gSc->GetXaxis()->SetLimits(xlo, xhi);
    gSc->SetMinimum(ylo);
    gSc->SetMaximum(yhi);

    // diagonal only within visible overlap
    double diagMin = std::max(xlo, ylo);
    double diagMax = std::min(xhi, yhi);

    TLine* diag = new TLine(diagMin, diagMin, diagMax, diagMax);
    diag->SetLineStyle(2);
    diag->SetLineWidth(2);
    diag->Draw("same");
  }

  c->SaveAs(outpng.c_str());
  std::cout << "[INFO] Wrote " << outpng << "\n";
  delete c;
}
