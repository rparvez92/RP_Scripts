// PlotLsFitOffsetsBySettingAndDpBin.C
//
// Read least-squares fit summary table from:
//   results/tables/LsFitOffsetsBySettingAndDpBin_summary.csv
//
// Make 4 PNGs:
//   dthe_vs_dpidx.png
//   dpe_vs_dpidx.png
//   dthp_vs_dpidx.png
//   dpp_vs_dpidx.png
//
// Each PNG has two pads:
//   top    -> Setting A
//   bottom -> Setting B
//
// X-axis uses numeric dp_idx.
//
// Output PNGs:
//   results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/*.png
//
// Usage:
//   root -l -b -q 'macros/PlotLsFitOffsetsBySettingAndDpBin.C("results/tables/LsFitOffsetsBySettingAndDpBin_summary.csv")'

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TPad.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct SummaryRow {
  std::string setting;
  int dp_idx = -999;
  std::string dp_label;

  int n_rows = 0;
  int n_runs = 0;

  double dthe = NAN, dthe_err = NAN;
  double dpe  = NAN, dpe_err  = NAN;
  double dthp = NAN, dthp_err = NAN;
  double dpp  = NAN, dpp_err  = NAN;

  double chi2 = NAN;
  int ndf = -999;
  double chi2_per_ndf = NAN;
};

std::string Trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inQuotes = false;

  for (char c : line) {
    if (c == '"') {
      inQuotes = !inQuotes;
    } else if (c == ',' && !inQuotes) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

double ToDouble(const std::string& s) {
  const std::string t = Trim(s);
  if (t.empty()) return NAN;
  return std::stod(t);
}

int ToInt(const std::string& s) {
  const std::string t = Trim(s);
  if (t.empty()) return -999;
  return std::stoi(t);
}

bool ReadSummaryCSV(const std::string& csvPath, std::vector<SummaryRow>& rows) {
  std::ifstream fin(csvPath);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open CSV: " << csvPath << std::endl;
    return false;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) {
    std::cerr << "[ERROR] Empty CSV: " << csvPath << std::endl;
    return false;
  }

  auto headers = SplitCSV(headerLine);
  std::map<std::string, int> col;
  for (int i = 0; i < (int)headers.size(); ++i) {
    col[Trim(headers[i])] = i;
  }

  const std::vector<std::string> need = {
    "setting", "dp_idx", "dp_label", "n_rows", "n_runs",
    "dthe", "dthe_err", "dpe", "dpe_err",
    "dthp", "dthp_err", "dpp", "dpp_err",
    "chi2", "ndf", "chi2_per_ndf"
  };

  for (const auto& k : need) {
    if (!col.count(k)) {
      std::cerr << "[ERROR] Missing required column: " << k << std::endl;
      return false;
    }
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    auto f = SplitCSV(line);
    if ((int)f.size() < (int)headers.size()) f.resize(headers.size());

    SummaryRow r;
    r.setting      = Trim(f[col["setting"]]);
    r.dp_idx       = ToInt(f[col["dp_idx"]]);
    r.dp_label     = Trim(f[col["dp_label"]]);
    r.n_rows       = ToInt(f[col["n_rows"]]);
    r.n_runs       = ToInt(f[col["n_runs"]]);

    r.dthe         = ToDouble(f[col["dthe"]]);
    r.dthe_err     = ToDouble(f[col["dthe_err"]]);
    r.dpe          = ToDouble(f[col["dpe"]]);
    r.dpe_err      = ToDouble(f[col["dpe_err"]]);
    r.dthp         = ToDouble(f[col["dthp"]]);
    r.dthp_err     = ToDouble(f[col["dthp_err"]]);
    r.dpp          = ToDouble(f[col["dpp"]]);
    r.dpp_err      = ToDouble(f[col["dpp_err"]]);

    r.chi2         = ToDouble(f[col["chi2"]]);
    r.ndf          = ToInt(f[col["ndf"]]);
    r.chi2_per_ndf = ToDouble(f[col["chi2_per_ndf"]]);

    rows.push_back(r);
  }

  return true;
}

bool SortByDpIdx(const SummaryRow& a, const SummaryRow& b) {
  return a.dp_idx < b.dp_idx;
}

double GetY(const SummaryRow& r, const std::string& parName) {
  if (parName == "dthe") return r.dthe;
  if (parName == "dpe")  return r.dpe;
  if (parName == "dthp") return r.dthp;
  if (parName == "dpp")  return r.dpp;
  return NAN;
}

double GetE(const SummaryRow& r, const std::string& parName) {
  if (parName == "dthe") return r.dthe_err;
  if (parName == "dpe")  return r.dpe_err;
  if (parName == "dthp") return r.dthp_err;
  if (parName == "dpp")  return r.dpp_err;
  return NAN;
}

std::vector<SummaryRow> FilterSetting(const std::vector<SummaryRow>& rows,
                                      const std::string& setting) {
  std::vector<SummaryRow> out;
  for (const auto& r : rows) {
    if (r.setting == setting) out.push_back(r);
  }
  std::sort(out.begin(), out.end(), SortByDpIdx);
  return out;
}

TGraphErrors* BuildGraph(const std::vector<SummaryRow>& rows,
                         const std::string& parName) {
  TGraphErrors* g = new TGraphErrors((int)rows.size());
  for (int i = 0; i < (int)rows.size(); ++i) {
    g->SetPoint(i, rows[i].dp_idx, GetY(rows[i], parName));
    g->SetPointError(i, 0.0, GetE(rows[i], parName));
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.2);
  g->SetLineWidth(2);
  return g;
}

void FindRanges(const std::vector<SummaryRow>& rows,
                const std::string& parName,
                double& xmin, double& xmax,
                double& ymin, double& ymax) {
  xmin =  1e99;
  xmax = -1e99;
  ymin =  1e99;
  ymax = -1e99;

  for (const auto& r : rows) {
    const double x = r.dp_idx;
    const double y = GetY(r, parName);
    const double e = GetE(r, parName);

    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y - e);
    ymax = std::max(ymax, y + e);
  }

  if (!std::isfinite(xmin) || !std::isfinite(xmax)) {
    xmin = -0.5;
    xmax =  5.5;
  }
  if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
    ymin = -1.0;
    ymax =  1.0;
  }

  if (xmin == xmax) {
    xmin -= 0.5;
    xmax += 0.5;
  }

  double yr = ymax - ymin;
  if (yr <= 0.0) yr = 1.0;
  double ypad = 0.18 * yr;
  ymin -= ypad;
  ymax += ypad;

  xmin -= 0.5;
  xmax += 0.5;
}

void DrawSettingPad(TPad* pad,
                    const std::vector<SummaryRow>& rows,
                    const std::string& setting,
                    const std::string& parName,
                    const std::string& yTitle,
                    bool drawXTitle) {
  pad->cd();
  pad->SetGrid();

  auto srows = FilterSetting(rows, setting);
  double xmin, xmax, ymin, ymax;
  FindRanges(srows, parName, xmin, xmax, ymin, ymax);

  TH1D* frame = new TH1D(("frame_" + setting + "_" + parName).c_str(), "",
                         100, xmin, xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetYaxis()->SetTitle(yTitle.c_str());
  frame->GetXaxis()->SetTitle(drawXTitle ? "dp_{idx}" : "");
  frame->GetXaxis()->SetNdivisions(510);
  frame->Draw();

  TLine* zero = new TLine(xmin, 0.0, xmax, 0.0);
  zero->SetLineStyle(2);
  zero->Draw("same");

  TGraphErrors* g = BuildGraph(srows, parName);
  g->Draw("P SAME");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.05);
  lat.DrawLatex(0.13, 0.88, ("Setting " + setting).c_str());
}

void DrawOnePlot(const std::vector<SummaryRow>& rows,
                 const std::string& parName,
                 const std::string& yTitle,
                 const std::string& unitText,
                 const std::string& outPng) {
  TCanvas* c = new TCanvas(("c_" + parName).c_str(), parName.c_str(), 1000, 900);
  c->Divide(1, 2, 0.0, 0.0);

  TPad* p1 = (TPad*)c->cd(1);
  p1->SetBottomMargin(0.04);
  p1->SetTopMargin(0.10);
  p1->SetLeftMargin(0.12);
  p1->SetRightMargin(0.05);

  TPad* p2 = (TPad*)c->cd(2);
  p2->SetTopMargin(0.04);
  p2->SetBottomMargin(0.12);
  p2->SetLeftMargin(0.12);
  p2->SetRightMargin(0.05);

  DrawSettingPad(p1, rows, "A", parName, yTitle, false);
  DrawSettingPad(p2, rows, "B", parName, yTitle, true);

  c->cd();
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.12, 0.975, ("Fitted " + parName + " vs dp_{idx}").c_str());
  lat.DrawLatex(0.70, 0.975, unitText.c_str());

  c->SaveAs(outPng.c_str());
  delete c;
}

} // namespace

void PlotLsFitOffsetsBySettingAndDpBin(
    const char* csvPath = "results/tables/LsFitOffsetsBySettingAndDpBin_summary.csv") {
  gStyle->SetOptStat(0);

  std::vector<SummaryRow> rows;
  if (!ReadSummaryCSV(csvPath, rows)) return;

  fs::create_directories("results/PNGs/PlotLsFitOffsetsBySettingAndDpBin");

  DrawOnePlot(rows,
              "dthe",
              "d#theta_{e} fitted offset",
              "Units: mrad",
              "results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/dthe_vs_dpidx.png");

  DrawOnePlot(rows,
              "dpe",
              "dp_{e} fitted offset",
              "Units: 0.1%",
              "results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/dpe_vs_dpidx.png");

  DrawOnePlot(rows,
              "dthp",
              "d#theta_{p} fitted offset",
              "Units: mrad",
              "results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/dthp_vs_dpidx.png");

  DrawOnePlot(rows,
              "dpp",
              "dp_{p} fitted offset",
              "Units: 0.1%",
              "results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/dpp_vs_dpidx.png");

  std::cout << "[INFO] Wrote PNGs to results/PNGs/PlotLsFitOffsetsBySettingAndDpBin/\n";
}
