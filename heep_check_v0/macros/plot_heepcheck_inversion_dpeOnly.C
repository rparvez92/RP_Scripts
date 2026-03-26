// File: heep_check_v0/macros/plot_heepcheck_inversion_dpeOnly.C
//
// Purpose:
//   Make plots from:
//     results/tables/heepcheck_inversion_dpeOnly.csv
//
// Produces:
//
//   First set (6 plots total, one per dp_idx):
//     results/PNGs/kin_offsets_dpeOnly/kin_offsets_vs_run_<dp_idx>.png
//
//   Second set (one plot per run):
//     results/PNGs/kin_offsets_dpeOnly/kin_offsets_vs_dpbinCenter_<run>.png
//
// Plot choices:
//   - same 4-panel format as full4
//   - markers only
//   - no legend
//   - grid on
//   - auto y-range
//   - second set uses narrow bins only as points
//   - full-bin (dp_idx=0) shown as dashed blue horizontal reference line
//   - annotate each pad with full-bin value

#include <TCanvas.h>
#include <TGraph.h>
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
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Row {
  std::string fit_model;
  int         de0_fixed = 0;
  int         dthe_fixed = 0;
  int         dthp_fixed = 0;
  int         dpp_fixed = 0;
  std::string sign_convention;

  int         run = -1;
  int         dp_idx = -999;
  std::string dp_label;
  double      dp_lo = 0.0;
  double      dp_hi = 0.0;

  double      dthe = 0.0;
  double      dpe  = 0.0;
  double      dthp = 0.0;
  double      dpp  = 0.0;

  double      sse = 0.0;
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

bool LoadCSV(const std::string& path, std::vector<Row>& rows) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open input CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty input CSV: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);

  const std::vector<std::string> required = {
    "fit_model","de0_fixed","dthe_fixed","dthp_fixed","dpp_fixed","sign_convention",
    "run","dp_idx","dp_label","dp_lo","dp_hi",
    "dthe","dpe","dthp","dpp","sse"
  };

  for (const auto& c : required) {
    if (!idx.count(c)) {
      std::cerr << "[ERROR] Missing required column: " << c << "\n";
      return false;
    }
  }

  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto tok = SplitCSVSimple(line);
    if ((int)tok.size() < (int)hdr.size()) continue;

    Row r;
    r.fit_model = tok[idx.at("fit_model")];
    ToInt(tok[idx.at("de0_fixed")], r.de0_fixed);
    ToInt(tok[idx.at("dthe_fixed")], r.dthe_fixed);
    ToInt(tok[idx.at("dthp_fixed")], r.dthp_fixed);
    ToInt(tok[idx.at("dpp_fixed")], r.dpp_fixed);
    r.sign_convention = tok[idx.at("sign_convention")];

    if (!ToInt(tok[idx.at("run")], r.run)) continue;
    if (!ToInt(tok[idx.at("dp_idx")], r.dp_idx)) continue;
    r.dp_label = tok[idx.at("dp_label")];
    if (!ToDouble(tok[idx.at("dp_lo")], r.dp_lo)) continue;
    if (!ToDouble(tok[idx.at("dp_hi")], r.dp_hi)) continue;
    if (!ToDouble(tok[idx.at("dthe")], r.dthe)) continue;
    if (!ToDouble(tok[idx.at("dpe")],  r.dpe))  continue;
    if (!ToDouble(tok[idx.at("dthp")], r.dthp)) continue;
    if (!ToDouble(tok[idx.at("dpp")],  r.dpp))  continue;
    if (!ToDouble(tok[idx.at("sse")],  r.sse))  continue;

    rows.push_back(r);
  }

  return true;
}

std::vector<int> UniqueRuns(const std::vector<Row>& rows) {
  std::set<int> s;
  for (const auto& r : rows) s.insert(r.run);
  return std::vector<int>(s.begin(), s.end());
}

double DpCenter(const Row& r) {
  return 0.5 * (r.dp_lo + r.dp_hi);
}

void SetPadStyle() {
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.12);
}

void ApplyGraphStyle(TGraph* g) {
  if (!g) return;
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
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
  lat.SetTextSize(0.055);
  lat.DrawLatex(0.14, 0.93, text.c_str());
}

void DrawFullBinAnnotation(double value) {
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.16, 0.84, Form("full-bin = %.3f", value));
}

template <typename Getter>
void DrawGraphVsRun(const std::vector<Row>& rowsForDp,
                    Getter getter,
                    const std::string& title,
                    const std::string& ytitle) {
  std::vector<Row> v = rowsForDp;
  std::sort(v.begin(), v.end(), [](const Row& a, const Row& b){ return a.run < b.run; });

  std::vector<double> x, y;
  for (const auto& r : v) {
    x.push_back(r.run);
    y.push_back(getter(r));
  }

  TGraph* g = new TGraph((int)x.size(), x.data(), y.data());
  ApplyGraphStyle(g);
  g->Draw("AP");
  SetAxisTitles(g, "Run", ytitle);
  DrawPadTitle(title);
}

template <typename Getter>
void DrawGraphVsDpCenter(const std::vector<Row>& rowsForRunNarrow,
                         double fullBinValue,
                         Getter getter,
                         const std::string& title,
                         const std::string& ytitle) {
  std::vector<Row> v = rowsForRunNarrow;
  std::sort(v.begin(), v.end(), [](const Row& a, const Row& b){ return DpCenter(a) < DpCenter(b); });

  std::vector<double> x, y;
  for (const auto& r : v) {
    x.push_back(DpCenter(r));
    y.push_back(getter(r));
  }

  TGraph* g = new TGraph((int)x.size(), x.data(), y.data());
  ApplyGraphStyle(g);
  g->Draw("AP");
  SetAxisTitles(g, "dp bin center", ytitle);
  DrawPadTitle(title);

  if (!x.empty()) {
    const double xmin = *std::min_element(x.begin(), x.end());
    const double xmax = *std::max_element(x.begin(), x.end());
    TLine* line = new TLine(xmin, fullBinValue, xmax, fullBinValue);
    line->SetLineStyle(2);
    line->SetLineColor(kBlue);
    line->SetLineWidth(2);
    line->Draw("same");
  }

  DrawFullBinAnnotation(fullBinValue);
}

bool PlotOneDpIdx(const std::vector<Row>& rows, int dp_idx, const std::string& outdir) {
  std::vector<Row> use;
  for (const auto& r : rows) if (r.dp_idx == dp_idx) use.push_back(r);

  if (use.empty()) {
    std::cerr << "[WARN] No rows found for dp_idx=" << dp_idx << "\n";
    return false;
  }

  TCanvas* c = new TCanvas(Form("c_run_dpeOnly_%d", dp_idx),
                           Form("kinematic offsets vs run (dp_idx = %d)", dp_idx),
                           1400, 1000);
  c->Divide(2,2);

  c->cd(1); SetPadStyle();
  DrawGraphVsRun(use, [](const Row& r){ return r.dthe; },
                 Form("kinematic offsets vs run (dp_idx = %d)", dp_idx), "dthe [1 mrad]");
  c->cd(2); SetPadStyle();
  DrawGraphVsRun(use, [](const Row& r){ return r.dpe; },
                 Form("kinematic offsets vs run (dp_idx = %d)", dp_idx), "dpe [0.1%]");
  c->cd(3); SetPadStyle();
  DrawGraphVsRun(use, [](const Row& r){ return r.dthp; },
                 Form("kinematic offsets vs run (dp_idx = %d)", dp_idx), "dthp [1 mrad]");
  c->cd(4); SetPadStyle();
  DrawGraphVsRun(use, [](const Row& r){ return r.dpp; },
                 Form("kinematic offsets vs run (dp_idx = %d)", dp_idx), "dpp [0.1%]");

  const std::string outpng = outdir + "/kin_offsets_vs_run_" + std::to_string(dp_idx) + ".png";
  c->SaveAs(outpng.c_str());
  std::cout << "[INFO] Wrote " << outpng << "\n";
  delete c;
  return true;
}

bool PlotOneRunDpCenter(const std::vector<Row>& rows, int run, const std::string& outdir) {
  std::vector<Row> narrow;
  const Row* full = nullptr;

  for (const auto& r : rows) {
    if (r.run != run) continue;
    if (r.dp_idx == 0) full = &r;
    else if (r.dp_idx >= 1 && r.dp_idx <= 5) narrow.push_back(r);
  }

  if (!full || narrow.empty()) {
    std::cerr << "[WARN] Missing full or narrow bins for run=" << run << "\n";
    return false;
  }

  TCanvas* c = new TCanvas(Form("c_dpcenter_dpeOnly_%d", run),
                           Form("kinematic offsets vs dp bin center (run %d)", run),
                           1400, 1000);
  c->Divide(2,2);

  c->cd(1); SetPadStyle();
  DrawGraphVsDpCenter(narrow, full->dthe, [](const Row& r){ return r.dthe; },
                      Form("kinematic offsets vs dp bin center (run %d)", run), "dthe [1 mrad]");
  c->cd(2); SetPadStyle();
  DrawGraphVsDpCenter(narrow, full->dpe, [](const Row& r){ return r.dpe; },
                      Form("kinematic offsets vs dp bin center (run %d)", run), "dpe [0.1%]");
  c->cd(3); SetPadStyle();
  DrawGraphVsDpCenter(narrow, full->dthp, [](const Row& r){ return r.dthp; },
                      Form("kinematic offsets vs dp bin center (run %d)", run), "dthp [1 mrad]");
  c->cd(4); SetPadStyle();
  DrawGraphVsDpCenter(narrow, full->dpp, [](const Row& r){ return r.dpp; },
                      Form("kinematic offsets vs dp bin center (run %d)", run), "dpp [0.1%]");

  const std::string outpng = outdir + "/kin_offsets_vs_dpbinCenter_" + std::to_string(run) + ".png";
  c->SaveAs(outpng.c_str());
  std::cout << "[INFO] Wrote " << outpng << "\n";
  delete c;
  return true;
}

} // namespace

void plot_heepcheck_inversion_dpeOnly(
    const char* input_csv = "results/tables/heepcheck_inversion_dpeOnly.csv",
    const char* out_dir   = "results/PNGs/kin_offsets_dpeOnly/")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::string incsv = input_csv ? input_csv : "results/tables/heepcheck_inversion_dpeOnly.csv";
  const std::string outdir = out_dir ? out_dir : "results/PNGs/kin_offsets_dpeOnly/";

  if (!EnsureDir(outdir)) {
    std::cerr << "[ERROR] Could not create output directory: " << outdir << "\n";
    return;
  }

  std::vector<Row> rows;
  if (!LoadCSV(incsv, rows)) return;

  if (rows.empty()) {
    std::cerr << "[ERROR] No valid rows loaded from: " << incsv << "\n";
    return;
  }

  std::cout << "[INFO] Loaded " << rows.size() << " rows from " << incsv << "\n";

  for (int dp_idx = 0; dp_idx <= 5; ++dp_idx) {
    PlotOneDpIdx(rows, dp_idx, outdir);
  }

  const auto runs = UniqueRuns(rows);
  for (int run : runs) {
    PlotOneRunDpCenter(rows, run, outdir);
  }

  std::cout << "[INFO] Done.\n";
}
