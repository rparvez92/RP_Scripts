// File: heep_check_v0/macros/plot_heepcheck_reference_from_dpeOnly.C
//
// Purpose:
//   Inspect the reference table built by:
//     results/tables/heepcheck_reference_from_dpeOnly.csv
//
// Produces:
//   - results/PNGs/heepcheck_reference_from_dpeOnly/reference_params_vs_run.png
//   - results/PNGs/heepcheck_reference_from_dpeOnly/reference_residuals_vs_run.png
//   - results/PNGs/heepcheck_reference_from_dpeOnly/reference_sse_vs_run.png
//
// Plot choices:
//   - markers only
//   - no legend
//   - grid on
//   - dashed horizontal zero line on residual pads and SSE pad

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
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Row {
  int run = -1;
  std::string setting;
  std::string coeff_tag;
  double dpe_ref_input = 0.0;
  double dthe_ref = 0.0;
  double dthp_ref = 0.0;
  double dpp_ref = 0.0;
  double r_W = 0.0;
  double r_Em = 0.0;
  double r_Pmz = 0.0;
  double r_Pmy = 0.0;
  double SSE_ref = 0.0;
  int large_dthe_warn = 0;
  int large_dthp_warn = 0;
  int large_dpp_warn = 0;
  int poor_ref_consistency_warn = 0;
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

bool EnsureDir(const std::string& path) {
  try {
    if (!std::filesystem::exists(path)) std::filesystem::create_directories(path);
  } catch (...) {
    return false;
  }
  return true;
}

bool LoadCSV(const std::string& path, std::vector<Row>& out) {
  std::ifstream fin(path);
  if (!fin) {
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
  const std::vector<std::string> req = {
    "run","setting","coeff_tag","dpe_ref_input","dthe_ref","dthp_ref","dpp_ref",
    "r_W","r_Em","r_Pmz","r_Pmy","SSE_ref",
    "large_dthe_warn","large_dthp_warn","large_dpp_warn","poor_ref_consistency_warn"
  };
  for (const auto& c : req) {
    if (!idx.count(c)) {
      std::cerr << "[ERROR] Missing required column '" << c << "' in " << path << "\n";
      return false;
    }
  }

  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto tok = SplitCSVSimple(line);
    auto getS = [&](const std::string& k)->std::string {
      auto it = idx.find(k);
      if (it == idx.end() || it->second >= (int)tok.size()) return "";
      return tok[it->second];
    };
    Row r;
    if (!ToInt(getS("run"), r.run)) continue;
    r.setting = getS("setting");
    r.coeff_tag = getS("coeff_tag");
    ToDouble(getS("dpe_ref_input"), r.dpe_ref_input);
    ToDouble(getS("dthe_ref"), r.dthe_ref);
    ToDouble(getS("dthp_ref"), r.dthp_ref);
    ToDouble(getS("dpp_ref"), r.dpp_ref);
    ToDouble(getS("r_W"), r.r_W);
    ToDouble(getS("r_Em"), r.r_Em);
    ToDouble(getS("r_Pmz"), r.r_Pmz);
    ToDouble(getS("r_Pmy"), r.r_Pmy);
    ToDouble(getS("SSE_ref"), r.SSE_ref);
    ToInt(getS("large_dthe_warn"), r.large_dthe_warn);
    ToInt(getS("large_dthp_warn"), r.large_dthp_warn);
    ToInt(getS("large_dpp_warn"), r.large_dpp_warn);
    ToInt(getS("poor_ref_consistency_warn"), r.poor_ref_consistency_warn);
    out.push_back(r);
  }
  std::sort(out.begin(), out.end(), [](const Row& a, const Row& b){ return a.run < b.run; });
  return true;
}

void ApplyGraphStyle(TGraph* g) {
  if (!g) return;
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->SetLineColor(kBlack);
  g->SetMarkerColor(kBlack);
}

void SetPadStyle() {
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.11);
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

void DrawZeroLine(const std::vector<double>& x) {
  if (x.empty()) return;
  const double xmin = *std::min_element(x.begin(), x.end());
  const double xmax = *std::max_element(x.begin(), x.end());
  TLine* line = new TLine(xmin, 0.0, xmax, 0.0);
  line->SetLineStyle(2);
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
  line->Draw("same");
}

void DrawWarningsAnnotation(const std::vector<Row>& rows) {
  int n1=0,n2=0,n3=0,n4=0;
  for (const auto& r : rows) {
    if (r.large_dthe_warn) ++n1;
    if (r.large_dthp_warn) ++n2;
    if (r.large_dpp_warn) ++n3;
    if (r.poor_ref_consistency_warn) ++n4;
  }
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.16, 0.84, Form("warn counts: dthe=%d dthp=%d dpp=%d poor=%d", n1,n2,n3,n4));
}

template <typename Getter>
void DrawGraphVsRun(const std::vector<Row>& rows,
                    Getter getter,
                    const std::string& title,
                    const std::string& ytitle,
                    bool drawZero = false) {
  std::vector<double> x, y;
  for (const auto& r : rows) {
    x.push_back(r.run);
    y.push_back(getter(r));
  }
  TGraph* g = new TGraph((int)x.size(), x.data(), y.data());
  ApplyGraphStyle(g);
  g->Draw("AP");
  SetAxisTitles(g, "Run", ytitle);
  DrawPadTitle(title);
  if (drawZero) DrawZeroLine(x);
}

} // namespace

void plot_heepcheck_reference_from_dpeOnly(
    const char* input_csv = "results/tables/heepcheck_reference_from_dpeOnly.csv",
    const char* out_dir   = "results/PNGs/heepcheck_reference_from_dpeOnly/")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const std::string incsv = input_csv ? input_csv : "results/tables/heepcheck_reference_from_dpeOnly.csv";
  const std::string outdir = out_dir ? out_dir : "results/PNGs/heepcheck_reference_from_dpeOnly/";

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

  {
    TCanvas* c = new TCanvas("c_ref_params","reference parameters vs run",1400,1000);
    c->Divide(2,2);
    c->cd(1); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.dpe_ref_input; }, "reference parameters vs run", "dpe_ref [0.1%]", false);
    c->cd(2); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.dthe_ref; },      "reference parameters vs run", "dthe_ref [mrad]", false);
    c->cd(3); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.dthp_ref; },      "reference parameters vs run", "dthp_ref [mrad]", false);
    c->cd(4); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.dpp_ref; },       "reference parameters vs run", "dpp_ref [0.1%]", false);
    c->SaveAs((outdir + "/reference_params_vs_run.png").c_str());
    std::cout << "[INFO] Wrote " << outdir + "/reference_params_vs_run.png" << "\n";
    delete c;
  }

  {
    TCanvas* c = new TCanvas("c_ref_resids","reference residuals vs run",1400,1000);
    c->Divide(2,2);
    c->cd(1); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.r_W; },   "reference residuals vs run", "r_W [MeV]", true);
    c->cd(2); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.r_Em; },  "reference residuals vs run", "r_Em [MeV]", true);
    c->cd(3); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.r_Pmz; }, "reference residuals vs run", "r_Pmz [MeV]", true);
    c->cd(4); SetPadStyle(); DrawGraphVsRun(rows, [](const Row& r){ return r.r_Pmy; }, "reference residuals vs run", "r_Pmy [MeV]", true);
    c->SaveAs((outdir + "/reference_residuals_vs_run.png").c_str());
    std::cout << "[INFO] Wrote " << outdir + "/reference_residuals_vs_run.png" << "\n";
    delete c;
  }

  {
    TCanvas* c = new TCanvas("c_ref_sse","reference SSE vs run",1000,700);
    c->cd(); SetPadStyle();
    DrawGraphVsRun(rows, [](const Row& r){ return r.SSE_ref; }, "reference SSE vs run", "SSE_ref", true);
    DrawWarningsAnnotation(rows);
    c->SaveAs((outdir + "/reference_sse_vs_run.png").c_str());
    std::cout << "[INFO] Wrote " << outdir + "/reference_sse_vs_run.png" << "\n";
    delete c;
  }

  std::cout << "[INFO] Done.\n";
}
