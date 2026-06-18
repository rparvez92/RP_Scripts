// PlotCombined4and5Pass.C
//
// Plot the combined 4+5-pass fit, toy errors, and residuals.
//
// Inputs:
//   results/tables/Combined4and5Pass/Combined4and5PassFitParams.csv
//   results/tables/Combined4and5Pass/Combined4and5PassErrors_params.csv
//   results/tables/Combined4and5Pass/Combined4and5PassResiduals_WithToyErrors.csv
//
// Outputs:
//   results/PNGs/Combined4and5Pass/chi2ndf.png
//   results/PNGs/Combined4and5Pass/kinOffsets.png
//   results/PNGs/Combined4and5Pass/residuals.png
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/PlotCombined4and5Pass.C'

#include <TBox.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {

const std::vector<std::string> kSettings = {"A", "B", "C", "E", "D"};
const std::vector<std::string> kVars = {"W", "Em", "Pmz", "Pmx"};
const std::map<std::string, double> kX = {{"A", 0.0}, {"B", 1.0}, {"C", 2.0}, {"E", 3.0}, {"D", 4.0}};
const std::map<std::string, std::string> kBandLabels = {
  {"A", "4pass, P_{SHMS}=-7.07"},
  {"B", "4pass, P_{SHMS}=-5.67"},
  {"C", "5pass, P_{SHMS}=-6.720"},
  {"E", "5pass, P_{SHMS}=-6.605"},
  {"D", "5pass, P_{SHMS}=-6.384"}
};

struct ParamErr {
  double central = NAN;
  double err = NAN;
};

struct ResidualRow {
  std::string pass;
  std::string setting;
  std::string variable;
  double residual = NAN;
  double residualErr = NAN;
};

struct FitSummary {
  double chi2ndf = NAN;
  std::map<std::string, double> params;
};

static std::string Trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
  size_t e = s.size();
  while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
  if (e > b && s[b] == '"' && s[e - 1] == '"') { ++b; --e; }
  return s.substr(b, e - b);
}

static std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
  return s;
}

static std::vector<std::string> SplitCsvLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inQuotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      inQuotes = !inQuotes;
    } else if (c == ',' && !inQuotes) {
      out.push_back(Trim(cur));
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(Trim(cur));
  return out;
}

static std::map<std::string, int> HeaderMap(const std::vector<std::string>& h) {
  std::map<std::string, int> out;
  for (int i = 0; i < (int)h.size(); ++i) out[Lower(h[i])] = i;
  return out;
}

static std::string Field(const std::vector<std::string>& row, const std::map<std::string, int>& col, const std::string& key) {
  auto it = col.find(Lower(key));
  if (it == col.end()) return "";
  const int i = it->second;
  return (i >= 0 && i < (int)row.size()) ? row[i] : "";
}

static double ToDouble(const std::string& s) {
  const std::string t = Trim(s);
  if (t.empty()) return NAN;
  char* end = nullptr;
  const double v = std::strtod(t.c_str(), &end);
  return (end && *end == '\0') ? v : NAN;
}

static FitSummary ReadFit(const TString& path) {
  FitSummary out;
  std::ifstream f(path.Data());
  if (!f) { std::cerr << "[ERROR] Cannot open " << path << "\n"; return out; }
  std::string line;
  if (!std::getline(f, line)) return out;
  const auto headers = SplitCsvLine(line);
  const auto col = HeaderMap(headers);
  if (!std::getline(f, line)) return out;
  const auto row = SplitCsvLine(line);
  out.chi2ndf = ToDouble(Field(row, col, "chi2_ndf"));
  for (const auto& h : headers) {
    if (h.rfind("dthe", 0) == 0 || h.rfind("dpe", 0) == 0 || h.rfind("dpp", 0) == 0 || h == "dthp") {
      out.params[h] = ToDouble(Field(row, col, h));
    }
  }
  return out;
}

static std::map<std::string, ParamErr> ReadParamErrors(const TString& path, const FitSummary& fit) {
  std::map<std::string, ParamErr> out;
  std::ifstream f(path.Data());
  if (!f) { std::cerr << "[ERROR] Cannot open " << path << "\n"; return out; }
  std::string line;
  if (!std::getline(f, line)) return out;
  const auto col = HeaderMap(SplitCsvLine(line));
  while (std::getline(f, line)) {
    if (Trim(line).empty()) continue;
    const auto row = SplitCsvLine(line);
    const std::string name = Field(row, col, "parameter");
    ParamErr e;
    e.central = ToDouble(Field(row, col, "central_value"));
    e.err = ToDouble(Field(row, col, "toy_std"));
    if (!std::isfinite(e.central) && fit.params.count(name)) e.central = fit.params.at(name);
    out[name] = e;
  }
  return out;
}

static std::vector<ResidualRow> ReadResiduals(const TString& path) {
  std::vector<ResidualRow> out;
  std::ifstream f(path.Data());
  if (!f) { std::cerr << "[ERROR] Cannot open " << path << "\n"; return out; }
  std::string line;
  if (!std::getline(f, line)) return out;
  const auto col = HeaderMap(SplitCsvLine(line));
  while (std::getline(f, line)) {
    if (Trim(line).empty()) continue;
    const auto row = SplitCsvLine(line);
    ResidualRow r;
    r.pass = Field(row, col, "pass");
    r.setting = Field(row, col, "setting");
    r.variable = Field(row, col, "variable");
    r.residual = ToDouble(Field(row, col, "residual"));
    r.residualErr = ToDouble(Field(row, col, "measured_err"));
    out.push_back(r);
  }
  return out;
}

static TH1F* MakeCategoryFrame(const char* name, const char* title, const char* ytitle, double ymin, double ymax) {
  TH1F* h = new TH1F(name, title, 5, -0.5, 4.5);
  for (int i = 0; i < 5; ++i) h->GetXaxis()->SetBinLabel(i + 1, kSettings[i].c_str());
  h->GetXaxis()->SetTitle("Setting");
  h->GetYaxis()->SetTitle(ytitle);
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->SetStats(0);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetXaxis()->SetTitleSize(0.052);
  h->GetYaxis()->SetTitleSize(0.052);
  return h;
}

static void Range(const std::vector<double>& y, const std::vector<double>& e, double& ymin, double& ymax) {
  ymin = std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < y.size(); ++i) {
    if (!std::isfinite(y[i])) continue;
    const double err = (i < e.size() && std::isfinite(e[i])) ? e[i] : 0.0;
    ymin = std::min(ymin, y[i] - err);
    ymax = std::max(ymax, y[i] + err);
  }
  ymin = std::min(ymin, 0.0);
  ymax = std::max(ymax, 0.0);
  if (!std::isfinite(ymin) || !std::isfinite(ymax) || ymin == ymax) { ymin = -1; ymax = 1; }
  const double pad = 0.18 * (ymax - ymin);
  ymin -= pad;
  ymax += pad;
}

static ParamErr GetParam(const std::map<std::string, ParamErr>& errs, const FitSummary& fit, const std::string& name) {
  auto it = errs.find(name);
  if (it != errs.end()) return it->second;
  ParamErr e;
  auto ip = fit.params.find(name);
  e.central = (ip == fit.params.end()) ? NAN : ip->second;
  e.err = 0.0;
  return e;
}

static void DrawBandLabels(double ymin, double ymax) {
  const int fills[5] = {kTeal - 9, kOrange - 9, kViolet - 9, kPink - 9, kGreen - 9};
  for (int i = 0; i < 5; ++i) {
    TBox* b = new TBox(i - 0.5, ymin, i + 0.5, ymax);
    b->SetFillColorAlpha(fills[i], 0.35);
    b->SetLineColor(0);
    b->Draw("same");
  }
  TLatex lat;
  lat.SetTextAlign(22);
  lat.SetTextSize(0.032);
  lat.SetTextFont(62);
  const double y = ymax - 0.08 * (ymax - ymin);
  for (int i = 0; i < 5; ++i) lat.DrawLatex(i, y, kBandLabels.at(kSettings[i]).c_str());
}

static void PlotChi2(const FitSummary& fit, const TString& out) {
  TCanvas c("c_comb45_chi2", "Combined 4+5 chi2", 1000, 720);
  c.SetGrid();
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);
  TH1F* frame = MakeCategoryFrame("h_chi2_frame", "Combined 4+5-pass fit", "#chi^{2}/ndf", 0, std::max(1.0, 1.35 * fit.chi2ndf));
  frame->Draw();
  TGraphErrors* g = new TGraphErrors();
  for (int i = 0; i < 5; ++i) {
    g->SetPoint(i, i, fit.chi2ndf);
    g->SetPointError(i, 0, 0);
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.25);
  g->SetMarkerColor(kBlack);
  g->Draw("P same");
  c.SaveAs(out);
}

static void DrawParamPad(const std::vector<std::string>& names,
                         const std::vector<double>& xs,
                         const std::map<std::string, ParamErr>& errs,
                         const FitSummary& fit,
                         const char* title,
                         const char* ytitle) {
  std::vector<double> y, ey;
  for (const auto& name : names) {
    const ParamErr e = GetParam(errs, fit, name);
    y.push_back(e.central);
    ey.push_back(std::isfinite(e.err) ? e.err : 0.0);
  }
  double ymin, ymax;
  Range(y, ey, ymin, ymax);
  TH1F* frame = MakeCategoryFrame(Form("h_%s", title), title, ytitle, ymin, ymax);
  frame->GetYaxis()->SetTitleOffset(0.95);
  frame->Draw();
  TGraphErrors* g = new TGraphErrors();
  for (int i = 0; i < (int)names.size(); ++i) {
    g->SetPoint(i, xs[i], y[i]);
    g->SetPointError(i, 0, ey[i]);
  }
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);
  g->Draw("P same");
}

static void PlotKin(const FitSummary& fit, const std::map<std::string, ParamErr>& errs, const TString& out) {
  TCanvas c("c_comb45_kin", "Combined 4+5 kinematic offsets", 1400, 1000);
  c.Divide(2, 2);
  const std::vector<double> x5 = {0, 1, 2, 3, 4};
  c.cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
  DrawParamPad({"dthe_A","dthe_B","dthe_C","dthe_E","dthe_D"}, x5, errs, fit,
               "Combined 4+5 Offsets: #Delta#theta_{e}", "#Delta#theta_{e} (mrad)");
  c.cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
  DrawParamPad({"dpe_A","dpe_B","dpe_C","dpe_E","dpe_D"}, x5, errs, fit,
               "Combined 4+5 Offsets: #Deltap_{e}", "#Deltap_{e} (0.1%)");
  c.cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
  DrawParamPad({"dthp"}, {2.0}, errs, fit,
               "Combined 4+5 Offsets: #Delta#theta_{p}", "#Delta#theta_{p} (mrad)");
  c.cd(4); gPad->SetGrid(); gPad->SetLeftMargin(0.19); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13);
  DrawParamPad({"dpp_A","dpp_B","dpp_C","dpp_E","dpp_D"}, x5, errs, fit,
               "Combined 4+5 Offsets: #Deltap_{p}", "#Deltap_{p} (0.1%)");
  c.SaveAs(out);
}

static void PlotResiduals(const std::vector<ResidualRow>& rows, const TString& out) {
  std::vector<double> y, ey;
  for (const auto& r : rows) { y.push_back(r.residual); ey.push_back(r.residualErr); }
  double ymin, ymax;
  Range(y, ey, ymin, ymax);

  TCanvas c("c_comb45_residuals", "Combined 4+5 residuals", 1500, 700);
  c.SetGrid();
  c.SetLeftMargin(0.09);
  c.SetRightMargin(0.03);
  c.SetBottomMargin(0.12);
  TH1F* frame = MakeCategoryFrame("h_resid_frame", "Combined 4+5-pass residuals", "Residual (MeV)", ymin, ymax);
  frame->GetYaxis()->SetTitleOffset(0.85);
  frame->Draw();
  DrawBandLabels(ymin, ymax);
  frame->Draw("axis same");

  TLine* zero = new TLine(-0.5, 0, 4.5, 0);
  zero->SetLineColor(kGray + 2);
  zero->SetLineStyle(2);
  zero->Draw("same");

  const std::map<std::string, int> colors = {{"W", kBlack}, {"Em", kBlue + 1}, {"Pmz", kRed + 1}, {"Pmx", kGreen + 2}};
  const std::map<std::string, int> markers = {{"W", 20}, {"Em", 21}, {"Pmz", 22}, {"Pmx", 33}};
  const std::map<std::string, double> dx = {{"W", -0.18}, {"Em", -0.06}, {"Pmz", 0.06}, {"Pmx", 0.18}};
  TLegend* leg = new TLegend(0.71, 0.12, 0.95, 0.24);
  leg->SetNColumns(4);
  leg->SetBorderSize(1);
  leg->SetFillColorAlpha(kWhite, 0.75);
  for (const auto& var : kVars) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for (const auto& r : rows) {
      if (r.variable != var) continue;
      g->SetPoint(ip, kX.at(r.setting) + dx.at(var), r.residual);
      g->SetPointError(ip, 0, r.residualErr);
      ++ip;
    }
    g->SetMarkerColor(colors.at(var));
    g->SetLineColor(colors.at(var));
    g->SetMarkerStyle(markers.at(var));
    g->SetMarkerSize(var == "Pmx" ? 1.5 : 1.1);
    g->SetLineWidth(2);
    g->Draw("P same");
    leg->AddEntry(g, var.c_str(), "p");
  }
  leg->Draw();
  c.SaveAs(out);
}

} // namespace

void PlotCombined4and5Pass() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  const TString inDir = "results/tables/Combined4and5Pass";
  const TString outDir = "results/PNGs/Combined4and5Pass";
  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);

  const FitSummary fit = ReadFit(inDir + "/Combined4and5PassFitParams.csv");
  const auto errs = ReadParamErrors(inDir + "/Combined4and5PassErrors_params.csv", fit);
  const auto residuals = ReadResiduals(inDir + "/Combined4and5PassResiduals_WithToyErrors.csv");

  PlotChi2(fit, outDir + "/chi2ndf.png");
  PlotKin(fit, errs, outDir + "/kinOffsets.png");
  PlotResiduals(residuals, outDir + "/residuals.png");
  std::cout << "[INFO] Wrote " << outDir << "/chi2ndf.png\n";
  std::cout << "[INFO] Wrote " << outDir << "/kinOffsets.png\n";
  std::cout << "[INFO] Wrote " << outDir << "/residuals.png\n";
}
