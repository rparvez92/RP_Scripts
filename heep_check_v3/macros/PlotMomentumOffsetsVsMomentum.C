// PlotMomentumOffsetsVsMomentum.C
//
// Plot fitted HMS/SHMS momentum offsets versus their central momenta for the
// standard combined 4+5-pass fit.
//
// Inputs:
//   results/tables/Combined4and5Pass/Combined4and5PassErrors_params.csv
//
// Outputs:
//   results/PNGs/Combined4and5Pass/dpp_vs_hms_p.png
//   results/PNGs/Combined4and5Pass/dpe_vs_shms_p.png
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/PlotMomentumOffsetsVsMomentum.C'

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
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
const std::map<std::string, double> kHmsP = {
  {"A", 2.28}, {"B", 3.75}, {"C", 4.72}, {"E", 4.72}, {"D", 4.72}
};
const std::map<std::string, double> kShmsP = {
  {"A", -7.07}, {"B", -5.67}, {"C", -6.72}, {"E", -6.6048}, {"D", -6.384}
};
const std::map<std::string, double> kHmsVisualDx = {
  {"A", 0.00}, {"B", 0.00}, {"C", -0.035}, {"E", 0.00}, {"D", 0.035}
};

struct ParamErr {
  double central = NAN;
  double err = NAN;
};

struct CaseDef {
  std::string label;
  TString table;
  TString outDir;
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

static std::string Field(const std::vector<std::string>& row,
                         const std::map<std::string, int>& col,
                         const std::string& key) {
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

static std::map<std::string, ParamErr> ReadParamErrors(const TString& path) {
  std::map<std::string, ParamErr> out;
  std::ifstream f(path.Data());
  if (!f) {
    std::cerr << "[ERROR] Cannot open " << path << "\n";
    return out;
  }
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
    out[name] = e;
  }
  return out;
}

static ParamErr GetParam(const std::map<std::string, ParamErr>& pars, const std::string& name) {
  auto it = pars.find(name);
  if (it == pars.end()) {
    std::cerr << "[WARN] Missing parameter " << name << "\n";
    return ParamErr();
  }
  return it->second;
}

static void Range(const std::vector<double>& y, const std::vector<double>& ey, double& ymin, double& ymax) {
  ymin = std::numeric_limits<double>::infinity();
  ymax = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < y.size(); ++i) {
    if (!std::isfinite(y[i])) continue;
    const double e = (i < ey.size() && std::isfinite(ey[i])) ? ey[i] : 0.0;
    ymin = std::min(ymin, y[i] - e);
    ymax = std::max(ymax, y[i] + e);
  }
  ymin = std::min(ymin, 0.0);
  ymax = std::max(ymax, 0.0);
  if (!std::isfinite(ymin) || !std::isfinite(ymax) || ymin == ymax) {
    ymin = -1.0;
    ymax = 1.0;
  }
  const double pad = 0.18 * (ymax - ymin);
  ymin -= pad;
  ymax += pad;
}

static void DrawSettingLabels(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& ey,
                              double ymin,
                              double ymax) {
  TLatex lat;
  lat.SetTextAlign(22);
  lat.SetTextSize(0.035);
  lat.SetTextFont(62);
  const double dy = 0.055 * (ymax - ymin);
  for (int i = 0; i < (int)kSettings.size(); ++i) {
    double ypos = y[i] + ((std::isfinite(ey[i]) ? ey[i] : 0.0) + dy);
    if (ypos > ymax - dy) ypos = y[i] - ((std::isfinite(ey[i]) ? ey[i] : 0.0) + dy);
    lat.DrawLatex(x[i], ypos, kSettings[i].c_str());
  }
}

static void PlotOne(const CaseDef& cdef,
                    const std::map<std::string, ParamErr>& pars,
                    bool hms,
                    const TString& outName) {
  std::vector<double> x, y, ey;
  x.reserve(kSettings.size());
  y.reserve(kSettings.size());
  ey.reserve(kSettings.size());

  for (const auto& setting : kSettings) {
    const std::string pname = std::string(hms ? "dpp_" : "dpe_") + setting;
    const ParamErr p = GetParam(pars, pname);
    const double x0 = hms ? kHmsP.at(setting) + kHmsVisualDx.at(setting) : kShmsP.at(setting);
    x.push_back(x0);
    y.push_back(p.central);
    ey.push_back(std::isfinite(p.err) ? p.err : 0.0);
  }

  double ymin, ymax;
  Range(y, ey, ymin, ymax);
  double xmin = *std::min_element(x.begin(), x.end());
  double xmax = *std::max_element(x.begin(), x.end());
  const double xpad = std::max(0.12, 0.12 * (xmax - xmin));
  xmin -= xpad;
  xmax += xpad;

  TCanvas canv(Form("c_%s_%s", cdef.label.c_str(), hms ? "hms" : "shms"),
               hms ? "HMS Momentum Offset" : "SHMS Momentum Offset",
               1000, 720);
  canv.SetGrid();
  canv.SetLeftMargin(0.13);
  canv.SetRightMargin(0.04);
  canv.SetBottomMargin(0.12);

  TH1F* frame = new TH1F(Form("h_%s_%s", cdef.label.c_str(), hms ? "hms" : "shms"),
                         hms ? "HMS Momentum Offset" : "SHMS Momentum Offset",
                         100, xmin, xmax);
  frame->SetStats(0);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->SetTitle("Momentum (GeV)");
  frame->GetYaxis()->SetTitle("Momentum Offset (0.1%)");
  frame->GetXaxis()->SetTitleSize(0.048);
  frame->GetYaxis()->SetTitleSize(0.048);
  frame->GetXaxis()->SetLabelSize(0.043);
  frame->GetYaxis()->SetLabelSize(0.043);
  frame->GetYaxis()->SetTitleOffset(1.05);
  frame->Draw();

  TGraphErrors* g = new TGraphErrors();
  for (int i = 0; i < (int)x.size(); ++i) {
    g->SetPoint(i, x[i], y[i]);
    g->SetPointError(i, 0.0, ey[i]);
  }
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.25);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);
  g->Draw("P same");
  DrawSettingLabels(x, y, ey, ymin, ymax);

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.16, 0.84, cdef.label.c_str());
  if (hms) {
    lat.SetTextSize(0.028);
    lat.DrawLatex(0.16, 0.79, "C/E/D are slightly x-offset for visibility");
  }

  canv.SaveAs(cdef.outDir + "/" + outName);
}

static void PlotCase(const CaseDef& cdef) {
  if (gSystem->AccessPathName(cdef.outDir)) gSystem->mkdir(cdef.outDir, kTRUE);
  const auto pars = ReadParamErrors(cdef.table);
  PlotOne(cdef, pars, true, "dpp_vs_hms_p.png");
  PlotOne(cdef, pars, false, "dpe_vs_shms_p.png");
  std::cout << "[INFO] Wrote " << cdef.outDir << "/dpp_vs_hms_p.png\n";
  std::cout << "[INFO] Wrote " << cdef.outDir << "/dpe_vs_shms_p.png\n";
}

} // namespace

void PlotMomentumOffsetsVsMomentum() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const CaseDef combined = {
    "Combined 4+5",
    "results/tables/Combined4and5Pass/Combined4and5PassErrors_params.csv",
    "results/PNGs/Combined4and5Pass"
  };
  PlotCase(combined);
}
