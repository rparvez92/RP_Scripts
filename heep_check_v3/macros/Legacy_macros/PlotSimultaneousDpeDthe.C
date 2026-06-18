// PlotSimultaneousDpeDthe.C
//
// Plot combined dpe/dthe fit results for either 4-pass or 5-pass.
//
// Inputs for <pass> = 4pass or 5pass:
//   results/tables/SimultaneousDpeDtheFitParams_<pass>.csv
//   results/tables/SimultaneousDpeDtheErrors_<pass>_params.csv
//   results/tables/SimultaneousDpeDtheResiduals_<pass>.csv
//
// Outputs:
//   results/PNGs/Combined/chi2ndf_<pass>.png
//   results/PNGs/Combined/kinOffsets_<pass>.png
//   results/PNGs/Combined/residuals_<pass>.png
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/PlotSimultaneousDpeDthe.C("4pass")'
//   root -l -b -q 'macros/PlotSimultaneousDpeDthe.C("5pass")'

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

struct PlotConfig {
  TString passTag;
  TString titlePass;
  std::vector<double> xCenters;
  std::vector<double> settingParamXCenters;
  std::vector<double> commonParamXCenters;
  std::vector<std::string> dtheNames;
  std::vector<std::string> dpeNames;
  std::vector<std::string> dppNames;
  std::string dthpName = "dthp";
};

struct ParamErr {
  double central = std::numeric_limits<double>::quiet_NaN();
  double minus = std::numeric_limits<double>::quiet_NaN();
  double plus = std::numeric_limits<double>::quiet_NaN();
  double toy_std = std::numeric_limits<double>::quiet_NaN();
};

struct FitSummary {
  double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
  std::map<std::string, double> params;
};

struct ResidualRow {
  std::string setting;
  std::string dp_label;
  double delta_center = std::numeric_limits<double>::quiet_NaN();
  std::string variable;
  double residual = std::numeric_limits<double>::quiet_NaN();
  double residual_err = std::numeric_limits<double>::quiet_NaN();
};

static std::string Trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  if (e > b && s[b] == '"' && s[e - 1] == '"') {
    ++b;
    --e;
  }
  return s.substr(b, e - b);
}

static std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

static TString NormalizePassTag(const char* passC) {
  TString tag(passC ? passC : "");
  tag.ToLower();
  if (tag == "4" || tag == "4pass" || tag == "pass4") return "4pass";
  if (tag == "5" || tag == "5pass" || tag == "pass5") return "5pass";
  std::cerr << "[ERROR] Unknown pass argument '" << tag << "'. Use 4pass or 5pass.\n";
  return "";
}

static PlotConfig MakeConfig(const TString& passTag) {
  PlotConfig cfg;
  cfg.passTag = passTag;
  if (passTag == "4pass") {
    cfg.titlePass = "4-pass";
    cfg.xCenters = {-0.5};
    cfg.settingParamXCenters = {-0.56, -0.44};
    cfg.commonParamXCenters = {-0.5};
    cfg.dtheNames = {"dthe_b1_A", "dthe_b1_B"};
    cfg.dpeNames = {"dpe_b1_A", "dpe_b1_B"};
    cfg.dppNames = {"dpp_b1_A", "dpp_b1_B"};
    return cfg;
  }
  cfg.titlePass = "5-pass";
  cfg.xCenters = {1.0, 3.0, 6.0};
  cfg.settingParamXCenters = cfg.xCenters;
  cfg.commonParamXCenters = {3.0};
  cfg.dtheNames = {"dthe_C", "dthe_E", "dthe_D"};
  cfg.dpeNames = {"dpe_C", "dpe_E", "dpe_D"};
  cfg.dppNames = {"dpp"};
  return cfg;
}

static std::vector<std::string> SplitCsvLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inQuotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
        cur.push_back('"');
        ++i;
      } else {
        inQuotes = !inQuotes;
      }
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

static double ToDouble(const std::string& s) {
  const std::string t = Trim(s);
  if (t.empty()) return std::numeric_limits<double>::quiet_NaN();
  char* end = nullptr;
  const double v = std::strtod(t.c_str(), &end);
  if (!(end && *end == '\0')) return std::numeric_limits<double>::quiet_NaN();
  return v;
}

static std::map<std::string, int> HeaderMap(const std::vector<std::string>& headers) {
  std::map<std::string, int> out;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) out[Lower(headers[i])] = i;
  return out;
}

static std::string GetField(const std::vector<std::string>& fields,
                            const std::map<std::string, int>& col,
                            const std::string& key) {
  const auto it = col.find(Lower(key));
  if (it == col.end()) return "";
  const int i = it->second;
  return (i >= 0 && i < static_cast<int>(fields.size())) ? fields[i] : "";
}

static FitSummary ReadFitSummary(const TString& path) {
  FitSummary out;
  std::ifstream fin(path.Data());
  if (!fin) {
    std::cerr << "[ERROR] Cannot open " << path << "\n";
    return out;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) return out;
  const auto headers = SplitCsvLine(headerLine);
  const auto col = HeaderMap(headers);

  std::string line;
  if (!std::getline(fin, line)) return out;
  const auto fields = SplitCsvLine(line);

  out.chi2_ndf = ToDouble(GetField(fields, col, "chi2_ndf"));
  for (const auto& h : headers) {
    if (h.rfind("dthe", 0) == 0 || h.rfind("dpe", 0) == 0 || h == "dthp" || h == "dpp") {
      out.params[h] = ToDouble(GetField(fields, col, h));
    }
  }
  return out;
}

static std::map<std::string, ParamErr> ReadParamErrors(const TString& path) {
  std::map<std::string, ParamErr> out;
  std::ifstream fin(path.Data());
  if (!fin) {
    std::cerr << "[ERROR] Cannot open " << path << "\n";
    return out;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) return out;
  const auto col = HeaderMap(SplitCsvLine(headerLine));

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    const std::string name = GetField(fields, col, "parameter");
    if (name.empty()) continue;
    ParamErr e;
    e.central = ToDouble(GetField(fields, col, "central_value"));
    e.minus = ToDouble(GetField(fields, col, "minus_err"));
    e.plus = ToDouble(GetField(fields, col, "plus_err"));
    e.toy_std = ToDouble(GetField(fields, col, "toy_std"));
    out[name] = e;
  }
  return out;
}

static std::vector<ResidualRow> ReadResiduals(const TString& path) {
  std::vector<ResidualRow> out;
  std::ifstream fin(path.Data());
  if (!fin) {
    std::cerr << "[ERROR] Cannot open " << path << "\n";
    return out;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) return out;
  const auto col = HeaderMap(SplitCsvLine(headerLine));

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    ResidualRow r;
    r.setting = GetField(fields, col, "setting");
    r.dp_label = GetField(fields, col, "dp_label");
    r.delta_center = ToDouble(GetField(fields, col, "delta_center"));
    r.variable = GetField(fields, col, "variable");
    r.residual = ToDouble(GetField(fields, col, "residual"));
    r.residual_err = ToDouble(GetField(fields, col, "residual_err_Total"));
    out.push_back(r);
  }
  return out;
}

static void SetGraphBase(TGraph* g, const char* title, const char* ytitle) {
  g->SetTitle(title);
  g->GetXaxis()->SetTitle("#delta bin center (%)");
  g->GetYaxis()->SetTitle(ytitle);
  g->GetXaxis()->SetTitleSize(0.052);
  g->GetYaxis()->SetTitleSize(0.052);
  g->GetXaxis()->SetLabelSize(0.045);
  g->GetYaxis()->SetLabelSize(0.045);
  g->GetYaxis()->SetTitleOffset(1.0);
}

static void RangeFromValues(const std::vector<double>& y,
                            const std::vector<double>& ey,
                            double& ymin,
                            double& ymax) {
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

static std::pair<double, double> XLimits(const std::vector<double>& x) {
  double xmin = *std::min_element(x.begin(), x.end());
  double xmax = *std::max_element(x.begin(), x.end());
  const double pad = std::max(0.5, 0.15 * (xmax - xmin));
  return {xmin - pad, xmax + pad};
}

static ParamErr ErrorFor(const std::map<std::string, ParamErr>& errs,
                         const FitSummary& fit,
                         const std::string& name) {
  auto it = errs.find(name);
  if (it != errs.end()) return it->second;
  ParamErr e;
  auto ip = fit.params.find(name);
  e.central = (ip == fit.params.end() ? std::numeric_limits<double>::quiet_NaN() : ip->second);
  e.toy_std = 0.0;
  return e;
}

static void PlotChi2(const PlotConfig& cfg, const FitSummary& fit, const TString& outPath) {
  const int n = static_cast<int>(cfg.xCenters.size());
  std::vector<double> y(n, fit.chi2_ndf), ex(n, 0.0), ey(n, 0.0);

  TCanvas c("c_chi2ndf_combined", "chi2/ndf combined", 1000, 750);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);
  c.SetGrid();

  TGraphErrors g(n, cfg.xCenters.data(), y.data(), ex.data(), ey.data());
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.35);
  g.SetMarkerColor(kBlack);
  g.SetLineColor(kBlack);
  SetGraphBase(&g, Form("Combined %s fit", cfg.titlePass.Data()), "#chi^{2}/ndf");
  const auto xlim = XLimits(cfg.xCenters);
  g.GetXaxis()->SetLimits(xlim.first, xlim.second);
  const double ymax = std::max(1.0, fit.chi2_ndf * 1.35);
  g.SetMinimum(0.0);
  g.SetMaximum(ymax);
  g.Draw("AP");
  c.SaveAs(outPath);
}

static void DrawParamPad(const PlotConfig& cfg,
                         const std::map<std::string, ParamErr>& errs,
                         const FitSummary& fit,
                         const std::vector<std::string>& names,
                         const std::vector<double>& xValues,
                         const char* title,
                         const char* ytitle) {
  const int n = static_cast<int>(names.size());
  std::vector<double> y(n), ey(n), ex(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const ParamErr pe = ErrorFor(errs, fit, names[i]);
    y[i] = pe.central;
    ey[i] = std::isfinite(pe.toy_std) ? pe.toy_std : 0.0;
  }

  double ymin = 0.0, ymax = 0.0;
  RangeFromValues(y, ey, ymin, ymax);

  TGraphErrors* g = new TGraphErrors(n, xValues.data(), y.data(), ex.data(), ey.data());
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);
  SetGraphBase(g, title, ytitle);
  const auto xlim = XLimits(cfg.xCenters);
  g->GetXaxis()->SetLimits(xlim.first, xlim.second);
  g->SetMinimum(ymin);
  g->SetMaximum(ymax);
  g->Draw("AP");
  g->Draw("P same");
}

static void PlotKinOffsets(const PlotConfig& cfg,
                           const FitSummary& fit,
                           const std::map<std::string, ParamErr>& errs,
                           const TString& outPath) {
  TCanvas c("c_kin_offsets_combined", "Combined kinematic offsets", 1400, 1000);
  c.Divide(2, 2);

  std::vector<std::string> commonDthp(1, cfg.dthpName);
  const std::vector<std::vector<std::string>> names = {cfg.dtheNames, cfg.dpeNames, commonDthp, cfg.dppNames};
  const std::vector<std::vector<double>> xValues = {
    cfg.settingParamXCenters,
    cfg.settingParamXCenters,
    cfg.commonParamXCenters,
    cfg.passTag == "4pass" ? cfg.settingParamXCenters : cfg.commonParamXCenters
  };
  const char* ytitles[4] = {"#Delta#theta_{e} (mrad)", "#Deltap_{e} (0.1%)", "#Delta#theta_{p} (mrad)", "#Deltap_{p} (0.1%)"};
  const char* labels[4] = {"#Delta#theta_{e}", "#Deltap_{e}", "#Delta#theta_{p}", "#Deltap_{p}"};

  for (int i = 0; i < 4; ++i) {
    c.cd(i + 1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.13);
    gPad->SetGrid();
    DrawParamPad(cfg, errs, fit, names[i],
                 xValues[i],
                 Form("Combined %s Offsets: %s", cfg.titlePass.Data(), labels[i]),
                 ytitles[i]);
  }
  c.SaveAs(outPath);
}

static void PlotResiduals(const PlotConfig& cfg, const std::vector<ResidualRow>& rows, const TString& outPath) {
  const std::vector<std::string> vars = {"W", "Em", "Pmz", "Pmx"};
  const std::map<std::string, int> colors = {
    {"W", kBlack},
    {"Em", kBlue + 1},
    {"Pmz", kRed + 1},
    {"Pmx", kGreen + 2}
  };
  const std::map<std::string, int> markers = {
    {"W", 20},
    {"Em", 21},
    {"Pmz", 22},
    {"Pmx", 33}
  };

  std::vector<double> allY, allE;
  for (const auto& r : rows) {
    allY.push_back(r.residual);
    allE.push_back(r.residual_err);
  }
  double ymin = 0.0, ymax = 0.0;
  RangeFromValues(allY, allE, ymin, ymax);

  TCanvas c("c_residuals_combined", "Combined residuals", 1200, 800);
  c.SetLeftMargin(0.11);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);
  c.SetGrid();

  TGraphErrors frame;
  frame.SetPoint(0, cfg.xCenters.front(), 0.0);
  frame.SetPoint(1, cfg.xCenters.back(), 0.0);
  SetGraphBase(&frame, Form("Combined %s residuals", cfg.titlePass.Data()), "Residual (MeV)");
  const auto xlim = XLimits(cfg.xCenters);
  frame.GetXaxis()->SetLimits(xlim.first, xlim.second);
  frame.SetMinimum(ymin);
  frame.SetMaximum(ymax);
  frame.Draw("AP");
  frame.GetXaxis()->SetNdivisions(507);

  TLegend leg(0.74, 0.72, 0.92, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.036);

  const std::map<std::string, double> xOffsets = {
    {"W", -0.09},
    {"Em", -0.03},
    {"Pmz", 0.03},
    {"Pmx", 0.09}
  };
  for (const auto& var : vars) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for (const auto& r : rows) {
      if (r.variable != var) continue;
      g->SetPoint(ip, r.delta_center + xOffsets.at(var), r.residual);
      g->SetPointError(ip, 0.0, r.residual_err);
      ++ip;
    }
    g->SetMarkerColor(colors.at(var));
    g->SetLineColor(colors.at(var));
    g->SetMarkerStyle(markers.at(var));
    g->SetMarkerSize(var == "Pmx" ? 1.7 : 1.2);
    g->SetLineWidth(2);
    g->Draw("P same");
    leg.AddEntry(g, var.c_str(), "p");
  }
  leg.Draw();
  c.SaveAs(outPath);
}

static void DrawResidualPad(const PlotConfig& cfg,
                            const std::vector<ResidualRow>& rows,
                            const std::string& setting,
                            const char* title,
                            const char* ytitle,
                            bool drawXTitle,
                            bool drawLegend,
                            double ymin,
                            double ymax) {
  const std::vector<std::string> vars = {"W", "Em", "Pmz", "Pmx"};
  const std::map<std::string, int> colors = {
    {"W", kBlack},
    {"Em", kBlue + 1},
    {"Pmz", kRed + 1},
    {"Pmx", kGreen + 2}
  };
  const std::map<std::string, int> markers = {
    {"W", 20},
    {"Em", 21},
    {"Pmz", 22},
    {"Pmx", 33}
  };
  const std::map<std::string, double> xOffsets = {
    {"W", -0.09},
    {"Em", -0.03},
    {"Pmz", 0.03},
    {"Pmx", 0.09}
  };

  TGraphErrors* frame = new TGraphErrors();
  frame->SetPoint(0, cfg.xCenters.front(), 0.0);
  frame->SetPoint(1, cfg.xCenters.back(), 0.0);
  SetGraphBase(frame, title, ytitle);
  const auto xlim = XLimits(cfg.xCenters);
  frame->GetXaxis()->SetLimits(xlim.first, xlim.second);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->SetNdivisions(507);
  if (!drawXTitle) {
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
  }
  frame->Draw("AP");

  TLegend* leg = nullptr;
  if (drawLegend) {
    leg = new TLegend(0.76, 0.66, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
  }

  for (const auto& var : vars) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for (const auto& r : rows) {
      if (r.setting != setting || r.variable != var) continue;
      g->SetPoint(ip, r.delta_center + xOffsets.at(var), r.residual);
      g->SetPointError(ip, 0.0, r.residual_err);
      ++ip;
    }
    g->SetMarkerColor(colors.at(var));
    g->SetLineColor(colors.at(var));
    g->SetMarkerStyle(markers.at(var));
    g->SetMarkerSize(var == "Pmx" ? 1.5 : 1.1);
    g->SetLineWidth(2);
    g->Draw("P same");
    if (leg) leg->AddEntry(g, var.c_str(), "p");
  }
  if (leg) leg->Draw();
}

static void PlotResiduals4Pass(const PlotConfig& cfg, const std::vector<ResidualRow>& rows, const TString& outPath) {
  std::vector<double> allY, allE;
  for (const auto& r : rows) {
    allY.push_back(r.residual);
    allE.push_back(r.residual_err);
  }
  double ymin = 0.0, ymax = 0.0;
  RangeFromValues(allY, allE, ymin, ymax);

  TCanvas c("c_residuals_4pass_combined", "Combined 4-pass residuals", 1200, 1000);
  c.Divide(1, 2);

  c.cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetGrid();
  DrawResidualPad(cfg, rows, "A", "Combined 4-pass residuals", "Residual_{A} (MeV)",
                  false, true, ymin, ymax);

  c.cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.16);
  gPad->SetTopMargin(0.04);
  gPad->SetGrid();
  DrawResidualPad(cfg, rows, "B", "", "Residual_{B} (MeV)",
                  true, false, ymin, ymax);

  c.SaveAs(outPath);
}

}  // namespace

void PlotSimultaneousDpeDthe(const char* passC, const char* outDirC = "results/PNGs/Combined") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const TString passTag = NormalizePassTag(passC);
  if (passTag.IsNull()) return;
  const PlotConfig cfg = MakeConfig(passTag);

  const TString fitCsv = Form("results/tables/SimultaneousDpeDtheFitParams_%s.csv", passTag.Data());
  const TString errCsv = Form("results/tables/SimultaneousDpeDtheErrors_%s_params.csv", passTag.Data());
  const TString residCsv = Form("results/tables/SimultaneousDpeDtheResiduals_%s.csv", passTag.Data());

  const TString outDir(outDirC);
  gSystem->mkdir(outDir, kTRUE);

  const FitSummary fit = ReadFitSummary(fitCsv);
  const auto errs = ReadParamErrors(errCsv);
  const auto residuals = ReadResiduals(residCsv);

  const TString chi2Path = outDir + Form("/chi2ndf_%s.png", passTag.Data());
  const TString kinPath = outDir + Form("/kinOffsets_%s.png", passTag.Data());
  const TString residPath = outDir + Form("/residuals_%s.png", passTag.Data());

  PlotChi2(cfg, fit, chi2Path);
  PlotKinOffsets(cfg, fit, errs, kinPath);
  if (passTag == "4pass") {
    PlotResiduals4Pass(cfg, residuals, residPath);
  } else {
    PlotResiduals(cfg, residuals, residPath);
  }

  std::cout << "[INFO] Wrote " << chi2Path << "\n";
  std::cout << "[INFO] Wrote " << kinPath << "\n";
  std::cout << "[INFO] Wrote " << residPath << "\n";
}
