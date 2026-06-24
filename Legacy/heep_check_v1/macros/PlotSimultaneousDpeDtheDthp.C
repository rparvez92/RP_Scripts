// PlotSimultaneousDpeDtheDthp.C
//
// Read a simultaneous-fit CSV produced by:
//   macros/SimultaneousDpeDtheDthp.C
//
// and produce, for each setting, three PNGs:
//
//   results/PNGs/PlotSimultaneousDpeDtheDthp/<Constrained,Unconstrained>/setting<A,B>_kinOffsets.png
//   results/PNGs/PlotSimultaneousDpeDtheDthp/<Constrained,Unconstrained>/setting<A,B>_residuals.png
//   results/PNGs/PlotSimultaneousDpeDtheDthp/<Constrained,Unconstrained>/setting<A,B>_chi2ndf.png
//
// Plot types:
//   1) kinematic offsets vs #delta bin center
//      - dthe (UL), dpe (UR), dthp (LL), dpp (LR)
//      - y-axis title: Offset (<unit>)
//      - narrow bins: black markers with error bars
//
//   2) residuals vs #delta bin center
//      - W (UL), Em (UR), Pmz (LL), Pmy (LR)
//      - y-axis title: Residuals (MeV)
//      - narrow bins: black markers with error bars
//
//   3) chi2/ndf vs #delta bin center
//      - single pad
//      - narrow bins: black markers
//
// Plotting policy:
//   - only narrow bins are plotted
//   - no full-bin representation is drawn
//   - no legend is drawn
//   - y-axis ranges are computed from plotted points and their errors
//   - rows with fit_valid != 1 are skipped
//   - residuals are interpreted as predicted - measured
//   - the CSV stores bin-specific dthe, dpe, and dthp columns
//     (for example dthe_b1_fit, dpe_b1_fit, dthp_b1_fit), which are mapped
//     here onto the plotted series row-by-row
//
// Default mode-specific paths:
//   constrained:
//     results/tables/SimultaneousDpeDtheDthp_constrained.csv
//     results/PNGs/PlotSimultaneousDpeDtheDthp/Constrained
//
//   unconstrained:
//     results/tables/SimultaneousDpeDtheDthp_unconstrained.csv
//     results/PNGs/PlotSimultaneousDpeDtheDthp/Unconstrained
//
// Example:
//   root -l -b -q 'macros/PlotSimultaneousDpeDtheDthp.C("constrained")'
//   root -l -b -q 'macros/PlotSimultaneousDpeDtheDthp.C("unconstrained")'
//   root -l -b -q 'macros/PlotSimultaneousDpeDtheDthp.C("constrained","results/tables/SimultaneousDpeDtheDthp_constrained.csv","results/PNGs/PlotSimultaneousDpeDtheDthp/Constrained")'

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Row {
  std::string setting;
  std::string mode;
  std::string bin;

  double dthe_fit = std::numeric_limits<double>::quiet_NaN();
  double dthe_fit_err = std::numeric_limits<double>::quiet_NaN();
  double dpe_fit = std::numeric_limits<double>::quiet_NaN();
  double dpe_fit_err = std::numeric_limits<double>::quiet_NaN();
  double dthp_fit = std::numeric_limits<double>::quiet_NaN();
  double dthp_fit_err = std::numeric_limits<double>::quiet_NaN();
  double dpp_fit = std::numeric_limits<double>::quiet_NaN();
  double dpp_fit_err = std::numeric_limits<double>::quiet_NaN();

  double dW_resid = std::numeric_limits<double>::quiet_NaN();
  double dW_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid_err = std::numeric_limits<double>::quiet_NaN();

  double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
  int fit_valid = 0;
};

struct SeriesDef {
  std::string key;
  std::string padTitle;
  std::string yTitle;
};

static std::string Trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  if (e > b && s[b] == '"' && s[e - 1] == '"') {
    ++b; --e;
  }
  return s.substr(b, e - b);
}

static std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

static std::string NormalizeMode(const std::string& s) {
  const std::string t = Lower(Trim(s));
  if (t == "constrained") return "constrained";
  if (t == "unconstrained" || t.empty()) return "unconstrained";
  return "";
}

static std::string NormalizeSetting(const std::string& s) {
  const std::string t = Lower(Trim(s));
  if (t == "a" || t == "setting a" || t == "settinga") return "A";
  if (t == "b" || t == "setting b" || t == "settingb") return "B";
  return Trim(s);
}

static std::string ModeDirLabel(const std::string& mode) {
  return mode == "constrained" ? "Constrained" : "Unconstrained";
}

static std::vector<std::string> SplitCsvLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inQuotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
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

static bool ToDouble(const std::string& s, double& v) {
  const std::string t = Trim(s);
  if (t.empty()) return false;
  char* end = nullptr;
  v = std::strtod(t.c_str(), &end);
  return end && *end == '\0';
}

static int ToIntOrDefault(const std::string& s, int defVal = 0) {
  const std::string t = Trim(s);
  if (t.empty()) return defVal;
  char* end = nullptr;
  long v = std::strtol(t.c_str(), &end, 10);
  if (!(end && *end == '\0')) return defVal;
  return static_cast<int>(v);
}

static int BinOrder(const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (b == "b1") return 1;
  if (b == "b2") return 2;
  if (b == "b3") return 3;
  if (b == "b4") return 4;
  if (b == "b5") return 5;
  return 99;
}

static double BinCenterFromLabel(const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (b == "b1") return -1.5;
  if (b == "b2") return -0.5;
  if (b == "b3") return  0.5;
  if (b == "b4") return  1.5;
  if (b == "b5") return  2.5;
  return std::numeric_limits<double>::quiet_NaN();
}

static double GetValue(const Row& r, const std::string& key) {
  if (key == "dthe") return r.dthe_fit;
  if (key == "dpe")  return r.dpe_fit;
  if (key == "dthp") return r.dthp_fit;
  if (key == "dpp")  return r.dpp_fit;
  if (key == "W")    return r.dW_resid;
  if (key == "Em")   return r.dEm_resid;
  if (key == "Pmz")  return r.dPmz_resid;
  if (key == "Pmy")  return r.dPmy_resid;
  if (key == "chi2ndf") return r.chi2_ndf;
  return std::numeric_limits<double>::quiet_NaN();
}

static double GetError(const Row& r, const std::string& key) {
  if (key == "dthe") return r.dthe_fit_err;
  if (key == "dpe")  return r.dpe_fit_err;
  if (key == "dthp") return r.dthp_fit_err;
  if (key == "dpp")  return r.dpp_fit_err;
  if (key == "W")    return r.dW_resid_err;
  if (key == "Em")   return r.dEm_resid_err;
  if (key == "Pmz")  return r.dPmz_resid_err;
  if (key == "Pmy")  return r.dPmy_resid_err;
  return 0.0;
}

static std::vector<Row> ReadSimultaneousCsv(const std::string& csvPath) {
  std::ifstream fin(csvPath);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input CSV: " << csvPath << std::endl;
    return {};
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) {
    std::cerr << "[ERROR] Empty CSV: " << csvPath << std::endl;
    return {};
  }

  const auto headers = SplitCsvLine(headerLine);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[Lower(headers[i])] = i;

  auto idx = [&](const std::string& name) -> int {
    auto it = col.find(Lower(name));
    return (it == col.end()) ? -1 : it->second;
  };

  const int iSetting   = idx("setting");
  const int iMode      = idx("mode");
  const int iBin       = idx("bin");
  const int iDpp       = idx("dpp_fit");
  const int iDppErr    = idx("dpp_fit_err");
  const int iWResid    = idx("dw_resid");
  const int iWResidErr = idx("dw_resid_err");
  const int iEmResid   = idx("dem_resid");
  const int iEmResidErr= idx("dem_resid_err");
  const int iPmzResid  = idx("dpmz_resid");
  const int iPmzResidErr = idx("dpmz_resid_err");
  const int iPmyResid  = idx("dpmy_resid");
  const int iPmyResidErr = idx("dpmy_resid_err");
  const int iChi2ndf   = idx("chi2_ndf");
  const int iFitValid  = idx("fit_valid");

  std::vector<std::string> missing;
  auto need = [&](int i, const std::string& name) {
    if (i < 0) missing.push_back(name);
  };

  need(iSetting, "setting");
  need(iMode, "mode");
  need(iBin, "bin");
  need(iDpp, "dpp_fit");
  need(iDppErr, "dpp_fit_err");
  need(iWResid, "dW_resid");
  need(iWResidErr, "dW_resid_err");
  need(iEmResid, "dEm_resid");
  need(iEmResidErr, "dEm_resid_err");
  need(iPmzResid, "dPmz_resid");
  need(iPmzResidErr, "dPmz_resid_err");
  need(iPmyResid, "dPmy_resid");
  need(iPmyResidErr, "dPmy_resid_err");
  need(iChi2ndf, "chi2_ndf");
  need(iFitValid, "fit_valid");

  if (!missing.empty()) {
    std::cerr << "[ERROR] Missing required columns in " << csvPath << ":\n";
    for (const auto& m : missing) std::cerr << "  - " << m << "\n";
    return {};
  }

  std::vector<Row> rows;
  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);

    auto getS = [&](int i) -> std::string {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    auto getD = [&](int i) -> double {
      double v = std::numeric_limits<double>::quiet_NaN();
      if (i >= 0 && i < static_cast<int>(fields.size())) ToDouble(fields[i], v);
      return v;
    };
    auto getI = [&](int i) -> int {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? ToIntOrDefault(fields[i], 0) : 0;
    };

    Row r;
    r.setting       = NormalizeSetting(getS(iSetting));
    r.mode          = NormalizeMode(getS(iMode));
    r.bin           = getS(iBin);
    r.dpp_fit       = getD(iDpp);
    r.dpp_fit_err   = getD(iDppErr);
    r.dW_resid      = getD(iWResid);
    r.dW_resid_err  = getD(iWResidErr);
    r.dEm_resid     = getD(iEmResid);
    r.dEm_resid_err = getD(iEmResidErr);
    r.dPmz_resid    = getD(iPmzResid);
    r.dPmz_resid_err= getD(iPmzResidErr);
    r.dPmy_resid    = getD(iPmyResid);
    r.dPmy_resid_err= getD(iPmyResidErr);
    r.chi2_ndf      = getD(iChi2ndf);
    r.fit_valid     = getI(iFitValid);

    const std::string binLower = Lower(Trim(r.bin));
    r.dthe_fit = getD(idx("dthe_" + binLower + "_fit"));
    r.dthe_fit_err = getD(idx("dthe_" + binLower + "_fit_err"));
    r.dpe_fit = getD(idx("dpe_" + binLower + "_fit"));
    r.dpe_fit_err = getD(idx("dpe_" + binLower + "_fit_err"));
    r.dthp_fit = getD(idx("dthp_" + binLower + "_fit"));
    r.dthp_fit_err = getD(idx("dthp_" + binLower + "_fit_err"));

    rows.push_back(r);
  }

  return rows;
}

static void ComputeYRange(const std::vector<Row>& rows,
                          const std::string& key,
                          bool useErrors,
                          double& yMin, double& yMax) {
  double mn = std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();

  for (const auto& r : rows) {
    const double y = GetValue(r, key);
    if (!std::isfinite(y)) continue;
    const double ey = useErrors ? GetError(r, key) : 0.0;
    const double lo = y - (std::isfinite(ey) ? ey : 0.0);
    const double hi = y + (std::isfinite(ey) ? ey : 0.0);
    mn = std::min(mn, lo);
    mx = std::max(mx, hi);
  }

  if (!std::isfinite(mn) || !std::isfinite(mx)) {
    yMin = -1.0;
    yMax = 1.0;
    return;
  }

  if (std::fabs(mx - mn) < 1e-12) {
    const double pad = (std::fabs(mx) > 0.0) ? 0.25 * std::fabs(mx) : 1.0;
    yMin = mn - pad;
    yMax = mx + pad;
    return;
  }

  const double span = mx - mn;
  const double pad = 0.18 * span;
  yMin = mn - pad;
  yMax = mx + pad;
}

static std::vector<Row> SortRows(const std::vector<Row>& rows) {
  std::vector<Row> out = rows;
  std::sort(out.begin(), out.end(),
            [](const Row& a, const Row& b) {
              return BinOrder(a.bin) < BinOrder(b.bin);
            });
  return out;
}

static void DrawSeriesPad(TPad* pad,
                          const std::vector<Row>& rows,
                          const SeriesDef& sdef,
                          bool useErrors) {
  pad->cd();
  pad->SetGrid(1,1);
  pad->SetTicks(1,1);
  pad->SetLeftMargin(0.14);
  pad->SetRightMargin(0.06);
  pad->SetBottomMargin(0.13);
  pad->SetTopMargin(0.08);

  std::vector<double> x, y, ex, ey;
  const auto sorted = SortRows(rows);

  for (const auto& r : sorted) {
    const double xc = BinCenterFromLabel(r.bin);
    const double yy = GetValue(r, sdef.key);
    const double ee = useErrors ? GetError(r, sdef.key) : 0.0;
    if (!std::isfinite(xc) || !std::isfinite(yy)) continue;
    x.push_back(xc);
    y.push_back(yy);
    ex.push_back(0.0);
    ey.push_back(std::isfinite(ee) ? ee : 0.0);
  }

  double yMin = 0.0, yMax = 0.0;
  ComputeYRange(sorted, sdef.key, useErrors, yMin, yMax);

  TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()));
  for (int i = 0; i < static_cast<int>(x.size()); ++i) {
    g->SetPoint(i, x[i], y[i]);
    g->SetPointError(i, ex[i], ey[i]);
  }

  g->SetTitle("");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.10);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);

  g->GetXaxis()->SetTitle("#delta bin center (%)");
  g->GetYaxis()->SetTitle(sdef.yTitle.c_str());
  g->GetXaxis()->SetLimits(-2.2, 3.2);
  g->GetYaxis()->SetRangeUser(yMin, yMax);
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  g->GetXaxis()->SetTitleSize(0.050);
  g->GetYaxis()->SetTitleSize(0.050);
  g->GetXaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetTitleOffset(1.2);
  g->Draw("AP");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.050);
  lat.DrawLatex(0.16, 0.92, sdef.padTitle.c_str());

  pad->Modified();
}

static void MakeOffsetsCanvas(const std::string& setting,
                              const std::vector<Row>& rows,
                              const std::string& outDir,
                              const std::string& modeLabel) {
  TCanvas* c = new TCanvas(Form("c_sim_dpdtp_kin_%s_%s", setting.c_str(), modeLabel.c_str()),
                           "SimultaneousDpeDtheDthp kinematic offsets", 1200, 900);
  c->Divide(2,2,0.002,0.002);

  const SeriesDef defs[4] = {
    {"dthe", "dthe", "Offset (mrad)"},
    {"dpe",  "dpe",  "Offset (0.1%)"},
    {"dthp", "dthp", "Offset (mrad)"},
    {"dpp",  "dpp",  "Offset (0.1%)"}
  };

  for (int i = 0; i < 4; ++i) {
    DrawSeriesPad((TPad*)c->cd(i + 1), rows, defs[i], true);
  }

  const std::string outName =
    outDir + "/setting" + setting + "_kinOffsets.png";
  c->SaveAs(outName.c_str());
  delete c;
}

static void MakeResidualsCanvas(const std::string& setting,
                                const std::vector<Row>& rows,
                                const std::string& outDir,
                                const std::string& modeLabel) {
  TCanvas* c = new TCanvas(Form("c_sim_dpdtp_res_%s_%s", setting.c_str(), modeLabel.c_str()),
                           "SimultaneousDpeDtheDthp residuals", 1200, 900);
  c->Divide(2,2,0.002,0.002);

  const SeriesDef defs[4] = {
    {"W",   "W",   "Residuals (MeV)"},
    {"Em",  "Em",  "Residuals (MeV)"},
    {"Pmz", "Pmz", "Residuals (MeV)"},
    {"Pmy", "Pmy", "Residuals (MeV)"}
  };

  for (int i = 0; i < 4; ++i) {
    DrawSeriesPad((TPad*)c->cd(i + 1), rows, defs[i], true);
  }

  const std::string outName =
    outDir + "/setting" + setting + "_residuals.png";
  c->SaveAs(outName.c_str());
  delete c;
}

static void MakeChi2Canvas(const std::string& setting,
                           const std::vector<Row>& rows,
                           const std::string& outDir,
                           const std::string& modeLabel) {
  TCanvas* c = new TCanvas(Form("c_sim_dpdtp_chi2_%s_%s", setting.c_str(), modeLabel.c_str()),
                           "SimultaneousDpeDtheDthp chi2/ndf", 900, 700);
  c->cd();
  gPad->SetGrid(1,1);
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.06);
  gPad->SetBottomMargin(0.13);
  gPad->SetTopMargin(0.08);

  std::vector<double> x, y;
  const auto sorted = SortRows(rows);
  for (const auto& r : sorted) {
    const double xc = BinCenterFromLabel(r.bin);
    const double yy = r.chi2_ndf;
    if (!std::isfinite(xc) || !std::isfinite(yy)) continue;
    x.push_back(xc);
    y.push_back(yy);
  }

  double yMin = 0.0, yMax = 0.0;
  ComputeYRange(sorted, "chi2ndf", false, yMin, yMax);

  TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()));
  for (int i = 0; i < static_cast<int>(x.size()); ++i) {
    g->SetPoint(i, x[i], y[i]);
    g->SetPointError(i, 0.0, 0.0);
  }

  g->SetTitle("");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.10);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetLineWidth(2);
  g->GetXaxis()->SetTitle("#delta bin center (%)");
  g->GetYaxis()->SetTitle("chi2/ndf");
  g->GetXaxis()->SetLimits(-2.2, 3.2);
  g->GetYaxis()->SetRangeUser(yMin, yMax);
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();
  g->GetXaxis()->SetTitleSize(0.050);
  g->GetYaxis()->SetTitleSize(0.050);
  g->GetXaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetLabelSize(0.042);
  g->GetYaxis()->SetTitleOffset(1.2);
  g->Draw("AP");

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.050);
  lat.DrawLatex(0.16, 0.92, "chi2/ndf");

  const std::string outName =
    outDir + "/setting" + setting + "_chi2ndf.png";
  c->SaveAs(outName.c_str());
  delete c;
}

} // namespace

void PlotSimultaneousDpeDtheDthp(const char* modeIn = "unconstrained",
                                 const char* inCsvIn = "",
                                 const char* outDirIn = "") {
  const std::string mode = NormalizeMode(modeIn ? modeIn : "");
  if (mode.empty()) {
    std::cerr << "[ERROR] Unknown mode. Use \"constrained\" or \"unconstrained\"." << std::endl;
    return;
  }

  const std::string modeDir = ModeDirLabel(mode);

  std::string inCsv = Trim(inCsvIn ? inCsvIn : "");
  if (inCsv.empty()) {
    inCsv = "results/tables/SimultaneousDpeDtheDthp_" + mode + ".csv";
  }

  std::string outDir = Trim(outDirIn ? outDirIn : "");
  if (outDir.empty()) {
    outDir = "results/PNGs/PlotSimultaneousDpeDtheDthp/" + modeDir;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir.c_str(), true);

  const auto allRows = ReadSimultaneousCsv(inCsv);
  if (allRows.empty()) {
    std::cerr << "[ERROR] No rows read from " << inCsv << std::endl;
    return;
  }

  std::map<std::string, std::vector<Row>> bySetting;
  int skippedInvalid = 0;
  for (const auto& r : allRows) {
    if (r.fit_valid != 1) {
      ++skippedInvalid;
      continue;
    }
    if (r.mode != mode) continue;
    bySetting[r.setting].push_back(r);
  }

  if (bySetting.empty()) {
    std::cerr << "[ERROR] No valid rows found for mode " << mode << " in " << inCsv << std::endl;
    return;
  }

  if (bySetting.count("A")) {
    MakeOffsetsCanvas("A", bySetting["A"], outDir, modeDir);
    MakeResidualsCanvas("A", bySetting["A"], outDir, modeDir);
    MakeChi2Canvas("A", bySetting["A"], outDir, modeDir);
  }
  if (bySetting.count("B")) {
    MakeOffsetsCanvas("B", bySetting["B"], outDir, modeDir);
    MakeResidualsCanvas("B", bySetting["B"], outDir, modeDir);
    MakeChi2Canvas("B", bySetting["B"], outDir, modeDir);
  }

  for (const auto& kv : bySetting) {
    if (kv.first == "A" || kv.first == "B") continue;
    MakeOffsetsCanvas(kv.first, kv.second, outDir, modeDir);
    MakeResidualsCanvas(kv.first, kv.second, outDir, modeDir);
    MakeChi2Canvas(kv.first, kv.second, outDir, modeDir);
  }

  std::cout << "[INFO] Wrote PlotSimultaneousDpeDtheDthp PNGs to " << outDir
            << " (skipped invalid rows: " << skippedInvalid << ")" << std::endl;
}
