// PlotGuessedOffsets.C
//
// Read guessed kinematic offsets, propagated prediction uncertainties,
// and residuals from:
//   results/tables/GuessedOffsets.csv
//
// and produce, for each setting, two PNGs with 4 pads each:
//
//   1) Guessed kinematic offsets
//      UL: dthe_guess  (mrad)
//      UR: dpe_guess   (0.1%)
//      LL: dthp_guess  (mrad)
//      LR: dpp_guess   (0.1%)
//
//   2) Residuals vs delta-bin center
//      UL: dW_resid    (MeV)
//      UR: dEm_resid   (MeV)
//      LL: dPmz_resid  (MeV)
//      LR: dPmy_resid  (MeV)
//
// Design choices:
//  - Guessed-offset pads draw narrow bins (b1..b5) as black markers with vertical
//    error bars.
//  - Residual pads draw narrow bins as black markers with vertical error bars,
//    using the propagated residual errors from the CSV.
//  - Prediction-uncertainty columns also exist in the CSV for downstream use,
//    but they are not plotted in this macro.
//  - Y-axis ranges are computed from the plotted content so they stay reasonably tight.
//  - No legend is drawn, since each pad shows only one class of data points.
//  - Full-bin rows remain in the CSV but are not plotted.
//  - The expected narrow-bin centers are for five equal bins in [-2,3]:
//      b1=-1.5, b2=-0.5, b3=0.5, b4=1.5, b5=2.5
//
// Example:
//   root -l -b -q 'macros/PlotGuessedOffsets.C("results/tables/GuessedOffsets.csv","results/PNGs/PlotGuessedOffsets")'

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TBox.h>
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
  std::string bin;

  double dthe_guess     = std::numeric_limits<double>::quiet_NaN();
  double dthe_guess_err = std::numeric_limits<double>::quiet_NaN();

  double dpe_guess      = std::numeric_limits<double>::quiet_NaN();
  double dpe_guess_err  = std::numeric_limits<double>::quiet_NaN();

  double dthp_guess     = std::numeric_limits<double>::quiet_NaN();
  double dthp_guess_err = std::numeric_limits<double>::quiet_NaN();

  double dpp_guess      = std::numeric_limits<double>::quiet_NaN();
  double dpp_guess_err  = std::numeric_limits<double>::quiet_NaN();

  double dW_resid       = std::numeric_limits<double>::quiet_NaN();
  double dW_resid_err   = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid      = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid_err  = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid     = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid     = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid_err = std::numeric_limits<double>::quiet_NaN();
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

static bool ToDouble(const std::string& s, double& v) {
  const std::string t = Trim(s);
  if (t.empty()) return false;
  char* end = nullptr;
  v = std::strtod(t.c_str(), &end);
  return end && *end == '\0';
}

static std::string Lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

static double BinCenterFromLabel(const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (b == "b1") return -1.5;
  if (b == "b2") return -0.5;
  if (b == "b3") return  0.5;
  if (b == "b4") return  1.5;
  if (b == "b5") return  2.5;
  return std::numeric_limits<double>::quiet_NaN(); // full has no center here
}

static int BinOrder(const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (b == "full") return 0;
  if (b == "b1")   return 1;
  if (b == "b2")   return 2;
  if (b == "b3")   return 3;
  if (b == "b4")   return 4;
  if (b == "b5")   return 5;
  return 99;
}

static std::string NormalizeSetting(const std::string& s) {
  const std::string t = Lower(Trim(s));
  if (t == "a" || t == "setting a" || t == "settinga") return "A";
  if (t == "b" || t == "setting b" || t == "settingb") return "B";
  return Trim(s);
}

static double GetValue(const Row& r, const std::string& key) {
  if (key == "dthe") return r.dthe_guess;
  if (key == "dpe")  return r.dpe_guess;
  if (key == "dthp") return r.dthp_guess;
  if (key == "dpp")  return r.dpp_guess;
  if (key == "W")    return r.dW_resid;
  if (key == "Em")   return r.dEm_resid;
  if (key == "Pmz")  return r.dPmz_resid;
  if (key == "Pmy")  return r.dPmy_resid;
  return std::numeric_limits<double>::quiet_NaN();
}

static double GetError(const Row& r, const std::string& key) {
  if (key == "dthe") return r.dthe_guess_err;
  if (key == "dpe")  return r.dpe_guess_err;
  if (key == "dthp") return r.dthp_guess_err;
  if (key == "dpp")  return r.dpp_guess_err;
  if (key == "W")    return r.dW_resid_err;
  if (key == "Em")   return r.dEm_resid_err;
  if (key == "Pmz")  return r.dPmz_resid_err;
  if (key == "Pmy")  return r.dPmy_resid_err;
  return std::numeric_limits<double>::quiet_NaN();
}

static std::vector<Row> ReadGuessedOffsets(const std::string& csvPath) {
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
  for (int i = 0; i < (int)headers.size(); ++i) {
    col[Lower(headers[i])] = i;
  }

  auto idx = [&](const std::string& name) -> int {
    auto it = col.find(Lower(name));
    return (it == col.end()) ? -1 : it->second;
  };

  const int iSetting       = idx("setting");
  const int iBin           = idx("bin");
  const int iDthe          = idx("dthe_guess");
  const int iDtheErr       = idx("dthe_guess_err");
  const int iDpe           = idx("dpe_guess");
  const int iDpeErr        = idx("dpe_guess_err");
  const int iDthp          = idx("dthp_guess");
  const int iDthpErr       = idx("dthp_guess_err");
  const int iDpp           = idx("dpp_guess");
  const int iDppErr        = idx("dpp_guess_err");
  const int iWResid        = idx("dw_resid");
  const int iWResidErr     = idx("dw_resid_err");
  const int iEmResid       = idx("dem_resid");
  const int iEmResidErr    = idx("dem_resid_err");
  const int iPmzResid      = idx("dpmz_resid");
  const int iPmzResidErr   = idx("dpmz_resid_err");
  const int iPmyResid      = idx("dpmy_resid");
  const int iPmyResidErr   = idx("dpmy_resid_err");

  if (iSetting < 0 || iBin < 0 ||
      iDthe < 0 || iDtheErr < 0 ||
      iDpe  < 0 || iDpeErr  < 0 ||
      iDthp < 0 || iDthpErr < 0 ||
      iDpp  < 0 || iDppErr  < 0 ||
      iWResid < 0 || iWResidErr < 0 ||
      iEmResid < 0 || iEmResidErr < 0 ||
      iPmzResid < 0 || iPmzResidErr < 0 ||
      iPmyResid < 0 || iPmyResidErr < 0) {
    std::cerr << "[ERROR] Missing one or more required columns in " << csvPath << std::endl;
    std::cerr << "Required: setting, bin, "
              << "dthe_guess, dthe_guess_err, dpe_guess, dpe_guess_err, "
              << "dthp_guess, dthp_guess_err, dpp_guess, dpp_guess_err, "
              << "dW_resid, dW_resid_err, dEm_resid, dEm_resid_err, "
              << "dPmz_resid, dPmz_resid_err, dPmy_resid, dPmy_resid_err"
              << std::endl;
    return {};
  }

  std::vector<Row> rows;
  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);

    auto getS = [&](int i) -> std::string {
      return (i >= 0 && i < (int)fields.size()) ? Trim(fields[i]) : "";
    };
    auto getD = [&](int i) -> double {
      double v = std::numeric_limits<double>::quiet_NaN();
      if (i >= 0 && i < (int)fields.size()) ToDouble(fields[i], v);
      return v;
    };

    Row r;
    r.setting        = NormalizeSetting(getS(iSetting));
    r.bin            = getS(iBin);
    r.dthe_guess     = getD(iDthe);
    r.dthe_guess_err = getD(iDtheErr);
    r.dpe_guess      = getD(iDpe);
    r.dpe_guess_err  = getD(iDpeErr);
    r.dthp_guess     = getD(iDthp);
    r.dthp_guess_err = getD(iDthpErr);
    r.dpp_guess      = getD(iDpp);
    r.dpp_guess_err  = getD(iDppErr);
    r.dW_resid       = getD(iWResid);
    r.dW_resid_err   = getD(iWResidErr);
    r.dEm_resid      = getD(iEmResid);
    r.dEm_resid_err  = getD(iEmResidErr);
    r.dPmz_resid     = getD(iPmzResid);
    r.dPmz_resid_err = getD(iPmzResidErr);
    r.dPmy_resid     = getD(iPmyResid);
    r.dPmy_resid_err = getD(iPmyResidErr);

    rows.push_back(r);
  }

  return rows;
}

static void ComputeYRange(const std::vector<Row>& rows,
                          const std::string& key,
                          bool useErrors,
                          double& yMin, double& yMax) {
  double mn =  std::numeric_limits<double>::infinity();
  double mx = -std::numeric_limits<double>::infinity();

  for (const auto& r : rows) {
    const std::string b = Lower(Trim(r.bin));
    if (b == "full") continue;

    const double y  = GetValue(r, key);
    const double ey = useErrors ? GetError(r, key) : 0.0;
    if (!std::isfinite(y)) continue;

    const double lo = y - (std::isfinite(ey) ? ey : 0.0);
    const double hi = y + (std::isfinite(ey) ? ey : 0.0);

    mn = std::min(mn, lo);
    mx = std::max(mx, hi);
  }

  if (!std::isfinite(mn) || !std::isfinite(mx)) {
    yMin = -1.0;
    yMax =  1.0;
    return;
  }

  if (std::fabs(mx - mn) < 1e-12) {
    const double pad = (std::fabs(mx) > 0.0) ? 0.25 * std::fabs(mx) : 1.0;
    yMin = mn - pad;
    yMax = mx + pad;
    return;
  }

  const double span = mx - mn;
  const double pad  = 0.18 * span;
  yMin = mn - pad;
  yMax = mx + pad;
}

static void DrawOnePad(TPad* pad,
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

  std::vector<Row> sorted = rows;
  std::sort(sorted.begin(), sorted.end(),
            [](const Row& a, const Row& b) {
              return BinOrder(a.bin) < BinOrder(b.bin);
            });

  for (const auto& r : sorted) {
    const std::string b = Lower(Trim(r.bin));
    if (b == "full") continue;

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

  TGraphErrors* g = new TGraphErrors((int)x.size());
  for (int i = 0; i < (int)x.size(); ++i) {
    g->SetPoint(i, x[i], y[i]);
    g->SetPointError(i, ex[i], ey[i]);
  }

  g->SetTitle("");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.15);
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

static void MakeSettingCanvas(const std::string& setting,
                              const std::vector<Row>& rows,
                              const std::string& outDir) {
  TCanvas* c = new TCanvas(Form("c_%s", setting.c_str()),
                           Form("Guessed offsets - Setting %s", setting.c_str()),
                           1200, 900);
  c->Divide(2,2,0.002,0.002);

  const SeriesDef defs[4] = {
    {"dthe", "dthe", "guess offsets (mrad)"},
    {"dpe",  "dpe",  "guess offsets (0.1%)"},
    {"dthp", "dthp", "guess offsets (mrad)"},
    {"dpp",  "dpp",  "guess offsets (0.1%)"}
  };

  for (int i = 0; i < 4; ++i) {
    DrawOnePad((TPad*)c->cd(i + 1), rows, defs[i], true);
  }

  c->SaveAs((outDir + "/GuessedOffsets_Setting" + setting + ".png").c_str());
  delete c;
}

static void MakeResidualCanvas(const std::string& setting,
                               const std::vector<Row>& rows,
                               const std::string& outDir) {
  TCanvas* c = new TCanvas(Form("c_res_%s", setting.c_str()),
                           Form("Residuals - Setting %s", setting.c_str()),
                           1200, 900);
  c->Divide(2,2,0.002,0.002);

  const SeriesDef defs[4] = {
    {"W",   "W",   "Residuals (MeV)"},
    {"Em",  "Em",  "Residuals (MeV)"},
    {"Pmz", "Pmz", "Residuals (MeV)"},
    {"Pmy", "Pmy", "Residuals (MeV)"}
  };

  for (int i = 0; i < 4; ++i) {
    DrawOnePad((TPad*)c->cd(i + 1), rows, defs[i], true);
  }

  c->SaveAs((outDir + "/GuessedResiduals_Setting" + setting + ".png").c_str());
  delete c;
}

} // end anonymous namespace

void PlotGuessedOffsets(const char* inCsv  = "results/tables/GuessedOffsets.csv",
                        const char* outDir = "results/PNGs/PlotGuessedOffsets") {
  gStyle->SetOptStat(0);  
  gSystem->mkdir(outDir, kTRUE);

  const std::vector<Row> rows = ReadGuessedOffsets(inCsv);
  if (rows.empty()) {
    std::cerr << "[ERROR] No rows read from " << inCsv << std::endl;
    return;
  }

  std::map<std::string, std::vector<Row>> bySetting;
  for (const auto& r : rows) {
    bySetting[r.setting].push_back(r);
  }

  if (bySetting.empty()) {
    std::cerr << "[ERROR] No settings found in input." << std::endl;
    return;
  }

  if (bySetting.count("A")) MakeSettingCanvas("A", bySetting["A"], outDir);
  if (bySetting.count("B")) MakeSettingCanvas("B", bySetting["B"], outDir);
  if (bySetting.count("A")) MakeResidualCanvas("A", bySetting["A"], outDir);
  if (bySetting.count("B")) MakeResidualCanvas("B", bySetting["B"], outDir);

  // In case settings are not normalized exactly as A/B, also handle remaining keys.
  for (const auto& kv : bySetting) {
    if (kv.first == "A" || kv.first == "B") continue;
    MakeSettingCanvas(kv.first, kv.second, outDir);
    MakeResidualCanvas(kv.first, kv.second, outDir);
  }

  std::cout << "[INFO] Wrote guessed-offset and residual PNGs to " << outDir << std::endl;
}
