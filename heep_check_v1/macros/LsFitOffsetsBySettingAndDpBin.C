// LsFitOffsetsBySettingAndDpBin.C
//
// Overconstrained weighted least-squares / chi2 minimization
// for HEEP kinematic offsets.
//
// For each (setting, dp_idx):
//   fit one common parameter vector
//
//      x = [ dthe, dpe, dthp, dpp ]^T
//
// using run-wise measured offsets from:
//   heep_check_dpBinnedOffsets.csv
//
// Observables:
//   y = [ W, Em, Pmz, Pmy ]
//
// Response matrices:
//   Setting A -> SHMS angle 12.465
//   Setting B -> SHMS angle 19.34
//
// Assumed run mapping:
//   Setting B runs = {23849, 23850, 23851}
//   all other runs in this CSV = Setting A
//
// Outputs:
//   1) LsFitOffsetsBySettingAndDpBin_summary.csv
//   2) LsFitOffsetsBySettingAndDpBin_residuals.csv
//   3) residual PNGs under results/PNGs/LsFitOffsetsBySettingAndDpBin/
//
// Usage:
//   root -l -b -q 'macros/LsFitOffsetsBySettingAndDpBin.C("results/tables/heep_check_dpBinnedOffsets.csv")'
//

#include <TCanvas.h>
#include <TDecompLU.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

namespace {

struct Row {
  int run = -1;
  int dp_idx = -999;
  std::string dp_label;
  std::string var;       // W, Em, Pmz, Pmy
  double dmu_MeV = NAN;
  double dmu_err_MeV = NAN;
  int fit_valid = 0;
};

struct FitResult {
  bool ok = false;
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

  double predW = NAN, predEm = NAN, predPmz = NAN, predPmy = NAN;
};

struct ResidRow {
  std::string setting;
  int dp_idx = -999;
  std::string dp_label;
  int run = -1;
  std::string var;
  double y_obs = NAN;
  double y_err = NAN;
  double y_pred = NAN;
  double resid = NAN;
  double pull = NAN;
};

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

std::string Trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
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

bool ReadDpBinnedCSV(const std::string& csvPath, std::vector<Row>& rows) {
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
  std::map<std::string,int> col;
  for (int i = 0; i < (int)headers.size(); ++i) {
    col[Trim(headers[i])] = i;
  }

  const std::vector<std::string> need = {
    "run", "dp_idx", "dp_label", "var", "dmu_MeV", "dmu_err_MeV", "fit_valid"
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

    Row r;
    r.run         = ToInt(f[col["run"]]);
    r.dp_idx      = ToInt(f[col["dp_idx"]]);
    r.dp_label    = Trim(f[col["dp_label"]]);
    r.var         = Trim(f[col["var"]]);
    r.dmu_MeV     = ToDouble(f[col["dmu_MeV"]]);
    r.dmu_err_MeV = ToDouble(f[col["dmu_err_MeV"]]);
    r.fit_valid   = ToInt(f[col["fit_valid"]]);
    rows.push_back(r);
  }

  return true;
}

std::string SettingFromRun(int run) {
  static const std::set<int> runsB = {23849, 23850, 23851};
  return runsB.count(run) ? "B" : "A";
}

TMatrixD BuildResponseMatrix(const std::string& setting) {
  TMatrixD A(4,4);

  // rows: W, Em, Pmz, Pmy
  // cols: dthe, dpe, dthp, dpp

  if (setting == "A") {
    A(0,0) = -14.08; A(0,1) = -8.62; A(0,2) =  0.00; A(0,3) =  0.00;
    A(1,0) =   0.00; A(1,1) = -7.06; A(1,2) =  0.00; A(1,3) = -2.10;
    A(2,0) =  -5.75; A(2,1) =  4.10; A(2,2) =  0.00; A(2,3) =  2.27;
    A(3,0) =   4.10; A(3,1) =  5.75; A(3,2) = -2.27; A(3,3) =  0.00;
  } else if (setting == "B") {
    A(0,0) = -17.33; A(0,1) = -8.62; A(0,2) =  0.00; A(0,3) =  0.00;
    A(1,0) =   0.00; A(1,1) = -5.66; A(1,2) =  0.00; A(1,3) = -3.63;
    A(2,0) =  -4.30; A(2,1) =  3.69; A(2,2) =  0.00; A(2,3) =  3.74;
    A(3,0) =   3.69; A(3,1) =  4.30; A(3,2) = -3.74; A(3,3) =  0.00;
  } else {
    std::cerr << "[ERROR] Unknown setting: " << setting << std::endl;
  }

  return A;
}

int VarIndex(const std::string& var) {
  if (var == "W")   return 0;
  if (var == "Em")  return 1;
  if (var == "Pmz") return 2;
  if (var == "Pmy") return 3;
  return -1;
}

double PredictForVar(const TVectorD& x, const TMatrixD& A4, const std::string& var) {
  int i = VarIndex(var);
  if (i < 0) return NAN;
  double y = 0.0;
  for (int j = 0; j < 4; ++j) y += A4(i,j) * x(j);
  return y;
}

bool SolveWeightedLS(const std::vector<Row>& rows,
                     const std::string& setting,
                     FitResult& fitOut,
                     std::vector<ResidRow>& residOut)
{
  if (rows.size() < 4) {
    std::cerr << "[WARN] Too few rows to fit setting " << setting << std::endl;
    return false;
  }

  const int M = (int)rows.size();
  TMatrixD A(M,4);
  TVectorD y(M);
  TVectorD sig(M);

  TMatrixD A4 = BuildResponseMatrix(setting);

  std::set<int> uniqueRuns;
  std::string dpLabel = rows.front().dp_label;
  int dpIdx = rows.front().dp_idx;

  for (int i = 0; i < M; ++i) {
    const auto& r = rows[i];
    uniqueRuns.insert(r.run);

    const int obs = VarIndex(r.var);
    if (obs < 0) {
      std::cerr << "[WARN] Unknown var = " << r.var << std::endl;
      return false;
    }

    for (int j = 0; j < 4; ++j) A(i,j) = A4(obs,j);
    y(i) = r.dmu_MeV;
    sig(i) = r.dmu_err_MeV;

    if (!(sig(i) > 0.0) || !std::isfinite(sig(i))) {
      std::cerr << "[WARN] Bad sigma for run " << r.run
                << " var " << r.var << " sigma=" << sig(i) << std::endl;
      return false;
    }
  }

  TMatrixD W(M,M);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) W(i,j) = 0.0;
    W(i,i) = 1.0 / (sig(i) * sig(i));
  }

  TMatrixD AT(TMatrixD::kTransposed, A);
  TMatrixD N = AT * W * A;   // normal matrix = A^T W A
  TMatrixD Cov = N;

  TDecompLU lu(Cov);
  Bool_t ok = lu.Invert(Cov);
  if (!ok) {
    std::cerr << "[WARN] Could not invert normal matrix for "
              << setting << " dp_idx=" << dpIdx << std::endl;
    return false;
  }

  TVectorD rhs = AT * (W * y);
  TVectorD x = Cov * rhs;

  TVectorD yfit = A * x;
  TVectorD resid = y - yfit;

  double chi2 = 0.0;
  for (int i = 0; i < M; ++i) {
    chi2 += resid(i)*resid(i)/(sig(i)*sig(i));
  }
  int ndf = M - 4;

  fitOut.ok = true;
  fitOut.setting = setting;
  fitOut.dp_idx = dpIdx;
  fitOut.dp_label = dpLabel;
  fitOut.n_rows = M;
  fitOut.n_runs = (int)uniqueRuns.size();

  fitOut.dthe = x(0); fitOut.dthe_err = std::sqrt(std::fabs(Cov(0,0)));
  fitOut.dpe  = x(1); fitOut.dpe_err  = std::sqrt(std::fabs(Cov(1,1)));
  fitOut.dthp = x(2); fitOut.dthp_err = std::sqrt(std::fabs(Cov(2,2)));
  fitOut.dpp  = x(3); fitOut.dpp_err  = std::sqrt(std::fabs(Cov(3,3)));

  fitOut.chi2 = chi2;
  fitOut.ndf = ndf;

  fitOut.predW   = PredictForVar(x, A4, "W");
  fitOut.predEm  = PredictForVar(x, A4, "Em");
  fitOut.predPmz = PredictForVar(x, A4, "Pmz");
  fitOut.predPmy = PredictForVar(x, A4, "Pmy");

  for (int i = 0; i < M; ++i) {
    ResidRow rr;
    rr.setting = setting;
    rr.dp_idx = dpIdx;
    rr.dp_label = dpLabel;
    rr.run = rows[i].run;
    rr.var = rows[i].var;
    rr.y_obs = y(i);
    rr.y_err = sig(i);
    rr.y_pred = yfit(i);
    rr.resid = resid(i);
    rr.pull = resid(i) / sig(i);
    residOut.push_back(rr);
  }

  return true;
}

void MakeResidualPlot(const std::vector<ResidRow>& residRows,
                      const std::string& outPng)
{
  if (residRows.empty()) return;

  std::vector<ResidRow> Wv, Emv, Pmzv, Pmyv;
  for (const auto& r : residRows) {
    if (r.var == "W") Wv.push_back(r);
    else if (r.var == "Em") Emv.push_back(r);
    else if (r.var == "Pmz") Pmzv.push_back(r);
    else if (r.var == "Pmy") Pmyv.push_back(r);
  }

  auto sorter = [](const ResidRow& a, const ResidRow& b){ return a.run < b.run; };
  std::sort(Wv.begin(),   Wv.end(), sorter);
  std::sort(Emv.begin(),  Emv.end(), sorter);
  std::sort(Pmzv.begin(), Pmzv.end(), sorter);
  std::sort(Pmyv.begin(), Pmyv.end(), sorter);

  auto makeGraph = [](const std::vector<ResidRow>& v, const char* title) {
    TGraphErrors* g = new TGraphErrors((int)v.size());
    for (int i = 0; i < (int)v.size(); ++i) {
      g->SetPoint(i, v[i].run, v[i].resid);
      g->SetPointError(i, 0.0, v[i].y_err);
    }
    g->SetTitle(title);
    g->SetMarkerStyle(20);
    g->GetXaxis()->SetTitle("Run");
    g->GetYaxis()->SetTitle("Residual (MeV)");
    return g;
  };

  TCanvas* c = new TCanvas("c_resid","Residuals",1400,900);
  c->Divide(2,2);

  struct Block { std::vector<ResidRow>* v; const char* title; };
  Block blocks[4] = {
    {&Wv,   "W residuals"},
    {&Emv,  "Em residuals"},
    {&Pmzv, "Pmz residuals"},
    {&Pmyv, "Pmy residuals"}
  };

  for (int ipad = 0; ipad < 4; ++ipad) {
    c->cd(ipad+1);
    gPad->SetGrid();

    auto* v = blocks[ipad].v;
    TGraphErrors* g = makeGraph(*v, blocks[ipad].title);

    double ymin =  1e9;
    double ymax = -1e9;
    double xmin =  1e9;
    double xmax = -1e9;

    for (const auto& r : *v) {
      ymin = std::min(ymin, r.resid - r.y_err);
      ymax = std::max(ymax, r.resid + r.y_err);
      xmin = std::min(xmin, (double)r.run);
      xmax = std::max(xmax, (double)r.run);
    }
    if (!v->empty()) {
      double pad = 0.15 * std::max(1.0, ymax - ymin);
      g->GetYaxis()->SetRangeUser(ymin - pad, ymax + pad);
      g->GetXaxis()->SetLimits(xmin - 0.8, xmax + 0.8);
    }

    g->Draw("AP");

    TLine* z = new TLine(xmin - 1.0, 0.0, xmax + 1.0, 0.0);
    z->SetLineStyle(2);
    z->Draw("same");
  }

  c->SaveAs(outPng.c_str());
  delete c;
}

} // namespace


void LsFitOffsetsBySettingAndDpBin(const char* csvPath = "heep_check_dpBinnedOffsets.csv") {
  gStyle->SetOptStat(0);

  std::vector<Row> allRows;
  if (!ReadDpBinnedCSV(csvPath, allRows)) return;

  std::map<std::pair<std::string,int>, std::vector<Row>> groups;

  // Build groups by (setting, dp_idx), using only valid rows
  for (const auto& r : allRows) {
    if (r.fit_valid != 1) continue;
    if (!(r.var == "W" || r.var == "Em" || r.var == "Pmz" || r.var == "Pmy")) continue;

    std::string setting = SettingFromRun(r.run);
    groups[{setting, r.dp_idx}].push_back(r);
  }

  fs::create_directories("results/PNGs/LsFitOffsetsBySettingAndDpBin");

  std::ofstream fsum("results/tables/LsFitOffsetsBySettingAndDpBin_summary.csv");
  fsum << "setting,dp_idx,dp_label,n_rows,n_runs,"
       << "dthe,dthe_err,dpe,dpe_err,dthp,dthp_err,dpp,dpp_err,"
       << "chi2,ndf,chi2_per_ndf,"
       << "predW,predEm,predPmz,predPmy\n";

  std::ofstream fres("results/tables/LsFitOffsetsBySettingAndDpBin_residuals.csv");
  fres << "setting,dp_idx,dp_label,run,var,y_obs,y_err,y_pred,resid,pull\n";

  for (const auto& kv : groups) {
    const std::string setting = kv.first.first;
    const int dp_idx = kv.first.second;
    const auto& rows = kv.second;

    FitResult fit;
    std::vector<ResidRow> residRows;
    bool ok = SolveWeightedLS(rows, setting, fit, residRows);

    if (!ok) {
      std::cerr << "[WARN] Fit failed for setting=" << setting
                << " dp_idx=" << dp_idx << std::endl;
      continue;
    }

    const double chi2ndf = (fit.ndf > 0) ? (fit.chi2 / fit.ndf) : NAN;

    std::cout << "\n====================================================\n";
    std::cout << "setting = " << fit.setting
              << "   dp_idx = " << fit.dp_idx
              << "   dp_label = " << fit.dp_label << "\n";
    std::cout << "n_rows = " << fit.n_rows
              << "   n_runs = " << fit.n_runs << "\n";
    std::cout << "dthe = " << fit.dthe << " +/- " << fit.dthe_err << "  [mrad]\n";
    std::cout << "dpe  = " << fit.dpe  << " +/- " << fit.dpe_err  << "  [0.1%]\n";
    std::cout << "dthp = " << fit.dthp << " +/- " << fit.dthp_err << "  [mrad]\n";
    std::cout << "dpp  = " << fit.dpp  << " +/- " << fit.dpp_err  << "  [0.1%]\n";
    std::cout << "chi2/ndf = " << fit.chi2 << " / " << fit.ndf;
    if (fit.ndf > 0) std::cout << " = " << chi2ndf;
    std::cout << "\n";
    std::cout << "predicted offsets: "
              << "W="   << fit.predW
              << ", Em="  << fit.predEm
              << ", Pmz=" << fit.predPmz
              << ", Pmy=" << fit.predPmy << "\n";

    fsum << fit.setting << ","
         << fit.dp_idx << ","
         << fit.dp_label << ","
         << fit.n_rows << ","
         << fit.n_runs << ","
         << fit.dthe << "," << fit.dthe_err << ","
         << fit.dpe  << "," << fit.dpe_err  << ","
         << fit.dthp << "," << fit.dthp_err << ","
         << fit.dpp  << "," << fit.dpp_err  << ","
         << fit.chi2 << "," << fit.ndf << "," << chi2ndf << ","
         << fit.predW << "," << fit.predEm << ","
         << fit.predPmz << "," << fit.predPmy << "\n";

    for (const auto& rr : residRows) {
      fres << rr.setting << ","
           << rr.dp_idx << ","
           << rr.dp_label << ","
           << rr.run << ","
           << rr.var << ","
           << rr.y_obs << ","
           << rr.y_err << ","
           << rr.y_pred << ","
           << rr.resid << ","
           << rr.pull << "\n";
    }

    std::ostringstream png;
    png << "results/PNGs/LsFitOffsetsBySettingAndDpBin/residuals_setting"
        << fit.setting << "_dp" << fit.dp_idx
        << "_" << fit.dp_label << ".png";
    MakeResidualPlot(residRows, png.str());
  }

  fsum.close();
  fres.close();

  std::cout << "\n[INFO] Wrote: results/tables/LsFitOffsetsBySettingAndDpBin_summary.csv\n";
  std::cout << "[INFO] Wrote: results/tables/LsFitOffsetsBySettingAndDpBin_residuals.csv\n";
  std::cout << "[INFO] Residual PNGs in results/PNGs/LsFitOffsetsBySettingAndDpBin/\n";
}
