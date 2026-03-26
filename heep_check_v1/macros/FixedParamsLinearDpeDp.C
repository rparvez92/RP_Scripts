// FixedParamsLinearDpeDp.C
//
// XEM2-style HEEP check macro:
//   * Fix (dthe, dthp, dpp) to anchor values per setting
//   * For each narrow dp bin, solve ONE independent floating dpe
//   * Use chi2 minimization for that 1-parameter solve
//   * Use Setting A and Setting B separately
//   * Use only narrow bins in the fit; keep full bin (dp_idx=0) diagnostic only
//
// Default input:
//   results/tables/heep_check_dpBinnedOffsets.csv
//
// Default outputs:
//   results/tables/FixedParamsLinearDpeDp_summary.csv
//   results/tables/FixedParamsLinearDpeDp_residuals.csv
//   results/PNGs/FixedParamsLinearDpeDp/*.png
//
// Observable mapping:
//   W   -> dW
//   Em  -> dE_m
//   Pmz -> dp_m(par)
//   Pmy -> dp_m(perp)
//
// Units:
//   dthe, dthp in mrad
//   dpe, dpp   in 0.1%
//   measured offsets and predictions in MeV

#include <TCanvas.h>
#include <TError.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TSystem.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace FixedParamsLinearDpeDpNS {

struct InputRow {
  int run = -1;
  int dp_idx = -999;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  std::string var;
  double y_obs = 0.0;   // dmu_MeV
  double y_err = 0.0;   // dmu_err_MeV
  int fit_valid = 0;
};

struct Coeffs {
  double c_dthe = 0.0;
  double c_dpe  = 0.0;
  double c_dthp = 0.0;
  double c_dpp  = 0.0;
};

struct FixedParams {
  double dthe = 0.0;
  double dthp = 0.0;
  double dpp  = 0.0;
};

struct ResidualRow {
  InputRow in;
  std::string setting;
  bool used_in_fit = false;
  bool used_in_summary = false;   // narrow bin + fit_valid==1
  double bin_center = 0.0;
  Coeffs coeff;
  double fixed_term = 0.0;
  double dpe_fit = 0.0;
  double dpe_fit_err = 0.0;
  double y_pred = 0.0;
  double resid = 0.0;
};

struct BinFitResult {
  std::string setting;
  int dp_idx = -999;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  double bin_center = 0.0;

  FixedParams fixed;
  double dpe = 0.0;
  double dpe_err = 0.0;

  double chi2 = 0.0;
  int ndf = 0;
  int n_rows_fit = 0;
  int n_runs_fit = 0;

  std::vector<ResidualRow> rows;
};

static std::vector<std::string> Split(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) out.push_back(item);
  return out;
}

static std::string Trim(const std::string& s) {
  const char* ws = " \t\r\n";
  const std::size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const std::size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

static bool ParseDouble(const std::string& s, double& v) {
  char* end = nullptr;
  errno = 0;
  v = std::strtod(s.c_str(), &end);
  return (errno == 0 && end != s.c_str());
}

static bool ParseInt(const std::string& s, int& v) {
  char* end = nullptr;
  errno = 0;
  long x = std::strtol(s.c_str(), &end, 10);
  if (!(errno == 0 && end != s.c_str())) return false;
  v = static_cast<int>(x);
  return true;
}

static void EnsureDir(const std::string& path) {
  if (path.empty()) return;
  const int rc = gSystem->mkdir(path.c_str(), kTRUE);
  if (rc != 0 && gSystem->AccessPathName(path.c_str())) {
    std::cerr << "[ERROR] Failed to create directory: " << path << std::endl;
  }
}

static std::string ParentDir(const std::string& path) {
  std::size_t pos = path.find_last_of('/');
  if (pos == std::string::npos) return "";
  return path.substr(0, pos);
}

static bool IsRunInSettingA(int run) {
  if (run == 23843) return false;
  return (run >= 23839 && run <= 23848);
}

static bool IsRunInSettingB(int run) {
  return (run >= 23849 && run <= 23851);
}

static std::string GetSettingFromRun(int run) {
  if (IsRunInSettingA(run)) return "A";
  if (IsRunInSettingB(run)) return "B";
  return "";
}

static bool AllowedFitBin(const std::string& setting, int dp_idx) {
  if (setting == "A") return (dp_idx >= 1 && dp_idx <= 4);
  if (setting == "B") return (dp_idx >= 1 && dp_idx <= 5);
  return false;
}

static bool IsSupportedVar(const std::string& v) {
  return (v == "W" || v == "Em" || v == "Pmz" || v == "Pmy");
}

static double BinCenter(const InputRow& r) {
  return 0.5 * (r.dp_lo + r.dp_hi);
}

static Coeffs GetCoeffs(const std::string& setting, const std::string& var) {
  Coeffs c;
  if (setting == "A") {
    if (var == "W")   { c.c_dthe = -14.08; c.c_dpe = -8.62; c.c_dthp =  0.00; c.c_dpp =  0.00; }
    if (var == "Em")  { c.c_dthe =   0.00; c.c_dpe = -7.06; c.c_dthp =  0.00; c.c_dpp = -2.10; }
    if (var == "Pmz") { c.c_dthe =  -5.75; c.c_dpe =  4.10; c.c_dthp =  0.00; c.c_dpp =  2.27; }
    if (var == "Pmy") { c.c_dthe =   4.10; c.c_dpe =  5.75; c.c_dthp = -2.27; c.c_dpp =  0.00; }
  } else if (setting == "B") {
    if (var == "W")   { c.c_dthe = -17.33; c.c_dpe = -8.62; c.c_dthp =  0.00; c.c_dpp =  0.00; }
    if (var == "Em")  { c.c_dthe =   0.00; c.c_dpe = -5.66; c.c_dthp =  0.00; c.c_dpp = -3.63; }
    if (var == "Pmz") { c.c_dthe =  -4.30; c.c_dpe =  3.69; c.c_dthp =  0.00; c.c_dpp =  3.74; }
    if (var == "Pmy") { c.c_dthe =   3.69; c.c_dpe =  4.30; c.c_dthp = -3.74; c.c_dpp =  0.00; }
  }
  return c;
}

static std::vector<InputRow> ReadMeasuredCsv(const std::string& csvPath) {
  std::ifstream fin(csvPath.c_str());
  std::vector<InputRow> rows;
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input CSV: " << csvPath << std::endl;
    return rows;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) {
    std::cerr << "[ERROR] Empty input CSV: " << csvPath << std::endl;
    return rows;
  }

  std::vector<std::string> headers = Split(headerLine, ',');
  std::map<std::string, int> idx;
  for (std::size_t i = 0; i < headers.size(); ++i) idx[Trim(headers[i])] = static_cast<int>(i);

  const std::vector<std::string> need = {
    "run","dp_idx","dp_label","dp_lo","dp_hi","var","dmu_MeV","dmu_err_MeV","fit_valid"
  };
  for (std::size_t i = 0; i < need.size(); ++i) {
    if (idx.find(need[i]) == idx.end()) {
      std::cerr << "[ERROR] Missing required column: " << need[i] << std::endl;
      return {};
    }
  }

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    std::vector<std::string> tok = Split(line, ',');
    if (tok.size() < headers.size()) tok.resize(headers.size());

    InputRow r;
    if (!ParseInt(Trim(tok[idx["run"]]), r.run)) continue;
    if (!ParseInt(Trim(tok[idx["dp_idx"]]), r.dp_idx)) continue;
    r.dp_label = Trim(tok[idx["dp_label"]]);
    ParseDouble(Trim(tok[idx["dp_lo"]]), r.dp_lo);
    ParseDouble(Trim(tok[idx["dp_hi"]]), r.dp_hi);
    r.var = Trim(tok[idx["var"]]);
    if (!ParseDouble(Trim(tok[idx["dmu_MeV"]]), r.y_obs)) continue;
    if (!ParseDouble(Trim(tok[idx["dmu_err_MeV"]]), r.y_err)) continue;
    if (!ParseInt(Trim(tok[idx["fit_valid"]]), r.fit_valid)) r.fit_valid = 0;

    rows.push_back(r);
  }

  return rows;
}

static bool SolveSingleDpe(const std::vector<ResidualRow>& rows,
                           double& dpe, double& dpe_err,
                           double& chi2, int& ndf) {
  double SAA = 0.0;
  double SAb = 0.0;
  int nuse = 0;

  for (std::size_t i = 0; i < rows.size(); ++i) {
    const ResidualRow& p = rows[i];
    if (!p.used_in_fit) continue;
    if (p.in.y_err <= 0.0) continue;

    const double w = 1.0 / (p.in.y_err * p.in.y_err);
    const double A = p.coeff.c_dpe;
    const double b = p.in.y_obs - p.fixed_term;

    SAA += w * A * A;
    SAb += w * A * b;
    ++nuse;
  }

  if (nuse < 1 || SAA <= 0.0) {
    dpe = 0.0;
    dpe_err = 0.0;
    chi2 = 0.0;
    ndf = 0;
    return false;
  }

  dpe = SAb / SAA;
  dpe_err = std::sqrt(1.0 / SAA);

  chi2 = 0.0;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const ResidualRow& p = rows[i];
    if (!p.used_in_fit) continue;
    if (p.in.y_err <= 0.0) continue;

    const double ypred = p.fixed_term + p.coeff.c_dpe * dpe;
    const double resid = p.in.y_obs - ypred;
    chi2 += (resid * resid) / (p.in.y_err * p.in.y_err);
  }

  ndf = nuse - 1;
  return true;
}

static BinFitResult FitOneSettingOneBin(const std::string& setting,
                                        int dp_idx,
                                        const std::vector<InputRow>& allRows,
                                        const FixedParams& fixed) {
  BinFitResult out;
  out.setting = setting;
  out.dp_idx = dp_idx;
  out.fixed = fixed;

  std::set<int> runset_fit;

  for (std::size_t i = 0; i < allRows.size(); ++i) {
    const InputRow& r = allRows[i];
    if (GetSettingFromRun(r.run) != setting) continue;
    if (!IsSupportedVar(r.var)) continue;
    if (r.dp_idx != dp_idx) continue;

    ResidualRow rr;
    rr.in = r;
    rr.setting = setting;
    rr.bin_center = BinCenter(r);
    rr.coeff = GetCoeffs(setting, r.var);
    rr.fixed_term = rr.coeff.c_dthe * fixed.dthe
                  + rr.coeff.c_dthp * fixed.dthp
                  + rr.coeff.c_dpp  * fixed.dpp;

    rr.used_in_fit = (r.fit_valid == 1 && AllowedFitBin(setting, r.dp_idx) && r.y_err > 0.0);
    rr.used_in_summary = (r.fit_valid == 1 && AllowedFitBin(setting, r.dp_idx));
    out.rows.push_back(rr);

    if (out.dp_label.empty()) out.dp_label = r.dp_label;
    out.dp_lo = r.dp_lo;
    out.dp_hi = r.dp_hi;
    out.bin_center = rr.bin_center;

    if (rr.used_in_fit) runset_fit.insert(r.run);
  }

  out.n_rows_fit = 0;
  for (std::size_t i = 0; i < out.rows.size(); ++i) {
    if (out.rows[i].used_in_fit) ++out.n_rows_fit;
  }
  out.n_runs_fit = static_cast<int>(runset_fit.size());

  bool ok = SolveSingleDpe(out.rows, out.dpe, out.dpe_err, out.chi2, out.ndf);
  if (!ok) {
    std::cerr << "[WARN] Could not solve dpe for Setting " << setting
              << ", dp_idx=" << dp_idx << std::endl;
  }

  for (std::size_t i = 0; i < out.rows.size(); ++i) {
    ResidualRow& rr = out.rows[i];
    rr.dpe_fit = out.dpe;
    rr.dpe_fit_err = out.dpe_err;
    rr.y_pred = rr.fixed_term + rr.coeff.c_dpe * out.dpe;
    rr.resid = rr.in.y_obs - rr.y_pred;
  }

  std::cout << "\n============================================================\n";
  std::cout << "FixedParamsLinearDpeDp : Setting " << setting
            << ", dp_idx=" << dp_idx
            << ", center=" << out.bin_center << "\n";
  std::cout << "  fixed dthe = " << fixed.dthe << " mrad\n";
  std::cout << "  fixed dthp = " << fixed.dthp << " mrad\n";
  std::cout << "  fixed dpp  = " << fixed.dpp  << " (0.1%)\n";
  std::cout << "  fit model  : one floating dpe per narrow bin\n";
  std::cout << "  dpe = " << out.dpe << " +/- " << out.dpe_err << " (0.1%)\n";
  std::cout << "  chi2 = " << out.chi2 << "\n";
  std::cout << "  ndf  = " << out.ndf << "\n";
  std::cout << "  chi2/ndf = " << (out.ndf > 0 ? out.chi2 / out.ndf : 0.0) << "\n";
  std::cout << "  n_rows_fit = " << out.n_rows_fit << "\n";
  std::cout << "  n_runs_fit = " << out.n_runs_fit << "\n";
  std::cout << "============================================================\n";

  return out;
}

static void WriteSummaryCsv(const std::string& outCsv,
                            const std::vector<BinFitResult>& results) {
  EnsureDir(ParentDir(outCsv));
  std::ofstream fout(outCsv.c_str());

  fout << "setting,dp_idx,dp_label,dp_lo,dp_hi,bin_center,"
          "fixed_dthe_mrad,fixed_dthp_mrad,fixed_dpp_tenthpct,"
          "dpe_tenthpct,dpe_err_tenthpct,chi2,ndf,chi2_ndf,n_rows_fit,n_runs_fit\n";

  for (std::size_t i = 0; i < results.size(); ++i) {
    const BinFitResult& r = results[i];
    fout << r.setting << ","
         << r.dp_idx << ","
         << r.dp_label << ","
         << r.dp_lo << ","
         << r.dp_hi << ","
         << r.bin_center << ","
         << r.fixed.dthe << ","
         << r.fixed.dthp << ","
         << r.fixed.dpp << ","
         << r.dpe << ","
         << r.dpe_err << ","
         << r.chi2 << ","
         << r.ndf << ","
         << (r.ndf > 0 ? r.chi2 / r.ndf : 0.0) << ","
         << r.n_rows_fit << ","
         << r.n_runs_fit << "\n";
  }

  std::cout << "[INFO] Wrote summary CSV: " << outCsv << std::endl;
}

static void WriteResidualCsv(const std::string& outCsv,
                             const std::vector<BinFitResult>& results) {
  EnsureDir(ParentDir(outCsv));
  std::ofstream fout(outCsv.c_str());

  fout << "setting,run,dp_idx,dp_label,dp_lo,dp_hi,bin_center,var,"
          "used_in_fit,used_in_summary,fit_valid,"
          "y_obs_MeV,y_err_MeV,"
          "coeff_dthe,coeff_dpe,coeff_dthp,coeff_dpp,"
          "fixed_dthe_mrad,fixed_dthp_mrad,fixed_dpp_tenthpct,"
          "fixed_term_MeV,dpe_fit_tenthpct,dpe_fit_err_tenthpct,y_pred_MeV,resid_MeV\n";

  for (std::size_t ir = 0; ir < results.size(); ++ir) {
    const BinFitResult& r = results[ir];
    for (std::size_t i = 0; i < r.rows.size(); ++i) {
      const ResidualRow& p = r.rows[i];
      fout << r.setting << ","
           << p.in.run << ","
           << p.in.dp_idx << ","
           << p.in.dp_label << ","
           << p.in.dp_lo << ","
           << p.in.dp_hi << ","
           << p.bin_center << ","
           << p.in.var << ","
           << (p.used_in_fit ? 1 : 0) << ","
           << (p.used_in_summary ? 1 : 0) << ","
           << p.in.fit_valid << ","
           << p.in.y_obs << ","
           << p.in.y_err << ","
           << p.coeff.c_dthe << ","
           << p.coeff.c_dpe << ","
           << p.coeff.c_dthp << ","
           << p.coeff.c_dpp << ","
           << r.fixed.dthe << ","
           << r.fixed.dthp << ","
           << r.fixed.dpp << ","
           << p.fixed_term << ","
           << p.dpe_fit << ","
           << p.dpe_fit_err << ","
           << p.y_pred << ","
           << p.resid << "\n";
    }
  }

  std::cout << "[INFO] Wrote residual CSV: " << outCsv << std::endl;
}

static void DrawDpeVsCenterPng(const std::string& outDir,
                               const std::vector<BinFitResult>& results) {
  EnsureDir(outDir);

  TCanvas* c = new TCanvas("c_dpe_vs_center", "c_dpe_vs_center", 900, 900);
  c->Divide(1, 2);

  const std::vector<std::string> settings = {"A", "B"};

  for (int is = 0; is < 2; ++is) {
    c->cd(is + 1);
    gPad->SetGrid();

    std::vector<double> xs, ys, exs, eys;
    for (std::size_t i = 0; i < results.size(); ++i) {
      const BinFitResult& r = results[i];
      if (r.setting != settings[is]) continue;
      xs.push_back(r.bin_center);
      ys.push_back(r.dpe);
      exs.push_back(0.0);
      eys.push_back(r.dpe_err);
    }

    if (xs.empty()) continue;

    TGraphErrors* g = new TGraphErrors(static_cast<int>(xs.size()),
                                       &xs[0], &ys[0], &exs[0], &eys[0]);
    g->SetTitle(Form("Setting %s;bin center;fitted dpe (0.1%%)", settings[is].c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetLineWidth(2);
    g->Draw("AP");

    TLine* l0 = new TLine(g->GetXaxis()->GetXmin(), 0.0,
                          g->GetXaxis()->GetXmax(), 0.0);
    l0->SetLineStyle(2);
    l0->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.16, 0.86, "One floating dpe solved independently per narrow bin");
  }

  const std::string outName = outDir + "/dpe_vs_bin_center.png";
  std::cout << "[INFO] Writing PNG: " << outName << std::endl;
  c->SaveAs(outName.c_str());
  delete c;
}

static bool ComputeSimpleMeanAndSEM(const std::vector<double>& vals,
                                    double& mean, double& sem) {
  if (vals.empty()) return false;
  const int N = static_cast<int>(vals.size());
  double sum = 0.0;
  for (int i = 0; i < N; ++i) sum += vals[i];
  mean = sum / N;

  if (N <= 1) {
    sem = 0.0;
    return true;
  }

  double s2 = 0.0;
  for (int i = 0; i < N; ++i) {
    const double d = vals[i] - mean;
    s2 += d * d;
  }
  const double rms = std::sqrt(s2 / (N - 1));
  sem = rms / std::sqrt(static_cast<double>(N));
  return true;
}

static void DrawResidualMeansPng(const std::string& outDir,
                                 const std::string& setting,
                                 const std::vector<BinFitResult>& results) {
  EnsureDir(outDir);

  TCanvas* c = new TCanvas(Form("c_resid_%s", setting.c_str()),
                           Form("c_resid_%s", setting.c_str()), 1100, 900);
  c->Divide(2, 2);

  const std::vector<std::string> vars = {"W", "Em", "Pmz", "Pmy"};

  for (int iv = 0; iv < 4; ++iv) {
    c->cd(iv + 1);
    gPad->SetGrid();

    std::vector<double> xs, ys, exs, eys;

    for (std::size_t ib = 0; ib < results.size(); ++ib) {
      const BinFitResult& bin = results[ib];
      if (bin.setting != setting) continue;

      std::vector<double> resids;
      for (std::size_t i = 0; i < bin.rows.size(); ++i) {
        const ResidualRow& rr = bin.rows[i];
        if (!rr.used_in_summary) continue;
        if (rr.in.var != vars[iv]) continue;
        resids.push_back(rr.resid);
      }

      double mean = 0.0, sem = 0.0;
      if (!ComputeSimpleMeanAndSEM(resids, mean, sem)) continue;

      xs.push_back(bin.bin_center);
      ys.push_back(mean);
      exs.push_back(0.0);
      eys.push_back(sem);
    }

    if (xs.empty()) continue;

    TGraphErrors* g = new TGraphErrors(static_cast<int>(xs.size()),
                                       &xs[0], &ys[0], &exs[0], &eys[0]);
    g->SetTitle(Form("Setting %s residual means: %s;bin center;mean residual over runs (MeV)",
                     setting.c_str(), vars[iv].c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetLineWidth(2);
    g->Draw("AP");

    TLine* l0 = new TLine(g->GetXaxis()->GetXmin(), 0.0,
                          g->GetXaxis()->GetXmax(), 0.0);
    l0->SetLineStyle(2);
    l0->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.14, 0.86, "Residual point = simple mean over runs");
    lat.DrawLatex(0.14, 0.80, "Error bar = SEM = RMS / #sqrt{N}");
  }

  const std::string outName = outDir + Form("/residual_means_setting%s.png", setting.c_str());
  std::cout << "[INFO] Writing PNG: " << outName << std::endl;
  c->SaveAs(outName.c_str());
  delete c;
}

} // namespace FixedParamsLinearDpeDpNS

void FixedParamsLinearDpeDp(
    const char* inputCsv = "results/tables/heep_check_dpBinnedOffsets.csv",
    const char* outSummaryCsv = "results/tables/FixedParamsLinearDpeDp_summary.csv",
    const char* outResidualCsv = "results/tables/FixedParamsLinearDpeDp_residuals.csv",
    const char* outPngDir = "results/PNGs/FixedParamsLinearDpeDp") {

  using namespace FixedParamsLinearDpeDpNS;

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::vector<InputRow> rows = ReadMeasuredCsv(inputCsv);
  if (rows.empty()) {
    std::cerr << "[ERROR] No input rows loaded. Aborting." << std::endl;
    return;
  }

  FixedParams fixedA;
  fixedA.dthe = -0.30;
  fixedA.dthp = -7.20;
  fixedA.dpp  = -3.30;

  FixedParams fixedB;
  fixedB.dthe = -1.15;
  fixedB.dthp = -4.30;
  fixedB.dpp  = -1.46;

  std::vector<BinFitResult> results;

  for (int dp = 1; dp <= 4; ++dp) {
    results.push_back(FitOneSettingOneBin("A", dp, rows, fixedA));
  }
  for (int dp = 1; dp <= 5; ++dp) {
    results.push_back(FitOneSettingOneBin("B", dp, rows, fixedB));
  }

  WriteSummaryCsv(outSummaryCsv, results);
  WriteResidualCsv(outResidualCsv, results);

  EnsureDir(outPngDir);

  DrawDpeVsCenterPng(outPngDir, results);
  DrawResidualMeansPng(outPngDir, "A", results);
  DrawResidualMeansPng(outPngDir, "B", results);

  std::cout << "[INFO] FixedParamsLinearDpeDp finished.\n"
            << "       Summary CSV   : " << outSummaryCsv << "\n"
            << "       Residual CSV  : " << outResidualCsv << "\n"
            << "       PNG directory : " << outPngDir << std::endl;
}
