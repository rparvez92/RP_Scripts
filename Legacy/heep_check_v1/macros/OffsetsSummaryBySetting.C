// OffsetsSummaryBySetting.C
//
// Purpose:
//   Read MeasuredOffsetsByRun.csv and summarize measured offsets
//   by (setting, dp bin, variable), using only fit_valid == 1.
//
//   Settings used here:
//     Setting A: runs 23839,23840,23841,23842,23844,23845,23846,23847,23848
//     Setting B: runs 23849,23850,23851
//
// Output:
//   1) summary CSV with one row per (setting, dp bin, variable)
//   2) terminal summary table
//
// Example:
//   root -l -b -q 'macros/OffsetsSummaryBySetting.C+("results/tables/MeasuredOffsetsByRun.csv", "results/tables/OffsetsSummaryBySetting.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "TString.h"
#include "TSystem.h"

struct Row {
  int run = -999;
  int dp_idx = -999;
  TString dp_label = "";
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  TString var = "";

  double dmu_MeV = -999.0;
  double dmu_err_MeV = -999.0;

  double entriesD = 0.0;
  double entriesS = 0.0;

  int statusD = -999;
  int statusS = -999;
  int fit_valid = 0;
};

struct Agg {
  TString setting = "";
  int dp_idx = -999;
  TString dp_label = "";
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  TString var = "";

  int nRowsAll = 0;
  int nValid = 0;

  std::set<int> runsAll;
  std::set<int> runsValid;

  std::vector<double> dmuVals;
  std::vector<double> dmuErrVals;

  double minVal =  1e300;
  double maxVal = -1e300;
};

static std::vector<TString> SplitCSV(const std::string& line) {
  std::vector<TString> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(TString(item));
  }
  return out;
}

static bool IsSettingA(int run) {
  static const std::set<int> sA = {
    23839, 23840, 23841, 23842, 23844, 23845, 23846, 23847, 23848
  };
  return (sA.count(run) > 0);
}

static bool IsSettingB(int run) {
  static const std::set<int> sB = {
    23849, 23850, 23851
  };
  return (sB.count(run) > 0);
}

static TString GetSettingLabel(int run) {
  if (IsSettingA(run)) return "A";
  if (IsSettingB(run)) return "B";
  return "UNKNOWN";
}

static TString Key(const TString& setting, int dp_idx, const TString& var) {
  return Form("%s__dp%d__%s", setting.Data(), dp_idx, var.Data());
}

static double Mean(const std::vector<double>& v) {
  if (v.empty()) return 0.0;
  double s = 0.0;
  for (double x : v) s += x;
  return s / v.size();
}

static double SampleRMS(const std::vector<double>& v) {
  int n = (int)v.size();
  if (n < 2) return 0.0;
  double mu = Mean(v);
  double s2 = 0.0;
  for (double x : v) s2 += (x - mu) * (x - mu);
  return std::sqrt(s2 / (n - 1));
}

static double WeightedMean(const std::vector<double>& vals,
                           const std::vector<double>& errs) {
  if (vals.empty() || vals.size() != errs.size()) return 0.0;
  double sw = 0.0;
  double swx = 0.0;
  for (size_t i = 0; i < vals.size(); ++i) {
    double w = 1.0 / (errs[i] * errs[i]);
    sw += w;
    swx += w * vals[i];
  }
  return (sw > 0.0 ? swx / sw : 0.0);
}

static double WeightedMeanErr(const std::vector<double>& errs) {
  if (errs.empty()) return 0.0;
  double sw = 0.0;
  for (double err : errs) {
    double w = 1.0 / (err * err);
    sw += w;
  }
  return (sw > 0.0 ? std::sqrt(1.0 / sw) : 0.0);
}

void OffsetsSummaryBySetting(
    const char* inCsvC  = "results/tables/MeasuredOffsetsByRun.csv",
    const char* outCsvC = "results/tables/OffsetsSummaryBySetting.csv") {

  TString inCsv(inCsvC);
  TString outCsv(outCsvC);

  std::ifstream fin(inCsv.Data());
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open input CSV: " << inCsv << "\n";
    return;
  }

  TString parent = outCsv;
  parent = gSystem->DirName(parent);
  if (!parent.IsNull() && gSystem->AccessPathName(parent)) {
    gSystem->mkdir(parent, kTRUE);
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty CSV: " << inCsv << "\n";
    return;
  }

  std::vector<TString> header = SplitCSV(line);
  std::map<std::string,int> col;
  for (int i = 0; i < (int)header.size(); ++i) {
    col[header[i].Data()] = i;
  }

  auto need = [&](const char* name) {
    if (!col.count(name)) {
      std::cerr << "[ERROR] Missing required column: " << name << "\n";
      return false;
    }
    return true;
  };

  if (!need("run")         ||
      !need("dp_idx")      ||
      !need("dp_label")    ||
      !need("dp_lo")       ||
      !need("dp_hi")       ||
      !need("var")         ||
      !need("dmu_MeV")     ||
      !need("dmu_err_MeV") ||
      !need("entriesD")    ||
      !need("entriesS")    ||
      !need("statusD")     ||
      !need("statusS")     ||
      !need("fit_valid")) {
    return;
  }

  std::map<TString, Agg> M;

  int nInputRows = 0;
  int nUnknownSetting = 0;

  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    std::vector<TString> tok = SplitCSV(line);
    if ((int)tok.size() < (int)header.size()) continue;

    Row r;
    r.run         = tok[col["run"]].Atoi();
    r.dp_idx      = tok[col["dp_idx"]].Atoi();
    r.dp_label    = tok[col["dp_label"]];
    r.dp_lo       = tok[col["dp_lo"]].Atof();
    r.dp_hi       = tok[col["dp_hi"]].Atof();
    r.var         = tok[col["var"]];
    r.dmu_MeV     = tok[col["dmu_MeV"]].Atof();
    r.dmu_err_MeV = tok[col["dmu_err_MeV"]].Atof();
    r.entriesD    = tok[col["entriesD"]].Atof();
    r.entriesS    = tok[col["entriesS"]].Atof();
    r.statusD     = tok[col["statusD"]].Atoi();
    r.statusS     = tok[col["statusS"]].Atoi();
    r.fit_valid   = tok[col["fit_valid"]].Atoi();

    TString setting = GetSettingLabel(r.run);
    if (setting == "UNKNOWN") {
      ++nUnknownSetting;
      continue;
    }

    TString k = Key(setting, r.dp_idx, r.var);
    if (!M.count(k)) {
      Agg a;
      a.setting  = setting;
      a.dp_idx   = r.dp_idx;
      a.dp_label = r.dp_label;
      a.dp_lo    = r.dp_lo;
      a.dp_hi    = r.dp_hi;
      a.var      = r.var;
      M[k] = a;
    }

    Agg& a = M[k];
    a.nRowsAll++;
    a.runsAll.insert(r.run);

    if (r.fit_valid == 1) {
      a.nValid++;
      a.runsValid.insert(r.run);
      a.dmuVals.push_back(r.dmu_MeV);
      a.dmuErrVals.push_back(r.dmu_err_MeV);
      a.minVal = std::min(a.minVal, r.dmu_MeV);
      a.maxVal = std::max(a.maxVal, r.dmu_MeV);
    }

    ++nInputRows;
  }
  fin.close();

  std::vector<Agg> V;
  for (const auto& kv : M) V.push_back(kv.second);

  std::sort(V.begin(), V.end(),
            [](const Agg& a, const Agg& b) {
              if (a.setting != b.setting) return a.setting < b.setting;
              if (a.dp_idx != b.dp_idx) return a.dp_idx < b.dp_idx;
              return a.var < b.var;
            });

  std::ofstream fout(outCsv.Data(), std::ios::out);
  if (!fout.is_open()) {
    std::cerr << "[ERROR] Cannot create output CSV: " << outCsv << "\n";
    return;
  }

  fout
    << "setting,dp_idx,dp_label,dp_lo,dp_hi,var,"
    << "n_runs_all,n_runs_valid,n_rows_all,n_valid,n_invalid,"
    << "mean_dmu_MeV,weighted_mean_dmu_MeV,weighted_mean_err_dmu_MeV,rms_dmu_MeV,sem_dmu_MeV,"
    << "min_dmu_MeV,max_dmu_MeV,"
    << "mean_input_dmuErr_MeV,"
    << "summary_status"
    << "\n";

  std::cout << "\n=== Measured offset summary by setting / dp bin / variable ===\n\n";
  std::cout
    << std::left
    << std::setw(8)  << "Setting"
    << std::setw(8)  << "dp_idx"
    << std::setw(10) << "dp_label"
    << std::setw(8)  << "var"
    << std::setw(12) << "valid/all"
    << std::setw(16) << "mean_dmu(MeV)"
    << std::setw(14) << "RMS"
    << std::setw(14) << "SEM"
    << "status"
    << "\n";

  for (const auto& a : V) {
    int nInvalid = a.nRowsAll - a.nValid;

    double meanVal = 0.0;
    double weightedMeanVal = 0.0;
    double weightedMeanErrVal = 0.0;
    double rmsVal = 0.0;
    double semVal = 0.0;
    double meanErr = 0.0;

    if (!a.dmuVals.empty()) {
      meanVal = Mean(a.dmuVals);
      weightedMeanVal = WeightedMean(a.dmuVals, a.dmuErrVals);
      weightedMeanErrVal = WeightedMeanErr(a.dmuErrVals);
      rmsVal  = SampleRMS(a.dmuVals);
      semVal  = (a.dmuVals.size() > 0 ? rmsVal / std::sqrt((double)a.dmuVals.size()) : 0.0);
      meanErr = Mean(a.dmuErrVals);
    }

    TString status = "ok";
    if (a.nValid == 0) {
      status = "no valid rows";
    } else if (nInvalid > 0) {
      status = "mixed validity";
    }

    double minVal = (a.nValid > 0 ? a.minVal : -999.0);
    double maxVal = (a.nValid > 0 ? a.maxVal : -999.0);

    fout
      << a.setting << ","
      << a.dp_idx << ","
      << a.dp_label << ","
      << a.dp_lo << ","
      << a.dp_hi << ","
      << a.var << ","
      << a.runsAll.size() << ","
      << a.runsValid.size() << ","
      << a.nRowsAll << ","
      << a.nValid << ","
      << nInvalid << ","
      << std::fixed << std::setprecision(6)
      << meanVal << ","
      << weightedMeanVal << ","
      << weightedMeanErrVal << ","
      << rmsVal << ","
      << semVal << ","
      << minVal << ","
      << maxVal << ","
      << meanErr << ","
      << "\"" << status << "\""
      << "\n";

    std::cout
      << std::left
      << std::setw(8)  << a.setting.Data()
      << std::setw(8)  << a.dp_idx
      << std::setw(10) << a.dp_label.Data()
      << std::setw(8)  << a.var.Data()
      << std::setw(12) << Form("%d/%d", a.nValid, a.nRowsAll)
      << std::setw(16) << Form("%.3f", meanVal)
      << std::setw(14) << Form("%.3f", rmsVal)
      << std::setw(14) << Form("%.3f", semVal)
      << status.Data()
      << "\n";
  }

  fout.close();

  std::cout << "\n[INFO] Read rows: " << nInputRows << "\n";
  if (nUnknownSetting > 0) {
    std::cout << "[WARN] Skipped rows with UNKNOWN setting: " << nUnknownSetting << "\n";
  }
  std::cout << "[INFO] Wrote summary CSV: " << outCsv << "\n\n";
}
