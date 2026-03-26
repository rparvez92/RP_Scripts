// DiagnoseDpBinnedOffsets.C
//
// Purpose:
//   Read heep_check_dpBinnedOffsets.csv and build a setting-wise diagnosis table.
//   This does NOT refit anything. It only summarizes existing measured-offset rows.
//
//   Settings used here:
//     Setting A: runs 23839,23840,23841,23842,23844,23845,23846,23847,23848
//     Setting B: runs 23849,23850,23851
//
// Output:
//   1) diagnosis CSV with one row per (setting, dp bin)
//   2) terminal summary table
//
// Example:
//   root -l -b -q 'DiagnoseDpBinnedOffsets.C+(
//      "results/tables/heep_check_dpBinnedOffsets.csv",
//      "results/tables/heep_check_diagnosis_bySetting.csv"
//   )'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iomanip>
#include <algorithm>

#include "TString.h"
#include "TSystem.h"

struct Row {
  int run = -999;
  int dp_idx = -999;
  TString dp_label = "";
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  TString var = "";

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

  int nRows = 0;
  int nValid = 0;

  std::set<int> runs;
  std::set<TString> vars;

  std::map<int,int> statusD_counts;
  std::map<int,int> statusS_counts;

  double entriesS_sum = 0.0;
  double entriesS_min =  1e300;
  double entriesS_max = -1e300;
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

static TString Key(const TString& setting, int dp_idx) {
  return Form("%s__dp%d", setting.Data(), dp_idx);
}

static TString DominantStatusString(const std::map<int,int>& counts) {
  if (counts.empty()) return "none";

  int bestStatus = 999999;
  int bestCount  = -1;
  for (const auto& kv : counts) {
    if (kv.second > bestCount) {
      bestCount  = kv.second;
      bestStatus = kv.first;
    }
  }

  // If there is only one unique status, print that directly.
  if ((int)counts.size() == 1) return Form("%d", bestStatus);

  // Otherwise print all counts compactly, e.g. "0(24);-999(12)"
  TString s = "";
  bool first = true;
  for (const auto& kv : counts) {
    if (!first) s += ";";
    s += Form("%d(%d)", kv.first, kv.second);
    first = false;
  }
  return s;
}

void DiagnoseDpBinnedOffsets(const char* inCsvC  = "results/tables/heep_check_dpBinnedOffsets.csv",
                             const char* outCsvC = "results/tables/heep_check_diagnosis_bySetting.csv") {
  TString inCsv(inCsvC);
  TString outCsv(outCsvC);

  std::ifstream fin(inCsv.Data());
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open input CSV: " << inCsv << "\n";
    return;
  }

  // Ensure output directory exists
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

  if (!need("run")      ||
      !need("dp_idx")   ||
      !need("dp_label") ||
      !need("dp_lo")    ||
      !need("dp_hi")    ||
      !need("var")      ||
      !need("entriesS") ||
      !need("statusD")  ||
      !need("statusS")  ||
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
    r.run       = tok[col["run"]].Atoi();
    r.dp_idx    = tok[col["dp_idx"]].Atoi();
    r.dp_label  = tok[col["dp_label"]];
    r.dp_lo     = tok[col["dp_lo"]].Atof();
    r.dp_hi     = tok[col["dp_hi"]].Atof();
    r.var       = tok[col["var"]];
    r.entriesS  = tok[col["entriesS"]].Atof();
    r.statusD   = tok[col["statusD"]].Atoi();
    r.statusS   = tok[col["statusS"]].Atoi();
    r.fit_valid = tok[col["fit_valid"]].Atoi();

    TString setting = GetSettingLabel(r.run);
    if (setting == "UNKNOWN") {
      ++nUnknownSetting;
      continue;
    }

    TString k = Key(setting, r.dp_idx);
    if (!M.count(k)) {
      Agg a;
      a.setting  = setting;
      a.dp_idx   = r.dp_idx;
      a.dp_label = r.dp_label;
      a.dp_lo    = r.dp_lo;
      a.dp_hi    = r.dp_hi;
      M[k] = a;
    }

    Agg& a = M[k];
    a.nRows++;
    if (r.fit_valid == 1) a.nValid++;

    a.runs.insert(r.run);
    a.vars.insert(r.var);

    a.statusD_counts[r.statusD]++;
    a.statusS_counts[r.statusS]++;

    a.entriesS_sum += r.entriesS;
    a.entriesS_min = std::min(a.entriesS_min, r.entriesS);
    a.entriesS_max = std::max(a.entriesS_max, r.entriesS);

    ++nInputRows;
  }
  fin.close();

  std::vector<Agg> V;
  for (const auto& kv : M) V.push_back(kv.second);

  std::sort(V.begin(), V.end(),
            [](const Agg& a, const Agg& b) {
              if (a.setting != b.setting) return a.setting < b.setting;
              return a.dp_idx < b.dp_idx;
            });

  std::ofstream fout(outCsv.Data(), std::ios::out);
  if (!fout.is_open()) {
    std::cerr << "[ERROR] Cannot create output CSV: " << outCsv << "\n";
    return;
  }

  fout
    << "setting,dp_idx,dp_label,dp_lo,dp_hi,"
    << "n_runs,n_vars,n_rows,n_valid,n_invalid,"
    << "statusD_summary,statusS_summary,"
    << "entriesS_mean,entriesS_min,entriesS_max,"
    << "diagnosis,usable_for_next_stage"
    << "\n";

  std::cout << "\n=== Dp-binned offset diagnosis by setting ===\n\n";
  std::cout
    << std::left
    << std::setw(8)  << "Setting"
    << std::setw(8)  << "dp_idx"
    << std::setw(10) << "dp_label"
    << std::setw(12) << "valid/rows"
    << std::setw(14) << "statusD"
    << std::setw(18) << "statusS"
    << std::setw(14) << "entriesS_mean"
    << "diagnosis"
    << "\n";

  for (const auto& a : V) {
    double entriesS_mean = (a.nRows > 0 ? a.entriesS_sum / a.nRows : 0.0);
    int nInvalid = a.nRows - a.nValid;

    TString statusD_summary = DominantStatusString(a.statusD_counts);
    TString statusS_summary = DominantStatusString(a.statusS_counts);

    TString diagnosis = "usable";
    int usable = 1;

    if (a.nValid == 0) {
      usable = 0;
      if (a.entriesS_max == 0.0) {
        diagnosis = "unusable: SIM has zero entries";
      } else {
        diagnosis = "unusable: no valid fits";
      }
    } else if (a.nValid < a.nRows) {
      usable = 0;
      diagnosis = "mixed validity: inspect row-level failures";
    }

    fout
      << a.setting << ","
      << a.dp_idx << ","
      << a.dp_label << ","
      << a.dp_lo << ","
      << a.dp_hi << ","
      << a.runs.size() << ","
      << a.vars.size() << ","
      << a.nRows << ","
      << a.nValid << ","
      << nInvalid << ","
      << "\"" << statusD_summary << "\","
      << "\"" << statusS_summary << "\","
      << std::fixed << std::setprecision(3) << entriesS_mean << ","
      << a.entriesS_min << ","
      << a.entriesS_max << ","
      << "\"" << diagnosis << "\","
      << usable
      << "\n";

    std::cout
      << std::left
      << std::setw(8)  << a.setting.Data()
      << std::setw(8)  << a.dp_idx
      << std::setw(10) << a.dp_label.Data()
      << std::setw(12) << Form("%d/%d", a.nValid, a.nRows)
      << std::setw(14) << statusD_summary.Data()
      << std::setw(18) << statusS_summary.Data()
      << std::setw(14) << Form("%.1f", entriesS_mean)
      << diagnosis.Data()
      << "\n";
  }

  fout.close();

  std::cout << "\n[INFO] Read rows: " << nInputRows << "\n";
  if (nUnknownSetting > 0) {
    std::cout << "[WARN] Skipped rows with UNKNOWN setting: " << nUnknownSetting << "\n";
  }
  std::cout << "[INFO] Wrote diagnosis CSV: " << outCsv << "\n\n";
}
