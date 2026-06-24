// File: heep_check_v0/macros/heepcheck_reference_from_dpeOnly.C
//
// Purpose:
//   Build a nonzero reference point for each run using:
//     de0_ref  = 0
//     dpe_ref  = full-bin dpe from heepcheck_inversion_dpeOnly.csv
//     dthe_ref = solved from observed W offset
//     dthp_ref = solved from observed Pmy offset
//     dpp_ref  = least-squares compromise between observed Em and Pmz
//
// Input:
//   - results/tables/heep_dp_binned_offsets.csv     (observed offsets; use dp_idx=0 only)
//   - results/tables/heepcheck_inversion_dpeOnly.csv (full-bin dpe_ref; use dp_idx=0 only)
//
// Output:
//   - results/tables/heepcheck_reference_from_dpeOnly.csv
//
// Conventions:
//   - observed offsets are data - sim, taken from dmu_MeV
//   - dpe_ref and dpp_ref are in heepcheck momentum units [0.1%]
//   - dthe_ref and dthp_ref are in [mrad]
//   - Pmy in the offsets CSV corresponds to the transverse component mislabeled as dPmx in
//     inversion CSVs; here it is consistently written as Pmy.
//
// Suggested usage from heep_check_v0/:
//   root -l -b -q 'macros/heepcheck_reference_from_dpeOnly.C("one",23840)'
//   root -l -b -q 'macros/heepcheck_reference_from_dpeOnly.C("all",-1)'

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct OffsetObsRow {
  int run = -1;
  int dp_idx = -999;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  std::string var;
  double dmu_MeV = 0.0;
  double dmu_err_MeV = 0.0;
};

struct DpeOnlyRow {
  int run = -1;
  int dp_idx = -999;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  double dpe = 0.0;
};

struct ObsRunData {
  bool hasW = false;
  bool hasEm = false;
  bool hasPmz = false;
  bool hasPmy = false;

  int run = -1;
  int dp_idx = 0;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;

  double W_obs = 0.0;
  double Em_obs = 0.0;
  double Pmz_obs = 0.0;
  double Pmy_obs = 0.0;

  double W_err = 0.0;
  double Em_err = 0.0;
  double Pmz_err = 0.0;
  double Pmy_err = 0.0;
};

struct CoeffSet {
  std::string tag;
  double aW = 0.0, bW = 0.0;
  double bE = 0.0, cE = 0.0;
  double aZ = 0.0, bZ = 0.0, cZ = 0.0;
  double aY = 0.0, bY = 0.0, cY = 0.0;
};

struct RefResult {
  bool ok = false;
  int run = -1;
  int dp_idx = 0;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  std::string setting;
  std::string coeff_tag;
  std::string notes;

  double W_obs = 0.0;
  double Em_obs = 0.0;
  double Pmz_obs = 0.0;
  double Pmy_obs = 0.0;

  double W_err = 0.0;
  double Em_err = 0.0;
  double Pmz_err = 0.0;
  double Pmy_err = 0.0;

  double de0_ref = 0.0;
  double dpe_ref_input = 0.0;
  double dthe_ref = 0.0;
  double dthp_ref = 0.0;
  double dpp_ref = 0.0;

  double W_pred_ref = 0.0;
  double Em_pred_ref = 0.0;
  double Pmz_pred_ref = 0.0;
  double Pmy_pred_ref = 0.0;

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

constexpr double kLargeDthe = 3.0;
constexpr double kLargeDthp = 3.0;
constexpr double kLargeDpp  = 5.0;
constexpr double kPoorResidual = 10.0;

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

bool EnsurePathForFile(const std::string& path) {
  try {
    std::filesystem::path p(path);
    auto parent = p.parent_path();
    if (!parent.empty() && !std::filesystem::exists(parent)) {
      std::filesystem::create_directories(parent);
    }
  } catch (...) {
    return false;
  }
  return true;
}

bool RequiredColumnsPresent(const std::map<std::string,int>& idx,
                            const std::vector<std::string>& req,
                            const std::string& path) {
  bool ok = true;
  for (const auto& c : req) {
    if (!idx.count(c)) {
      std::cerr << "[ERROR] Missing required column '" << c << "' in " << path << "\n";
      ok = false;
    }
  }
  return ok;
}

bool IsSettingA(int run) {
  return (run >= 23839 && run <= 23848 && run != 23843);
}

bool IsSettingB(int run) {
  return (run >= 23849 && run <= 23851);
}

bool GetCoeffSetForRun(int run, CoeffSet& c) {
  if (IsSettingA(run)) {
    c.tag = "A";
    c.aW = -14.08; c.bW = -8.62;
    c.bE = -7.06;  c.cE = -2.10;
    c.aZ = -5.75;  c.bZ =  4.10; c.cZ =  2.27;
    c.aY =  4.10;  c.bY =  5.75; c.cY = -2.27;
    return true;
  }
  if (IsSettingB(run)) {
    c.tag = "B";
    c.aW = -17.33; c.bW = -8.62;
    c.bE =  -5.66; c.cE = -3.63;
    c.aZ =  -4.30; c.bZ =  3.69; c.cZ =  3.74;
    c.aY =   3.69; c.bY =  4.30; c.cY = -3.74;
    return true;
  }
  return false;
}

bool LoadObservedOffsets(const std::string& path, std::map<int, ObsRunData>& out) {
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open observed-offset CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty observed-offset CSV: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);
  if (!RequiredColumnsPresent(idx, {"run","dp_idx","dp_label","dp_lo","dp_hi","var","dmu_MeV","dmu_err_MeV"}, path)) {
    return false;
  }

  int nread = 0;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto tok = SplitCSVSimple(line);
    auto getS = [&](const std::string& k)->std::string {
      auto it = idx.find(k);
      if (it == idx.end() || it->second >= (int)tok.size()) return "";
      return tok[it->second];
    };

    OffsetObsRow r;
    if (!ToInt(getS("run"), r.run)) continue;
    if (!ToInt(getS("dp_idx"), r.dp_idx)) continue;
    if (r.dp_idx != 0) continue; // full bin only
    r.dp_label = getS("dp_label");
    ToDouble(getS("dp_lo"), r.dp_lo);
    ToDouble(getS("dp_hi"), r.dp_hi);
    r.var = getS("var");
    ToDouble(getS("dmu_MeV"), r.dmu_MeV);
    ToDouble(getS("dmu_err_MeV"), r.dmu_err_MeV);

    ObsRunData& d = out[r.run];
    d.run = r.run;
    d.dp_idx = r.dp_idx;
    d.dp_label = r.dp_label;
    d.dp_lo = r.dp_lo;
    d.dp_hi = r.dp_hi;

    if (r.var == "W") {
      d.W_obs = r.dmu_MeV; d.W_err = r.dmu_err_MeV; d.hasW = true;
    } else if (r.var == "Em") {
      d.Em_obs = r.dmu_MeV; d.Em_err = r.dmu_err_MeV; d.hasEm = true;
    } else if (r.var == "Pmz") {
      d.Pmz_obs = r.dmu_MeV; d.Pmz_err = r.dmu_err_MeV; d.hasPmz = true;
    } else if (r.var == "Pmy") {
      d.Pmy_obs = r.dmu_MeV; d.Pmy_err = r.dmu_err_MeV; d.hasPmy = true;
    }
    ++nread;
  }

  std::cout << "[INFO] Parsed " << nread << " full-bin observed rows from " << path << "\n";
  return true;
}

bool LoadDpeOnly(const std::string& path, std::map<int, DpeOnlyRow>& out) {
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open dpeOnly inversion CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty dpeOnly inversion CSV: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);
  if (!RequiredColumnsPresent(idx, {"run","dp_idx","dp_label","dp_lo","dp_hi","dpe"}, path)) {
    return false;
  }

  int nread = 0;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto tok = SplitCSVSimple(line);
    auto getS = [&](const std::string& k)->std::string {
      auto it = idx.find(k);
      if (it == idx.end() || it->second >= (int)tok.size()) return "";
      return tok[it->second];
    };

    DpeOnlyRow r;
    if (!ToInt(getS("run"), r.run)) continue;
    if (!ToInt(getS("dp_idx"), r.dp_idx)) continue;
    if (r.dp_idx != 0) continue; // full bin only
    r.dp_label = getS("dp_label");
    ToDouble(getS("dp_lo"), r.dp_lo);
    ToDouble(getS("dp_hi"), r.dp_hi);
    if (!ToDouble(getS("dpe"), r.dpe)) continue;
    out[r.run] = r;
    ++nread;
  }

  std::cout << "[INFO] Parsed " << nread << " full-bin dpeOnly rows from " << path << "\n";
  return true;
}

bool HasAllObservables(const ObsRunData& d) {
  return d.hasW && d.hasEm && d.hasPmz && d.hasPmy;
}

RefResult BuildReferenceForRun(int run,
                               const ObsRunData& obs,
                               const DpeOnlyRow& dpeRow) {
  RefResult r;
  r.run = run;
  r.dp_idx = obs.dp_idx;
  r.dp_label = obs.dp_label;
  r.dp_lo = obs.dp_lo;
  r.dp_hi = obs.dp_hi;
  r.W_obs = obs.W_obs;
  r.Em_obs = obs.Em_obs;
  r.Pmz_obs = obs.Pmz_obs;
  r.Pmy_obs = obs.Pmy_obs;
  r.W_err = obs.W_err;
  r.Em_err = obs.Em_err;
  r.Pmz_err = obs.Pmz_err;
  r.Pmy_err = obs.Pmy_err;
  r.de0_ref = 0.0;
  r.dpe_ref_input = dpeRow.dpe;

  CoeffSet c;
  if (!GetCoeffSetForRun(run, c)) {
    r.notes = "run_outside_setting_map";
    return r;
  }
  r.setting = std::string("setting_") + c.tag;
  r.coeff_tag = c.tag;

  // Step 1: dthe from W
  if (std::fabs(c.aW) < 1e-12) {
    r.notes = "aW_zero";
    return r;
  }
  r.dthe_ref = (r.W_obs - c.bW * r.dpe_ref_input) / c.aW;

  // Step 2: dthp from Pmy
  if (std::fabs(c.cY) < 1e-12) {
    r.notes = "cY_zero";
    return r;
  }
  r.dthp_ref = (r.Pmy_obs - c.aY * r.dthe_ref - c.bY * r.dpe_ref_input) / c.cY;

  // Step 3: dpp from LS compromise between Em and Pmz
  const double AE = c.bE * r.dpe_ref_input;
  const double AZ = c.aZ * r.dthe_ref + c.bZ * r.dpe_ref_input;
  const double denom = c.cE * c.cE + c.cZ * c.cZ;
  if (denom <= 0.0) {
    r.notes = "lsq_denom_nonpositive";
    return r;
  }
  r.dpp_ref = (c.cE * (r.Em_obs - AE) + c.cZ * (r.Pmz_obs - AZ)) / denom;

  // Predictions at reference point
  r.W_pred_ref   = c.aW * r.dthe_ref + c.bW * r.dpe_ref_input;
  r.Em_pred_ref  = c.bE * r.dpe_ref_input + c.cE * r.dpp_ref;
  r.Pmz_pred_ref = c.aZ * r.dthe_ref + c.bZ * r.dpe_ref_input + c.cZ * r.dpp_ref;
  r.Pmy_pred_ref = c.aY * r.dthe_ref + c.bY * r.dpe_ref_input + c.cY * r.dthp_ref;

  // Residuals and SSE
  r.r_W   = r.W_pred_ref   - r.W_obs;
  r.r_Em  = r.Em_pred_ref  - r.Em_obs;
  r.r_Pmz = r.Pmz_pred_ref - r.Pmz_obs;
  r.r_Pmy = r.Pmy_pred_ref - r.Pmy_obs;
  r.SSE_ref = r.r_W*r.r_W + r.r_Em*r.r_Em + r.r_Pmz*r.r_Pmz + r.r_Pmy*r.r_Pmy;

  // Warnings
  r.large_dthe_warn = (std::fabs(r.dthe_ref) > kLargeDthe) ? 1 : 0;
  r.large_dthp_warn = (std::fabs(r.dthp_ref) > kLargeDthp) ? 1 : 0;
  r.large_dpp_warn  = (std::fabs(r.dpp_ref)  > kLargeDpp ) ? 1 : 0;
  r.poor_ref_consistency_warn = (std::fabs(r.r_Em) > kPoorResidual || std::fabs(r.r_Pmz) > kPoorResidual) ? 1 : 0;

  std::ostringstream note;
  bool first = true;
  auto addNote = [&](const std::string& s){ if (!first) note << ";"; note << s; first = false; };
  if (r.large_dthe_warn) addNote("large_dthe");
  if (r.large_dthp_warn) addNote("large_dthp");
  if (r.large_dpp_warn)  addNote("large_dpp");
  if (r.poor_ref_consistency_warn) addNote("poor_ref_consistency");
  r.notes = note.str();

  r.ok = true;
  return r;
}

bool WriteCSVHeaderIfNeeded(const std::string& out_csv) {
  if (!EnsurePathForFile(out_csv)) {
    std::cerr << "[ERROR] Could not create parent directory for: " << out_csv << "\n";
    return false;
  }
  const bool exists = std::filesystem::exists(out_csv);
  const bool empty  = (!exists || std::filesystem::file_size(out_csv) == 0);
  if (!empty) return true;

  std::ofstream fout(out_csv, std::ios::out);
  if (!fout) {
    std::cerr << "[ERROR] Cannot create output CSV: " << out_csv << "\n";
    return false;
  }

  fout
    << "run,dp_idx,dp_label,dp_lo,dp_hi,setting,coeff_tag,"
    << "W_obs,Em_obs,Pmz_obs,Pmy_obs,"
    << "W_err,Em_err,Pmz_err,Pmy_err,"
    << "de0_ref,dpe_ref_input,dthe_ref,dthp_ref,dpp_ref,"
    << "W_pred_ref,Em_pred_ref,Pmz_pred_ref,Pmy_pred_ref,"
    << "r_W,r_Em,r_Pmz,r_Pmy,SSE_ref,"
    << "large_dthe_warn,large_dthp_warn,large_dpp_warn,poor_ref_consistency_warn,notes\n";
  return true;
}

bool AppendCSVRow(const RefResult& r, const std::string& out_csv) {
  if (!WriteCSVHeaderIfNeeded(out_csv)) return false;
  std::ofstream fout(out_csv, std::ios::app);
  if (!fout) {
    std::cerr << "[ERROR] Cannot append to output CSV: " << out_csv << "\n";
    return false;
  }

  fout << r.run << ","
       << r.dp_idx << ","
       << r.dp_label << ","
       << std::fixed << std::setprecision(3)
       << r.dp_lo << ","
       << r.dp_hi << ","
       << r.setting << ","
       << r.coeff_tag << ","
       << std::setprecision(6)
       << r.W_obs << ","
       << r.Em_obs << ","
       << r.Pmz_obs << ","
       << r.Pmy_obs << ","
       << r.W_err << ","
       << r.Em_err << ","
       << r.Pmz_err << ","
       << r.Pmy_err << ","
       << r.de0_ref << ","
       << r.dpe_ref_input << ","
       << r.dthe_ref << ","
       << r.dthp_ref << ","
       << r.dpp_ref << ","
       << r.W_pred_ref << ","
       << r.Em_pred_ref << ","
       << r.Pmz_pred_ref << ","
       << r.Pmy_pred_ref << ","
       << r.r_W << ","
       << r.r_Em << ","
       << r.r_Pmz << ","
       << r.r_Pmy << ","
       << r.SSE_ref << ","
       << r.large_dthe_warn << ","
       << r.large_dthp_warn << ","
       << r.large_dpp_warn << ","
       << r.poor_ref_consistency_warn << ","
       << r.notes << "\n";
  return true;
}

void PrintSummary(const RefResult& r) {
  std::cout << "run             : " << r.run << "\n";
  std::cout << "setting         : " << r.setting << " (coeff_tag=" << r.coeff_tag << ")\n";
  std::cout << "dp_idx          : " << r.dp_idx << "\n";
  std::cout << "dp_label        : " << r.dp_label << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "dp_lo, dp_hi    : " << r.dp_lo << ", " << r.dp_hi << "\n";
  std::cout << std::setprecision(6);
  std::cout << "W_obs           : " << r.W_obs << " MeV\n";
  std::cout << "Em_obs          : " << r.Em_obs << " MeV\n";
  std::cout << "Pmz_obs         : " << r.Pmz_obs << " MeV\n";
  std::cout << "Pmy_obs         : " << r.Pmy_obs << " MeV\n";
  std::cout << "dpe_ref_input   : " << r.dpe_ref_input << " [0.1%]\n";
  std::cout << "dthe_ref        : " << r.dthe_ref << " [mrad]\n";
  std::cout << "dthp_ref        : " << r.dthp_ref << " [mrad]\n";
  std::cout << "dpp_ref         : " << r.dpp_ref << " [0.1%]\n";
  std::cout << "W_pred_ref      : " << r.W_pred_ref << " MeV\n";
  std::cout << "Em_pred_ref     : " << r.Em_pred_ref << " MeV\n";
  std::cout << "Pmz_pred_ref    : " << r.Pmz_pred_ref << " MeV\n";
  std::cout << "Pmy_pred_ref    : " << r.Pmy_pred_ref << " MeV\n";
  std::cout << "r_W             : " << r.r_W << " MeV\n";
  std::cout << "r_Em            : " << r.r_Em << " MeV\n";
  std::cout << "r_Pmz           : " << r.r_Pmz << " MeV\n";
  std::cout << "r_Pmy           : " << r.r_Pmy << " MeV\n";
  std::cout << "SSE_ref         : " << r.SSE_ref << "\n";
  std::cout << "warn_large_dthe : " << r.large_dthe_warn << "\n";
  std::cout << "warn_large_dthp : " << r.large_dthp_warn << "\n";
  std::cout << "warn_large_dpp  : " << r.large_dpp_warn << "\n";
  std::cout << "warn_poor_ref   : " << r.poor_ref_consistency_warn << "\n";
  if (!r.notes.empty()) std::cout << "notes           : " << r.notes << "\n";
}

std::vector<int> IntersectRuns(const std::map<int, ObsRunData>& obs,
                               const std::map<int, DpeOnlyRow>& dpe) {
  std::vector<int> runs;
  for (const auto& kv : obs) {
    if (dpe.count(kv.first)) runs.push_back(kv.first);
  }
  std::sort(runs.begin(), runs.end());
  return runs;
}

bool BuildReferenceOne(int run,
                       const std::string& offsets_csv,
                       const std::string& dpeonly_csv,
                       const std::string& out_csv,
                       bool append_csv = true) {
  std::map<int, ObsRunData> obs;
  std::map<int, DpeOnlyRow> dpe;
  if (!LoadObservedOffsets(offsets_csv, obs)) return false;
  if (!LoadDpeOnly(dpeonly_csv, dpe)) return false;

  if (!obs.count(run)) {
    std::cerr << "[ERROR] Run " << run << " not found in observed-offset table.\n";
    return false;
  }
  if (!dpe.count(run)) {
    std::cerr << "[ERROR] Run " << run << " not found in dpeOnly table.\n";
    return false;
  }
  if (!HasAllObservables(obs[run])) {
    std::cerr << "[ERROR] Run " << run << " is missing one or more full-bin observables among W, Em, Pmz, Pmy.\n";
    return false;
  }

  std::cout << "\n============================================================\n";
  std::cout << "HEEPCHECK REFERENCE FROM DPEONLY\n";
  std::cout << "reference model  : de0_ref=0, dpe_ref from dpeOnly(full-bin)\n";
  std::cout << "construction     : dthe from W, dthp from Pmy, dpp from LSQ(Em,Pmz)\n";

  RefResult r = BuildReferenceForRun(run, obs[run], dpe[run]);
  if (!r.ok) {
    std::cerr << "[ERROR] Failed to build reference for run=" << run
              << " notes='" << r.notes << "'\n";
    return false;
  }
  PrintSummary(r);

  if (append_csv) {
    if (!AppendCSVRow(r, out_csv)) return false;
    std::cout << "[INFO] Appended result to: " << out_csv << "\n";
  }

  std::cout << "============================================================\n\n";
  return true;
}

bool BuildReferenceAll(const std::string& offsets_csv,
                       const std::string& dpeonly_csv,
                       const std::string& out_csv) {
  std::map<int, ObsRunData> obs;
  std::map<int, DpeOnlyRow> dpe;
  if (!LoadObservedOffsets(offsets_csv, obs)) return false;
  if (!LoadDpeOnly(dpeonly_csv, dpe)) return false;

  const auto runs = IntersectRuns(obs, dpe);
  std::cout << "[INFO] Found " << runs.size() << " runs in the intersection of observed-offset and dpeOnly tables.\n";
  bool all_ok = true;
  for (size_t i = 0; i < runs.size(); ++i) {
    const int run = runs[i];
    std::cout << "[INFO] Processing " << (i+1) << "/" << runs.size() << " : run=" << run << "\n";
    if (!HasAllObservables(obs[run])) {
      all_ok = false;
      std::cerr << "[WARN] Skipping run=" << run << " because one or more of W, Em, Pmz, Pmy is missing in full-bin observed offsets.\n";
      continue;
    }
    RefResult r = BuildReferenceForRun(run, obs[run], dpe[run]);
    if (!r.ok) {
      all_ok = false;
      std::cerr << "[WARN] Failed for run=" << run << " notes='" << r.notes << "'\n";
      continue;
    }
    if (!AppendCSVRow(r, out_csv)) {
      all_ok = false;
      std::cerr << "[WARN] Failed to write output for run=" << run << "\n";
    }
  }
  return all_ok;
}

} // namespace

void heepcheck_reference_from_dpeOnly(
    const char* mode = "one",
    int run = -1,
    const char* offsets_csv = "results/tables/heep_dp_binned_offsets.csv",
    const char* dpeonly_csv = "results/tables/heepcheck_inversion_dpeOnly.csv",
    const char* out_csv     = "results/tables/heepcheck_reference_from_dpeOnly.csv")
{
  const std::string smode = Trim(mode ? mode : "one");

  if (smode == "one") {
    if (run < 0) {
      std::cerr << "[ERROR] For mode='one', please provide a valid run.\n";
      return;
    }
    BuildReferenceOne(run, offsets_csv, dpeonly_csv, out_csv, true);
    return;
  }

  if (smode == "all") {
    BuildReferenceAll(offsets_csv, dpeonly_csv, out_csv);
    return;
  }

  std::cerr << "[ERROR] Unknown mode: " << smode << "\n"
            << "        Use mode='one' or mode='all'.\n";
}
