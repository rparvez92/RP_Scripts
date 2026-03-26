// File: heep_check_v0/macros/heepcheck_inversion_full4.C
//
// Purpose:
//   Invert the HEEP-check linear response (with de0 fixed to 0) to solve for
//     dthe, dpe, dthp, dpp
//   from measured offsets in
//     W, Em, Pmz, Pmx
//
// Conventions:
//   - measured offsets are data - sim  (from dmu_MeV in the offsets CSV)
//   - e0 is taken from ebeam in the bigtable, converted GeV -> MeV
//   - the0 is taken from shms_th in degrees
//   - parameter units match heepcheck_sign.f exactly:
//       dthe : 1 mrad units
//       dpe  : 0.1% units
//       dthp : 1 mrad units
//       dpp  : 0.1% units
//
// First interface:
//   HeepCheckSolveOne(run, dp_idx)
//
// Optional batch wrapper:
//   HeepCheckSolveAll()
//
// Suggested usage from heep_check_v0/ directory:
//   root -l -b -q 'macros/heepcheck_inversion_full4.C("one",23839,0)'
//   root -l -b -q 'macros/heepcheck_inversion_full4.C("all",-1,-1)'

#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h>
#include <TMath.h>

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
#include <tuple>
#include <vector>

namespace {

constexpr double kAm     = 938.30;       // MeV
constexpr double kRadFac = 0.017453293;  // deg -> rad

struct BigRow {
  int    run = -1;
  double e0_MeV = 0.0;   // from ebeam [GeV] * 1000
  double the0_deg = 0.0; // from shms_th [deg]
};

struct OffsetRow {
  int         run = -1;
  int         dp_idx = -999;
  std::string dp_label;
  double      dp_lo = 0.0;
  double      dp_hi = 0.0;
  std::string var;
  double      dmu_MeV = 0.0;
  double      dmu_err_MeV = 0.0;
};

struct SolveInput {
  int    run = -1;
  int    dp_idx = -999;
  double e0_MeV = 0.0;
  double the0_deg = 0.0;

  std::string dp_label;
  double      dp_lo = 0.0;
  double      dp_hi = 0.0;

  double dW_MeV = 0.0;
  double dW_err_MeV = 0.0;

  double dEm_MeV = 0.0;
  double dEm_err_MeV = 0.0;

  double dPmz_MeV = 0.0;
  double dPmz_err_MeV = 0.0;

  double dPmx_MeV = 0.0;
  double dPmx_err_MeV = 0.0;
};

struct SolveResult {
  bool ok = false;

  int    run = -1;
  int    dp_idx = -999;
  std::string dp_label;
  double dp_lo = 0.0;
  double dp_hi = 0.0;

  double e0_MeV = 0.0;
  double the0_deg = 0.0;

  double pe0_MeV = 0.0;
  double q0_MeV  = 0.0;
  double thq0_deg = 0.0;
  double W0_MeV   = 0.0;

  double dW_MeV = 0.0;
  double dW_err_MeV = 0.0;
  double dEm_MeV = 0.0;
  double dEm_err_MeV = 0.0;
  double dPmz_MeV = 0.0;
  double dPmz_err_MeV = 0.0;
  double dPmx_MeV = 0.0;
  double dPmx_err_MeV = 0.0;

  // solved parameters in heepcheck units
  double dthe = 0.0; // 1 mrad units
  double dpe  = 0.0; // 0.1% units
  double dthp = 0.0; // 1 mrad units
  double dpp  = 0.0; // 0.1% units

  // closure
  double pred_dW_MeV   = 0.0;
  double pred_dEm_MeV  = 0.0;
  double pred_dPmz_MeV = 0.0;
  double pred_dPmx_MeV = 0.0;

  double res_dW_MeV   = 0.0;
  double res_dEm_MeV  = 0.0;
  double res_dPmz_MeV = 0.0;
  double res_dPmx_MeV = 0.0;

  double det = 0.0;
  int    matrix_ok = 0;

  TMatrixD D;     // 4x4
  TVectorD y;     // measured
  TVectorD x;     // solved
  TVectorD ypred; // predicted

  SolveResult() : D(4,4), y(4), x(4), ypred(4) {}
};

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

bool IsInteger(const std::string& s) {
  if (s.empty()) return false;
  char* end = nullptr;
  errno = 0;
  std::strtol(s.c_str(), &end, 10);
  return errno == 0 && end && *end == '\0';
}

bool ToInt(const std::string& s, int& val) {
  if (!IsInteger(s)) return false;
  val = std::stoi(s);
  return true;
}

bool ToDouble(const std::string& s, double& val) {
  if (s.empty()) return false;
  char* end = nullptr;
  errno = 0;
  val = std::strtod(s.c_str(), &end);
  return errno == 0 && end && *end == '\0';
}

double side3(double x, double y, double t_deg) {
  double arg = x*x + y*y - 2.0*x*y*std::cos(t_deg * kRadFac);
  if (arg < 0.0) arg = 0.0;
  return std::sqrt(arg);
}

double angleDeg(double x, double y, double z) {
  double arg = (y*y + z*z - x*x) / (2.0*y*z);
  if (arg < -1.0) arg = -1.0;
  if (arg >  1.0) arg =  1.0;
  return std::acos(arg) / kRadFac;
}

void heepkin(double e0_MeV, double the0_deg,
             double& pe0_MeV, double& thq0_deg, double& q0_MeV) {
  pe0_MeV = e0_MeV /
            (1.0 + 2.0 * e0_MeV / kAm *
             std::pow(std::sin(the0_deg * kRadFac / 2.0), 2));
  q0_MeV   = side3(e0_MeV, pe0_MeV, the0_deg);
  thq0_deg = angleDeg(pe0_MeV, q0_MeV, e0_MeV);
}

double Wfunc(double e, double e1, double th_deg) {
  return std::sqrt(2.0 * kAm * (e - e1) + kAm*kAm
                   - 4.0 * e * e1 * std::pow(std::sin(th_deg * kRadFac / 2.0), 2));
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

bool LoadBigtable(const std::string& path, std::map<int, BigRow>& out) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open bigtable CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Bigtable CSV is empty: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);

  for (const auto& need : {"run", "ebeam", "shms_th"}) {
    if (!idx.count(need)) {
      std::cerr << "[ERROR] Missing required column in bigtable: " << need << "\n";
      return false;
    }
  }

  int lineno = 1;
  while (std::getline(fin, line)) {
    ++lineno;
    if (Trim(line).empty()) continue;

    const auto tok = SplitCSVSimple(line);
    if ((int)tok.size() < (int)hdr.size()) continue;

    int run = -1;
    double ebeam_GeV = 0.0;
    double shms_th_deg = 0.0;

    if (!ToInt(tok[idx.at("run")], run)) continue;
    if (!ToDouble(tok[idx.at("ebeam")], ebeam_GeV)) continue;
    if (!ToDouble(tok[idx.at("shms_th")], shms_th_deg)) continue;

    if (out.count(run)) {
      std::cerr << "[ERROR] Duplicate run in bigtable: " << run << "\n";
      return false;
    }

    BigRow r;
    r.run = run;
    r.e0_MeV = ebeam_GeV * 1000.0;
    r.the0_deg = shms_th_deg;
    out[run] = r;
  }

  return true;
}

bool LoadOffsets(const std::string& path, std::vector<OffsetRow>& out) {
  std::ifstream fin(path);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open offsets CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Offsets CSV is empty: " << path << "\n";
    return false;
  }

  const auto hdr = SplitCSVSimple(line);
  const auto idx = HeaderIndex(hdr);

  for (const auto& need : {"run", "dp_idx", "dp_label", "dp_lo", "dp_hi",
                           "var", "dmu_MeV", "dmu_err_MeV"}) {
    if (!idx.count(need)) {
      std::cerr << "[ERROR] Missing required column in offsets CSV: " << need << "\n";
      return false;
    }
  }

  int lineno = 1;
  while (std::getline(fin, line)) {
    ++lineno;
    if (Trim(line).empty()) continue;

    const auto tok = SplitCSVSimple(line);
    if ((int)tok.size() < (int)hdr.size()) continue;

    OffsetRow r;

    if (!ToInt(tok[idx.at("run")], r.run)) continue;
    if (!ToInt(tok[idx.at("dp_idx")], r.dp_idx)) continue;
    r.dp_label = tok[idx.at("dp_label")];
    if (!ToDouble(tok[idx.at("dp_lo")], r.dp_lo)) continue;
    if (!ToDouble(tok[idx.at("dp_hi")], r.dp_hi)) continue;
    r.var = tok[idx.at("var")];
    if (!ToDouble(tok[idx.at("dmu_MeV")], r.dmu_MeV)) continue;
    if (!ToDouble(tok[idx.at("dmu_err_MeV")], r.dmu_err_MeV)) continue;

    out.push_back(r);
  }

  return true;
}

bool BuildSolveInput(int run, int dp_idx,
                     const std::map<int, BigRow>& big,
                     const std::vector<OffsetRow>& offs,
                     SolveInput& in) {
  auto itB = big.find(run);
  if (itB == big.end()) {
    std::cerr << "[ERROR] Run " << run << " not found in bigtable.\n";
    return false;
  }

  std::map<std::string, OffsetRow> m;
  int nrows = 0;

  for (const auto& r : offs) {
    if (r.run != run || r.dp_idx != dp_idx) continue;
    ++nrows;

    if (m.count(r.var)) {
      std::cerr << "[ERROR] Duplicate offsets row for run=" << run
                << ", dp_idx=" << dp_idx << ", var=" << r.var << "\n";
      return false;
    }
    m[r.var] = r;
  }

  for (const auto& v : {"W", "Em", "Pmz", "Pmx"}) {
    if (!m.count(v)) {
      std::cerr << "[ERROR] Missing var=" << v
                << " for run=" << run << ", dp_idx=" << dp_idx << "\n";
      return false;
    }
  }

  in.run = run;
  in.dp_idx = dp_idx;
  in.e0_MeV = itB->second.e0_MeV;
  in.the0_deg = itB->second.the0_deg;

  in.dp_label = m["W"].dp_label;
  in.dp_lo = m["W"].dp_lo;
  in.dp_hi = m["W"].dp_hi;

  auto sameBinMeta = [&](const OffsetRow& r) {
    return (r.dp_label == in.dp_label &&
            std::abs(r.dp_lo - in.dp_lo) < 1e-12 &&
            std::abs(r.dp_hi - in.dp_hi) < 1e-12);
  };

  for (const auto& v : {"Em", "Pmz", "Pmx"}) {
    if (!sameBinMeta(m[v])) {
      std::cerr << "[ERROR] Inconsistent dp metadata across vars for run=" << run
                << ", dp_idx=" << dp_idx << "\n";
      return false;
    }
  }

  in.dW_MeV      = m["W"].dmu_MeV;
  in.dW_err_MeV  = m["W"].dmu_err_MeV;
  in.dEm_MeV     = m["Em"].dmu_MeV;
  in.dEm_err_MeV = m["Em"].dmu_err_MeV;
  in.dPmz_MeV    = m["Pmz"].dmu_MeV;
  in.dPmz_err_MeV= m["Pmz"].dmu_err_MeV;
  in.dPmx_MeV    = m["Pmx"].dmu_MeV;
  in.dPmx_err_MeV= m["Pmx"].dmu_err_MeV;

  return true;
}

void FillMatrixFull4(double e0, double the0,
                     double& pe0, double& q0, double& thq0, double& W0,
                     TMatrixD& D) {
  heepkin(e0, the0, pe0, thq0, q0);
  W0 = Wfunc(e0, pe0, the0);

  D.ResizeTo(4,4);
  D.Zero();

  // Column ordering:
  // 0: dthe [1 mrad]
  // 1: dpe  [0.1%]
  // 2: dthp [1 mrad]
  // 3: dpp  [0.1%]

  // Row ordering:
  // 0: W
  // 1: Em
  // 2: Pmz
  // 3: Pmx

  // --- dthe column (variation of +1 mrad in electron angle)
  {
    const double dth = 0.001; // rad = 1 mrad
    D(0,0) = Wfunc(e0, pe0, the0 + dth / kRadFac) - W0;
    D(1,0) = 0.0;
    D(2,0) = -pe0 * dth * std::sin((the0 + thq0) * kRadFac);
    D(3,0) =  pe0 * dth * std::cos((the0 + thq0) * kRadFac);
  }

  // --- dpe column (variation of +0.1% in scattered electron energy)
  {
    const double dpe_phys = 0.001 * pe0;
    D(0,1) = Wfunc(e0, pe0 + dpe_phys, the0) - W0;
    D(1,1) = -dpe_phys;
    D(2,1) =  dpe_phys * std::cos((the0 + thq0) * kRadFac);
    D(3,1) =  dpe_phys * std::sin((the0 + thq0) * kRadFac);
  }

  // --- dthp column (variation of +1 mrad in proton angle)
  {
    const double dthp = 0.001; // rad = 1 mrad
    D(0,2) = 0.0;
    D(1,2) = 0.0;
    D(2,2) = 0.0;
    D(3,2) = -q0 * dthp;
  }

  // --- dpp column (variation of +0.1% in proton momentum)
  {
    const double dpp_phys = 0.001 * q0;
    D(0,3) = 0.0;
    D(1,3) = -q0 / std::sqrt(kAm*kAm + q0*q0) * dpp_phys;
    D(2,3) = dpp_phys;
    D(3,3) = 0.0;
  }
}

void PrintMatrix(const TMatrixD& D) {
  std::cout << "\nMatrix D (rows: W, Em, Pmz, Pmx ; cols: dthe, dpe, dthp, dpp)\n";
  std::cout << "Units: dthe,dthp in 1 mrad ; dpe,dpp in 0.1%\n\n";

  std::cout << std::setw(12) << ""
            << std::setw(14) << "dthe"
            << std::setw(14) << "dpe"
            << std::setw(14) << "dthp"
            << std::setw(14) << "dpp" << "\n";

  const char* rows[4] = {"W", "Em", "Pmz", "Pmx"};
  std::cout << std::fixed << std::setprecision(6);
  for (int i = 0; i < 4; ++i) {
    std::cout << std::setw(12) << rows[i];
    for (int j = 0; j < 4; ++j) {
      std::cout << std::setw(14) << D(i,j);
    }
    std::cout << "\n";
  }
}

void PrintVectorMeasured(const SolveInput& in) {
  std::cout << "\nMeasured offsets used in inversion (data - sim) [MeV]\n";
  std::cout << std::fixed << std::setprecision(6)
            << "  dW   = " << in.dW_MeV   << " +/- " << in.dW_err_MeV   << "\n"
            << "  dEm  = " << in.dEm_MeV  << " +/- " << in.dEm_err_MeV  << "\n"
            << "  dPmz = " << in.dPmz_MeV << " +/- " << in.dPmz_err_MeV << "\n"
            << "  dPmx = " << in.dPmx_MeV << " +/- " << in.dPmx_err_MeV << "\n";
}

void PrintSolveSummary(const SolveResult& r) {
  std::cout << "\nSolved parameters (heepcheck units)\n";
  std::cout << std::fixed << std::setprecision(6)
            << "  dthe = " << r.dthe << "   [1 mrad]\n"
            << "  dpe  = " << r.dpe  << "   [0.1%]\n"
            << "  dthp = " << r.dthp << "   [1 mrad]\n"
            << "  dpp  = " << r.dpp  << "   [0.1%]\n";

  std::cout << "\nClosure check [MeV]\n";
  std::cout << std::fixed << std::setprecision(6)
            << "               measured      predicted      residual\n"
            << "  dW   : " << std::setw(12) << r.dW_MeV
            << "  "     << std::setw(12) << r.pred_dW_MeV
            << "  "     << std::setw(12) << r.res_dW_MeV   << "\n"
            << "  dEm  : " << std::setw(12) << r.dEm_MeV
            << "  "     << std::setw(12) << r.pred_dEm_MeV
            << "  "     << std::setw(12) << r.res_dEm_MeV  << "\n"
            << "  dPmz : " << std::setw(12) << r.dPmz_MeV
            << "  "     << std::setw(12) << r.pred_dPmz_MeV
            << "  "     << std::setw(12) << r.res_dPmz_MeV << "\n"
            << "  dPmx : " << std::setw(12) << r.dPmx_MeV
            << "  "     << std::setw(12) << r.pred_dPmx_MeV
            << "  "     << std::setw(12) << r.res_dPmx_MeV << "\n";

  std::cout << "\nMatrix diagnostics\n"
            << "  determinant = " << r.det << "\n"
            << "  matrix_ok   = " << r.matrix_ok << "\n";
}

bool WriteCSVRow(const SolveResult& r, const std::string& outcsv) {
  if (!EnsurePathForFile(outcsv)) {
    std::cerr << "[ERROR] Could not create parent directory for: " << outcsv << "\n";
    return false;
  }

  const bool needHeader = !std::filesystem::exists(outcsv) ||
                          std::filesystem::file_size(outcsv) == 0;

  std::ofstream fout(outcsv, std::ios::app);
  if (!fout.is_open()) {
    std::cerr << "[ERROR] Cannot open output CSV for append: " << outcsv << "\n";
    return false;
  }

  if (needHeader) {
    fout << "fit_model,de0_fixed,sign_convention,"
         << "run,dp_idx,dp_label,dp_lo,dp_hi,"
         << "e0_MeV,the0_deg,pe0_MeV,q0_MeV,thq0_deg,W0_MeV,"
         << "dW_MeV,dW_err_MeV,dEm_MeV,dEm_err_MeV,dPmz_MeV,dPmz_err_MeV,dPmx_MeV,dPmx_err_MeV,"
         << "dthe,dpe,dthp,dpp,"
         << "pred_dW_MeV,pred_dEm_MeV,pred_dPmz_MeV,pred_dPmx_MeV,"
         << "res_dW_MeV,res_dEm_MeV,res_dPmz_MeV,res_dPmx_MeV,"
         << "det,matrix_ok\n";
  }

  fout << "full4,1,data_minus_sim,"
       << r.run << ","
       << r.dp_idx << ","
       << r.dp_label << ","
       << std::setprecision(15) << r.dp_lo << ","
       << std::setprecision(15) << r.dp_hi << ","
       << std::setprecision(15) << r.e0_MeV << ","
       << std::setprecision(15) << r.the0_deg << ","
       << std::setprecision(15) << r.pe0_MeV << ","
       << std::setprecision(15) << r.q0_MeV << ","
       << std::setprecision(15) << r.thq0_deg << ","
       << std::setprecision(15) << r.W0_MeV << ","
       << std::setprecision(15) << r.dW_MeV << ","
       << std::setprecision(15) << r.dW_err_MeV << ","
       << std::setprecision(15) << r.dEm_MeV << ","
       << std::setprecision(15) << r.dEm_err_MeV << ","
       << std::setprecision(15) << r.dPmz_MeV << ","
       << std::setprecision(15) << r.dPmz_err_MeV << ","
       << std::setprecision(15) << r.dPmx_MeV << ","
       << std::setprecision(15) << r.dPmx_err_MeV << ","
       << std::setprecision(15) << r.dthe << ","
       << std::setprecision(15) << r.dpe << ","
       << std::setprecision(15) << r.dthp << ","
       << std::setprecision(15) << r.dpp << ","
       << std::setprecision(15) << r.pred_dW_MeV << ","
       << std::setprecision(15) << r.pred_dEm_MeV << ","
       << std::setprecision(15) << r.pred_dPmz_MeV << ","
       << std::setprecision(15) << r.pred_dPmx_MeV << ","
       << std::setprecision(15) << r.res_dW_MeV << ","
       << std::setprecision(15) << r.res_dEm_MeV << ","
       << std::setprecision(15) << r.res_dPmz_MeV << ","
       << std::setprecision(15) << r.res_dPmx_MeV << ","
       << std::setprecision(15) << r.det << ","
       << r.matrix_ok << "\n";

  return true;
}

SolveResult SolveOneInternal(const SolveInput& in) {
  SolveResult r;
  r.run = in.run;
  r.dp_idx = in.dp_idx;
  r.dp_label = in.dp_label;
  r.dp_lo = in.dp_lo;
  r.dp_hi = in.dp_hi;
  r.e0_MeV = in.e0_MeV;
  r.the0_deg = in.the0_deg;

  r.dW_MeV = in.dW_MeV;
  r.dW_err_MeV = in.dW_err_MeV;
  r.dEm_MeV = in.dEm_MeV;
  r.dEm_err_MeV = in.dEm_err_MeV;
  r.dPmz_MeV = in.dPmz_MeV;
  r.dPmz_err_MeV = in.dPmz_err_MeV;
  r.dPmx_MeV = in.dPmx_MeV;
  r.dPmx_err_MeV = in.dPmx_err_MeV;

  FillMatrixFull4(in.e0_MeV, in.the0_deg,
                  r.pe0_MeV, r.q0_MeV, r.thq0_deg, r.W0_MeV, r.D);

  r.y(0) = in.dW_MeV;
  r.y(1) = in.dEm_MeV;
  r.y(2) = in.dPmz_MeV;
  r.y(3) = in.dPmx_MeV;

  TDecompLU lu(r.D);
  Double_t d1 = 0.0, d2 = 0.0;
  r.matrix_ok = lu.Decompose() ? 1 : 0;
  lu.Det(d1, d2);
  r.det = d1 * std::pow(2.0, d2);

  if (!r.matrix_ok) {
    std::cerr << "[ERROR] Matrix decomposition failed for run=" << in.run
              << ", dp_idx=" << in.dp_idx << "\n";
    r.ok = false;
    return r;
  }

  Bool_t okSolve = kFALSE;
  r.x = lu.Solve(r.y, okSolve);
  if (!okSolve) {
    std::cerr << "[ERROR] Matrix solve failed for run=" << in.run
              << ", dp_idx=" << in.dp_idx << "\n";
    r.ok = false;
    return r;
  }

  r.dthe = r.x(0);
  r.dpe  = r.x(1);
  r.dthp = r.x(2);
  r.dpp  = r.x(3);

  r.ypred = r.D * r.x;

  r.pred_dW_MeV   = r.ypred(0);
  r.pred_dEm_MeV  = r.ypred(1);
  r.pred_dPmz_MeV = r.ypred(2);
  r.pred_dPmx_MeV = r.ypred(3);

  r.res_dW_MeV   = r.dW_MeV   - r.pred_dW_MeV;
  r.res_dEm_MeV  = r.dEm_MeV  - r.pred_dEm_MeV;
  r.res_dPmz_MeV = r.dPmz_MeV - r.pred_dPmz_MeV;
  r.res_dPmx_MeV = r.dPmx_MeV - r.pred_dPmx_MeV;

  r.ok = true;
  return r;
}

std::vector<std::pair<int,int>> UniqueRunDpPairs(const std::vector<OffsetRow>& offs) {
  std::set<std::pair<int,int>> s;
  for (const auto& r : offs) s.insert({r.run, r.dp_idx});
  return std::vector<std::pair<int,int>>(s.begin(), s.end());
}

} // namespace

bool HeepCheckSolveOne(
    int run,
    int dp_idx,
    const std::string& bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const std::string& offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const std::string& out_csv      = "results/tables/heepcheck_inversion_full4.csv",
    bool append_output_csv          = true)
{
  std::map<int, BigRow> big;
  std::vector<OffsetRow> offs;

  if (!LoadBigtable(bigtable_csv, big)) return false;
  if (!LoadOffsets(offsets_csv, offs)) return false;

  SolveInput in;
  if (!BuildSolveInput(run, dp_idx, big, offs, in)) return false;

  std::cout << "\n============================================================\n";
  std::cout << "HEEPCHECK INVERSION FULL4\n";
  std::cout << "fit_model       : full4\n";
  std::cout << "de0_fixed       : 1\n";
  std::cout << "sign_convention : data_minus_sim\n";
  std::cout << "run             : " << in.run << "\n";
  std::cout << "dp_idx          : " << in.dp_idx << "\n";
  std::cout << "dp_label        : " << in.dp_label << "\n";
  std::cout << "dp_lo, dp_hi    : " << in.dp_lo << ", " << in.dp_hi << "\n";
  std::cout << "e0              : " << in.e0_MeV << " MeV\n";
  std::cout << "the0            : " << in.the0_deg << " deg\n";

  SolveResult r = SolveOneInternal(in);

  if (!r.ok) return false;

  std::cout << "pe0             : " << r.pe0_MeV << " MeV\n";
  std::cout << "q0              : " << r.q0_MeV << " MeV\n";
  std::cout << "thq0            : " << r.thq0_deg << " deg\n";
  std::cout << "W0              : " << r.W0_MeV << " MeV\n";

  PrintMatrix(r.D);
  PrintVectorMeasured(in);
  PrintSolveSummary(r);

  if (append_output_csv) {
    if (!WriteCSVRow(r, out_csv)) return false;
    std::cout << "\n[INFO] Appended result to: " << out_csv << "\n";
  }

  std::cout << "============================================================\n\n";
  return true;
}

bool HeepCheckSolveAll(
    const std::string& bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const std::string& offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const std::string& out_csv      = "results/tables/heepcheck_inversion_full4.csv")
{
  std::map<int, BigRow> big;
  std::vector<OffsetRow> offs;

  if (!LoadBigtable(bigtable_csv, big)) return false;
  if (!LoadOffsets(offsets_csv, offs)) return false;

  const auto pairs = UniqueRunDpPairs(offs);

  std::cout << "[INFO] Found " << pairs.size() << " unique (run, dp_idx) pairs.\n";

  bool all_ok = true;
  for (size_t i = 0; i < pairs.size(); ++i) {
    const int run    = pairs[i].first;
    const int dp_idx = pairs[i].second;

    std::cout << "\n[INFO] Processing " << (i+1) << "/" << pairs.size()
              << " : run=" << run << ", dp_idx=" << dp_idx << "\n";

    bool ok = HeepCheckSolveOne(run, dp_idx, bigtable_csv, offsets_csv, out_csv, true);
    if (!ok) {
      all_ok = false;
      std::cerr << "[WARN] Failed for run=" << run << ", dp_idx=" << dp_idx << "\n";
    }
  }

  return all_ok;
}

// Dispatcher for convenient ROOT calling.
// Examples:
//   root -l -b -q 'macros/heepcheck_inversion_full4.C("one",23839,0)'
//   root -l -b -q 'macros/heepcheck_inversion_full4.C("all",-1,-1)'
void heepcheck_inversion_full4(
    const char* mode = "one",
    int run = -1,
    int dp_idx = 0,
    const char* bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const char* offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const char* out_csv      = "results/tables/heepcheck_inversion_full4.csv")
{
  const std::string smode = Trim(mode ? mode : "one");

  if (smode == "one") {
    if (run < 0) {
      std::cerr << "[ERROR] For mode='one', please provide a valid run.\n";
      return;
    }
    HeepCheckSolveOne(run, dp_idx, bigtable_csv, offsets_csv, out_csv, true);
    return;
  }

  if (smode == "all") {
    HeepCheckSolveAll(bigtable_csv, offsets_csv, out_csv);
    return;
  }

  std::cerr << "[ERROR] Unknown mode: " << smode << "\n"
            << "        Use mode='one' or mode='all'.\n";
}
