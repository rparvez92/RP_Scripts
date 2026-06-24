// File: heep_check_v0/macros/heepcheck_inversion_dpeDppDthp.C
//
// Purpose:
//   Fit dpe, dthp, and dpp (with de0=0, dthe=0 fixed) using
//   ordinary unweighted least squares with all 4 observables:
//
//     W, Em, Pmz, Pmx
//
// Input:
//   - bigtable/rsidis_bigtable_pass0.csv
//   - results/tables/heep_dp_binned_offsets.csv
//
// Output:
//   - results/tables/heepcheck_inversion_dpeDppDthp.csv
//
// Conventions:
//   - measured offsets are data - sim
//   - e0 = ebeam * 1000 [MeV]
//   - the0 = shms_th [deg]
//   - dpe, dpp are reported in heepcheck units: [0.1%]
//   - dthp is reported in heepcheck units: [1 mrad]
//
// Suggested usage from heep_check_v0/:
//   root -l -b -q 'macros/heepcheck_inversion_dpeDppDthp.C("one",23839,0)'
//   root -l -b -q 'macros/heepcheck_inversion_dpeDppDthp.C("all",-1,-1)'

#include <TDecompLU.h>
#include <TMatrixD.h>
#include <TVectorD.h>

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

constexpr double kAm     = 938.30;
constexpr double kRadFac = 0.017453293;

// Warning thresholds only; they do NOT constrain the solve.
constexpr double kNearSingularRelTol = 1e-12;
constexpr double kWarnAbsDpe         = 25.0;  // [0.1%]
constexpr double kWarnAbsDpp         = 25.0;  // [0.1%]
constexpr double kWarnAbsDthp        = 10.0;  // [1 mrad]

struct BigRow {
  int    run = -1;
  double e0_MeV = 0.0;
  double the0_deg = 0.0;
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
  double q0_MeV = 0.0;
  double thq0_deg = 0.0;
  double W0_MeV = 0.0;

  double dW_MeV = 0.0;
  double dW_err_MeV = 0.0;
  double dEm_MeV = 0.0;
  double dEm_err_MeV = 0.0;
  double dPmz_MeV = 0.0;
  double dPmz_err_MeV = 0.0;
  double dPmx_MeV = 0.0;
  double dPmx_err_MeV = 0.0;

  double dthe = 0.0;
  double dpe  = 0.0;
  double dthp = 0.0;
  double dpp  = 0.0;

  double pred_dW_MeV = 0.0;
  double pred_dEm_MeV = 0.0;
  double pred_dPmz_MeV = 0.0;
  double pred_dPmx_MeV = 0.0;

  double res_dW_MeV = 0.0;
  double res_dEm_MeV = 0.0;
  double res_dPmz_MeV = 0.0;
  double res_dPmx_MeV = 0.0;

  double sse = 0.0;
  double det_AtA = 0.0;
  double maxAbsAtA = 0.0;
  int near_singular_warn = 0;
  int large_param_warn = 0;

  TVectorD y;
  TVectorD c_dpe;
  TVectorD c_dthp;
  TVectorD c_dpp;
  TVectorD ypred;

  TMatrixD A;      // design matrix 4x3
  TMatrixD AtA;    // 3x3 normal matrix
  TVectorD Aty;    // 3-vector RHS
  TVectorD xfit;   // fitted (dpe,dthp,dpp)

  SolveResult()
      : y(4), c_dpe(4), c_dthp(4), c_dpp(4), ypred(4), A(4,3), AtA(3,3), Aty(3), xfit(3) {}
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

  while (std::getline(fin, line)) {
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

  for (const auto& need : {"run","dp_idx","dp_label","dp_lo","dp_hi","var","dmu_MeV","dmu_err_MeV"}) {
    if (!idx.count(need)) {
      std::cerr << "[ERROR] Missing required column in offsets CSV: " << need << "\n";
      return false;
    }
  }

  while (std::getline(fin, line)) {
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

bool SameOffsetRow(const OffsetRow& a, const OffsetRow& b) {
  auto close = [](double x, double y) {
    return std::fabs(x - y) <= 1e-12 * std::max(1.0, std::max(std::fabs(x), std::fabs(y)));
  };
  return a.run == b.run &&
         a.dp_idx == b.dp_idx &&
         a.dp_label == b.dp_label &&
         close(a.dp_lo, b.dp_lo) &&
         close(a.dp_hi, b.dp_hi) &&
         a.var == b.var &&
         close(a.dmu_MeV, b.dmu_MeV) &&
         close(a.dmu_err_MeV, b.dmu_err_MeV);
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
  for (const auto& r : offs) {
    if (r.run != run || r.dp_idx != dp_idx) continue;
    if (!m.count(r.var)) {
      m[r.var] = r;
    } else {
      if (!SameOffsetRow(m[r.var], r)) {
        std::cerr << "[ERROR] Non-identical duplicate offsets row for run=" << run
                  << ", dp_idx=" << dp_idx << ", var=" << r.var << "\n";
        return false;
      }
    }
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
    return r.dp_label == in.dp_label &&
           std::fabs(r.dp_lo - in.dp_lo) < 1e-12 &&
           std::fabs(r.dp_hi - in.dp_hi) < 1e-12;
  };

  for (const auto& v : {"Em", "Pmz", "Pmx"}) {
    if (!sameBinMeta(m[v])) {
      std::cerr << "[ERROR] Inconsistent dp metadata across vars for run=" << run
                << ", dp_idx=" << dp_idx << "\n";
      return false;
    }
  }

  in.dW_MeV = m["W"].dmu_MeV;
  in.dW_err_MeV = m["W"].dmu_err_MeV;
  in.dEm_MeV = m["Em"].dmu_MeV;
  in.dEm_err_MeV = m["Em"].dmu_err_MeV;
  in.dPmz_MeV = m["Pmz"].dmu_MeV;
  in.dPmz_err_MeV = m["Pmz"].dmu_err_MeV;
  in.dPmx_MeV = m["Pmx"].dmu_MeV;
  in.dPmx_err_MeV = m["Pmx"].dmu_err_MeV;

  return true;
}

void FillColumnsDpeDppDthp(double e0, double the0,
                           double& pe0, double& q0, double& thq0, double& W0,
                           TVectorD& c_dpe, TVectorD& c_dthp, TVectorD& c_dpp) {
  heepkin(e0, the0, pe0, thq0, q0);
  W0 = Wfunc(e0, pe0, the0);

  c_dpe.ResizeTo(4);
  c_dthp.ResizeTo(4);
  c_dpp.ResizeTo(4);
  c_dpe.Zero();
  c_dthp.Zero();
  c_dpp.Zero();

  // dpe column: +0.1% scattered electron momentum.
  const double dpe_phys = 0.001 * pe0;
  c_dpe(0) = Wfunc(e0, pe0 + dpe_phys, the0) - W0;
  c_dpe(1) = -dpe_phys;
  c_dpe(2) =  dpe_phys * std::cos((the0 + thq0) * kRadFac);
  c_dpe(3) =  dpe_phys * std::sin((the0 + thq0) * kRadFac);

  // dthp column: +1 mrad proton-angle offset.
  // First-order small-angle model in the q-frame:
  //   - W and Em are unchanged at first order.
  //   - Pmz is unchanged at first order.
  //   - Pmx gains q0 * dtheta.
  const double dthp_phys_rad = 0.001; // 1 mrad
  c_dthp(0) = 0.0;
  c_dthp(1) = 0.0;
  c_dthp(2) = 0.0;
  c_dthp(3) = q0 * dthp_phys_rad;

  // dpp column: +0.1% proton momentum.
  const double dpp_phys = 0.001 * q0;
  c_dpp(0) = 0.0;
  c_dpp(1) = -q0 / std::sqrt(kAm*kAm + q0*q0) * dpp_phys;
  c_dpp(2) = dpp_phys;
  c_dpp(3) = 0.0;
}

double Determinant3x3(const TMatrixD& M) {
  return M(0,0) * (M(1,1)*M(2,2) - M(1,2)*M(2,1))
       - M(0,1) * (M(1,0)*M(2,2) - M(1,2)*M(2,0))
       + M(0,2) * (M(1,0)*M(2,1) - M(1,1)*M(2,0));
}

double MaxAbsElement(const TMatrixD& M) {
  double out = 0.0;
  for (int i = 0; i < M.GetNrows(); ++i) {
    for (int j = 0; j < M.GetNcols(); ++j) {
      out = std::max(out, std::fabs(M(i,j)));
    }
  }
  return out;
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

  r.y(0) = in.dW_MeV;
  r.y(1) = in.dEm_MeV;
  r.y(2) = in.dPmz_MeV;
  r.y(3) = in.dPmx_MeV;

  FillColumnsDpeDppDthp(in.e0_MeV, in.the0_deg,
                        r.pe0_MeV, r.q0_MeV, r.thq0_deg, r.W0_MeV,
                        r.c_dpe, r.c_dthp, r.c_dpp);

  // Design matrix A = [c_dpe | c_dthp | c_dpp]
  r.A.ResizeTo(4,3);
  for (int i = 0; i < 4; ++i) {
    r.A(i,0) = r.c_dpe(i);
    r.A(i,1) = r.c_dthp(i);
    r.A(i,2) = r.c_dpp(i);
  }

  // Normal equations: (A^T A) x = A^T y
  r.AtA.ResizeTo(3,3);
  r.Aty.ResizeTo(3);
  r.AtA = TMatrixD(TMatrixD::kAtA, r.A);
  r.Aty = r.A.T() * r.y;

  r.det_AtA = Determinant3x3(r.AtA);
  r.maxAbsAtA = MaxAbsElement(r.AtA);
  const double scale3 = std::pow(std::max(1.0, r.maxAbsAtA), 3);
  if (std::fabs(r.det_AtA) < kNearSingularRelTol * scale3) {
    r.near_singular_warn = 1;
    std::cerr << "[WARN] Near-singular normal matrix for run=" << in.run
              << ", dp_idx=" << in.dp_idx
              << " : det(AtA)=" << r.det_AtA
              << ", max|AtA|=" << r.maxAbsAtA << "\n";
  }

  TDecompLU lu(r.AtA);
  if (!lu.Decompose()) {
    std::cerr << "[ERROR] Normal-equation matrix decomposition failed for run=" << in.run
              << ", dp_idx=" << in.dp_idx << "\n";
    r.ok = false;
    return r;
  }

  Bool_t okSolve = kFALSE;
  r.xfit = lu.Solve(r.Aty, okSolve);
  if (!okSolve) {
    std::cerr << "[ERROR] Normal-equation solve failed for run=" << in.run
              << ", dp_idx=" << in.dp_idx << "\n";
    r.ok = false;
    return r;
  }

  r.dthe = 0.0;
  r.dpe  = r.xfit(0);
  r.dthp = r.xfit(1);
  r.dpp  = r.xfit(2);

  if (std::fabs(r.dpe) > kWarnAbsDpe ||
      std::fabs(r.dpp) > kWarnAbsDpp ||
      std::fabs(r.dthp) > kWarnAbsDthp) {
    r.large_param_warn = 1;
    std::cerr << "[WARN] Suspiciously large fitted parameter(s) for run=" << in.run
              << ", dp_idx=" << in.dp_idx
              << " : dpe=" << r.dpe
              << ", dthp=" << r.dthp
              << ", dpp=" << r.dpp << "\n";
  }

  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(r.xfit(i))) {
      std::cerr << "[ERROR] Non-finite fitted parameter for run=" << in.run
                << ", dp_idx=" << in.dp_idx << "\n";
      r.ok = false;
      return r;
    }
  }

  r.ypred.ResizeTo(4);
  for (int i = 0; i < 4; ++i) {
    r.ypred(i) = r.c_dpe(i) * r.dpe + r.c_dthp(i) * r.dthp + r.c_dpp(i) * r.dpp;
  }

  r.pred_dW_MeV   = r.ypred(0);
  r.pred_dEm_MeV  = r.ypred(1);
  r.pred_dPmz_MeV = r.ypred(2);
  r.pred_dPmx_MeV = r.ypred(3);

  r.res_dW_MeV   = r.dW_MeV   - r.pred_dW_MeV;
  r.res_dEm_MeV  = r.dEm_MeV  - r.pred_dEm_MeV;
  r.res_dPmz_MeV = r.dPmz_MeV - r.pred_dPmz_MeV;
  r.res_dPmx_MeV = r.dPmx_MeV - r.pred_dPmx_MeV;

  r.sse = r.res_dW_MeV * r.res_dW_MeV
        + r.res_dEm_MeV * r.res_dEm_MeV
        + r.res_dPmz_MeV * r.res_dPmz_MeV
        + r.res_dPmx_MeV * r.res_dPmx_MeV;

  r.ok = true;
  return r;
}

void PrintColumn(const TVectorD& c, const std::string& name, const std::string& unitText) {
  std::cout << "\n" << name << " response column (rows: W, Em, Pmz, Pmx)\n";
  std::cout << "Units correspond to one unit of " << name << " = " << unitText << "\n";
  std::cout << std::fixed << std::setprecision(6)
            << "  W   : " << c(0) << "\n"
            << "  Em  : " << c(1) << "\n"
            << "  Pmz : " << c(2) << "\n"
            << "  Pmx : " << c(3) << "\n";
}

void PrintMeasured(const SolveInput& in) {
  std::cout << "\nMeasured offsets used in fit (data - sim) [MeV]\n";
  std::cout << std::fixed << std::setprecision(6)
            << "  dW   = " << in.dW_MeV   << " +/- " << in.dW_err_MeV   << "\n"
            << "  dEm  = " << in.dEm_MeV  << " +/- " << in.dEm_err_MeV  << "\n"
            << "  dPmz = " << in.dPmz_MeV << " +/- " << in.dPmz_err_MeV << "\n"
            << "  dPmx = " << in.dPmx_MeV << " +/- " << in.dPmx_err_MeV << "\n";
}

void PrintNormalSystem(const TMatrixD& AtA, const TVectorD& Aty) {
  std::cout << "\nNormal-equation system for (dpe, dthp, dpp)\n";
  std::cout << "  (A^T A) x = A^T y\n\n";
  std::cout << std::fixed << std::setprecision(6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      std::cout << "  AtA(" << i << "," << j << ") = " << AtA(i,j) << "\n";
    }
  }
  for (int i = 0; i < 3; ++i) {
    std::cout << "  Aty(" << i << ")   = " << Aty(i) << "\n";
  }
}

void PrintSummary(const SolveResult& r) {
  std::cout << "\nSolved parameters (heepcheck units)\n";
  std::cout << std::fixed << std::setprecision(6)
            << "  dthe = " << r.dthe << "   [1 mrad] (fixed)\n"
            << "  dpe  = " << r.dpe  << "   [0.1%]   (fit)\n"
            << "  dthp = " << r.dthp << "   [1 mrad] (fit)\n"
            << "  dpp  = " << r.dpp  << "   [0.1%]   (fit)\n";

  std::cout << "\nClosure check [MeV]\n";
  std::cout << std::fixed << std::setprecision(6)
            << "               measured      predicted      residual\n"
            << "  dW   : " << std::setw(12) << r.dW_MeV
            << "  "     << std::setw(12) << r.pred_dW_MeV
            << "  "     << std::setw(12) << r.res_dW_MeV << "\n"
            << "  dEm  : " << std::setw(12) << r.dEm_MeV
            << "  "     << std::setw(12) << r.pred_dEm_MeV
            << "  "     << std::setw(12) << r.res_dEm_MeV << "\n"
            << "  dPmz : " << std::setw(12) << r.dPmz_MeV
            << "  "     << std::setw(12) << r.pred_dPmz_MeV
            << "  "     << std::setw(12) << r.res_dPmz_MeV << "\n"
            << "  dPmx : " << std::setw(12) << r.dPmx_MeV
            << "  "     << std::setw(12) << r.pred_dPmx_MeV
            << "  "     << std::setw(12) << r.res_dPmx_MeV << "\n";

  std::cout << "\nFit diagnostics\n"
            << "  SSE                = " << r.sse << "\n"
            << "  det(AtA)           = " << r.det_AtA << "\n"
            << "  max|AtA|           = " << r.maxAbsAtA << "\n"
            << "  near_singular_warn = " << r.near_singular_warn << "\n"
            << "  large_param_warn   = " << r.large_param_warn << "\n";
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
    std::cerr << "[ERROR] Cannot open output CSV: " << outcsv << "\n";
    return false;
  }

  if (needHeader) {
    fout << "fit_model,de0_fixed,dthe_fixed,dthp_fixed,dpe_fixed,dpp_fixed,sign_convention,"
         << "run,dp_idx,dp_label,dp_lo,dp_hi,"
         << "e0_MeV,the0_deg,pe0_MeV,q0_MeV,thq0_deg,W0_MeV,"
         << "dW_MeV,dW_err_MeV,dEm_MeV,dEm_err_MeV,dPmz_MeV,dPmz_err_MeV,dPmx_MeV,dPmx_err_MeV,"
         << "dthe,dpe,dthp,dpp,"
         << "pred_dW_MeV,pred_dEm_MeV,pred_dPmz_MeV,pred_dPmx_MeV,"
         << "res_dW_MeV,res_dEm_MeV,res_dPmz_MeV,res_dPmx_MeV,"
         << "sse,det_AtA,maxAbsAtA,near_singular_warn,large_param_warn\n";
  }

  fout << "dpeDppDthp,1,1,0,0,0,data_minus_sim,"
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
       << std::setprecision(15) << r.sse << ","
       << std::setprecision(15) << r.det_AtA << ","
       << std::setprecision(15) << r.maxAbsAtA << ","
       << r.near_singular_warn << ","
       << r.large_param_warn << "\n";

  return true;
}

std::vector<std::pair<int,int>> UniqueRunDpPairs(const std::vector<OffsetRow>& offs) {
  std::set<std::pair<int,int>> s;
  for (const auto& r : offs) s.insert({r.run, r.dp_idx});
  return std::vector<std::pair<int,int>>(s.begin(), s.end());
}

} // namespace

bool HeepCheckSolveOneDpeDppDthp(
    int run,
    int dp_idx,
    const std::string& bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const std::string& offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const std::string& out_csv      = "results/tables/heepcheck_inversion_dpeDppDthp.csv",
    bool append_output_csv          = true)
{
  std::map<int, BigRow> big;
  std::vector<OffsetRow> offs;

  if (!LoadBigtable(bigtable_csv, big)) return false;
  if (!LoadOffsets(offsets_csv, offs)) return false;

  SolveInput in;
  if (!BuildSolveInput(run, dp_idx, big, offs, in)) return false;

  std::cout << "\n============================================================\n";
  std::cout << "HEEPCHECK INVERSION DPE-DPP-DTHP\n";
  std::cout << "fit_model       : dpeDppDthp\n";
  std::cout << "de0_fixed       : 1\n";
  std::cout << "dthe_fixed      : 1\n";
  std::cout << "dthp_fixed      : 0\n";
  std::cout << "dpe_fixed       : 0\n";
  std::cout << "dpp_fixed       : 0\n";
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

  PrintColumn(r.c_dpe,  "dpe",  "0.1%");
  PrintColumn(r.c_dthp, "dthp", "1 mrad");
  PrintColumn(r.c_dpp,  "dpp",  "0.1%");
  PrintMeasured(in);
  PrintNormalSystem(r.AtA, r.Aty);
  PrintSummary(r);

  if (append_output_csv) {
    if (!WriteCSVRow(r, out_csv)) return false;
    std::cout << "\n[INFO] Appended result to: " << out_csv << "\n";
  }

  std::cout << "============================================================\n\n";
  return true;
}

bool HeepCheckSolveAllDpeDppDthp(
    const std::string& bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const std::string& offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const std::string& out_csv      = "results/tables/heepcheck_inversion_dpeDppDthp.csv")
{
  std::map<int, BigRow> big;
  std::vector<OffsetRow> offs;

  if (!LoadBigtable(bigtable_csv, big)) return false;
  if (!LoadOffsets(offsets_csv, offs)) return false;

  const auto pairs = UniqueRunDpPairs(offs);
  std::cout << "[INFO] Found " << pairs.size() << " unique (run, dp_idx) pairs.\n";

  bool all_ok = true;
  for (size_t i = 0; i < pairs.size(); ++i) {
    const int run = pairs[i].first;
    const int dp_idx = pairs[i].second;

    std::cout << "\n[INFO] Processing " << (i+1) << "/" << pairs.size()
              << " : run=" << run << ", dp_idx=" << dp_idx << "\n";

    bool ok = HeepCheckSolveOneDpeDppDthp(run, dp_idx, bigtable_csv, offsets_csv, out_csv, true);
    if (!ok) {
      all_ok = false;
      std::cerr << "[WARN] Failed for run=" << run << ", dp_idx=" << dp_idx << "\n";
    }
  }

  return all_ok;
}

void heepcheck_inversion_dpeDppDthp(
    const char* mode = "one",
    int run = -1,
    int dp_idx = 0,
    const char* bigtable_csv = "bigtable/rsidis_bigtable_pass0.csv",
    const char* offsets_csv  = "results/tables/heep_dp_binned_offsets.csv",
    const char* out_csv      = "results/tables/heepcheck_inversion_dpeDppDthp.csv")
{
  const std::string smode = Trim(mode ? mode : "one");

  if (smode == "one") {
    if (run < 0) {
      std::cerr << "[ERROR] For mode='one', please provide a valid run.\n";
      return;
    }
    HeepCheckSolveOneDpeDppDthp(run, dp_idx, bigtable_csv, offsets_csv, out_csv, true);
    return;
  }

  if (smode == "all") {
    HeepCheckSolveAllDpeDppDthp(bigtable_csv, offsets_csv, out_csv);
    return;
  }

  std::cerr << "[ERROR] Unknown mode: " << smode << "\n"
            << "        Use mode='one' or mode='all'.\n";
}
