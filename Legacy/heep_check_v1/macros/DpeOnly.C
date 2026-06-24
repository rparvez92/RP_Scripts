// DpeOnly.C
//
// Purpose:
//   Perform "DpeOnly" fits using GuessedOffsets.csv as input.
//
// Workflow:
//   For each setting (A/B), and for each source bin among {full,b1,b2,b3,b4,b5},
//   we freeze the three guessed offsets from the source bin:
//
//      dthe = dthe_guess(source_bin)
//      dthp = dthp_guess(source_bin)
//      dpp  = dpp_guess(source_bin)
//
//   Then for each target bin among {full,b1,b2,b3,b4,b5}, we fit only dpe by
//   minimizing a weighted chi-square between predicted and measured observable
//   offsets:
//
//      observables = { dW, dEm, dPmz, dPmy }
//
// Important note on the minimization:
//   We fit one free parameter, dpe, by direct chi2 minimization.
//   For each observable i, with fixed (dthe,dthp,dpp), the prediction is:
//
//      pred_i(dpe) = a_i + b_i * dpe
//
//   The weighted chi-square is:
//
//      chi2(dpe) = sum_i [ (y_i - (a_i + b_i dpe))^2 / sigma_i_tot^2 ]
//   and dpe_fit / dpe_fit_err are taken from the minimizer result.
//
// Error model used:
//   The denominator uses a total observable variance:
//
//      sigma_i_tot^2 = sigma_i_meas^2 + sigma_i_pred^2
//
//   where:
//     - sigma_i_meas is the measured observable error from the target bin
//     - sigma_i_pred is the propagated uncertainty from the frozen source-bin
//       guessed offsets (dthe, dthp, dpp), using the fitted source-bin
//       covariance matrix
//
// Assumptions:
//   - Response is linear.
//   - Frozen guessed offsets are treated as uncertain with covariance taken
//     from GuessedOffsets.csv.
//   - No covariance matrix between observables is modeled.
//   - If an observable has invalid / non-positive total uncertainty, it is
//     skipped from that fit.
//   - A fit is declared invalid if fewer than 2 observables remain.
//   - (Setting A, b5) is excluded entirely as both source and target because
//     there is no faithful offset measurement there.
//   - Residuals are reported as predicted - measured to match GuessedOffsets.C.
//
// Input:
//   results/tables/GuessedOffsets.csv
//
// Output:
//   results/tables/DpeOnly.csv
//
// Example:
//   root -l -b -q 'macros/DpeOnly.C("results/tables/GuessedOffsets.csv","results/tables/DpeOnly.csv")'

#include <TMath.h>
#include <TMinuit.h>

#include <algorithm>
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

namespace {

struct Row {
  std::string setting;
  std::string bin;

  double dW_measured   = std::numeric_limits<double>::quiet_NaN();
  double dW_err        = std::numeric_limits<double>::quiet_NaN();
  double dEm_measured  = std::numeric_limits<double>::quiet_NaN();
  double dEm_err       = std::numeric_limits<double>::quiet_NaN();
  double dPmz_measured = std::numeric_limits<double>::quiet_NaN();
  double dPmz_err      = std::numeric_limits<double>::quiet_NaN();
  double dPmy_measured = std::numeric_limits<double>::quiet_NaN();
  double dPmy_err      = std::numeric_limits<double>::quiet_NaN();

  double dthe_guess     = std::numeric_limits<double>::quiet_NaN();
  double dthe_guess_err = std::numeric_limits<double>::quiet_NaN();
  double dpe_guess      = std::numeric_limits<double>::quiet_NaN();
  double dpe_guess_err  = std::numeric_limits<double>::quiet_NaN();
  double dthp_guess     = std::numeric_limits<double>::quiet_NaN();
  double dthp_guess_err = std::numeric_limits<double>::quiet_NaN();
  double dpp_guess      = std::numeric_limits<double>::quiet_NaN();
  double dpp_guess_err  = std::numeric_limits<double>::quiet_NaN();

  double cov_dthe_dthe = std::numeric_limits<double>::quiet_NaN();
  double cov_dthe_dthp = std::numeric_limits<double>::quiet_NaN();
  double cov_dthe_dpp  = std::numeric_limits<double>::quiet_NaN();
  double cov_dthp_dthp = std::numeric_limits<double>::quiet_NaN();
  double cov_dthp_dpp  = std::numeric_limits<double>::quiet_NaN();
  double cov_dpp_dpp   = std::numeric_limits<double>::quiet_NaN();

  int guess_valid = 0;
  int source_guess_trivial = 0;
};

struct ObsModel {
  double a = 0.0;          // fixed part: pred = a + b*dpe
  double b = 0.0;          // coefficient of dpe
  double sigma_meas = std::numeric_limits<double>::quiet_NaN();
  double sigma_pred = std::numeric_limits<double>::quiet_NaN();
  double sigma_tot  = std::numeric_limits<double>::quiet_NaN();
  double y_meas     = std::numeric_limits<double>::quiet_NaN();
  bool used         = false;
};

struct FitContext {
  const ObsModel* obs[4] = {nullptr, nullptr, nullptr, nullptr};
};

static FitContext gFitContext;

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

static std::string NormalizeSetting(const std::string& s) {
  const std::string t = Lower(Trim(s));
  if (t == "a" || t == "setting a" || t == "settinga") return "A";
  if (t == "b" || t == "setting b" || t == "settingb") return "B";
  return Trim(s);
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

static bool IsFinitePos(double x) {
  return std::isfinite(x) && x > 0.0;
}

static bool IsNearlyZero(double x, double eps=1e-14) {
  return std::isfinite(x) && std::fabs(x) < eps;
}

static bool IsExcludedSettingBin(const std::string& setting, const std::string& bin) {
  return setting == "A" && Lower(Trim(bin)) == "b5";
}

static void DpeChi2Fcn(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval,
                       Double_t* par, Int_t /*iflag*/) {
  const double dpe = par[0];
  double chi2 = 0.0;

  for (const ObsModel* M : gFitContext.obs) {
    if (!M || !M->used) continue;
    const double pred = M->a + M->b * dpe;
    const double resid = M->y_meas - pred;
    chi2 += resid * resid / (M->sigma_tot * M->sigma_tot);
  }

  fval = chi2;
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

  const int iSetting  = idx("setting");
  const int iBin      = idx("bin");

  const int iDW       = idx("dw_measured");
  const int iDWErr    = idx("dw_err");
  const int iDEm      = idx("dem_measured");
  const int iDEmErr   = idx("dem_err");
  const int iDPmz     = idx("dpmz_measured");
  const int iDPmzErr  = idx("dpmz_err");
  const int iDPmy     = idx("dpmy_measured");
  const int iDPmyErr  = idx("dpmy_err");

  const int iDthe     = idx("dthe_guess");
  const int iDtheErr  = idx("dthe_guess_err");
  const int iDpe      = idx("dpe_guess");
  const int iDpeErr   = idx("dpe_guess_err");
  const int iDthp     = idx("dthp_guess");
  const int iDthpErr  = idx("dthp_guess_err");
  const int iDpp      = idx("dpp_guess");
  const int iDppErr   = idx("dpp_guess_err");
  const int iCovDtheDthe = idx("cov_dthe_dthe");
  const int iCovDtheDthp = idx("cov_dthe_dthp");
  const int iCovDtheDpp  = idx("cov_dthe_dpp");
  const int iCovDthpDthp = idx("cov_dthp_dthp");
  const int iCovDthpDpp  = idx("cov_dthp_dpp");
  const int iCovDppDpp   = idx("cov_dpp_dpp");
  const int iGuessValid  = idx("guess_valid");

  int iTrivial        = idx("source_guess_trivial");
  if (iTrivial < 0) iTrivial = idx("guess_trivial");

  std::vector<std::string> missing;
  auto need = [&](int i, const std::string& name) {
    if (i < 0) missing.push_back(name);
  };

  need(iSetting, "setting");
  need(iBin, "bin");
  need(iDW, "dW_measured");
  need(iDWErr, "dW_err");
  need(iDEm, "dEm_measured");
  need(iDEmErr, "dEm_err");
  need(iDPmz, "dPmz_measured");
  need(iDPmzErr, "dPmz_err");
  need(iDPmy, "dPmy_measured");
  need(iDPmyErr, "dPmy_err");
  need(iDthe, "dthe_guess");
  need(iDtheErr, "dthe_guess_err");
  need(iDpe, "dpe_guess");
  need(iDpeErr, "dpe_guess_err");
  need(iDthp, "dthp_guess");
  need(iDthpErr, "dthp_guess_err");
  need(iDpp, "dpp_guess");
  need(iDppErr, "dpp_guess_err");
  need(iCovDtheDthe, "cov_dthe_dthe");
  need(iCovDtheDthp, "cov_dthe_dthp");
  need(iCovDtheDpp, "cov_dthe_dpp");
  need(iCovDthpDthp, "cov_dthp_dthp");
  need(iCovDthpDpp, "cov_dthp_dpp");
  need(iCovDppDpp, "cov_dpp_dpp");
  need(iGuessValid, "guess_valid");

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

    r.dW_measured    = getD(iDW);
    r.dW_err         = getD(iDWErr);
    r.dEm_measured   = getD(iDEm);
    r.dEm_err        = getD(iDEmErr);
    r.dPmz_measured  = getD(iDPmz);
    r.dPmz_err       = getD(iDPmzErr);
    r.dPmy_measured  = getD(iDPmy);
    r.dPmy_err       = getD(iDPmyErr);

    r.dthe_guess     = getD(iDthe);
    r.dthe_guess_err = getD(iDtheErr);
    r.dpe_guess      = getD(iDpe);
    r.dpe_guess_err  = getD(iDpeErr);
    r.dthp_guess     = getD(iDthp);
    r.dthp_guess_err = getD(iDthpErr);
    r.dpp_guess      = getD(iDpp);
    r.dpp_guess_err  = getD(iDppErr);
    r.cov_dthe_dthe  = getD(iCovDtheDthe);
    r.cov_dthe_dthp  = getD(iCovDtheDthp);
    r.cov_dthe_dpp   = getD(iCovDtheDpp);
    r.cov_dthp_dthp  = getD(iCovDthpDthp);
    r.cov_dthp_dpp   = getD(iCovDthpDpp);
    r.cov_dpp_dpp    = getD(iCovDppDpp);
    r.guess_valid    = ToIntOrDefault(getS(iGuessValid), 0);

    r.source_guess_trivial = (iTrivial >= 0 && iTrivial < (int)fields.size())
      ? ToIntOrDefault(fields[iTrivial], 0)
      : 0;

    rows.push_back(r);
  }

  return rows;
}

static std::string CsvQuoteIfNeeded(const std::string& s) {
  if (s.find(',') == std::string::npos && s.find('"') == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

static double FrozenVariance(const Row& source,
                             double c_dthe,
                             double c_dthp,
                             double c_dpp) {
  if (!std::isfinite(source.cov_dthe_dthe) ||
      !std::isfinite(source.cov_dthe_dthp) ||
      !std::isfinite(source.cov_dthe_dpp) ||
      !std::isfinite(source.cov_dthp_dthp) ||
      !std::isfinite(source.cov_dthp_dpp) ||
      !std::isfinite(source.cov_dpp_dpp)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double var =
      c_dthe * c_dthe * source.cov_dthe_dthe +
      c_dthp * c_dthp * source.cov_dthp_dthp +
      c_dpp * c_dpp * source.cov_dpp_dpp +
      2.0 * c_dthe * c_dthp * source.cov_dthe_dthp +
      2.0 * c_dthe * c_dpp  * source.cov_dthe_dpp +
      2.0 * c_dthp * c_dpp  * source.cov_dthp_dpp;

  if (!std::isfinite(var)) return std::numeric_limits<double>::quiet_NaN();
  return (var >= 0.0) ? var : (std::fabs(var) < 1e-12 ? 0.0 : std::numeric_limits<double>::quiet_NaN());
}

static void BuildObservableModels(const std::string& setting,
                                  const Row& source,
                                  const Row& target,
                                  ObsModel& W,
                                  ObsModel& Em,
                                  ObsModel& Pmz,
                                  ObsModel& Pmy) {
  const double dthe    = source.dthe_guess;
  const double dthp    = source.dthp_guess;
  const double dpp     = source.dpp_guess;

  if (setting == "A") {
    // Setting A:
    //   dW   = -14.08*dthe - 8.62*dpe
    //   dEm  = -7.06*dpe - 2.10*dpp
    //   dPmz = -5.75*dthe + 4.10*dpe + 2.27*dpp
    //   dPmy =  4.10*dthe + 5.75*dpe - 2.27*dthp

    W.a = -14.08 * dthe;
    W.b = -8.62;
    W.y_meas = target.dW_measured;
    W.sigma_meas = target.dW_err;
    {
      const double var = FrozenVariance(source, -14.08, 0.0, 0.0);
      W.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Em.a = -2.10 * dpp;
    Em.b = -7.06;
    Em.y_meas = target.dEm_measured;
    Em.sigma_meas = target.dEm_err;
    {
      const double var = FrozenVariance(source, 0.0, 0.0, -2.10);
      Em.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Pmz.a = -5.75 * dthe + 2.27 * dpp;
    Pmz.b = 4.10;
    Pmz.y_meas = target.dPmz_measured;
    Pmz.sigma_meas = target.dPmz_err;
    {
      const double var = FrozenVariance(source, -5.75, 0.0, 2.27);
      Pmz.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Pmy.a = 4.10 * dthe - 2.27 * dthp;
    Pmy.b = 5.75;
    Pmy.y_meas = target.dPmy_measured;
    Pmy.sigma_meas = target.dPmy_err;
    {
      const double var = FrozenVariance(source, 4.10, -2.27, 0.0);
      Pmy.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }
  } else {
    // Setting B:
    //   dW   = -17.33*dthe - 8.62*dpe
    //   dEm  = -5.66*dpe - 3.63*dpp
    //   dPmz = -4.30*dthe + 3.69*dpe + 3.74*dpp
    //   dPmy =  3.69*dthe + 4.30*dpe - 3.74*dthp

    W.a = -17.33 * dthe;
    W.b = -8.62;
    W.y_meas = target.dW_measured;
    W.sigma_meas = target.dW_err;
    {
      const double var = FrozenVariance(source, -17.33, 0.0, 0.0);
      W.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Em.a = -3.63 * dpp;
    Em.b = -5.66;
    Em.y_meas = target.dEm_measured;
    Em.sigma_meas = target.dEm_err;
    {
      const double var = FrozenVariance(source, 0.0, 0.0, -3.63);
      Em.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Pmz.a = -4.30 * dthe + 3.74 * dpp;
    Pmz.b = 3.69;
    Pmz.y_meas = target.dPmz_measured;
    Pmz.sigma_meas = target.dPmz_err;
    {
      const double var = FrozenVariance(source, -4.30, 0.0, 3.74);
      Pmz.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }

    Pmy.a = 3.69 * dthe - 3.74 * dthp;
    Pmy.b = 4.30;
    Pmy.y_meas = target.dPmy_measured;
    Pmy.sigma_meas = target.dPmy_err;
    {
      const double var = FrozenVariance(source, 3.69, -3.74, 0.0);
      Pmy.sigma_pred = (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }
  }

  auto finalize = [](ObsModel& M) {
    const double sMeas2 = (std::isfinite(M.sigma_meas) && M.sigma_meas >= 0.0) ? M.sigma_meas * M.sigma_meas : 0.0;
    const double sPred2 = (std::isfinite(M.sigma_pred) && M.sigma_pred >= 0.0) ? M.sigma_pred * M.sigma_pred : 0.0;
    M.sigma_tot = std::sqrt(sMeas2 + sPred2);
    M.used = std::isfinite(M.y_meas) && IsFinitePos(M.sigma_tot) && std::isfinite(M.a) && std::isfinite(M.b);
  };

  finalize(W);
  finalize(Em);
  finalize(Pmz);
  finalize(Pmy);
}

static void WriteHeader(std::ofstream& fout) {
  fout
    << "setting,source_bin,target_bin,"
    << "dthe_fix,dthe_fix_err,dthp_fix,dthp_fix_err,dpp_fix,dpp_fix_err,"
    << "dW_measured,dW_err,dEm_measured,dEm_err,dPmz_measured,dPmz_err,dPmy_measured,dPmy_err,"
    << "dW_pred_err_from_fixed,dEm_pred_err_from_fixed,dPmz_pred_err_from_fixed,dPmy_pred_err_from_fixed,"
    << "dW_sigma_total,dEm_sigma_total,dPmz_sigma_total,dPmy_sigma_total,"
    << "dpe_fit,dpe_fit_err,"
    << "dW_pred,dEm_pred,dPmz_pred,dPmy_pred,"
    << "dW_resid,dW_resid_err,dEm_resid,dEm_resid_err,dPmz_resid,dPmz_resid_err,dPmy_resid,dPmy_resid_err,"
    << "chi2_min,n_obs_used,ndf,chi2_ndf,fit_prob,sse_unweighted,"
    << "source_guess_trivial,fit_valid,fit_note"
    << "\n";
}

static void WriteRow(std::ofstream& fout,
                     const Row& source,
                     const Row& target,
                     const ObsModel& W,
                     const ObsModel& Em,
                     const ObsModel& Pmz,
                     const ObsModel& Pmy,
                     double dpe_fit,
                     double dpe_fit_err,
                     double dW_pred,
                     double dEm_pred,
                     double dPmz_pred,
                     double dPmy_pred,
                     double dW_resid,
                     double dW_resid_err,
                     double dEm_resid,
                     double dEm_resid_err,
                     double dPmz_resid,
                     double dPmz_resid_err,
                     double dPmy_resid,
                     double dPmy_resid_err,
                     double chi2_min,
                     int n_obs_used,
                     int ndf,
                     double chi2_ndf,
                     double fit_prob,
                     double sse_unweighted,
                     int fit_valid,
                     const std::string& fit_note) {
  fout << CsvQuoteIfNeeded(source.setting) << ","
       << CsvQuoteIfNeeded(source.bin) << ","
       << CsvQuoteIfNeeded(target.bin) << ","
       << source.dthe_guess << ","
       << source.dthe_guess_err << ","
       << source.dthp_guess << ","
       << source.dthp_guess_err << ","
       << source.dpp_guess << ","
       << source.dpp_guess_err << ","
       << target.dW_measured << ","
       << target.dW_err << ","
       << target.dEm_measured << ","
       << target.dEm_err << ","
       << target.dPmz_measured << ","
       << target.dPmz_err << ","
       << target.dPmy_measured << ","
       << target.dPmy_err << ","
       << W.sigma_pred << ","
       << Em.sigma_pred << ","
       << Pmz.sigma_pred << ","
       << Pmy.sigma_pred << ","
       << W.sigma_tot << ","
       << Em.sigma_tot << ","
       << Pmz.sigma_tot << ","
       << Pmy.sigma_tot << ","
       << dpe_fit << ","
       << dpe_fit_err << ","
       << dW_pred << ","
       << dEm_pred << ","
       << dPmz_pred << ","
       << dPmy_pred << ","
       << dW_resid << ","
       << dW_resid_err << ","
       << dEm_resid << ","
       << dEm_resid_err << ","
       << dPmz_resid << ","
       << dPmz_resid_err << ","
       << dPmy_resid << ","
       << dPmy_resid_err << ","
       << chi2_min << ","
       << n_obs_used << ","
       << ndf << ","
       << chi2_ndf << ","
       << fit_prob << ","
       << sse_unweighted << ","
       << source.source_guess_trivial << ","
       << fit_valid << ","
       << CsvQuoteIfNeeded(fit_note)
       << "\n";
}

} // namespace

void DpeOnly(const char* inCsv  = "results/tables/GuessedOffsets.csv",
             const char* outCsv = "results/tables/DpeOnly.csv") {
  const auto rows = ReadGuessedOffsets(inCsv);
  if (rows.empty()) {
    std::cerr << "[ERROR] No rows read from " << inCsv << std::endl;
    return;
  }

  std::map<std::string, std::map<std::string, Row>> bySettingBin;
  std::set<std::string> settings;
  std::set<std::string> binsSeen;

  for (const auto& r : rows) {
    bySettingBin[r.setting][r.bin] = r;
    settings.insert(r.setting);
    binsSeen.insert(r.bin);
  }

  const std::vector<std::string> binOrder = {"full","b1","b2","b3","b4","b5"};

  std::ofstream fout(outCsv);
  if (!fout) {
    std::cerr << "[ERROR] Cannot open output CSV for writing: " << outCsv << std::endl;
    return;
  }

  fout << std::setprecision(10);
  WriteHeader(fout);

  int nRowsWritten = 0;
  int nValidFits   = 0;
  int nInvalidFits = 0;

  std::vector<std::string> orderedSettings;
  if (settings.count("A")) orderedSettings.push_back("A");
  if (settings.count("B")) orderedSettings.push_back("B");
  for (const auto& s : settings) {
    if (s != "A" && s != "B") orderedSettings.push_back(s);
  }

  for (const auto& setting : orderedSettings) {
    auto itS = bySettingBin.find(setting);
    if (itS == bySettingBin.end()) continue;

    for (const auto& sourceBin : binOrder) {
      auto itSource = itS->second.find(sourceBin);
      if (itSource == itS->second.end()) continue;
      const Row& source = itSource->second;

      for (const auto& targetBin : binOrder) {
        auto itTarget = itS->second.find(targetBin);
        if (itTarget == itS->second.end()) continue;
        const Row& target = itTarget->second;

        if (IsExcludedSettingBin(setting, source.bin) ||
            IsExcludedSettingBin(setting, target.bin)) {
          continue;
        }

        if (source.guess_valid != 1) {
          continue;
        }

        ObsModel W, Em, Pmz, Pmy;
        BuildObservableModels(setting, source, target, W, Em, Pmz, Pmy);

        std::vector<ObsModel*> obs = {&W, &Em, &Pmz, &Pmy};

        int n_obs_used = 0;
        int n_skipped_invalid_sigma = 0;
        int n_skipped_invalid_y = 0;

        for (auto* M : obs) {
          if (!std::isfinite(M->y_meas)) {
            M->used = false;
            ++n_skipped_invalid_y;
            continue;
          }
          if (!IsFinitePos(M->sigma_tot)) {
            M->used = false;
            ++n_skipped_invalid_sigma;
            continue;
          }
          M->used = true;
          ++n_obs_used;
        }

        double dpe_fit      = std::numeric_limits<double>::quiet_NaN();
        double dpe_fit_err  = std::numeric_limits<double>::quiet_NaN();
        double chi2_min     = std::numeric_limits<double>::quiet_NaN();
        double chi2_ndf     = std::numeric_limits<double>::quiet_NaN();
        double fit_prob     = std::numeric_limits<double>::quiet_NaN();
        double sse_unw      = std::numeric_limits<double>::quiet_NaN();

        double dW_pred      = std::numeric_limits<double>::quiet_NaN();
        double dEm_pred     = std::numeric_limits<double>::quiet_NaN();
        double dPmz_pred    = std::numeric_limits<double>::quiet_NaN();
        double dPmy_pred    = std::numeric_limits<double>::quiet_NaN();

        double dW_resid      = std::numeric_limits<double>::quiet_NaN();
        double dW_resid_err  = std::numeric_limits<double>::quiet_NaN();
        double dEm_resid     = std::numeric_limits<double>::quiet_NaN();
        double dEm_resid_err = std::numeric_limits<double>::quiet_NaN();
        double dPmz_resid    = std::numeric_limits<double>::quiet_NaN();
        double dPmz_resid_err= std::numeric_limits<double>::quiet_NaN();
        double dPmy_resid    = std::numeric_limits<double>::quiet_NaN();
        double dPmy_resid_err= std::numeric_limits<double>::quiet_NaN();

        int ndf = -1;
        int fit_valid = 0;
        std::string fit_note;

        if (n_obs_used < 2) {
          fit_valid = 0;
          fit_note = "fewer than 2 observables with valid total uncertainty";
        } else {
          for (int i = 0; i < 4; ++i) gFitContext.obs[i] = obs[i];
          bool minuit_ok = true;
          {
            TMinuit minuit(1);
            minuit.SetPrintLevel(-1);
            minuit.SetFCN(DpeChi2Fcn);

            Double_t arglist[2];
            Int_t ierflg = 0;
            arglist[0] = 1.0;
            minuit.mnexcm("SET ERR", arglist, 1, ierflg);

            minuit.mnparm(0, "dpe", 0.0, 0.01, 0.0, 0.0, ierflg);
            if (ierflg != 0) {
              fit_valid = 0;
              fit_note = "minuit_parameter_setup_failed";
              minuit_ok = false;
            }

            if (minuit_ok) {
              arglist[0] = 1000.0;
              arglist[1] = 1.0;
              minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
              if (ierflg != 0) {
                fit_valid = 0;
                fit_note = "migrad_failed";
                minuit_ok = false;
              }
            }

            if (minuit_ok) {
              arglist[0] = 200.0;
              minuit.mnexcm("HESSE", arglist, 1, ierflg);
              if (ierflg != 0) {
                fit_valid = 0;
                fit_note = "hesse_failed";
                minuit_ok = false;
              }
            }

            if (minuit_ok) {
              Double_t val = 0.0, err = 0.0;
              minuit.GetParameter(0, val, err);
              dpe_fit = val;
              dpe_fit_err = err;
            }
          }

          if (minuit_ok) {
            dW_pred   = W.a   + W.b   * dpe_fit;
            dEm_pred  = Em.a  + Em.b  * dpe_fit;
            dPmz_pred = Pmz.a + Pmz.b * dpe_fit;
            dPmy_pred = Pmy.a + Pmy.b * dpe_fit;

            if (std::isfinite(target.dW_measured)) {
              dW_resid = dW_pred - target.dW_measured;
              dW_resid_err = W.sigma_tot;
            }
            if (std::isfinite(target.dEm_measured)) {
              dEm_resid = dEm_pred - target.dEm_measured;
              dEm_resid_err = Em.sigma_tot;
            }
            if (std::isfinite(target.dPmz_measured)) {
              dPmz_resid = dPmz_pred - target.dPmz_measured;
              dPmz_resid_err = Pmz.sigma_tot;
            }
            if (std::isfinite(target.dPmy_measured)) {
              dPmy_resid = dPmy_pred - target.dPmy_measured;
              dPmy_resid_err = Pmy.sigma_tot;
            }

            chi2_min = 0.0;
            sse_unw  = 0.0;

            for (auto* M : obs) {
              if (!M->used) continue;
              const double pred = M->a + M->b * dpe_fit;
              const double resid = M->y_meas - pred;
              chi2_min += resid * resid / (M->sigma_tot * M->sigma_tot);
              sse_unw  += resid * resid;
            }

            ndf = n_obs_used - 1;
            if (ndf > 0) {
              chi2_ndf = chi2_min / ndf;
              fit_prob = TMath::Prob(chi2_min, ndf);
            }

            fit_valid = 1;
            std::ostringstream note;
            note << "ok";
            if (n_skipped_invalid_sigma > 0) note << "; skipped_invalid_sigma=" << n_skipped_invalid_sigma;
            if (n_skipped_invalid_y > 0)     note << "; skipped_invalid_y=" << n_skipped_invalid_y;
            if (source.source_guess_trivial) note << "; source_guess_trivial=1";
            fit_note = note.str();
          }
        }
        WriteRow(fout, source, target, W, Em, Pmz, Pmy,
                dpe_fit, dpe_fit_err,
                dW_pred, dEm_pred, dPmz_pred, dPmy_pred,
                dW_resid, dW_resid_err,
                dEm_resid, dEm_resid_err,
                dPmz_resid, dPmz_resid_err,
                dPmy_resid, dPmy_resid_err,
                chi2_min, n_obs_used, ndf, chi2_ndf, fit_prob, sse_unw,
                fit_valid, fit_note);

        ++nRowsWritten;
        if (fit_valid) ++nValidFits;
        else ++nInvalidFits;
      }
    }
  }

  fout.close();

  std::cout << "[INFO] Wrote: " << outCsv << "\n"
            << "[INFO] Rows written: " << nRowsWritten << "\n"
            << "[INFO] Valid fits : " << nValidFits << "\n"
            << "[INFO] Invalid fits: " << nInvalidFits << std::endl;
}
