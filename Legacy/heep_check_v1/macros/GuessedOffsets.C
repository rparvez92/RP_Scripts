/*
  GuessedOffsets.C

  Purpose
  -------
  Read per-setting / per-dp-bin weighted mean observable offsets from
    results/tables/OffsetsSummaryBySetting.csv
  and produce a chi2-minimized table of kinematic offsets in
    results/tables/GuessedOffsets.csv

  Philosophy of the guess
  -----------------------
  For each (setting, dp bin), solve for the four kinematic offsets
    (dthe, dpe, dthp, dpp)
  by minimizing the weighted chi2 built from all four observables
  (W, Em, Pmz, Pmy)
  and their measured uncertainties from the input CSV, using a bounded
  Minuit fit.

  The objective function is

    chi2 = Sum_i [ (y_i(measured) - y_i(predicted))^2 / sigma_i^2 ]

  where y_i(predicted) is given by the setting-dependent response matrix.
  Each fit parameter is constrained to lie in [-5, +5].

  The fitted parameter errors are taken from the local curvature at the minimum,
  i.e. from the covariance matrix returned after the Minuit/HESSE evaluation.

  Response matrices used
  ----------------------
  Setting A:
    dW   = -14.08*dthe - 8.62*dpe
    dEm  =              - 7.06*dpe - 2.10*dpp
    dPmz = - 5.75*dthe + 4.10*dpe + 2.27*dpp
    dPmy = + 4.10*dthe + 5.75*dpe - 2.27*dthp

  Setting B:
    dW   = -17.33*dthe - 8.62*dpe
    dEm  =              - 5.66*dpe - 3.63*dpp
    dPmz = - 4.30*dthe + 3.69*dpe + 3.74*dpp
    dPmy = + 3.69*dthe + 4.30*dpe - 3.74*dthp

  Notes
  -----
  * (Setting A, bin b5) is skipped on purpose because there is no faithful
    offset measurement there.
  * Residuals are stored as predicted - measured.
  * The output CSV stores the full fitted covariance matrix for
    (dthe, dpe, dthp, dpp) so downstream steps can propagate source-bin
    uncertainty consistently.
  * The output CSV also stores propagated prediction uncertainties for
    (dW, dEm, dPmz, dPmy) to make later stages simpler and more transparent.
  * Residual errors are computed from the stored prediction uncertainties and
    the measured observable uncertainties, assuming prediction and measurement
    are independent.

  Run command:
    root -l -b -q 'macros/GuessedOffsets.C("results/tables/OffsetsSummaryBySetting.csv","results/tables/GuessedOffsets.csv")'
*/

#include <TString.h>
#include <TMinuit.h>

#include <algorithm>
#include <cmath>
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

struct ObsRow {
  double value = std::numeric_limits<double>::quiet_NaN();
  double err   = std::numeric_limits<double>::quiet_NaN();
  bool   ok    = false;
};

struct BinData {
  std::string setting;
  int         dp_idx = -1;
  std::string dp_label;
  double      dp_lo = std::numeric_limits<double>::quiet_NaN();
  double      dp_hi = std::numeric_limits<double>::quiet_NaN();

  ObsRow W, Em, Pmz, Pmy;
  std::string summary_status = "";
};

struct GuessOut {
  double dthe = std::numeric_limits<double>::quiet_NaN();
  double dpe  = std::numeric_limits<double>::quiet_NaN();
  double dthp = std::numeric_limits<double>::quiet_NaN();
  double dpp  = std::numeric_limits<double>::quiet_NaN();

  double dthe_err = std::numeric_limits<double>::quiet_NaN();
  double dpe_err  = std::numeric_limits<double>::quiet_NaN();
  double dthp_err = std::numeric_limits<double>::quiet_NaN();
  double dpp_err  = std::numeric_limits<double>::quiet_NaN();

  double cov_dthe_dthe = std::numeric_limits<double>::quiet_NaN();
  double cov_dthe_dpe  = std::numeric_limits<double>::quiet_NaN();
  double cov_dthe_dthp = std::numeric_limits<double>::quiet_NaN();
  double cov_dthe_dpp  = std::numeric_limits<double>::quiet_NaN();
  double cov_dpe_dpe   = std::numeric_limits<double>::quiet_NaN();
  double cov_dpe_dthp  = std::numeric_limits<double>::quiet_NaN();
  double cov_dpe_dpp   = std::numeric_limits<double>::quiet_NaN();
  double cov_dthp_dthp = std::numeric_limits<double>::quiet_NaN();
  double cov_dthp_dpp  = std::numeric_limits<double>::quiet_NaN();
  double cov_dpp_dpp   = std::numeric_limits<double>::quiet_NaN();

  double predW   = std::numeric_limits<double>::quiet_NaN();
  double predW_err = std::numeric_limits<double>::quiet_NaN();
  double predEm  = std::numeric_limits<double>::quiet_NaN();
  double predEm_err = std::numeric_limits<double>::quiet_NaN();
  double predPmz = std::numeric_limits<double>::quiet_NaN();
  double predPmz_err = std::numeric_limits<double>::quiet_NaN();
  double predPmy = std::numeric_limits<double>::quiet_NaN();
  double predPmy_err = std::numeric_limits<double>::quiet_NaN();

  double residW   = std::numeric_limits<double>::quiet_NaN();
  double residEm  = std::numeric_limits<double>::quiet_NaN();
  double residPmz = std::numeric_limits<double>::quiet_NaN();
  double residPmy = std::numeric_limits<double>::quiet_NaN();

  double residW_err   = std::numeric_limits<double>::quiet_NaN();
  double residEm_err  = std::numeric_limits<double>::quiet_NaN();
  double residPmz_err = std::numeric_limits<double>::quiet_NaN();
  double residPmy_err = std::numeric_limits<double>::quiet_NaN();

  std::string guess_strategy = "bounded_minuit_chi2_all_observables";
  int guess_valid = 0;
  int source_guess_trivial = 0;
  std::string guess_note;
};

static std::string Trim(const std::string& s) {
  const char* ws = " \t\n\r";
  const auto b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
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
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  for (auto& s : out) s = Trim(s);
  return out;
}

static bool ToInt(const std::string& s, int& x) {
  try {
    size_t pos = 0;
    x = std::stoi(s, &pos);
    return pos == s.size();
  } catch (...) { return false; }
}

static bool ToDouble(const std::string& s, double& x) {
  try {
    size_t pos = 0;
    x = std::stod(s, &pos);
    return pos == s.size();
  } catch (...) { return false; }
}

static bool NearlyZero(double x, double eps = 1e-12) {
  return std::fabs(x) < eps;
}

static bool IsFinite(double x) {
  return std::isfinite(x);
}

static bool IsExcludedBin(const BinData& b) {
  return b.setting == "A" && b.dp_label == "b5";
}

static void Predict(const std::string& setting,
                    double dthe, double dpe, double dthp, double dpp,
                    double& W, double& Em, double& Pmz, double& Pmy) {
  if (setting == "A") {
    W   = -14.08 * dthe - 8.62 * dpe;
    Em  = -7.06 * dpe  - 2.10 * dpp;
    Pmz = - 5.75 * dthe + 4.10 * dpe + 2.27 * dpp;
    Pmy = + 4.10 * dthe + 5.75 * dpe - 2.27 * dthp;
  } else {
    W   = -17.33 * dthe - 8.62 * dpe;
    Em  = -5.66 * dpe  - 3.63 * dpp;
    Pmz = - 4.30 * dthe + 3.69 * dpe + 3.74 * dpp;
    Pmy = + 3.69 * dthe + 4.30 * dpe - 3.74 * dthp;
  }
}

static bool BuildResponseMatrix(const std::string& setting, double A[4][4]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      A[i][j] = 0.0;
    }
  }

  if (setting == "A") {
    A[0][0] = -14.08; A[0][1] = -8.62;
    A[1][1] =  -7.06; A[1][3] = -2.10;
    A[2][0] =  -5.75; A[2][1] =  4.10; A[2][3] =  2.27;
    A[3][0] =   4.10; A[3][1] =  5.75; A[3][2] = -2.27;
    return true;
  }

  if (setting == "B") {
    A[0][0] = -17.33; A[0][1] = -8.62;
    A[1][1] =  -5.66; A[1][3] = -3.63;
    A[2][0] =  -4.30; A[2][1] =  3.69; A[2][3] =  3.74;
    A[3][0] =   3.69; A[3][1] =  4.30; A[3][2] = -3.74;
    return true;
  }

  return false;
}

struct FitContext {
  double A[4][4];
  double y[4];
  double sigma[4];
  bool valid = false;
};

static FitContext gFitContext;

static void Chi2Fcn(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval,
                    Double_t* par, Int_t /*iflag*/) {
  if (!gFitContext.valid) {
    fval = 1e30;
    return;
  }

  double chi2 = 0.0;
  for (int i = 0; i < 4; ++i) {
    double pred = 0.0;
    for (int j = 0; j < 4; ++j) pred += gFitContext.A[i][j] * par[j];
    const double resid = pred - gFitContext.y[i];
    chi2 += resid * resid / (gFitContext.sigma[i] * gFitContext.sigma[i]);
  }

  fval = chi2;
}

static bool SolveChi2Guess(const BinData& b, GuessOut& g) {
  if (!(b.W.ok && b.Em.ok && b.Pmz.ok && b.Pmy.ok)) {
    g.guess_note = "missing_one_or_more_observables";
    return false;
  }

  if (IsExcludedBin(b)) {
    g.guess_note = "excluded_settingA_b5_no_faithful_measurement";
    return false;
  }

  const ObsRow obs[4] = {b.W, b.Em, b.Pmz, b.Pmy};
  for (const auto& o : obs) {
    if (!IsFinite(o.err) || o.err <= 0.0) {
      g.guess_note = "nonpositive_or_invalid_observable_error";
      return false;
    }
  }

  double A[4][4];
  if (!BuildResponseMatrix(b.setting, A)) {
    g.guess_note = "unknown_setting";
    return false;
  }

  for (int i = 0; i < 4; ++i) {
    gFitContext.y[i] = (i == 0 ? b.W.value : i == 1 ? b.Em.value : i == 2 ? b.Pmz.value : b.Pmy.value);
    gFitContext.sigma[i] = (i == 0 ? b.W.err : i == 1 ? b.Em.err : i == 2 ? b.Pmz.err : b.Pmy.err);
    for (int j = 0; j < 4; ++j) {
      gFitContext.A[i][j] = A[i][j];
    }
  }
  gFitContext.valid = true;

  TMinuit minuit(4);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(Chi2Fcn);

  Double_t arglist[2];
  Int_t ierflg = 0;
  arglist[0] = 1.0;
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  const char* names[4] = {"dthe", "dpe", "dthp", "dpp"};
  const double starts[4] = {0.0, 0.0, 0.0, 0.0};
  const double steps[4]  = {0.01, 0.01, 0.01, 0.01};
  const double lo[4]     = {-5.0, -5.0, -5.0, -5.0};
  const double hi[4]     = { 5.0,  5.0,  5.0,  5.0};

  for (int i = 0; i < 4; ++i) {
    minuit.mnparm(i, names[i], starts[i], steps[i], lo[i], hi[i], ierflg);
    if (ierflg != 0) {
      gFitContext.valid = false;
      g.guess_note = "minuit_parameter_setup_failed";
      return false;
    }
  }

  arglist[0] = 2000.0;
  arglist[1] = 1.0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  if (ierflg != 0) {
    gFitContext.valid = false;
    g.guess_note = "migrad_failed";
    return false;
  }

  arglist[0] = 500.0;
  minuit.mnexcm("HESSE", arglist, 1, ierflg);
  if (ierflg != 0) {
    gFitContext.valid = false;
    g.guess_note = "hesse_failed";
    return false;
  }

  Double_t val = 0.0, err = 0.0;
  minuit.GetParameter(0, val, err);
  g.dthe = val; g.dthe_err = err;
  minuit.GetParameter(1, val, err);
  g.dpe = val; g.dpe_err = err;
  minuit.GetParameter(2, val, err);
  g.dthp = val; g.dthp_err = err;
  minuit.GetParameter(3, val, err);
  g.dpp = val; g.dpp_err = err;

  Double_t amin = 0.0, edm = 0.0, errdef = 0.0;
  Int_t nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  if (icstat < 2) {
    gFitContext.valid = false;
    g.guess_note = "minuit_covariance_not_reliable";
    return false;
  }

  Double_t covPacked[16];
  minuit.mnemat(covPacked, 4);
  double cov[4][4];
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      cov[i][j] = covPacked[4 * i + j];
    }
  }

  g.cov_dthe_dthe = cov[0][0];
  g.cov_dthe_dpe  = cov[0][1];
  g.cov_dthe_dthp = cov[0][2];
  g.cov_dthe_dpp  = cov[0][3];
  g.cov_dpe_dpe   = cov[1][1];
  g.cov_dpe_dthp  = cov[1][2];
  g.cov_dpe_dpp   = cov[1][3];
  g.cov_dthp_dthp = cov[2][2];
  g.cov_dthp_dpp  = cov[2][3];
  g.cov_dpp_dpp   = cov[3][3];

  Predict(b.setting, g.dthe, g.dpe, g.dthp, g.dpp,
          g.predW, g.predEm, g.predPmz, g.predPmy);
  g.residW   = g.predW   - b.W.value;
  g.residEm  = g.predEm  - b.Em.value;
  g.residPmz = g.predPmz - b.Pmz.value;
  g.residPmy = g.predPmy - b.Pmy.value;

  const double jacW[4]   = {A[0][0], A[0][1], A[0][2], A[0][3]};
  const double jacEm[4]  = {A[1][0], A[1][1], A[1][2], A[1][3]};
  const double jacPmz[4] = {A[2][0], A[2][1], A[2][2], A[2][3]};
  const double jacPmy[4] = {A[3][0], A[3][1], A[3][2], A[3][3]};

  auto PredVariance = [&](const double jac[4]) {
    double var = 0.0;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        var += jac[i] * cov[i][j] * jac[j];
      }
    }
    return var;
  };

  const double predWVar   = PredVariance(jacW);
  const double predEmVar  = PredVariance(jacEm);
  const double predPmzVar = PredVariance(jacPmz);
  const double predPmyVar = PredVariance(jacPmy);

  g.predW_err   = (predWVar   >= 0.0) ? std::sqrt(predWVar)   : std::numeric_limits<double>::quiet_NaN();
  g.predEm_err  = (predEmVar  >= 0.0) ? std::sqrt(predEmVar)  : std::numeric_limits<double>::quiet_NaN();
  g.predPmz_err = (predPmzVar >= 0.0) ? std::sqrt(predPmzVar) : std::numeric_limits<double>::quiet_NaN();
  g.predPmy_err = (predPmyVar >= 0.0) ? std::sqrt(predPmyVar) : std::numeric_limits<double>::quiet_NaN();

  const double measWVar   = (IsFinite(b.W.err)   && b.W.err   >= 0.0) ? b.W.err   * b.W.err   : std::numeric_limits<double>::quiet_NaN();
  const double measEmVar  = (IsFinite(b.Em.err)  && b.Em.err  >= 0.0) ? b.Em.err  * b.Em.err  : std::numeric_limits<double>::quiet_NaN();
  const double measPmzVar = (IsFinite(b.Pmz.err) && b.Pmz.err >= 0.0) ? b.Pmz.err * b.Pmz.err : std::numeric_limits<double>::quiet_NaN();
  const double measPmyVar = (IsFinite(b.Pmy.err) && b.Pmy.err >= 0.0) ? b.Pmy.err * b.Pmy.err : std::numeric_limits<double>::quiet_NaN();

  g.residW_err   = (IsFinite(g.predW_err)   && IsFinite(measWVar))   ? std::sqrt(g.predW_err   * g.predW_err   + measWVar)   : std::numeric_limits<double>::quiet_NaN();
  g.residEm_err  = (IsFinite(g.predEm_err)  && IsFinite(measEmVar))  ? std::sqrt(g.predEm_err  * g.predEm_err  + measEmVar)  : std::numeric_limits<double>::quiet_NaN();
  g.residPmz_err = (IsFinite(g.predPmz_err) && IsFinite(measPmzVar)) ? std::sqrt(g.predPmz_err * g.predPmz_err + measPmzVar) : std::numeric_limits<double>::quiet_NaN();
  g.residPmy_err = (IsFinite(g.predPmy_err) && IsFinite(measPmyVar)) ? std::sqrt(g.predPmy_err * g.predPmy_err + measPmyVar) : std::numeric_limits<double>::quiet_NaN();

  g.guess_valid = 1;
  g.guess_note = "bounded_minuit_chi2_applied";
  auto atBound = [](double x) {
    return std::fabs(std::fabs(x) - 5.0) < 1e-3;
  };
  if (atBound(g.dthe) || atBound(g.dpe) || atBound(g.dthp) || atBound(g.dpp)) {
    g.guess_note = "bounded_minuit_chi2_applied_parameter_at_limit";
  }
  gFitContext.valid = false;
  return true;
}

static std::string CsvEscape(const std::string& s) {
  if (s.find_first_of(",\"\n") == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

} // namespace

void GuessedOffsets(const char* inCsv  = "results/tables/OffsetsSummaryBySetting.csv",
                    const char* outCsv = "results/tables/GuessedOffsets.csv") {
  std::ifstream fin(inCsv);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input CSV: " << inCsv << std::endl;
    return;
  }

  std::string headerLine;
  if (!std::getline(fin, headerLine)) {
    std::cerr << "[ERROR] Empty input CSV: " << inCsv << std::endl;
    return;
  }

  const auto headers = SplitCsvLine(headerLine);
  std::map<std::string, size_t> col;
  for (size_t i = 0; i < headers.size(); ++i) col[headers[i]] = i;

  const std::vector<std::string> required = {
    "setting", "dp_idx", "dp_label", "dp_lo", "dp_hi", "var",
    "weighted_mean_dmu_MeV", "weighted_mean_err_dmu_MeV", "summary_status"
  };
  for (const auto& key : required) {
    if (!col.count(key)) {
      std::cerr << "[ERROR] Missing required column: " << key << std::endl;
      return;
    }
  }

  std::map<std::string, BinData> bins;
  auto makeKey = [](const std::string& setting, int dp_idx) {
    std::ostringstream os;
    os << setting << "__" << dp_idx;
    return os.str();
  };

  std::string line;
  int nLine = 1;
  while (std::getline(fin, line)) {
    ++nLine;
    if (Trim(line).empty()) continue;
    const auto tok = SplitCsvLine(line);
    if (tok.size() < headers.size()) continue;

    std::string setting = tok[col["setting"]];
    int dp_idx = -1;
    double dp_lo = 0.0, dp_hi = 0.0, value = 0.0, err = 0.0;
    if (!ToInt(tok[col["dp_idx"]], dp_idx)) {
      std::cerr << "[WARN] Bad dp_idx at line " << nLine << std::endl;
      continue;
    }
    ToDouble(tok[col["dp_lo"]], dp_lo);
    ToDouble(tok[col["dp_hi"]], dp_hi);
    if (!ToDouble(tok[col["weighted_mean_dmu_MeV"]], value)) {
      std::cerr << "[WARN] Bad weighted_mean_dmu_MeV at line " << nLine << std::endl;
      continue;
    }
    if (!ToDouble(tok[col["weighted_mean_err_dmu_MeV"]], err)) {
      std::cerr << "[WARN] Bad weighted_mean_err_dmu_MeV at line " << nLine << std::endl;
      continue;
    }
    const std::string dp_label = tok[col["dp_label"]];
    const std::string var = tok[col["var"]];
    const std::string summary_status = tok[col["summary_status"]];

    BinData& b = bins[makeKey(setting, dp_idx)];
    b.setting = setting;
    b.dp_idx = dp_idx;
    b.dp_label = dp_label;
    b.dp_lo = dp_lo;
    b.dp_hi = dp_hi;
    b.summary_status = summary_status;

    ObsRow row;
    row.value = value;
    row.err = err;
    row.ok = true;

    if      (var == "W")   b.W = row;
    else if (var == "Em")  b.Em = row;
    else if (var == "Pmz") b.Pmz = row;
    else if (var == "Pmy") b.Pmy = row;
  }
  fin.close();

  std::vector<BinData> rows;
  rows.reserve(bins.size());
  for (auto& kv : bins) rows.push_back(kv.second);
  std::sort(rows.begin(), rows.end(), [](const BinData& a, const BinData& b) {
    if (a.setting != b.setting) return a.setting < b.setting;
    return a.dp_idx < b.dp_idx;
  });

  std::ofstream fout(outCsv);
  if (!fout) {
    std::cerr << "[ERROR] Cannot open output CSV: " << outCsv << std::endl;
    return;
  }

  fout << std::fixed << std::setprecision(6);
  fout << "setting,dp_idx,bin,dp_lo,dp_hi,summary_status,";
  fout << "dW_measured,dW_err,dEm_measured,dEm_err,dPmz_measured,dPmz_err,dPmy_measured,dPmy_err,";
  fout << "dthe_guess,dthe_guess_err,dpe_guess,dpe_guess_err,dthp_guess,dthp_guess_err,dpp_guess,dpp_guess_err,";
  fout << "cov_dthe_dthe,cov_dthe_dpe,cov_dthe_dthp,cov_dthe_dpp,"
       << "cov_dpe_dpe,cov_dpe_dthp,cov_dpe_dpp,"
       << "cov_dthp_dthp,cov_dthp_dpp,cov_dpp_dpp,";
  fout << "dW_pred,dW_pred_err,dEm_pred,dEm_pred_err,dPmz_pred,dPmz_pred_err,dPmy_pred,dPmy_pred_err,";
  fout << "dW_resid,dW_resid_err,dEm_resid,dEm_resid_err,dPmz_resid,dPmz_resid_err,dPmy_resid,dPmy_resid_err,";
  fout << "guess_strategy,guess_valid,source_guess_trivial,guess_note\n";

  int nGood = 0;
  int nBad  = 0;
  for (const auto& b : rows) {
    GuessOut g;
    const bool ok = SolveChi2Guess(b, g);
    if (ok && g.guess_valid) {
      ++nGood;
    } else {
      ++nBad;
      // Keep prediction / residual fields as NaN if guess failed.
    }

    fout << CsvEscape(b.setting) << ","
         << b.dp_idx << ","
         << CsvEscape(b.dp_label) << ","
         << b.dp_lo << ","
         << b.dp_hi << ","
         << CsvEscape(b.summary_status) << ","
         << b.W.value << ","   << b.W.err << ","
         << b.Em.value << ","  << b.Em.err << ","
         << b.Pmz.value << "," << b.Pmz.err << ","
         << b.Pmy.value << "," << b.Pmy.err << ","
         << g.dthe << "," << g.dthe_err << ","
         << g.dpe  << "," << g.dpe_err  << ","
         << g.dthp << "," << g.dthp_err << ","
         << g.dpp  << "," << g.dpp_err  << ","
         << g.cov_dthe_dthe << "," << g.cov_dthe_dpe << "," << g.cov_dthe_dthp << "," << g.cov_dthe_dpp << ","
         << g.cov_dpe_dpe << "," << g.cov_dpe_dthp << "," << g.cov_dpe_dpp << ","
         << g.cov_dthp_dthp << "," << g.cov_dthp_dpp << "," << g.cov_dpp_dpp << ","
         << g.predW << "," << g.predW_err << ","
         << g.predEm << "," << g.predEm_err << ","
         << g.predPmz << "," << g.predPmz_err << ","
         << g.predPmy << "," << g.predPmy_err << ","
         << g.residW << "," << g.residW_err << ","
         << g.residEm << "," << g.residEm_err << ","
         << g.residPmz << "," << g.residPmz_err << ","
         << g.residPmy << "," << g.residPmy_err << ","
         << CsvEscape(g.guess_strategy) << ","
         << g.guess_valid << ","
         << g.source_guess_trivial << ","
         << CsvEscape(g.guess_note)
         << "\n";
  }
  fout.close();

  std::cout << "[INFO] Wrote: " << outCsv << std::endl;
  std::cout << "[INFO] Rows written: " << rows.size() << std::endl;
  std::cout << "[INFO] Valid guesses: " << nGood << std::endl;
  std::cout << "[INFO] Invalid guesses: " << nBad  << std::endl;
  std::cout << "[INFO] Run with:" << std::endl;
  std::cout << "  root -l -b -q 'macros/GuessedOffsets.C(\"results/tables/OffsetsSummaryBySetting.csv\",\"results/tables/GuessedOffsets.csv\")'" << std::endl;
}
