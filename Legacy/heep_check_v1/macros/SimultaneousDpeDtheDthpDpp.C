/*
  SimultaneousDpeDtheDthpDpp.C

  Purpose
  -------
  Read setting-wise measured observable offsets from
    results/tables/OffsetsSummaryBySetting.csv
  and perform a simultaneous fit, per setting, under the hypothesis that:

    - dthe is bin-dependent
    - dpe   is bin-dependent
    - dthp  is bin-dependent
    - dpp   is bin-dependent

  Parameter ordering
  ------------------
  For consistency in code, comments, and covariance output, the fit parameter
  ordering is:

    dthe_b1, dthe_b2, dthe_b3, dthe_b4, [dthe_b5 when present],
    dpe_b1,  dpe_b2,  dpe_b3,  dpe_b4,  [dpe_b5 when present],
    dthp_b1, dthp_b2, dthp_b3, dthp_b4, [dthp_b5 when present],
    dpp_b1,  dpp_b2,  dpp_b3,  dpp_b4,  [dpp_b5 when present]

  Bin usage
  ---------
    Setting A: b1, b2, b3, b4
    Setting B: b1, b2, b3, b4, b5

  We do not use:
    - full bins
    - any row with summary_status != "ok"
    - Setting A, b5

  Fit modes
  ---------
    unconstrained: no parameter bounds
    constrained  : every fit parameter is limited to [-5, +5]

  Objective function
  ------------------
    chi2 = sum_i [ (y_i(measured) - y_i(predicted))^2 / sigma_i^2 ]

  where:
    - y_i(measured) is weighted_mean_dmu_MeV
    - sigma_i       is weighted_mean_err_dmu_MeV
    - y_i(predicted) is given by the setting-dependent linear response model

  Output
  ------
  One row per fitted narrow bin. The output also includes the full covariance
  matrix of the simultaneous fit, as well as predicted values, predicted
  errors, residuals, and residual errors. The CSV schema is kept fixed across
  settings by always reserving columns for b1 through b5 for each parameter
  family; fields not used by a given setting are written as nan.

  Default output names:
    results/tables/SimultaneousDpeDtheDthpDpp_unconstrained.csv
    results/tables/SimultaneousDpeDtheDthpDpp_constrained.csv

  Run examples:
    root -l -b -q 'macros/SimultaneousDpeDtheDthpDpp.C("results/tables/OffsetsSummaryBySetting.csv","","unconstrained")'
    root -l -b -q 'macros/SimultaneousDpeDtheDthpDpp.C("results/tables/OffsetsSummaryBySetting.csv","","constrained")'
*/

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
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ObsSummary {
  double value = std::numeric_limits<double>::quiet_NaN();
  double err   = std::numeric_limits<double>::quiet_NaN();
  bool ok      = false;
};

struct BinData {
  std::string setting;
  int dp_idx = -1;
  std::string bin;
  double dp_lo = std::numeric_limits<double>::quiet_NaN();
  double dp_hi = std::numeric_limits<double>::quiet_NaN();
  std::string summary_status;

  ObsSummary W;
  ObsSummary Em;
  ObsSummary Pmz;
  ObsSummary Pmy;
};

struct FitBinResult {
  std::string bin;
  double dthe = std::numeric_limits<double>::quiet_NaN();
  double dthe_err = std::numeric_limits<double>::quiet_NaN();
  double dpe = std::numeric_limits<double>::quiet_NaN();
  double dpe_err = std::numeric_limits<double>::quiet_NaN();
  double dthp = std::numeric_limits<double>::quiet_NaN();
  double dthp_err = std::numeric_limits<double>::quiet_NaN();
  double dpp = std::numeric_limits<double>::quiet_NaN();
  double dpp_err = std::numeric_limits<double>::quiet_NaN();

  double dW_measured = std::numeric_limits<double>::quiet_NaN();
  double dW_err = std::numeric_limits<double>::quiet_NaN();
  double dEm_measured = std::numeric_limits<double>::quiet_NaN();
  double dEm_err = std::numeric_limits<double>::quiet_NaN();
  double dPmz_measured = std::numeric_limits<double>::quiet_NaN();
  double dPmz_err = std::numeric_limits<double>::quiet_NaN();
  double dPmy_measured = std::numeric_limits<double>::quiet_NaN();
  double dPmy_err = std::numeric_limits<double>::quiet_NaN();

  double dW_pred = std::numeric_limits<double>::quiet_NaN();
  double dW_pred_err = std::numeric_limits<double>::quiet_NaN();
  double dEm_pred = std::numeric_limits<double>::quiet_NaN();
  double dEm_pred_err = std::numeric_limits<double>::quiet_NaN();
  double dPmz_pred = std::numeric_limits<double>::quiet_NaN();
  double dPmz_pred_err = std::numeric_limits<double>::quiet_NaN();
  double dPmy_pred = std::numeric_limits<double>::quiet_NaN();
  double dPmy_pred_err = std::numeric_limits<double>::quiet_NaN();

  double dW_resid = std::numeric_limits<double>::quiet_NaN();
  double dW_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid = std::numeric_limits<double>::quiet_NaN();
  double dEm_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid = std::numeric_limits<double>::quiet_NaN();
  double dPmz_resid_err = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid = std::numeric_limits<double>::quiet_NaN();
  double dPmy_resid_err = std::numeric_limits<double>::quiet_NaN();
};

struct SettingFitResult {
  std::string setting;
  std::string mode;

  double chi2_min = std::numeric_limits<double>::quiet_NaN();
  int n_obs_used = 0;
  int n_params = 0;
  int ndf = -1;
  double chi2_ndf = std::numeric_limits<double>::quiet_NaN();
  double fit_prob = std::numeric_limits<double>::quiet_NaN();
  int fit_valid = 0;
  std::string fit_note;

  std::vector<std::string> param_names;
  std::vector<double> params;
  std::vector<double> param_errs;
  std::vector<double> covariance;
  std::vector<FitBinResult> bins;
};

struct ObsEquation {
  int bin_index = -1;
  double coeff_dthe = 0.0;
  double coeff_dpe = 0.0;
  double coeff_dthp = 0.0;
  double coeff_dpp = 0.0;
  double measured = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
};

struct FitContext {
  std::vector<ObsEquation> equations;
  std::vector<int> dthe_param_index_by_bin;
  std::vector<int> dpe_param_index_by_bin;
  std::vector<int> dthp_param_index_by_bin;
  std::vector<int> dpp_param_index_by_bin;
  bool valid = false;
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

static bool ToInt(const std::string& s, int& x) {
  try {
    size_t pos = 0;
    x = std::stoi(Trim(s), &pos);
    return pos == Trim(s).size();
  } catch (...) {
    return false;
  }
}

static bool ToDouble(const std::string& s, double& x) {
  try {
    size_t pos = 0;
    const std::string t = Trim(s);
    x = std::stod(t, &pos);
    return pos == t.size();
  } catch (...) {
    return false;
  }
}

static std::string NormalizeSetting(const std::string& s) {
  const std::string t = Lower(Trim(s));
  if (t == "a" || t == "setting a" || t == "settinga") return "A";
  if (t == "b" || t == "setting b" || t == "settingb") return "B";
  return Trim(s);
}

static bool IsAllowedBinForSetting(const std::string& setting, const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (setting == "A") return b == "b1" || b == "b2" || b == "b3" || b == "b4";
  if (setting == "B") return b == "b1" || b == "b2" || b == "b3" || b == "b4" || b == "b5";
  return false;
}

static std::vector<std::string> BinOrderForSetting(const std::string& setting) {
  if (setting == "A") return {"b1", "b2", "b3", "b4"};
  if (setting == "B") return {"b1", "b2", "b3", "b4", "b5"};
  return {};
}

static std::vector<std::string> GlobalParamNames() {
  return {"dthe_b1", "dthe_b2", "dthe_b3", "dthe_b4", "dthe_b5",
          "dpe_b1", "dpe_b2", "dpe_b3", "dpe_b4", "dpe_b5",
          "dthp_b1", "dthp_b2", "dthp_b3", "dthp_b4", "dthp_b5",
          "dpp_b1", "dpp_b2", "dpp_b3", "dpp_b4", "dpp_b5"};
}

static std::string CsvQuoteIfNeeded(const std::string& s) {
  if (s.find_first_of(",\"\n") == std::string::npos) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

static void SimultaneousChi2Fcn(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval,
                                Double_t* par, Int_t /*iflag*/) {
  if (!gFitContext.valid) {
    fval = 1e30;
    return;
  }

  double chi2 = 0.0;
  for (const auto& eq : gFitContext.equations) {
    const int idx_dthe = gFitContext.dthe_param_index_by_bin[eq.bin_index];
    const int idx_dpe = gFitContext.dpe_param_index_by_bin[eq.bin_index];
    const int idx_dthp = gFitContext.dthp_param_index_by_bin[eq.bin_index];
    const int idx_dpp = gFitContext.dpp_param_index_by_bin[eq.bin_index];
    const double pred =
        eq.coeff_dthe * par[idx_dthe] +
        eq.coeff_dpe  * par[idx_dpe] +
        eq.coeff_dthp * par[idx_dthp] +
        eq.coeff_dpp  * par[idx_dpp];
    const double resid = eq.measured - pred;
    chi2 += resid * resid / (eq.sigma * eq.sigma);
  }

  fval = chi2;
}

static std::map<std::string, std::map<std::string, BinData>>
ReadOffsetsSummary(const std::string& csvPath) {
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
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[Lower(headers[i])] = i;

  auto idx = [&](const std::string& name) -> int {
    auto it = col.find(Lower(name));
    return (it == col.end()) ? -1 : it->second;
  };

  const int iSetting = idx("setting");
  const int iDpIdx = idx("dp_idx");
  const int iBin = idx("dp_label");
  const int iDpLo = idx("dp_lo");
  const int iDpHi = idx("dp_hi");
  const int iVar = idx("var");
  const int iValue = idx("weighted_mean_dmu_mev");
  const int iErr = idx("weighted_mean_err_dmu_mev");
  const int iStatus = idx("summary_status");

  if (iSetting < 0 || iDpIdx < 0 || iBin < 0 || iDpLo < 0 || iDpHi < 0 ||
      iVar < 0 || iValue < 0 || iErr < 0 || iStatus < 0) {
    std::cerr << "[ERROR] Missing one or more required columns in " << csvPath << std::endl;
    return {};
  }

  std::map<std::string, std::map<std::string, BinData>> out;
  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);

    auto getS = [&](int i) -> std::string {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    auto getD = [&](int i) -> double {
      double v = std::numeric_limits<double>::quiet_NaN();
      if (i >= 0 && i < static_cast<int>(fields.size())) ToDouble(fields[i], v);
      return v;
    };

    const std::string setting = NormalizeSetting(getS(iSetting));
    const std::string bin = getS(iBin);
    const std::string status = Lower(getS(iStatus));
    const std::string var = getS(iVar);

    if (Lower(bin) == "full") continue;
    if (!IsAllowedBinForSetting(setting, bin)) continue;
    if (status != "ok") continue;

    BinData& b = out[setting][bin];
    b.setting = setting;
    b.bin = bin;
    b.summary_status = getS(iStatus);
    ToInt(getS(iDpIdx), b.dp_idx);
    b.dp_lo = getD(iDpLo);
    b.dp_hi = getD(iDpHi);

    ObsSummary obs;
    obs.value = getD(iValue);
    obs.err = getD(iErr);
    obs.ok = std::isfinite(obs.value) && std::isfinite(obs.err) && obs.err > 0.0;

    if (var == "W") b.W = obs;
    else if (var == "Em") b.Em = obs;
    else if (var == "Pmz") b.Pmz = obs;
    else if (var == "Pmy") b.Pmy = obs;
  }

  return out;
}

static bool BuildSettingEquations(const std::string& setting,
                                  const std::map<std::string, BinData>& byBin,
                                  std::vector<std::string>& bins,
                                  std::vector<ObsEquation>& eqs,
                                  std::string& note) {
  bins.clear();
  eqs.clear();

  const auto desired = BinOrderForSetting(setting);
  if (desired.empty()) {
    note = "unsupported_setting";
    return false;
  }

  for (const auto& bin : desired) {
    auto it = byBin.find(bin);
    if (it == byBin.end()) {
      note = "missing_required_bin_" + bin;
      return false;
    }
    const BinData& b = it->second;
    if (!(b.W.ok && b.Em.ok && b.Pmz.ok && b.Pmy.ok)) {
      note = "missing_one_or_more_valid_observables_in_" + bin;
      return false;
    }
    bins.push_back(bin);
  }

  for (int ibin = 0; ibin < static_cast<int>(bins.size()); ++ibin) {
    const BinData& b = byBin.at(bins[ibin]);

    if (setting == "A") {
      eqs.push_back({ibin, -14.08, -8.62,  0.0,   0.0,  b.W.value,   b.W.err});
      eqs.push_back({ibin,   0.0,  -7.06,  0.0,  -2.10, b.Em.value,  b.Em.err});
      eqs.push_back({ibin, -5.75,   4.10,  0.0,   2.27, b.Pmz.value, b.Pmz.err});
      eqs.push_back({ibin,  4.10,   5.75, -2.27,  0.0,  b.Pmy.value, b.Pmy.err});
    } else if (setting == "B") {
      eqs.push_back({ibin, -17.33, -8.62,  0.0,   0.0,  b.W.value,   b.W.err});
      eqs.push_back({ibin,   0.0,  -5.66,  0.0,  -3.63, b.Em.value,  b.Em.err});
      eqs.push_back({ibin, -4.30,   3.69,  0.0,   3.74, b.Pmz.value, b.Pmz.err});
      eqs.push_back({ibin,  3.69,   4.30, -3.74,  0.0,  b.Pmy.value, b.Pmy.err});
    } else {
      note = "unsupported_setting";
      return false;
    }
  }

  return true;
}

static void BuildParamNames(const std::vector<std::string>& bins,
                            std::vector<std::string>& names) {
  names.clear();
  for (const auto& bin : bins) names.push_back("dthe_" + bin);
  for (const auto& bin : bins) names.push_back("dpe_" + bin);
  for (const auto& bin : bins) names.push_back("dthp_" + bin);
  for (const auto& bin : bins) names.push_back("dpp_" + bin);
}

static double PredVarianceForBin(const std::vector<double>& cov,
                                 int npar,
                                 int idx_dthe,
                                 int idx_dpe,
                                 int idx_dthp,
                                 int idx_dpp,
                                 double c_dthe,
                                 double c_dpe,
                                 double c_dthp,
                                 double c_dpp) {
  const double coeffs[4] = {c_dthe, c_dpe, c_dthp, c_dpp};
  const int idxs[4] = {idx_dthe, idx_dpe, idx_dthp, idx_dpp};

  double var = 0.0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      var += coeffs[i] * cov[idxs[i] * npar + idxs[j]] * coeffs[j];
    }
  }
  if (!std::isfinite(var)) return std::numeric_limits<double>::quiet_NaN();
  if (var < 0.0 && std::fabs(var) < 1e-12) return 0.0;
  return (var >= 0.0) ? var : std::numeric_limits<double>::quiet_NaN();
}

static bool RunSettingFit(const std::string& setting,
                          const std::map<std::string, BinData>& byBin,
                          bool constrained,
                          SettingFitResult& res) {
  res = SettingFitResult{};
  res.setting = setting;
  res.mode = constrained ? "constrained" : "unconstrained";

  std::vector<std::string> bins;
  std::vector<ObsEquation> eqs;
  std::string note;
  if (!BuildSettingEquations(setting, byBin, bins, eqs, note)) {
    res.fit_valid = 0;
    res.fit_note = note;
    return false;
  }

  BuildParamNames(bins, res.param_names);
  const int nbin = static_cast<int>(bins.size());
  const int npar = static_cast<int>(res.param_names.size());
  res.n_params = npar;

  gFitContext.equations = eqs;
  gFitContext.dthe_param_index_by_bin.resize(nbin);
  gFitContext.dpe_param_index_by_bin.resize(nbin);
  gFitContext.dthp_param_index_by_bin.resize(nbin);
  gFitContext.dpp_param_index_by_bin.resize(nbin);
  for (int i = 0; i < nbin; ++i) {
    gFitContext.dthe_param_index_by_bin[i] = i;
    gFitContext.dpe_param_index_by_bin[i] = nbin + i;
    gFitContext.dthp_param_index_by_bin[i] = 2 * nbin + i;
    gFitContext.dpp_param_index_by_bin[i] = 3 * nbin + i;
  }
  gFitContext.valid = true;

  TMinuit minuit(npar);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(SimultaneousChi2Fcn);

  Double_t arglist[2];
  Int_t ierflg = 0;
  arglist[0] = 1.0;
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  for (int i = 0; i < npar; ++i) {
    const char* pname = res.param_names[i].c_str();
    const double start = 0.0;
    const double step = 0.01;
    const double low = constrained ? -5.0 : 0.0;
    const double high = constrained ? 5.0 : 0.0;
    minuit.mnparm(i, pname, start, step, low, high, ierflg);
    if (ierflg != 0) {
      gFitContext.valid = false;
      res.fit_valid = 0;
      res.fit_note = "minuit_parameter_setup_failed";
      return false;
    }
  }

  arglist[0] = 5000.0;
  arglist[1] = 1.0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  if (ierflg != 0) {
    gFitContext.valid = false;
    res.fit_valid = 0;
    res.fit_note = "migrad_failed";
    return false;
  }

  arglist[0] = 1000.0;
  minuit.mnexcm("HESSE", arglist, 1, ierflg);
  if (ierflg != 0) {
    gFitContext.valid = false;
    res.fit_valid = 0;
    res.fit_note = "hesse_failed";
    return false;
  }

  Double_t amin = 0.0, edm = 0.0, errdef = 0.0;
  Int_t nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  if (icstat < 2) {
    gFitContext.valid = false;
    res.fit_valid = 0;
    res.fit_note = "minuit_covariance_not_reliable";
    return false;
  }

  res.params.resize(npar, std::numeric_limits<double>::quiet_NaN());
  res.param_errs.resize(npar, std::numeric_limits<double>::quiet_NaN());
  for (int i = 0; i < npar; ++i) {
    Double_t val = 0.0, err = 0.0;
    minuit.GetParameter(i, val, err);
    res.params[i] = val;
    res.param_errs[i] = err;
  }

  res.covariance.resize(npar * npar, std::numeric_limits<double>::quiet_NaN());
  std::vector<Double_t> covPacked(npar * npar, 0.0);
  minuit.mnemat(covPacked.data(), npar);
  for (int i = 0; i < npar * npar; ++i) res.covariance[i] = covPacked[i];

  res.chi2_min = amin;
  res.n_obs_used = static_cast<int>(eqs.size());
  res.ndf = res.n_obs_used - res.n_params;
  if (res.ndf > 0) {
    res.chi2_ndf = res.chi2_min / res.ndf;
    res.fit_prob = TMath::Prob(res.chi2_min, res.ndf);
  }
  res.fit_valid = 1;
  res.fit_note = constrained ? "ok_constrained" : "ok_unconstrained";

  auto atBound = [](double x) { return std::fabs(std::fabs(x) - 5.0) < 1e-3; };
  if (constrained) {
    bool onLimit = false;
    for (int i = 0; i < npar; ++i) onLimit = onLimit || atBound(res.params[i]);
    if (onLimit) res.fit_note += "_parameter_at_limit";
  }

  for (int ibin = 0; ibin < nbin; ++ibin) {
    const BinData& b = byBin.at(bins[ibin]);
    const int idx_dthe = ibin;
    const int idx_dpe = nbin + ibin;
    const int idx_dthp = 2 * nbin + ibin;
    const int idx_dpp = 3 * nbin + ibin;

    FitBinResult bout;
    bout.bin = bins[ibin];
    bout.dthe = res.params[idx_dthe];
    bout.dthe_err = res.param_errs[idx_dthe];
    bout.dpe = res.params[idx_dpe];
    bout.dpe_err = res.param_errs[idx_dpe];
    bout.dthp = res.params[idx_dthp];
    bout.dthp_err = res.param_errs[idx_dthp];
    bout.dpp = res.params[idx_dpp];
    bout.dpp_err = res.param_errs[idx_dpp];

    bout.dW_measured = b.W.value;      bout.dW_err = b.W.err;
    bout.dEm_measured = b.Em.value;    bout.dEm_err = b.Em.err;
    bout.dPmz_measured = b.Pmz.value;  bout.dPmz_err = b.Pmz.err;
    bout.dPmy_measured = b.Pmy.value;  bout.dPmy_err = b.Pmy.err;

    if (setting == "A") {
      bout.dW_pred   = -14.08 * bout.dthe + -8.62 * bout.dpe;
      bout.dEm_pred  = -7.06  * bout.dpe + -2.10 * bout.dpp;
      bout.dPmz_pred = -5.75  * bout.dthe +  4.10 * bout.dpe + 2.27 * bout.dpp;
      bout.dPmy_pred =  4.10  * bout.dthe +  5.75 * bout.dpe - 2.27 * bout.dthp;

      const double vW   = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -14.08, -8.62,  0.0,   0.0);
      const double vEm  = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp,   0.0,  -7.06,  0.0,  -2.10);
      const double vPmz = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -5.75,   4.10,  0.0,   2.27);
      const double vPmy = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp,  4.10,   5.75, -2.27,  0.0);
      bout.dW_pred_err   = std::isfinite(vW)   ? std::sqrt(vW)   : std::numeric_limits<double>::quiet_NaN();
      bout.dEm_pred_err  = std::isfinite(vEm)  ? std::sqrt(vEm)  : std::numeric_limits<double>::quiet_NaN();
      bout.dPmz_pred_err = std::isfinite(vPmz) ? std::sqrt(vPmz) : std::numeric_limits<double>::quiet_NaN();
      bout.dPmy_pred_err = std::isfinite(vPmy) ? std::sqrt(vPmy) : std::numeric_limits<double>::quiet_NaN();
    } else {
      bout.dW_pred   = -17.33 * bout.dthe + -8.62 * bout.dpe;
      bout.dEm_pred  = -5.66  * bout.dpe + -3.63 * bout.dpp;
      bout.dPmz_pred = -4.30  * bout.dthe +  3.69 * bout.dpe + 3.74 * bout.dpp;
      bout.dPmy_pred =  3.69  * bout.dthe +  4.30 * bout.dpe - 3.74 * bout.dthp;

      const double vW   = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -17.33, -8.62,  0.0,   0.0);
      const double vEm  = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp,   0.0,  -5.66,  0.0,  -3.63);
      const double vPmz = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, -4.30,   3.69,  0.0,   3.74);
      const double vPmy = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp,  3.69,   4.30, -3.74,  0.0);
      bout.dW_pred_err   = std::isfinite(vW)   ? std::sqrt(vW)   : std::numeric_limits<double>::quiet_NaN();
      bout.dEm_pred_err  = std::isfinite(vEm)  ? std::sqrt(vEm)  : std::numeric_limits<double>::quiet_NaN();
      bout.dPmz_pred_err = std::isfinite(vPmz) ? std::sqrt(vPmz) : std::numeric_limits<double>::quiet_NaN();
      bout.dPmy_pred_err = std::isfinite(vPmy) ? std::sqrt(vPmy) : std::numeric_limits<double>::quiet_NaN();
    }

    bout.dW_resid = bout.dW_pred - bout.dW_measured;
    bout.dEm_resid = bout.dEm_pred - bout.dEm_measured;
    bout.dPmz_resid = bout.dPmz_pred - bout.dPmz_measured;
    bout.dPmy_resid = bout.dPmy_pred - bout.dPmy_measured;

    bout.dW_resid_err   = (std::isfinite(bout.dW_pred_err)   && std::isfinite(bout.dW_err))   ? std::sqrt(bout.dW_pred_err   * bout.dW_pred_err   + bout.dW_err   * bout.dW_err)   : std::numeric_limits<double>::quiet_NaN();
    bout.dEm_resid_err  = (std::isfinite(bout.dEm_pred_err)  && std::isfinite(bout.dEm_err))  ? std::sqrt(bout.dEm_pred_err  * bout.dEm_pred_err  + bout.dEm_err  * bout.dEm_err)  : std::numeric_limits<double>::quiet_NaN();
    bout.dPmz_resid_err = (std::isfinite(bout.dPmz_pred_err) && std::isfinite(bout.dPmz_err)) ? std::sqrt(bout.dPmz_pred_err * bout.dPmz_pred_err + bout.dPmz_err * bout.dPmz_err) : std::numeric_limits<double>::quiet_NaN();
    bout.dPmy_resid_err = (std::isfinite(bout.dPmy_pred_err) && std::isfinite(bout.dPmy_err)) ? std::sqrt(bout.dPmy_pred_err * bout.dPmy_pred_err + bout.dPmy_err * bout.dPmy_err) : std::numeric_limits<double>::quiet_NaN();

    res.bins.push_back(bout);
  }

  gFitContext.valid = false;
  return true;
}

static void WriteHeader(std::ofstream& fout, const std::vector<std::string>& param_names) {
  fout << "setting,mode,bin,dp_idx,dp_lo,dp_hi,";
  for (const auto& name : param_names) {
    fout << name << "_fit," << name << "_fit_err,";
  }
  fout << "dW_measured,dW_err,dEm_measured,dEm_err,dPmz_measured,dPmz_err,dPmy_measured,dPmy_err,";
  fout << "dW_pred,dW_pred_err,dEm_pred,dEm_pred_err,dPmz_pred,dPmz_pred_err,dPmy_pred,dPmy_pred_err,";
  fout << "dW_resid,dW_resid_err,dEm_resid,dEm_resid_err,dPmz_resid,dPmz_resid_err,dPmy_resid,dPmy_resid_err,";
  fout << "chi2_min,n_obs_used,n_params,ndf,chi2_ndf,fit_prob,fit_valid,fit_note,";
  for (const auto& pi : param_names) {
    for (const auto& pj : param_names) {
      fout << "cov_" << pi << "_" << pj << ",";
    }
  }
  fout << "end_marker\n";
}

static int FindParamIndex(const std::vector<std::string>& names, const std::string& target) {
  for (int i = 0; i < static_cast<int>(names.size()); ++i) {
    if (names[i] == target) return i;
  }
  return -1;
}

static double ParamValueOrNaN(const SettingFitResult& res, const std::string& name, bool wantError) {
  const int idx = FindParamIndex(res.param_names, name);
  if (idx < 0) return std::numeric_limits<double>::quiet_NaN();
  if (wantError) {
    return (idx < static_cast<int>(res.param_errs.size())) ? res.param_errs[idx]
                                                           : std::numeric_limits<double>::quiet_NaN();
  }
  return (idx < static_cast<int>(res.params.size())) ? res.params[idx]
                                                     : std::numeric_limits<double>::quiet_NaN();
}

static double CovarianceOrNaN(const SettingFitResult& res,
                              const std::string& pi,
                              const std::string& pj) {
  const int i = FindParamIndex(res.param_names, pi);
  const int j = FindParamIndex(res.param_names, pj);
  if (i < 0 || j < 0) return std::numeric_limits<double>::quiet_NaN();
  if (res.n_params <= 0) return std::numeric_limits<double>::quiet_NaN();
  const int flat = i * res.n_params + j;
  return (flat >= 0 && flat < static_cast<int>(res.covariance.size()))
             ? res.covariance[flat]
             : std::numeric_limits<double>::quiet_NaN();
}

static void WriteSettingRows(std::ofstream& fout,
                             const SettingFitResult& res,
                             const std::map<std::string, BinData>& byBin,
                             const std::vector<std::string>& global_param_names) {
  for (const auto& bres : res.bins) {
    const BinData& b = byBin.at(bres.bin);
    fout << CsvQuoteIfNeeded(res.setting) << ","
         << CsvQuoteIfNeeded(res.mode) << ","
         << CsvQuoteIfNeeded(bres.bin) << ","
         << b.dp_idx << ","
         << b.dp_lo << ","
         << b.dp_hi << ",";

    for (const auto& name : global_param_names) {
      fout << ParamValueOrNaN(res, name, false) << ","
           << ParamValueOrNaN(res, name, true) << ",";
    }

    fout << bres.dW_measured << "," << bres.dW_err << ","
         << bres.dEm_measured << "," << bres.dEm_err << ","
         << bres.dPmz_measured << "," << bres.dPmz_err << ","
         << bres.dPmy_measured << "," << bres.dPmy_err << ","
         << bres.dW_pred << "," << bres.dW_pred_err << ","
         << bres.dEm_pred << "," << bres.dEm_pred_err << ","
         << bres.dPmz_pred << "," << bres.dPmz_pred_err << ","
         << bres.dPmy_pred << "," << bres.dPmy_pred_err << ","
         << bres.dW_resid << "," << bres.dW_resid_err << ","
         << bres.dEm_resid << "," << bres.dEm_resid_err << ","
         << bres.dPmz_resid << "," << bres.dPmz_resid_err << ","
         << bres.dPmy_resid << "," << bres.dPmy_resid_err << ","
         << res.chi2_min << "," << res.n_obs_used << "," << res.n_params << ","
         << res.ndf << "," << res.chi2_ndf << "," << res.fit_prob << ","
         << res.fit_valid << "," << CsvQuoteIfNeeded(res.fit_note) << ",";

    for (const auto& pi : global_param_names) {
      for (const auto& pj : global_param_names) {
        fout << CovarianceOrNaN(res, pi, pj) << ",";
      }
    }
    fout << "done\n";
  }
}

static std::string ResolveMode(const std::string& modeRaw) {
  const std::string mode = Lower(Trim(modeRaw));
  if (mode == "constrained") return "constrained";
  if (mode == "unconstrained" || mode.empty()) return "unconstrained";
  return "";
}

} // namespace

void SimultaneousDpeDtheDthpDpp(const char* inCsv = "results/tables/OffsetsSummaryBySetting.csv",
                                const char* outCsv = "",
                                const char* modeIn = "unconstrained") {
  const std::string mode = ResolveMode(modeIn ? modeIn : "");
  if (mode.empty()) {
    std::cerr << "[ERROR] Unknown mode. Use \"constrained\" or \"unconstrained\"." << std::endl;
    return;
  }
  const bool constrained = (mode == "constrained");

  const auto bySettingBin = ReadOffsetsSummary(inCsv);
  if (bySettingBin.empty()) {
    std::cerr << "[ERROR] No valid rows read from " << inCsv << std::endl;
    return;
  }

  std::string outPath = Trim(outCsv ? outCsv : "");
  if (outPath.empty()) {
    outPath = "results/tables/SimultaneousDpeDtheDthpDpp_" + mode + ".csv";
  }

  std::ofstream fout(outPath);
  if (!fout) {
    std::cerr << "[ERROR] Cannot open output CSV for writing: " << outPath << std::endl;
    return;
  }
  fout << std::fixed << std::setprecision(6);

  int rowsWritten = 0;
  int validSettings = 0;
  int invalidSettings = 0;
  bool wroteHeader = false;
  const std::vector<std::string> globalParamNames = GlobalParamNames();

  const std::vector<std::string> orderedSettings = {"A", "B"};
  for (const auto& setting : orderedSettings) {
    auto it = bySettingBin.find(setting);
    if (it == bySettingBin.end()) continue;

    SettingFitResult res;
    RunSettingFit(setting, it->second, constrained, res);

    if (!wroteHeader) {
      WriteHeader(fout, globalParamNames);
      wroteHeader = true;
    }

    if (res.fit_valid != 1) {
      ++invalidSettings;
      std::cerr << "[WARN] Setting " << setting << " fit invalid: " << res.fit_note << std::endl;
      continue;
    }

    WriteSettingRows(fout, res, it->second, globalParamNames);
    rowsWritten += static_cast<int>(res.bins.size());
    ++validSettings;
  }

  fout.close();

  std::cout << "[INFO] Wrote: " << outPath << "\n"
            << "[INFO] Mode: " << mode << "\n"
            << "[INFO] Rows written: " << rowsWritten << "\n"
            << "[INFO] Valid setting fits: " << validSettings << "\n"
            << "[INFO] Invalid setting fits: " << invalidSettings << std::endl;
}
