/*
  SimultaneousDpeDthe.C

  Purpose
  -------
  Read setting-wise measured observable offsets from a pass/study-specific CSV,
  for example:
    results/tables/MeasuredOffsetsBySetting_5pass.csv
  and perform a simultaneous fit, per setting, under the hypothesis that:

    - dthe is bin-dependent
    - dpe   is bin-dependent
    - dthp  is common to all fitted narrow bins in that setting
    - dpp   is common to all fitted narrow bins in that setting

  Parameter ordering
  ------------------
  For consistency in code, comments, and covariance output, the fit parameter
  ordering is:

    dthe_b1, dthe_b2, dthe_b3, dthe_b4, [dthe_b5 when present],
    dpe_b1,  dpe_b2,  dpe_b3,  dpe_b4,  [dpe_b5 when present],
    dthp,
    dpp

  Bin usage
  ---------
    Setting A: b1, b2, b3, b4
    Setting B: b1, b2, b3, b4, b5

  We do not use:
    - full bins
    - any row with fit_valid != 1
    - Setting A, b5

  Objective function
  ------------------
    chi2 = sum_i [ (y_i(measured) - y_i(predicted))^2 / sigma_i^2 ]

  where:
    - y_i(measured) is offset_MeV
    - sigma_i       is offset_err_MeV
    - y_i(predicted) is given by the setting-dependent linear response model

  Output
  ------
  One row per fitted narrow bin, with the shared fit results repeated on each row.
  The output also includes the full covariance matrix of the simultaneous fit,
  as well as predicted values, predicted errors, residuals, and residual errors.
  The CSV schema is kept fixed across settings by always reserving columns for
  dthe_b1 through dthe_b5, dpe_b1 through dpe_b5, and the corresponding full
  covariance block; fields not used by a given setting are written as nan.

  Default output name:
    results/tables/SimultaneousDpeDthe_5pass.csv

  Run examples:
    root -l -b -q 'macros/SimultaneousDpeDthe.C("results/tables/MeasuredOffsetsBySetting_5pass.csv","")'
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
  int fit_valid = 0;

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
  double dthp = std::numeric_limits<double>::quiet_NaN();
  double dthp_err = std::numeric_limits<double>::quiet_NaN();
  double dpp = std::numeric_limits<double>::quiet_NaN();
  double dpp_err = std::numeric_limits<double>::quiet_NaN();

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
  bool valid = false;
};

static FitContext gFitContext;

static double EvaluateChi2(const std::vector<double>& par);

static void MinuitFCN(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval,
                      Double_t* par, Int_t /*iflag*/) {
  const size_t nparLocal = static_cast<size_t>(gFitContext.dpe_param_index_by_bin.size()) * 2U + 2U;
  std::vector<double> pars(nparLocal, 0.0);
  for (size_t i = 0; i < nparLocal; ++i) pars[i] = par[i];
  fval = EvaluateChi2(pars);
}

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
  for (const std::string& label : {"a","b","c","d","e","f","g","h","i"}) {
    if (t == label || t == "setting " + label || t == "setting" + label) {
      std::string out = label;
      std::transform(out.begin(), out.end(), out.begin(),
                     [](unsigned char c){ return std::toupper(c); });
      return out;
    }
  }
  return Trim(s);
}

static bool IsKnownSetting(const std::string& setting) {
  return setting == "A" || setting == "B" || setting == "C" || setting == "D" ||
         setting == "E" || setting == "F" || setting == "G" || setting == "H" ||
         setting == "I";
}

static bool IsAllowedBinForSetting(const std::string& setting, const std::string& bin) {
  const std::string b = Lower(Trim(bin));
  if (setting == "A") return b == "b1" || b == "b2" || b == "b3" || b == "b4";
  if (IsKnownSetting(setting)) return b == "b1" || b == "b2" || b == "b3" || b == "b4" || b == "b5";
  return false;
}

static std::vector<std::string> BinOrderForSetting(const std::string& setting) {
  if (setting == "A") return {"b1", "b2", "b3", "b4"};
  if (IsKnownSetting(setting)) return {"b1", "b2", "b3", "b4", "b5"};
  return {};
}

static std::vector<std::string> OrderedSettings() {
  return {"A", "B", "C", "D", "E", "F", "G", "H", "I"};
}

struct ResponseCoefficients {
  double W_dthe = 0.0, W_dpe = 0.0, W_dthp = 0.0, W_dpp = 0.0;
  double Em_dthe = 0.0, Em_dpe = 0.0, Em_dthp = 0.0, Em_dpp = 0.0;
  double Pmz_dthe = 0.0, Pmz_dpe = 0.0, Pmz_dthp = 0.0, Pmz_dpp = 0.0;
  double Pmy_dthe = 0.0, Pmy_dpe = 0.0, Pmy_dthp = 0.0, Pmy_dpp = 0.0;
};

static bool GetResponseCoefficients(const std::string& setting, ResponseCoefficients& c) {
  if (setting == "A") {
    c = {-14.08, -8.62,  0.0,   0.0,
           0.0,  -7.06,  0.0,  -2.10,
          -5.75,  4.10,  0.0,   2.27,
           4.10,  5.75, -2.27,  0.0};
    return true;
  }
  if (setting == "B") {
    c = {-17.33, -8.62,  0.0,   0.0,
           0.0,  -5.66,  0.0,  -3.63,
          -4.30,  3.69,  0.0,   3.74,
           3.69,  4.30, -3.74,  0.0};
    return true;
  }
  if (setting == "C" || setting == "D" || setting == "E" || setting == "F" ||
      setting == "G" || setting == "H" || setting == "I") {
    c = {-24.63, -10.73,  0.0,   0.0,
           0.0,   -6.71,  0.0,  -4.73,
          -4.73,   4.75,  0.0,   4.81,
           4.75,   4.73, -4.81,  0.0};
    return true;
  }
  return false;
}

static std::vector<std::string> GlobalParamNames() {
  return {"dthe_b1", "dthe_b2", "dthe_b3", "dthe_b4", "dthe_b5",
          "dpe_b1", "dpe_b2", "dpe_b3", "dpe_b4", "dpe_b5",
          "dthp", "dpp"};
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

static double EvaluateChi2(const std::vector<double>& par) {
  if (!gFitContext.valid) return 1e30;
  double chi2 = 0.0;
  const int idx_dthp = static_cast<int>(gFitContext.dpe_param_index_by_bin.size()) * 2;
  const int idx_dpp  = idx_dthp + 1;

  for (const auto& eq : gFitContext.equations) {
    const int idx_dthe = gFitContext.dthe_param_index_by_bin[eq.bin_index];
    const int idx_dpe = gFitContext.dpe_param_index_by_bin[eq.bin_index];
    const double pred =
        eq.coeff_dthe * par[idx_dthe] +
        eq.coeff_dpe  * par[idx_dpe] +
        eq.coeff_dthp * par[idx_dthp] +
        eq.coeff_dpp  * par[idx_dpp];
    const double resid = eq.measured - pred;
    chi2 += resid * resid / (eq.sigma * eq.sigma);
  }

  return chi2;
}

static int ToIntOrDefault(const std::string& s, int defVal = 0) {
  const std::string t = Trim(s);
  if (t.empty()) return defVal;
  char* end = nullptr;
  long v = std::strtol(t.c_str(), &end, 10);
  if (!(end && *end == '\0')) return defVal;
  return static_cast<int>(v);
}

static std::map<std::string, std::map<std::string, BinData>>
ReadMeasuredOffsetsBySetting(const std::string& csvPath) {
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
  const int iValue = idx("offset_mev");
  const int iErr = idx("offset_err_mev");
  const int iFitValid = idx("fit_valid");

  if (iSetting < 0 || iDpIdx < 0 || iBin < 0 || iDpLo < 0 || iDpHi < 0 ||
      iVar < 0 || iValue < 0 || iErr < 0 || iFitValid < 0) {
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
    const std::string var = getS(iVar);
    const int fitValid = ToIntOrDefault(getS(iFitValid), 0);

    if (Lower(bin) == "full") continue;
    if (!IsAllowedBinForSetting(setting, bin)) continue;
    if (fitValid != 1) continue;

    BinData& b = out[setting][bin];
    b.setting = setting;
    b.bin = bin;
    b.fit_valid = fitValid;
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

    ResponseCoefficients c;
    if (!GetResponseCoefficients(setting, c)) {
      note = "missing_response_coefficients_for_setting_" + setting;
      return false;
    }

    eqs.push_back({ibin, c.W_dthe,   c.W_dpe,   c.W_dthp,   c.W_dpp,   b.W.value,   b.W.err});
    eqs.push_back({ibin, c.Em_dthe,  c.Em_dpe,  c.Em_dthp,  c.Em_dpp,  b.Em.value,  b.Em.err});
    eqs.push_back({ibin, c.Pmz_dthe, c.Pmz_dpe, c.Pmz_dthp, c.Pmz_dpp, b.Pmz.value, b.Pmz.err});
    eqs.push_back({ibin, c.Pmy_dthe, c.Pmy_dpe, c.Pmy_dthp, c.Pmy_dpp, b.Pmy.value, b.Pmy.err});
	  }

  return true;
}

static void BuildParamNames(const std::vector<std::string>& bins,
                            std::vector<std::string>& names) {
  names.clear();
  for (const auto& bin : bins) names.push_back("dthe_" + bin);
  for (const auto& bin : bins) names.push_back("dpe_" + bin);
  names.push_back("dthp");
  names.push_back("dpp");
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
                          SettingFitResult& res) {
  res = SettingFitResult{};
  res.setting = setting;

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
  for (int i = 0; i < nbin; ++i) {
    gFitContext.dthe_param_index_by_bin[i] = i;
    gFitContext.dpe_param_index_by_bin[i] = nbin + i;
  }
  gFitContext.valid = true;

  const int idx_dthp = 2 * nbin;
  const int idx_dpp = idx_dthp + 1;
  TMinuit minuit(npar);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(MinuitFCN);

  int ierflg = 0;
  double arglist[2] = {1.0, 0.0};
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);
  for (int i = 0; i < npar; ++i) {
    minuit.mnparm(i, res.param_names[i].c_str(), 0.0, 0.1, 0.0, 0.0, ierflg);
  }

  arglist[0] = 20000.0;
  arglist[1] = 1e-6;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  if (ierflg != 0) {
    gFitContext.valid = false;
    res.fit_valid = 0;
    res.fit_note = "tminuit_migrad_failed";
    return false;
  }

  minuit.mnexcm("HESSE", arglist, 0, ierflg);

  res.params.assign(npar, std::numeric_limits<double>::quiet_NaN());
  res.param_errs.assign(npar, std::numeric_limits<double>::quiet_NaN());
  for (int i = 0; i < npar; ++i) {
    double val = std::numeric_limits<double>::quiet_NaN();
    double err = std::numeric_limits<double>::quiet_NaN();
    minuit.GetParameter(i, val, err);
    res.params[i] = val;
    res.param_errs[i] = err;
  }

  res.covariance.assign(npar * npar, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> covFlat(npar * npar, 0.0);
  minuit.mnemat(covFlat.data(), npar);
  res.covariance = covFlat;

  res.dthp = res.params[idx_dthp];
  res.dthp_err = res.param_errs[idx_dthp];
  res.dpp = res.params[idx_dpp];
  res.dpp_err = res.param_errs[idx_dpp];

  double amin = std::numeric_limits<double>::quiet_NaN();
  double edm = std::numeric_limits<double>::quiet_NaN();
  double errdef = std::numeric_limits<double>::quiet_NaN();
  int nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  res.chi2_min = amin;
  res.n_obs_used = static_cast<int>(eqs.size());
  res.ndf = res.n_obs_used - res.n_params;
  if (res.ndf > 0) {
    res.chi2_ndf = res.chi2_min / res.ndf;
    res.fit_prob = TMath::Prob(res.chi2_min, res.ndf);
  }
  res.fit_valid = (icstat >= 1) ? 1 : 0;
  res.fit_note = (res.fit_valid == 1) ? "ok_tminuit" : "tminuit_covariance_unavailable";

  for (int ibin = 0; ibin < nbin; ++ibin) {
    const BinData& b = byBin.at(bins[ibin]);
    const int idx_dthe = ibin;
    const int idx_dpe = nbin + ibin;

    FitBinResult bout;
    bout.bin = bins[ibin];
    bout.dthe = res.params[idx_dthe];
    bout.dthe_err = res.param_errs[idx_dthe];
    bout.dpe = res.params[idx_dpe];
    bout.dpe_err = res.param_errs[idx_dpe];

    bout.dW_measured = b.W.value;      bout.dW_err = b.W.err;
    bout.dEm_measured = b.Em.value;    bout.dEm_err = b.Em.err;
    bout.dPmz_measured = b.Pmz.value;  bout.dPmz_err = b.Pmz.err;
    bout.dPmy_measured = b.Pmy.value;  bout.dPmy_err = b.Pmy.err;

    ResponseCoefficients c;
    GetResponseCoefficients(setting, c);

    bout.dW_pred   = c.W_dthe   * bout.dthe + c.W_dpe   * bout.dpe + c.W_dthp   * res.dthp + c.W_dpp   * res.dpp;
    bout.dEm_pred  = c.Em_dthe  * bout.dthe + c.Em_dpe  * bout.dpe + c.Em_dthp  * res.dthp + c.Em_dpp  * res.dpp;
    bout.dPmz_pred = c.Pmz_dthe * bout.dthe + c.Pmz_dpe * bout.dpe + c.Pmz_dthp * res.dthp + c.Pmz_dpp * res.dpp;
    bout.dPmy_pred = c.Pmy_dthe * bout.dthe + c.Pmy_dpe * bout.dpe + c.Pmy_dthp * res.dthp + c.Pmy_dpp * res.dpp;

    const double vW   = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, c.W_dthe,   c.W_dpe,   c.W_dthp,   c.W_dpp);
    const double vEm  = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, c.Em_dthe,  c.Em_dpe,  c.Em_dthp,  c.Em_dpp);
    const double vPmz = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, c.Pmz_dthe, c.Pmz_dpe, c.Pmz_dthp, c.Pmz_dpp);
    const double vPmy = PredVarianceForBin(res.covariance, npar, idx_dthe, idx_dpe, idx_dthp, idx_dpp, c.Pmy_dthe, c.Pmy_dpe, c.Pmy_dthp, c.Pmy_dpp);
    bout.dW_pred_err   = std::isfinite(vW)   ? std::sqrt(vW)   : std::numeric_limits<double>::quiet_NaN();
    bout.dEm_pred_err  = std::isfinite(vEm)  ? std::sqrt(vEm)  : std::numeric_limits<double>::quiet_NaN();
    bout.dPmz_pred_err = std::isfinite(vPmz) ? std::sqrt(vPmz) : std::numeric_limits<double>::quiet_NaN();
    bout.dPmy_pred_err = std::isfinite(vPmy) ? std::sqrt(vPmy) : std::numeric_limits<double>::quiet_NaN();

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
  fout << "setting,bin,dp_idx,dp_lo,dp_hi,";
  for (const auto& name : param_names) {
    if (name == "dthp" || name == "dpp") continue;
    fout << name << "_fit," << name << "_fit_err,";
  }
  fout << "dthp_fit,dthp_fit_err,dpp_fit,dpp_fit_err,";
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
         << CsvQuoteIfNeeded(bres.bin) << ","
         << b.dp_idx << ","
         << b.dp_lo << ","
         << b.dp_hi << ",";

    for (const auto& name : global_param_names) {
      if (name == "dthp" || name == "dpp") continue;
      fout << ParamValueOrNaN(res, name, false) << ","
           << ParamValueOrNaN(res, name, true) << ",";
    }

    fout << res.dthp << "," << res.dthp_err << ","
         << res.dpp << "," << res.dpp_err << ","
         << bres.dW_measured << "," << bres.dW_err << ","
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

} // namespace

void SimultaneousDpeDthe(const char* inCsv = "results/tables/MeasuredOffsetsBySetting_5pass.csv",
                         const char* outCsv = "") {
  const auto bySettingBin = ReadMeasuredOffsetsBySetting(inCsv);
  if (bySettingBin.empty()) {
    std::cerr << "[ERROR] No valid rows read from " << inCsv << std::endl;
    return;
  }

  std::string outPath = Trim(outCsv ? outCsv : "");
  if (outPath.empty()) {
    outPath = "results/tables/SimultaneousDpeDthe_5pass.csv";
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

	  const std::vector<std::string> orderedSettings = OrderedSettings();
  for (const auto& setting : orderedSettings) {
    auto it = bySettingBin.find(setting);
    if (it == bySettingBin.end()) continue;

    SettingFitResult res;
    RunSettingFit(setting, it->second, res);

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
            << "[INFO] Rows written: " << rowsWritten << "\n"
            << "[INFO] Valid setting fits: " << validSettings << "\n"
            << "[INFO] Invalid setting fits: " << invalidSettings << std::endl;
}
