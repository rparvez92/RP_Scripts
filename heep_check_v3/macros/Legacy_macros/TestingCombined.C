// TestingCombined.C
//
// TMinuit combined-mode response-matrix test for 5-pass C/E/D.
//
// Combined parameters:
//   dpe_C, dpe_E, dpe_D,
//   dtheta_e_C, dtheta_e_E, dtheta_e_D,
//   dtheta_p_common,
//   dpp_common
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/TestingCombined.C()'

#include <TMinuit.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::vector<std::string> kSettings = {"C", "E", "D"};
const std::vector<std::string> kVars = {"W", "Em", "Pmz", "Pmy"};
const std::vector<std::string> kParamNames = {
  "dpe_C", "dpe_E", "dpe_D",
  "dtheta_e_C", "dtheta_e_E", "dtheta_e_D",
  "dtheta_p_common", "dpp_common"
};

const double kResponse[4][4] = {
  {-24.63, -10.73,  0.00,  0.00},
  {  0.00,  -6.71,  0.00, -4.73},
  { -4.73,   4.75,  0.00,  4.81},
  {  4.75,   4.73, -4.81,  0.00}
};

const std::vector<double> kZeroInit(8, 0.0);
const std::vector<double> kJulioLikeInit = {
  0.9869, 1.0048, 1.4619,
  -0.8619, -0.6471, -1.1178,
  0.0, 0.0
};

struct ObsRow {
  std::string dataset;
  std::string setting;
  double shms_p_GeV = std::numeric_limits<double>::quiet_NaN();
  double delta_center = std::numeric_limits<double>::quiet_NaN();
  double dp_lo = std::numeric_limits<double>::quiet_NaN();
  double dp_hi = std::numeric_limits<double>::quiet_NaN();
  std::string var;
  double mu = std::numeric_limits<double>::quiet_NaN();
  double mu_err = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double sigma_err = std::numeric_limits<double>::quiet_NaN();
};

struct Dataset {
  std::string name;
  std::vector<ObsRow> rows;
};

struct FitResult {
  std::string method;
  std::string init;
  std::vector<double> params = std::vector<double>(8, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> errors = std::vector<double>(8, std::numeric_limits<double>::quiet_NaN());
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int minuit_status = -999;
  int cov_status = -999;
  double edm = std::numeric_limits<double>::quiet_NaN();
  std::string message;
};

struct Context {
  std::vector<ObsRow> rows;
  bool valid = false;
};

Context gContext;

static std::string Trim(const std::string& s) {
  size_t b = 0;
  while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  size_t e = s.size();
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  if (e > b && s[b] == '"' && s[e - 1] == '"') {
    ++b;
    --e;
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
  bool in_quotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        cur.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (c == ',' && !in_quotes) {
      out.push_back(Trim(cur));
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(Trim(cur));
  return out;
}

static double ToDouble(const std::string& s) {
  try {
    size_t pos = 0;
    const std::string t = Trim(s);
    const double x = std::stod(t, &pos);
    return pos == t.size() ? x : std::numeric_limits<double>::quiet_NaN();
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

static int SettingIndex(const std::string& setting) {
  for (int i = 0; i < static_cast<int>(kSettings.size()); ++i) {
    if (kSettings[i] == setting) return i;
  }
  return -1;
}

static int VarIndex(const std::string& var) {
  for (int i = 0; i < static_cast<int>(kVars.size()); ++i) {
    if (kVars[i] == var) return i;
  }
  return -1;
}

static double Predict(const std::vector<double>& p, const std::string& setting, const std::string& var) {
  const int iset = SettingIndex(setting);
  const int ivar = VarIndex(var);
  if (iset < 0 || ivar < 0) return std::numeric_limits<double>::quiet_NaN();

  const double dpe = p[iset];
  const double dtheta_e = p[3 + iset];
  const double dtheta_p = p[6];
  const double dpp = p[7];
  const double pars[4] = {dtheta_e, dpe, dtheta_p, dpp};

  double pred = 0.0;
  for (int j = 0; j < 4; ++j) pred += kResponse[ivar][j] * pars[j];
  return pred;
}

static double Chi2(const std::vector<ObsRow>& rows, const std::vector<double>& p) {
  double chi2 = 0.0;
  for (const auto& row : rows) {
    const double pred = Predict(p, row.setting, row.var);
    const double resid = row.mu - pred;
    chi2 += resid * resid / (row.mu_err * row.mu_err);
  }
  return chi2;
}

static void FCN(Int_t&, Double_t*, Double_t& fval, Double_t* par, Int_t) {
  if (!gContext.valid) {
    fval = 1e30;
    return;
  }
  std::vector<double> p(8, 0.0);
  for (int i = 0; i < 8; ++i) p[i] = par[i];
  fval = Chi2(gContext.rows, p);
}

static bool HasAllRows(const std::vector<ObsRow>& rows) {
  std::map<std::pair<std::string, std::string>, int> seen;
  for (const auto& row : rows) seen[{row.setting, row.var}]++;
  for (const auto& setting : kSettings) {
    for (const auto& var : kVars) {
      if (seen[{setting, var}] != 1) return false;
    }
  }
  return rows.size() == 12;
}

static std::vector<ObsRow> ReadStandardInput(const std::string& path, const std::string& dataset) {
  std::ifstream fin(path);
  std::vector<ObsRow> rows;
  if (!fin) {
    std::cerr << "[WARN] Cannot open " << path << "\n";
    return rows;
  }
  std::string header_line;
  if (!std::getline(fin, header_line)) return rows;
  const auto headers = SplitCsvLine(header_line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[Lower(headers[i])] = i;
  auto idx = [&](const std::string& name) {
    auto it = col.find(Lower(name));
    return it == col.end() ? -1 : it->second;
  };

  const int i_setting = idx("setting");
  const int i_shms = idx("shms_p_GeV");
  const int i_delta = idx("delta_center");
  const int i_dp_lo = idx("dp_lo");
  const int i_dp_hi = idx("dp_hi");
  const int i_var = idx("variable");
  const int i_mu = idx("mu");
  const int i_mu_err = idx("mu_err");
  const int i_sigma = idx("sigma");
  const int i_sigma_err = idx("sigma_err");
  if (i_setting < 0 || i_shms < 0 || i_delta < 0 || i_dp_lo < 0 ||
      i_dp_hi < 0 || i_var < 0 || i_mu < 0 || i_mu_err < 0 ||
      i_sigma < 0 || i_sigma_err < 0) return rows;

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    auto get = [&](int i) {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    ObsRow row;
    row.dataset = dataset;
    row.setting = get(i_setting);
    row.shms_p_GeV = ToDouble(get(i_shms));
    row.delta_center = ToDouble(get(i_delta));
    row.dp_lo = ToDouble(get(i_dp_lo));
    row.dp_hi = ToDouble(get(i_dp_hi));
    row.var = get(i_var);
    row.mu = ToDouble(get(i_mu));
    row.mu_err = ToDouble(get(i_mu_err));
    row.sigma = ToDouble(get(i_sigma));
    row.sigma_err = ToDouble(get(i_sigma_err));
    if (SettingIndex(row.setting) >= 0 && VarIndex(row.var) >= 0) rows.push_back(row);
  }
  if (!HasAllRows(rows)) {
    std::cerr << "[WARN] " << dataset << " is incomplete; skipping.\n";
    rows.clear();
  }
  return rows;
}

static std::vector<ObsRow> ReadRangeScanInput(const std::string& path, double kSigma, const std::string& dataset) {
  std::ifstream fin(path);
  std::vector<ObsRow> rows;
  if (!fin) {
    std::cerr << "[WARN] Cannot open " << path << "\n";
    return rows;
  }
  std::string header_line;
  if (!std::getline(fin, header_line)) return rows;
  const auto headers = SplitCsvLine(header_line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[Lower(headers[i])] = i;
  auto idx = [&](const std::string& name) {
    auto it = col.find(Lower(name));
    return it == col.end() ? -1 : it->second;
  };

  const int i_k = idx("fit_k_sigma");
  const int i_setting = idx("setting");
  const int i_shms = idx("shms_p_GeV");
  const int i_delta = idx("delta_center");
  const int i_dp_lo = idx("dp_lo");
  const int i_dp_hi = idx("dp_hi");
  const int i_var = idx("variable");
  const int i_mu = idx("mu");
  const int i_mu_err = idx("mu_err");
  const int i_sigma = idx("sigma");
  const int i_sigma_err = idx("sigma_err");
  const int i_valid = idx("fit_valid");
  if (i_k < 0 || i_setting < 0 || i_shms < 0 || i_delta < 0 ||
      i_dp_lo < 0 || i_dp_hi < 0 || i_var < 0 || i_mu < 0 ||
      i_mu_err < 0 || i_sigma < 0 || i_sigma_err < 0 || i_valid < 0) return rows;

  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    auto get = [&](int i) {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    if (get(i_valid) != "1") continue;
    if (std::abs(ToDouble(get(i_k)) - kSigma) > 1e-9) continue;
    ObsRow row;
    row.dataset = dataset;
    row.setting = get(i_setting);
    row.shms_p_GeV = ToDouble(get(i_shms));
    row.delta_center = ToDouble(get(i_delta));
    row.dp_lo = ToDouble(get(i_dp_lo));
    row.dp_hi = ToDouble(get(i_dp_hi));
    row.var = get(i_var);
    row.mu = ToDouble(get(i_mu));
    row.mu_err = ToDouble(get(i_mu_err));
    row.sigma = ToDouble(get(i_sigma));
    row.sigma_err = ToDouble(get(i_sigma_err));
    if (SettingIndex(row.setting) >= 0 && VarIndex(row.var) >= 0) rows.push_back(row);
  }
  if (!HasAllRows(rows)) {
    std::cerr << "[WARN] " << dataset << " is incomplete; skipping.\n";
    rows.clear();
  }
  return rows;
}

static std::vector<Dataset> BuildDatasets() {
  std::vector<Dataset> out;
  const auto julio = ReadStandardInput("results/tables/Julio/TestingJulio_input.csv", "julio_input");
  if (!julio.empty()) out.push_back({"julio_input", julio});
  const auto radwan = ReadStandardInput("results/tables/Julio/TestingRadwan_input.csv", "radwan_input");
  if (!radwan.empty()) out.push_back({"radwan_input", radwan});
  const auto rs15 = ReadRangeScanInput("results/tables/RangeScan/TestingRadwan_RangeScan_input.csv", 1.5, "radwan_rangescan_k1p5");
  if (!rs15.empty()) out.push_back({"radwan_rangescan_k1p5", rs15});
  const auto rs20 = ReadRangeScanInput("results/tables/RangeScan/TestingRadwan_RangeScan_input.csv", 2.0, "radwan_rangescan_k2p0");
  if (!rs20.empty()) out.push_back({"radwan_rangescan_k2p0", rs20});
  return out;
}

static FitResult RunMinuit(const Dataset& dataset,
                           const std::string& algorithm,
                           const std::string& initName,
                           const std::vector<double>& init) {
  gContext.rows = dataset.rows;
  gContext.valid = true;

  TMinuit minuit(8);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(FCN);

  int ierflg = 0;
  double arglist[2] = {1.0, 0.0};
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  for (int i = 0; i < 8; ++i) {
    minuit.mnparm(i, kParamNames[i].c_str(), init[i], 0.1, 0.0, 0.0, ierflg);
  }

  arglist[0] = 30000.0;
  arglist[1] = 1e-6;
  if (algorithm == "simplex") {
    minuit.mnexcm("SIMPLEX", arglist, 2, ierflg);
  } else {
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  }

  FitResult out;
  out.method = "tminuit_" + algorithm;
  out.init = initName;
  out.minuit_status = ierflg;
  out.message = (ierflg == 0 ? "ok" : "minuit_command_failed");

  int hesseStatus = 0;
  minuit.mnexcm("HESSE", arglist, 0, hesseStatus);

  for (int i = 0; i < 8; ++i) {
    double val = std::numeric_limits<double>::quiet_NaN();
    double err = std::numeric_limits<double>::quiet_NaN();
    minuit.GetParameter(i, val, err);
    out.params[i] = val;
    out.errors[i] = err;
  }

  double amin = std::numeric_limits<double>::quiet_NaN();
  double edm = std::numeric_limits<double>::quiet_NaN();
  double errdef = std::numeric_limits<double>::quiet_NaN();
  int nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  out.chi2 = amin;
  out.edm = edm;
  out.cov_status = icstat;
  (void)nparx;
  gContext.valid = false;
  return out;
}

static bool ExactLike(const std::vector<double>& p, double chi2) {
  double maxAbs = 0.0;
  for (double x : p) maxAbs = std::max(maxAbs, std::abs(x));
  return maxAbs > 20.0 || (chi2 < 1e-6 && maxAbs > 5.0);
}

static void WriteResult(std::ofstream& params,
                        std::ofstream& residuals,
                        std::ofstream& summary,
                        const Dataset& dataset,
                        const FitResult& fit) {
  const int nObs = static_cast<int>(dataset.rows.size());
  const int nParams = 8;
  const int ndf = nObs - nParams;
  const bool exact = ExactLike(fit.params, fit.chi2);

  double maxAbsPull = 0.0;
  for (const auto& row : dataset.rows) {
    const double pred = Predict(fit.params, row.setting, row.var);
    const double resid = row.mu - pred;
    const double pull = resid / row.mu_err;
    maxAbsPull = std::max(maxAbsPull, std::abs(pull));
    residuals
      << dataset.name << ","
      << fit.method << ","
      << fit.init << ","
      << row.setting << ","
      << row.shms_p_GeV << ","
      << row.delta_center << ","
      << row.dp_lo << ","
      << row.dp_hi << ","
      << row.var << ","
      << row.mu << ","
      << row.mu_err << ","
      << row.sigma << ","
      << row.sigma_err << ","
      << pred << ","
      << resid << ","
      << pull << "\n";
  }

  for (int i = 0; i < 8; ++i) {
    params
      << dataset.name << ","
      << fit.method << ","
      << fit.init << ","
      << kParamNames[i] << ","
      << fit.params[i] << ","
      << fit.errors[i] << "\n";
  }

  summary
    << dataset.name << ","
    << fit.method << ","
    << fit.init << ","
    << fit.chi2 << ","
    << nObs << ","
    << nParams << ","
    << ndf << ","
    << (ndf > 0 ? fit.chi2 / ndf : std::numeric_limits<double>::quiet_NaN()) << ","
    << maxAbsPull << ","
    << (exact ? 1 : 0) << ","
    << (fit.minuit_status == 0 ? "True" : "False") << ","
    << "\"" << fit.message << "\","
    << fit.minuit_status << ","
    << fit.cov_status << ","
    << fit.edm << "\n";
}

}  // namespace

void TestingCombined(
    const char* outDirC = "results/tables/CombinedTest") {
  const TString outDir(outDirC);
  gSystem->mkdir(outDir, kTRUE);

  const TString paramsPath = outDir + "/CombinedTest_TMinuit_params.csv";
  const TString residualsPath = outDir + "/CombinedTest_TMinuit_residuals.csv";
  const TString summaryPath = outDir + "/CombinedTest_TMinuit_summary.csv";

  std::ofstream params(paramsPath.Data(), std::ios::out);
  std::ofstream residuals(residualsPath.Data(), std::ios::out);
  std::ofstream summary(summaryPath.Data(), std::ios::out);
  params << std::fixed << std::setprecision(6);
  residuals << std::fixed << std::setprecision(6);
  summary << std::fixed << std::setprecision(6);

  params << "dataset,method,init,parameter,value,error\n";
  residuals
    << "dataset,method,init,setting,shms_p_GeV,delta_center,dp_lo,dp_hi,variable,"
    << "mu,mu_err,sigma,sigma_err,predicted,residual,pull\n";
  summary
    << "dataset,method,init,chi2,n_obs,n_params,ndf,chi2_ndf,max_abs_pull,"
    << "exact_inversion_flag,success,message,minuit_status,cov_status,edm\n";

  const auto datasets = BuildDatasets();
  for (const auto& dataset : datasets) {
    for (const auto& initItem : std::vector<std::pair<std::string, std::vector<double>>>{
           {"zero_init", kZeroInit},
           {"julio_like_init", kJulioLikeInit}
         }) {
      WriteResult(params, residuals, summary, dataset,
                  RunMinuit(dataset, "simplex", initItem.first, initItem.second));
      WriteResult(params, residuals, summary, dataset,
                  RunMinuit(dataset, "migrad", initItem.first, initItem.second));
    }
  }

  std::cout << "[INFO] Wrote " << paramsPath << "\n";
  std::cout << "[INFO] Wrote " << residualsPath << "\n";
  std::cout << "[INFO] Wrote " << summaryPath << "\n";
}
