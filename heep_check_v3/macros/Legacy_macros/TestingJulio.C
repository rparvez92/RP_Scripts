/*
  TestingJulio.C

  Single-setting 5-pass response-matrix test using Julio's input CSV and
  TMinuit. This mirrors the scipy test in TestingJulio.py, but runs two Minuit
  algorithms independently per setting:

    - tminuit_simplex
    - tminuit_migrad

  It also writes Julio's reported notebook parameters as a reference block.

  Run from RSIDIS/heep_check_v3:
    root -l -b -q 'macros/TestingJulio.C()'
*/

#include <TMinuit.h>
#include <TSystem.h>
#include <TString.h>

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

const std::vector<std::string> kSettings = {"C", "E", "D"};
const std::vector<std::string> kVars = {"W", "Em", "Pmz", "Pmy"};
const std::vector<std::string> kParamNames = {"dtheta_e", "dpe", "dtheta_p", "dpp"};

const double kResponse[4][4] = {
  {-24.63, -10.73,  0.00,  0.00},
  {  0.00,  -6.71,  0.00, -4.73},
  { -4.73,   4.75,  0.00,  4.81},
  {  4.75,   4.73, -4.81,  0.00}
};

const std::map<std::string, std::vector<double>> kJulioReported = {
  {"C", {-0.8619, 0.9869, -2.4567, 0.4436}},
  {"E", {-0.6471, 1.0048, -1.5902, 1.9413}},
  {"D", {-1.1178, 1.4619, -1.3684, 0.9741}}
};

struct ObsRow {
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

struct SettingData {
  std::string setting;
  double shms_p_GeV = std::numeric_limits<double>::quiet_NaN();
  double delta_center = std::numeric_limits<double>::quiet_NaN();
  double dp_lo = std::numeric_limits<double>::quiet_NaN();
  double dp_hi = std::numeric_limits<double>::quiet_NaN();
  std::map<std::string, ObsRow> by_var;
};

struct FitResult {
  std::string fit_type;
  std::vector<double> params = std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> errors = std::vector<double>(4, std::numeric_limits<double>::quiet_NaN());
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int minuit_status = -999;
  int cov_status = -999;
  double edm = std::numeric_limits<double>::quiet_NaN();
  int n_calls = -1;
  std::string message;
};

struct FitContext {
  SettingData data;
  bool valid = false;
};

FitContext gContext;

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

static double PredictOne(int ivar, const std::vector<double>& p) {
  double pred = 0.0;
  for (int ip = 0; ip < 4; ++ip) pred += kResponse[ivar][ip] * p[ip];
  return pred;
}

static double EvaluateChi2(const SettingData& data, const std::vector<double>& p) {
  double chi2 = 0.0;
  for (int iv = 0; iv < static_cast<int>(kVars.size()); ++iv) {
    const ObsRow& row = data.by_var.at(kVars[iv]);
    const double pred = PredictOne(iv, p);
    const double resid = row.mu - pred;
    chi2 += resid * resid / (row.mu_err * row.mu_err);
  }
  return chi2;
}

static void MinuitFCN(Int_t& /*npar*/, Double_t* /*grad*/, Double_t& fval,
                      Double_t* par, Int_t /*iflag*/) {
  if (!gContext.valid) {
    fval = 1e30;
    return;
  }
  std::vector<double> p(4, 0.0);
  for (int i = 0; i < 4; ++i) p[i] = par[i];
  fval = EvaluateChi2(gContext.data, p);
}

static std::map<std::string, SettingData> ReadInput(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input CSV: " << path << "\n";
    return {};
  }

  std::string header_line;
  if (!std::getline(fin, header_line)) {
    std::cerr << "[ERROR] Empty input CSV: " << path << "\n";
    return {};
  }

  const auto headers = SplitCsvLine(header_line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[Lower(headers[i])] = i;

  auto idx = [&](const std::string& name) -> int {
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
      i_sigma < 0 || i_sigma_err < 0) {
    std::cerr << "[ERROR] Missing required column in " << path << "\n";
    return {};
  }

  std::map<std::string, SettingData> out;
  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    auto get = [&](int i) -> std::string {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    auto getD = [&](int i) -> double {
      double x = std::numeric_limits<double>::quiet_NaN();
      ToDouble(get(i), x);
      return x;
    };

    ObsRow row;
    row.setting = get(i_setting);
    row.shms_p_GeV = getD(i_shms);
    row.delta_center = getD(i_delta);
    row.dp_lo = getD(i_dp_lo);
    row.dp_hi = getD(i_dp_hi);
    row.var = get(i_var);
    row.mu = getD(i_mu);
    row.mu_err = getD(i_mu_err);
    row.sigma = getD(i_sigma);
    row.sigma_err = getD(i_sigma_err);

    SettingData& data = out[row.setting];
    data.setting = row.setting;
    data.shms_p_GeV = row.shms_p_GeV;
    data.delta_center = row.delta_center;
    data.dp_lo = row.dp_lo;
    data.dp_hi = row.dp_hi;
    data.by_var[row.var] = row;
  }

  for (const auto& setting : kSettings) {
    auto it = out.find(setting);
    if (it == out.end()) {
      std::cerr << "[ERROR] Missing setting " << setting << " in " << path << "\n";
      return {};
    }
    for (const auto& var : kVars) {
      if (it->second.by_var.find(var) == it->second.by_var.end()) {
        std::cerr << "[ERROR] Missing " << setting << " " << var << " in " << path << "\n";
        return {};
      }
    }
  }

  return out;
}

static FitResult RunMinuitFit(const SettingData& data, const std::string& fit_type) {
  FitResult res;
  res.fit_type = fit_type;
  gContext.data = data;
  gContext.valid = true;

  TMinuit minuit(4);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(MinuitFCN);

  int ierflg = 0;
  double arglist[2] = {1.0, 0.0};
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  for (int i = 0; i < 4; ++i) {
    minuit.mnparm(i, kParamNames[i].c_str(), 0.0, 0.1, 0.0, 0.0, ierflg);
  }

  if (fit_type == "tminuit_simplex") {
    arglist[0] = 20000.0;
    arglist[1] = 1e-6;
    minuit.mnexcm("SIMPLEX", arglist, 2, ierflg);
  } else {
    arglist[0] = 20000.0;
    arglist[1] = 1e-6;
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  }
  res.minuit_status = ierflg;
  res.message = (ierflg == 0) ? "ok" : "minuit_command_failed";

  int hesse_status = 0;
  minuit.mnexcm("HESSE", arglist, 0, hesse_status);

  for (int i = 0; i < 4; ++i) {
    double val = std::numeric_limits<double>::quiet_NaN();
    double err = std::numeric_limits<double>::quiet_NaN();
    minuit.GetParameter(i, val, err);
    res.params[i] = val;
    res.errors[i] = err;
  }

  double amin = std::numeric_limits<double>::quiet_NaN();
  double edm = std::numeric_limits<double>::quiet_NaN();
  double errdef = std::numeric_limits<double>::quiet_NaN();
  int nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  res.chi2 = amin;
  res.edm = edm;
  res.cov_status = icstat;
  (void)nparx;
  res.n_calls = -1;
  gContext.valid = false;
  return res;
}

static FitResult BuildReferenceFit(const SettingData& data, const std::string& setting) {
  FitResult res;
  res.fit_type = "julio_reported_reference";
  res.params = kJulioReported.at(setting);
  res.errors.assign(4, std::numeric_limits<double>::quiet_NaN());
  res.chi2 = EvaluateChi2(data, res.params);
  res.minuit_status = -1;
  res.cov_status = -1;
  res.edm = std::numeric_limits<double>::quiet_NaN();
  res.n_calls = -1;
  res.message = "copied from Julio notebook; not refit here";
  return res;
}

static void WriteHeader(std::ofstream& fout) {
  fout
    << "fit_type,setting,shms_p_GeV,delta_center,dp_lo,dp_hi,variable,"
    << "mu,mu_err,sigma,sigma_err,predicted,residual,pull,"
    << "dtheta_e,dpe,dtheta_p,dpp,chi2,n_obs,n_params,ndf,success,message,nit,nfev,"
    << "dtheta_e_err,dpe_err,dtheta_p_err,dpp_err,minuit_status,cov_status,edm\n";
}

static void WriteRows(std::ofstream& fout, const SettingData& data, const FitResult& fit) {
  for (int iv = 0; iv < static_cast<int>(kVars.size()); ++iv) {
    const ObsRow& row = data.by_var.at(kVars[iv]);
    const double pred = PredictOne(iv, fit.params);
    const double resid = row.mu - pred;
    const double pull = resid / row.mu_err;

    fout
      << fit.fit_type << ","
      << data.setting << ","
      << data.shms_p_GeV << ","
      << data.delta_center << ","
      << data.dp_lo << ","
      << data.dp_hi << ","
      << row.var << ","
      << row.mu << ","
      << row.mu_err << ","
      << row.sigma << ","
      << row.sigma_err << ","
      << pred << ","
      << resid << ","
      << pull << ","
      << fit.params[0] << ","
      << fit.params[1] << ","
      << fit.params[2] << ","
      << fit.params[3] << ","
      << fit.chi2 << ","
      << 4 << ","
      << 4 << ","
      << 0 << ",";

    if (fit.fit_type == "julio_reported_reference") {
      fout << "reference,";
    } else {
      fout << (fit.minuit_status == 0 ? "True" : "False") << ",";
    }

    fout
      << "\"" << fit.message << "\","
      << fit.n_calls << ","
      << fit.n_calls << ","
      << fit.errors[0] << ","
      << fit.errors[1] << ","
      << fit.errors[2] << ","
      << fit.errors[3] << ","
      << fit.minuit_status << ","
      << fit.cov_status << ","
      << fit.edm << "\n";
  }
}

static void PrintSummary(const std::string& outCsv,
                         const std::map<std::string, std::vector<FitResult>>& results) {
  std::cout << "[INFO] Wrote output: " << outCsv << "\n";
  for (const auto& fit_type : {"tminuit_simplex", "tminuit_migrad", "julio_reported_reference"}) {
    std::cout << "\n" << fit_type << "\n";
    for (const auto& setting : kSettings) {
      const auto& fits = results.at(setting);
      auto it = std::find_if(fits.begin(), fits.end(), [&](const FitResult& r) {
        return r.fit_type == fit_type;
      });
      if (it == fits.end()) continue;
      double max_pull = 0.0;
      const SettingData& data = gContext.data;
      (void)data;
      std::cout
        << "  " << setting
        << ": dtheta_e=" << it->params[0]
        << ", dpe=" << it->params[1]
        << ", dtheta_p=" << it->params[2]
        << ", dpp=" << it->params[3]
        << ", chi2=" << it->chi2
        << ", status=" << it->minuit_status
        << ", cov_status=" << it->cov_status
        << "\n";
    }
  }
}

} // namespace

void TestingJulio(const char* inCsvC = "results/tables/Julio/TestingJulio_input.csv",
                  const char* outCsvC = "results/tables/Julio/TestingJulioTMinuit_output.csv") {
  const std::string inCsv = inCsvC ? inCsvC : "";
  const std::string outCsv = outCsvC ? outCsvC : "";
  const auto bySetting = ReadInput(inCsv);
  if (bySetting.empty()) return;

  TString outPath(outCsv.c_str());
  TString parent = gSystem->DirName(outPath);
  if (!parent.IsNull() && gSystem->AccessPathName(parent)) {
    gSystem->mkdir(parent, kTRUE);
  }

  std::ofstream fout(outCsv);
  if (!fout) {
    std::cerr << "[ERROR] Cannot write output CSV: " << outCsv << "\n";
    return;
  }
  fout.setf(std::ios::fixed);
  fout << std::setprecision(3);
  WriteHeader(fout);

  std::map<std::string, std::vector<FitResult>> allResults;
  for (const auto& setting : kSettings) {
    const SettingData& data = bySetting.at(setting);
    std::vector<FitResult> fits;
    fits.push_back(RunMinuitFit(data, "tminuit_simplex"));
    fits.push_back(RunMinuitFit(data, "tminuit_migrad"));
    fits.push_back(BuildReferenceFit(data, setting));

    for (const auto& fit : fits) WriteRows(fout, data, fit);
    allResults[setting] = fits;
  }

  PrintSummary(outCsv, allResults);
}
