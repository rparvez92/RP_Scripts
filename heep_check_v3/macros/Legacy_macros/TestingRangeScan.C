// TestingRangeScan.C
//
// TMinuit SIMPLEX/MIGRAD response-matrix fits for RangeScan measured offsets.
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/TestingRangeScan.C()'

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
#include <tuple>
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

struct ObsRow {
  double fit_k_sigma = std::numeric_limits<double>::quiet_NaN();
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

struct GroupData {
  double fit_k_sigma = std::numeric_limits<double>::quiet_NaN();
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
  std::string message;
};

struct FitContext {
  GroupData data;
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

static double PredictOne(int ivar, const std::vector<double>& p) {
  double pred = 0.0;
  for (int ip = 0; ip < 4; ++ip) pred += kResponse[ivar][ip] * p[ip];
  return pred;
}

static double EvaluateChi2(const GroupData& data, const std::vector<double>& p) {
  double chi2 = 0.0;
  for (int iv = 0; iv < static_cast<int>(kVars.size()); ++iv) {
    const ObsRow& row = data.by_var.at(kVars[iv]);
    const double resid = row.mu - PredictOne(iv, p);
    chi2 += resid * resid / (row.mu_err * row.mu_err);
  }
  return chi2;
}

static void MinuitFCN(Int_t&, Double_t*, Double_t& fval, Double_t* par, Int_t) {
  if (!gContext.valid) {
    fval = 1e30;
    return;
  }
  std::vector<double> p(4, 0.0);
  for (int i = 0; i < 4; ++i) p[i] = par[i];
  fval = EvaluateChi2(gContext.data, p);
}

static std::map<std::pair<double, std::string>, GroupData> ReadInput(const std::string& path) {
  std::ifstream fin(path);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input CSV: " << path << "\n";
    return {};
  }

  std::string header_line;
  if (!std::getline(fin, header_line)) return {};
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
  int i_delta = idx("bin_center");
  if (i_delta < 0) i_delta = idx("delta_center");
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
      i_mu_err < 0 || i_sigma < 0 || i_sigma_err < 0) {
    std::cerr << "[ERROR] Missing required column in " << path << "\n";
    return {};
  }

  std::map<std::pair<double, std::string>, GroupData> out;
  std::string line;
  while (std::getline(fin, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCsvLine(line);
    auto get = [&](int i) {
      return (i >= 0 && i < static_cast<int>(fields.size())) ? Trim(fields[i]) : "";
    };
    if (i_valid >= 0 && get(i_valid) != "1") continue;

    ObsRow row;
    row.fit_k_sigma = ToDouble(get(i_k));
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

    const auto key = std::make_pair(row.fit_k_sigma, row.setting);
    GroupData& data = out[key];
    data.fit_k_sigma = row.fit_k_sigma;
    data.setting = row.setting;
    data.shms_p_GeV = row.shms_p_GeV;
    data.delta_center = row.delta_center;
    data.dp_lo = row.dp_lo;
    data.dp_hi = row.dp_hi;
    data.by_var[row.var] = row;
  }

  return out;
}

static FitResult RunMinuit(const GroupData& data, const std::string& algorithm) {
  gContext.data = data;
  gContext.valid = true;

  TMinuit minuit(4);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(MinuitFCN);
  double arglist[10] = {0};
  int ierflg = 0;
  arglist[0] = 1.0;
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  for (int i = 0; i < 4; ++i) {
    minuit.DefineParameter(i, kParamNames[i].c_str(), 0.0, 0.1, 0.0, 0.0);
  }

  arglist[0] = 20000;
  arglist[1] = 1e-6;
  if (algorithm == "simplex") {
    minuit.mnexcm("SIMPLEX", arglist, 2, ierflg);
  } else {
    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  }

  FitResult out;
  out.fit_type = "tminuit_" + algorithm;
  out.minuit_status = ierflg;
  out.message = (ierflg == 0 ? "ok" : "minuit command returned nonzero status");
  int hesse_status = 0;
  minuit.mnexcm("HESSE", arglist, 0, hesse_status);
  for (int i = 0; i < 4; ++i) {
    double val = 0.0, err = 0.0;
    minuit.GetParameter(i, val, err);
    out.params[i] = val;
    out.errors[i] = err;
  }
  double amin = 0.0, edm = 0.0, errdef = 0.0;
  int nvpar = 0, nparx = 0, icstat = 0;
  minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  out.chi2 = amin;
  out.edm = edm;
  out.cov_status = icstat;
  gContext.valid = false;
  return out;
}

static bool ExactLike(const std::vector<double>& p, double chi2) {
  double maxAbs = 0.0;
  for (double x : p) maxAbs = std::max(maxAbs, std::abs(x));
  return maxAbs > 20.0 || (chi2 < 1e-6 && maxAbs > 5.0);
}

static void WriteRows(std::ofstream& out,
                      std::ofstream& summary,
                      const GroupData& data,
                      const FitResult& fit) {
  const bool exact = ExactLike(fit.params, fit.chi2);
  double maxAbsPull = 0.0;
  for (int iv = 0; iv < static_cast<int>(kVars.size()); ++iv) {
    const ObsRow& row = data.by_var.at(kVars[iv]);
    const double pred = PredictOne(iv, fit.params);
    const double resid = row.mu - pred;
    const double pull = resid / row.mu_err;
    maxAbsPull = std::max(maxAbsPull, std::abs(pull));
    out
      << fit.fit_type << ","
      << data.fit_k_sigma << ","
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
      << 0 << ","
      << (fit.minuit_status == 0 ? "True" : "False") << ","
      << "\"" << fit.message << "\","
      << fit.errors[0] << ","
      << fit.errors[1] << ","
      << fit.errors[2] << ","
      << fit.errors[3] << ","
      << fit.minuit_status << ","
      << fit.cov_status << ","
      << fit.edm << ","
      << (exact ? 1 : 0)
      << "\n";
  }

  summary
    << fit.fit_type << ","
    << data.fit_k_sigma << ","
    << data.setting << ","
    << fit.params[0] << ","
    << fit.params[1] << ","
    << fit.params[2] << ","
    << fit.params[3] << ","
    << fit.chi2 << ","
    << maxAbsPull << ","
    << (exact ? 1 : 0) << ","
    << fit.minuit_status << ","
    << fit.cov_status << ","
    << fit.edm << "\n";
}

}  // namespace

void TestingRangeScan(
    const char* inputCsvC = "results/tables/RangeScan/TestingRadwan_RangeScan_input.csv",
    const char* outputCsvC = "results/tables/RangeScan/TestingRadwan_RangeScan_TMinuit_output.csv",
    const char* summaryCsvC = "results/tables/RangeScan/TestingRadwan_RangeScan_TMinuit_summary.csv") {
  const auto groups = ReadInput(inputCsvC);
  if (groups.empty()) return;

  TString outputCsv(outputCsvC);
  TString summaryCsv(summaryCsvC);
  gSystem->mkdir(gSystem->DirName(outputCsv), kTRUE);
  gSystem->mkdir(gSystem->DirName(summaryCsv), kTRUE);

  std::ofstream out(outputCsv.Data(), std::ios::out);
  std::ofstream summary(summaryCsv.Data(), std::ios::out);
  out << std::fixed << std::setprecision(3);
  summary << std::fixed << std::setprecision(3);

  out
    << "fit_type,fit_k_sigma,setting,shms_p_GeV,delta_center,dp_lo,dp_hi,variable,"
    << "mu,mu_err,sigma,sigma_err,predicted,residual,pull,"
    << "dtheta_e,dpe,dtheta_p,dpp,chi2,n_obs,n_params,ndf,success,message,"
    << "dtheta_e_err,dpe_err,dtheta_p_err,dpp_err,minuit_status,cov_status,edm,exact_inversion_flag\n";
  summary
    << "method,fit_k_sigma,setting,dtheta_e,dpe,dtheta_p,dpp,chi2,max_abs_pull,"
    << "exact_inversion_flag,minuit_status,cov_status,edm\n";

  for (const auto& item : groups) {
    const GroupData& data = item.second;
    bool complete = true;
    for (const auto& var : kVars) {
      if (data.by_var.find(var) == data.by_var.end()) complete = false;
    }
    if (!complete) continue;

    WriteRows(out, summary, data, RunMinuit(data, "simplex"));
    WriteRows(out, summary, data, RunMinuit(data, "migrad"));
  }

  std::cout << "[INFO] Wrote " << outputCsv << "\n";
  std::cout << "[INFO] Wrote " << summaryCsv << "\n";
}
