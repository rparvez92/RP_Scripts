// YieldStability.C
//
// QA-replay coincidence-yield stability check.  Run from rate_dependence_v1:
//   root -l -b -q 'macros/YieldStability.C+()'

// The output is deliberately one row per requested QA run.  Unsupported run
// types and failed runs stay in the table, but only OK rows enter the plot.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace {

constexpr const char *kTreeName = "T";
constexpr const char *kCTBranch = "CTime.ePiCoinTime_ROC2";
constexpr int kNBins = 400;
constexpr double kHistLow = 0.0;
constexpr double kHistHigh = 100.0;
constexpr const char *kDeltaCuts =
    "(P.gtr.p<=2.7 || P.hgcer.npeSum>1) && P.aero.npeSum>2 && "
    "P.cal.etottracknorm<0.8 && P.gtr.dp>-10 && P.gtr.dp<22 && "
    "H.cer.npeSum>2 && H.cal.etottracknorm>0.8 && abs(H.gtr.dp)<8";
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const std::vector<int> kRandomOffsets = {-4, -3, -2, 2, 3, 4};

using CSVRow = std::map<std::string, std::string>;

struct CSVTable {
  std::vector<std::string> header;
  std::vector<CSVRow> rows;
};

struct RunResult {
  int run = 0;
  std::string requestedRunType;
  std::string phase;
  std::string bigtablePath;
  std::string runType;
  std::string target;
  std::string rootFile;
  double ebeam = kNaN;
  double hmsP = kNaN;
  double shmsP = kNaN;
  Long64_t rootEntries = -1;
  Long64_t selectedEntries = -1;
  std::string ctMethod;
  std::string referenceRuns;
  double ctMean = kNaN;
  double ctMeanErr = kNaN;
  double ctSigma = kNaN;
  double ctSigmaErr = kNaN;
  double rfPeriod = kNaN;
  double ctLow = kNaN;
  double ctHigh = kNaN;
  double randomLeftLow = kNaN;
  double randomLeftHigh = kNaN;
  double randomRightLow = kNaN;
  double randomRightHigh = kNaN;
  double nCoin = kNaN;
  int nRandomPeaks = 0;
  double nRandomTotal = kNaN;
  double randomMean = kNaN;
  double goodCoin = kNaN;
  double goodCoinErr = kNaN;
  double bcm2Q = kNaN;
  double ps5 = kNaN;
  double ps6 = kNaN;
  std::string activePrescale;
  double prescaleFactor = kNaN;
  double compLivetime = kNaN;
  double hmsTrackingEff = kNaN;
  double shmsTrackingEff = kNaN;
  double triggerEff = 1.0;
  double pidEff = 1.0;
  double boilCorr = kNaN;
  double normFactor = kNaN;
  double normalizedYield = kNaN;
  double normalizedYieldErr = kNaN;
  double bigtableNormalizedYield = kNaN;
  double bigtableNormalizedYieldErr = kNaN;
  std::string status = "NOT_PROCESSED";
  std::string reason;
};

std::string Trim(const std::string &value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

std::vector<std::string> SplitCSV(const std::string &line, bool &valid) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  valid = true;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field += '"';
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (c == ',' && !quoted) {
      fields.push_back(Trim(field));
      field.clear();
    } else {
      field += c;
    }
  }
  valid = !quoted;
  fields.push_back(Trim(field));
  return fields;
}

bool ReadCSV(const std::string &path, CSVTable &table) {
  std::ifstream input(path);
  if (!input) {
    std::cerr << "[ERROR] Cannot open CSV: " << path << '\n';
    return false;
  }
  std::string line;
  bool valid = false;
  if (!std::getline(input, line)) {
    std::cerr << "[ERROR] Empty CSV: " << path << '\n';
    return false;
  }
  table.header = SplitCSV(line, valid);
  if (!valid) return false;
  size_t lineNumber = 1;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (Trim(line).empty()) continue;
    const auto fields = SplitCSV(line, valid);
    if (!valid) {
      std::cerr << "[ERROR] Unclosed quote in " << path << ':' << lineNumber
                << '\n';
      return false;
    }
    CSVRow row;
    for (size_t i = 0; i < table.header.size(); ++i)
      row[Lower(table.header[i])] = i < fields.size() ? fields[i] : "";
    table.rows.push_back(std::move(row));
  }
  return true;
}

std::string Get(const CSVRow &row, const std::string &name) {
  const auto found = row.find(Lower(name));
  return found == row.end() ? "" : found->second;
}

bool HasColumns(const CSVTable &table,
                const std::vector<std::string> &required,
                const std::string &label) {
  std::set<std::string> columns;
  for (const auto &name : table.header) columns.insert(Lower(name));
  bool ok = true;
  for (const auto &name : required) {
    if (!columns.count(Lower(name))) {
      std::cerr << "[ERROR] " << label << " is missing column " << name
                << '\n';
      ok = false;
    }
  }
  return ok;
}

double Number(const std::string &text) {
  const std::string value = Trim(text);
  if (value.empty()) return kNaN;
  char *end = nullptr;
  const double parsed = std::strtod(value.c_str(), &end);
  return end != value.c_str() && end && *end == '\0' ? parsed : kNaN;
}

int RunNumber(const std::string &text) {
  const double value = Number(text);
  return std::isfinite(value) ? static_cast<int>(std::lround(value)) : 0;
}

bool ValidMeasured(double value) {
  return std::isfinite(value) && value != -999.0;
}

bool IsCoinType(const std::string &type) {
  const std::string value = Lower(Trim(type));
  return value == "pi-sidis" || value == "pi+sidis" ||
         value == "piminus" || value == "piplus";
}

bool IsPositron(const CSVRow &row) {
  return Number(Get(row, "hms_p")) > 0.0;
}

void Fail(RunResult &result, const std::string &reason) {
  result.status = "FAILED";
  result.reason = reason;
}

double RFPeriod(double ebeam) {
  // Phase-2 bigtable energies are rounded (for example, 8.581 and 10.676),
  // so match the established nominal energies within 10 MeV.
  constexpr double kEnergyTolerance = 0.010;
  if (std::abs(ebeam - 6.4490) < kEnergyTolerance) return 2.0;
  if (std::abs(ebeam - 8.5831) < kEnergyTolerance ||
      std::abs(ebeam - 10.6716) < kEnergyTolerance) return 4.0;
  return kNaN;
}

std::string JoinRuns(const std::vector<int> &runs) {
  std::ostringstream output;
  for (size_t i = 0; i < runs.size(); ++i) {
    if (i) output << ';';
    output << runs[i];
  }
  return output.str();
}

double SampleStddev(const std::vector<double> &values) {
  if (values.size() < 2) return kNaN;
  double mean = 0.0;
  for (double value : values) mean += value;
  mean /= values.size();
  double sum = 0.0;
  for (double value : values) sum += (value - mean) * (value - mean);
  return std::sqrt(sum / (values.size() - 1));
}

bool AssignPriorElectronReference(const CSVTable &phaseTable,
                                  RunResult &result) {
  const int run = result.run;
  const double ebeam = result.ebeam;
  std::vector<const CSVRow *> candidates;
  for (const auto &row : phaseTable.rows) {
    const int candidateRun = RunNumber(Get(row, "run"));
    const double candidateEbeam = Number(Get(row, "ebeam"));
    const double mean = Number(Get(row, "ctmean"));
    const double sigma = Number(Get(row, "ctsigma"));
    if (candidateRun <= 0 || candidateRun >= run || IsPositron(row) ||
        !IsCoinType(Get(row, "run_type")) ||
        std::abs(candidateEbeam - ebeam) > 1e-9 ||
        !std::isfinite(mean) || mean < 47.0 || mean > 55.0 ||
        !std::isfinite(sigma) || sigma <= 0.0)
      continue;
    candidates.push_back(&row);
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const CSVRow *a, const CSVRow *b) {
              return RunNumber(Get(*a, "run")) > RunNumber(Get(*b, "run"));
            });
  if (candidates.size() < 5) {
    Fail(result, "INSUFFICIENT_PRIOR_ELECTRON_REFERENCES");
    return false;
  }
  candidates.resize(5);
  std::vector<int> runs;
  std::vector<double> means, sigmas;
  for (const auto *row : candidates) {
    runs.push_back(RunNumber(Get(*row, "run")));
    means.push_back(Number(Get(*row, "ctmean")));
    sigmas.push_back(Number(Get(*row, "ctsigma")));
  }
  result.referenceRuns = JoinRuns(runs);
  result.ctMean = 0.0;
  result.ctSigma = 0.0;
  for (double value : means) result.ctMean += value;
  for (double value : sigmas) result.ctSigma += value;
  result.ctMean /= means.size();
  result.ctSigma /= sigmas.size();
  result.ctMeanErr = SampleStddev(means) / std::sqrt(5.0);
  result.ctSigmaErr = SampleStddev(sigmas) / std::sqrt(5.0);
  result.ctMethod = "PRIOR_ELEC_AVERAGE";
  return true;
}

bool FitPeakTwice(TH1D *histogram, RunResult &result) {
  if (!histogram || histogram->GetEntries() <= 0) {
    Fail(result, "EMPTY_CTIME_HISTOGRAM");
    return false;
  }
  const double peak = histogram->GetXaxis()->GetBinCenter(
      histogram->GetMaximumBin());
  TF1 first(Form("ys_ct_first_%d", result.run), "gaus", peak - 1.5,
            peak + 1.5);
  TFitResultPtr firstFit = histogram->Fit(&first, "SQNR0");
  if (static_cast<int>(firstFit) != 0 || !(first.GetParameter(0) > 0.0) ||
      !std::isfinite(first.GetParameter(1)) ||
      !(std::abs(first.GetParameter(2)) > 0.0)) {
    Fail(result, "FIRST_CTIME_FIT_FAILED");
    return false;
  }
  const double mean1 = first.GetParameter(1);
  const double sigma1 = std::abs(first.GetParameter(2));
  TF1 second(Form("ys_ct_second_%d", result.run), "gaus",
             mean1 - 2.0 * sigma1, mean1 + 2.0 * sigma1);
  second.SetParameters(first.GetParameter(0), mean1, sigma1);
  TFitResultPtr secondFit = histogram->Fit(&second, "SQNR0");
  if (static_cast<int>(secondFit) != 0 || !(second.GetParameter(0) > 0.0)) {
    Fail(result, "SECOND_CTIME_FIT_FAILED");
    return false;
  }
  result.ctMean = second.GetParameter(1);
  result.ctMeanErr = second.GetParError(1);
  result.ctSigma = std::abs(second.GetParameter(2));
  result.ctSigmaErr = second.GetParError(2);
  if (!std::isfinite(result.ctMean) || !std::isfinite(result.ctMeanErr) ||
      !std::isfinite(result.ctSigma) || !std::isfinite(result.ctSigmaErr) ||
      result.ctSigma <= 0.0) {
    Fail(result, "SECOND_CTIME_FIT_INVALID");
    return false;
  }
  result.ctMethod = "TWO_STAGE_GAUSSIAN_FIT";
  return true;
}

double IntegralHalfOpen(TH1D *histogram, double low, double high,
                        double &error) {
  const int lowBin = histogram->GetXaxis()->FindFixBin(low);
  const int highBin = histogram->GetXaxis()->FindFixBin(high) - 1;
  if (highBin < lowBin) {
    error = 0.0;
    return 0.0;
  }
  return histogram->IntegralAndError(lowBin, highBin, error);
}

void CountCoincidenceWindows(TH1D *histogram, RunResult &result) {
  const double halfWidth = result.rfPeriod / 2.0;
  result.ctLow = result.ctMean - halfWidth;
  result.ctHigh = result.ctMean + halfWidth;
  double coinError = 0.0;
  result.nCoin = IntegralHalfOpen(histogram, result.ctLow, result.ctHigh,
                                 coinError);
  result.nRandomPeaks = static_cast<int>(kRandomOffsets.size());
  result.nRandomTotal = 0.0;
  double randomVariance = 0.0;
  std::vector<std::pair<double, double>> windows;
  for (int offset : kRandomOffsets) {
    const double center = result.ctMean + offset * result.rfPeriod;
    const double low = center - halfWidth;
    const double high = center + halfWidth;
    double error = 0.0;
    result.nRandomTotal += IntegralHalfOpen(histogram, low, high, error);
    randomVariance += error * error;
    windows.emplace_back(low, high);
  }
  result.randomLeftLow = windows[0].first;
  result.randomLeftHigh = windows[2].second;
  result.randomRightLow = windows[3].first;
  result.randomRightHigh = windows[5].second;
  result.randomMean = result.nRandomTotal / result.nRandomPeaks;
  result.goodCoin = result.nCoin - result.randomMean;
  result.goodCoinErr = std::hypot(
      coinError, std::sqrt(randomVariance) / result.nRandomPeaks);
}

bool LoadNormalization(const CSVRow &row, RunResult &result) {
  result.bcm2Q = Number(Get(row, "bcm2_q"));
  result.ps5 = Number(Get(row, "ps5"));
  result.ps6 = Number(Get(row, "ps6"));
  result.compLivetime = Number(Get(row, "comp_livetime"));
  result.hmsTrackingEff = Number(Get(row, "h_esing_eff"));
  result.shmsTrackingEff = Number(Get(row, "p_hadron_eff"));
  result.boilCorr = Number(Get(row, "boil_corr"));
  result.bigtableNormalizedYield = Number(Get(row, "normyield"));
  result.bigtableNormalizedYieldErr = Number(Get(row, "normyield_err"));
  const bool ps5Positive = ValidMeasured(result.ps5) && result.ps5 > 0.0;
  const bool ps6Positive = ValidMeasured(result.ps6) && result.ps6 > 0.0;
  if (ps5Positive == ps6Positive) {
    Fail(result, "ACTIVE_PRESCALE_NOT_UNIQUE");
    return false;
  }
  result.activePrescale = ps5Positive ? "ps5" : "ps6";
  result.prescaleFactor = ps5Positive ? result.ps5 : result.ps6;
  if (!(ValidMeasured(result.bcm2Q) && result.bcm2Q > 0.0))
    Fail(result, "INVALID_BCM2_Q");
  else if (!(ValidMeasured(result.compLivetime) &&
             result.compLivetime > 0.0))
    Fail(result, "INVALID_COMP_LIVETIME");
  else if (!(ValidMeasured(result.hmsTrackingEff) &&
             result.hmsTrackingEff > 0.0 && result.hmsTrackingEff <= 1.0))
    Fail(result, "INVALID_HMS_TRACKING_EFF");
  else if (!(ValidMeasured(result.shmsTrackingEff) &&
             result.shmsTrackingEff > 0.0 && result.shmsTrackingEff <= 1.0))
    Fail(result, "INVALID_SHMS_TRACKING_EFF");
  else if (!(ValidMeasured(result.boilCorr) && result.boilCorr > 0.0))
    Fail(result, "INVALID_BOIL_CORR");
  if (result.status == "FAILED") return false;
  result.normFactor = result.prescaleFactor * result.boilCorr /
      (result.bcm2Q * result.compLivetime * result.hmsTrackingEff *
       result.shmsTrackingEff * result.triggerEff * result.pidEff);
  return true;
}

std::string CSVQuote(const std::string &value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
  std::string escaped = "\"";
  for (char c : value) {
    if (c == '"') escaped += '"';
    escaped += c;
  }
  return escaped + '"';
}

void WriteNumber(std::ostream &output, double value) {
  if (std::isfinite(value)) output << std::setprecision(10) << value;
  else output << "nan";
}

bool WriteCSV(const std::string &path, const std::vector<RunResult> &rows) {
  std::ofstream output(path);
  if (!output) return false;
  output << "run,requested_run_type,phase,bigtable_source,run_type,target,"
            "ebeam,hms_p,shms_p,root_file,root_entries,selected_entries,"
            "ct_method,ct_reference_runs,ctmean,ctmean_err,ctsigma,ctsigma_err,"
            "rf_period_ns,ct_low,ct_high,random_left_low,random_left_high,"
            "random_right_low,random_right_high,n_coin,n_random_peaks,"
            "n_random_total,random_mean,good_coin,good_coin_err,BCM2_Q,ps5,"
            "ps6,active_prescale,prescale_factor,comp_livetime,h_esing_Eff,"
            "p_hadron_Eff,trigger_eff,pid_eff,boil_corr,norm_factor,"
            "normalized_yield,normalized_yield_err,bigtable_normyield,"
            "bigtable_normyield_err,status,reason\n";
  for (const auto &r : rows) {
    output << r.run << ',' << CSVQuote(r.requestedRunType) << ','
           << r.phase << ',' << CSVQuote(r.bigtablePath) << ','
           << CSVQuote(r.runType) << ',' << CSVQuote(r.target) << ',';
    WriteNumber(output, r.ebeam); output << ',';
    WriteNumber(output, r.hmsP); output << ',';
    WriteNumber(output, r.shmsP); output << ',' << CSVQuote(r.rootFile) << ','
      << r.rootEntries << ',' << r.selectedEntries << ',' << r.ctMethod << ','
      << CSVQuote(r.referenceRuns) << ',';
    for (double value : {r.ctMean, r.ctMeanErr, r.ctSigma, r.ctSigmaErr,
                         r.rfPeriod, r.ctLow, r.ctHigh, r.randomLeftLow,
                         r.randomLeftHigh, r.randomRightLow,
                         r.randomRightHigh, r.nCoin}) {
      WriteNumber(output, value); output << ',';
    }
    output << r.nRandomPeaks << ',';
    for (double value : {r.nRandomTotal, r.randomMean, r.goodCoin,
                         r.goodCoinErr, r.bcm2Q, r.ps5, r.ps6}) {
      WriteNumber(output, value); output << ',';
    }
    output << r.activePrescale << ',';
    for (double value : {r.prescaleFactor, r.compLivetime,
                         r.hmsTrackingEff, r.shmsTrackingEff, r.triggerEff,
                         r.pidEff, r.boilCorr, r.normFactor,
                         r.normalizedYield, r.normalizedYieldErr,
                         r.bigtableNormalizedYield,
                         r.bigtableNormalizedYieldErr}) {
      WriteNumber(output, value); output << ',';
    }
    output << r.status << ',' << CSVQuote(r.reason) << '\n';
  }
  return output.good();
}

double NiceRunSpacing(double span) {
  if (!(span > 0.0)) return 1.0;
  const double raw = span / 8.0;
  const double magnitude = std::pow(10.0, std::floor(std::log10(raw)));
  const double fraction = raw / magnitude;
  const double nice = fraction <= 1.0 ? 1.0 : fraction <= 2.0 ? 2.0
                       : fraction <= 5.0 ? 5.0 : 10.0;
  return std::max(1.0, nice * magnitude);
}

struct RunAxisRange {
  double low = 0.0;
  double high = 1.0;
  int divisions = 1;
};

RunAxisRange IntegerRunAxis(const std::vector<double> &runs) {
  RunAxisRange axis;
  if (runs.empty()) return axis;
  const auto limits = std::minmax_element(runs.begin(), runs.end());
  const double span = *limits.second - *limits.first;
  const double spacing = NiceRunSpacing(span);
  const double padding = std::max({1.0, 0.05 * span, 0.5 * spacing});
  axis.low = std::floor((*limits.first - padding) / spacing) * spacing;
  axis.high = std::ceil((*limits.second + padding) / spacing) * spacing;
  if (axis.low == axis.high) {
    axis.low -= spacing;
    axis.high += spacing;
  }
  axis.divisions = std::max(1, static_cast<int>(std::lround(
      (axis.high - axis.low) / spacing)));
  return axis;
}

void AddExtent(std::vector<double> &extents, double value, double error) {
  if (!std::isfinite(value)) return;
  const double safeError = std::isfinite(error) && error >= 0.0 ? error : 0.0;
  extents.push_back(value - safeError);
  extents.push_back(value + safeError);
}

bool DrawPlot(const std::string &path, const std::vector<RunResult> &rows) {
  std::vector<double> calculatedRuns, calculatedYields, calculatedErrors;
  std::vector<double> bigtableRuns, bigtableYields, bigtableErrors;
  std::vector<double> ratioRuns, ratios, ratioErrors;
  std::vector<double> allRuns, extents, ratioExtents;
  for (const auto &row : rows) {
    if (row.status == "OK") {
      calculatedRuns.push_back(row.run);
      calculatedYields.push_back(row.normalizedYield);
      calculatedErrors.push_back(row.normalizedYieldErr);
      allRuns.push_back(row.run);
      AddExtent(extents, row.normalizedYield, row.normalizedYieldErr);
    }
    if (IsCoinType(row.requestedRunType) &&
        ValidMeasured(row.bigtableNormalizedYield) &&
        ValidMeasured(row.bigtableNormalizedYieldErr) &&
        row.bigtableNormalizedYieldErr >= 0.0) {
      bigtableRuns.push_back(row.run);
      bigtableYields.push_back(row.bigtableNormalizedYield);
      bigtableErrors.push_back(row.bigtableNormalizedYieldErr);
      allRuns.push_back(row.run);
      AddExtent(extents, row.bigtableNormalizedYield,
                row.bigtableNormalizedYieldErr);
    }
    if (row.status == "OK" &&
        ValidMeasured(row.bigtableNormalizedYield) &&
        row.bigtableNormalizedYield != 0.0 &&
        ValidMeasured(row.bigtableNormalizedYieldErr) &&
        row.bigtableNormalizedYieldErr >= 0.0) {
      const double ratio =
          row.normalizedYield / row.bigtableNormalizedYield;
      const double ratioError = std::hypot(
          row.normalizedYieldErr / row.bigtableNormalizedYield,
          row.normalizedYield * row.bigtableNormalizedYieldErr /
              (row.bigtableNormalizedYield * row.bigtableNormalizedYield));
      if (std::isfinite(ratio) && std::isfinite(ratioError)) {
        ratioRuns.push_back(row.run);
        ratios.push_back(ratio);
        ratioErrors.push_back(ratioError);
        AddExtent(ratioExtents, ratio, ratioError);
      }
    }
  }
  if (calculatedRuns.empty() || allRuns.empty() || extents.empty() ||
      ratioRuns.empty() || ratioExtents.empty())
    return false;
  const RunAxisRange runAxis = IntegerRunAxis(allRuns);
  const auto yLimits = std::minmax_element(extents.begin(), extents.end());
  const double ySpan = *yLimits.second - *yLimits.first;
  double yPadding = std::max(0.15 * ySpan,
      0.05 * std::max(std::abs(*yLimits.first), std::abs(*yLimits.second)));
  if (!(yPadding > 0.0)) yPadding = 0.1;
  const auto ratioLimits =
      std::minmax_element(ratioExtents.begin(), ratioExtents.end());
  const double ratioSpan = *ratioLimits.second - *ratioLimits.first;
  double ratioPadding = std::max(0.10 * ratioSpan, 0.05);
  if (!(ratioPadding > 0.0)) ratioPadding = 0.05;

  TCanvas canvas("c_yield_stability", "QA yield stability", 1400, 1050);
  TPad upperPad("p_yield_stability", "Normalized yields", 0.0, 0.30, 1.0,
                1.0);
  TPad lowerPad("p_yield_ratio", "Calculated / bigtable", 0.0, 0.0, 1.0,
                0.30);
  upperPad.SetLeftMargin(0.11);
  upperPad.SetRightMargin(0.04);
  upperPad.SetTopMargin(0.08);
  upperPad.SetBottomMargin(0.02);
  upperPad.SetGridx();
  upperPad.SetGridy();
  lowerPad.SetLeftMargin(0.11);
  lowerPad.SetRightMargin(0.04);
  lowerPad.SetTopMargin(0.02);
  lowerPad.SetBottomMargin(0.30);
  lowerPad.SetGridx();
  lowerPad.SetGridy();
  upperPad.Draw();
  lowerPad.Draw();

  upperPad.cd();
  TH1D frame("h_yield_stability_frame",
             "Normalized Yield vs Run;Run;Normalized yield (mC^{-1})", 1,
             runAxis.low, runAxis.high);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(*yLimits.first - yPadding);
  frame.SetMaximum(*yLimits.second + yPadding);
  frame.GetXaxis()->SetNdivisions(runAxis.divisions, false);
  frame.GetXaxis()->SetNoExponent(true);
  frame.GetXaxis()->SetDecimals(false);
  frame.GetXaxis()->SetLabelSize(0.0);
  frame.GetXaxis()->SetTitleSize(0.0);
  frame.GetYaxis()->SetTitleSize(0.050);
  frame.GetYaxis()->SetLabelSize(0.042);
  frame.GetYaxis()->SetTitleOffset(1.00);
  frame.Draw();
  TGraphErrors bigtableGraph(bigtableRuns.size(), bigtableRuns.data(),
                             bigtableYields.data(), nullptr,
                             bigtableErrors.data());
  bigtableGraph.SetMarkerStyle(20);
  bigtableGraph.SetMarkerSize(0.75);
  bigtableGraph.SetMarkerColor(kRed + 1);
  bigtableGraph.SetLineColor(kRed + 1);
  bigtableGraph.SetLineWidth(1);
  if (!bigtableRuns.empty()) bigtableGraph.Draw("P SAME");
  TGraphErrors calculatedGraph(calculatedRuns.size(), calculatedRuns.data(),
                               calculatedYields.data(), nullptr,
                               calculatedErrors.data());
  calculatedGraph.SetMarkerStyle(24);
  calculatedGraph.SetMarkerSize(0.85);
  calculatedGraph.SetMarkerColor(kBlue + 1);
  calculatedGraph.SetLineColor(kBlue + 1);
  calculatedGraph.SetLineWidth(1);
  calculatedGraph.Draw("P SAME");
  TLegend legend(0.62, 0.76, 0.89, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&calculatedGraph, "Calculated normyield", "lep");
  if (!bigtableRuns.empty())
    legend.AddEntry(&bigtableGraph, "Bigtable normyield", "lep");
  legend.Draw();

  lowerPad.cd();
  TH1D ratioFrame(
      "h_yield_stability_ratio_frame",
      ";Run;Calculated / bigtable", 1, runAxis.low, runAxis.high);
  ratioFrame.SetDirectory(nullptr);
  ratioFrame.SetStats(false);
  ratioFrame.SetMinimum(*ratioLimits.first - ratioPadding);
  ratioFrame.SetMaximum(*ratioLimits.second + ratioPadding);
  ratioFrame.GetXaxis()->SetNdivisions(runAxis.divisions, false);
  ratioFrame.GetXaxis()->SetNoExponent(true);
  ratioFrame.GetXaxis()->SetDecimals(false);
  ratioFrame.GetXaxis()->SetTitleSize(0.12);
  ratioFrame.GetXaxis()->SetLabelSize(0.10);
  ratioFrame.GetXaxis()->SetTitleOffset(1.05);
  ratioFrame.GetYaxis()->SetTitleSize(0.10);
  ratioFrame.GetYaxis()->SetLabelSize(0.085);
  ratioFrame.GetYaxis()->SetTitleOffset(0.48);
  ratioFrame.GetYaxis()->SetNdivisions(505);
  ratioFrame.Draw();
  TLine unity(runAxis.low, 1.0, runAxis.high, 1.0);
  unity.SetLineColor(kBlack);
  unity.SetLineStyle(2);
  unity.SetLineWidth(2);
  unity.Draw("SAME");
  TGraphErrors ratioGraph(ratioRuns.size(), ratioRuns.data(), ratios.data(),
                          nullptr, ratioErrors.data());
  ratioGraph.SetMarkerStyle(20);
  ratioGraph.SetMarkerSize(0.65);
  ratioGraph.SetMarkerColor(kBlack);
  ratioGraph.SetLineColor(kBlack);
  ratioGraph.SetLineWidth(1);
  ratioGraph.Draw("P SAME");

  canvas.cd();
  canvas.SaveAs(path.c_str());
  return !gSystem->AccessPathName(path.c_str());
}

std::string TemporaryPath(const std::string &path) {
  const auto dot = path.find_last_of('.');
  return dot == std::string::npos ? path + ".tmp"
                                  : path.substr(0, dot) + ".tmp" + path.substr(dot);
}

}  // namespace

int YieldStability(
    const char *rootDir =
        "/lustre24/expphy/volatile/hallc/c-rsidis/pdbforce/replay/ROOTfiles",
    const char *runList = "bigtable/final_qa_runlist.csv",
    const char *phase1Bigtable = "bigtable/rsidis_bigtable_pass0p1.csv",
    const char *phase2Bigtable = "bigtable/rsidis_bigtable_phase2.csv") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  CSVTable requests, phase1, phase2;
  if (!ReadCSV(runList, requests) || !ReadCSV(phase1Bigtable, phase1) ||
      !ReadCSV(phase2Bigtable, phase2))
    return 1;
  if (!HasColumns(requests, {"run", "run_type"}, "QA run list")) return 1;
  const std::vector<std::string> required = {
      "run", "ebeam", "target", "hms_p", "shms_p", "run_type",
      "BCM2_Q", "ps5", "ps6", "comp_livetime", "h_esing_Eff",
      "p_hadron_Eff", "boil_corr", "ctmean", "ctsigma", "normyield",
      "normyield_err"};
  if (!HasColumns(phase1, required, "phase-1 bigtable") ||
      !HasColumns(phase2, required, "phase-2 bigtable"))
    return 1;

  struct Match { const CSVRow *row; const CSVTable *table; std::string phase; std::string path; };
  std::map<int, std::vector<Match>> byRun;
  for (const auto &row : phase1.rows)
    byRun[RunNumber(Get(row, "run"))].push_back(
        {&row, &phase1, "phase1", phase1Bigtable});
  for (const auto &row : phase2.rows)
    byRun[RunNumber(Get(row, "run"))].push_back(
        {&row, &phase2, "phase2", phase2Bigtable});

  std::set<int> requestedRuns;
  std::vector<RunResult> results;
  results.reserve(requests.rows.size());
  size_t eligible = 0, skipped = 0;
  for (const auto &request : requests.rows) {
    RunResult result;
    result.run = RunNumber(Get(request, "run"));
    result.requestedRunType = Get(request, "run_type");
    if (result.run <= 0 || !requestedRuns.insert(result.run).second) {
      Fail(result, result.run <= 0 ? "INVALID_REQUESTED_RUN"
                                  : "DUPLICATE_REQUESTED_RUN");
      results.push_back(std::move(result));
      continue;
    }
    const auto found = byRun.find(result.run);
    if (found == byRun.end() || found->second.size() != 1) {
      Fail(result, found == byRun.end() ? "BIGTABLE_RUN_NOT_FOUND"
                                       : "BIGTABLE_RUN_DUPLICATED");
      results.push_back(std::move(result));
      continue;
    }
    const Match &match = found->second.front();
    const CSVRow &metadata = *match.row;
    result.phase = match.phase;
    result.bigtablePath = match.path;
    result.runType = Get(metadata, "run_type");
    result.target = Get(metadata, "target");
    result.ebeam = Number(Get(metadata, "ebeam"));
    result.hmsP = Number(Get(metadata, "hms_p"));
    result.shmsP = Number(Get(metadata, "shms_p"));
    if (!IsCoinType(result.requestedRunType)) {
      result.status = "SKIPPED_UNSUPPORTED_RUN_TYPE";
      result.reason = "ONLY_PI_SIDIS_AND_PI_PLUS_SIDIS_ARE_SUPPORTED";
      ++skipped;
      results.push_back(std::move(result));
      continue;
    }
    ++eligible;
    result.rfPeriod = RFPeriod(result.ebeam);
    result.rootFile = std::string(rootDir) + "/coin_replay_production_" +
                      std::to_string(result.run) + "_-1.root";
    if (!IsCoinType(result.runType))
      Fail(result, "BIGTABLE_RUN_TYPE_NOT_COINCIDENCE");
    else if (!std::isfinite(result.rfPeriod))
      Fail(result, "UNSUPPORTED_BEAM_ENERGY");
    else if (!LoadNormalization(metadata, result)) {
      // Stable reason assigned by LoadNormalization.
    } else if (gSystem->AccessPathName(result.rootFile.c_str())) {
      Fail(result, "ROOT_FILE_MISSING");
    } else {
      std::unique_ptr<TFile> file(TFile::Open(result.rootFile.c_str(), "READ"));
      if (!file || file->IsZombie()) {
        Fail(result, "ROOT_FILE_OPEN_FAILED");
      } else {
        TTree *tree = dynamic_cast<TTree *>(file->Get(kTreeName));
        if (!tree) {
          Fail(result, "TREE_MISSING");
        } else if (!tree->GetBranch(kCTBranch)) {
          Fail(result, "CTIME_BRANCH_MISSING");
        } else {
          result.rootEntries = tree->GetEntries();
          gROOT->cd();
          TH1D histogram(Form("h_qa_yield_%d", result.run), "", kNBins,
                         kHistLow, kHistHigh);
          histogram.Sumw2();
          result.selectedEntries =
              tree->Project(histogram.GetName(), kCTBranch, kDeltaCuts);
          if (result.selectedEntries < 0) {
            Fail(result, "CUT_PROJECTION_FAILED");
          } else {
            const bool timingOK = IsPositron(metadata)
                ? AssignPriorElectronReference(*match.table, result)
                : FitPeakTwice(&histogram, result);
            if (timingOK) {
              CountCoincidenceWindows(&histogram, result);
              result.normalizedYield = result.goodCoin * result.normFactor;
              result.normalizedYieldErr = result.goodCoinErr * result.normFactor;
              if (std::isfinite(result.normalizedYield) &&
                  std::isfinite(result.normalizedYieldErr)) {
                result.status = "OK";
                result.reason.clear();
              } else {
                Fail(result, "NORMALIZED_YIELD_INVALID");
              }
            }
          }
        }
      }
    }
    results.push_back(std::move(result));
  }

  std::sort(results.begin(), results.end(),
            [](const RunResult &a, const RunResult &b) { return a.run < b.run; });
  const std::string finalCSV = "results/tables/yield_stability_QA.csv";
  const std::string finalPDF = "results/PDFs/yield_stability_QA.pdf";
  const std::string tempCSV = TemporaryPath(finalCSV);
  const std::string tempPDF = TemporaryPath(finalPDF);
  gSystem->mkdir(gSystem->DirName(finalCSV.c_str()), true);
  gSystem->mkdir(gSystem->DirName(finalPDF.c_str()), true);
  gSystem->Unlink(tempCSV.c_str());
  gSystem->Unlink(tempPDF.c_str());
  if (!WriteCSV(tempCSV, results)) {
    std::cerr << "[ERROR] Cannot write temporary CSV: " << tempCSV << '\n';
    return 1;
  }
  if (!DrawPlot(tempPDF, results)) {
    std::cerr << "[ERROR] No successful runs to plot or PDF write failed.\n";
    gSystem->Unlink(tempCSV.c_str());
    gSystem->Unlink(tempPDF.c_str());
    return 2;
  }
  if (gSystem->Rename(tempCSV.c_str(), finalCSV.c_str()) != 0 ||
      gSystem->Rename(tempPDF.c_str(), finalPDF.c_str()) != 0) {
    std::cerr << "[ERROR] Failed to publish final output files.\n";
    return 1;
  }

  size_t plotted = 0, failed = 0;
  std::map<std::string, size_t> phases;
  for (const auto &result : results) {
    if (result.status == "OK") ++plotted;
    else if (result.status == "FAILED") ++failed;
    if (!result.phase.empty()) ++phases[result.phase];
  }
  std::cout << "[SUMMARY] Requested: " << results.size()
            << ", eligible: " << eligible << ", plotted: " << plotted
            << ", skipped: " << skipped << ", failed: " << failed << '\n'
            << "[SUMMARY] Matched phase1: " << phases["phase1"]
            << ", phase2: " << phases["phase2"] << '\n'
            << "[OUTPUT] " << finalCSV << '\n'
            << "[OUTPUT] " << finalPDF << '\n';
  return plotted == 0 ? 2 : 0;
}
