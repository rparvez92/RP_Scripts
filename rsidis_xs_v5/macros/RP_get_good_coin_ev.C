// RP_get_good_coin_ev.C
//
// Determine run-by-run coincidence-time windows and random-subtracted counts.
// Run from rsidis_xs_v5 with:
//   root -l -q 'macros/RP_get_good_coin_ev.C()'
// or
//   root -l -b -q 'macros/RP_get_good_coin_ev.C("bigtable/.../input.csv")'

// The central and random integrations are equal-width RF cells.  The fitted
// sigma is diagnostic; it is not used to set the counting-window width.


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
#include <sstream>
#include <string>
#include <set>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TText.h"
#include "TTree.h"

namespace {

constexpr int kNBins = 400;
constexpr double kHistLow = 0.0;
constexpr double kHistHigh = 100.0;
constexpr const char* kCTBranch = "CTime_ePiCoinTime_ROC2";
constexpr const char* kTreeName = "T";
constexpr const char* kDeltaCoinCuts =
  "(P_gtr_p<=2.7 || P_hgcer_npeSum>1) && P_aero_npeSum>2 && "
  "P_cal_etottracknorm<0.8 && P_gtr_dp>-10 && P_gtr_dp<22 && "
  "H_cer_npeSum>2 && H_cal_etottracknorm>0.8 && abs(H_gtr_dp)<8";
constexpr const char* kFullCoinCuts =
  "(P_gtr_p<=2.7 || P_hgcer_npeSum>1) && P_aero_npeSum>2 && "
  "P_cal_etottracknorm<0.8 && P_gtr_dp>-10 && P_gtr_dp<22 && "
  "H_cer_npeSum>2 && H_cal_etottracknorm>0.8 && abs(H_gtr_dp)<8 && "
  "H_gtr_th>-0.15 && H_gtr_th<0.15 && "
  "H_gtr_ph>-0.10 && H_gtr_ph<0.10 && "
  "P_gtr_th>-0.15 && P_gtr_th<0.15 && "
  "P_gtr_ph>-0.10 && P_gtr_ph<0.10";

const double kNaN = std::numeric_limits<double>::quiet_NaN();
const std::vector<int> kRandomOffsets = {-4, -3, -2, 2, 3, 4};
constexpr const char* kPhase1RootDir =
  "/Volumes/T7/RSIDIS/Phase1/Pass0p1/ROOTfiles/SkimmedData";
constexpr const char* kPhase2RootDir =
  "/Volumes/T7/RSIDIS/Phase2/Passm1/ROOTfiles/SkimmedData";

struct InputRow {
  int run = 0;
  double ebeam = kNaN;
  std::string target;
  std::string runType;
  double hmsP = kNaN;
  double hmsTh = kNaN;
  double shmsP = kNaN;
  double shmsTh = kNaN;
  double x = kNaN;
  double q2 = kNaN;
  double z = kNaN;
  double thpq = kNaN;
  double replayCtmean = kNaN;
  double replayCtsigma = kNaN;
  double ransubcoin = kNaN;
};

struct FitInfo {
  bool ok = false;
  std::string status = "NOT_PROCESSED";
  double mean = kNaN;
  double meanErr = kNaN;
  double sigma = kNaN;
  double sigmaErr = kNaN;
  std::unique_ptr<TF1> function;
};

struct TierCounts {
  Long64_t selectedEntries = 0;
  double nCoin = kNaN;
  int nRandomPeaks = 0;
  double nRandomTotal = kNaN;
  double randomMean = kNaN;
  double goodCoin = kNaN;
  double goodCoinErr = kNaN;
};

struct OutputRow {
  InputRow input;
  std::string ctMethod = "NOT_ASSIGNED";
  std::string ctPosRefStatus = "NOT_APPLICABLE";
  std::vector<int> ctPosRefRuns;
  double ctPosRefCtmeanStddev = kNaN;
  double ctPosRefCtmeanRange = kNaN;
  double rfPeriod = kNaN;
  Long64_t rootEntries = 0;
  TierCounts delta;
  TierCounts full;
  FitInfo fit;
  double ctLow = kNaN;
  double ctHigh = kNaN;
  double leftLow = kNaN;
  double leftHigh = kNaN;
  double rightLow = kNaN;
  double rightHigh = kNaN;
  double fullByDelta = kNaN;
  double ctmeanResidual = kNaN;
  double ransubcoinByRPGoodcoin = kNaN;
  std::vector<std::pair<double, double>> randomWindows;
};

struct DatasetMetadata {
  std::string phase = "UNKNOWN";
  std::string passName = "UNKNOWN";
  std::string runType = "UNKNOWN";
  std::string target = "UNKNOWN";
  std::string charge = "UNKNOWN";
  std::string setting = "UNKNOWN";
  std::string rootDir;
  std::vector<std::string> warnings;

  bool IsValid() const { return warnings.empty(); }
};

std::string Trim(const std::string& value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) return "";
  const auto end = value.find_last_not_of(" \t\r\n");
  return value.substr(begin, end - begin + 1);
}

std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (char c : line) {
    if (c == '"') quoted = !quoted;
    else if (c == ',' && !quoted) {
      fields.push_back(Trim(field));
      field.clear();
    } else field.push_back(c);
  }
  fields.push_back(Trim(field));
  return fields;
}

double ParseDouble(const std::string& text) {
  const std::string trimmed = Trim(text);
  if (trimmed.empty()) return kNaN;
  char* end = nullptr;
  const double value = std::strtod(trimmed.c_str(), &end);
  return end != trimmed.c_str() && *end == '\0' ? value : kNaN;
}

bool ReadInputCSV(const std::string& path, std::vector<InputRow>& rows,
                  bool warnSkippedRuns = true) {
  std::ifstream input(path);
  if (!input) {
    std::cerr << "ERROR: cannot open input CSV: " << path << "\n";
    return false;
  }

  std::string line;
  if (!std::getline(input, line)) return false;
  const auto header = SplitCSV(line);
  std::map<std::string, size_t> column;
  for (size_t i = 0; i < header.size(); ++i) column[ToLower(header[i])] = i;
  for (const char* required : {"run", "ebeam", "target", "run_type",
                               "hms_p", "hms_th", "shms_p", "shms_th",
                               "x", "q2", "z", "thpq", "ctmean", "ctsigma",
                               "ransubcoin"}) {
    if (!column.count(required)) {
      std::cerr << "ERROR: missing input column: " << required << "\n";
      return false;
    }
  }

  while (std::getline(input, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCSV(line);
    auto get = [&](const char* name) -> std::string {
      const size_t i = column[ToLower(name)];
      return i < fields.size() ? fields[i] : "";
    };
    InputRow row;
    row.run = std::stoi(get("run"));
    row.ebeam = ParseDouble(get("ebeam"));
    row.target = get("target");
    row.runType = get("run_type");
    row.hmsP = ParseDouble(get("hms_p"));
    row.hmsTh = ParseDouble(get("hms_th"));
    row.shmsP = ParseDouble(get("shms_p"));
    row.shmsTh = ParseDouble(get("shms_th"));
    row.x = ParseDouble(get("x"));
    row.q2 = ParseDouble(get("q2"));
    row.z = ParseDouble(get("z"));
    row.thpq = ParseDouble(get("thpq"));
    row.replayCtmean = ParseDouble(get("ctmean"));
    row.replayCtsigma = ParseDouble(get("ctsigma"));
    row.ransubcoin = ParseDouble(get("ransubcoin"));
    if (row.runType == "PI-SIDIS") row.runType = "PIMINUS";
    if (row.runType == "PI+SIDIS") row.runType = "PIPLUS";
    if (row.runType != "PIMINUS" && row.runType != "PIPLUS") {
      if (warnSkippedRuns)
        std::cerr << "WARNING: skipping non-PI-SIDIS run " << row.run << "\n";
      continue;
    }
    rows.push_back(row);
  }
  return !rows.empty();
}

bool IsPositron(const InputRow& row) {
  return std::isfinite(row.hmsP) && row.hmsP > 0.0;
}

std::string PhaseBigtablePath(const std::string& phase) {
  if (phase == "phase1") return "bigtable/rsidis_bigtable_phase1.csv";
  if (phase == "phase2") return "bigtable/rsidis_bigtable_phase2.csv";
  return "";
}

double SampleStddev(const std::vector<double>& values) {
  if (values.size() < 2) return kNaN;
  double mean = 0.0;
  for (double value : values) mean += value;
  mean /= values.size();
  double sumSquares = 0.0;
  for (double value : values) {
    const double delta = value - mean;
    sumSquares += delta * delta;
  }
  return std::sqrt(sumSquares / (values.size() - 1));
}

void AssignPositronReference(const InputRow& positron,
                             const std::vector<InputRow>& referenceTable,
                             OutputRow& output) {
  output.ctMethod = "PRIOR_ELEC_AVERAGE";
  output.fit.status = "NOT_APPLICABLE";

  std::vector<const InputRow*> candidates;
  for (const auto& row : referenceTable) {
    if (row.run >= positron.run || !(row.hmsP < 0.0)) continue;
    if (row.ebeam != positron.ebeam) continue;
    if (!std::isfinite(row.replayCtmean) || row.replayCtmean < 47.0 ||
        row.replayCtmean > 55.0) continue;
    if (!std::isfinite(row.replayCtsigma) || !(row.replayCtsigma > 0.0))
      continue;
    candidates.push_back(&row);
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const InputRow* a, const InputRow* b) {
              return a->run > b->run;
            });

  constexpr size_t kRequiredReferences = 5;
  const size_t count = std::min(kRequiredReferences, candidates.size());
  std::vector<double> means;
  std::vector<double> sigmas;
  for (size_t i = 0; i < count; ++i) {
    output.ctPosRefRuns.push_back(candidates[i]->run);
    means.push_back(candidates[i]->replayCtmean);
    sigmas.push_back(candidates[i]->replayCtsigma);
  }
  std::sort(output.ctPosRefRuns.begin(), output.ctPosRefRuns.end());

  if (count != kRequiredReferences) {
    output.ctPosRefStatus = "INSUFFICIENT_PRIOR_ELEC_REFERENCES";
    return;
  }

  const auto meanLimits = std::minmax_element(means.begin(), means.end());
  output.ctPosRefCtmeanRange = *meanLimits.second - *meanLimits.first;
  output.ctPosRefCtmeanStddev = SampleStddev(means);

  output.fit.mean = 0.0;
  output.fit.sigma = 0.0;
  for (double value : means) output.fit.mean += value;
  for (double value : sigmas) output.fit.sigma += value;
  output.fit.mean /= count;
  output.fit.sigma /= count;
  output.fit.meanErr = output.ctPosRefCtmeanStddev / std::sqrt(count);
  output.fit.sigmaErr = SampleStddev(sigmas) / std::sqrt(count);
  output.fit.ok = std::isfinite(output.fit.mean) &&
                  std::isfinite(output.fit.meanErr) &&
                  std::isfinite(output.fit.sigma) &&
                  std::isfinite(output.fit.sigmaErr) &&
                  output.fit.sigma > 0.0;
  if (!output.fit.ok) {
    output.ctPosRefStatus = "INVALID_REFERENCE_STATISTICS";
    return;
  }

  std::vector<std::string> warnings;
  if (output.ctPosRefCtmeanRange > 1.0)
    warnings.push_back("WARNING_VERY_LARGE_CTMEAN_SPREAD");
  else if (output.ctPosRefCtmeanRange > 0.5)
    warnings.push_back("WARNING_LARGE_CTMEAN_SPREAD");
  const int oldestGap = positron.run - output.ctPosRefRuns.front();
  if (oldestGap > 30) warnings.push_back("WARNING_STALE_REFERENCE");

  if (warnings.empty()) {
    output.ctPosRefStatus = "OK";
  } else {
    output.ctPosRefStatus = warnings.front();
    for (size_t i = 1; i < warnings.size(); ++i)
      output.ctPosRefStatus += "+" + warnings[i];
  }
}

double RFPeriod(double ebeam) {
  if (ebeam == 6.4490) return 2.0;
  if (ebeam == 8.5831 || ebeam == 10.6716) return 4.0;
  return kNaN;
}

std::string PassFromEnergy(double ebeam) {
  if (ebeam == 6.4490) return "pass3";
  if (ebeam == 8.5831) return "pass4";
  if (ebeam == 10.6716) return "pass5";
  return "UNKNOWN";
}

std::vector<std::string> SplitPath(const std::string& path) {
  std::vector<std::string> parts;
  std::string part;
  std::istringstream stream(path);
  while (std::getline(stream, part, '/')) {
    if (!part.empty() && part != ".") parts.push_back(part);
  }
  return parts;
}

std::vector<std::string> WrapPath(const std::string& path,
                                  size_t maxCharacters = 72) {
  std::vector<std::string> lines;
  std::string current;
  for (const auto& part : SplitPath(path)) {
    const std::string piece = (current.empty() ? "" : "/") + part;
    if (!current.empty() && current.size() + piece.size() > maxCharacters) {
      lines.push_back(current + "/");
      current = part;
    } else {
      current += piece;
    }
  }
  if (!current.empty()) lines.push_back(current);
  if (lines.empty()) lines.push_back(path);
  return lines;
}

std::string InputStem(const std::string& path) {
  std::string name = gSystem->BaseName(path.c_str());
  const auto dot = name.rfind('.');
  if (dot != std::string::npos) name.erase(dot);
  return name;
}

std::string RootPathForRun(const std::string& rootDir, int run) {
  return rootDir + "/skimmed_coin_replay_production_" +
         std::to_string(run) + "_-1.root";
}

std::string SettingNumber(double value) {
  if (value == -999.0) return "neg999";
  std::ostringstream out;
  out << std::fixed << std::setprecision(10) << value;
  std::string text = out.str();
  while (text.size() > 2 && text.back() == '0') text.pop_back();
  if (!text.empty() && text.back() == '.') text.push_back('0');
  if (!text.empty() && text.front() == '-') text.replace(0, 1, "neg");
  std::replace(text.begin(), text.end(), '.', 'p');
  return text;
}

std::string SettingName(const InputRow& row) {
  return "x" + SettingNumber(row.x) +
         "Q2" + SettingNumber(row.q2) +
         "z" + SettingNumber(row.z) +
         "thpq" + SettingNumber(row.thpq);
}

bool SameNumber(double a, double b) {
  return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) < 1.0e-9;
}

void AddUniqueWarning(std::vector<std::string>& warnings,
                      const std::string& warning) {
  if (std::find(warnings.begin(), warnings.end(), warning) == warnings.end())
    warnings.push_back(warning);
}

DatasetMetadata ValidateDataset(const std::string& inputPath,
                                const std::vector<InputRow>& rows) {
  DatasetMetadata meta;
  const auto parts = SplitPath(inputPath);
  size_t phaseIndex = parts.size();
  for (size_t i = 0; i < parts.size(); ++i) {
    if (parts[i] == "phase1" || parts[i] == "phase2") {
      meta.phase = parts[i];
      phaseIndex = i;
      break;
    }
  }

  if (phaseIndex == parts.size()) {
    AddUniqueWarning(meta.warnings, "Input path does not contain phase1 or phase2");
  } else {
    if (phaseIndex + 1 < parts.size()) meta.passName = parts[phaseIndex + 1];
    if (phaseIndex + 2 < parts.size()) meta.runType = parts[phaseIndex + 2];
    if (phaseIndex + 3 < parts.size()) meta.target = parts[phaseIndex + 3];
    if (phaseIndex + 4 < parts.size()) meta.charge = parts[phaseIndex + 4];
    if (phaseIndex + 5 < parts.size()) meta.setting = parts[phaseIndex + 5];
  }

  if (meta.phase == "phase1") meta.rootDir = kPhase1RootDir;
  else if (meta.phase == "phase2") meta.rootDir = kPhase2RootDir;
  else meta.rootDir = "";

  if (rows.empty()) {
    AddUniqueWarning(meta.warnings, "Input CSV contains no data rows");
    return meta;
  }

  const InputRow& first = rows.front();
  const std::string firstPass = PassFromEnergy(first.ebeam);
  const std::string firstCharge = first.hmsP < 0.0 ? "Elec" : "Pos";
  const std::string firstSetting = SettingName(first);

  if (meta.passName != firstPass)
    AddUniqueWarning(meta.warnings, "Path pass disagrees with Ebeam");
  if (meta.runType != first.runType)
    AddUniqueWarning(meta.warnings, "Path run type disagrees with CSV");
  if (meta.target != first.target)
    AddUniqueWarning(meta.warnings, "Path target disagrees with CSV");
  if (meta.charge != firstCharge)
    AddUniqueWarning(meta.warnings, "Path Elec/Pos disagrees with hms_p sign");
  if (meta.setting != firstSetting)
    AddUniqueWarning(meta.warnings, "Path setting disagrees with x-Q2-z-thpq");

  for (const auto& row : rows) {
    if (PassFromEnergy(row.ebeam) != firstPass)
      AddUniqueWarning(meta.warnings, "Rows contain more than one pass");
    if (row.runType != first.runType)
      AddUniqueWarning(meta.warnings, "Rows contain more than one run type");
    if (row.target != first.target)
      AddUniqueWarning(meta.warnings, "Rows contain more than one target");
    if ((row.hmsP < 0.0 ? "Elec" : "Pos") != firstCharge)
      AddUniqueWarning(meta.warnings, "Rows contain both Elec and Pos");
    if (!SameNumber(row.x, first.x))
      AddUniqueWarning(meta.warnings, "Rows contain more than one x value");
    if (!SameNumber(row.q2, first.q2))
      AddUniqueWarning(meta.warnings, "Rows contain more than one Q2 value");
    if (!SameNumber(row.z, first.z))
      AddUniqueWarning(meta.warnings, "Rows contain more than one z value");
    if (!SameNumber(row.thpq, first.thpq))
      AddUniqueWarning(meta.warnings, "Rows contain more than one thpq value");
  }
  return meta;
}

FitInfo FitPeakTwice(TH1D* hist, int run) {
  FitInfo info;
  if (!hist || hist->GetEntries() <= 0) {
    info.status = "EMPTY_HISTOGRAM";
    return info;
  }

  const double peak = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
  TF1 first(Form("ct_first_%d", run), "gaus", peak - 1.5, peak + 1.5);
  TFitResultPtr firstResult = hist->Fit(&first, "SQNR0");
  if (static_cast<int>(firstResult) != 0) {
    info.status = "FIRST_FIT_FAILED";
    return info;
  }

  const double amplitude1 = first.GetParameter(0);
  const double mean1 = first.GetParameter(1);
  const double sigma1 = std::fabs(first.GetParameter(2));
  if (!(amplitude1 > 0.0) || !std::isfinite(mean1) ||
      !std::isfinite(sigma1) || !(sigma1 > 0.0)) {
    info.status = "FIRST_FIT_INVALID";
    return info;
  }

  info.function.reset(new TF1(Form("ct_second_%d", run), "gaus",
                              mean1 - 2.0 * sigma1, mean1 + 2.0 * sigma1));
  info.function->SetParameters(amplitude1, mean1, sigma1);
  TFitResultPtr secondResult = hist->Fit(info.function.get(), "SQNR0");
  if (static_cast<int>(secondResult) != 0) {
    info.status = "SECOND_FIT_FAILED";
    info.function.reset();
    return info;
  }

  const double amplitude2 = info.function->GetParameter(0);
  info.mean = info.function->GetParameter(1);
  info.meanErr = info.function->GetParError(1);
  info.sigma = std::fabs(info.function->GetParameter(2));
  info.sigmaErr = info.function->GetParError(2);
  if (!(amplitude2 > 0.0) || !std::isfinite(info.mean) ||
      !std::isfinite(info.meanErr) || !std::isfinite(info.sigma) ||
      !std::isfinite(info.sigmaErr) || !(info.sigma > 0.0)) {
    info.status = "SECOND_FIT_INVALID";
    info.function.reset();
    return info;
  }

  info.ok = true;
  info.status = "OK";
  info.function->SetLineColor(kRed + 1);
  info.function->SetLineWidth(3);
  return info;
}

double IntegralHalfOpen(TH1D* hist, double low, double high, double& error) {
  const int lowBin = hist->GetXaxis()->FindFixBin(low);
  const int highBin = hist->GetXaxis()->FindFixBin(high) - 1;
  if (highBin < lowBin) {
    error = 0.0;
    return 0.0;
  }
  return hist->IntegralAndError(lowBin, highBin, error);
}

void CountWindows(TH1D* hist, OutputRow& row, TierCounts& counts,
                  bool saveWindows) {
  const double halfWidth = row.rfPeriod / 2.0;
  row.ctLow = row.fit.mean - halfWidth;
  row.ctHigh = row.fit.mean + halfWidth;

  double coinError = 0.0;
  counts.nCoin = IntegralHalfOpen(hist, row.ctLow, row.ctHigh, coinError);

  double randomVariance = 0.0;
  counts.nRandomTotal = 0.0;
  for (int offset : kRandomOffsets) {
    const double center = row.fit.mean + offset * row.rfPeriod;
    const double low = center - halfWidth;
    const double high = center + halfWidth;
    double error = 0.0;
    const double count = IntegralHalfOpen(hist, low, high, error);
    counts.nRandomTotal += count;
    randomVariance += error * error;
    if (saveWindows) row.randomWindows.emplace_back(low, high);
  }

  counts.nRandomPeaks = static_cast<int>(kRandomOffsets.size());
  if (saveWindows) {
    row.leftLow = row.randomWindows[0].first;
    row.leftHigh = row.randomWindows[2].second;
    row.rightLow = row.randomWindows[3].first;
    row.rightHigh = row.randomWindows[5].second;
  }
  counts.randomMean = counts.nRandomTotal / counts.nRandomPeaks;
  const double randomMeanError =
    std::sqrt(randomVariance) / counts.nRandomPeaks;
  counts.goodCoin = counts.nCoin - counts.randomMean;
  counts.goodCoinErr = std::hypot(coinError, randomMeanError);
}

void FinalizeCounts(OutputRow& row) {
  if (std::isfinite(row.delta.goodCoin) && row.delta.goodCoin != 0.0 &&
      std::isfinite(row.full.goodCoin))
    row.fullByDelta = row.full.goodCoin / row.delta.goodCoin;
  if (row.ctMethod != "PRIOR_ELEC_AVERAGE" &&
      row.input.replayCtmean != -999.0 &&
      std::isfinite(row.input.replayCtmean) && std::isfinite(row.fit.mean)) {
    row.ctmeanResidual = row.input.replayCtmean - row.fit.mean;
  }
  if (row.ctMethod != "PRIOR_ELEC_AVERAGE" &&
      row.input.ransubcoin != -999.0 && std::isfinite(row.input.ransubcoin) &&
      std::isfinite(row.delta.goodCoin) && row.delta.goodCoin != 0.0) {
    row.ransubcoinByRPGoodcoin =
      row.input.ransubcoin / row.delta.goodCoin;
  }
}

void DrawBand(double low, double high, double ymax, Color_t color) {
  TBox* band = new TBox(low, 0.0, high, ymax);
  band->SetFillColorAlpha(color, 0.25);
  band->SetLineColor(color);
  band->Draw("same");
}

TText* AddMetadataLine(TPaveText& box, const std::string& text,
                       double size = 0.025, Color_t color = kBlack,
                       Font_t font = 42) {
  TText* line = box.AddText(text.c_str());
  line->SetTextSize(size);
  line->SetTextColor(color);
  line->SetTextFont(font);
  return line;
}

void PrepareMetadataCanvas(TCanvas* canvas) {
  canvas->Clear();
  canvas->SetMargin(0.03, 0.03, 0.03, 0.03);
  canvas->SetGrid(0, 0);
}

void DrawMetadataPage1(TCanvas* canvas, const DatasetMetadata& meta,
                       const std::vector<InputRow>& rows,
                       const std::string& inputPath,
                       const std::string& pdfPath) {
  PrepareMetadataCanvas(canvas);
  TPaveText box(0.05, 0.03, 0.95, 0.95, "NDC");
  box.SetFillColor(kWhite);
  box.SetBorderSize(1);
  box.SetTextAlign(12);
  box.SetMargin(0.04);

  AddMetadataLine(box, "RP_get_good_coin_ev: Dataset and Selection Metadata",
                  0.034, kBlue + 2, 62);
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Input dataset", 0.029, kBlack, 62);
  AddMetadataLine(box, "Input CSV:", 0.025, kBlack, 42);
  for (const auto& pathLine : WrapPath(inputPath))
    AddMetadataLine(box, "  " + pathLine, 0.020);
  AddMetadataLine(box, Form("Number of runs: %zu", rows.size()));
  AddMetadataLine(box, Form("Phase: %s, Pass: %s, Run Type: %s",
                            meta.phase.c_str(), meta.passName.c_str(),
                            meta.runType.c_str()));
  AddMetadataLine(box, Form("Target: %s, Charge category: %s",
                            meta.target.c_str(), meta.charge.c_str()));
  size_t electronRuns = 0;
  for (const auto& row : rows)
    if (!IsPositron(row)) ++electronRuns;
  AddMetadataLine(box, Form("Timing methods: %zu electron fit run(s), "
                            "%zu positron reference run(s)",
                            electronRuns, rows.size() - electronRuns));
  if (!rows.empty()) {
    const auto& first = rows.front();
    AddMetadataLine(box, Form("x = %.2f, Q^{2} = %.2f GeV^{2}, z = %.2f, "
                              "#theta_{pq} = %.2f#circ",
                              first.x, first.q2, first.z, first.thpq));
  }
  AddMetadataLine(box, "Setting directory: " + meta.setting);
  AddMetadataLine(box, "ROOT directory:", 0.025, kBlack, 42);
  AddMetadataLine(box, "  " + meta.rootDir, 0.021);
  AddMetadataLine(box, "ROOT filename: "
                       "skimmed_coin_replay_production_<run>_-1.root");
  AddMetadataLine(box, Form("Coincidence-time branch: %s", kCTBranch));
  AddMetadataLine(box, " ");

  if (meta.IsValid()) {
    AddMetadataLine(box, "Setting validation: PASS", 0.031, kGreen + 2, 62);
  } else {
    AddMetadataLine(box, "Setting validation: WARNING", 0.031, kRed + 1, 62);
    for (const auto& warning : meta.warnings)
      AddMetadataLine(box, "  - " + warning, 0.022, kRed + 1);
  }
  AddMetadataLine(box, " ");

  AddMetadataLine(box, "Histogram configuration", 0.029, kBlack, 62);
  AddMetadataLine(box, Form("%d bins over [%.0f, %.0f] ns; bin width = %.2f ns",
                            kNBins, kHistLow, kHistHigh,
                            (kHistHigh - kHistLow) / kNBins));
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Baseline PID + Delta-only acceptance", 0.029, kBlack, 62);
  AddMetadataLine(box, "  (P_gtr_p <= 2.7 || P_hgcer_npeSum > 1)", 0.022);
  AddMetadataLine(box, "  && P_aero_npeSum > 2 && P_cal_etottracknorm < 0.8", 0.022);
  AddMetadataLine(box, "  && P_gtr_dp > -10 && P_gtr_dp < 22", 0.022);
  AddMetadataLine(box, "  && H_cer_npeSum > 2", 0.022);
  AddMetadataLine(box, "  && H_cal_etottracknorm > 0.8", 0.022);
  AddMetadataLine(box, "  && abs(H_gtr_dp) < 8", 0.022);
  AddMetadataLine(box, "Full-cut adds target-angle acceptance:", 0.024,
                  kBlack, 62);
  AddMetadataLine(box, "  -0.15 < H_gtr_th < 0.15; "
                       "-0.10 < H_gtr_ph < 0.10", 0.021);
  AddMetadataLine(box, "  -0.15 < P_gtr_th < 0.15; "
                       "-0.10 < P_gtr_ph < 0.10", 0.021);

  box.Draw();
  canvas->Print(pdfPath.c_str());
}

void DrawMetadataPage2(TCanvas* canvas, const std::string& pdfPath) {
  PrepareMetadataCanvas(canvas);
  TPaveText box(0.05, 0.03, 0.95, 0.95, "NDC");
  box.SetFillColor(kWhite);
  box.SetBorderSize(1);
  box.SetTextAlign(12);
  box.SetMargin(0.04);

  AddMetadataLine(box, "RP_get_good_coin_ev: Coincidence-Time Methods",
                  0.034, kBlue + 2, 62);
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Electron runs: two-stage Gaussian fit",
                  0.029, kBlue + 1, 62);
  AddMetadataLine(box, "1. Locate the maximum coincidence-time histogram bin.");
  AddMetadataLine(box, "2. First Gaussian fit: maximum-bin center #pm 1.5 ns.");
  AddMetadataLine(box, "3. Second Gaussian fit: first mean #pm 2 #times first sigma.");
  AddMetadataLine(box, "4. CTmean and CTsigma are taken from the second fit.");
  AddMetadataLine(box, "ROOT file, T tree, and coincidence-time branch must exist.");
  AddMetadataLine(box, "Both fits must succeed; amplitude must be positive.");
  AddMetadataLine(box, "Mean, sigma, and their errors must be finite; sigma > 0.");
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Positron runs: five-prior-electron reference",
                  0.029, kRed + 1, 62);
  AddMetadataLine(box, "Search backward in the phase bigtable and require exactly 5 runs.");
  AddMetadataLine(box, "References: earlier PI+SIDIS/PI-SIDIS electron runs, same exact Ebeam.");
  AddMetadataLine(box, "Each reference requires 47 #leq ctmean #leq 55 ns and ctsigma > 0.");
  AddMetadataLine(box, "CTmean/CTsigma = arithmetic means of the 5 reference values.");
  AddMetadataLine(box, "Errors = sample standard deviation / #sqrt{5}.");
  AddMetadataLine(box, "ctmean range > 0.5/1.0 ns gives large/very-large-spread warning.");
  AddMetadataLine(box, "Positron run - oldest reference run > 30 gives stale warning.");
  AddMetadataLine(box, "Warnings remain usable; fewer than 5 references is invalid.");
  AddMetadataLine(box, "Fit_status = NOT_APPLICABLE for a valid positron reference.");
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Method provenance", 0.029, kBlack, 62);
  AddMetadataLine(box, "CT_method: TWO_STAGE_GAUSSIAN_FIT or PRIOR_ELEC_AVERAGE.");
  AddMetadataLine(box, "CT_POS_ref_* records status, count, runs, ctmean spread and range.");
  AddMetadataLine(box, "Failed or insufficient runs remain in the PDF and CSV.");

  box.Draw();
  canvas->Print(pdfPath.c_str());
}

void DrawMetadataPage3(TCanvas* canvas, const std::string& pdfPath) {
  PrepareMetadataCanvas(canvas);
  TPaveText box(0.05, 0.03, 0.95, 0.95, "NDC");
  box.SetFillColor(kWhite);
  box.SetBorderSize(1);
  box.SetTextAlign(12);
  box.SetMargin(0.04);

  AddMetadataLine(box, "RP_get_good_coin_ev: Counting and Output",
                  0.034, kBlue + 2, 62);
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "RF-cell counting", 0.029, kBlack, 62);
  AddMetadataLine(box, "RF period: 6.4490 GeV #rightarrow 2 ns; "
                       "8.5831/10.6716 GeV #rightarrow 4 ns.");
  AddMetadataLine(box, "Central window: CTmean #pm RFperiodNs/2.");
  AddMetadataLine(box, "Random centers: CTmean + k #times RFperiodNs,");
  AddMetadataLine(box, "  k = {-4, -3, -2, +2, +3, +4}.");
  AddMetadataLine(box, "Each random window has half-width RFperiodNs/2.");
  AddMetadataLine(box, "All integration windows are half-open: [low, high).");
  AddMetadataLine(box, "CTmean, CTsigma, and all RF-window positions are "
                       "calibrated from Delta-only.");
  AddMetadataLine(box, "The same fixed windows are integrated separately for "
                       "Delta-only and Full-cut.");
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Random subtraction", 0.029, kBlack, 62);
  AddMetadataLine(box, "N_rndm_peak = 6");
  AddMetadataLine(box, "Rndm_mean = N_rndm_total / N_rndm_peak");
  AddMetadataLine(box, "RP_Goodcoin_<tier> = N_coin_<tier> - "
                       "Rndm_mean_<tier>");
  AddMetadataLine(box, "Goodcoin_err = sqrt(N_coin + "
                       "N_rndm_total / N_rndm_peak^{2})");
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Initial-replay comparison", 0.029, kBlack, 62);
  AddMetadataLine(box, "ransubcoin = initial replay-script good-coin result.");
  AddMetadataLine(box, "Replay comparison = ransubcoin / "
                       "RP_Goodcoin_delta.");
  AddMetadataLine(box, "Ratio is nan when ransubcoin = -999 or "
                       "RP_Goodcoin_delta = 0.");
  AddMetadataLine(box, "ctmean = initial replay-script fitted coincidence mean.");
  AddMetadataLine(box, "CTmean_residual = ctmean - CTmean (this macro).");
  AddMetadataLine(box, "Positron residual and replay ratio are nan/not applicable.");
  AddMetadataLine(box, " ");
  AddMetadataLine(box, "Summary plots and output", 0.029, kBlack, 62);
  AddMetadataLine(box, "CTmean axis: [40, 60] ns; CTsigma axis: [0, 0.7] ns.");
  AddMetadataLine(box, "Run axes use integer ticks and padded limits.");
  AddMetadataLine(box, "Electron: blue value/red error; positron: red value/blue error.");
  AddMetadataLine(box, "CSV floating-point precision: four decimal places.");

  box.Draw();
  canvas->Print(pdfPath.c_str());
}

std::string JoinRuns(const std::vector<int>& runs) {
  std::ostringstream out;
  for (size_t i = 0; i < runs.size(); ++i) {
    if (i) out << ";";
    out << runs[i];
  }
  return out.str();
}

void DrawRunPage(TCanvas* canvas, TH1D* hist, const OutputRow& row,
                 const std::string& pdfPath) {
  canvas->Clear();
  // Leave enough room for four-/five-digit count labels without clipping the
  // vertical axis title. Canvas pixel dimensions do not change NDC margins.
  canvas->SetLeftMargin(0.15);
  canvas->SetRightMargin(0.04);
  canvas->SetBottomMargin(0.10);
  canvas->SetTopMargin(0.20);
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("CTime_ePiCoinTime_ROC2 (ns)");
  hist->GetYaxis()->SetTitle("Counts / 0.25 ns");
  hist->GetYaxis()->SetTitleOffset(1.65);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->SetStats(false);
  hist->SetMaximum(std::max(1.0, hist->GetMaximum() * 1.25));
  hist->Draw("hist");

  TLatex title;
  title.SetNDC(true);
  title.SetTextAlign(22);
  title.SetTextFont(42);
  title.SetTextSize(0.028);
  title.DrawLatex(0.52, 0.975,
    Form("Run = %d, E_{beam} = %.4f GeV, Target = %s, Run Type = %s",
         row.input.run, row.input.ebeam, row.input.target.c_str(),
         row.input.runType.c_str()));
  title.DrawLatex(0.52, 0.935,
    Form("hms_p = %.4f GeV/c, hms_th = %.4f#circ",
         row.input.hmsP, row.input.hmsTh));
  title.DrawLatex(0.52, 0.895,
    Form("shms_p = %.4f GeV/c, shms_th = %.4f#circ",
         row.input.shmsP, row.input.shmsTh));
  title.DrawLatex(0.52, 0.855,
    Form("x = %.2f, Q^{2} = %.2f GeV^{2}, z = %.2f, #theta_{pq} = %.2f#circ",
         row.input.x, row.input.q2, row.input.z, row.input.thpq));

  TLine* referenceLine = nullptr;
  if (row.fit.ok) {
    const double ymax = hist->GetMaximum();
    DrawBand(row.ctLow, row.ctHigh, ymax, kGreen + 1);
    for (const auto& window : row.randomWindows)
      DrawBand(window.first, window.second, ymax, kBlue);
    hist->Draw("hist same");
    if (row.fit.function) {
      row.fit.function->Draw("same");
    } else if (row.ctMethod == "PRIOR_ELEC_AVERAGE") {
      referenceLine = new TLine(row.fit.mean, 0.0, row.fit.mean, ymax);
      referenceLine->SetLineColor(kRed + 1);
      referenceLine->SetLineWidth(3);
      referenceLine->Draw("same");
    }
  }

  TPaveText info(0.16, 0.34, 0.60, 0.76, "NDC");
  info.SetFillStyle(0);
  info.SetBorderSize(0);
  info.SetTextAlign(12);
  info.SetTextFont(42);
  info.SetTextSize(0.021);
  info.AddText(Form("CT method: %s", row.ctMethod.c_str()));
  info.AddText(Form("Fit status: %s", row.fit.status.c_str()));
  if (row.ctMethod == "PRIOR_ELEC_AVERAGE") {
    info.AddText("CT POS ref status:");
    info.AddText(("  " + row.ctPosRefStatus).c_str());
    info.AddText(Form("CT POS ref count: %zu", row.ctPosRefRuns.size()));
    info.AddText("CT POS ref runs:");
    info.AddText(("  " + JoinRuns(row.ctPosRefRuns)).c_str());
    if (std::isfinite(row.ctPosRefCtmeanRange))
      info.AddText(Form("Ref ctmean stddev/range: %.4f / %.4f ns",
                        row.ctPosRefCtmeanStddev, row.ctPosRefCtmeanRange));
  }
  info.AddText(Form("Tree / Delta / Full entries: %lld / %lld / %lld",
                    row.rootEntries, row.delta.selectedEntries,
                    row.full.selectedEntries));
  if (row.fit.ok) {
    info.AddText(Form("CTmean = %.4f #pm %.4f ns", row.fit.mean, row.fit.meanErr));
    info.AddText(Form("CTsigma = %.4f #pm %.4f ns", row.fit.sigma, row.fit.sigmaErr));
    if (row.ctMethod == "PRIOR_ELEC_AVERAGE") {
      info.AddText("ctmean = untrusted for positron");
      info.AddText("ctmean - CTmean = not applicable");
    } else {
      info.AddText(Form("ctmean = %.4f ns", row.input.replayCtmean));
      if (std::isfinite(row.ctmeanResidual))
        info.AddText(Form("ctmean - CTmean = %.4f ns", row.ctmeanResidual));
      else
        info.AddText("ctmean - CTmean = nan");
    }
    info.AddText(Form("Central: [%.4f, %.4f) ns", row.ctLow, row.ctHigh));
    info.AddText(Form("Delta: N_{coin}=%.0f, N_{rndm,total}=%.0f",
                      row.delta.nCoin, row.delta.nRandomTotal));
    info.AddText(Form("Delta RP_Goodcoin = %.4f #pm %.4f",
                      row.delta.goodCoin, row.delta.goodCoinErr));
    info.AddText(Form("Full: N_{coin}=%.0f, N_{rndm,total}=%.0f",
                      row.full.nCoin, row.full.nRandomTotal));
    info.AddText(Form("Full RP_Goodcoin = %.4f #pm %.4f",
                      row.full.goodCoin, row.full.goodCoinErr));
    if (std::isfinite(row.fullByDelta))
      info.AddText(Form("Full / Delta Goodcoin = %.4f", row.fullByDelta));
  }
  info.AddText(Form("ransubcoin = %.4f", row.input.ransubcoin));
  if (row.ctMethod == "PRIOR_ELEC_AVERAGE")
    info.AddText("ransubcoin / Delta RP_Goodcoin = not applicable");
  else if (std::isfinite(row.ransubcoinByRPGoodcoin))
    info.AddText(Form("ransubcoin / Delta RP_Goodcoin = %.4f",
                      row.ransubcoinByRPGoodcoin));
  else
    info.AddText("ransubcoin / Delta RP_Goodcoin = nan");
  info.Draw();

  TLegend legend(0.61, 0.58, 0.95, 0.70);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  TBox centralLegend(0, 0, 0, 0);
  centralLegend.SetFillColorAlpha(kGreen + 1, 0.25);
  TBox randomLegend(0, 0, 0, 0);
  randomLegend.SetFillColorAlpha(kBlue, 0.25);
  legend.AddEntry(&centralLegend, "Central RF cell", "f");
  legend.AddEntry(&randomLegend, "Random RF cells", "f");
  if (row.fit.function)
    legend.AddEntry(row.fit.function.get(), "Two-stage Gaussian fit", "l");
  else if (referenceLine)
    legend.AddEntry(referenceLine, "Prior-electron CTmean", "l");
  legend.Draw();
  canvas->Print(pdfPath.c_str());
}

double NiceRunSpacing(double span, int maxIntervals = 6) {
  if (!(span > 0.0)) return 1.0;
  const double raw = span / maxIntervals;
  const double magnitude = std::pow(10.0, std::floor(std::log10(raw)));
  const double fraction = raw / magnitude;
  double niceFraction = 1.0;
  if (fraction <= 1.0) niceFraction = 1.0;
  else if (fraction <= 2.0) niceFraction = 2.0;
  else if (fraction <= 5.0) niceFraction = 5.0;
  else niceFraction = 10.0;
  // Run numbers are integers, so sub-run tick spacing is never meaningful.
  return std::max(1.0, niceFraction * magnitude);
}

void DrawSummaryPage(TCanvas* canvas, const std::vector<OutputRow>& rows,
                     bool drawMean, const std::string& pdfPath) {
  canvas->Clear();
  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.12);
  canvas->SetTopMargin(0.10);
  canvas->SetGrid(1, 1);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  std::vector<double> allRuns;
  std::vector<double> electronRun, electronValue, electronError;
  std::vector<double> positronRun, positronValue, positronError;
  for (const auto& row : rows) {
    if (!row.fit.ok) continue;
    const double value = drawMean ? row.fit.mean : row.fit.sigma;
    const double error = drawMean ? row.fit.meanErr : row.fit.sigmaErr;
    allRuns.push_back(row.input.run);
    if (row.ctMethod == "PRIOR_ELEC_AVERAGE") {
      positronRun.push_back(row.input.run);
      positronValue.push_back(value);
      positronError.push_back(error);
    } else {
      electronRun.push_back(row.input.run);
      electronValue.push_back(value);
      electronError.push_back(error);
    }
  }

  double axisMin = 0.0;
  double axisMax = 1.0;
  double spacing = 1.0;
  int intervals = 1;
  if (!allRuns.empty()) {
    const auto limits = std::minmax_element(allRuns.begin(), allRuns.end());
    const double runMin = *limits.first;
    const double runMax = *limits.second;
    spacing = NiceRunSpacing(runMax - runMin);
    const double span = runMax - runMin;
    const double padding = std::max({1.0, 0.05 * span, 0.5 * spacing});
    axisMin = std::floor((runMin - padding) / spacing) * spacing;
    axisMax = std::ceil((runMax + padding) / spacing) * spacing;
    if (axisMin == axisMax) {
      axisMin -= spacing;
      axisMax += spacing;
    }
    intervals = std::max(
      1, static_cast<int>(std::lround((axisMax - axisMin) / spacing)));
  }

  TH1D frame(drawMean ? "h_ctmean_summary_frame" : "h_ctsigma_summary_frame",
             drawMean ? "CTmean VS Run;Run;CTmean (ns)"
                      : "CTsigma VS Run;Run;CTsigma (ns)",
             1, axisMin, axisMax);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(drawMean ? 40.0 : 0.0);
  frame.SetMaximum(drawMean ? 60.0 : 0.7);
  frame.GetXaxis()->SetNdivisions(intervals, false);
  frame.GetXaxis()->SetNoExponent(true);
  frame.GetXaxis()->SetDecimals(false);
  frame.GetXaxis()->SetLabelSize(0.028);
  frame.SetLineColor(kWhite);
  frame.Draw();

  TGraphErrors electronGraph(electronRun.size(), electronRun.data(),
                             electronValue.data(), nullptr,
                             electronError.data());
  electronGraph.SetMarkerStyle(20);
  electronGraph.SetMarkerSize(1.20);
  electronGraph.SetMarkerColor(kBlue + 1);
  electronGraph.SetLineColor(kRed + 1);
  electronGraph.SetLineWidth(2);

  TGraphErrors positronGraph(positronRun.size(), positronRun.data(),
                             positronValue.data(), nullptr,
                             positronError.data());
  positronGraph.SetMarkerStyle(20);
  positronGraph.SetMarkerSize(1.20);
  positronGraph.SetMarkerColor(kRed + 1);
  positronGraph.SetLineColor(kBlue + 1);
  positronGraph.SetLineWidth(2);
  gStyle->SetEndErrorSize(6);
  if (!electronRun.empty()) electronGraph.Draw("P same");
  if (!positronRun.empty()) positronGraph.Draw("P same");

  canvas->Modified();
  canvas->Update();
  canvas->Print(pdfPath.c_str());
}

void WriteNumber(std::ostream& out, double value) {
  if (std::isfinite(value)) out << std::fixed << std::setprecision(4) << value;
  else out << "nan";
}

void WriteOutputCSV(const std::string& path, const std::vector<OutputRow>& rows) {
  std::ofstream out(path);
  out << "Run,Ebeam,Target,Run_type,CT_method,CT_POS_ref_status,"
         "CT_POS_ref_count,CT_POS_ref_runs,CT_POS_ref_ctmean_stddev,"
         "CT_POS_ref_ctmean_range,hms_p,hms_th,shms_p,shms_th,x,Q2,z,thpq,"
         "RFperiodNs,CTmean,CTmean_err,CTsigma,CTsigma_err,ctmean,"
         "CTmean_residual,CT_low,CT_high,N_rndm_peak,Left_rndm_low,"
         "Left_rndm_high,Right_rndm_low,Right_rndm_high,"
         "N_coin_delta,N_rndm_total_delta,Rndm_mean_delta,"
         "RP_Goodcoin_delta,RP_Goodcoin_err_delta,CoinCuts_entries_delta,"
         "N_coin_full,N_rndm_total_full,Rndm_mean_full,"
         "RP_Goodcoin_full,RP_Goodcoin_err_full,CoinCuts_entries_full,"
         "RP_Goodcoin_full_by_delta,ransubcoin,"
         "ransubcoin_by_RP_Goodcoin_delta,Fit_status,Root_entries\n";
  for (const auto& row : rows) {
    out << row.input.run << ',';
    WriteNumber(out, row.input.ebeam); out << ',' << row.input.target << ','
        << row.input.runType << ',';
    out << row.ctMethod << ',' << row.ctPosRefStatus << ','
        << row.ctPosRefRuns.size() << ',' << JoinRuns(row.ctPosRefRuns) << ',';
    WriteNumber(out, row.ctPosRefCtmeanStddev); out << ',';
    WriteNumber(out, row.ctPosRefCtmeanRange); out << ',';
    WriteNumber(out, row.input.hmsP); out << ',';
    WriteNumber(out, row.input.hmsTh); out << ',';
    WriteNumber(out, row.input.shmsP); out << ',';
    WriteNumber(out, row.input.shmsTh); out << ',';
    WriteNumber(out, row.input.x); out << ',';
    WriteNumber(out, row.input.q2); out << ',';
    WriteNumber(out, row.input.z); out << ',';
    WriteNumber(out, row.input.thpq); out << ',';
    WriteNumber(out, row.rfPeriod); out << ',';
    WriteNumber(out, row.fit.mean); out << ',';
    WriteNumber(out, row.fit.meanErr); out << ',';
    WriteNumber(out, row.fit.sigma); out << ',';
    WriteNumber(out, row.fit.sigmaErr); out << ',';
    WriteNumber(out, row.input.replayCtmean); out << ',';
    WriteNumber(out, row.ctmeanResidual); out << ',';
    WriteNumber(out, row.ctLow); out << ',';
    WriteNumber(out, row.ctHigh); out << ',';
    out << row.delta.nRandomPeaks << ',';
    WriteNumber(out, row.leftLow); out << ',';
    WriteNumber(out, row.leftHigh); out << ',';
    WriteNumber(out, row.rightLow); out << ',';
    WriteNumber(out, row.rightHigh); out << ',';
    WriteNumber(out, row.delta.nCoin); out << ',';
    WriteNumber(out, row.delta.nRandomTotal); out << ',';
    WriteNumber(out, row.delta.randomMean); out << ',';
    WriteNumber(out, row.delta.goodCoin); out << ',';
    WriteNumber(out, row.delta.goodCoinErr); out << ','
        << row.delta.selectedEntries << ',';
    WriteNumber(out, row.full.nCoin); out << ',';
    WriteNumber(out, row.full.nRandomTotal); out << ',';
    WriteNumber(out, row.full.randomMean); out << ',';
    WriteNumber(out, row.full.goodCoin); out << ',';
    WriteNumber(out, row.full.goodCoinErr); out << ','
        << row.full.selectedEntries << ',';
    WriteNumber(out, row.fullByDelta); out << ',';
    WriteNumber(out, row.input.ransubcoin); out << ',';
    WriteNumber(out, row.ransubcoinByRPGoodcoin); out << ',' << row.fit.status << ','
        << row.rootEntries << '\n';
  }
}

} // namespace

int RP_get_good_coin_ev(const std::string& inputCSVArgument = "") {
  gStyle->SetOptStat(0);

  std::string inputCSV = Trim(inputCSVArgument);
  if (inputCSV.empty()) {
    std::cout << "Enter input CSV path relative to rsidis_xs_v5:\n> ";
    std::getline(std::cin, inputCSV);
    inputCSV = Trim(inputCSV);
  }
  if (inputCSV.empty()) {
    std::cerr << "ERROR: no input CSV path was provided\n";
    return 1;
  }

  std::vector<InputRow> inputs;
  if (!ReadInputCSV(inputCSV, inputs)) return 1;
  const DatasetMetadata metadata = ValidateDataset(inputCSV, inputs);
  const bool hasPositron = std::any_of(
    inputs.begin(), inputs.end(),
    [](const InputRow& row) { return IsPositron(row); });
  std::vector<InputRow> positronReferenceTable;
  if (hasPositron) {
    const std::string referencePath = PhaseBigtablePath(metadata.phase);
    if (referencePath.empty() ||
        !ReadInputCSV(referencePath, positronReferenceTable, false)) {
      std::cerr << "ERROR: cannot load positron reference bigtable for "
                << metadata.phase << "\n";
      return 1;
    }
  }

  const std::string stem = InputStem(inputCSV);
  const std::string outputPDF =
    "results/PDFs/RP_get_good_coin_ev/" + stem + ".pdf";
  const std::string outputCSV =
    "results/Tables/RP_get_good_coin_ev/" + stem + ".csv";
  gSystem->mkdir(gSystem->DirName(outputPDF.c_str()), true);
  gSystem->mkdir(gSystem->DirName(outputCSV.c_str()), true);

  std::vector<OutputRow> outputs;
  std::vector<std::unique_ptr<TH1D>> histograms;
  outputs.reserve(inputs.size());
  histograms.reserve(inputs.size());

  for (const auto& input : inputs) {
    OutputRow row;
    row.input = input;
    row.ctMethod = IsPositron(input) ? "PRIOR_ELEC_AVERAGE"
                                     : "TWO_STAGE_GAUSSIAN_FIT";
    if (IsPositron(input)) {
      row.ctPosRefStatus = "NOT_EVALUATED";
      AssignPositronReference(input, positronReferenceTable, row);
    }
    row.rfPeriod = RFPeriod(input.ebeam);
    const std::string rootPath = RootPathForRun(metadata.rootDir, input.run);
    std::unique_ptr<TH1D> hist(new TH1D(Form("h_ct_%d", input.run), "",
                                        kNBins, kHistLow, kHistHigh));
    hist->Sumw2();
    std::unique_ptr<TH1D> fullHist(
      new TH1D(Form("h_ct_full_%d", input.run), "",
               kNBins, kHistLow, kHistHigh));
    fullHist->Sumw2();

    if (!std::isfinite(row.rfPeriod)) {
      row.fit.ok = false;
      row.fit.status = "UNKNOWN_EBEAM";
    } else if (gSystem->AccessPathName(rootPath.c_str())) {
      row.fit.ok = false;
      row.fit.status = "FILE_MISSING";
    } else {
      std::unique_ptr<TFile> file(TFile::Open(rootPath.c_str(), "READ"));
      if (!file || file->IsZombie()) {
        row.fit.ok = false;
        row.fit.status = "FILE_OPEN_FAILED";
      } else {
        TTree* tree = dynamic_cast<TTree*>(file->Get(kTreeName));
        if (!tree) {
          row.fit.ok = false;
          row.fit.status = "TREE_MISSING";
        } else if (!tree->GetBranch(kCTBranch)) {
          row.fit.ok = false;
          row.fit.status = "CTIME_BRANCH_MISSING";
        } else if (!tree->GetBranch("H_gtr_th") ||
                   !tree->GetBranch("H_gtr_ph") ||
                   !tree->GetBranch("P_gtr_th") ||
                   !tree->GetBranch("P_gtr_ph")) {
          row.fit.ok = false;
          row.fit.status = "FULL_CUT_BRANCH_MISSING";
        } else {
          row.rootEntries = tree->GetEntries();
          gROOT->cd();
          hist->SetDirectory(gROOT);
          fullHist->SetDirectory(gROOT);
          row.delta.selectedEntries =
            tree->Project(hist->GetName(), kCTBranch, kDeltaCoinCuts);
          row.full.selectedEntries =
            tree->Project(fullHist->GetName(), kCTBranch, kFullCoinCuts);
          if (row.delta.selectedEntries < 0 || row.full.selectedEntries < 0) {
            row.fit.ok = false;
            row.fit.status = "CUT_PROJECTION_FAILED";
          } else {
            if (!IsPositron(input))
              row.fit = FitPeakTwice(hist.get(), input.run);
            if (row.fit.ok) {
              CountWindows(hist.get(), row, row.delta, true);
              CountWindows(fullHist.get(), row, row.full, false);
              FinalizeCounts(row);
            }
          }
        }
      }
    }
    histograms.push_back(std::move(hist));
    outputs.push_back(std::move(row));
  }

  // A square drawable canvas keeps every PDF page matched to the PNG aspect.
  TCanvas canvas("c_rp_good_coin", "RP good coincidence", 1200, 1200);
  canvas.SetCanvasSize(1200, 1200);
  canvas.Print((outputPDF + "[").c_str());
  DrawMetadataPage1(&canvas, metadata, inputs, inputCSV, outputPDF);
  DrawMetadataPage2(&canvas, outputPDF);
  DrawMetadataPage3(&canvas, outputPDF);
  for (size_t i = 0; i < outputs.size(); ++i)
    DrawRunPage(&canvas, histograms[i].get(), outputs[i], outputPDF);
  DrawSummaryPage(&canvas, outputs, true, outputPDF);
  DrawSummaryPage(&canvas, outputs, false, outputPDF);
  canvas.Print((outputPDF + "]").c_str());

  WriteOutputCSV(outputCSV, outputs);
  std::cout << "Setting validation: "
            << (metadata.IsValid() ? "PASS" : "WARNING") << "\n";
  for (const auto& warning : metadata.warnings)
    std::cout << "  - " << warning << "\n";
  std::cout << "Wrote " << outputPDF << "\nWrote " << outputCSV << "\n";
  return 0;
}
