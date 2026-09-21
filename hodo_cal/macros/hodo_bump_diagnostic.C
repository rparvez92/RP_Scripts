// hodo_bump_diagnostic.C
//
// Diagnose localized structure in HMS beta versus focal-plane x by tracing
// track-associated hodoscope hits through every timing-correction stage.
//
// Baseline-only usage:
//   root -l -b -q 'macros/hodo_bump_diagnostic.C+("/path/to/coin_replay_production_28425_-1.root")'
//
// Baseline/trial comparison:
//   root -l -b -q 'macros/hodo_bump_diagnostic.C+("baseline.root","trial.root")'

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace hodo_bump {

constexpr int kRun = 28425;
constexpr double kBumpLo = 12.0;
constexpr double kBumpHi = 18.0;
constexpr double kControlLo = 6.0;
constexpr double kControlHi = 24.0;
constexpr double kInvalid = 1.0e30;

const char *kCoinCut =
    "(H.gtr.dp>-8)&&(H.gtr.dp<8)"
    "&&(H.cal.etottracknorm>0.7)&&(H.cer.npeSum>2.0)"
    "&&(P.gtr.dp>-10)&&(P.gtr.dp<22)"
    "&&(P.cal.etottracknorm<0.8)"
    "&&(((P.gtr.p<2.7)&&(P.aero.npeSum>2))"
    "||((P.gtr.p>=2.7)&&(P.hgcer.npeSum>1)&&(P.aero.npeSum>2)))";

struct RunningStats {
  long long n = 0;
  double sum = 0.0;
  double sum2 = 0.0;
  void Fill(double value) {
    if (!std::isfinite(value) || std::abs(value) >= kInvalid)
      return;
    ++n;
    sum += value;
    sum2 += value * value;
  }
  double Mean() const {
    return n ? sum / static_cast<double>(n)
             : std::numeric_limits<double>::quiet_NaN();
  }
  double Sigma() const {
    if (n < 2)
      return std::numeric_limits<double>::quiet_NaN();
    const double variance =
        (sum2 - sum * sum / static_cast<double>(n)) / (n - 1.0);
    return std::sqrt(std::max(0.0, variance));
  }
};

enum Stage { kRaw, kUncorr, kWalk, kTof, kCorr, kNStages };
const std::array<const char *, kNStages> kStageNames = {
    "raw", "uncorrected", "time_walk", "tof", "corrected"};

struct SideSummary {
  std::array<RunningStats, kNStages> bump;
  std::array<RunningStats, kNStages> control;
  RunningStats bumpAdcAmp;
  RunningStats controlAdcAmp;
  RunningStats bumpAdcInt;
  RunningStats controlAdcInt;
  RunningStats bumpRawAdcAmp;
  RunningStats controlRawAdcAmp;
};

struct PaddleSummary {
  long long selectedBump = 0;
  long long selectedControl = 0;
  long long acceptedBump = 0;
  long long acceptedControl = 0;
  long long bothSidesBump = 0;
  long long bothSidesControl = 0;
  long long posValidBump = 0;
  long long posValidControl = 0;
  long long negValidBump = 0;
  long long negValidControl = 0;
  RunningStats betaBump;
  RunningStats betaControl;
  RunningStats posMinusNegBump;
  RunningStats posMinusNegControl;
  RunningStats posNegMeanBump;
  RunningStats posNegMeanControl;
  RunningStats alongBump;
  RunningStats alongControl;
  SideSummary pos;
  SideSummary neg;
  double score = 0.0;
  std::string candidateStage;
};

struct SideReaders {
  std::unique_ptr<TTreeReaderArray<Double_t>> used;
  std::unique_ptr<TTreeReaderArray<Double_t>> uncorr;
  std::unique_ptr<TTreeReaderArray<Double_t>> walk;
  std::unique_ptr<TTreeReaderArray<Double_t>> tof;
  std::unique_ptr<TTreeReaderArray<Double_t>> corr;
  std::unique_ptr<TTreeReaderArray<Double_t>> adcAmp;
  std::unique_ptr<TTreeReaderArray<Double_t>> adcInt;
  std::unique_ptr<TTreeReaderArray<Double_t>> rawCounter;
  std::unique_ptr<TTreeReaderArray<Double_t>> rawTime;
  std::unique_ptr<TTreeReaderArray<Double_t>> rawAdcCounter;
  std::unique_ptr<TTreeReaderArray<Double_t>> rawAdcAmp;
};

struct PlaneReaders {
  std::string name;
  int paddles = 0;
  SideReaders pos;
  SideReaders neg;
  std::unique_ptr<TTreeReaderValue<Double_t>> trackX;
  std::unique_ptr<TTreeReaderValue<Double_t>> trackY;
  std::unique_ptr<TTreeReaderValue<Double_t>> fpTime;
};

struct PlaneOutput {
  std::string name;
  int paddles = 0;
  std::vector<PaddleSummary> summary;
  std::unique_ptr<TH2D> paddleVsXfp;
  std::unique_ptr<TH2D> betaVsPaddle;
  std::unique_ptr<TH2D> fpTimeVsXfp;
  std::unique_ptr<TH2D> residualVsXfp;
  std::vector<std::unique_ptr<TH1D>> betaBump;
  std::vector<std::unique_ptr<TH1D>> betaControl;
};

struct AnalysisResult {
  std::string label;
  long long allEvents = 0;
  long long selectedEvents = 0;
  long long bumpEvents = 0;
  long long controlEvents = 0;
  std::unique_ptr<TH2D> betaVsXfp;
  std::unique_ptr<TH1D> betaBump;
  std::unique_ptr<TH1D> betaControl;
  std::array<PlaneOutput, 4> planes;
};

std::string Unique(const std::string &base, const std::string &label) {
  std::string result = base + "_" + label;
  std::replace_if(result.begin(), result.end(),
                  [](char c) { return !std::isalnum(c); }, '_');
  return result;
}

template <typename H> void Detach(H *hist) {
  hist->SetDirectory(nullptr);
  hist->SetStats(false);
}

bool HasBranch(TTree *tree, const std::string &name) {
  return tree && tree->GetBranch(name.c_str());
}

std::vector<std::string> RequiredBranches() {
  std::vector<std::string> names = {
      "H.gtr.beta", "H.gtr.dp", "H.dc.x_fp", "H.cal.etottracknorm",
      "H.cer.npeSum", "P.gtr.dp", "P.gtr.p", "P.cal.etottracknorm",
      "P.aero.npeSum", "P.hgcer.npeSum"};
  for (const std::string plane : {"1x", "1y", "2x", "2y"}) {
    const std::string prefix = "H.hod." + plane + ".";
    names.push_back(prefix + "TrackXPos");
    names.push_back(prefix + "TrackYPos");
    names.push_back(prefix + "fptime");
    for (const std::string side : {"Pos", "Neg"}) {
      names.push_back(prefix + "Good" + side + "AdcHitUsed");
      names.push_back(prefix + "Good" + side + "TdcTimeUnCorr");
      names.push_back(prefix + "Good" + side + "TdcTimeWalkCorr");
      names.push_back(prefix + "Good" + side + "TdcTimeTOFCorr");
      names.push_back(prefix + "Good" + side + "TdcTimeCorr");
      names.push_back(prefix + "Good" + side + "AdcPulseAmp");
      names.push_back(prefix + "Good" + side + "AdcPulseInt");
    }
    names.push_back(prefix + "posTdcCounter");
    names.push_back(prefix + "posTdcTimeRaw");
    names.push_back(prefix + "posAdcCounter");
    names.push_back(prefix + "posAdcPulseAmpRaw");
    names.push_back(prefix + "negTdcCounter");
    names.push_back(prefix + "negTdcTimeRaw");
    names.push_back(prefix + "negAdcCounter");
    names.push_back(prefix + "negAdcPulseAmpRaw");
  }
  return names;
}

std::unique_ptr<TTreeReaderArray<Double_t>> Array(TTreeReader &reader,
                                                  const std::string &name) {
  return std::make_unique<TTreeReaderArray<Double_t>>(reader, name.c_str());
}

SideReaders MakeSide(TTreeReader &reader, const std::string &prefix,
                     const std::string &side) {
  SideReaders output;
  output.used = Array(reader, prefix + "Good" + side + "AdcHitUsed");
  output.uncorr = Array(reader, prefix + "Good" + side + "TdcTimeUnCorr");
  output.walk = Array(reader, prefix + "Good" + side + "TdcTimeWalkCorr");
  output.tof = Array(reader, prefix + "Good" + side + "TdcTimeTOFCorr");
  output.corr = Array(reader, prefix + "Good" + side + "TdcTimeCorr");
  output.adcAmp = Array(reader, prefix + "Good" + side + "AdcPulseAmp");
  output.adcInt = Array(reader, prefix + "Good" + side + "AdcPulseInt");
  const std::string lower = side == "Pos" ? "pos" : "neg";
  output.rawCounter = Array(reader, prefix + lower + "TdcCounter");
  output.rawTime = Array(reader, prefix + lower + "TdcTimeRaw");
  output.rawAdcCounter = Array(reader, prefix + lower + "AdcCounter");
  output.rawAdcAmp = Array(reader, prefix + lower + "AdcPulseAmpRaw");
  return output;
}

PlaneReaders MakePlaneReaders(TTreeReader &reader, const std::string &name,
                              int paddles) {
  PlaneReaders output;
  output.name = name;
  output.paddles = paddles;
  const std::string prefix = "H.hod." + name + ".";
  output.pos = MakeSide(reader, prefix, "Pos");
  output.neg = MakeSide(reader, prefix, "Neg");
  output.trackX = std::make_unique<TTreeReaderValue<Double_t>>(
      reader, (prefix + "TrackXPos").c_str());
  output.trackY = std::make_unique<TTreeReaderValue<Double_t>>(
      reader, (prefix + "TrackYPos").c_str());
  output.fpTime = std::make_unique<TTreeReaderValue<Double_t>>(
      reader, (prefix + "fptime").c_str());
  return output;
}

PlaneOutput MakePlaneOutput(const std::string &name, int paddles,
                            const std::string &label) {
  PlaneOutput output;
  output.name = name;
  output.paddles = paddles;
  output.summary.resize(paddles);
  const std::string tag = Unique("hodo_" + name, label);
  output.paddleVsXfp = std::make_unique<TH2D>(
      (tag + "_paddle_xfp").c_str(),
      ("HMS " + name + " accepted paddle vs xfp;H.dc.x_fp [cm];Paddle").c_str(),
      80, -45, 45, paddles, 0.5, paddles + 0.5);
  output.betaVsPaddle = std::make_unique<TH2D>(
      (tag + "_beta_paddle").c_str(),
      ("HMS beta vs " + name + " accepted paddle;Paddle;H.gtr.beta").c_str(),
      paddles, 0.5, paddles + 0.5, 160, 0.2, 1.2);
  output.fpTimeVsXfp = std::make_unique<TH2D>(
      (tag + "_fptime_xfp").c_str(),
      ("HMS " + name + " fptime vs xfp;H.dc.x_fp [cm];H.hod." + name +
       ".fptime [ns]").c_str(),
      80, -45, 45, 200, 45, 80);
  output.residualVsXfp = std::make_unique<TH2D>(
      (tag + "_residual_xfp").c_str(),
      ("HMS " + name + " time residual vs xfp;H.dc.x_fp [cm];t_{" + name +
       "}-mean(other planes) [ns]").c_str(),
      80, -45, 45, 200, -10, 10);
  for (int paddle = 1; paddle <= paddles; ++paddle) {
    auto bump = std::make_unique<TH1D>(
        (tag + "_beta_bump_p" + std::to_string(paddle)).c_str(),
        (name + " paddle " + std::to_string(paddle) +
         ";H.gtr.beta;Counts").c_str(),
        120, 0.2, 1.2);
    auto control = std::make_unique<TH1D>(
        (tag + "_beta_control_p" + std::to_string(paddle)).c_str(),
        (name + " paddle " + std::to_string(paddle) +
         ";H.gtr.beta;Counts").c_str(),
        120, 0.2, 1.2);
    Detach(bump.get());
    Detach(control.get());
    output.betaBump.push_back(std::move(bump));
    output.betaControl.push_back(std::move(control));
  }
  Detach(output.paddleVsXfp.get());
  Detach(output.betaVsPaddle.get());
  Detach(output.fpTimeVsXfp.get());
  Detach(output.residualVsXfp.get());
  return output;
}

double At(const TTreeReaderArray<Double_t> &values, int index) {
  return index >= 0 && index < static_cast<int>(values.GetSize())
             ? values[index]
             : std::numeric_limits<double>::quiet_NaN();
}

bool Used(const TTreeReaderArray<Double_t> &values, int index) {
  const double value = At(values, index);
  return std::isfinite(value) && std::abs(value) < kInvalid && value > 0.5;
}

double RawForPaddle(const TTreeReaderArray<Double_t> &counter,
                    const TTreeReaderArray<Double_t> &value, int paddle) {
  const auto count = std::min(counter.GetSize(), value.GetSize());
  for (std::size_t index = 0; index < count; ++index) {
    if (std::lround(counter[index]) == paddle)
      return value[index];
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void FillSide(SideSummary &summary, const SideReaders &readers, int index,
              bool bump) {
  auto &stages = bump ? summary.bump : summary.control;
  const std::array<double, 4> corrected = {
      At(*readers.uncorr, index), At(*readers.walk, index),
      At(*readers.tof, index), At(*readers.corr, index)};
  for (int stage = kUncorr; stage <= kCorr; ++stage) {
    const double value = corrected[stage - kUncorr];
    // hcana uses large negative sentinels for a missing/invalid timing side.
    if (std::isfinite(value) && value > -100.0 && value < 200.0)
      stages[stage].Fill(value);
  }
  (bump ? summary.bumpAdcAmp : summary.controlAdcAmp)
      .Fill(At(*readers.adcAmp, index));
  (bump ? summary.bumpAdcInt : summary.controlAdcInt)
      .Fill(At(*readers.adcInt, index));
}

void FillRaw(SideSummary &summary, const SideReaders &readers, int paddle,
             bool bump) {
  auto &stages = bump ? summary.bump : summary.control;
  stages[kRaw].Fill(
      RawForPaddle(*readers.rawCounter, *readers.rawTime, paddle));
  (bump ? summary.bumpRawAdcAmp : summary.controlRawAdcAmp)
      .Fill(RawForPaddle(*readers.rawAdcCounter, *readers.rawAdcAmp, paddle));
}

double Effect(const RunningStats &a, const RunningStats &b) {
  if (a.n < 20 || b.n < 20)
    return 0.0;
  const double sa = a.Sigma();
  const double sb = b.Sigma();
  const double pooled = std::sqrt(0.5 * (sa * sa + sb * sb));
  return pooled > 1.0e-12 ? std::abs(a.Mean() - b.Mean()) / pooled : 0.0;
}

double Shift(const RunningStats &a, const RunningStats &b) {
  return a.n && b.n ? a.Mean() - b.Mean()
                    : std::numeric_limits<double>::quiet_NaN();
}

void Score(PaddleSummary &summary) {
  double best = Effect(summary.betaBump, summary.betaControl);
  double largestTimingEffect = 0.0;
  int firstStage = -1;
  for (int side = 0; side < 2; ++side) {
    const SideSummary &data = side == 0 ? summary.pos : summary.neg;
    for (int stage = 0; stage < kNStages; ++stage) {
      const double effect = Effect(data.bump[stage], data.control[stage]);
      best = std::max(best, effect);
      if (effect > largestTimingEffect) {
        largestTimingEffect = effect;
        firstStage = stage;
      }
    }
  }
  best = std::max(best,
                  Effect(summary.posMinusNegBump, summary.posMinusNegControl));
  summary.score = best;
  summary.candidateStage =
      firstStage >= 0 ? kStageNames[firstStage] : "insufficient_statistics";
}

bool ValidCorrectedSide(const SideReaders &side, int index) {
  const double value = At(*side.corr, index);
  return Used(*side.used, index) && std::isfinite(value) && value > -100.0 &&
         value < 200.0;
}

void AddPlaneAnomalyScore(PlaneOutput &plane) {
  std::vector<double> twoSidedFractions;
  std::vector<double> betaSigmas;
  for (const auto &p : plane.summary) {
    const long long accepted = p.acceptedBump + p.acceptedControl;
    const long long both = p.bothSidesBump + p.bothSidesControl;
    if (accepted >= 500) {
      twoSidedFractions.push_back(static_cast<double>(both) / accepted);
      if (std::isfinite(p.betaBump.Sigma()))
        betaSigmas.push_back(p.betaBump.Sigma());
    }
  }
  auto median = [](std::vector<double> values) {
    if (values.empty())
      return 0.0;
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
  };
  const double medianGood = median(twoSidedFractions);
  const double medianBetaSigma = median(betaSigmas);
  for (auto &p : plane.summary) {
    const long long accepted = p.acceptedBump + p.acceptedControl;
    const long long both = p.bothSidesBump + p.bothSidesControl;
    if (accepted < 500)
      continue;
    const double coverage =
        std::sqrt(std::min(1.0, static_cast<double>(p.acceptedBump) / 5000.0));
    const double goodFraction = static_cast<double>(both) / accepted;
    const double missingPenalty =
        3.0 * std::max(0.0, medianGood - goodFraction) * coverage;
    const double widthPenalty =
        medianBetaSigma > 0.0 && std::isfinite(p.betaBump.Sigma())
            ? std::max(0.0, p.betaBump.Sigma() / medianBetaSigma - 1.0) *
                  coverage
            : 0.0;
    p.score += missingPenalty + widthPenalty;
    if (goodFraction < 0.8 && p.acceptedBump >= 1000) {
      const bool posWorse = p.posValidBump < p.negValidBump;
      const SideSummary &side = posWorse ? p.pos : p.neg;
      const long long used = side.bumpAdcAmp.n;
      const long long rawAdc = side.bumpRawAdcAmp.n;
      const long long rawTdc = side.bump[kRaw].n;
      const long long corrected = side.bump[kCorr].n;
      std::string cause = "missing_or_invalid_";
      cause += posWorse ? "pos" : "neg";
      if (rawAdc < 0.8 * p.acceptedBump)
        cause += "_raw_adc";
      else if (used < 0.8 * rawAdc)
        cause += "_adc_selection";
      if (rawTdc < 0.8 * p.acceptedBump)
        cause += "_raw_tdc";
      if (used > 0 && corrected < 0.8 * used)
        cause += "_timing";
      p.candidateStage = cause;
    }
  }
}

bool Analyze(const char *path, const std::string &label, AnalysisResult &result) {
  std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "[ERROR] Cannot open ROOT file: " << path << '\n';
    return false;
  }
  auto *tree = dynamic_cast<TTree *>(file->Get("T"));
  if (!tree) {
    std::cerr << "[ERROR] Tree 'T' is missing from " << path << '\n';
    return false;
  }
  std::vector<std::string> missing;
  for (const auto &branch : RequiredBranches())
    if (!HasBranch(tree, branch))
      missing.push_back(branch);
  if (!missing.empty()) {
    std::cerr << "[ERROR] Missing " << missing.size() << " required branches:\n";
    for (const auto &branch : missing)
      std::cerr << "  " << branch << '\n';
    return false;
  }

  result.label = label;
  const std::string tag = Unique("hodo_bump", label);
  result.betaVsXfp = std::make_unique<TH2D>(
      (tag + "_beta_xfp").c_str(),
      "Phase 2 run 28425: HMS beta vs HMS xfp;H.dc.x_fp [cm];H.gtr.beta",
      80, -45, 45, 120, 0.2, 1.2);
  result.betaBump = std::make_unique<TH1D>(
      (tag + "_beta_bump").c_str(),
      "Bump and control beta spectra;H.gtr.beta;Normalized counts", 120, 0.2,
      1.2);
  result.betaControl = std::make_unique<TH1D>(
      (tag + "_beta_control").c_str(),
      "Bump and control beta spectra;H.gtr.beta;Normalized counts", 120, 0.2,
      1.2);
  Detach(result.betaVsXfp.get());
  Detach(result.betaBump.get());
  Detach(result.betaControl.get());
  const std::array<std::string, 4> names = {"1x", "1y", "2x", "2y"};
  const std::array<int, 4> paddles = {16, 10, 16, 10};
  for (int plane = 0; plane < 4; ++plane)
    result.planes[plane] = MakePlaneOutput(names[plane], paddles[plane], label);

  TTreeReader reader(tree);
  TTreeReaderValue<Double_t> beta(reader, "H.gtr.beta");
  TTreeReaderValue<Double_t> xfp(reader, "H.dc.x_fp");
  std::array<PlaneReaders, 4> inputs = {
      MakePlaneReaders(reader, "1x", 16), MakePlaneReaders(reader, "1y", 10),
      MakePlaneReaders(reader, "2x", 16), MakePlaneReaders(reader, "2y", 10)};
  TTreeFormula selection(Unique("coin_selection", label).c_str(), kCoinCut,
                         tree);
  result.allEvents = tree->GetEntries();
  long long reportAt = 500000;
  while (reader.Next()) {
    const auto entry = reader.GetCurrentEntry();
    if (entry >= reportAt) {
      std::cout << "[INFO] " << label << ": processed " << entry << "/"
                << result.allEvents << " entries\n";
      reportAt += 500000;
    }
    selection.GetNdata();
    if (selection.EvalInstance() == 0.0)
      continue;
    const double x = *xfp;
    const double b = *beta;
    if (!std::isfinite(x) || !std::isfinite(b))
      continue;
    ++result.selectedEvents;
    result.betaVsXfp->Fill(x, b);
    const bool inBump = x >= kBumpLo && x <= kBumpHi;
    const bool inControl =
        (x >= kControlLo && x < kBumpLo) || (x > kBumpHi && x <= kControlHi);
    if (inBump) {
      ++result.bumpEvents;
      result.betaBump->Fill(b);
    }
    if (inControl) {
      ++result.controlEvents;
      result.betaControl->Fill(b);
    }

    std::array<double, 4> fpTimes;
    for (int plane = 0; plane < 4; ++plane)
      fpTimes[plane] = **inputs[plane].fpTime;
    for (int plane = 0; plane < 4; ++plane) {
      auto &input = inputs[plane];
      auto &output = result.planes[plane];
      const double fp = fpTimes[plane];
      if (std::isfinite(fp) && std::abs(fp) < kInvalid)
        output.fpTimeVsXfp->Fill(x, fp);
      double otherSum = 0.0;
      int otherCount = 0;
      for (int other = 0; other < 4; ++other) {
        if (other == plane || !std::isfinite(fpTimes[other]) ||
            std::abs(fpTimes[other]) >= kInvalid)
          continue;
        otherSum += fpTimes[other];
        ++otherCount;
      }
      if (otherCount && std::isfinite(fp) && std::abs(fp) < kInvalid)
        output.residualVsXfp->Fill(x, fp - otherSum / otherCount);

      for (int index = 0; index < input.paddles; ++index) {
        auto &summary = output.summary[index];
        if (inBump)
          ++summary.selectedBump;
        if (inControl)
          ++summary.selectedControl;
        const bool posUsed = Used(*input.pos.used, index);
        const bool negUsed = Used(*input.neg.used, index);
        if (!posUsed && !negUsed)
          continue;
        const int paddle = index + 1;
        output.paddleVsXfp->Fill(x, paddle);
        output.betaVsPaddle->Fill(paddle, b);
        if (inBump) {
          ++summary.acceptedBump;
          summary.betaBump.Fill(b);
          output.betaBump[index]->Fill(b);
        }
        if (inControl) {
          ++summary.acceptedControl;
          summary.betaControl.Fill(b);
          output.betaControl[index]->Fill(b);
        }
        if (inBump || inControl) {
          FillRaw(summary.pos, input.pos, paddle, inBump);
          FillRaw(summary.neg, input.neg, paddle, inBump);
        }
        if (posUsed && (inBump || inControl))
          FillSide(summary.pos, input.pos, index, inBump);
        if (negUsed && (inBump || inControl))
          FillSide(summary.neg, input.neg, index, inBump);
        const bool posValid = ValidCorrectedSide(input.pos, index);
        const bool negValid = ValidCorrectedSide(input.neg, index);
        if (inBump) {
          summary.posValidBump += posValid;
          summary.negValidBump += negValid;
        }
        if (inControl) {
          summary.posValidControl += posValid;
          summary.negValidControl += negValid;
        }
        if (posValid && negValid) {
          if (inBump)
            ++summary.bothSidesBump;
          if (inControl)
            ++summary.bothSidesControl;
          const double pos = At(*input.pos.corr, index);
          const double neg = At(*input.neg.corr, index);
          if (inBump) {
            summary.posMinusNegBump.Fill(pos - neg);
            summary.posNegMeanBump.Fill(0.5 * (pos + neg));
          }
          if (inControl) {
            summary.posMinusNegControl.Fill(pos - neg);
            summary.posNegMeanControl.Fill(0.5 * (pos + neg));
          }
        }
        const double along =
            input.name.back() == 'x' ? **input.trackY : **input.trackX;
        if (inBump)
          summary.alongBump.Fill(along);
        if (inControl)
          summary.alongControl.Fill(along);
      }
    }
  }
  for (auto &plane : result.planes) {
    for (auto &paddle : plane.summary)
      Score(paddle);
    AddPlaneAnomalyScore(plane);
  }
  std::cout << "[INFO] " << label << ": all=" << result.allEvents
            << " selected=" << result.selectedEvents
            << " bump=" << result.bumpEvents
            << " control=" << result.controlEvents << '\n';
  return true;
}

std::string CsvNumber(double value) {
  if (!std::isfinite(value))
    return "";
  std::ostringstream stream;
  stream << std::setprecision(10) << value;
  return stream.str();
}

bool WriteCsv(const AnalysisResult &result, const std::string &path) {
  std::ofstream output(path);
  if (!output) {
    std::cerr << "[ERROR] Cannot create " << path << '\n';
    return false;
  }
  output << "run,label,plane,paddle,side,selected_bump,selected_control,"
            "accepted_bump,accepted_control,two_sided_bump,two_sided_control,"
            "pos_valid_bump,pos_valid_control,neg_valid_bump,neg_valid_control,"
            "bump_occupancy,control_occupancy,bump_two_sided_fraction,"
            "control_two_sided_fraction,beta_bump_n,beta_bump_mean,"
            "beta_bump_sigma,beta_control_n,beta_control_mean,"
            "beta_control_sigma,beta_shift,pos_minus_neg_shift,"
            "pos_neg_mean_shift,track_along_bump_mean,"
            "track_along_control_mean,raw_bump_mean,raw_control_mean,"
            "uncorrected_bump_mean,uncorrected_control_mean,"
            "time_walk_bump_mean,time_walk_control_mean,tof_bump_mean,"
            "tof_control_mean,corrected_bump_mean,corrected_control_mean,"
            "side_adc_used_bump,side_adc_used_control,raw_adc_bump_n,"
            "raw_adc_control_n,raw_tdc_bump_n,raw_tdc_control_n,"
            "uncorrected_bump_n,uncorrected_control_n,time_walk_bump_n,"
            "time_walk_control_n,tof_bump_n,tof_control_n,corrected_bump_n,"
            "corrected_control_n,raw_adc_amp_bump_mean,"
            "raw_adc_amp_control_mean,adc_amp_bump_mean,adc_amp_control_mean,"
            "adc_int_bump_mean,adc_int_control_mean,score,candidate_stage\n";
  for (const auto &plane : result.planes) {
    for (int index = 0; index < plane.paddles; ++index) {
      const auto &p = plane.summary[index];
      for (int sideIndex = 0; sideIndex < 2; ++sideIndex) {
        const auto &side = sideIndex == 0 ? p.pos : p.neg;
        output << kRun << ',' << result.label << ',' << plane.name << ','
               << index + 1 << ',' << (sideIndex == 0 ? "pos" : "neg") << ','
               << p.selectedBump << ',' << p.selectedControl << ','
               << p.acceptedBump << ',' << p.acceptedControl << ','
               << p.bothSidesBump << ',' << p.bothSidesControl << ','
               << p.posValidBump << ',' << p.posValidControl << ','
               << p.negValidBump << ',' << p.negValidControl << ','
               << CsvNumber(p.selectedBump
                                ? static_cast<double>(p.acceptedBump) /
                                      p.selectedBump
                                : NAN)
               << ','
               << CsvNumber(p.selectedControl
                                ? static_cast<double>(p.acceptedControl) /
                                      p.selectedControl
                                : NAN)
               << ','
               << CsvNumber(p.acceptedBump
                                ? static_cast<double>(p.bothSidesBump) /
                                      p.acceptedBump
                                : NAN)
               << ','
               << CsvNumber(p.acceptedControl
                                ? static_cast<double>(p.bothSidesControl) /
                                      p.acceptedControl
                                : NAN)
               << ',' << p.betaBump.n << ',' << CsvNumber(p.betaBump.Mean())
               << ',' << CsvNumber(p.betaBump.Sigma()) << ',' << p.betaControl.n
               << ',' << CsvNumber(p.betaControl.Mean()) << ','
               << CsvNumber(p.betaControl.Sigma()) << ','
               << CsvNumber(Shift(p.betaBump, p.betaControl)) << ','
               << CsvNumber(Shift(p.posMinusNegBump, p.posMinusNegControl))
               << ',' << CsvNumber(Shift(p.posNegMeanBump, p.posNegMeanControl))
               << ',' << CsvNumber(p.alongBump.Mean()) << ','
               << CsvNumber(p.alongControl.Mean());
        for (int stage = 0; stage < kNStages; ++stage)
          output << ',' << CsvNumber(side.bump[stage].Mean()) << ','
                 << CsvNumber(side.control[stage].Mean());
        output << ',' << side.bumpAdcAmp.n << ',' << side.controlAdcAmp.n
               << ',' << side.bumpRawAdcAmp.n << ','
               << side.controlRawAdcAmp.n;
        for (int stage = 0; stage < kNStages; ++stage)
          output << ',' << side.bump[stage].n << ',' << side.control[stage].n;
        output << ',' << CsvNumber(side.bumpRawAdcAmp.Mean()) << ','
               << CsvNumber(side.controlRawAdcAmp.Mean()) << ','
               << CsvNumber(side.bumpAdcAmp.Mean()) << ','
               << CsvNumber(side.controlAdcAmp.Mean()) << ','
               << CsvNumber(side.bumpAdcInt.Mean()) << ','
               << CsvNumber(side.controlAdcInt.Mean()) << ','
               << CsvNumber(p.score) << ',' << p.candidateStage << '\n';
      }
    }
  }
  return output.good();
}

void DrawRegionLines(double ymin, double ymax) {
  for (double x : {kBumpLo, kBumpHi}) {
    auto *line = new TLine(x, ymin, x, ymax);
    line->SetLineColor(kRed + 1);
    line->SetLineWidth(2);
    line->Draw();
  }
}

void PrintPage(TCanvas &canvas, const std::string &pdf) {
  canvas.Print(pdf.c_str());
  canvas.Clear();
}

struct Candidate {
  std::string plane;
  int paddle = 0;
  double score = 0.0;
  std::string stage;
};

std::vector<Candidate> Candidates(const AnalysisResult &result);

bool WritePdf(const AnalysisResult &result, const std::string &path) {
  TCanvas canvas("hodo_bump_canvas", "Hodoscope bump diagnostic", 1500, 1000);
  canvas.Print((path + "[").c_str());
  canvas.Divide(2, 1);
  canvas.cd(1);
  result.betaVsXfp->Draw("COLZ");
  DrawRegionLines(0.2, 1.2);
  canvas.cd(2);
  auto bump = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(result.betaBump->Clone("beta_bump_pdf")));
  auto control = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(result.betaControl->Clone("beta_control_pdf")));
  if (bump->Integral())
    bump->Scale(1.0 / bump->Integral());
  if (control->Integral())
    control->Scale(1.0 / control->Integral());
  bump->SetLineColor(kRed + 1);
  bump->SetLineWidth(2);
  control->SetLineColor(kBlue + 1);
  control->SetLineWidth(2);
  bump->SetMaximum(1.2 * std::max(bump->GetMaximum(), control->GetMaximum()));
  bump->Draw("HIST");
  control->Draw("HIST SAME");
  TLegend legend(0.58, 0.75, 0.88, 0.88);
  legend.AddEntry(bump.get(), "xfp 12-18 cm", "l");
  legend.AddEntry(control.get(), "xfp 6-12, 18-24 cm", "l");
  legend.Draw();
  PrintPage(canvas, path);

  for (const auto &plane : result.planes) {
    canvas.Divide(2, 2);
    canvas.cd(1);
    plane.paddleVsXfp->Draw("COLZ");
    DrawRegionLines(0.5, plane.paddles + 0.5);
    canvas.cd(2);
    plane.betaVsPaddle->Draw("COLZ");
    canvas.cd(3);
    plane.fpTimeVsXfp->Draw("COLZ");
    DrawRegionLines(45, 80);
    canvas.cd(4);
    plane.residualVsXfp->Draw("COLZ");
    DrawRegionLines(-10, 10);
    PrintPage(canvas, path);

    const int columns = plane.paddles == 16 ? 4 : 5;
    const int rows = plane.paddles == 16 ? 4 : 2;
    canvas.Divide(columns, rows);
    for (int index = 0; index < plane.paddles; ++index) {
      canvas.cd(index + 1);
      auto *b = plane.betaBump[index].get();
      auto *c = plane.betaControl[index].get();
      b->SetLineColor(kRed + 1);
      c->SetLineColor(kBlue + 1);
      b->SetLineWidth(2);
      c->SetLineWidth(2);
      const double bIntegral = b->Integral();
      const double cIntegral = c->Integral();
      if (bIntegral)
        b->Scale(1.0 / bIntegral);
      if (cIntegral)
        c->Scale(1.0 / cIntegral);
      b->SetMaximum(1.2 * std::max(b->GetMaximum(), c->GetMaximum()));
      b->GetYaxis()->SetTitle("Normalized counts");
      b->Draw("HIST");
      c->Draw("HIST SAME");
    }
    PrintPage(canvas, path);

    std::vector<double> paddle(plane.paddles), score(plane.paddles),
        bumpOcc(plane.paddles), controlOcc(plane.paddles), goodBump(plane.paddles),
        goodControl(plane.paddles);
    for (int index = 0; index < plane.paddles; ++index) {
      const auto &p = plane.summary[index];
      paddle[index] = index + 1;
      score[index] = p.score;
      bumpOcc[index] = p.selectedBump
                           ? static_cast<double>(p.acceptedBump) / p.selectedBump
                           : 0.0;
      controlOcc[index] =
          p.selectedControl
              ? static_cast<double>(p.acceptedControl) / p.selectedControl
              : 0.0;
      goodBump[index] = p.acceptedBump
                            ? static_cast<double>(p.bothSidesBump) / p.acceptedBump
                            : 0.0;
      goodControl[index] =
          p.acceptedControl
              ? static_cast<double>(p.bothSidesControl) / p.acceptedControl
              : 0.0;
    }
    canvas.Divide(2, 2);
    canvas.cd(1);
    TGraph bumpGraph(plane.paddles, paddle.data(), bumpOcc.data());
    TGraph controlGraph(plane.paddles, paddle.data(), controlOcc.data());
    bumpGraph.SetTitle(("HMS " + plane.name +
                        " occupancy;Paddle;Accepted / selected events").c_str());
    bumpGraph.SetMarkerStyle(20);
    bumpGraph.SetMarkerColor(kRed + 1);
    bumpGraph.SetLineColor(kRed + 1);
    controlGraph.SetMarkerStyle(22);
    controlGraph.SetMarkerColor(kBlue + 1);
    controlGraph.SetLineColor(kBlue + 1);
    bumpGraph.Draw("APL");
    controlGraph.Draw("PL SAME");
    TLegend occLegend(0.58, 0.75, 0.88, 0.88);
    occLegend.AddEntry(&bumpGraph, "Bump region", "lp");
    occLegend.AddEntry(&controlGraph, "Control region", "lp");
    occLegend.Draw();
    canvas.cd(2);
    TGraph goodBumpGraph(plane.paddles, paddle.data(), goodBump.data());
    TGraph goodControlGraph(plane.paddles, paddle.data(), goodControl.data());
    goodBumpGraph.SetTitle(("HMS " + plane.name +
                            " two-sided good-hit fraction;Paddle;Both PMTs / accepted events")
                               .c_str());
    goodBumpGraph.SetMarkerStyle(20);
    goodBumpGraph.SetMarkerColor(kRed + 1);
    goodBumpGraph.SetLineColor(kRed + 1);
    goodControlGraph.SetMarkerStyle(22);
    goodControlGraph.SetMarkerColor(kBlue + 1);
    goodControlGraph.SetLineColor(kBlue + 1);
    goodBumpGraph.Draw("APL");
    goodControlGraph.Draw("PL SAME");
    TLegend goodLegend(0.58, 0.75, 0.88, 0.88);
    goodLegend.AddEntry(&goodBumpGraph, "Bump region", "lp");
    goodLegend.AddEntry(&goodControlGraph, "Control region", "lp");
    goodLegend.Draw();
    canvas.cd(3);
    TGraph scoreGraph(plane.paddles, paddle.data(), score.data());
    scoreGraph.SetTitle(("HMS " + plane.name +
                         " candidate ranking;Paddle;Maximum standardized effect")
                            .c_str());
    scoreGraph.SetMarkerStyle(20);
    scoreGraph.Draw("APL");
    canvas.cd(4);
    std::vector<double> betaShift(plane.paddles);
    for (int index = 0; index < plane.paddles; ++index)
      betaShift[index] = Shift(plane.summary[index].betaBump,
                               plane.summary[index].betaControl);
    TGraph betaShiftGraph(plane.paddles, paddle.data(), betaShift.data());
    betaShiftGraph.SetTitle(("HMS " + plane.name +
                             " beta displacement;Paddle;Bump mean - control mean")
                                .c_str());
    betaShiftGraph.SetMarkerStyle(20);
    betaShiftGraph.Draw("APL");
    PrintPage(canvas, path);
  }

  const auto candidates = Candidates(result);
  const std::size_t candidateCount = std::min<std::size_t>(10, candidates.size());
  for (std::size_t rank = 0; rank < candidateCount; ++rank) {
    const auto &candidate = candidates[rank];
    const PlaneOutput *plane = nullptr;
    for (const auto &entry : result.planes)
      if (entry.name == candidate.plane)
        plane = &entry;
    if (!plane)
      continue;
    const auto &p = plane->summary[candidate.paddle - 1];
    std::array<double, kNStages> stageIndex{}, posBump{}, posControl{},
        negBump{}, negControl{}, posShift{}, negShift{};
    for (int stage = 0; stage < kNStages; ++stage) {
      stageIndex[stage] = stage + 1;
      posBump[stage] = p.pos.bump[stage].Mean();
      posControl[stage] = p.pos.control[stage].Mean();
      negBump[stage] = p.neg.bump[stage].Mean();
      negControl[stage] = p.neg.control[stage].Mean();
      posShift[stage] = Shift(p.pos.bump[stage], p.pos.control[stage]);
      negShift[stage] = Shift(p.neg.bump[stage], p.neg.control[stage]);
    }
    canvas.Divide(2, 2);
    canvas.cd(1);
    TGraph posBumpGraph(kNStages, stageIndex.data(), posBump.data());
    TGraph posControlGraph(kNStages, stageIndex.data(), posControl.data());
    posBumpGraph.SetTitle(("Rank " + std::to_string(rank + 1) + ": " +
                           candidate.plane + " paddle " +
                           std::to_string(candidate.paddle) +
                           " positive PMT;Correction stage;Mean time")
                              .c_str());
    posBumpGraph.SetMarkerStyle(20);
    posBumpGraph.SetMarkerColor(kRed + 1);
    posBumpGraph.SetLineColor(kRed + 1);
    posControlGraph.SetMarkerStyle(22);
    posControlGraph.SetMarkerColor(kBlue + 1);
    posControlGraph.SetLineColor(kBlue + 1);
    posBumpGraph.Draw("APL");
    posControlGraph.Draw("PL SAME");
    canvas.cd(2);
    TGraph negBumpGraph(kNStages, stageIndex.data(), negBump.data());
    TGraph negControlGraph(kNStages, stageIndex.data(), negControl.data());
    negBumpGraph.SetTitle((candidate.plane + " paddle " +
                           std::to_string(candidate.paddle) +
                           " negative PMT;Correction stage;Mean time")
                              .c_str());
    negBumpGraph.SetMarkerStyle(20);
    negBumpGraph.SetMarkerColor(kRed + 1);
    negBumpGraph.SetLineColor(kRed + 1);
    negControlGraph.SetMarkerStyle(22);
    negControlGraph.SetMarkerColor(kBlue + 1);
    negControlGraph.SetLineColor(kBlue + 1);
    negBumpGraph.Draw("APL");
    negControlGraph.Draw("PL SAME");
    canvas.cd(3);
    TGraph posShiftGraph(kNStages, stageIndex.data(), posShift.data());
    TGraph negShiftGraph(kNStages, stageIndex.data(), negShift.data());
    posShiftGraph.SetTitle("Bump-control timing shift;Correction stage;Mean shift");
    posShiftGraph.SetMarkerStyle(20);
    posShiftGraph.SetMarkerColor(kRed + 1);
    posShiftGraph.SetLineColor(kRed + 1);
    negShiftGraph.SetMarkerStyle(22);
    negShiftGraph.SetMarkerColor(kBlue + 1);
    negShiftGraph.SetLineColor(kBlue + 1);
    posShiftGraph.Draw("APL");
    negShiftGraph.Draw("PL SAME");
    TLegend sideLegend(0.58, 0.75, 0.88, 0.88);
    sideLegend.AddEntry(&posShiftGraph, "Positive PMT", "lp");
    sideLegend.AddEntry(&negShiftGraph, "Negative PMT", "lp");
    sideLegend.Draw();
    canvas.cd(4);
    const double categories[] = {1, 2, 3, 4};
    const double values[] = {p.pos.bumpAdcAmp.Mean(), p.pos.controlAdcAmp.Mean(),
                             p.neg.bumpAdcAmp.Mean(), p.neg.controlAdcAmp.Mean()};
    TGraph adcGraph(4, categories, values);
    adcGraph.SetTitle("ADC amplitude means;1=pos bump, 2=pos control, 3=neg bump, 4=neg control;Mean pulse amplitude");
    adcGraph.SetMarkerStyle(20);
    adcGraph.Draw("AP");
    PrintPage(canvas, path);
  }
  canvas.Print((path + "]").c_str());
  return !gSystem->AccessPathName(path.c_str());
}

std::vector<Candidate> Candidates(const AnalysisResult &result) {
  std::vector<Candidate> candidates;
  for (const auto &plane : result.planes)
    for (int index = 0; index < plane.paddles; ++index)
      candidates.push_back(
          {plane.name, index + 1, plane.summary[index].score,
           plane.summary[index].candidateStage});
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &a, const Candidate &b) {
              return a.score > b.score;
            });
  return candidates;
}

bool WriteConclusion(const AnalysisResult &result, const std::string &path) {
  std::ofstream output(path);
  if (!output)
    return false;
  output << "Run " << kRun << " HMS beta-vs-xfp bump diagnostic\n\n"
         << "Selection: " << kCoinCut << "\n"
         << "All events: " << result.allEvents << "\n"
         << "Selected events: " << result.selectedEvents << "\n"
         << "Bump-region events (12 <= xfp <= 18 cm): " << result.bumpEvents
         << "\nControl-region events: " << result.controlEvents << "\n\n"
         << "Highest-ranked plane/paddle candidates:\n";
  const auto candidates = Candidates(result);
  const std::size_t count = std::min<std::size_t>(10, candidates.size());
  for (std::size_t index = 0; index < count; ++index)
  {
    const auto &candidate = candidates[index];
    const PlaneOutput *plane = nullptr;
    for (const auto &entry : result.planes)
      if (entry.name == candidate.plane)
        plane = &entry;
    const auto &p = plane->summary[candidate.paddle - 1];
    output << index + 1 << ". " << candidate.plane << " paddle "
           << candidate.paddle << ", score=" << candidate.score
           << ", candidate stage=" << candidate.stage
           << ", bump accepted=" << p.acceptedBump
           << ", pos valid=" << p.posValidBump
           << ", neg valid=" << p.negValidBump << '\n';
  }
  if (!candidates.empty() &&
      candidates.front().stage.find("raw_adc_raw_tdc") != std::string::npos) {
    output << "\nConclusion: the leading candidate loses both raw ADC and raw "
              "TDC information on one PMT side. This is a raw-level "
              "channel/readout problem, not a time-walk, propagation-velocity, "
              "cable-offset, or LCoeff calibration effect. Do not change those "
              "calibration parameters; inspect the PMT/HV, signal path, "
              "discriminator/TDC channel, and run-period hardware records first.\n";
  }
  output << "\nThe ranking is evidence triage, not an automatic calibration "
            "decision. Inspect the PDF and CSV correction-stage shifts before "
            "changing exactly one parameter.\n";
  return output.good();
}

bool Publish(const std::string &temporary, const std::string &finalPath) {
  if (gSystem->Rename(temporary.c_str(), finalPath.c_str()) != 0) {
    std::cerr << "[ERROR] Cannot publish " << finalPath << '\n';
    return false;
  }
  return true;
}

bool WriteComparison(const AnalysisResult &baseline, const AnalysisResult &trial,
                     const std::string &pdfPath, const std::string &csvPath) {
  const std::string pdfTmp = pdfPath + ".tmp.pdf";
  TCanvas canvas("hodo_compare_canvas", "Before/after comparison", 1500, 1000);
  canvas.Print((pdfTmp + "[").c_str());
  canvas.Divide(2, 1);
  canvas.cd(1);
  baseline.betaVsXfp->SetTitle("Baseline;H.dc.x_fp [cm];H.gtr.beta");
  baseline.betaVsXfp->Draw("COLZ");
  DrawRegionLines(0.2, 1.2);
  canvas.cd(2);
  trial.betaVsXfp->SetTitle("Trial;H.dc.x_fp [cm];H.gtr.beta");
  trial.betaVsXfp->Draw("COLZ");
  DrawRegionLines(0.2, 1.2);
  PrintPage(canvas, pdfTmp);
  for (int plane = 0; plane < 4; ++plane) {
    canvas.Divide(2, 1);
    canvas.cd(1);
    baseline.planes[plane].residualVsXfp->SetTitle(
        ("Baseline " + baseline.planes[plane].name +
         ";H.dc.x_fp [cm];Plane residual [ns]").c_str());
    baseline.planes[plane].residualVsXfp->Draw("COLZ");
    canvas.cd(2);
    trial.planes[plane].residualVsXfp->SetTitle(
        ("Trial " + trial.planes[plane].name +
         ";H.dc.x_fp [cm];Plane residual [ns]").c_str());
    trial.planes[plane].residualVsXfp->Draw("COLZ");
    PrintPage(canvas, pdfTmp);
  }
  canvas.Print((pdfTmp + "]").c_str());

  const std::string csvTmp = csvPath + ".tmp";
  std::ofstream csv(csvTmp);
  csv << "run,plane,paddle,baseline_score,trial_score,score_change,"
         "baseline_beta_shift,trial_beta_shift,beta_shift_change,"
         "baseline_stage,trial_stage\n";
  for (int plane = 0; plane < 4; ++plane) {
    const auto &before = baseline.planes[plane];
    const auto &after = trial.planes[plane];
    for (int index = 0; index < before.paddles; ++index) {
      const auto &b = before.summary[index];
      const auto &a = after.summary[index];
      const double bShift = Shift(b.betaBump, b.betaControl);
      const double aShift = Shift(a.betaBump, a.betaControl);
      csv << kRun << ',' << before.name << ',' << index + 1 << ',' << b.score
          << ',' << a.score << ',' << a.score - b.score << ','
          << CsvNumber(bShift) << ',' << CsvNumber(aShift) << ','
          << CsvNumber(aShift - bShift) << ',' << b.candidateStage << ','
          << a.candidateStage << '\n';
    }
  }
  csv.close();
  return Publish(pdfTmp, pdfPath) && Publish(csvTmp, csvPath);
}

} // namespace hodo_bump

void hodo_bump_diagnostic(
    const char *BaselineFile =
        "/Volumes/T7/RSIDIS/SampleFiles/RootFiles/"
        "coin_replay_production_28425_-1.root",
    const char *TrialFile = "", int RunNumber = 28425,
    const char *OutputDir = "results/Diagnostics/run28425") {
  using namespace hodo_bump;
  if (RunNumber != kRun) {
    std::cerr << "[ERROR] This diagnostic is fixed to run " << kRun
              << "; received " << RunNumber << ".\n";
    return;
  }
  if (!BaselineFile || !*BaselineFile || !OutputDir || !*OutputDir) {
    std::cerr << "[ERROR] BaselineFile and OutputDir must be non-empty.\n";
    return;
  }
  if (gSystem->mkdir(OutputDir, true) != 0 &&
      gSystem->AccessPathName(OutputDir)) {
    std::cerr << "[ERROR] Cannot create output directory " << OutputDir << '\n';
    return;
  }

  AnalysisResult baseline;
  if (!Analyze(BaselineFile, "baseline", baseline))
    return;
  const std::string base(OutputDir);
  const std::string pdf = base + "/hodo_bump_diagnostic_run28425.pdf";
  const std::string csv = base + "/hodo_bump_paddle_summary_run28425.csv";
  const std::string conclusion = base + "/hodo_bump_conclusion_run28425.txt";
  const std::string pdfTmp = pdf + ".tmp.pdf";
  const std::string csvTmp = csv + ".tmp";
  const std::string conclusionTmp = conclusion + ".tmp";
  if (!WritePdf(baseline, pdfTmp) || !WriteCsv(baseline, csvTmp) ||
      !WriteConclusion(baseline, conclusionTmp)) {
    std::cerr << "[ERROR] Failed to create baseline diagnostic outputs.\n";
    return;
  }
  if (!Publish(pdfTmp, pdf) || !Publish(csvTmp, csv) ||
      !Publish(conclusionTmp, conclusion))
    return;
  std::cout << "[INFO] Published baseline diagnostic outputs under " << base
            << '\n';

  if (TrialFile && *TrialFile) {
    AnalysisResult trial;
    if (!Analyze(TrialFile, "trial", trial))
      return;
    const std::string trialCsv =
        base + "/hodo_bump_paddle_summary_run28425_trial.csv";
    const std::string trialTmp = trialCsv + ".tmp";
    if (!WriteCsv(trial, trialTmp) || !Publish(trialTmp, trialCsv))
      return;
    if (!WriteComparison(
            baseline, trial,
            base + "/hodo_bump_before_after_run28425.pdf",
            base + "/hodo_bump_before_after_run28425.csv"))
      return;
    std::cout << "[INFO] Published baseline/trial comparison outputs.\n";
  }
}
