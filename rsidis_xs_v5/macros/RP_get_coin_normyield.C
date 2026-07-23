// RP_get_coin_normyield.C
//
// Calculate run-by-run normalized coincidence yields from an
// RP_get_good_coin_ev per-setting CSV and the corresponding phase bigtable.
// Run from rsidis_xs_v5 with:
//   root -l -b -q 'macros/RP_get_coin_normyield.C("results/Tables/RP_get_good_coin_ev/<setting>.csv")'


#include <algorithm>
#include <cctype>
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

#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kTriggerEff = 1.0;
constexpr double kPIDEff = 1.0;
constexpr double kProvisionalCompLivetime = 1.0;

using CSVRow = std::map<std::string, std::string>;

struct CSVTable {
  std::vector<std::string> header;
  std::vector<CSVRow> rows;
};

struct NormRow {
  CSVRow good;
  CSVRow big;
  int run = 0;
  bool isPositron = false;
  std::string psTrigger = "INVALID";
  double psFactor = kNaN;
  double charge = kNaN;
  double bigtableCompLivetime = kNaN;
  double compLivetime = kProvisionalCompLivetime;
  double hmsEff = kNaN;
  double shmsEff = kNaN;
  double boilCorr = kNaN;
  double goodcoin = kNaN;
  double goodcoinErr = kNaN;
  double replayNormyield = kNaN;
  double replayNormyieldErr = kNaN;
  double normFactor = kNaN;
  double rpNormyield = kNaN;
  double rpNormyieldErr = kNaN;
  double replayByRP = kNaN;
  std::string status = "NOT_PROCESSED";
  std::string reason;
};

std::string Trim(const std::string& value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
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
  fields.push_back(Trim(field));
  return fields;
}

std::string CSVQuote(const std::string& value) {
  if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
  std::string escaped;
  for (char c : value) {
    escaped += c;
    if (c == '"') escaped += '"';
  }
  return '"' + escaped + '"';
}

double ParseDouble(const std::string& value) {
  const std::string text = Trim(value);
  if (text.empty()) return kNaN;
  char* end = nullptr;
  const double number = std::strtod(text.c_str(), &end);
  return end != text.c_str() && *end == '\0' ? number : kNaN;
}

int ParseRun(const std::string& value) {
  const double number = ParseDouble(value);
  return std::isfinite(number) ? static_cast<int>(std::lround(number)) : 0;
}

std::string Get(const CSVRow& row, const std::string& name) {
  const auto exact = row.find(name);
  if (exact != row.end()) return exact->second;
  const std::string lowerName = ToLower(name);
  for (const auto& item : row)
    if (ToLower(item.first) == lowerName) return item.second;
  return "";
}

bool ReadCSV(const std::string& path, CSVTable& table) {
  std::ifstream input(path);
  if (!input) {
    std::cerr << "ERROR: cannot open CSV: " << path << "\n";
    return false;
  }
  std::string line;
  if (!std::getline(input, line)) {
    std::cerr << "ERROR: empty CSV: " << path << "\n";
    return false;
  }
  table.header = SplitCSV(line);
  while (std::getline(input, line)) {
    if (Trim(line).empty()) continue;
    const auto fields = SplitCSV(line);
    CSVRow row;
    for (size_t i = 0; i < table.header.size(); ++i)
      // Preserve exact spelling: the Goodcoin schema contains both CTmean
      // (this analysis) and ctmean (legacy replay), which are distinct.
      row[table.header[i]] = i < fields.size() ? fields[i] : "";
    table.rows.push_back(std::move(row));
  }
  return true;
}

bool HasColumns(const CSVTable& table,
                const std::vector<std::string>& required,
                const std::string& label) {
  std::set<std::string> available;
  for (const auto& name : table.header) available.insert(ToLower(name));
  bool ok = true;
  for (const auto& name : required) {
    if (!available.count(ToLower(name))) {
      std::cerr << "ERROR: " << label << " is missing column " << name << "\n";
      ok = false;
    }
  }
  return ok;
}

std::string InputStem(const std::string& path) {
  std::string name = gSystem->BaseName(path.c_str());
  const auto dot = name.find_last_of('.');
  return dot == std::string::npos ? name : name.substr(0, dot);
}

std::string PhaseFromPathAndRows(const std::string& path,
                                 const CSVTable& table) {
  const std::string lower = ToLower(path);
  if (lower.find("phase1") != std::string::npos) return "phase1";
  if (lower.find("phase2") != std::string::npos) return "phase2";
  if (!table.rows.empty()) {
    const double run = ParseDouble(Get(table.rows.front(), "run"));
    if (std::isfinite(run)) {
      // Phase-1 production runs precede the phase-2 run block.
      if (run < 30000) return "phase1";
      return "phase2";
    }
  }
  return "UNKNOWN";
}

std::string LeafBigtablePath(const std::string& inputPath) {
  const std::string stem = InputStem(inputPath);
  // SplitCSV is intentionally not used for underscores; retain setting names
  // verbatim while extracting the five categorical path components.
  std::vector<std::string> tokens;
  std::stringstream stream(stem);
  std::string token;
  while (std::getline(stream, token, '_')) tokens.push_back(token);
  if (tokens.size() != 6) return "";
  const std::string& phase = tokens[0];
  const std::string& passName = tokens[1];
  const std::string& runType = tokens[2];
  const std::string& target = tokens[3];
  const std::string& charge = tokens[4];
  const std::string& setting = tokens[5];
  return "bigtable/" + phase + "/" + passName + "/" + runType + "/" +
         target + "/" + charge + "/" + setting + "/" + stem + ".csv";
}

void AddReason(NormRow& row, const std::string& reason) {
  if (!row.reason.empty()) row.reason += ';';
  row.reason += reason;
}

bool ValidMeasured(double value) {
  return std::isfinite(value) && value != -999.0;
}

NormRow BuildNormRow(const CSVRow& good, const std::vector<const CSVRow*>& matches) {
  NormRow out;
  out.good = good;
  out.run = ParseRun(Get(good, "Run"));
  out.goodcoin = ParseDouble(Get(good, "RP_Goodcoin"));
  out.goodcoinErr = ParseDouble(Get(good, "RP_Goodcoin_err"));
  out.isPositron = ParseDouble(Get(good, "hms_p")) > 0.0;

  if (matches.size() != 1) {
    out.status = "FAILED";
    AddReason(out, matches.empty() ? "BIGTABLE_RUN_NOT_FOUND"
                                   : "BIGTABLE_RUN_DUPLICATED");
    return out;
  }
  out.big = *matches.front();

  const std::string goodTarget = Get(good, "Target");
  const std::string goodRunType = Get(good, "Run_type");
  const bool categoryMatches =
    goodTarget == Get(out.big, "target") &&
    goodRunType == Get(out.big, "run_type") &&
    ((ParseDouble(Get(good, "hms_p")) < 0.0) ==
     (ParseDouble(Get(out.big, "hms_p")) < 0.0));
  bool settingMatches = true;
  for (const char* field : {"ebeam", "x", "q2", "z", "thpq"}) {
    const double goodValue = ParseDouble(Get(good, field));
    const double leafValue = ParseDouble(Get(out.big, field));
    if (!std::isfinite(goodValue) || !std::isfinite(leafValue) ||
        std::abs(goodValue - leafValue) > 5e-4) {
      settingMatches = false;
      break;
    }
  }
  if (!categoryMatches) AddReason(out, "LEAF_CATEGORY_MISMATCH");
  if (!settingMatches) AddReason(out, "LEAF_SETTING_MISMATCH");

  out.charge = ParseDouble(Get(out.big, "BCM2_Q"));
  out.bigtableCompLivetime = ParseDouble(Get(out.big, "comp_livetime"));
  out.hmsEff = ParseDouble(Get(out.big, "h_esing_Eff"));
  out.shmsEff = ParseDouble(Get(out.big, "p_hadron_Eff"));
  out.boilCorr = ParseDouble(Get(out.big, "boil_corr"));
  out.replayNormyield = ParseDouble(Get(out.big, "normyield"));
  out.replayNormyieldErr = ParseDouble(Get(out.big, "normyield_err"));

  const double ps5 = ParseDouble(Get(out.big, "ps5"));
  const double ps6 = ParseDouble(Get(out.big, "ps6"));
  const bool ps5Positive = ValidMeasured(ps5) && ps5 > 0.0;
  const bool ps6Positive = ValidMeasured(ps6) && ps6 > 0.0;
  if (ps5Positive != ps6Positive) {
    out.psTrigger = ps5Positive ? "ps5" : "ps6";
    out.psFactor = ps5Positive ? ps5 : ps6;
  } else {
    AddReason(out, "ACTIVE_PRESCALE_NOT_UNIQUE");
  }
  if (std::isfinite(out.psFactor) && std::abs(out.psFactor - 1.0) > 1e-9)
    AddReason(out, "ACTIVE_PRESCALE_NOT_ONE");
  if (!(ValidMeasured(out.charge) && out.charge > 0.0))
    AddReason(out, "INVALID_BCM2_Q");
  if (!(ValidMeasured(out.hmsEff) && out.hmsEff > 0.0 && out.hmsEff <= 1.0))
    AddReason(out, "INVALID_HMS_TRACKING_EFF");
  if (!(ValidMeasured(out.shmsEff) && out.shmsEff > 0.0 && out.shmsEff <= 1.0))
    AddReason(out, "INVALID_SHMS_TRACKING_EFF");
  if (!(ValidMeasured(out.boilCorr) && out.boilCorr > 0.0))
    AddReason(out, "INVALID_BOIL_CORR");
  if (!(ValidMeasured(out.goodcoin))) AddReason(out, "INVALID_RP_GOODCOIN");
  if (!(ValidMeasured(out.goodcoinErr) && out.goodcoinErr >= 0.0))
    AddReason(out, "INVALID_RP_GOODCOIN_ERR");

  if (!out.reason.empty()) {
    out.status = "FAILED";
    return out;
  }

  out.normFactor = out.psFactor * out.boilCorr /
    (out.charge * out.compLivetime * out.hmsEff * out.shmsEff *
     kTriggerEff * kPIDEff);
  out.rpNormyield = out.goodcoin * out.normFactor;
  out.rpNormyieldErr = out.goodcoinErr * out.normFactor;
  if (ValidMeasured(out.replayNormyield) && out.rpNormyield != 0.0)
    out.replayByRP = out.replayNormyield / out.rpNormyield;
  out.status = "OK";
  out.reason = "PROVISIONAL_COMP_LIVETIME_1P0";
  return out;
}

void WriteNumber(std::ostream& output, double value) {
  if (std::isfinite(value)) output << std::fixed << std::setprecision(4) << value;
  else output << "nan";
}

void WriteOutputCSV(const std::string& path, const CSVTable& input,
                    const std::vector<NormRow>& rows) {
  std::ofstream output(path);
  for (size_t i = 0; i < input.header.size(); ++i) {
    if (i) output << ',';
    output << CSVQuote(input.header[i]);
  }
  output << ",BCM2_Q,ps5,ps6,PS_trigger,PS_factor,bigtable_comp_livetime,"
            "comp_livetime,comp_livetime_source,h_esing_Eff,p_hadron_Eff,"
            "Trigger_eff,PID_eff,boil_corr,Norm_factor,RP_Normyield,"
            "RP_Normyield_err,normyield,normyield_err,"
            "normyield_by_RP_Normyield,Norm_status,Norm_reason\n";

  for (const auto& row : rows) {
    for (size_t i = 0; i < input.header.size(); ++i) {
      if (i) output << ',';
      output << CSVQuote(Get(row.good, input.header[i]));
    }
    output << ','; WriteNumber(output, row.charge);
    output << ','; WriteNumber(output, ParseDouble(Get(row.big, "ps5")));
    output << ','; WriteNumber(output, ParseDouble(Get(row.big, "ps6")));
    output << ',' << row.psTrigger << ','; WriteNumber(output, row.psFactor);
    output << ','; WriteNumber(output, row.bigtableCompLivetime);
    output << ','; WriteNumber(output, row.compLivetime);
    output << ",PROVISIONAL_CONSTANT,";
    WriteNumber(output, row.hmsEff); output << ',';
    WriteNumber(output, row.shmsEff); output << ',';
    WriteNumber(output, kTriggerEff); output << ',';
    WriteNumber(output, kPIDEff); output << ',';
    WriteNumber(output, row.boilCorr); output << ',';
    WriteNumber(output, row.normFactor); output << ',';
    WriteNumber(output, row.rpNormyield); output << ',';
    WriteNumber(output, row.rpNormyieldErr); output << ',';
    WriteNumber(output, row.replayNormyield); output << ',';
    WriteNumber(output, row.replayNormyieldErr); output << ',';
    WriteNumber(output, row.replayByRP);
    output << ',' << row.status << ',' << row.reason << '\n';
  }
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

struct AxisRange {
  double low = 0.0;
  double high = 1.0;
  int divisions = 1;
};

AxisRange RunAxis(const std::vector<double>& runs) {
  AxisRange result;
  if (runs.empty()) return result;
  const auto limits = std::minmax_element(runs.begin(), runs.end());
  const double span = *limits.second - *limits.first;
  const double spacing = NiceRunSpacing(span);
  const double padding = std::max({1.0, 0.05 * span, 0.5 * spacing});
  result.low = std::floor((*limits.first - padding) / spacing) * spacing;
  result.high = std::ceil((*limits.second + padding) / spacing) * spacing;
  if (result.low == result.high) {
    result.low -= spacing;
    result.high += spacing;
  }
  result.divisions = std::max(1, static_cast<int>(std::lround(
    (result.high - result.low) / spacing)));
  return result;
}

void ConfigureCanvas(TCanvas* canvas) {
  canvas->Clear();
  canvas->SetLeftMargin(0.18);
  canvas->SetRightMargin(0.06);
  canvas->SetTopMargin(0.10);
  canvas->SetBottomMargin(0.12);
  canvas->SetGrid(1, 1);
}

void DrawMetadataPage(TCanvas* canvas, const std::string& inputPath,
                      const std::string& phase,
                      const std::vector<NormRow>& rows,
                      const std::string& pdfPath) {
  canvas->Clear();
  canvas->SetMargin(0.07, 0.05, 0.06, 0.06);
  canvas->SetGrid(0, 0);
  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextAlign(13);
  text.SetTextSize(0.024);
  double y = 0.94;
  auto line = [&](const std::string& value, double gap = 0.032) {
    text.DrawLatex(0.09, y, value.c_str());
    y -= gap;
  };
  text.SetTextFont(62);
  text.SetTextSize(0.046);
  line("RP coincidence normalized-yield analysis", 0.052);
  text.SetTextFont(42);
  text.SetTextSize(0.024);
  const std::string inputStem = InputStem(inputPath);
  line("Input setting:");
  if (inputStem.size() <= 58) {
    line("  " + inputStem);
  } else {
    line("  " + inputStem.substr(0, 58));
    line("  " + inputStem.substr(58));
  }
  line("Normalization lookup: matching setting-wise bigtable leaf CSV");
  line("Phase: " + phase + ", runs: " + std::to_string(rows.size()), 0.042);
  text.SetTextFont(62); line("Normalization", 0.038);
  text.SetTextFont(42);
  line("C_{norm} = PS #times boil_corr / (BCM2_Q #times comp_livetime #times #epsilon_{HMS} #times #epsilon_{SHMS}");
  line("                                      #times #epsilon_{trigger} #times #epsilon_{PID})");
  line("RP_Normyield = RP_Goodcoin #times C_{norm}");
  line("RP_Normyield_err = RP_Goodcoin_err #times C_{norm}", 0.042);
  line("BCM2_Q is in mC; #epsilon_{HMS}=h_esing_Eff; #epsilon_{SHMS}=p_hadron_Eff.");
  line("Exactly one of ps5 and ps6 must be positive, and it must equal 1.0.");
  line("Trigger efficiency = 1.0 and PID efficiency = 1.0 (provisional constants).", 0.042);
  text.SetTextColor(kRed + 1);
  text.SetTextFont(62);
  line("WARNING: comp_livetime = 1.0 is provisional for every run.", 0.042);
  text.SetTextColor(kBlack);
  text.SetTextFont(42);
  line("The bigtable comp_livetime is retained in the CSV but is not used yet.");
  line("RP_Normyield_err propagates RP_Goodcoin_err only.");
  line("Normalization systematic uncertainties are not included.", 0.042);
  line("normyield / RP_Normyield compares replay with this analysis.");
  line("The ratio combines counting and normalization differences.");
  line("The same comparison rule is applied to electron and positron runs.");

  size_t ok = 0, failed = 0, comparable = 0;
  for (const auto& row : rows) {
    ok += row.status == "OK";
    failed += row.status != "OK";
    comparable += std::isfinite(row.replayByRP);
  }
  y -= 0.006;
  text.SetTextFont(62);
  line("Validation summary", 0.038);
  text.SetTextFont(42);
  line("Valid: " + std::to_string(ok) + ", failed: " +
       std::to_string(failed) + ", replay comparisons: " +
       std::to_string(comparable));
  canvas->Print(pdfPath.c_str());
}

enum class PlotKind { kNormyield, kNormFactor, kRatio };

void AddYExtent(std::vector<double>& extents, double value, double error = 0.0) {
  if (!std::isfinite(value)) return;
  const double safeError = std::isfinite(error) && error >= 0.0 ? error : 0.0;
  extents.push_back(value - safeError);
  extents.push_back(value + safeError);
}

std::pair<double, double> PaddedYRange(const std::vector<double>& extents) {
  if (extents.empty()) return {0.0, 1.0};

  const auto limits = std::minmax_element(extents.begin(), extents.end());
  const double low = *limits.first;
  const double high = *limits.second;
  const double span = high - low;
  const double scale = std::max(std::abs(low), std::abs(high));
  double padding = std::max(0.15 * span, 0.05 * scale);
  if (!(padding > 0.0)) padding = 0.1;
  return {low - padding, high + padding};
}

void DrawDiagnosticPage(TCanvas* canvas, const std::vector<NormRow>& rows,
                        PlotKind kind, const std::string& pdfPath) {
  ConfigureCanvas(canvas);
  std::vector<double> allRuns, elecRun, elecValue, elecError;
  std::vector<double> posRun, posValue, posError;
  std::vector<double> replayRun, replayValue, replayError;
  std::vector<double> yExtents;
  for (const auto& row : rows) {
    if (row.status != "OK") continue;
    double value = kNaN, error = 0.0;
    if (kind == PlotKind::kNormyield) {
      value = row.rpNormyield;
      error = row.rpNormyieldErr;
      if (ValidMeasured(row.replayNormyield)) {
        replayRun.push_back(row.run);
        replayValue.push_back(row.replayNormyield);
        replayError.push_back(ValidMeasured(row.replayNormyieldErr)
                              ? row.replayNormyieldErr : 0.0);
        AddYExtent(yExtents, replayValue.back(), replayError.back());
      }
    } else if (kind == PlotKind::kNormFactor) {
      value = row.normFactor;
    } else {
      value = row.replayByRP;
    }
    if (!std::isfinite(value)) continue;
    allRuns.push_back(row.run);
    AddYExtent(yExtents, value, error);
    if (row.isPositron) {
      posRun.push_back(row.run); posValue.push_back(value); posError.push_back(error);
    } else {
      elecRun.push_back(row.run); elecValue.push_back(value); elecError.push_back(error);
    }
  }
  if (kind == PlotKind::kNormyield)
    allRuns.insert(allRuns.end(), replayRun.begin(), replayRun.end());
  const AxisRange axis = RunAxis(allRuns);
  const auto yRange = PaddedYRange(yExtents);
  const double yLow = yRange.first;
  const double yHigh = yRange.second;

  std::string title, yTitle;
  if (kind == PlotKind::kNormyield) {
    title = "Normalized Yield VS Run"; yTitle = "Normalized yield (mC^{-1})";
  } else if (kind == PlotKind::kNormFactor) {
    title = "Normalization Factor VS Run"; yTitle = "Normalization factor (mC^{-1})";
  } else {
    title = "normyield / RP_Normyield VS Run"; yTitle = "normyield / RP_Normyield";
  }
  TH1D frame("h_norm_diag_frame", (title + ";Run;" + yTitle).c_str(),
             1, axis.low, axis.high);
  frame.SetDirectory(nullptr);
  frame.SetStats(false);
  frame.SetMinimum(yLow);
  frame.SetMaximum(yHigh);
  frame.GetXaxis()->SetNdivisions(axis.divisions, false);
  frame.GetXaxis()->SetNoExponent(true);
  frame.GetXaxis()->SetDecimals(false);
  frame.GetXaxis()->SetLabelSize(0.028);
  frame.Draw();

  TGraphErrors electron(elecRun.size(), elecRun.data(), elecValue.data(),
                        nullptr, elecError.data());
  electron.SetMarkerStyle(20); electron.SetMarkerSize(1.3);
  electron.SetMarkerColor(kBlue + 1); electron.SetLineColor(kRed + 1);
  electron.SetLineWidth(2);
  TGraphErrors positron(posRun.size(), posRun.data(), posValue.data(),
                        nullptr, posError.data());
  positron.SetMarkerStyle(20); positron.SetMarkerSize(1.3);
  positron.SetMarkerColor(kRed + 1); positron.SetLineColor(kBlue + 1);
  positron.SetLineWidth(2);
  TGraphErrors replay(replayRun.size(), replayRun.data(), replayValue.data(),
                      nullptr, replayError.data());
  replay.SetMarkerStyle(20); replay.SetMarkerSize(1.2);
  replay.SetMarkerColor(kBlack); replay.SetLineColor(kBlack);
  replay.SetLineWidth(2);
  gStyle->SetEndErrorSize(6);
  if (!elecRun.empty()) electron.Draw("P same");
  if (!posRun.empty()) positron.Draw("P same");
  if (kind == PlotKind::kNormyield && !replayRun.empty()) replay.Draw("P same");
  if (kind == PlotKind::kNormyield) {
    TLegend legend(0.57, 0.74, 0.91, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.025);
    if (!elecRun.empty())
      legend.AddEntry(&electron, "RP_Normyield: electron", "lep");
    if (!posRun.empty())
      legend.AddEntry(&positron, "RP_Normyield: positron", "lep");
    if (!replayRun.empty())
      legend.AddEntry(&replay, "normyield: bigtable", "lep");
    legend.DrawClone("same");
  }
  if (kind == PlotKind::kRatio && axis.low < axis.high) {
    TLine reference(axis.low, 1.0, axis.high, 1.0);
    reference.SetLineColor(kBlack);
    reference.SetLineStyle(2);
    reference.SetLineWidth(2);
    reference.Draw();
  }
  canvas->Print(pdfPath.c_str());
}

} // namespace

int RP_get_coin_normyield(const std::string& inputCSVArgument = "") {
  gStyle->SetOptStat(0);
  gStyle->SetGridColor(kGray + 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  std::string inputPath = Trim(inputCSVArgument);
  if (inputPath.empty()) {
    std::cout << "Enter RP_get_good_coin_ev CSV path relative to rsidis_xs_v5:\n> ";
    std::getline(std::cin, inputPath);
    inputPath = Trim(inputPath);
  }
  if (inputPath.empty()) {
    std::cerr << "ERROR: no input CSV was provided\n";
    return 1;
  }

  CSVTable good;
  if (!ReadCSV(inputPath, good)) return 1;
  if (!HasColumns(good, {"Run", "RP_Goodcoin", "RP_Goodcoin_err", "hms_p"},
                  "Goodcoin table")) return 1;
  const std::string phase = PhaseFromPathAndRows(inputPath, good);
  const std::string bigtablePath = LeafBigtablePath(inputPath);
  if (bigtablePath.empty()) {
    std::cerr << "ERROR: cannot derive setting leaf CSV from input filename\n";
    return 1;
  }

  CSVTable bigtable;
  if (!ReadCSV(bigtablePath, bigtable)) return 1;
  if (!HasColumns(bigtable,
                  {"run", "BCM2_Q", "ps5", "ps6", "comp_livetime",
                   "h_esing_Eff", "p_hadron_Eff", "boil_corr",
                   "normyield", "normyield_err"}, "Setting leaf bigtable")) return 1;

  std::map<int, std::vector<const CSVRow*>> byRun;
  for (const auto& row : bigtable.rows) byRun[ParseRun(Get(row, "run"))].push_back(&row);
  std::vector<NormRow> outputs;
  outputs.reserve(good.rows.size());
  for (const auto& row : good.rows) {
    const int run = ParseRun(Get(row, "Run"));
    outputs.push_back(BuildNormRow(row, byRun[run]));
  }

  const std::string stem = InputStem(inputPath);
  const std::string outputCSV =
    "results/Tables/RP_get_coin_normyield/" + stem + ".csv";
  const std::string outputPDF =
    "results/PDFs/RP_get_coin_normyield/" + stem + ".pdf";
  gSystem->mkdir(gSystem->DirName(outputCSV.c_str()), true);
  gSystem->mkdir(gSystem->DirName(outputPDF.c_str()), true);

  WriteOutputCSV(outputCSV, good, outputs);
  TCanvas canvas("c_rp_coin_normyield", "RP coin normalized yield", 1200, 1200);
  canvas.SetCanvasSize(1200, 1200);
  canvas.Print((outputPDF + "[").c_str());
  DrawMetadataPage(&canvas, inputPath, phase, outputs, outputPDF);
  DrawDiagnosticPage(&canvas, outputs, PlotKind::kNormyield, outputPDF);
  DrawDiagnosticPage(&canvas, outputs, PlotKind::kNormFactor, outputPDF);
  DrawDiagnosticPage(&canvas, outputs, PlotKind::kRatio, outputPDF);
  canvas.Print((outputPDF + "]").c_str());

  size_t failed = 0;
  for (const auto& row : outputs) failed += row.status != "OK";
  std::cout << "Runs: " << outputs.size() << ", normalization failures: "
            << failed << "\nWrote " << outputCSV << "\nWrote " << outputPDF
            << "\n";
  return failed == outputs.size() && !outputs.empty() ? 2 : 0;
}
