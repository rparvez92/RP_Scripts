// hodo_calib_qc_batch.C
//
// Phase-1/Phase-2 hodoscope calibration QC using unskimmed replay ROOT files.
//
// Supported modes:
//   hms           : HMS electron selection, fit H.gtr.beta
//   shms          : SHMS electron selection, fit P.gtr.beta
//   coin          : HMS electron + SHMS pion selection, fit CTime.ePiCoinTime_ROC2
//
// Input filename conventions:
//   hms           : hms_coin_replay_production_<run>_-1.root
//   shms_*        : shms_coin_replay_production_<run>_-1.root
//   coin          : coin_replay_production_<run>_-1.root
//
// Standard full-dataset output:
//   results/Phase<phase>/tables/hodo_qc_<spec>_summary_phase<phase>.csv
//   results/Phase<phase>/pdfs/hodo_qc_<spec>_by_run_phase<phase>.pdf
// A nonempty OutputSuffix replaces the standard _phase<phase> tag, preserving
// established names such as hodo_qc_<spec>_by_run_QA.pdf.
//
// Compile only:
//   root -l -b -q -e '.L macros/hodo_calib_qc_batch.C+'
//
// Process a complete phase (SHMSDIS, HMSDIS, then COIN):
//   root -l -b -q 'macros/hodo_calib_qc_batch.C+(1)'
//   root -l -b -q 'macros/hodo_calib_qc_batch.C+(2)'
// Process an explicit sample while keeping the standard output names:
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch.C+(1, "", "", "", false, false, "23853,23856,23861")'
//
// Process only files available in an alternate ROOT directory, writing _QA
// outputs (set the final argument true first for a selection-only preflight):
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch.C+(2, "bigtable/rsidis_bigtable_phase2.csv", "/path/to/ROOTfiles", "_QA", true, false)'

#include <TCanvas.h>
#include <TCut.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

const char *kDefaultRootDirPhase1 =
    "/cache/hallc/c-rsidis/analysis/replays/pass0p1";
const char *kDefaultRootDirPhase2 =
    "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles";
const char *kDefaultBigtablePhase1 =
    "bigtable/rsidis_bigtable_pass0p1.csv";
const char *kDefaultBigtablePhase2 =
    "bigtable/rsidis_bigtable_phase2.csv";

const int kCanvasWidth = 1250;
const int kCanvasHeight = 850;

struct RunMetadata {
  int run = 0;
  TString runType;
  TString target;
  double hmsP = std::nan("");
  double shmsP = std::nan("");
};

struct RunGroups {
  std::vector<RunMetadata> coin;
  std::vector<RunMetadata> hms;
  std::vector<RunMetadata> shms;
  int excludedPolarity = 0;
  int excludedRunType = 0;
};

struct FitDiagnostics {
  bool fitAttempted = false;
  int minimizerStatus = -1;
  int fitValid = -1;
  int covarianceStatus = -1;
  double chi2 = std::nan("");
  double ndf = std::nan("");
  double chi2Ndf = std::nan("");
  double edm = std::nan("");
  double candidateMean = std::nan("");
  double candidateSigma = std::nan("");
  double fitLow = std::nan("");
  double fitHigh = std::nan("");
  TString failureReason;
};

struct RunSummary {
  RunMetadata metadata;
  TString spec;
  TString filePath;
  TString fitVariable;
  Long64_t allEvents = -1;
  Long64_t selectedEvents = -1;
  double fitMean = std::nan("");
  double fitSigma = std::nan("");
  double fitEntries = std::nan("");
  FitDiagnostics fitDiagnostics;
  TString status = "NOT_RUN";
};

TString MakeFileName(const TString &spec, int run) {
  if (spec == "hms")
    return TString::Format("hms_coin_replay_production_%d_-1.root", run);
  if (spec == "shms")
    return TString::Format("shms_coin_replay_production_%d_-1.root", run);
  if (spec == "coin")
    return TString::Format("coin_replay_production_%d_-1.root", run);
  return "";
}

bool ValidateOutputSuffix(const TString &suffix) {
  for (Ssiz_t index = 0; index < suffix.Length(); ++index) {
    const unsigned char ch =
        static_cast<unsigned char>(suffix[index]);
    if (!std::isalnum(ch) && ch != '_' && ch != '-') {
      std::cerr << "[ERROR] Unsafe OutputSuffix '" << suffix
                << "'. Use only letters, digits, '_' or '-'.\n";
      return false;
    }
  }
  return true;
}

std::vector<RunMetadata> KeepAvailableRuns(
    const TString &spec, const TString &rootDir,
    const std::vector<RunMetadata> &runs) {
  std::vector<RunMetadata> available;
  available.reserve(runs.size());
  for (const RunMetadata &run : runs) {
    const TString path = TString::Format(
        "%s/%s", rootDir.Data(), MakeFileName(spec, run.run).Data());
    if (!gSystem->AccessPathName(path))
      available.push_back(run);
  }
  return available;
}

std::string Trim(const std::string &text) {
  const std::string whitespace = " \t\r\n";
  const auto first = text.find_first_not_of(whitespace);
  if (first == std::string::npos)
    return "";
  return text.substr(first, text.find_last_not_of(whitespace) - first + 1);
}

std::vector<std::string> ParseCsvLine(const std::string &line, bool &valid) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  valid = true;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];
    if (ch == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field += '"';
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (ch == ',' && !quoted) {
      fields.push_back(Trim(field));
      field.clear();
    } else {
      field += ch;
    }
  }
  if (quoted)
    valid = false;
  fields.push_back(Trim(field));
  return fields;
}

bool ParseIntStrict(const std::string &text, int &value) {
  char *end = nullptr;
  const long parsed = std::strtol(text.c_str(), &end, 10);
  if (text.empty() || !end || *end != '\0')
    return false;
  value = static_cast<int>(parsed);
  return true;
}

bool ParseDoubleStrict(const std::string &text, double &value) {
  char *end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return !text.empty() && end && *end == '\0' && std::isfinite(value);
}

bool ParseRunsList(const TString &text, std::set<int> &runs) {
  std::stringstream input(text.Data());
  std::string token;
  while (std::getline(input, token, ',')) {
    token = Trim(token);
    if (token.empty()) {
      std::cerr << "[ERROR] Empty item in RunsList '" << text << "'.\n";
      return false;
    }
    const std::size_t dash = token.find('-');
    int first = 0;
    int last = 0;
    if (dash == std::string::npos) {
      if (!ParseIntStrict(token, first)) {
        std::cerr << "[ERROR] Invalid run '" << token << "' in RunsList.\n";
        return false;
      }
      last = first;
    } else {
      if (token.find('-', dash + 1) != std::string::npos ||
          !ParseIntStrict(Trim(token.substr(0, dash)), first) ||
          !ParseIntStrict(Trim(token.substr(dash + 1)), last) ||
          last < first || last - first > 100000) {
        std::cerr << "[ERROR] Invalid run range '" << token
                  << "' in RunsList.\n";
        return false;
      }
    }
    for (int run = first; run <= last; ++run) {
      if (run <= 0 || !runs.insert(run).second) {
        std::cerr << "[ERROR] Invalid or duplicate requested run " << run
                  << " in RunsList.\n";
        return false;
      }
    }
  }
  if (runs.empty()) {
    std::cerr << "[ERROR] RunsList did not contain any runs.\n";
    return false;
  }
  return true;
}

bool ReadBigtable(const TString &path, RunGroups &groups) {
  std::ifstream input(path.Data());
  if (!input) {
    std::cerr << "[ERROR] Cannot open bigtable: " << path << '\n';
    return false;
  }

  std::string line;
  if (!std::getline(input, line)) {
    std::cerr << "[ERROR] Bigtable is empty: " << path << '\n';
    return false;
  }
  bool valid = true;
  const auto headers = ParseCsvLine(line, valid);
  if (!valid) {
    std::cerr << "[ERROR] Malformed quoted CSV header in " << path << '\n';
    return false;
  }
  std::map<std::string, std::size_t> columns;
  for (std::size_t i = 0; i < headers.size(); ++i)
    columns[Trim(headers[i])] = i;
  const std::vector<std::string> required = {
      "run", "run_type", "target", "hms_p", "shms_p"};
  for (const auto &name : required) {
    if (!columns.count(name)) {
      std::cerr << "[ERROR] Bigtable is missing required column '"
                << name << "'.\n";
      return false;
    }
  }

  std::set<int> selectedRuns;
  int lineNumber = 1;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (Trim(line).empty())
      continue;
    const auto fields = ParseCsvLine(line, valid);
    if (!valid || fields.size() != headers.size()) {
      std::cerr << "[ERROR] Malformed CSV row at line " << lineNumber
                << ": expected " << headers.size() << " fields, got "
                << fields.size() << ".\n";
      return false;
    }

    RunMetadata row;
    const std::string runText = fields[columns["run"]];
    const std::string hmsText = fields[columns["hms_p"]];
    const std::string shmsText = fields[columns["shms_p"]];
    if (!ParseIntStrict(runText, row.run) ||
        !ParseDoubleStrict(hmsText, row.hmsP) ||
        !ParseDoubleStrict(shmsText, row.shmsP)) {
      std::cerr << "[ERROR] Invalid run/hms_p/shms_p value at bigtable line "
                << lineNumber << ".\n";
      return false;
    }
    row.runType = fields[columns["run_type"]];
    row.target = fields[columns["target"]];

    TString spec;
    if (row.runType == "PI+SIDIS" || row.runType == "PI-SIDIS")
      spec = "coin";
    else if (row.runType == "HMSDIS")
      spec = "hms";
    else if (row.runType == "SHMSDIS")
      spec = "shms";
    else {
      ++groups.excludedRunType;
      continue;
    }
    if (row.hmsP >= 0.0) {
      ++groups.excludedPolarity;
      continue;
    }
    if (!selectedRuns.insert(row.run).second) {
      std::cerr << "[ERROR] Duplicate selected run " << row.run
                << " at bigtable line " << lineNumber << ".\n";
      return false;
    }
    if (spec == "coin")
      groups.coin.push_back(row);
    else if (spec == "hms")
      groups.hms.push_back(row);
    else
      groups.shms.push_back(row);
  }

  const auto byRun = [](const RunMetadata &a, const RunMetadata &b) {
    return a.run < b.run;
  };
  std::sort(groups.coin.begin(), groups.coin.end(), byRun);
  std::sort(groups.hms.begin(), groups.hms.end(), byRun);
  std::sort(groups.shms.begin(), groups.shms.end(), byRun);
  return true;
}

bool ApplyRunFilter(const TString &text, RunGroups &groups) {
  std::set<int> requested;
  if (!ParseRunsList(text, requested))
    return false;

  std::set<int> matched;
  const auto filter = [&](std::vector<RunMetadata> &runs) {
    runs.erase(std::remove_if(runs.begin(), runs.end(), [&](const RunMetadata &row) {
                 if (!requested.count(row.run))
                   return true;
                 matched.insert(row.run);
                 return false;
               }),
               runs.end());
  };
  filter(groups.coin);
  filter(groups.hms);
  filter(groups.shms);

  if (matched != requested) {
    std::cerr << "[ERROR] RunsList contains runs not selected from the "
              << "requested phase bigtable:";
    for (int run : requested) {
      if (!matched.count(run))
        std::cerr << ' ' << run;
    }
    std::cerr << "\n";
    return false;
  }
  return true;
}

TCut HmsElectronCuts() {
  return TCut("(H.gtr.dp>-8) && (H.gtr.dp<8)"
              " && (H.cal.etottracknorm>0.7)"
              " && (H.cer.npeSum>2.0)");
}

TCut ShmsElectronCuts() {
  return TCut("(P.gtr.dp>-10) && (P.gtr.dp<22)"
              " && (P.cal.etottracknorm>0.7)");
}

TCut ShmsPionCuts() {
  const TCut shmsBase =
      "(P.gtr.dp>-10) && (P.gtr.dp<22)"
      " && (P.cal.etottracknorm<0.8)";
  const TCut shmsAero = "(P.gtr.p<2.7) && (P.aero.npeSum>2)";
  const TCut shmsHgc =
      "(P.gtr.p>=2.7) && (P.hgcer.npeSum>1)"
      " && (P.aero.npeSum>2)";
  return shmsBase && (shmsAero || shmsHgc);
}

TCut BuildCuts(const TString &spec) {
  if (spec == "hms")
    return HmsElectronCuts();
  if (spec == "shms")
    return ShmsElectronCuts();
  if (spec == "coin")
    return HmsElectronCuts() && ShmsPionCuts();
  return TCut("");
}

void PrintPhysicsLogic(int phase) {
  std::cout
      << "\n===== PHASE-" << phase
      << " HODOSCOPE QC PHYSICS LOGIC =====\n"
      << "HMS electron selection:\n"
      << "  -8 < H.gtr.dp < 8\n"
      << "  H.cal.etottracknorm > 0.7\n"
      << "  H.cer.npeSum > 2.0\n\n"
      << "SHMS electron selection:\n"
      << "  -10 < P.gtr.dp < 22\n"
      << "  P.cal.etottracknorm > 0.7\n"
      << "  No P.ngcer/P.hgcer electron cuts are applied.\n\n"
      << "SHMS pion selection:\n"
      << "  -10 < P.gtr.dp < 22\n"
      << "  P.cal.etottracknorm < 0.8\n"
      << "  P.gtr.p < 2.7: P.aero.npeSum > 2\n"
      << "  P.gtr.p >= 2.7: P.hgcer.npeSum > 1 AND "
         "P.aero.npeSum > 2\n\n"
      << "COIN selection: HMS electron AND SHMS pion selections above.\n"
      << "No coincidence-time gate is applied.\n"
      << "Beta Gaussian fit: +/-0.03 around peak, bounded to [0.9,1.1].\n"
      << "Coin-time Gaussian + constant-background fit: +/-0.75 ns around "
         "peak, bounded to [0,100].\n"
      << "Fits require at least 50 selected entries.\n"
      << "Beta lines at 0.95 and 1.05 are visual guides only.\n"
      << "================================================\n\n";
}

std::vector<const char *> RequiredBranches(const TString &spec) {
  const std::vector<const char *> hms = {
      "H.gtr.dp", "H.gtr.beta", "H.cal.etottracknorm",
      "H.cer.npeSum", "H.dc.x_fp"};
  const std::vector<const char *> shms = {
      "P.gtr.dp", "P.gtr.beta", "P.cal.etottracknorm", "P.dc.x_fp"};
  const std::vector<const char *> shmsPion = {
      "P.gtr.dp", "P.gtr.beta", "P.cal.etottracknorm",
      "P.hgcer.npeSum", "P.aero.npeSum", "P.gtr.p", "P.dc.x_fp"};

  if (spec == "hms")
    return hms;
  if (spec == "shms")
    return shms;

  std::vector<const char *> coin = hms;
  coin.insert(coin.end(), shmsPion.begin(), shmsPion.end());
  coin.push_back("CTime.ePiCoinTime_ROC2");
  return coin;
}

bool ValidateAndEnableBranches(TTree *tree, const TString &spec, int run) {
  bool valid = true;
  const std::vector<const char *> branches = RequiredBranches(spec);
  for (const char *name : branches) {
    if (!tree->GetBranch(name)) {
      std::cerr << "[ERROR] Run " << run << " is missing required branch "
                << name << '\n';
      valid = false;
    }
  }
  if (!valid)
    return false;

  tree->SetBranchStatus("*", 0);
  for (const char *name : branches)
    tree->SetBranchStatus(name, 1);
  return true;
}

TString OutputTag(int phase, const TString &suffix) {
  if (!suffix.IsNull())
    return suffix;
  return TString::Format("_phase%d", phase);
}

TString PdfPath(const TString &spec, int phase, const TString &suffix = "") {
  return TString::Format(
      "results/Phase%d/pdfs/hodo_qc_%s_by_run%s.pdf", phase,
      spec.Data(), OutputTag(phase, suffix).Data());
}

TString CsvPath(const TString &spec, int phase, const TString &suffix = "") {
  return TString::Format(
      "results/Phase%d/tables/hodo_qc_%s_summary%s.csv", phase,
      spec.Data(), OutputTag(phase, suffix).Data());
}

void SaveCanvas(TCanvas &canvas, const TString &pdfPath) {
  canvas.Print(pdfPath);
}

void OpenPdf(const TString &pdfPath) {
  TCanvas opener("c_pdf_open", "", 1, 1);
  opener.Print(TString::Format("%s[", pdfPath.Data()));
}

void ClosePdf(const TString &pdfPath) {
  TCanvas closer("c_pdf_close", "", 1, 1);
  closer.Print(TString::Format("%s]", pdfPath.Data()));
}

void DrawBetaVsXfp(TTree *tree, const TString &selectionSpec,
                   const TString &viewSpec, int phase, int run,
                   const TString &pdfPath) {
  const bool hmsView = viewSpec == "hms";
  const TString expression =
      hmsView ? "H.gtr.beta:H.dc.x_fp" : "P.gtr.beta:P.dc.x_fp";
  const TString histName =
      TString::Format("h_beta_xfp_%s_%d", viewSpec.Data(), run);
  const TString arm = hmsView ? "HMS" : "SHMS";
  const TString xVariable = hmsView ? "H.dc.x_fp" : "P.dc.x_fp";
  const TString betaVariable = hmsView ? "H.gtr.beta" : "P.gtr.beta";
  const TString title = TString::Format(
      "Phase %d run %d: %s beta vs %s xfp;%s [cm];%s", phase,
      run, arm.Data(), arm.Data(), xVariable.Data(), betaVariable.Data());

  TH2D hist(histName, title, 80, -45, 45, 120, 0.2, 1.2);
  hist.Sumw2();
  tree->Project(histName, expression, BuildCuts(selectionSpec));

  TCanvas canvas(TString::Format("c_beta_xfp_%s_%d", viewSpec.Data(), run),
                 "", kCanvasWidth, kCanvasHeight);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.18);
  canvas.SetBottomMargin(0.13);
  canvas.SetTopMargin(0.10);
  gStyle->SetOptStat(0);
  hist.GetXaxis()->SetTitleOffset(1.15);
  hist.GetYaxis()->SetTitleOffset(1.15);
  hist.GetXaxis()->SetLabelSize(0.035);
  hist.GetYaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetTitleOffset(1.25);
  hist.Draw("COLZ");
  const double xmin = hist.GetXaxis()->GetXmin();
  const double xmax = hist.GetXaxis()->GetXmax();
  TLine low(xmin, 0.95, xmax, 0.95);
  TLine high(xmin, 1.05, xmax, 1.05);
  low.SetLineStyle(2);
  high.SetLineStyle(2);
  low.SetLineColor(kBlack);
  high.SetLineColor(kBlack);
  low.SetLineWidth(5);
  high.SetLineWidth(5);
  low.Draw("SAME");
  high.Draw("SAME");

  SaveCanvas(canvas, pdfPath);
}

void DrawCoinTimeVsHmsXfp(TTree *tree, int phase, int run,
                          const TString &pdfPath) {
  const TString histName = TString::Format("h_ctime_hms_xfp_%d", run);
  TH2D hist(
      histName,
      TString::Format(
          "Phase %d run %d: CTime (ROC2) vs HMS xfp;"
          "H.dc.x_fp [cm];CTime.ePiCoinTime_ROC2 [ns]",
          phase, run),
      80, -45, 45, 400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2:H.dc.x_fp",
                BuildCuts("coin"));

  TCanvas canvas(TString::Format("c_ctime_hms_xfp_%d", run), "",
                 kCanvasWidth, kCanvasHeight);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.18);
  canvas.SetBottomMargin(0.13);
  canvas.SetTopMargin(0.10);
  gStyle->SetOptStat(0);
  hist.GetXaxis()->SetTitleOffset(1.15);
  hist.GetYaxis()->SetTitleOffset(1.15);
  hist.GetXaxis()->SetLabelSize(0.035);
  hist.GetYaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetTitleOffset(1.25);
  hist.Draw("COLZ");
  SaveCanvas(canvas, pdfPath);
}

void DrawCoinTimeVsShmsXfp(TTree *tree, int phase, int run,
                           const TString &pdfPath) {
  const TString histName = TString::Format("h_ctime_shms_xfp_%d", run);
  TH2D hist(
      histName,
      TString::Format(
          "Phase %d run %d: CTime (ROC2) vs SHMS xfp;"
          "P.dc.x_fp [cm];CTime.ePiCoinTime_ROC2 [ns]",
          phase, run),
      80, -45, 45, 400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2:P.dc.x_fp",
                BuildCuts("coin"));

  TCanvas canvas(TString::Format("c_ctime_shms_xfp_%d", run), "",
                 kCanvasWidth, kCanvasHeight);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.18);
  canvas.SetBottomMargin(0.13);
  canvas.SetTopMargin(0.10);
  gStyle->SetOptStat(0);
  hist.GetXaxis()->SetTitleOffset(1.15);
  hist.GetYaxis()->SetTitleOffset(1.15);
  hist.GetXaxis()->SetLabelSize(0.035);
  hist.GetYaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetLabelSize(0.035);
  hist.GetZaxis()->SetTitleOffset(1.25);
  hist.Draw("COLZ");
  SaveCanvas(canvas, pdfPath);
}

bool CaptureFitDiagnostics(const TFitResultPtr &result, const TF1 &fit,
                           bool enforceBounds, FitDiagnostics &diagnostics) {
  diagnostics.fitAttempted = true;
  diagnostics.minimizerStatus = static_cast<int>(result);
  diagnostics.candidateMean = fit.GetParameter(1);
  diagnostics.candidateSigma = std::abs(fit.GetParameter(2));
  diagnostics.chi2 = fit.GetChisquare();
  diagnostics.ndf = fit.GetNDF();
  if (diagnostics.ndf > 0.0)
    diagnostics.chi2Ndf = diagnostics.chi2 / diagnostics.ndf;

  if (result.Get()) {
    diagnostics.fitValid = result->IsValid() ? 1 : 0;
    diagnostics.covarianceStatus = result->CovMatrixStatus();
    diagnostics.edm = result->Edm();
  } else {
    diagnostics.fitValid = 0;
  }

  if (diagnostics.minimizerStatus != 0) {
    diagnostics.failureReason = TString::Format(
        "MINIMIZER_STATUS_%d", diagnostics.minimizerStatus);
  } else if (diagnostics.fitValid != 1) {
    diagnostics.failureReason = "INVALID_FIT_RESULT";
  } else if (!std::isfinite(diagnostics.candidateMean) ||
             !std::isfinite(diagnostics.candidateSigma)) {
    diagnostics.failureReason = "NONFINITE_PARAMETERS";
  } else if (enforceBounds &&
             (diagnostics.candidateMean < diagnostics.fitLow ||
              diagnostics.candidateMean > diagnostics.fitHigh)) {
    diagnostics.failureReason = "MEAN_OUT_OF_RANGE";
  } else if (enforceBounds &&
             (diagnostics.candidateSigma < 0.05 ||
              diagnostics.candidateSigma > 2.0)) {
    diagnostics.failureReason = "SIGMA_OUT_OF_RANGE";
  } else if (enforceBounds &&
             (std::abs(diagnostics.candidateMean - diagnostics.fitLow) < 1e-4 ||
              std::abs(diagnostics.candidateMean - diagnostics.fitHigh) < 1e-4)) {
    diagnostics.failureReason = "MEAN_AT_LIMIT";
  } else if (enforceBounds &&
             (std::abs(diagnostics.candidateSigma - 0.05) < 1e-4 ||
              std::abs(diagnostics.candidateSigma - 2.0) < 1e-4)) {
    diagnostics.failureReason = "SIGMA_AT_LIMIT";
  } else {
    diagnostics.failureReason = "";
  }
  return diagnostics.failureReason.IsNull();
}

void PrintFitFailure(int run, const TString &variable, double entries,
                     const FitDiagnostics &diagnostics) {
  std::cerr << "[WARN] Fit rejected: run=" << run
            << ", variable=" << variable
            << ", entries=" << entries
            << ", window=[" << diagnostics.fitLow << ','
            << diagnostics.fitHigh << ']'
            << ", attempted=" << (diagnostics.fitAttempted ? 1 : 0)
            << ", status=" << diagnostics.minimizerStatus
            << ", valid=" << diagnostics.fitValid
            << ", covariance_status=" << diagnostics.covarianceStatus
            << ", candidate_mean=" << diagnostics.candidateMean
            << ", candidate_sigma=" << diagnostics.candidateSigma
            << ", chi2/ndf=" << diagnostics.chi2Ndf
            << ", edm=" << diagnostics.edm
            << ", reason=" << diagnostics.failureReason << '\n';
}

bool ComputeBetaMetrics(TTree *tree, const TString &spec, int run,
                        double &mean, double &sigma, double &entries,
                        FitDiagnostics &diagnostics) {
  mean = sigma = entries = std::nan("");
  const TString variable = spec == "hms" ? "H.gtr.beta" : "P.gtr.beta";
  const TString histName = TString::Format("h_beta_fit_%s_%d", spec.Data(), run);
  TH1D hist(histName, TString::Format(";%s;Counts", variable.Data()),
            200, 0.2, 1.2);
  hist.Sumw2();
  tree->Project(histName, variable, BuildCuts(spec));
  entries = hist.GetEntries();
  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  diagnostics.fitLow = std::max(0.9, peak - 0.03);
  diagnostics.fitHigh = std::min(1.1, peak + 0.03);
  if (entries < 50) {
    diagnostics.failureReason = "LOW_STATISTICS";
    return false;
  }
  if (diagnostics.fitHigh <= diagnostics.fitLow) {
    diagnostics.fitValid = 0;
    diagnostics.failureReason = "INVALID_FIT_WINDOW";
    return false;
  }

  TF1 fit(TString::Format("f_beta_%s_%d", spec.Data(), run),
          "gaus", diagnostics.fitLow, diagnostics.fitHigh);
  const TFitResultPtr result = hist.Fit(&fit, "QNRS");
  if (!CaptureFitDiagnostics(result, fit, false, diagnostics))
    return false;
  mean = diagnostics.candidateMean;
  sigma = diagnostics.candidateSigma;
  return true;
}

bool FitCoinTimePeak(TH1D &hist, TF1 &fit,
                     double fitLow, double peak, double fitHigh,
                     FitDiagnostics &diagnostics) {
  diagnostics.fitLow = fitLow;
  diagnostics.fitHigh = fitHigh;
  const double maximum = hist.GetMaximum();
  double background = maximum;
  const int firstBin = hist.GetXaxis()->FindFixBin(fitLow);
  const int lastBin = hist.GetXaxis()->FindFixBin(fitHigh);
  for (int bin = firstBin; bin <= lastBin; ++bin)
    background = std::min(background, hist.GetBinContent(bin));
  background = std::max(0.0, background);
  fit.SetParameters(std::max(1.0, maximum - background), peak, 0.35,
                    background);
  fit.SetParLimits(0, 0.0, std::max(1.0, maximum * 10.0));
  fit.SetParLimits(1, fitLow, fitHigh);
  fit.SetParLimits(2, 0.05, 2.0);
  fit.SetParLimits(3, 0.0, std::max(1.0, maximum * 2.0));

  const TFitResultPtr result = hist.Fit(&fit, "QNRSB");
  return CaptureFitDiagnostics(result, fit, true, diagnostics);
}

void DrawCoinTime1D(TTree *tree, int phase, int run,
                    const TString &pdfPath) {
  const TString histName = TString::Format("h_ctime_%d", run);
  TH1D hist(histName,
            TString::Format(
                "Phase %d run %d: Coincidence Time (ROC2);"
                "CTime.ePiCoinTime_ROC2 [ns];Counts",
                phase, run),
            400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2", BuildCuts("coin"));

  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  const double fitLow = std::max(0.0, peak - 0.75);
  const double fitHigh = std::min(100.0, peak + 0.75);
  TF1 fit(TString::Format("f_ctime_%d", run), "gaus(0)+pol0(3)",
          fitLow, fitHigh);
  FitDiagnostics displayDiagnostics;
  const bool fitValid = FitCoinTimePeak(
      hist, fit, fitLow, peak, fitHigh, displayDiagnostics);

  TCanvas canvas(TString::Format("c_ctime_%d", run), "", kCanvasWidth,
                 kCanvasHeight);
  canvas.SetLeftMargin(0.13);
  canvas.SetRightMargin(0.05);
  canvas.SetBottomMargin(0.13);
  canvas.SetTopMargin(0.10);
  gStyle->SetOptStat(0);
  hist.GetXaxis()->SetTitleOffset(1.15);
  hist.GetYaxis()->SetTitleOffset(1.15);
  hist.GetXaxis()->SetLabelSize(0.035);
  hist.GetYaxis()->SetLabelSize(0.035);
  hist.Draw("HIST");
  if (fitValid) {
    fit.SetLineColor(kRed);
    fit.SetLineWidth(3);
    fit.Draw("SAME");
  } else {
    std::cerr << "[WARN] Coin-time display fit rejected for run "
              << run << '\n';
  }

  SaveCanvas(canvas, pdfPath);
}

bool ComputeCoinTimeMetrics(TTree *tree, int run, double &mean,
                            double &sigma, double &entries,
                            FitDiagnostics &diagnostics) {
  mean = sigma = entries = std::nan("");
  const TString histName = TString::Format("h_ctime_fit_%d", run);
  TH1D hist(histName, ";CTime.ePiCoinTime_ROC2 [ns];Counts", 400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2", BuildCuts("coin"));
  entries = hist.GetEntries();
  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  const double fitLow = std::max(0.0, peak - 0.75);
  const double fitHigh = std::min(100.0, peak + 0.75);
  diagnostics.fitLow = fitLow;
  diagnostics.fitHigh = fitHigh;
  if (entries < 50) {
    diagnostics.failureReason = "LOW_STATISTICS";
    return false;
  }

  TF1 fit(TString::Format("f_ctime_metric_%d", run),
          "gaus(0)+pol0(3)", fitLow, fitHigh);
  if (!FitCoinTimePeak(hist, fit, fitLow, peak, fitHigh, diagnostics))
    return false;
  mean = diagnostics.candidateMean;
  sigma = diagnostics.candidateSigma;
  return true;
}

void DrawDualTrend(const std::vector<int> &runs,
                   const std::vector<double> &means,
                   const std::vector<double> &sigmas,
                   const TString &spec, int phase, bool coinTime,
                   const TString &pdfPath) {
  if (runs.empty())
    return;

  const int count = static_cast<int>(runs.size());
  const double meanLow = coinTime ? 45.0 : 0.5;
  const double meanHigh = coinTime ? 55.0 : 1.5;
  const double sigmaLow = 0.0;
  const double sigmaHigh = coinTime ? 1.0 : 0.1;
  const TString frameName =
      TString::Format("h_trend_%s", coinTime ? "ctime" : spec.Data());
  const TString title =
      coinTime ? TString::Format(
                     "Phase %d COIN: CTime (ROC2) mean / sigma vs run;Run;"
                     "CTime.ePiCoinTime_ROC2 mean [ns]",
                     phase)
               : TString::Format(
                     "Phase %d %s: beta mean / sigma vs run;Run;%s mean",
                     phase, spec.Data(),
                     spec == "hms" ? "H.gtr.beta" : "P.gtr.beta");

  TCanvas canvas(TString::Format("c_trend_%s", spec.Data()), "",
                 kCanvasWidth, kCanvasHeight);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.13);
  canvas.SetBottomMargin(0.18);
  canvas.SetTopMargin(0.10);
  TH1F frame(frameName, title, count, 0.0, static_cast<double>(count));
  const int labelStep = std::max(1, (count + 19) / 20);
  for (int index = 0; index < count; ++index) {
    if (index % labelStep == 0)
      frame.GetXaxis()->SetBinLabel(index + 1,
                                    TString::Format("%d", runs[index]));
  }
  frame.SetMinimum(meanLow);
  frame.SetMaximum(meanHigh);
  frame.GetXaxis()->LabelsOption("v");
  frame.GetXaxis()->SetLabelSize(0.030);
  frame.GetYaxis()->SetTitleOffset(1.20);
  frame.GetYaxis()->SetLabelSize(0.035);
  frame.Draw("HIST");

  std::vector<double> x(count), scaledSigma(count);
  for (int index = 0; index < count; ++index) {
    x[index] = index + 0.5;
    scaledSigma[index] =
        meanLow + (sigmas[index] - sigmaLow) *
                      (meanHigh - meanLow) / (sigmaHigh - sigmaLow);
  }

  TGraph meanGraph(count, x.data(), means.data());
  TGraph sigmaGraph(count, x.data(), scaledSigma.data());
  meanGraph.SetMarkerStyle(20);
  meanGraph.SetMarkerSize(1.1);
  meanGraph.SetMarkerColor(kBlack);
  sigmaGraph.SetMarkerStyle(22);
  sigmaGraph.SetMarkerSize(1.1);
  sigmaGraph.SetMarkerColor(kBlue + 1);
  meanGraph.Draw("P SAME");
  sigmaGraph.Draw("P SAME");

  TGaxis rightAxis(count, meanLow, count, meanHigh,
                   sigmaLow, sigmaHigh, 510, "+L");
  rightAxis.SetTitle(
      coinTime ? "CTime.ePiCoinTime_ROC2 sigma [ns]"
               : (spec == "hms" ? "H.gtr.beta sigma"
                                 : "P.gtr.beta sigma"));
  rightAxis.SetTitleOffset(1.20);
  rightAxis.SetLabelSize(0.035);
  rightAxis.Draw();

  TLegend legend(0.12, 0.84, 0.24, 0.92);
  legend.AddEntry(&meanGraph, coinTime ? "CTime mean" : "beta mean", "p");
  legend.AddEntry(&sigmaGraph, coinTime ? "CTime sigma" : "beta sigma", "p");
  legend.Draw();

  SaveCanvas(canvas, pdfPath);
}

RunSummary ProcessOneRun(const TString &spec, const TString &rootDir,
                         int phase, const RunMetadata &metadata,
                         const TString &pdfPath,
                         std::vector<int> &trendRuns,
                         std::vector<double> &means,
                         std::vector<double> &sigmas) {
  RunSummary summary;
  summary.metadata = metadata;
  const int run = metadata.run;
  summary.spec = spec;
  summary.fitVariable =
      spec == "coin" ? "CTime.ePiCoinTime_ROC2"
                     : (spec == "hms" ? "H.gtr.beta" : "P.gtr.beta");
  summary.filePath =
      TString::Format("%s/%s", rootDir.Data(), MakeFileName(spec, run).Data());

  if (gSystem->AccessPathName(summary.filePath)) {
    summary.status = "MISSING_FILE";
    std::cerr << "[WARN] Missing file for run " << run << ": "
              << summary.filePath << '\n';
    return summary;
  }

  TFile *file = TFile::Open(summary.filePath, "READ");
  if (!file || file->IsZombie()) {
    summary.status = "ZOMBIE_FILE";
    std::cerr << "[WARN] Could not open run " << run << ": "
              << summary.filePath << '\n';
    if (file)
      delete file;
    return summary;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("T"));
  if (!tree) {
    summary.status = "MISSING_TREE";
    std::cerr << "[WARN] Tree 'T' missing for run " << run << ": "
              << summary.filePath << '\n';
    file->Close();
    delete file;
    return summary;
  }
  if (!ValidateAndEnableBranches(tree, spec, run)) {
    summary.status = "MISSING_BRANCH";
    file->Close();
    delete file;
    return summary;
  }

  summary.allEvents = tree->GetEntries();
  summary.selectedEvents = tree->GetEntries(BuildCuts(spec));
  std::cout << "[RUN " << run << "] all=" << summary.allEvents
            << ", selected=" << summary.selectedEvents << '\n';

  if (spec == "hms" || spec == "shms") {
    const TString viewSpec = spec == "hms" ? "hms" : "shms";
    DrawBetaVsXfp(tree, spec, viewSpec, phase, run, pdfPath);
    if (ComputeBetaMetrics(tree, spec, run, summary.fitMean,
                           summary.fitSigma, summary.fitEntries,
                           summary.fitDiagnostics)) {
      summary.status = "OK";
      trendRuns.push_back(run);
      means.push_back(summary.fitMean);
      sigmas.push_back(summary.fitSigma);
      std::cout << "[RUN " << run << "] beta mean=" << summary.fitMean
                << ", sigma=" << summary.fitSigma
                << ", fit entries=" << summary.fitEntries << '\n';
    } else {
      summary.status = summary.selectedEvents < 50 ? "LOW_STATISTICS"
                                                   : "FIT_FAILED";
      PrintFitFailure(run, summary.fitVariable, summary.fitEntries,
                      summary.fitDiagnostics);
    }
  } else {
    DrawBetaVsXfp(tree, "coin", "hms", phase, run, pdfPath);
    DrawBetaVsXfp(tree, "coin", "shms", phase, run, pdfPath);
    DrawCoinTimeVsHmsXfp(tree, phase, run, pdfPath);
    DrawCoinTimeVsShmsXfp(tree, phase, run, pdfPath);
    DrawCoinTime1D(tree, phase, run, pdfPath);

    if (ComputeCoinTimeMetrics(tree, run, summary.fitMean,
                               summary.fitSigma, summary.fitEntries,
                               summary.fitDiagnostics)) {
      summary.status = "OK";
      trendRuns.push_back(run);
      means.push_back(summary.fitMean);
      sigmas.push_back(summary.fitSigma);
      std::cout << "[RUN " << run << "] CTime mean=" << summary.fitMean
                << " ns, sigma=" << summary.fitSigma
                << " ns, fit entries=" << summary.fitEntries << '\n';
    } else {
      summary.status = summary.selectedEvents < 50 ? "LOW_STATISTICS"
                                                   : "FIT_FAILED";
      PrintFitFailure(run, summary.fitVariable, summary.fitEntries,
                      summary.fitDiagnostics);
    }
  }

  file->Close();
  delete file;
  return summary;
}

std::string CsvEscape(const TString &value) {
  std::string text = value.Data();
  if (text.find_first_of(",\"\r\n") == std::string::npos)
    return text;
  std::string escaped = "\"";
  for (char ch : text) {
    if (ch == '"')
      escaped += '"';
    escaped += ch;
  }
  return escaped + '"';
}

bool WriteCsv(const TString &path, const std::vector<RunSummary> &summaries) {
  std::ofstream out(path.Data());
  if (!out) {
    std::cerr << "[ERROR] Cannot create CSV: " << path << '\n';
    return false;
  }
  out << "run,spec,run_type,target,hms_p,shms_p,fit_variable,"
      << "all_events,selected_events,"
      << "fit_mean,fit_sigma,fit_entries,fit_attempted,fit_status,fit_valid,"
      << "covariance_status,fit_chi2,fit_ndf,fit_chi2_ndf,fit_edm,"
      << "candidate_mean,candidate_sigma,failure_reason,status\n";
  for (const RunSummary &row : summaries) {
    out << row.metadata.run << ','
        << row.spec << ','
        << CsvEscape(row.metadata.runType) << ','
        << CsvEscape(row.metadata.target) << ','
        << row.metadata.hmsP << ','
        << row.metadata.shmsP << ','
        << row.fitVariable << ','
        << row.allEvents << ','
        << row.selectedEvents << ',';
    if (std::isfinite(row.fitMean))
      out << row.fitMean;
    out << ',';
    if (std::isfinite(row.fitSigma))
      out << row.fitSigma;
    out << ',';
    if (std::isfinite(row.fitEntries))
      out << row.fitEntries;
    out << ',';
    out << (row.fitDiagnostics.fitAttempted ? 1 : 0) << ',';
    if (row.fitDiagnostics.fitAttempted)
      out << row.fitDiagnostics.minimizerStatus;
    out << ',';
    if (row.fitDiagnostics.fitAttempted || row.fitDiagnostics.fitValid >= 0)
      out << row.fitDiagnostics.fitValid;
    out << ',';
    if (row.fitDiagnostics.fitAttempted)
      out << row.fitDiagnostics.covarianceStatus;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.chi2))
      out << row.fitDiagnostics.chi2;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.ndf))
      out << row.fitDiagnostics.ndf;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.chi2Ndf))
      out << row.fitDiagnostics.chi2Ndf;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.edm))
      out << row.fitDiagnostics.edm;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.candidateMean))
      out << row.fitDiagnostics.candidateMean;
    out << ',';
    if (std::isfinite(row.fitDiagnostics.candidateSigma))
      out << row.fitDiagnostics.candidateSigma;
    out << ',' << CsvEscape(row.fitDiagnostics.failureReason)
        << ',' << row.status << '\n';
  }
  out.close();
  if (!out) {
    std::cerr << "[ERROR] Failed while writing CSV: " << path << '\n';
    return false;
  }
  return true;
}

TString TemporaryPath(const TString &path) {
  const Ssiz_t dot = path.Last('.');
  if (dot == kNPOS)
    return path + ".tmp";
  TString result = path;
  result.Insert(dot, ".tmp");
  return result;
}

bool PublishFile(const TString &temporary, const TString &finalPath) {
  if (gSystem->Rename(temporary, finalPath) != 0) {
    std::cerr << "[ERROR] Cannot publish " << temporary << " as "
              << finalPath << ".\n";
    return false;
  }
  return true;
}

bool ProcessCategory(const TString &spec, const TString &rootDir,
                     int phase, const TString &outputSuffix,
                     const std::vector<RunMetadata> &runs) {
  const TString finalPdf = PdfPath(spec, phase, outputSuffix);
  const TString finalCsv = CsvPath(spec, phase, outputSuffix);
  const TString temporaryPdf = TemporaryPath(finalPdf);
  const TString temporaryCsv = TemporaryPath(finalCsv);
  gSystem->Unlink(temporaryPdf);
  gSystem->Unlink(temporaryCsv);

  std::cout << "\n[INFO] Starting " << spec << " calibration check for "
            << runs.size() << " runs.\n"
            << "[INFO] PDF output: " << finalPdf << '\n'
            << "[INFO] CSV output: " << finalCsv << '\n';

  std::vector<int> trendRuns;
  std::vector<double> means;
  std::vector<double> sigmas;
  std::vector<RunSummary> summaries;
  trendRuns.reserve(runs.size());
  means.reserve(runs.size());
  sigmas.reserve(runs.size());
  summaries.reserve(runs.size());

  OpenPdf(temporaryPdf);
  int processed = 0;
  int ok = 0;
  for (const RunMetadata &run : runs) {
    RunSummary summary = ProcessOneRun(spec, rootDir, phase, run, temporaryPdf,
                                       trendRuns, means, sigmas);
    if (summary.status != "MISSING_FILE" &&
        summary.status != "ZOMBIE_FILE" &&
        summary.status != "MISSING_TREE" &&
        summary.status != "MISSING_BRANCH")
      ++processed;
    if (summary.status == "OK")
      ++ok;
    summaries.push_back(summary);
  }

  DrawDualTrend(trendRuns, means, sigmas, spec, phase, spec == "coin",
                temporaryPdf);
  ClosePdf(temporaryPdf);
  if (!WriteCsv(temporaryCsv, summaries))
    return false;
  if (!PublishFile(temporaryPdf, finalPdf) ||
      !PublishFile(temporaryCsv, finalCsv))
    return false;

  std::cout << "[SUMMARY] " << spec << ": successfully opened/processed "
            << processed << " of " << runs.size() << " requested runs.\n"
            << "[SUMMARY] " << spec << ": OK fits " << ok << " of "
            << runs.size() << " requested runs.\n";
  return true;
}

} // namespace

void hodo_calib_qc_batch(
    int Phase = 2, const char *BigtablePath = "", const char *RootDir = "",
    const char *OutputSuffix = "", bool AvailableFilesOnly = false,
    bool SelectionOnly = false, const char *RunsList = "") {
  gROOT->SetBatch(kTRUE);

  if (Phase != 1 && Phase != 2) {
    std::cerr << "[ERROR] Phase must be 1 or 2.\n";
    return;
  }
  TString bigtablePath = BigtablePath ? BigtablePath : "";
  if (bigtablePath.IsNull())
    bigtablePath = Phase == 1 ? kDefaultBigtablePhase1
                              : kDefaultBigtablePhase2;
  TString rootDir = RootDir ? RootDir : "";
  if (rootDir.IsNull())
    rootDir = Phase == 1 ? kDefaultRootDirPhase1 : kDefaultRootDirPhase2;
  const TString outputSuffix = OutputSuffix ? OutputSuffix : "";
  if (!ValidateOutputSuffix(outputSuffix))
    return;

  RunGroups groups;
  if (!ReadBigtable(bigtablePath, groups))
    return;
  const std::size_t fullCoin = groups.coin.size();
  const std::size_t fullHms = groups.hms.size();
  const std::size_t fullShms = groups.shms.size();
  const TString runsList = RunsList ? RunsList : "";
  if (!runsList.IsNull() && !ApplyRunFilter(runsList, groups))
    return;
  if (groups.coin.empty() || groups.hms.empty() || groups.shms.empty()) {
    std::cerr << "[ERROR] Run selection produced an empty required "
              << "category; no outputs were changed.\n";
    return;
  }

  const std::size_t bigtableCoin = groups.coin.size();
  const std::size_t bigtableHms = groups.hms.size();
  const std::size_t bigtableShms = groups.shms.size();

  std::cout << "[INFO] Bigtable: " << bigtablePath << '\n'
            << "[INFO] ROOT directory: " << rootDir << '\n'
            << "[INFO] Output suffix: '" << outputSuffix << "'\n"
            << "[SELECTION] Full bigtable COIN (PI+SIDIS/PI-SIDIS, hms_p < 0): "
            << fullCoin << '\n'
            << "[SELECTION] Full bigtable HMS (HMSDIS, hms_p < 0): "
            << fullHms << '\n'
            << "[SELECTION] Full bigtable SHMS (SHMSDIS, hms_p < 0): "
            << fullShms << '\n'
            << "[SELECTION] Excluded selected-type rows with hms_p >= 0: "
            << groups.excludedPolarity << '\n'
            << "[SELECTION] Excluded other run types: "
            << groups.excludedRunType << '\n';
  if (!runsList.IsNull())
    std::cout << "[SELECTION] Explicit sample COIN: " << groups.coin.size()
              << '\n'
              << "[SELECTION] Explicit sample HMS: " << groups.hms.size()
              << '\n'
              << "[SELECTION] Explicit sample SHMS: " << groups.shms.size()
              << '\n';

  if (AvailableFilesOnly) {
    groups.coin = KeepAvailableRuns("coin", rootDir, groups.coin);
    groups.hms = KeepAvailableRuns("hms", rootDir, groups.hms);
    groups.shms = KeepAvailableRuns("shms", rootDir, groups.shms);
    std::cout << "[SELECTION] Available-file COIN: " << groups.coin.size()
              << " (omitted " << bigtableCoin - groups.coin.size() << ")\n"
              << "[SELECTION] Available-file HMS: " << groups.hms.size()
              << " (omitted " << bigtableHms - groups.hms.size() << ")\n"
              << "[SELECTION] Available-file SHMS: " << groups.shms.size()
              << " (omitted " << bigtableShms - groups.shms.size() << ")\n";
    if (groups.coin.empty() || groups.hms.empty() || groups.shms.empty()) {
      std::cerr << "[ERROR] Available-file selection produced an empty "
                << "required category; no outputs were changed.\n";
      return;
    }
  }

  if (SelectionOnly) {
    std::cout << "[INFO] Selection-only mode complete; no outputs were "
              << "changed.\n";
    return;
  }

  const TString resultsDir = TString::Format("results/Phase%d", Phase);
  gSystem->mkdir(resultsDir, true);
  gSystem->mkdir(resultsDir + "/pdfs", true);
  gSystem->mkdir(resultsDir + "/tables", true);
  PrintPhysicsLogic(Phase);

  if (!ProcessCategory("shms", rootDir, Phase, outputSuffix, groups.shms))
    return;
  if (!ProcessCategory("hms", rootDir, Phase, outputSuffix, groups.hms))
    return;
  if (!ProcessCategory("coin", rootDir, Phase, outputSuffix, groups.coin))
    return;
}
