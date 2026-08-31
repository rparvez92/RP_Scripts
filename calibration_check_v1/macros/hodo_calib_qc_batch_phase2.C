// hodo_calib_qc_batch_phase2.C
//
// Phase-2 hodoscope calibration QC using unskimmed/main replay ROOT files.
//
// Supported modes:
//   hms           : HMS electron selection, fit H.gtr.beta
//   shms_electron : SHMS electron selection, fit P.gtr.beta
//   coin          : HMS electron + SHMS pion selection, fit CTime.ePiCoinTime_ROC2
//
// Input filename conventions:
//   hms           : hms_coin_replay_production_<run>_-1.root
//   shms_*        : shms_coin_replay_production_<run>_-1.root
//   coin          : coin_replay_production_<run>_-1.root
//
// Standard output:
//   results/Phase2/tables/hodo_qc_<spec>_summary.csv
//   results/Phase2/pdfs/hodo_qc_<spec>_by_run.pdf
//
// Compile only:
//   root -l -b -q -e '.L macros/hodo_calib_qc_batch_phase2.C+'
//
// Process the complete Phase-2 bigtable (SHMSDIS, HMSDIS, then COIN):
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch_phase2.C+()'

#include <TCanvas.h>
#include <TCut.h>
#include <TF1.h>
#include <TFile.h>
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

const char *kDefaultRootDir =
    "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles";

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
  std::vector<RunMetadata> shmsElectron;
  int excludedPolarity = 0;
  int excludedRunType = 0;
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
  TString status = "NOT_RUN";
};

TString MakeFileName(const TString &spec, int run) {
  if (spec == "hms")
    return TString::Format("hms_coin_replay_production_%d_-1.root", run);
  if (spec == "shms_electron")
    return TString::Format("shms_coin_replay_production_%d_-1.root", run);
  if (spec == "coin")
    return TString::Format("coin_replay_production_%d_-1.root", run);
  return "";
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
      spec = "shms_electron";
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
      groups.shmsElectron.push_back(row);
  }

  const auto byRun = [](const RunMetadata &a, const RunMetadata &b) {
    return a.run < b.run;
  };
  std::sort(groups.coin.begin(), groups.coin.end(), byRun);
  std::sort(groups.hms.begin(), groups.hms.end(), byRun);
  std::sort(groups.shmsElectron.begin(), groups.shmsElectron.end(), byRun);
  return true;
}

TCut HmsElectronCuts() {
  return TCut("(H.gtr.dp>-8) && (H.gtr.dp<8)"
              " && (H.gtr.beta>0) && (H.gtr.beta<1.2)"
              " && (H.cal.etottracknorm>0.7)"
              " && (H.cer.npeSum>2.0)");
}

TCut ShmsElectronCuts() {
  return TCut("(P.gtr.dp>-10) && (P.gtr.dp<22)"
              " && (P.gtr.beta>0) && (P.gtr.beta<1.2)"
              " && (P.cal.etottracknorm>0.7)");
}

TCut ShmsPionCuts() {
  const TCut shmsBase =
      "(P.gtr.dp>-10) && (P.gtr.dp<22)"
      " && (P.gtr.beta>0) && (P.gtr.beta<1.2)"
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
  if (spec == "shms_electron")
    return ShmsElectronCuts();
  if (spec == "coin")
    return HmsElectronCuts() && ShmsPionCuts();
  return TCut("");
}

void PrintPhysicsLogic() {
  std::cout
      << "\n===== PHASE-2 HODOSCOPE QC PHYSICS LOGIC =====\n"
      << "HMS electron selection:\n"
      << "  -8 < H.gtr.dp < 8\n"
      << "  0 < H.gtr.beta < 1.2\n"
      << "  H.cal.etottracknorm > 0.7\n"
      << "  H.cer.npeSum > 2.0\n\n"
      << "SHMS electron selection:\n"
      << "  -10 < P.gtr.dp < 22\n"
      << "  0 < P.gtr.beta < 1.2\n"
      << "  P.cal.etottracknorm > 0.7\n"
      << "  No P.ngcer/P.hgcer electron cuts are applied.\n\n"
      << "SHMS pion selection:\n"
      << "  -10 < P.gtr.dp < 22\n"
      << "  0 < P.gtr.beta < 1.2\n"
      << "  P.cal.etottracknorm < 0.8\n"
      << "  P.gtr.p < 2.7: P.aero.npeSum > 2\n"
      << "  P.gtr.p >= 2.7: P.hgcer.npeSum > 1 AND "
         "P.aero.npeSum > 2\n\n"
      << "COIN selection: HMS electron AND SHMS pion selections above.\n"
      << "No coincidence-time gate is applied.\n"
      << "Beta Gaussian fit: +/-0.03 around peak, bounded to [0.9,1.1].\n"
      << "Coin-time Gaussian fit: +/-2 ns around peak, bounded to [0,100].\n"
      << "Fits require at least 50 selected entries.\n"
      << "Beta lines at 0.95 and 1.05 are visual guides only.\n"
      << "================================================\n\n";
}

std::vector<const char *> RequiredBranches(const TString &spec) {
  const std::vector<const char *> hms = {
      "H.gtr.dp", "H.gtr.beta", "H.cal.etottracknorm",
      "H.cer.npeSum", "H.dc.x_fp"};
  const std::vector<const char *> shmsElectron = {
      "P.gtr.dp", "P.gtr.beta", "P.cal.etottracknorm", "P.dc.x_fp"};
  const std::vector<const char *> shmsPion = {
      "P.gtr.dp", "P.gtr.beta", "P.cal.etottracknorm",
      "P.hgcer.npeSum", "P.aero.npeSum", "P.gtr.p", "P.dc.x_fp"};

  if (spec == "hms")
    return hms;
  if (spec == "shms_electron")
    return shmsElectron;

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

TString PdfPath(const TString &spec) {
  return TString::Format("results/Phase2/pdfs/hodo_qc_%s_by_run.pdf",
                         spec.Data());
}

TString CsvPath(const TString &spec) {
  return TString::Format("results/Phase2/tables/hodo_qc_%s_summary.csv",
                         spec.Data());
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
                   const TString &viewSpec, int run,
                   const TString &pdfPath) {
  const bool hmsView = viewSpec == "hms";
  const TString expression =
      hmsView ? "H.gtr.beta:H.dc.x_fp" : "P.gtr.beta:P.dc.x_fp";
  const TString histName =
      TString::Format("h_beta_xfp_%s_%d", viewSpec.Data(), run);
  const TString title = TString::Format(
      "Phase 2 run %d: beta vs xfp (%s);xfp (cm);beta",
      run, hmsView ? "HMS" : "SHMS");

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

bool ComputeBetaMetrics(TTree *tree, const TString &spec, int run,
                        double &mean, double &sigma, double &entries) {
  mean = sigma = entries = std::nan("");
  const TString variable = spec == "hms" ? "H.gtr.beta" : "P.gtr.beta";
  const TString histName = TString::Format("h_beta_fit_%s_%d", spec.Data(), run);
  TH1D hist(histName, ";beta;Counts", 200, 0.2, 1.2);
  hist.Sumw2();
  tree->Project(histName, variable, BuildCuts(spec));
  entries = hist.GetEntries();
  if (entries < 50)
    return false;

  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  const double fitLow = std::max(0.9, peak - 0.03);
  const double fitHigh = std::min(1.1, peak + 0.03);
  TF1 fit(TString::Format("f_beta_%s_%d", spec.Data(), run),
          "gaus", fitLow, fitHigh);
  if (hist.Fit(&fit, "QNR") != 0)
    return false;
  mean = fit.GetParameter(1);
  sigma = std::abs(fit.GetParameter(2));
  return std::isfinite(mean) && std::isfinite(sigma);
}

bool FitCoinTimePeak(TH1D &hist, TF1 &fit,
                     double fitLow, double peak, double fitHigh) {
  const double maximum = hist.GetMaximum();
  fit.SetParameters(maximum, peak, 0.5);
  fit.SetParLimits(0, 0.0, std::max(1.0, maximum * 10.0));
  fit.SetParLimits(1, fitLow, fitHigh);
  fit.SetParLimits(2, 0.05, 2.0);

  if (hist.Fit(&fit, "QNR") != 0)
    return false;

  const double mean = fit.GetParameter(1);
  const double sigma = std::abs(fit.GetParameter(2));
  return std::isfinite(mean) && std::isfinite(sigma) &&
         mean >= fitLow && mean <= fitHigh &&
         sigma >= 0.05 && sigma <= 2.0;
}

void DrawCoinTime1D(TTree *tree, int run, const TString &pdfPath) {
  const TString histName = TString::Format("h_ctime_%d", run);
  TH1D hist(histName,
            TString::Format(
                "Phase 2 run %d: Coincidence Time (ROC2);"
                "CTime.ePiCoinTime_ROC2 (ns);Counts",
                run),
            400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2", BuildCuts("coin"));

  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  const double fitLow = std::max(0.0, peak - 2.0);
  const double fitHigh = std::min(100.0, peak + 2.0);
  TF1 fit(TString::Format("f_ctime_%d", run), "gaus", fitLow, fitHigh);
  const bool fitValid = FitCoinTimePeak(hist, fit, fitLow, peak, fitHigh);

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
                            double &sigma, double &entries) {
  mean = sigma = entries = std::nan("");
  const TString histName = TString::Format("h_ctime_fit_%d", run);
  TH1D hist(histName, ";CTime.ePiCoinTime_ROC2 (ns);Counts", 400, 0, 100);
  hist.Sumw2();
  tree->Project(histName, "CTime.ePiCoinTime_ROC2", BuildCuts("coin"));
  entries = hist.GetEntries();
  if (entries < 50)
    return false;

  const double peak = hist.GetBinCenter(hist.GetMaximumBin());
  const double fitLow = std::max(0.0, peak - 2.0);
  const double fitHigh = std::min(100.0, peak + 2.0);
  TF1 fit(TString::Format("f_ctime_metric_%d", run),
          "gaus", fitLow, fitHigh);
  if (!FitCoinTimePeak(hist, fit, fitLow, peak, fitHigh))
    return false;
  mean = fit.GetParameter(1);
  sigma = std::abs(fit.GetParameter(2));
  return true;
}

void DrawDualTrend(const std::vector<int> &runs,
                   const std::vector<double> &means,
                   const std::vector<double> &sigmas,
                   const TString &spec, bool coinTime,
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
      coinTime
          ? "Phase 2 COIN: CTime (ROC2) mean / sigma vs run;Run;"
            "CTime mean (ns)"
          : TString::Format(
                "Phase 2 %s: beta mean / sigma vs run;Run;beta mean",
                spec.Data());

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
  sigmaGraph.SetMarkerStyle(3);
  sigmaGraph.SetMarkerSize(1.1);
  meanGraph.Draw("P SAME");
  sigmaGraph.Draw("P SAME");

  TGaxis rightAxis(count, meanLow, count, meanHigh,
                   sigmaLow, sigmaHigh, 510, "+L");
  rightAxis.SetTitle(coinTime ? "CTime sigma (ns)" : "beta sigma");
  rightAxis.SetTitleOffset(1.20);
  rightAxis.SetLabelSize(0.035);
  rightAxis.Draw();

  TLegend legend(0.12, 0.84, 0.24, 0.92);
  legend.AddEntry(&meanGraph, "mean", "p");
  legend.AddEntry(&sigmaGraph, "sigma", "p");
  legend.Draw();

  SaveCanvas(canvas, pdfPath);
}

RunSummary ProcessOneRun(const TString &spec, const TString &rootDir,
                         const RunMetadata &metadata, const TString &pdfPath,
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

  if (spec == "hms" || spec == "shms_electron") {
    const TString viewSpec = spec == "hms" ? "hms" : "shms";
    DrawBetaVsXfp(tree, spec, viewSpec, run, pdfPath);
    if (ComputeBetaMetrics(tree, spec, run, summary.fitMean,
                           summary.fitSigma, summary.fitEntries)) {
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
      std::cerr << "[WARN] Beta fit failed or had fewer than 50 entries for run "
                << run << '\n';
    }
  } else {
    DrawBetaVsXfp(tree, "coin", "hms", run, pdfPath);
    DrawBetaVsXfp(tree, "coin", "shms", run, pdfPath);
    DrawCoinTime1D(tree, run, pdfPath);

    if (ComputeCoinTimeMetrics(tree, run, summary.fitMean,
                               summary.fitSigma, summary.fitEntries)) {
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
      std::cerr
          << "[WARN] Coin-time fit failed or had fewer than 50 entries for run "
          << run << '\n';
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
      << "fit_mean,fit_sigma,fit_entries,status\n";
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
    out << ','
        << row.status << '\n';
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
                     const std::vector<RunMetadata> &runs) {
  const TString finalPdf = PdfPath(spec);
  const TString finalCsv = CsvPath(spec);
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
    RunSummary summary = ProcessOneRun(spec, rootDir, run, temporaryPdf,
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

  DrawDualTrend(trendRuns, means, sigmas, spec, spec == "coin", temporaryPdf);
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

void hodo_calib_qc_batch_phase2(
    const char *BigtablePath = "bigtable/rsidis_bigtable_phase2.csv",
    const char *RootDir =
        "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles") {
  gROOT->SetBatch(kTRUE);

  const TString bigtablePath = BigtablePath ? BigtablePath : "";
  if (bigtablePath.IsNull()) {
    std::cerr << "[ERROR] BigtablePath must not be empty.\n";
    return;
  }
  TString rootDir = RootDir ? RootDir : "";
  if (rootDir.IsNull())
    rootDir = kDefaultRootDir;

  RunGroups groups;
  if (!ReadBigtable(bigtablePath, groups))
    return;
  if (groups.coin.empty() || groups.hms.empty() || groups.shmsElectron.empty()) {
    std::cerr << "[ERROR] Bigtable selection produced an empty required "
              << "category; no outputs were changed.\n";
    return;
  }

  std::cout << "[INFO] Bigtable: " << bigtablePath << '\n'
            << "[INFO] ROOT directory: " << rootDir << '\n'
            << "[SELECTION] COIN (PI+SIDIS/PI-SIDIS, hms_p < 0): "
            << groups.coin.size() << '\n'
            << "[SELECTION] HMS (HMSDIS, hms_p < 0): "
            << groups.hms.size() << '\n'
            << "[SELECTION] SHMS electron (SHMSDIS, hms_p < 0): "
            << groups.shmsElectron.size() << '\n'
            << "[SELECTION] Excluded selected-type rows with hms_p >= 0: "
            << groups.excludedPolarity << '\n'
            << "[SELECTION] Excluded other run types: "
            << groups.excludedRunType << '\n';

  gSystem->mkdir("results/Phase2", true);
  gSystem->mkdir("results/Phase2/pdfs", true);
  gSystem->mkdir("results/Phase2/tables", true);
  PrintPhysicsLogic();

  if (!ProcessCategory("shms_electron", rootDir, groups.shmsElectron))
    return;
  if (!ProcessCategory("hms", rootDir, groups.hms))
    return;
  if (!ProcessCategory("coin", rootDir, groups.coin))
    return;
}
