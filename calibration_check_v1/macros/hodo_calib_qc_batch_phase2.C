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
// Optional PNG output:
//   results/Phase2/<spec>PNGs/
//
// Compile only:
//   root -l -b -q -e '.L macros/hodo_calib_qc_batch_phase2.C+'
//
// Standard compact output:
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch_phase2.C+("coin","/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles","27122-27132")'
//
// Compact output plus PNGs:
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch_phase2.C+("coin","/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles","27122-27132","pdf,csv,png")'

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
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const char *kDefaultRootDir =
    "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles";

// These defaults are intentionally conservative. During live Phase-2 running,
// prefer passing explicit run lists generated from rsidis_runlist_phaseII.dat.
const std::vector<int> kPhase2HmsRuns;
const std::vector<int> kPhase2ShmsElectronRuns;
const std::vector<int> kPhase2CoinRuns = {
    27122, 27123, 27124, 27125, 27126, 27127,
    27128, 27129, 27130, 27131, 27132};

const int kCanvasWidth = 1250;
const int kCanvasHeight = 850;

struct OutputOptions {
  bool pdf = true;
  bool csv = true;
  bool png = false;
};

struct RunSummary {
  int run = 0;
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

TString NormalizeSpec(const TString &rawSpec) {
  TString spec = rawSpec;
  spec.ToLower();
  spec.ReplaceAll("-", "_");
  return spec;
}

bool IsValidSpec(const TString &spec) {
  return spec == "hms" || spec == "shms_electron" || spec == "coin";
}

TString MakeFileName(const TString &spec, int run) {
  if (spec == "hms")
    return TString::Format("hms_coin_replay_production_%d_-1.root", run);
  if (spec == "shms_electron")
    return TString::Format("shms_coin_replay_production_%d_-1.root", run);
  if (spec == "coin")
    return TString::Format("coin_replay_production_%d_-1.root", run);
  return "";
}

std::vector<int> ParseRunsList(const std::string &text) {
  std::vector<int> runs;
  if (text.empty())
    return runs;

  std::stringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
                token.end());
    if (token.empty())
      continue;

    const std::string::size_type dash = token.find('-');
    if (dash == std::string::npos) {
      runs.push_back(std::stoi(token));
      continue;
    }

    int first = std::stoi(token.substr(0, dash));
    int last = std::stoi(token.substr(dash + 1));
    if (last < first)
      std::swap(first, last);
    for (int run = first; run <= last; ++run)
      runs.push_back(run);
  }
  return runs;
}

OutputOptions ParseOutputOptions(const TString &modeText) {
  OutputOptions options;
  options.pdf = false;
  options.csv = false;
  options.png = false;

  TString text = modeText;
  text.ToLower();
  text.ReplaceAll(" ", "");
  if (text.IsNull())
    text = "pdf,csv";

  std::stringstream stream(text.Data());
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (token == "pdf")
      options.pdf = true;
    else if (token == "csv")
      options.csv = true;
    else if (token == "png")
      options.png = true;
    else if (token == "all") {
      options.pdf = true;
      options.csv = true;
      options.png = true;
    } else if (!token.empty()) {
      std::cerr << "[WARN] Unknown output token '" << token
                << "' ignored. Valid tokens: pdf,csv,png,all\n";
    }
  }

  if (!options.pdf && !options.csv && !options.png) {
    std::cerr << "[WARN] No valid output modes requested; using pdf,csv.\n";
    options.pdf = true;
    options.csv = true;
  }
  return options;
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

TString PngDir(const TString &spec) {
  return TString::Format("results/Phase2/%sPNGs", spec.Data());
}

TString PdfPath(const TString &spec) {
  return TString::Format("results/Phase2/pdfs/hodo_qc_%s_by_run.pdf",
                         spec.Data());
}

TString CsvPath(const TString &spec) {
  return TString::Format("results/Phase2/tables/hodo_qc_%s_summary.csv",
                         spec.Data());
}

void SaveCanvas(TCanvas &canvas, const TString &pdfPath,
                const TString &pngPath, const OutputOptions &options) {
  if (options.pdf)
    canvas.Print(pdfPath);
  if (options.png)
    canvas.SaveAs(pngPath);
}

void OpenPdf(const TString &pdfPath, const OutputOptions &options) {
  if (options.pdf) {
    TCanvas opener("c_pdf_open", "", 1, 1);
    opener.Print(TString::Format("%s[", pdfPath.Data()));
  }
}

void ClosePdf(const TString &pdfPath, const OutputOptions &options) {
  if (options.pdf) {
    TCanvas closer("c_pdf_close", "", 1, 1);
    closer.Print(TString::Format("%s]", pdfPath.Data()));
  }
}

void DrawBetaVsXfp(TTree *tree, const TString &selectionSpec,
                   const TString &viewSpec, int run,
                   const TString &pdfPath, const OutputOptions &options) {
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

  TString pngPath;
  if (selectionSpec == "coin") {
    pngPath = TString::Format("%s/coin_%s_run%d_beta_vs_xfp.png",
                              PngDir("coin").Data(), viewSpec.Data(), run);
  } else {
    pngPath = TString::Format("%s/%s_run%d_beta_vs_xfp.png",
                              PngDir(selectionSpec).Data(),
                              selectionSpec.Data(), run);
  }
  SaveCanvas(canvas, pdfPath, pngPath, options);
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

void DrawCoinTime1D(TTree *tree, int run, const TString &pdfPath,
                    const OutputOptions &options) {
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

  const TString pngPath = TString::Format(
      "%s/coin_run%d_ctime1D_ROC2.png", PngDir("coin").Data(), run);
  SaveCanvas(canvas, pdfPath, pngPath, options);
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
                   const TString &pdfPath, const OutputOptions &options) {
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
  const int labelStep = count > 20 ? 5 : 1;
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

  const TString pngPath =
      coinTime
          ? TString::Format("%s/coin_ctime_trends_ROC2.png",
                            PngDir("coin").Data())
          : TString::Format("%s/%s_beta_trends.png",
                            PngDir(spec).Data(), spec.Data());
  SaveCanvas(canvas, pdfPath, pngPath, options);
}

RunSummary ProcessOneRun(const TString &spec, const TString &rootDir, int run,
                         const TString &pdfPath,
                         const OutputOptions &options,
                         std::vector<int> &trendRuns,
                         std::vector<double> &means,
                         std::vector<double> &sigmas) {
  RunSummary summary;
  summary.run = run;
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
    DrawBetaVsXfp(tree, spec, viewSpec, run, pdfPath, options);
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
    DrawBetaVsXfp(tree, "coin", "hms", run, pdfPath, options);
    DrawBetaVsXfp(tree, "coin", "shms", run, pdfPath, options);
    DrawCoinTime1D(tree, run, pdfPath, options);

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

void WriteCsv(const TString &path, const std::vector<RunSummary> &summaries) {
  std::ofstream out(path.Data());
  out << "run,spec,fit_variable,all_events,selected_events,"
      << "fit_mean,fit_sigma,fit_entries,status\n";
  for (const RunSummary &row : summaries) {
    out << row.run << ','
        << row.spec << ','
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
}

} // namespace

void hodo_calib_qc_batch_phase2(const char *Spec = "",
                                const char *RootDir = "",
                                const char *RunsList = "",
                                const char *OutputMode = "pdf,csv") {
  gROOT->SetBatch(kTRUE);

  const TString spec = NormalizeSpec(Spec ? Spec : "");
  if (!IsValidSpec(spec)) {
    std::cerr << "[ERROR] Spec must be 'hms', 'shms_electron', "
              << "or 'coin'.\n";
    return;
  }

  TString rootDir = RootDir ? RootDir : "";
  if (rootDir.IsNull())
    rootDir = kDefaultRootDir;

  const OutputOptions options =
      ParseOutputOptions(OutputMode ? OutputMode : "pdf,csv");

  std::vector<int> runs = ParseRunsList(RunsList ? RunsList : "");
  if (runs.empty()) {
    if (spec == "hms")
      runs = kPhase2HmsRuns;
    else if (spec == "shms_electron")
      runs = kPhase2ShmsElectronRuns;
    else
      runs = kPhase2CoinRuns;
  }
  if (runs.empty()) {
    std::cerr << "[ERROR] No Phase-2 " << spec
              << " runs are configured. Supply RunsList explicitly.\n";
    return;
  }

  PrintPhysicsLogic();
  std::cout << "[INFO] Spec: " << spec << '\n'
            << "[INFO] ROOT directory: " << rootDir << '\n'
            << "[INFO] Requested runs: " << runs.size() << '\n'
            << "[INFO] Output modes: "
            << (options.pdf ? "pdf " : "")
            << (options.csv ? "csv " : "")
            << (options.png ? "png " : "") << '\n';

  gSystem->mkdir("results/Phase2", true);
  gSystem->mkdir("results/Phase2/pdfs", true);
  gSystem->mkdir("results/Phase2/tables", true);
  if (options.png)
    gSystem->mkdir(PngDir(spec), true);

  const TString pdfPath = PdfPath(spec);
  const TString csvPath = CsvPath(spec);
  if (options.pdf)
    std::cout << "[INFO] PDF output: " << pdfPath << '\n';
  if (options.csv)
    std::cout << "[INFO] CSV output: " << csvPath << '\n';
  if (options.png)
    std::cout << "[INFO] PNG output directory: " << PngDir(spec) << '\n';

  std::vector<int> trendRuns;
  std::vector<double> means;
  std::vector<double> sigmas;
  std::vector<RunSummary> summaries;
  trendRuns.reserve(runs.size());
  means.reserve(runs.size());
  sigmas.reserve(runs.size());
  summaries.reserve(runs.size());

  OpenPdf(pdfPath, options);
  int processed = 0;
  int ok = 0;
  for (int run : runs) {
    RunSummary summary = ProcessOneRun(spec, rootDir, run, pdfPath, options,
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

  if (spec == "coin")
    DrawDualTrend(trendRuns, means, sigmas, spec, true, pdfPath, options);
  else
    DrawDualTrend(trendRuns, means, sigmas, spec, false, pdfPath, options);
  ClosePdf(pdfPath, options);

  if (options.csv)
    WriteCsv(csvPath, summaries);

  std::cout << "[SUMMARY] Successfully opened/processed " << processed
            << " of " << runs.size() << " requested runs.\n"
            << "[SUMMARY] OK fits: " << ok << " of " << runs.size()
            << " requested runs.\n";
}
