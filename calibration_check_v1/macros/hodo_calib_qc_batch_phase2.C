// hodo_calib_qc_batch_phase2.C
//
// Phase-2 hodoscope calibration QC using unskimmed/main replay ROOT files.
// The physics cuts and fit windows intentionally match hodo_calib_qc_batch.C.
//
// Input filename conventions:
//   hms  : hms_coin_replay_production_<run>_-1.root
//   shms : shms_coin_replay_production_<run>_-1.root
//   coin : coin_replay_production_<run>_-1.root
//
// Output:
//   results/Phase2/hmsPNGs/
//   results/Phase2/shmsPNGs/
//   results/Phase2/coinPNGs/
//
// Compile only:
//   root -l -b -q -e '.L macros/hodo_calib_qc_batch_phase2.C+'
//
// First COIN test (run only after reviewing the printed physics cuts):
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch_phase2.C+("coin","/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles","27128")'
//
// Current Phase-2 COIN set:
//   root -l -b -q \
//     'macros/hodo_calib_qc_batch_phase2.C+("coin","/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles","27122-27133")'

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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const char *kDefaultRootDir =
    "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles";

// Dedicated HMS/SHMS Phase-2 runs are not available yet.
const std::vector<int> kPhase2HmsRuns;
const std::vector<int> kPhase2ShmsRuns;
const std::vector<int> kPhase2CoinRuns = {
    27122, 27123, 27124, 27125, 27126, 27127,
    27128, 27129, 27130, 27131, 27132, 27133};

TString MakeFileName(const TString &spec, int run) {
  if (spec == "hms")
    return TString::Format("hms_coin_replay_production_%d_-1.root", run);
  if (spec == "shms")
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

TCut BuildCuts(const TString &spec) {
  const TCut hmsPid =
      "(H.gtr.dp>-8) && (H.gtr.dp<8)"
      " && (H.gtr.beta>0) && (H.gtr.beta<1.2)"
      " && (H.cal.etottracknorm>0.7)"
      " && (H.cer.npeSum>2.0)";

  const TCut shmsBase =
      "(P.gtr.dp>-10) && (P.gtr.dp<22)"
      " && (P.gtr.beta>0) && (P.gtr.beta<1.2)"
      " && (P.cal.etottracknorm<0.8)";

  const TCut shmsAero =
      "(P.gtr.p<2.7) && (P.aero.npeSum>2)";
  const TCut shmsHgc =
      "(P.gtr.p>=2.7) && (P.hgcer.npeSum>1)"
      " && (P.aero.npeSum>2)";
  const TCut shmsPid = shmsBase && (shmsAero || shmsHgc);

  if (spec == "hms")
    return hmsPid;
  if (spec == "shms")
    return shmsPid;
  if (spec == "coin")
    return hmsPid && shmsPid;
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
      << "SHMS pion selection:\n"
      << "  -10 < P.gtr.dp < 22\n"
      << "  0 < P.gtr.beta < 1.2\n"
      << "  P.cal.etottracknorm < 0.8\n"
      << "  P.gtr.p < 2.7: P.aero.npeSum > 2\n"
      << "  P.gtr.p >= 2.7: P.hgcer.npeSum > 1 AND "
         "P.aero.npeSum > 2\n\n"
      << "COIN selection: HMS AND SHMS selections above.\n"
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
  const std::vector<const char *> shms = {
      "P.gtr.dp", "P.gtr.beta", "P.cal.etottracknorm",
      "P.ngcer.npeSum", "P.hgcer.npeSum", "P.aero.npeSum",
      "P.gtr.p", "P.dc.x_fp"};

  if (spec == "hms")
    return hms;
  if (spec == "shms")
    return shms;

  std::vector<const char *> coin = hms;
  coin.insert(coin.end(), shms.begin(), shms.end());
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

TString OutputDir(const TString &spec) {
  return TString::Format("results/Phase2/%sPNGs", spec.Data());
}

void DrawBetaVsXfp(TTree *tree, const TString &selectionSpec,
                   const TString &viewSpec, int run) {
  const bool hmsView = viewSpec == "hms";
  const TString expression =
      hmsView ? "H.gtr.beta:H.dc.x_fp" : "P.gtr.beta:P.dc.x_fp";
  const TString histName =
      TString::Format("h_beta_xfp_%s_%d", viewSpec.Data(), run);
  const TString title = TString::Format(
      "Phase 2 run %d: #beta vs x_{fp} (%s);x_{fp} (cm);#beta",
      run, hmsView ? "HMS" : "SHMS");

  TH2D hist(histName, title, 80, -45, 45, 120, 0.2, 1.2);
  hist.Sumw2();
  tree->Project(histName, expression, BuildCuts(selectionSpec));

  TCanvas canvas(TString::Format("c_beta_xfp_%s_%d", viewSpec.Data(), run),
                 "", 900, 700);
  canvas.SetRightMargin(0.12);
  gStyle->SetOptStat(0);
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

  TString output;
  if (selectionSpec == "coin") {
    output = TString::Format("%s/coin_%s_run%d_beta_vs_xfp.png",
                             OutputDir("coin").Data(), viewSpec.Data(), run);
  } else {
    output = TString::Format("%s/%s_run%d_beta_vs_xfp.png",
                             OutputDir(selectionSpec).Data(),
                             selectionSpec.Data(), run);
  }
  canvas.SaveAs(output);
}

bool ComputeBetaMetrics(TTree *tree, const TString &spec, int run,
                        double &mean, double &sigma, double &entries) {
  mean = sigma = entries = std::nan("");
  const TString variable = spec == "shms" ? "P.gtr.beta" : "H.gtr.beta";
  const TString histName = TString::Format("h_beta_fit_%s_%d", spec.Data(), run);
  TH1D hist(histName, ";#beta;Counts", 200, 0.2, 1.2);
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
  sigma = fit.GetParameter(2);
  return true;
}

void DrawCoinTime1D(TTree *tree, int run) {
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
  hist.Fit(&fit, "QNR");

  TCanvas canvas(TString::Format("c_ctime_%d", run), "", 800, 600);
  canvas.SetLeftMargin(0.12);
  gStyle->SetOptStat(0);
  hist.Draw("HIST");
  fit.SetLineColor(kRed);
  fit.SetLineWidth(3);
  fit.Draw("SAME");
  canvas.SaveAs(TString::Format(
      "%s/coin_run%d_ctime1D_ROC2.png", OutputDir("coin").Data(), run));
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
  if (hist.Fit(&fit, "QNR") != 0)
    return false;
  mean = fit.GetParameter(1);
  sigma = fit.GetParameter(2);
  return true;
}

void DrawDualTrend(const std::vector<int> &runs,
                   const std::vector<double> &means,
                   const std::vector<double> &sigmas,
                   const TString &spec, bool coinTime) {
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
                "Phase 2 %s: #beta mean / sigma vs run;Run;#beta mean",
                spec.Data());

  TCanvas canvas(TString::Format("c_trend_%s", spec.Data()), "", 1000, 600);
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
  frame.GetXaxis()->SetLabelSize(0.035);
  frame.GetYaxis()->SetTitleOffset(1.2);
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
  rightAxis.SetTitle(coinTime ? "CTime sigma (ns)" : "#beta sigma");
  rightAxis.SetTitleOffset(1.2);
  rightAxis.Draw();

  TLegend legend(0.12, 0.84, 0.24, 0.92);
  legend.AddEntry(&meanGraph, "mean", "p");
  legend.AddEntry(&sigmaGraph, "sigma", "p");
  legend.Draw();

  const TString output =
      coinTime
          ? TString::Format("%s/coin_ctime_trends_ROC2.png",
                            OutputDir("coin").Data())
          : TString::Format("%s/%s_beta_trends.png",
                            OutputDir(spec).Data(), spec.Data());
  canvas.SaveAs(output);
}

bool ProcessOneRun(const TString &spec, const TString &rootDir, int run,
                   std::vector<int> &trendRuns,
                   std::vector<double> &means,
                   std::vector<double> &sigmas) {
  const TString path =
      TString::Format("%s/%s", rootDir.Data(), MakeFileName(spec, run).Data());
  TFile *file = TFile::Open(path, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "[WARN] Could not open run " << run << ": " << path << '\n';
    if (file)
      delete file;
    return false;
  }

  TTree *tree = dynamic_cast<TTree *>(file->Get("T"));
  if (!tree) {
    std::cerr << "[WARN] Tree 'T' missing for run " << run << ": "
              << path << '\n';
    file->Close();
    delete file;
    return false;
  }
  if (!ValidateAndEnableBranches(tree, spec, run)) {
    file->Close();
    delete file;
    return false;
  }

  const Long64_t allEvents = tree->GetEntries();
  const Long64_t selectedEvents = tree->GetEntries(BuildCuts(spec));
  std::cout << "[RUN " << run << "] all=" << allEvents
            << ", selected=" << selectedEvents << '\n';

  if (spec == "hms" || spec == "shms") {
    DrawBetaVsXfp(tree, spec, spec, run);
    double mean, sigma, entries;
    if (ComputeBetaMetrics(tree, spec, run, mean, sigma, entries)) {
      trendRuns.push_back(run);
      means.push_back(mean);
      sigmas.push_back(sigma);
      std::cout << "[RUN " << run << "] beta mean=" << mean
                << ", sigma=" << sigma
                << ", fit entries=" << entries << '\n';
    } else {
      std::cerr << "[WARN] Beta fit failed or had fewer than 50 entries for run "
                << run << '\n';
    }
  } else {
    DrawBetaVsXfp(tree, "coin", "hms", run);
    DrawBetaVsXfp(tree, "coin", "shms", run);
    DrawCoinTime1D(tree, run);

    double mean, sigma, entries;
    if (ComputeCoinTimeMetrics(tree, run, mean, sigma, entries)) {
      trendRuns.push_back(run);
      means.push_back(mean);
      sigmas.push_back(sigma);
      std::cout << "[RUN " << run << "] CTime mean=" << mean
                << " ns, sigma=" << sigma
                << " ns, fit entries=" << entries << '\n';
    } else {
      std::cerr
          << "[WARN] Coin-time fit failed or had fewer than 50 entries for run "
          << run << '\n';
    }
  }

  file->Close();
  delete file;
  return true;
}

} // namespace

void hodo_calib_qc_batch_phase2(const char *Spec = "",
                                const char *RootDir = "",
                                const char *RunsList = "") {
  gROOT->SetBatch(kTRUE);

  const TString spec = Spec ? Spec : "";
  if (!(spec == "hms" || spec == "shms" || spec == "coin")) {
    std::cerr << "[ERROR] Spec must be 'hms', 'shms', or 'coin'\n";
    return;
  }

  TString rootDir = RootDir ? RootDir : "";
  if (rootDir.IsNull())
    rootDir = kDefaultRootDir;

  std::vector<int> runs = ParseRunsList(RunsList ? RunsList : "");
  if (runs.empty()) {
    if (spec == "hms")
      runs = kPhase2HmsRuns;
    else if (spec == "shms")
      runs = kPhase2ShmsRuns;
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
            << "[INFO] Output directory: " << OutputDir(spec) << '\n';

  gSystem->mkdir("results/Phase2", true);
  gSystem->mkdir(OutputDir(spec), true);

  std::vector<int> trendRuns;
  std::vector<double> means;
  std::vector<double> sigmas;
  trendRuns.reserve(runs.size());
  means.reserve(runs.size());
  sigmas.reserve(runs.size());

  int processed = 0;
  for (int run : runs) {
    if (ProcessOneRun(spec, rootDir, run, trendRuns, means, sigmas))
      ++processed;
  }

  if (spec == "coin")
    DrawDualTrend(trendRuns, means, sigmas, spec, true);
  else
    DrawDualTrend(trendRuns, means, sigmas, spec, false);

  std::cout << "[SUMMARY] Successfully opened and processed " << processed
            << " of " << runs.size() << " requested runs.\n";
}
