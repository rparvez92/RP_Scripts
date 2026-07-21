// File: rate_dependence_v1/macros/TauVsRateAndCurrent.C
//
// UPDATED (axis fix):
// - For each pad, x-axis limits are computed from BOTH pi+ and pi- points,
//   with a small margin, so no series gets clipped.
// - Applied to BOTH tau_vs_shms34 and tau_vs_current.
//
// Filters remain:
//   (1) N_valid_for_fit2 >= 6   (implemented as nFitPoints from use_in_fit==1 rows)
//   (2) tau_ns_err / |tau_ns| < 1.0   (requires tau_ns_err present & finite)
//
// Run from rate_dependence_v1/:
//   root -l -b -q 'macros/TauVsRateAndCurrent.C()'
// or:
//   root -l -b -q 'macros/TauVsRateAndCurrent.C("results")'
//
// Produces:
//   <resultsDir>/tau/tau_vs_shms34.png
//   <resultsDir>/tau/tau_vs_current.png

#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>

namespace {

  // ---------- Small helpers ----------
  static inline std::string Trim(const std::string& s) {
    const char* ws = " \t\r\n";
    const auto b = s.find_first_not_of(ws);
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
  }

  static inline std::vector<std::string> SplitCSV(const std::string& line) {
    std::vector<std::string> out;
    std::string token;
    std::stringstream ss(line);
    while (std::getline(ss, token, ',')) out.push_back(Trim(token));
    return out;
  }

  static inline bool ToInt(const std::string& s, int& v) {
    try {
      size_t idx = 0;
      long x = std::stol(s, &idx);
      if (idx != s.size()) return false;
      v = (int)x;
      return true;
    } catch (...) { return false; }
  }

  static inline bool ToDouble(const std::string& s, double& v) {
    try {
      size_t idx = 0;
      double x = std::stod(s, &idx);
      if (idx != s.size()) return false;
      v = x;
      return true;
    } catch (...) { return false; }
  }

  static inline double Mean(const std::vector<double>& a) {
    if (a.empty()) return std::nan("");
    double s = 0;
    for (double x : a) s += x;
    return s / (double)a.size();
  }

  static inline double RMS(const std::vector<double>& a, double mean) {
    if (a.size() < 2) return 0.0;
    double s2 = 0;
    for (double x : a) {
      const double d = x - mean;
      s2 += d * d;
    }
    return std::sqrt(s2 / (double)a.size());
  }

  // Compute padded x-range across two buckets.
  // Returns false if no finite points exist.
  bool ComputeXRange(const std::vector<double>& x1, const std::vector<double>& x2,
                     double& xmin, double& xmax)
  {
    bool have = false;
    double mn = 1e300, mx = -1e300;

    auto scan = [&](const std::vector<double>& xs) {
      for (double v : xs) {
        if (!std::isfinite(v)) continue;
        mn = std::min(mn, v);
        mx = std::max(mx, v);
        have = true;
      }
    };
    scan(x1);
    scan(x2);

    if (!have) return false;

    // Add margin
    double span = mx - mn;
    if (!(span > 0) || !std::isfinite(span)) span = std::max(1.0, std::abs(mx) * 0.05);
    const double margin = 0.08 * span; // 8% padding
    xmin = mn - margin;
    xmax = mx + margin;

    // If xmin==xmax after padding, expand slightly
    if (!(xmax > xmin)) {
      xmin = mn - 1.0;
      xmax = mx + 1.0;
    }
    return true;
  }

  // ---------- Recursively collect files ----------
  void CollectFilesRec(const TString& dirPath,
                       const TString& targetFileName,
                       std::vector<TString>& outFiles)
  {
    TSystemDirectory dir(dirPath, dirPath);
    TList* files = dir.GetListOfFiles();
    if (!files) return;

    TIter next(files);
    while (TSystemFile* f = (TSystemFile*) next()) {
      TString name = f->GetName();
      if (name == "." || name == "..") continue;

      TString fullPath = dirPath + "/" + name;

      if (f->IsDirectory()) {
        CollectFilesRec(fullPath, targetFileName, outFiles);
      } else {
        if (name == targetFileName) outFiles.push_back(fullPath);
      }
    }
  }

  // ---------- Parse tau_ns and tau_ns_err from header ----------
  bool ParseTauHeader(const TString& csvPath, double& tauNs, double& tauErrNs, bool& hasTauErr)
  {
    std::ifstream in(csvPath.Data());
    if (!in.is_open()) return false;

    bool foundTau = false;
    bool foundErr = false;
    double tau = std::nan("");
    double terr = std::nan("");

    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] != '#') break;

      auto posTau = line.find("tau_ns:");
      if (posTau != std::string::npos) {
        std::string rhs = Trim(line.substr(posTau + std::string("tau_ns:").size()));
        double v;
        if (ToDouble(rhs, v)) { tau = v; foundTau = true; }
      }

      auto posErr = line.find("tau_ns_err:");
      if (posErr != std::string::npos) {
        std::string rhs = Trim(line.substr(posErr + std::string("tau_ns_err:").size()));
        double v;
        if (ToDouble(rhs, v)) { terr = v; foundErr = true; }
      }
    }

    if (!foundTau) return false;
    if (!std::isfinite(tau)) return false;
    if (tau <= 0.0) return false;
    if (tau > 1e6) return false;

    tauNs = tau;

    hasTauErr = (foundErr && std::isfinite(terr) && terr > 0.0 && terr < 1e6);
    tauErrNs  = hasTauErr ? terr : 0.0;

    return true;
  }

  // ---------- Extract pass & target ----------
  bool ExtractPassAndTarget(const TString& path, int& pass, TString& target)
  {
    pass = 0; target = "";
    if (path.Contains("/pass4/")) pass = 4;
    else if (path.Contains("/pass5/")) pass = 5;

    if (path.Contains("/LH2/")) target = "LH2";
    else if (path.Contains("/LD2/")) target = "LD2";

    return (pass != 0 && (target == "LH2" || target == "LD2"));
  }

  // ---------- Extract pion sign ----------
  TString ExtractPionSign(const TString& path)
  {
    if (path.Contains("/pi+")) return "pi+";
    if (path.Contains("/pi-")) return "pi-";
    return "unk";
  }

  // ---------- Pad index ----------
  int PadIndex(int pass, const TString& target)
  {
    if (pass == 4 && target == "LH2") return 0;
    if (pass == 4 && target == "LD2") return 1;
    if (pass == 5 && target == "LH2") return 2;
    if (pass == 5 && target == "LD2") return 3;
    return -1;
  }

  // ---------- Read yield_ratio CSV to get fit-runs and fit-rate stats ----------
  bool ReadFitRunsAndRate(const TString& yieldRatioCsv,
                          std::vector<int>& fitRuns,
                          double& meanRate,
                          double& rmsRate,
                          int& nFitPoints)
  {
    std::ifstream in(yieldRatioCsv.Data());
    if (!in.is_open()) return false;

    std::string line;

    while (std::getline(in, line)) {
      if (!line.empty() && line[0] != '#') break;
    }
    if (!in.good()) return false;

    const auto hdr = SplitCSV(line);

    auto findCol = [&](const std::string& name) -> int {
      for (size_t i = 0; i < hdr.size(); i++) {
        if (hdr[i] == name) return (int)i;
      }
      return -1;
    };

    const int iRun   = findCol("run");
    const int iRate  = findCol("shms34_rate_kHz");
    const int iUse   = findCol("use_in_fit");
    if (iRun < 0 || iRate < 0 || iUse < 0) return false;

    std::vector<double> rates;
    std::vector<int> runs;

    while (std::getline(in, line)) {
      if (line.empty()) continue;
      const auto row = SplitCSV(line);
      const int ncol = (int)row.size();
      if (iRun >= ncol || iRate >= ncol || iUse >= ncol) continue;

      int use = 0;
      if (!ToInt(row[iUse], use)) continue;
      if (use != 1) continue;

      int run = 0;
      double rate = 0.0;
      if (!ToInt(row[iRun], run)) continue;
      if (!ToDouble(row[iRate], rate)) continue;
      if (!std::isfinite(rate)) continue;

      runs.push_back(run);
      rates.push_back(rate);
    }

    nFitPoints = (int)rates.size();
    if (nFitPoints < 2) return false;

    std::sort(runs.begin(), runs.end());
    runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
    fitRuns = runs;

    meanRate = Mean(rates);
    rmsRate  = RMS(rates, meanRate);

    return std::isfinite(meanRate);
  }

  // ---------- Read current CSV into map(run -> BCM2_I) ----------
  bool ReadRunToCurrentMap(const TString& currentCsv,
                           std::unordered_map<int,double>& run2I)
  {
    std::ifstream in(currentCsv.Data());
    if (!in.is_open()) return false;

    std::string line;
    if (!std::getline(in, line)) return false;
    const auto hdr = SplitCSV(line);

    auto findCol = [&](const std::string& name) -> int {
      for (size_t i = 0; i < hdr.size(); i++) {
        if (hdr[i] == name) return (int)i;
      }
      return -1;
    };

    const int iRun = findCol("run");
    const int iI   = findCol("BCM2_I");
    if (iRun < 0 || iI < 0) return false;

    while (std::getline(in, line)) {
      if (line.empty()) continue;
      const auto row = SplitCSV(line);
      const int ncol = (int)row.size();
      if (iRun >= ncol || iI >= ncol) continue;

      int run = 0;
      double I = 0.0;
      if (!ToInt(row[iRun], run)) continue;
      if (!ToDouble(row[iI], I)) continue;
      if (!std::isfinite(I)) continue;

      run2I[run] = I;
    }
    return !run2I.empty();
  }

  // ---------- Compute mean/rms current over fitRuns ----------
  bool ComputeMeanCurrentOverRuns(const std::vector<int>& fitRuns,
                                  const std::unordered_map<int,double>& run2I,
                                  double& meanI,
                                  double& rmsI,
                                  int& nMatched)
  {
    std::vector<double> currents;
    for (int r : fitRuns) {
      auto it = run2I.find(r);
      if (it == run2I.end()) continue;
      currents.push_back(it->second);
    }
    nMatched = (int)currents.size();
    if (nMatched < 1) return false;

    meanI = Mean(currents);
    rmsI  = RMS(currents, meanI);
    return std::isfinite(meanI);
  }

  // ---------- Graph bucket ----------
  struct Bucket {
    std::vector<double> x, y, ex, ey;
  };

  // ---------- Style helper ----------
  void StyleGraph(TGraphErrors* g, const TString& sign)
  {
    g->SetLineWidth(2);
    if (sign == "pi+") {
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.0);
      g->SetMarkerColor(kBlue + 1);
      g->SetLineColor(kBlue + 1);
    } else if (sign == "pi-") {
      g->SetMarkerStyle(24);
      g->SetMarkerSize(1.0);
      g->SetMarkerColor(kRed + 1);
      g->SetLineColor(kRed + 1);
    } else {
      g->SetMarkerStyle(21);
      g->SetMarkerSize(1.0);
      g->SetMarkerColor(kBlack);
      g->SetLineColor(kBlack);
    }
  }

} // end anon namespace

// ============================================================================
// Main macro
// ============================================================================
void TauVsRateAndCurrent(const char* resultsDir = "results")
{
  gStyle->SetOptStat(0);

  // ---- Filters ----
  const int    kMinNFit2 = 6;
  const double kMaxRelTauErr = 1.0;

  const TString kResultsDir = resultsDir;
  const TString kYieldRatioName = "yield_ratio_vs_trigger_shms34.csv";
  const TString kCurrentName    = "hms_elclean_chargeNorm_vs_current.csv";

  // Collect all yield_ratio CSV files
  std::vector<TString> yieldRatioFiles;
  CollectFilesRec(kResultsDir, kYieldRatioName, yieldRatioFiles);

  std::cout << "[TauVsRateAndCurrent] Found " << yieldRatioFiles.size()
            << " candidate '" << kYieldRatioName << "' files under '" << kResultsDir << "'.\n";

  // Ensure output directory exists
  const TString outDir = kResultsDir + "/tau";
  if (gSystem->AccessPathName(outDir)) {
    int mk = gSystem->mkdir(outDir, kTRUE);
    if (mk != 0) {
      std::cerr << "[ERROR] Could not create output directory: " << outDir << "\n";
      return;
    }
  }

  const TString outPngRate = outDir + "/tau_vs_shms34.png";
  const TString outPngCurr = outDir + "/tau_vs_current.png";

  // Buckets: pad 0..3 and sign 0=pi+, 1=pi-
  Bucket rateBuckets[4][2];
  Bucket currBuckets[4][2];

  auto signIndex = [](const TString& s) -> int {
    if (s == "pi+") return 0;
    if (s == "pi-") return 1;
    return -1;
  };

  // Counters
  int nProcessed = 0;
  int nSkipBadCat = 0, nSkipSign = 0;
  int nSkipTau = 0, nSkipNoTauErr = 0, nSkipRelErr = 0;
  int nSkipFitPts = 0, nSkipMinFitPts = 0;
  int nSkipNoCurrentFile = 0, nSkipCurrentParse = 0, nSkipNoRunMatch = 0;

  struct LogLine {
    TString shortPath;
    int pass;
    TString target;
    TString sign;
    double tau, tauerr, relerr;
    int nFit;
    double meanRate, rmsRate;
    int nRunMatch;
    double meanI, rmsI;
  };
  std::vector<LogLine> logs;

  for (const auto& yrPath : yieldRatioFiles) {

    int pass = 0;
    TString target;
    if (!ExtractPassAndTarget(yrPath, pass, target)) {
      std::cerr << "[WARN] Cannot categorize pass/target from path, skipping: " << yrPath << "\n";
      nSkipBadCat++;
      continue;
    }
    const int pad = PadIndex(pass, target);
    if (pad < 0) { nSkipBadCat++; continue; }

    const TString sign = ExtractPionSign(yrPath);
    const int sidx = signIndex(sign);
    if (sidx < 0) {
      std::cerr << "[WARN] Cannot identify pi sign (pi+/pi-) from path, skipping: " << yrPath << "\n";
      nSkipSign++;
      continue;
    }

    // tau + err
    double tau = 0.0, tauErr = 0.0;
    bool hasTauErr = false;
    if (!ParseTauHeader(yrPath, tau, tauErr, hasTauErr)) {
      std::cerr << "[WARN] Missing/invalid tau_ns in header, skipping: " << yrPath << "\n";
      nSkipTau++;
      continue;
    }
    if (!hasTauErr) {
      std::cerr << "[WARN] tau_ns_err missing/invalid; skipping (needed for rel-err cut): " << yrPath << "\n";
      nSkipNoTauErr++;
      continue;
    }
    const double relErr = tauErr / std::abs(tau);
    if (!std::isfinite(relErr) || relErr >= kMaxRelTauErr) {
      std::cerr << "[WARN] tau rel-err too large (tauErr/|tau|=" << relErr
                << " >= " << kMaxRelTauErr << "), skipping: " << yrPath << "\n";
      nSkipRelErr++;
      continue;
    }

    // fit runs + mean rate
    std::vector<int> fitRuns;
    double meanRate = std::nan(""), rmsRate = 0.0;
    int nFit = 0;
    if (!ReadFitRunsAndRate(yrPath, fitRuns, meanRate, rmsRate, nFit)) {
      std::cerr << "[WARN] <2 fit points or could not compute mean shms34_rate_kHz; skipping: " << yrPath << "\n";
      nSkipFitPts++;
      continue;
    }
    if (nFit < kMinNFit2) {
      std::cerr << "[WARN] Too few fit points (N_valid_for_fit2=" << nFit
                << " < " << kMinNFit2 << "), skipping: " << yrPath << "\n";
      nSkipMinFitPts++;
      continue;
    }

    // current file
    TString curPath = yrPath;
    curPath.ReplaceAll(kYieldRatioName, kCurrentName);

    if (gSystem->AccessPathName(curPath)) {
      std::cerr << "[WARN] Current CSV not found, skipping: " << curPath << "\n"
                << "       (paired with) " << yrPath << "\n";
      nSkipNoCurrentFile++;
      continue;
    }

    std::unordered_map<int,double> run2I;
    if (!ReadRunToCurrentMap(curPath, run2I)) {
      std::cerr << "[WARN] Could not parse run->BCM2_I map, skipping: " << curPath << "\n"
                << "       (paired with) " << yrPath << "\n";
      nSkipCurrentParse++;
      continue;
    }

    double meanI = std::nan(""), rmsI = 0.0;
    int nRunMatch = 0;
    if (!ComputeMeanCurrentOverRuns(fitRuns, run2I, meanI, rmsI, nRunMatch)) {
      std::cerr << "[WARN] No matching run(s) for current among fit runs; skipping: " << curPath << "\n"
                << "       (paired with) " << yrPath << "\n";
      nSkipNoRunMatch++;
      continue;
    }

    // store
    rateBuckets[pad][sidx].x.push_back(meanRate);
    rateBuckets[pad][sidx].ex.push_back(0.0);
    rateBuckets[pad][sidx].y.push_back(tau);
    rateBuckets[pad][sidx].ey.push_back(tauErr);

    currBuckets[pad][sidx].x.push_back(meanI);
    currBuckets[pad][sidx].ex.push_back(0.0);
    currBuckets[pad][sidx].y.push_back(tau);
    currBuckets[pad][sidx].ey.push_back(tauErr);

    TString sp = yrPath;
    sp.ReplaceAll(kResultsDir + "/", "");
    logs.push_back({sp, pass, target, sign, tau, tauErr, relErr, nFit, meanRate, rmsRate, nRunMatch, meanI, rmsI});

    nProcessed++;
  }

  // Summary
  std::cout << "\n[TauVsRateAndCurrent] Filters applied:\n"
            << "  Require N_valid_for_fit2 >= " << kMinNFit2 << "\n"
            << "  Require tau_err/|tau| < " << kMaxRelTauErr << "\n\n";

  std::cout << "[TauVsRateAndCurrent] Summary\n"
            << "  Processed (plotted):           " << nProcessed << "\n"
            << "  Skipped (bad pass/target):     " << nSkipBadCat << "\n"
            << "  Skipped (pi sign unknown):     " << nSkipSign << "\n"
            << "  Skipped (missing/invalid tau): " << nSkipTau << "\n"
            << "  Skipped (missing/invalid terr):" << nSkipNoTauErr << "\n"
            << "  Skipped (tau rel-err cut):     " << nSkipRelErr << "\n"
            << "  Skipped (<2 fit points):       " << nSkipFitPts << "\n"
            << "  Skipped (Nfit < min):          " << nSkipMinFitPts << "\n"
            << "  Skipped (missing current CSV): " << nSkipNoCurrentFile << "\n"
            << "  Skipped (current parse fail):  " << nSkipCurrentParse << "\n"
            << "  Skipped (no run match):        " << nSkipNoRunMatch << "\n\n";

  std::cout << "[TauVsRateAndCurrent] Per-setting fit spread log (rate RMS & current RMS)\n";
  std::cout << "  pass target sign  tau(ns)  tauErr  relErr  nFit  meanRate(kHz)  rmsRate  nRunMatch  meanI(uA)  rmsI   | path\n";
  for (const auto& L : logs) {
    std::cout
      << "  " << L.pass
      << "   " << L.target
      << "    " << L.sign
      << "   " << Form("%7.2f", L.tau)
      << "  " << Form("%6.2f", L.tauerr)
      << "  " << Form("%6.3f", L.relerr)
      << "   " << Form("%3d", L.nFit)
      << "     " << Form("%10.3f", L.meanRate)
      << "   " << Form("%7.3f", L.rmsRate)
      << "     " << Form("%3d", L.nRunMatch)
      << "     " << Form("%8.3f", L.meanI)
      << "   " << Form("%7.3f", L.rmsI)
      << "   | " << L.shortPath
      << "\n";
  }
  std::cout << "\n";

  // Plot helpers
  auto MakeCanvas = [&](const TString& name, const TString& title) -> TCanvas* {
    TCanvas* c = new TCanvas(name, title, 1200, 900);
    c->Divide(2, 2);
    return c;
  };

  auto DrawPad = [&](const TString& padTitle,
                     Bucket& bPlus, Bucket& bMinus,
                     const TString& xTitle, const TString& yTitle) {

    gPad->SetTicks(1, 1);

    // Compute x-range from both series
    double xmin = 0, xmax = 1;
    bool haveX = ComputeXRange(bPlus.x, bMinus.x, xmin, xmax);

    // Build graphs
    TGraphErrors* gPlus  = new TGraphErrors((int)bPlus.x.size());
    TGraphErrors* gMinus = new TGraphErrors((int)bMinus.x.size());

    for (int i = 0; i < (int)bPlus.x.size(); i++) {
      gPlus->SetPoint(i, bPlus.x[i], bPlus.y[i]);
      gPlus->SetPointError(i, bPlus.ex[i], bPlus.ey[i]);
    }
    for (int i = 0; i < (int)bMinus.x.size(); i++) {
      gMinus->SetPoint(i, bMinus.x[i], bMinus.y[i]);
      gMinus->SetPointError(i, bMinus.ex[i], bMinus.ey[i]);
    }

    StyleGraph(gPlus, "pi+");
    StyleGraph(gMinus, "pi-");

    // y-range from data across both series
    double yMin = 1e300, yMax = -1e300;
    auto updateY = [&](const Bucket& b) {
      for (size_t i = 0; i < b.y.size(); i++) {
        yMin = std::min(yMin, b.y[i] - b.ey[i]);
        yMax = std::max(yMax, b.y[i] + b.ey[i]);
      }
    };
    updateY(bPlus);
    updateY(bMinus);

    if (!(yMin < yMax) || !std::isfinite(yMin) || !std::isfinite(yMax)) {
      yMin = 0.0; yMax = 300.0;
    } else {
      const double pad = 0.10 * (yMax - yMin);
      yMin -= pad;
      yMax += pad;
      if (yMin < 0) yMin = 0;
    }

    // Draw with pi+ as frame, but force x-limits to include both series
    gPlus->SetTitle(padTitle);
    gPlus->GetXaxis()->SetTitle(xTitle);
    gPlus->GetYaxis()->SetTitle(yTitle);
    gPlus->GetYaxis()->SetRangeUser(yMin, yMax);
    if (haveX) gPlus->GetXaxis()->SetLimits(xmin, xmax);

    gPlus->Draw("AP");
    gMinus->Draw("P SAME");

    // Legend
    TLegend* leg = new TLegend(0.12, 0.78, 0.36, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gPlus, "pi+", "p");
    leg->AddEntry(gMinus, "pi-", "p");
    leg->Draw("same");

    // N text
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.12, 0.72, Form("N(pi+)=%zu, N(pi-)=%zu", bPlus.x.size(), bMinus.x.size()));
  };

  const char* padTitles[4] = {"Pass4, LH2", "Pass4, LD2", "Pass5, LH2", "Pass5, LD2"};

  // Plot 1: tau vs shms34 rate
  {
    TCanvas* cR = MakeCanvas("c_tau_vs_shms34", "Tau vs SHMS 3/4 trigger rate");
    for (int i = 0; i < 4; i++) {
      cR->cd(i + 1);
      DrawPad(padTitles[i],
              rateBuckets[i][0], rateBuckets[i][1],
              "Mean SHMS 3/4 trigger rate (kHz) [use_in_fit==1]",
              "#tau (ns)");
    }
    cR->SaveAs(outPngRate);
    std::cout << "[TauVsRateAndCurrent] Saved: " << outPngRate << "\n";
  }

  // Plot 2: tau vs current
  {
    TCanvas* cI = MakeCanvas("c_tau_vs_current", "Tau vs beam current");
    for (int i = 0; i < 4; i++) {
      cI->cd(i + 1);
      DrawPad(padTitles[i],
              currBuckets[i][0], currBuckets[i][1],
              "Mean beam current BCM2_I (#muA) [fit runs]",
              "#tau (ns)");
    }
    cI->SaveAs(outPngCurr);
    std::cout << "[TauVsRateAndCurrent] Saved: " << outPngCurr << "\n";
  }
}
