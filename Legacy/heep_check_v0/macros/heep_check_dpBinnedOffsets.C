// heep_check_dpBinnedOffsets.C
//
// Purpose:
//   For a single (or list of) 4-pass HEEP run(s), make 6 separate PNGs per run:
//     - 1 over full SHMS dp range: [-10, 20]
//     - 5 over user-provided focus dp range: [dpFocusLo, dpFocusHi] split into 5 equal bins
//   Each PNG is a 2x2 canvas with overlays (Data vs Sim) for:
//     { W, Em, Pmz, Pmy }
//   Sim is filled with event-by-event Weight from h10 tree.
//   Sim is area-normalized to Data in a restricted "peak window" (adaptive).
//   Both Data and Sim are fit with a Gaussian when possible; report Δμ ± σ(Δμ)
//   on the plot in red for valid fits.
//   Append results into a CSV with newline-safe row writing.
//
// Run from: RP_Scripts/heep_check_v0/
// Example:
//   root -l -b -q 'heep_check_dpBinnedOffsets.C+(23839,-2,3)'
// or loop all 4-pass runs hardcoded below:
//   root -l -b -q 'heep_check_dpBinnedOffsets.C+(0,-2,3,5,true)'

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TPad.h"

// -----------------------------------------------------------------------------
// Directories (edit if your paths change)
// -----------------------------------------------------------------------------
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check_v0/Pass0_ROOTfiles/";
static const TString kSimDir  = "/home/cdaq/users/rparvez/simc_gfortran/worksim/heep_check_v0/";

// File patterns
static TString GetDataFileName(int run) {
  return kDataDir + Form("coin_replay_production_%d_-1.root", run);
}
static TString GetSimFileName(int run) {
  return kSimDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
}

// Trees
static const TString kDataTreeName = "T";
static const TString kSimTreeName  = "h10";

// Baseline cuts (same spirit as your original macro, but with dp-range moved to dp-bin)
static const char* kDataCutBase = "(H.gtr.dp>-8) && (H.gtr.dp<8)";
static const char* kSimCutBase  = "(hsdelta>-8) && (hsdelta<8)";

// Sim weight branch
static const char* kSimWeightBranch = "Weight";

// -----------------------------------------------------------------------------
// Variable definitions
// -----------------------------------------------------------------------------
struct VarInfo {
  TString name;
  TString data_expr;
  TString sim_expr;
  TString x_title;
  int nbins;
  double xmin;
  double xmax;
};

static std::vector<VarInfo> GetVariables() {
  std::vector<VarInfo> v;

  // Keep same binning as your previous macro for consistency
  v.push_back({"W",
               "P.kin.primary.W",
               "W_recon",
               "W (GeV)",
               100, 0.6, 1.3});

  v.push_back({"Em",
               "H.kin.secondary.emiss",
               "Em_recon",
               "Em (GeV)",
               100, -0.05, 0.15});

  v.push_back({"Pmz",
               "H.kin.secondary.pmiss_z",
               "Pmz_recon",
               "Pmz (GeV/c)",
               100, -0.3, 0.3});

  v.push_back({"Pmy",
               "H.kin.secondary.pmiss_y",
               "Pmy_recon",
               "Pmy (GeV/c)",
               100, -0.3, 0.3});

  return v;
}

// -----------------------------------------------------------------------------
// DP bin definitions
// -----------------------------------------------------------------------------
struct DpBin {
  int idx;          // 0..5
  TString label;    // "full" or "b1".."b5"
  double lo;
  double hi;
};

static std::vector<DpBin> BuildDpBins(double dpFocusLo, double dpFocusHi, int nFocusBins) {
  std::vector<DpBin> bins;

  // Full-range bin always first (hard-coded)
  bins.push_back({0, "full", -10.0, 20.0});

  // Focus bins
  if (nFocusBins < 1) nFocusBins = 1;
  double step = (dpFocusHi - dpFocusLo) / nFocusBins;
  for (int i = 0; i < nFocusBins; i++) {
    double lo = dpFocusLo + i * step;
    double hi = dpFocusLo + (i + 1) * step;
    bins.push_back({i + 1, Form("b%d", i + 1), lo, hi});
  }
  return bins;
}

// -----------------------------------------------------------------------------
// Styling
// -----------------------------------------------------------------------------
static void StyleHists(TH1D* hData, TH1D* hSim) {
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);

  hSim->SetMarkerStyle(24);
  hSim->SetMarkerSize(0.7);
  hSim->SetLineColor(kRed+1);

  hData->GetYaxis()->SetTitle("Counts");
  hData->GetYaxis()->SetTitleOffset(1.15);
  hData->GetXaxis()->SetTitleOffset(1.05);

  hData->GetXaxis()->SetTitleSize(0.058);
  hData->GetXaxis()->SetLabelSize(0.048);
  hData->GetYaxis()->SetTitleSize(0.058);
  hData->GetYaxis()->SetLabelSize(0.048);

}

// -----------------------------------------------------------------------------
// Peak window estimation utilities
// -----------------------------------------------------------------------------
struct PeakWindow {
  double xLo = 0;
  double xHi = 0;
  double muSeed = 0;
  double sigSeed = 0;
  bool ok = false;
};

static PeakWindow EstimatePeakWindow(const TH1D* h,
                                     int seedHalfWindowBins = 10, // "±10 bins" around max
                                     double kSigma = 2.0,
                                     int minHalfWidthBins = 4,
                                     int maxHalfWidthBins = 25) {
  PeakWindow pw;
  if (!h) return pw;

  const int nb = h->GetNbinsX();
  if (nb < 5) return pw;

  int bMax = h->GetMaximumBin();
  if (bMax < 1 || bMax > nb) return pw;

  int bLo = std::max(1, bMax - seedHalfWindowBins);
  int bHi = std::min(nb, bMax + seedHalfWindowBins);

  double sumW = 0.0, sumWX = 0.0, sumWXX = 0.0;
  for (int b = bLo; b <= bHi; b++) {
    double w = h->GetBinContent(b);
    double x = h->GetXaxis()->GetBinCenter(b);
    if (w <= 0) continue;
    sumW   += w;
    sumWX  += w * x;
    sumWXX += w * x * x;
  }
  if (sumW <= 0) return pw;

  double mu = sumWX / sumW;
  double var = sumWXX / sumW - mu * mu;
  double sig = (var > 0 ? std::sqrt(var) : 0.0);

  // Fallback for pathological sigma
  double binW = h->GetXaxis()->GetBinWidth(1);
  if (!(sig > 0)) sig = seedHalfWindowBins * binW * 0.35;

  // Build normalization window ±kSigma*sig around seed mean
  double halfW = kSigma * sig;

  // Clamp half-width in bins
  double minHalfW = minHalfWidthBins * binW;
  double maxHalfW = maxHalfWidthBins * binW;
  halfW = std::max(halfW, minHalfW);
  halfW = std::min(halfW, maxHalfW);

  double xLo = mu - halfW;
  double xHi = mu + halfW;

  // Clamp to axis range
  double axLo = h->GetXaxis()->GetXmin();
  double axHi = h->GetXaxis()->GetXmax();
  xLo = std::max(xLo, axLo);
  xHi = std::min(xHi, axHi);

  pw.xLo = xLo;
  pw.xHi = xHi;
  pw.muSeed = mu;
  pw.sigSeed = sig;
  pw.ok = (xHi > xLo);
  return pw;
}

static double IntegralInXRange(const TH1D* h, double xLo, double xHi) {
  if (!h) return 0.0;
  if (xHi <= xLo) return 0.0;

  int b1 = h->GetXaxis()->FindBin(xLo);
  int b2 = h->GetXaxis()->FindBin(xHi);

  b1 = std::max(1, std::min(b1, h->GetNbinsX()));
  b2 = std::max(1, std::min(b2, h->GetNbinsX()));
  if (b2 < b1) std::swap(b1, b2);

  return h->Integral(b1, b2);
}

// -----------------------------------------------------------------------------
// Fit utilities
// -----------------------------------------------------------------------------
struct FitOut {
  double mu = 0.0;
  double muErr = 0.0;
  double sig = 0.0;
  double sigErr = 0.0;
  int status = -999;
  double chi2 = 0.0;
  int ndf = 0;
  bool ok = false;
  TF1* f = nullptr;
};

static FitOut FitGaus(TH1D* h, const TString& fname,
                      double xLo, double xHi,
                      Color_t lineColor) {
  FitOut out;
  if (!h) return out;
  if (!(xHi > xLo)) return out;
  if (h->GetEntries() < 30) return out;

  TF1* f = new TF1(fname, "gaus", xLo, xHi);
  f->SetLineColor(lineColor);
  f->SetLineWidth(2);

  // Q: quiet, R: use range, S: return fit result, N: don't draw automatically
  TFitResultPtr r = h->Fit(f, "QRSN");

  out.f = f;
  out.status = r.Get() ? r->Status() : -999;

  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf  = r->Ndf();
  }

  out.mu     = f->GetParameter(1);
  out.muErr  = f->GetParError(1);
  out.sig    = f->GetParameter(2);
  out.sigErr = f->GetParError(2);

  out.ok = (out.status == 0);
  return out;
}

// -----------------------------------------------------------------------------
// CSV writer (append-safe, newline-safe)
// -----------------------------------------------------------------------------
static void EnsureCsvHeader(const TString& csvPath) {
  // If file doesn't exist, create with header
  if (!gSystem->AccessPathName(csvPath)) return; // exists

  // Ensure parent dir exists
  TString parent = csvPath;
  parent = gSystem->DirName(parent);
  if (!parent.IsNull() && gSystem->AccessPathName(parent)) {
    gSystem->mkdir(parent, kTRUE);
  }

  std::ofstream csv(csvPath.Data(), std::ios::out);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot create CSV: " << csvPath << "\n";
    return;
  }

  csv
    << "run,dp_idx,dp_label,dp_lo,dp_hi,"
    << "var,"
    << "norm_win_lo,norm_win_hi,"
    << "intD_win,intS_win,sim_scale,"
    << "entriesD,entriesS,"
    << "muD,muD_err,sigD,sigD_err,statusD,chi2D,ndfD,"
    << "muS,muS_err,sigS,sigS_err,statusS,chi2S,ndfS,"
    << "dmu_GeV,dmu_err_GeV,dmu_MeV,dmu_err_MeV,fit_valid"
    << "\n";

  csv.close();
}

static void AppendCsvRow(const TString& csvPath,
                         int run, const DpBin& dp, const VarInfo& V,
                         const PeakWindow& pw,
                         double intDwin, double intSwin, double simScale,
                         double entriesD, double entriesS,
                         const FitOut& fD, const FitOut& fS,
                         double dmu, double dmuErr, double dmuMeV, double dmuErrMeV, bool fitValid) {
  std::ofstream csv(csvPath.Data(), std::ios::out | std::ios::app);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot append CSV: " << csvPath << "\n";
    return;
  }

  csv.setf(std::ios::fixed);
  csv << std::setprecision(6);

  csv
    << run << ","
    << dp.idx << ","
    << dp.label.Data() << ","
    << dp.lo << ","
    << dp.hi << ","
    << V.name.Data() << ","
    << pw.xLo << ","
    << pw.xHi << ","
    << intDwin << ","
    << intSwin << ","
    << simScale << ","
    << entriesD << ","
    << entriesS << ","
    << fD.mu << ","
    << fD.muErr << ","
    << fD.sig << ","
    << fD.sigErr << ","
    << fD.status << ","
    << fD.chi2 << ","
    << fD.ndf << ","
    << fS.mu << ","
    << fS.muErr << ","
    << fS.sig << ","
    << fS.sigErr << ","
    << fS.status << ","
    << fS.chi2 << ","
    << fS.ndf << ","
    << dmu << ","
    << dmuErr << ","
    << dmuMeV << ","
    << dmuErrMeV << ","
    << (fitValid ? 1 : 0)
    << "\n";

  csv.close();
}

// -----------------------------------------------------------------------------
// Main worker: one run
// -----------------------------------------------------------------------------
static void ProcessOneRun(int run,
                          double dpFocusLo, double dpFocusHi, int nFocusBins,
                          const TString& outDir,
                          const TString& csvPath) {
  TString dataFile = GetDataFileName(run);
  TString simFile  = GetSimFileName(run);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) { std::cerr << "[ERROR] Cannot open data file: " << dataFile << "\n"; return; }

  TFile* fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) { std::cerr << "[ERROR] Cannot open sim file:  " << simFile  << "\n"; return; }

  TTree* tData = (TTree*) fData->Get(kDataTreeName);
  TTree* tSim  = (TTree*) fSim->Get(kSimTreeName);
  if (!tData || !tSim) { std::cerr << "[ERROR] Missing trees in files for run " << run << "\n"; return; }

  std::vector<VarInfo> vars = GetVariables();
  std::vector<DpBin> dpBins = BuildDpBins(dpFocusLo, dpFocusHi, nFocusBins);

  for (const auto& dp : dpBins) {
    // One PNG per dp bin
    TCanvas* c = new TCanvas(Form("c_heep_%d_%s", run, dp.label.Data()),
                             Form("HEEP offsets run %d dp %s", run, dp.label.Data()),
                             1600, 1000);
    c->Divide(2, 2);

    // Build dp-cuts
    TString dataDpCut = Form("(P.gtr.dp>%g) && (P.gtr.dp<%g)", dp.lo, dp.hi);
    TString simDpCut  = Form("(ssdelta>%g) && (ssdelta<%g)", dp.lo, dp.hi);

    TString dataCut = Form("(%s) && (%s)", kDataCutBase, dataDpCut.Data());
    TString simCut  = Form("(%s) && (%s)", kSimCutBase,  simDpCut.Data());

    for (int i = 0; i < (int)vars.size(); i++) {
      const VarInfo& V = vars[i];
      c->cd(i + 1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.08);

      TString hDname = Form("hD_%s_r%d_%s", V.name.Data(), run, dp.label.Data());
      TString hSname = Form("hS_%s_r%d_%s", V.name.Data(), run, dp.label.Data());

      TH1D* hD = new TH1D(hDname, "", V.nbins, V.xmin, V.xmax);
      TH1D* hS = new TH1D(hSname, "", V.nbins, V.xmin, V.xmax);
      hD->Sumw2();
      hS->Sumw2();

      // Fill
      tData->Draw(Form("%s>>%s", V.data_expr.Data(), hDname.Data()), dataCut.Data(), "goff");

      // Sim selection = Weight * (cuts)
      TString simSel = Form("%s*(%s)", kSimWeightBranch, simCut.Data());
      tSim->Draw(Form("%s>>%s", V.sim_expr.Data(), hSname.Data()), simSel.Data(), "goff");

      // Determine restricted normalization window from DATA (adaptive)
      PeakWindow pw = EstimatePeakWindow(hD, /*seedHalfWindowBins=*/10, /*kSigma=*/2.0,
                                         /*minHalfWidthBins=*/4, /*maxHalfWidthBins=*/25);

      double intDwin = 0.0, intSwin = 0.0, simScale = 1.0;
      if (pw.ok) {
        intDwin = IntegralInXRange(hD, pw.xLo, pw.xHi);
        intSwin = IntegralInXRange(hS, pw.xLo, pw.xHi);
        if (intSwin > 0 && intDwin > 0) {
          simScale = intDwin / intSwin;
          hS->Scale(simScale);
        }
      }

      // Fit range for DATA from DATA peak window
      double fitLoD = pw.ok ? pw.xLo : V.xmin;
      double fitHiD = pw.ok ? pw.xHi : V.xmax;

      // After scaling, estimate a separate peak window for SIM
      PeakWindow pwS = EstimatePeakWindow(hS, /*seedHalfWindowBins=*/10, /*kSigma=*/2.0,
                                          /*minHalfWidthBins=*/4, /*maxHalfWidthBins=*/25);

      double fitLoS = pwS.ok ? pwS.xLo : fitLoD;
      double fitHiS = pwS.ok ? pwS.xHi : fitHiD;

      FitOut fD = FitGaus(hD, Form("fD_%s_r%d_%s", V.name.Data(), run, dp.label.Data()),
                          fitLoD, fitHiD, kBlack);

      FitOut fS = FitGaus(hS, Form("fS_%s_r%d_%s", V.name.Data(), run, dp.label.Data()),
                          fitLoS, fitHiS, kRed+1);

      // Δμ and error: only valid if BOTH fits are good
      double dmu       = -999.0;
      double dmuErr    = -999.0;
      double dmuMeV    = -999.0;
      double dmuErrMeV = -999.0;

      if (fD.ok && fS.ok) {
        dmu       = fD.mu - fS.mu;
        dmuErr    = std::sqrt(fD.muErr * fD.muErr + fS.muErr * fS.muErr);
        dmuMeV    = dmu * 1000.0;
        dmuErrMeV = dmuErr * 1000.0;
      }

      // Draw
      StyleHists(hD, hS);
      hD->GetXaxis()->SetTitle(V.x_title);

      double maxY = std::max(hD->GetMaximum(), hS->GetMaximum());
      hD->SetMaximum(1.25 * maxY);

      hD->Draw("E1");
      hS->Draw("E1 SAME");

      if (fD.f) fD.f->Draw("SAME");
      if (fS.f) fS.f->Draw("SAME");

      // Legend
      TLegend* L = new TLegend(0.58, 0.72, 0.90, 0.88);
      L->SetBorderSize(0);
      L->SetFillStyle(0);
      L->SetTextSize(0.042);
      L->AddEntry(hD, "Data", "lp");
      L->AddEntry(hS, "Sim (weighted, scaled)", "lp");
      L->Draw();

      // Text annotations
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.058);

      // dp bin label
      lat.SetTextColor(kBlack);
      lat.DrawLatex(0.50, 0.92, Form("Run %d, dp[%g,%g] (%s)", run, dp.lo, dp.hi, dp.label.Data()));

      // Fit means
      lat.SetTextSize(0.050);
      lat.DrawLatex(0.15, 0.84, Form("Data: #mu=%.5f #pm %.5f", fD.mu, fD.muErr));
      lat.DrawLatex(0.15, 0.79, Form("Sim:  #mu=%.5f #pm %.5f", fS.mu, fS.muErr));

      // Delta mu in red
      lat.SetTextColor(kRed+1);
      lat.SetTextSize(0.062);
      if (fD.ok && fS.ok) {
        lat.DrawLatex(0.15, 0.72, Form("#Delta#mu = %.2f #pm %.2f MeV", dmuMeV, dmuErrMeV));
      } else {
        lat.DrawLatex(0.15, 0.72, Form("#Delta#mu invalid  (statusD=%d, statusS=%d)", fD.status, fS.status));
      }
      // Optional: show normalization window (small)
      lat.SetTextColor(kGray+2);
      lat.SetTextSize(0.042);
      if (pw.ok) lat.DrawLatex(0.15, 0.66, Form("NormWin: [%.5f, %.5f], SimScale=%.4g", pw.xLo, pw.xHi, simScale));

      // CSV row
      AppendCsvRow(csvPath, run, dp, V, pw,
                   intDwin, intSwin, simScale,
                   hD->GetEntries(), hS->GetEntries(),
                   fD, fS, dmu, dmuErr, dmuMeV, dmuErrMeV, (fD.ok && fS.ok));
    }

    // Save PNG
    TString outName = outDir + Form("/heep_offsets_run%d_dp%s.png", run, dp.label.Data());
    c->SaveAs(outName);
    std::cout << "[INFO] Saved: " << outName << "\n";

    delete c;
  }

  fData->Close();
  fSim->Close();
  delete fData;
  delete fSim;
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------
void heep_check_dpBinnedOffsets(int run = 0,
                                double dpFocusLo = -2.0,
                                double dpFocusHi =  3.0,
                                int nFocusBins = 5,
                                bool loopRuns = false,
                                const char* outDirC = "results/PNGs/dp_binned",
                                const char* csvPathC = "results/tables/heep_check_dpBinnedOffsets.csv") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString outDir(outDirC);
  TString csvPath(csvPathC);

  // Ensure output dirs exist
  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);
  EnsureCsvHeader(csvPath);

  // Runs: either single or a hardcoded list (your 12 four-pass runs)
  std::vector<int> runs = {
    23839, 23840, 23841, 23842, 23844, 23845,
    23846, 23847, 23848, 23849, 23850, 23851
  };

  if (!loopRuns && run != 0) {
    ProcessOneRun(run, dpFocusLo, dpFocusHi, nFocusBins, outDir, csvPath);
    return;
  }

  // Loop mode
  for (int r : runs) {
    ProcessOneRun(r, dpFocusLo, dpFocusHi, nFocusBins, outDir, csvPath);
  }
}
