// WandEmPeakVsDelta.C
//
// Purpose:
//   For the 12 hardcoded 4-pass HEEP data runs, measure W_peak and Em_peak
//   versus SHMS delta bin using DATA only.
//   For each run and each dp bin (full + 5 narrow bins), make a histogram,
//   fit a Gaussian, and write the fitted peak information to a CSV.
//
// Output CSV:
//   heep_check_v1/results/tables/WandEmPeakVsDelta.csv
//
// Notes:
//   - Tree name: T
//   - W branch:  P.kin.primary.W
//   - Em branch: H.kin.secondary.emiss
//   - SHMS dp:   P.gtr.dp
//   - Base data cut: (H.gtr.dp>-8) && (H.gtr.dp<8)
//   - Units in CSV:
//       peak_GeV, peak_err_GeV, sigma_GeV, sigma_err_GeV
//       peak_MeV, peak_err_MeV, sigma_MeV, sigma_err_MeV
//   - For failed fits, fit-derived quantities are written as NaN
//   - For full bin, dp_center is written as -999.0 sentinel
//
// Example:
//   root -l -b -q 'macros/WandEmPeakVsDelta.C+(0,-2,3,5,true,"results/tables/WandEmPeakVsDelta.csv",true)'

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

// -----------------------------------------------------------------------------
// Paths
// -----------------------------------------------------------------------------
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check_v1/Pass0_ROOTfiles/";
static const TString kDataTreeName = "T";

// Base data cut
static const char* kDataCutBase = "(H.gtr.dp>-8) && (H.gtr.dp<8)";

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static TString GetDataFileName(int run) {
  return kDataDir + Form("coin_replay_production_%d_-1.root", run);
}

static TString GetSettingLabel(int run) {
  if (run >= 23839 && run <= 23848 && run != 23843) return "A";
  if (run >= 23849 && run <= 23851) return "B";
  return "UNKNOWN";
}

// -----------------------------------------------------------------------------
// Variable definitions
// -----------------------------------------------------------------------------
struct VarInfo {
  TString name;      // "W" or "Em"
  TString expr;      // TTree expression
  TString x_title;   // axis title
  int nbins;
  double xmin;
  double xmax;
};

static std::vector<VarInfo> GetVariables() {
  std::vector<VarInfo> v;

  // Reuse same histogram choices as MeasuredOffsetsByRun.C
  v.push_back({"W",
               "P.kin.primary.W",
               "W (GeV)",
               100, 0.6, 1.3});

  v.push_back({"Em",
               "H.kin.secondary.emiss",
               "Em (GeV)",
               100, -0.05, 0.15});

  return v;
}

// -----------------------------------------------------------------------------
// DP bins
// -----------------------------------------------------------------------------
struct DpBin {
  int idx;          // 0..5
  TString label;    // "full", "b1"... "b5"
  double lo;
  double hi;
  double center;    // -999.0 for full
};

static std::vector<DpBin> BuildDpBins(double dpFocusLo, double dpFocusHi, int nFocusBins) {
  std::vector<DpBin> bins;

  bins.push_back({0, "full", -10.0, 20.0, -999.0});

  if (nFocusBins < 1) nFocusBins = 1;
  double step = (dpFocusHi - dpFocusLo) / nFocusBins;
  for (int i = 0; i < nFocusBins; i++) {
    double lo = dpFocusLo + i * step;
    double hi = dpFocusLo + (i + 1) * step;
    double ctr = 0.5 * (lo + hi);
    bins.push_back({i + 1, Form("b%d", i + 1), lo, hi, ctr});
  }
  return bins;
}

// -----------------------------------------------------------------------------
// Peak-window estimation
// -----------------------------------------------------------------------------
struct PeakWindow {
  double xLo = 0.0;
  double xHi = 0.0;
  double muSeed = 0.0;
  double sigSeed = 0.0;
  bool ok = false;
};

static PeakWindow EstimatePeakWindow(const TH1D* h,
                                     int seedHalfWindowBins = 10,
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

  double binW = h->GetXaxis()->GetBinWidth(1);
  if (!(sig > 0.0)) sig = seedHalfWindowBins * binW * 0.35;

  double halfW = kSigma * sig;
  double minHalfW = minHalfWidthBins * binW;
  double maxHalfW = maxHalfWidthBins * binW;

  halfW = std::max(halfW, minHalfW);
  halfW = std::min(halfW, maxHalfW);

  double xLo = mu - halfW;
  double xHi = mu + halfW;

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

// -----------------------------------------------------------------------------
// Fit utilities
// -----------------------------------------------------------------------------
struct FitOut {
  double mu = std::numeric_limits<double>::quiet_NaN();
  double muErr = std::numeric_limits<double>::quiet_NaN();
  double sig = std::numeric_limits<double>::quiet_NaN();
  double sigErr = std::numeric_limits<double>::quiet_NaN();
  int status = -999;
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = -999;
  bool ok = false;
};

static FitOut FitGaus(TH1D* h, const TString& fname,
                      double xLo, double xHi) {
  FitOut out;
  if (!h) return out;
  if (!(xHi > xLo)) return out;
  if (h->GetEntries() < 30) return out;

  TF1 f(fname, "gaus", xLo, xHi);

  TFitResultPtr r = h->Fit(&f, "QRSN");
  out.status = r.Get() ? r->Status() : -999;

  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf  = r->Ndf();
  }

  out.mu     = f.GetParameter(1);
  out.muErr  = f.GetParError(1);
  out.sig    = f.GetParameter(2);
  out.sigErr = f.GetParError(2);

  bool finiteVals = std::isfinite(out.mu) && std::isfinite(out.muErr) &&
                    std::isfinite(out.sig) && std::isfinite(out.sigErr);
  bool saneErrors = (out.muErr >= 0.0) && (out.sigErr >= 0.0);
  bool saneWidth  = (out.sig > 0.0);

  out.ok = (out.status == 0) && finiteVals && saneErrors && saneWidth;
  return out;
}

// -----------------------------------------------------------------------------
// CSV writing
// -----------------------------------------------------------------------------
static void EnsureCsvHeader(const TString& csvPath) {
  if (!gSystem->AccessPathName(csvPath)) return;

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
    << "run,setting,ebeam_MeV,"
    << "dp_idx,dp_label,dp_lo,dp_hi,dp_center,"
    << "var,"
    << "entries,"
    << "fit_lo_GeV,fit_hi_GeV,mu_seed_GeV,sig_seed_GeV,"
    << "peak_GeV,peak_err_GeV,sigma_GeV,sigma_err_GeV,"
    << "peak_MeV,peak_err_MeV,sigma_MeV,sigma_err_MeV,"
    << "fit_status,chi2,ndf,chi2_ndf,fit_valid"
    << "\n";

  csv.close();
}

static void AppendCsvRow(const TString& csvPath,
                         int run,
                         const TString& setting,
                         double ebeamMeV,
                         const DpBin& dp,
                         const VarInfo& V,
                         double entries,
                         const PeakWindow& pw,
                         const FitOut& f) {
  std::ofstream csv(csvPath.Data(), std::ios::out | std::ios::app);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot append CSV: " << csvPath << "\n";
    return;
  }

  const double NaN = std::numeric_limits<double>::quiet_NaN();

  double peakGeV     = f.ok ? f.mu     : NaN;
  double peakErrGeV  = f.ok ? f.muErr  : NaN;
  double sigmaGeV    = f.ok ? f.sig    : NaN;
  double sigmaErrGeV = f.ok ? f.sigErr : NaN;

  double peakMeV     = f.ok ? 1000.0 * f.mu     : NaN;
  double peakErrMeV  = f.ok ? 1000.0 * f.muErr  : NaN;
  double sigmaMeV    = f.ok ? 1000.0 * f.sig    : NaN;
  double sigmaErrMeV = f.ok ? 1000.0 * f.sigErr : NaN;

  double chi2ndf = (f.ok && f.ndf > 0 && std::isfinite(f.chi2)) ? (f.chi2 / f.ndf) : NaN;

  csv.setf(std::ios::fixed);
  csv << std::setprecision(6);

  csv
    << run << ","
    << setting.Data() << ","
    << ebeamMeV << ","
    << dp.idx << ","
    << dp.label.Data() << ","
    << dp.lo << ","
    << dp.hi << ","
    << dp.center << ","
    << V.name.Data() << ","
    << entries << ","
    << (pw.ok ? pw.xLo     : NaN) << ","
    << (pw.ok ? pw.xHi     : NaN) << ","
    << (pw.ok ? pw.muSeed  : NaN) << ","
    << (pw.ok ? pw.sigSeed : NaN) << ","
    << peakGeV << ","
    << peakErrGeV << ","
    << sigmaGeV << ","
    << sigmaErrGeV << ","
    << peakMeV << ","
    << peakErrMeV << ","
    << sigmaMeV << ","
    << sigmaErrMeV << ","
    << f.status << ","
    << f.chi2 << ","
    << f.ndf << ","
    << chi2ndf << ","
    << (f.ok ? 1 : 0)
    << "\n";

  csv.close();
}

// -----------------------------------------------------------------------------
// One-run worker
// -----------------------------------------------------------------------------
static void ProcessOneRun(int run,
                          double dpFocusLo, double dpFocusHi, int nFocusBins,
                          const TString& csvPath,
                          bool verbose = true) {
  const TString setting = GetSettingLabel(run);
  const double ebeamMeV = 8583.1;

  TString dataFile = GetDataFileName(run);
  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[ERROR] Cannot open data file: " << dataFile << "\n";
    return;
  }

  TTree* tData = (TTree*) fData->Get(kDataTreeName);
  if (!tData) {
    std::cerr << "[ERROR] Missing tree '" << kDataTreeName << "' in file: " << dataFile << "\n";
    fData->Close();
    delete fData;
    return;
  }

  std::vector<VarInfo> vars = GetVariables();
  std::vector<DpBin> dpBins = BuildDpBins(dpFocusLo, dpFocusHi, nFocusBins);

  for (const auto& dp : dpBins) {
    TString dpCut = Form("(P.gtr.dp>%g) && (P.gtr.dp<%g)", dp.lo, dp.hi);
    TString dataCut = Form("(%s) && (%s)", kDataCutBase, dpCut.Data());

    for (const auto& V : vars) {
      TString hname = Form("h_%s_run%d_%s", V.name.Data(), run, dp.label.Data());
      TH1D* h = new TH1D(hname, "", V.nbins, V.xmin, V.xmax);
      h->Sumw2();

      tData->Draw(Form("%s>>%s", V.expr.Data(), hname.Data()), dataCut.Data(), "goff");

      PeakWindow pw = EstimatePeakWindow(h, 10, 2.0, 4, 25);
      double fitLo = pw.ok ? pw.xLo : V.xmin;
      double fitHi = pw.ok ? pw.xHi : V.xmax;

      FitOut f = FitGaus(h, Form("f_%s_run%d_%s", V.name.Data(), run, dp.label.Data()),
                         fitLo, fitHi);

      AppendCsvRow(csvPath, run, setting, ebeamMeV, dp, V, h->GetEntries(), pw, f);

      if (verbose) {
        std::cout
          << "[INFO] run=" << run
          << " setting=" << setting
          << " dp=" << dp.label
          << " var=" << V.name
          << " entries=" << h->GetEntries()
          << " fit_valid=" << (f.ok ? 1 : 0);

        if (f.ok) {
          std::cout
            << " peak=" << 1000.0 * f.mu
            << " +/- " << 1000.0 * f.muErr
            << " MeV";
        }
        std::cout << "\n";
      }

      delete h;
    }
  }

  fData->Close();
  delete fData;
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------
void WandEmPeakVsDelta(int run = 0,
                       double dpFocusLo = -2.0,
                       double dpFocusHi =  3.0,
                       int nFocusBins = 5,
                       bool loopRuns = true,
                       const char* csvPathC = "heep_check_v1/results/tables/WandEmPeakVsDelta.csv",
                       bool verbose = true) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString csvPath(csvPathC);
  EnsureCsvHeader(csvPath);

  std::vector<int> runs = {
    23839, 23840, 23841, 23842, 23844, 23845,
    23846, 23847, 23848, 23849, 23850, 23851
  };

  if (!loopRuns && run != 0) {
    ProcessOneRun(run, dpFocusLo, dpFocusHi, nFocusBins, csvPath, verbose);
    return;
  }

  for (int r : runs) {
    ProcessOneRun(r, dpFocusLo, dpFocusHi, nFocusBins, csvPath, verbose);
  }
}
