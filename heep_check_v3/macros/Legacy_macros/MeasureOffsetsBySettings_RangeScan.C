// MeasureOffsetsBySettings_RangeScan.C
//
// Isolated 5-pass C/D/E range-scan version of MeasuredOffsetsBySetting.C.
// It scans Gaussian fit windows of local_core_mu +/- k * local_core_sigma
// with the same local-core clamp and sim-normalization convention used by
// MeasuredOffsetsBySetting.C.
//
// Run from heep_check_v3/:
//   root -l -b -q 'macros/MeasureOffsetsBySettings_RangeScan.C()'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TChain.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

struct Setting {
  TString label;
  int run = 0;
  double eBeam = 10.6716;
  double hmsP = 4.72;
  double hmsAngleDeg = 26.395;
  double shmsP = std::numeric_limits<double>::quiet_NaN();
  double shmsAngleDeg = 18.56;
  double dpLo = std::numeric_limits<double>::quiet_NaN();
  double dpHi = std::numeric_limits<double>::quiet_NaN();
};

struct VarInfo {
  TString name;
  TString dataExpr;
  TString simExpr;
  int nbins = 0;
  double xmin = std::numeric_limits<double>::quiet_NaN();
  double xmax = std::numeric_limits<double>::quiet_NaN();
};

struct LocalCore {
  int peakBin = -1;
  int seedBinLo = -1;
  int seedBinHi = -1;
  double peakX = std::numeric_limits<double>::quiet_NaN();
  double mu = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double seedLo = std::numeric_limits<double>::quiet_NaN();
  double seedHi = std::numeric_limits<double>::quiet_NaN();
  bool ok = false;
};

struct FitResult {
  double lo = std::numeric_limits<double>::quiet_NaN();
  double hi = std::numeric_limits<double>::quiet_NaN();
  double entriesInRange = 0.0;
  double mu = std::numeric_limits<double>::quiet_NaN();
  double muErr = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double sigmaErr = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = -1;
  int status = -999;
  bool ok = false;
};

struct FitWindow {
  double lo = std::numeric_limits<double>::quiet_NaN();
  double hi = std::numeric_limits<double>::quiet_NaN();
  bool ok = false;
};

static std::vector<Setting> GetSettings() {
  return {
    {"C", 25476, 10.6716, 4.72, 26.395, -6.7200, 18.56, -0.5, 2.5},
    {"E", 25478, 10.6716, 4.72, 26.395, -6.6048, 18.56,  1.5, 4.5},
    {"D", 25477, 10.6716, 4.72, 26.395, -6.3840, 18.56,  4.5, 7.5}
  };
}

static std::vector<VarInfo> GetVars() {
  return {
    {"W",   "P.kin.primary.W",         "W_recon",   100,  0.6,   1.3},
    {"Em",  "H.kin.secondary.emiss",   "Em_recon",  100, -0.05,  0.15},
    {"Pmz", "H.kin.secondary.pmiss_z", "Pmz_recon", 100, -0.3,   0.3},
    {"Pmy", "H.kin.secondary.pmiss_y", "Pmy_recon", 100, -0.3,   0.3}
  };
}

static std::vector<double> GetKSigmas() {
  return {0.50, 0.75, 1.00, 1.50, 2.00};
}

static TString NormalizeDir(TString dir) {
  if (!dir.EndsWith("/")) dir += "/";
  return dir;
}

static LocalCore EstimateLocalCore(const TH1D* h, int seedHalfWindowBins = 10) {
  LocalCore out;
  if (!h || h->GetEntries() <= 0) return out;

  const int nb = h->GetNbinsX();
  const int bMax = h->GetMaximumBin();
  if (bMax < 1 || bMax > nb) return out;

  const int bLo = std::max(1, bMax - seedHalfWindowBins);
  const int bHi = std::min(nb, bMax + seedHalfWindowBins);

  double sumW = 0.0;
  double sumWX = 0.0;
  double sumWXX = 0.0;
  for (int b = bLo; b <= bHi; ++b) {
    const double w = h->GetBinContent(b);
    if (w <= 0.0) continue;
    const double x = h->GetXaxis()->GetBinCenter(b);
    sumW += w;
    sumWX += w * x;
    sumWXX += w * x * x;
  }
  if (sumW <= 0.0) return out;

  const double mu = sumWX / sumW;
  const double var = sumWXX / sumW - mu * mu;
  if (!(var > 0.0)) return out;

  out.peakBin = bMax;
  out.seedBinLo = bLo;
  out.seedBinHi = bHi;
  out.peakX = h->GetXaxis()->GetBinCenter(bMax);
  out.mu = mu;
  out.sigma = std::sqrt(var);
  out.seedLo = h->GetXaxis()->GetBinLowEdge(bLo);
  out.seedHi = h->GetXaxis()->GetBinUpEdge(bHi);
  out.ok = (out.sigma > 0.0);
  return out;
}

static double IntegralInRange(const TH1D* h, double xLo, double xHi) {
  if (!h || !(xHi > xLo)) return 0.0;
  int bLo = h->GetXaxis()->FindBin(xLo);
  int bHi = h->GetXaxis()->FindBin(xHi);
  bLo = std::max(1, std::min(bLo, h->GetNbinsX()));
  bHi = std::max(1, std::min(bHi, h->GetNbinsX()));
  if (bHi < bLo) std::swap(bLo, bHi);
  return h->Integral(bLo, bHi);
}

static FitWindow MakeClampedCoreWindow(const TH1D* h,
                                       const LocalCore& core,
                                       double kSigma,
                                       int minHalfWidthBins = 2,
                                       int maxHalfWidthBins = 20) {
  FitWindow out;
  if (!h || !core.ok || !(kSigma > 0.0)) return out;

  const double binW = h->GetXaxis()->GetBinWidth(1);
  double halfW = kSigma * core.sigma;
  halfW = std::max(halfW, minHalfWidthBins * binW);
  halfW = std::min(halfW, maxHalfWidthBins * binW);

  out.lo = std::max(core.mu - halfW, h->GetXaxis()->GetXmin());
  out.hi = std::min(core.mu + halfW, h->GetXaxis()->GetXmax());
  out.ok = (out.hi > out.lo);
  return out;
}

static FitResult FitWithCoreRange(TH1D* h,
                                  const TString& fName,
                                  const LocalCore& core,
                                  double kSigma,
                                  Color_t color) {
  FitResult out;
  if (!h || !core.ok) return out;

  const FitWindow win = MakeClampedCoreWindow(h, core, kSigma, 2, 20);
  if (!win.ok) return out;

  out.lo = win.lo;
  out.hi = win.hi;
  out.entriesInRange = IntegralInRange(h, out.lo, out.hi);
  if (!(out.hi > out.lo) || h->GetEntries() < 30) return out;

  TF1 f(fName, "gaus", out.lo, out.hi);
  f.SetParameters(h->GetMaximum(), core.mu, core.sigma);
  f.SetParLimits(2, 0.05 * core.sigma, 5.0 * core.sigma);
  f.SetLineColor(color);

  TFitResultPtr r = h->Fit(&f, "QRSN");
  out.status = r.Get() ? r->Status() : -999;
  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf = r->Ndf();
  }
  out.mu = f.GetParameter(1);
  out.muErr = f.GetParError(1);
  out.sigma = std::abs(f.GetParameter(2));
  out.sigmaErr = f.GetParError(2);
  out.ok = (out.status == 0);
  return out;
}

static void WriteMeasuredHeader(std::ofstream& out) {
  out
    << "fit_k_sigma,setting,run,e_beam_GeV,hms_p_GeV,hms_angle_deg,"
    << "shms_p_GeV,shms_angle_deg,dp_label,dp_lo,dp_hi,bin_center,var,"
    << "norm_win_lo_MeV,norm_win_hi_MeV,intD_win,intS_win,sim_scale,"
    << "entriesD,entriesS,entriesD_fit_range,entriesS_fit_range,"
    << "seed_muD_MeV,seed_sigD_MeV,seed_loD_MeV,seed_hiD_MeV,"
    << "fit_loD_MeV,fit_hiD_MeV,muD_MeV,muD_err_MeV,sigD_MeV,sigD_err_MeV,statusD,chi2D,ndfD,"
    << "seed_muS_MeV,seed_sigS_MeV,seed_loS_MeV,seed_hiS_MeV,"
    << "fit_loS_MeV,fit_hiS_MeV,muS_MeV,muS_err_MeV,sigS_MeV,sigS_err_MeV,statusS,chi2S,ndfS,"
    << "offset_MeV,offset_err_MeV,fit_valid"
    << "\n";
}

static void WriteInputHeader(std::ofstream& out) {
  out
    << "fit_k_sigma,setting,shms_p_GeV,bin_center,dp_lo,dp_hi,variable,"
    << "mu,mu_err,sigma,sigma_err,fit_valid,statusD,statusS"
    << "\n";
}

}  // namespace

void MeasureOffsetsBySettings_RangeScan(
    const char* dataDirC = "Pass0p1_DataROOTfiles",
    const char* simDirC = "Pass0p1_SimROOTfiles",
    const char* outMeasuredCsvC = "results/tables/RangeScan/MeasuredOffsetsBySetting_RangeScan.csv",
    const char* outInputCsvC = "results/tables/RangeScan/TestingRadwan_RangeScan_input.csv") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  const TString dataDir = NormalizeDir(TString(dataDirC));
  const TString simDir = NormalizeDir(TString(simDirC));
  const TString outMeasuredCsv(outMeasuredCsvC);
  const TString outInputCsv(outInputCsvC);
  gSystem->mkdir(gSystem->DirName(outMeasuredCsv), kTRUE);
  gSystem->mkdir(gSystem->DirName(outInputCsv), kTRUE);

  std::ofstream measured(outMeasuredCsv.Data(), std::ios::out);
  std::ofstream input(outInputCsv.Data(), std::ios::out);
  if (!measured.is_open() || !input.is_open()) {
    std::cerr << "[ERROR] Could not open RangeScan output CSVs.\n";
    return;
  }
  measured << std::fixed << std::setprecision(6);
  input << std::fixed << std::setprecision(6);
  WriteMeasuredHeader(measured);
  WriteInputHeader(input);

  const auto settings = GetSettings();
  const auto vars = GetVars();
  const auto kSigmas = GetKSigmas();

  for (const auto& setting : settings) {
    const TString dataFile = dataDir + Form("coin_replay_production_%d_-1.root", setting.run);
    const TString simFile = simDir + Form("recon_hcana_coin_replay_production_%d_-1.root", setting.run);
    if (gSystem->AccessPathName(dataFile)) {
      std::cerr << "[WARN] Missing data file: " << dataFile << "\n";
      continue;
    }
    if (gSystem->AccessPathName(simFile)) {
      std::cerr << "[WARN] Missing sim file: " << simFile << "\n";
      continue;
    }

    TChain data("T");
    TChain sim("h10");
    data.Add(dataFile);
    sim.Add(simFile);

    const TString dataCut = Form("(H.gtr.dp>-1) && (H.gtr.dp<1) && (P.gtr.dp>%g) && (P.gtr.dp<%g)",
                                 setting.dpLo, setting.dpHi);
    const TString simCut = Form("(hsdelta>-1) && (hsdelta<1) && (ssdelta>%g) && (ssdelta<%g)",
                                setting.dpLo, setting.dpHi);

    for (const auto& var : vars) {
      const TString hDName = Form("hD_rs_%s_%s", setting.label.Data(), var.name.Data());
      const TString hSName = Form("hS_rs_%s_%s", setting.label.Data(), var.name.Data());
      TH1D* hD = new TH1D(hDName, "", var.nbins, var.xmin, var.xmax);
      TH1D* hS = new TH1D(hSName, "", var.nbins, var.xmin, var.xmax);
      hD->Sumw2();
      hS->Sumw2();

      data.Draw(Form("%s>>%s", var.dataExpr.Data(), hDName.Data()), dataCut.Data(), "goff");
      sim.Draw(Form("%s>>%s", var.simExpr.Data(), hSName.Data()),
               Form("Weight*(%s)", simCut.Data()), "goff");

      const LocalCore coreD = EstimateLocalCore(hD, 10);

      for (double kSigma : kSigmas) {
        const FitWindow dataWin = MakeClampedCoreWindow(hD, coreD, kSigma, 2, 20);
        double intDwin = 0.0;
        double intSwin = 0.0;
        double simScale = 1.0;

        TH1D* hSScaled = static_cast<TH1D*>(hS->Clone(Form("%s_scaled_k%.2f", hSName.Data(), kSigma)));
        hSScaled->SetDirectory(nullptr);
        if (dataWin.ok) {
          intDwin = IntegralInRange(hD, dataWin.lo, dataWin.hi);
          intSwin = IntegralInRange(hSScaled, dataWin.lo, dataWin.hi);
          if (intDwin > 0.0 && intSwin > 0.0) {
            simScale = intDwin / intSwin;
            hSScaled->Scale(simScale);
          }
        }

        const LocalCore coreS = EstimateLocalCore(hSScaled, 10);
        FitResult fD = FitWithCoreRange(hD,
                                        Form("fD_rs_%s_%s_k%.2f", setting.label.Data(), var.name.Data(), kSigma),
                                        coreD, kSigma, kBlack);
        FitResult fS = FitWithCoreRange(hSScaled,
                                        Form("fS_rs_%s_%s_k%.2f", setting.label.Data(), var.name.Data(), kSigma),
                                        coreS, kSigma, kRed + 1);
        const bool fitValid = fD.ok && fS.ok;
        const double offset = fitValid ? (fD.mu - fS.mu) * 1000.0 : -999.0;
        const double offsetErr = fitValid
          ? std::sqrt(fD.muErr * fD.muErr + fS.muErr * fS.muErr) * 1000.0
          : -999.0;

        measured
          << kSigma << ","
          << setting.label.Data() << ","
          << setting.run << ","
          << setting.eBeam << ","
          << setting.hmsP << ","
          << setting.hmsAngleDeg << ","
          << setting.shmsP << ","
          << setting.shmsAngleDeg << ","
          << "b1,"
          << setting.dpLo << ","
          << setting.dpHi << ","
          << 0.5 * (setting.dpLo + setting.dpHi) << ","
          << var.name.Data() << ","
          << (dataWin.ok ? dataWin.lo * 1000.0 : -999.0) << ","
          << (dataWin.ok ? dataWin.hi * 1000.0 : -999.0) << ","
          << intDwin << ","
          << intSwin << ","
          << simScale << ","
          << hD->GetEntries() << ","
          << hS->GetEntries() << ","
          << fD.entriesInRange << ","
          << fS.entriesInRange << ","
          << coreD.mu * 1000.0 << ","
          << coreD.sigma * 1000.0 << ","
          << coreD.seedLo * 1000.0 << ","
          << coreD.seedHi * 1000.0 << ","
          << fD.lo * 1000.0 << ","
          << fD.hi * 1000.0 << ","
          << fD.mu * 1000.0 << ","
          << fD.muErr * 1000.0 << ","
          << fD.sigma * 1000.0 << ","
          << fD.sigmaErr * 1000.0 << ","
          << fD.status << ","
          << fD.chi2 << ","
          << fD.ndf << ","
          << coreS.mu * 1000.0 << ","
          << coreS.sigma * 1000.0 << ","
          << coreS.seedLo * 1000.0 << ","
          << coreS.seedHi * 1000.0 << ","
          << fS.lo * 1000.0 << ","
          << fS.hi * 1000.0 << ","
          << fS.mu * 1000.0 << ","
          << fS.muErr * 1000.0 << ","
          << fS.sigma * 1000.0 << ","
          << fS.sigmaErr * 1000.0 << ","
          << fS.status << ","
          << fS.chi2 << ","
          << fS.ndf << ","
          << offset << ","
          << offsetErr << ","
          << (fitValid ? 1 : 0)
          << "\n";

        input
          << kSigma << ","
          << setting.label.Data() << ","
          << setting.shmsP << ","
          << 0.5 * (setting.dpLo + setting.dpHi) << ","
          << setting.dpLo << ","
          << setting.dpHi << ","
          << var.name.Data() << ","
          << offset << ","
          << offsetErr << ","
          << (fitValid ? fD.sigma * 1000.0 : -999.0) << ","
          << (fitValid ? fD.sigmaErr * 1000.0 : -999.0) << ","
          << (fitValid ? 1 : 0) << ","
          << fD.status << ","
          << fS.status
          << "\n";

        delete hSScaled;
      }

      delete hD;
      delete hS;
    }
  }

  std::cout << "[INFO] Wrote " << outMeasuredCsv << "\n";
  std::cout << "[INFO] Wrote " << outInputCsv << "\n";
}
