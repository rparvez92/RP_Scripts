// Dry4PassTest_MeasuredOffsets.C
//
// Dry 4-pass offset measurement with tight HMS dp and one common SHMS dp bin.
//
// Run from heep_check_v3:
//   root -l -b -q 'macros/Dry4PassTest_MeasuredOffsets.C()'

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
  std::vector<int> runs;
  double eBeam = 8.5831;
  double hmsP = std::numeric_limits<double>::quiet_NaN();
  double hmsAngleDeg = std::numeric_limits<double>::quiet_NaN();
  double shmsP = std::numeric_limits<double>::quiet_NaN();
  double shmsAngleDeg = std::numeric_limits<double>::quiet_NaN();
};

struct VarInfo {
  TString name;
  TString dataExpr;
  TString simExpr;
  int nbins = 0;
  double xmin = 0.0;
  double xmax = 0.0;
};

struct PeakWindow {
  double lo = 0.0;
  double hi = 0.0;
  double mu = 0.0;
  double sigma = 0.0;
  bool ok = false;
};

struct FitOut {
  double mu = std::numeric_limits<double>::quiet_NaN();
  double muErr = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double sigmaErr = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = -1;
  int status = -999;
  bool ok = false;
};

static std::vector<Setting> GetSettings() {
  return {
    {"A", {23839, 23840, 23841, 23842, 23844, 23845, 23846, 23847, 23848},
     8.5831, 2.28, 42.0, -7.07, 12.465},
    {"B", {23849, 23850, 23851},
     8.5831, 3.75, 29.99, -5.67, 19.34}
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

static TString NormalizeDir(TString dir) {
  if (!dir.EndsWith("/")) dir += "/";
  return dir;
}

static PeakWindow EstimatePeakWindow(const TH1D* h,
                                     int seedHalfWindowBins = 10,
                                     double kSigma = 2.0,
                                     int minHalfWidthBins = 4,
                                     int maxHalfWidthBins = 25) {
  PeakWindow out;
  if (!h || h->GetEntries() <= 0) return out;
  const int nb = h->GetNbinsX();
  const int bMax = h->GetMaximumBin();
  const int bLo = std::max(1, bMax - seedHalfWindowBins);
  const int bHi = std::min(nb, bMax + seedHalfWindowBins);

  double sumW = 0.0, sumWX = 0.0, sumWXX = 0.0;
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
  double sigma = var > 0.0 ? std::sqrt(var) : 0.0;
  const double binW = h->GetXaxis()->GetBinWidth(1);
  if (!(sigma > 0.0)) sigma = seedHalfWindowBins * binW * 0.35;

  double halfW = kSigma * sigma;
  halfW = std::max(halfW, minHalfWidthBins * binW);
  halfW = std::min(halfW, maxHalfWidthBins * binW);
  out.lo = std::max(mu - halfW, h->GetXaxis()->GetXmin());
  out.hi = std::min(mu + halfW, h->GetXaxis()->GetXmax());
  out.mu = mu;
  out.sigma = sigma;
  out.ok = out.hi > out.lo;
  return out;
}

static double IntegralInRange(const TH1D* h, double lo, double hi) {
  if (!h || !(hi > lo)) return 0.0;
  int bLo = h->GetXaxis()->FindBin(lo);
  int bHi = h->GetXaxis()->FindBin(hi);
  bLo = std::max(1, std::min(bLo, h->GetNbinsX()));
  bHi = std::max(1, std::min(bHi, h->GetNbinsX()));
  if (bHi < bLo) std::swap(bLo, bHi);
  return h->Integral(bLo, bHi);
}

static FitOut FitGaus(TH1D* h, const TString& name, double lo, double hi) {
  FitOut out;
  if (!h || !(hi > lo) || h->GetEntries() < 30) return out;
  TF1 f(name, "gaus", lo, hi);
  f.SetParameters(h->GetMaximum(), 0.5 * (lo + hi), 0.25 * (hi - lo));
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
  out.ok = out.status == 0;
  return out;
}

}  // namespace

void Dry4PassTest_MeasuredOffsets(
    const char* dataDirC = "Pass0p1_DataROOTfiles",
    const char* simDirC = "Pass0p1_SimROOTfiles",
    const char* outCsvC = "results/tables/Dry4PassTest/MeasuredOffsets.csv") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  const TString dataDir = NormalizeDir(TString(dataDirC));
  const TString simDir = NormalizeDir(TString(simDirC));
  const TString outCsv(outCsvC);
  gSystem->mkdir(gSystem->DirName(outCsv), kTRUE);

  std::ofstream csv(outCsv.Data(), std::ios::out);
  if (!csv) {
    std::cerr << "[ERROR] Cannot write " << outCsv << "\n";
    return;
  }
  csv << std::fixed << std::setprecision(6);
  csv
    << "setting,run_group,e_beam_GeV,hms_p_GeV,hms_angle_deg,shms_p_GeV,shms_angle_deg,"
    << "dp_label,dp_lo,dp_hi,delta_center,var,"
    << "entriesD,entriesS,norm_win_lo,norm_win_hi,intD_win,intS_win,sim_scale,"
    << "muD_MeV,muD_err_MeV,sigD_MeV,sigD_err_MeV,statusD,chi2D,ndfD,"
    << "muS_MeV,muS_err_MeV,sigS_MeV,sigS_err_MeV,statusS,chi2S,ndfS,"
    << "offset_MeV,offset_err_MeV,fit_valid\n";

  const double dpLo = -2.0;
  const double dpHi = 1.0;
  const double deltaCenter = 0.5 * (dpLo + dpHi);
  const TString dataBaseCut = "(H.gtr.dp>-1) && (H.gtr.dp<1) && (P.gtr.dp>-2) && (P.gtr.dp<1)";
  const TString simBaseCut = "(hsdelta>-1) && (hsdelta<1) && (ssdelta>-2) && (ssdelta<1)";

  for (const auto& setting : GetSettings()) {
    TChain data("T");
    TChain sim("h10");
    TString runGroup;
    for (int run : setting.runs) {
      if (!runGroup.IsNull()) runGroup += ";";
      runGroup += Form("%d", run);
      const TString dataFile = dataDir + Form("coin_replay_production_%d_-1.root", run);
      const TString simFile = simDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
      if (!gSystem->AccessPathName(dataFile)) data.Add(dataFile);
      else std::cerr << "[WARN] Missing data " << dataFile << "\n";
      if (!gSystem->AccessPathName(simFile)) sim.Add(simFile);
      else std::cerr << "[WARN] Missing sim " << simFile << "\n";
    }

    for (const auto& var : GetVars()) {
      const TString hDName = Form("hD_dry4_%s_%s", setting.label.Data(), var.name.Data());
      const TString hSName = Form("hS_dry4_%s_%s", setting.label.Data(), var.name.Data());
      TH1D* hD = new TH1D(hDName, "", var.nbins, var.xmin, var.xmax);
      TH1D* hS = new TH1D(hSName, "", var.nbins, var.xmin, var.xmax);
      hD->Sumw2();
      hS->Sumw2();

      data.Draw(Form("%s>>%s", var.dataExpr.Data(), hDName.Data()), dataBaseCut.Data(), "goff");
      sim.Draw(Form("%s>>%s", var.simExpr.Data(), hSName.Data()),
               Form("Weight*(%s)", simBaseCut.Data()), "goff");

      const PeakWindow pwD = EstimatePeakWindow(hD, 10, 1.0, 4, 25);
      double intD = 0.0, intS = 0.0, simScale = 1.0;
      if (pwD.ok) {
        intD = IntegralInRange(hD, pwD.lo, pwD.hi);
        intS = IntegralInRange(hS, pwD.lo, pwD.hi);
        if (intD > 0.0 && intS > 0.0) {
          simScale = intD / intS;
          hS->Scale(simScale);
        }
      }
      const PeakWindow pwS = EstimatePeakWindow(hS, 10, 1.0, 4, 25);

      const FitOut fD = FitGaus(hD, Form("fD_dry4_%s_%s", setting.label.Data(), var.name.Data()),
                                pwD.ok ? pwD.lo : var.xmin, pwD.ok ? pwD.hi : var.xmax);
      const FitOut fS = FitGaus(hS, Form("fS_dry4_%s_%s", setting.label.Data(), var.name.Data()),
                                pwS.ok ? pwS.lo : (pwD.ok ? pwD.lo : var.xmin),
                                pwS.ok ? pwS.hi : (pwD.ok ? pwD.hi : var.xmax));

      const bool valid = fD.ok && fS.ok;
      const double muD = valid ? fD.mu * 1000.0 : -999.0;
      const double muDErr = valid ? fD.muErr * 1000.0 : -999.0;
      const double sigD = valid ? fD.sigma * 1000.0 : -999.0;
      const double sigDErr = valid ? fD.sigmaErr * 1000.0 : -999.0;
      const double muS = valid ? fS.mu * 1000.0 : -999.0;
      const double muSErr = valid ? fS.muErr * 1000.0 : -999.0;
      const double sigS = valid ? fS.sigma * 1000.0 : -999.0;
      const double sigSErr = valid ? fS.sigmaErr * 1000.0 : -999.0;
      const double offset = valid ? muD - muS : -999.0;
      const double offsetErr = valid ? std::sqrt(muDErr * muDErr + muSErr * muSErr) : -999.0;

      csv
        << setting.label.Data() << ","
        << "\"" << runGroup.Data() << "\","
        << setting.eBeam << ","
        << setting.hmsP << ","
        << setting.hmsAngleDeg << ","
        << setting.shmsP << ","
        << setting.shmsAngleDeg << ","
        << "b1,"
        << dpLo << ","
        << dpHi << ","
        << deltaCenter << ","
        << var.name.Data() << ","
        << hD->GetEntries() << ","
        << hS->GetEntries() << ","
        << (pwD.ok ? pwD.lo : -999.0) << ","
        << (pwD.ok ? pwD.hi : -999.0) << ","
        << intD << ","
        << intS << ","
        << simScale << ","
        << muD << ","
        << muDErr << ","
        << sigD << ","
        << sigDErr << ","
        << fD.status << ","
        << fD.chi2 << ","
        << fD.ndf << ","
        << muS << ","
        << muSErr << ","
        << sigS << ","
        << sigSErr << ","
        << fS.status << ","
        << fS.chi2 << ","
        << fS.ndf << ","
        << offset << ","
        << offsetErr << ","
        << (valid ? 1 : 0) << "\n";

      delete hD;
      delete hS;
    }
  }

  std::cout << "[INFO] Wrote " << outCsv << "\n";
}
