// FitRangeDeterminationData.C
//
// Data-only diagnostic for choosing tight Gaussian fit ranges for 5-pass
// settings C/D/E.  The macro applies the same HMS and SHMS dp cuts used in the
// current 5-pass comparison study, estimates the local peak core, fits only
// mu +/- 1 sigma, and saves both a CSV summary and per-variable PNGs.
//
// Run from heep_check_v3/:
//   root -l -b -q 'macros/FitRangeDeterminationData.C()'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

struct FRSetting {
  TString label;
  int run = 0;
  double shmsP = NAN;
  double dpLo = NAN;
  double dpHi = NAN;
};

struct FRVar {
  TString name;
  TString expr;
  TString title;
  int nbins = 0;
  double xmin = NAN;
  double xmax = NAN;
};

struct FRSeed {
  double peakX = NAN;
  double mu = NAN;
  double sigma = NAN;
  double lo = NAN;
  double hi = NAN;
  int peakBin = -1;
  int seedBinLo = -1;
  int seedBinHi = -1;
  bool ok = false;
};

struct FRFit {
  double mu = NAN;
  double muErr = NAN;
  double sigma = NAN;
  double sigmaErr = NAN;
  double chi2 = NAN;
  int ndf = -1;
  int status = -999;
  bool ok = false;
  TF1* func = nullptr;
};

static std::vector<FRSetting> GetFRSettings() {
  return {
    {"C", 25476, -6.7200, -0.5, 2.5},
    {"E", 25478, -6.6048,  1.5, 4.5},
    {"D", 25477, -6.3840,  4.5, 7.5}
  };
}

static std::vector<FRVar> GetFRVars() {
  return {
    {"W",   "P.kin.primary.W",          "W (GeV)",       100,  0.6,   1.3},
    {"Em",  "H.kin.secondary.emiss",    "Em (GeV)",      100, -0.05,  0.15},
    {"Pmz", "H.kin.secondary.pmiss_z",  "Pmz (GeV/c)",   100, -0.3,   0.3},
    {"Pmy", "H.kin.secondary.pmiss_y",  "Pmy (GeV/c)",   100, -0.3,   0.3}
  };
}

static FRSeed EstimateCoreWindow(const TH1D* h, int seedHalfWindowBins = 10) {
  FRSeed out;
  if (!h || h->GetEntries() <= 0) return out;

  const int nbins = h->GetNbinsX();
  const int peakBin = h->GetMaximumBin();
  if (peakBin < 1 || peakBin > nbins) return out;

  const int bLo = std::max(1, peakBin - seedHalfWindowBins);
  const int bHi = std::min(nbins, peakBin + seedHalfWindowBins);

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
  const double sigma = (var > 0.0 ? std::sqrt(var) : 0.0);
  if (!(sigma > 0.0)) return out;

  out.peakBin = peakBin;
  out.peakX = h->GetXaxis()->GetBinCenter(peakBin);
  out.seedBinLo = bLo;
  out.seedBinHi = bHi;
  out.mu = mu;
  out.sigma = sigma;
  out.lo = std::max(mu - sigma, h->GetXaxis()->GetXmin());
  out.hi = std::min(mu + sigma, h->GetXaxis()->GetXmax());
  out.ok = (out.hi > out.lo);
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

static FRFit FitCore(TH1D* h, const TString& name, const FRSeed& seed) {
  FRFit out;
  if (!h || !seed.ok || h->GetEntries() < 20) return out;

  TF1* f = new TF1(name, "gaus", seed.lo, seed.hi);
  f->SetParameters(h->GetMaximum(), seed.mu, seed.sigma);
  f->SetParLimits(2, 0.2 * seed.sigma, 3.0 * seed.sigma);
  f->SetLineColor(kRed + 1);
  f->SetLineWidth(3);

  TFitResultPtr r = h->Fit(f, "QRSN");
  out.func = f;
  out.status = r.Get() ? r->Status() : -999;
  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf = r->Ndf();
  }
  out.mu = f->GetParameter(1);
  out.muErr = f->GetParError(1);
  out.sigma = std::abs(f->GetParameter(2));
  out.sigmaErr = f->GetParError(2);
  out.ok = (out.status == 0);
  return out;
}

static void DrawFitRangePng(TH1D* h,
                            const FRSetting& setting,
                            const FRVar& var,
                            const FRSeed& seed,
                            const FRFit& fit,
                            const TString& outPng) {
  TCanvas c(Form("c_fr_%d_%s", setting.run, var.name.Data()),
            Form("FR Data %d %s", setting.run, var.name.Data()),
            1100, 800);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.12);
  c.SetTopMargin(0.09);

  h->SetStats(0);
  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
  h->SetMarkerColor(kBlue + 1);
  h->GetXaxis()->SetTitle(var.title);
  h->GetYaxis()->SetTitle("Counts");
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->Draw("E");

  const double yMax = std::max(1.0, h->GetMaximum() * 1.18);
  h->SetMaximum(yMax);

  if (seed.ok) {
    TBox box(seed.lo, 0.0, seed.hi, yMax * 0.96);
    box.SetFillColorAlpha(kOrange + 1, 0.22);
    box.SetLineColor(kOrange + 7);
    box.SetLineWidth(2);
    box.Draw("same");

    TLine lLo(seed.lo, 0.0, seed.lo, yMax * 0.96);
    TLine lHi(seed.hi, 0.0, seed.hi, yMax * 0.96);
    lLo.SetLineColor(kOrange + 7);
    lHi.SetLineColor(kOrange + 7);
    lLo.SetLineStyle(2);
    lHi.SetLineStyle(2);
    lLo.SetLineWidth(2);
    lHi.SetLineWidth(2);
    lLo.Draw("same");
    lHi.Draw("same");
  }

  h->Draw("E same");
  if (fit.func) fit.func->Draw("same");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.034);
  latex.DrawLatex(0.13, 0.945,
                  Form("Data run %d, setting %s, %s;  H.gtr.dp #in [-1,1], P.gtr.dp #in [%.1f,%.1f]",
                       setting.run, setting.label.Data(), var.name.Data(), setting.dpLo, setting.dpHi));

  latex.SetTextSize(0.031);
  if (seed.ok) {
    latex.DrawLatex(0.58, 0.84,
                    Form("seed #mu = %.4f, #sigma = %.4f", seed.mu, seed.sigma));
    latex.DrawLatex(0.58, 0.79,
                    Form("fit range: %.4f to %.4f", seed.lo, seed.hi));
  }
  if (fit.func) {
    latex.DrawLatex(0.58, 0.72,
                    Form("fit #mu = %.4f #pm %.4f", fit.mu, fit.muErr));
    latex.DrawLatex(0.58, 0.67,
                    Form("fit #sigma = %.4f #pm %.4f", fit.sigma, fit.sigmaErr));
    latex.DrawLatex(0.58, 0.62,
                    Form("status = %d, #chi^{2}/ndf = %.2f/%d", fit.status, fit.chi2, fit.ndf));
  }

  c.RedrawAxis();
  c.SaveAs(outPng);
}

void FitRangeDeterminationData(const char* dataDir = "Pass0p1_DataROOTfiles",
                               const char* outCsv = "results/tables/FitRangeDeterminationData_5pass.csv",
                               const char* outPngDir = "results/PNGs/FitRangeDetermination") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(4);

  TString dataDirS(dataDir);
  if (!dataDirS.EndsWith("/")) dataDirS += "/";

  const TString csvPath(outCsv);
  const TString pngDir(outPngDir);
  const TString csvParent = gSystem->DirName(csvPath);
  if (!csvParent.IsNull()) gSystem->mkdir(csvParent, kTRUE);
  gSystem->mkdir(pngDir, kTRUE);

  std::ofstream csv(csvPath.Data(), std::ios::out);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot write " << csvPath << "\n";
    return;
  }
  csv
    << "setting,run,shms_p_GeV,dp_lo,dp_hi,var,"
    << "hist_nbins,hist_xmin,hist_xmax,bin_width,"
    << "entries_total,entries_fit_range,"
    << "peak_x,seed_mu,seed_sigma,fit_range_lo,fit_range_hi,"
    << "seed_mu_MeV,seed_sigma_MeV,fit_range_lo_MeV,fit_range_hi_MeV,"
    << "fit_mu,fit_mu_err,fit_sigma,fit_sigma_err,"
    << "fit_mu_MeV,fit_mu_err_MeV,fit_sigma_MeV,fit_sigma_err_MeV,"
    << "status,chi2,ndf,chi2_ndf,fit_valid,png"
    << "\n";
  csv << std::fixed << std::setprecision(6);

  const auto settings = GetFRSettings();
  const auto vars = GetFRVars();
  const TString baseCut = "(H.gtr.dp>-1) && (H.gtr.dp<1)";

  for (const auto& setting : settings) {
    const TString fileName = dataDirS + Form("coin_replay_production_%d_-1.root", setting.run);
    if (gSystem->AccessPathName(fileName)) {
      std::cerr << "[WARN] Missing data file: " << fileName << "\n";
      continue;
    }

    TChain chain("T");
    chain.Add(fileName);
    const TString cut = Form("(%s) && (P.gtr.dp>%g) && (P.gtr.dp<%g)",
                             baseCut.Data(), setting.dpLo, setting.dpHi);

    for (const auto& var : vars) {
      const TString hName = Form("h_fr_%d_%s", setting.run, var.name.Data());
      TH1D* h = new TH1D(hName, "", var.nbins, var.xmin, var.xmax);
      h->Sumw2();

      chain.Draw(Form("%s>>%s", var.expr.Data(), hName.Data()), cut.Data(), "goff");

      const FRSeed seed = EstimateCoreWindow(h, 10);
      const FRFit fit = FitCore(h,
                                Form("f_fr_%d_%s", setting.run, var.name.Data()),
                                seed);
      const double entriesRange = seed.ok ? IntegralInRange(h, seed.lo, seed.hi) : 0.0;
      const double chi2Ndf = (fit.ndf > 0 ? fit.chi2 / fit.ndf : NAN);
      const double binWidth = h->GetXaxis()->GetBinWidth(1);
      const TString pngPath = pngDir + Form("/FR_Data_%d_%s.png", setting.run, var.name.Data());

      DrawFitRangePng(h, setting, var, seed, fit, pngPath);

      csv
        << setting.label.Data() << ","
        << setting.run << ","
        << setting.shmsP << ","
        << setting.dpLo << ","
        << setting.dpHi << ","
        << var.name.Data() << ","
        << var.nbins << ","
        << var.xmin << ","
        << var.xmax << ","
        << binWidth << ","
        << h->GetEntries() << ","
        << entriesRange << ","
        << seed.peakX << ","
        << seed.mu << ","
        << seed.sigma << ","
        << seed.lo << ","
        << seed.hi << ","
        << seed.mu * 1000.0 << ","
        << seed.sigma * 1000.0 << ","
        << seed.lo * 1000.0 << ","
        << seed.hi * 1000.0 << ","
        << fit.mu << ","
        << fit.muErr << ","
        << fit.sigma << ","
        << fit.sigmaErr << ","
        << fit.mu * 1000.0 << ","
        << fit.muErr * 1000.0 << ","
        << fit.sigma * 1000.0 << ","
        << fit.sigmaErr * 1000.0 << ","
        << fit.status << ","
        << fit.chi2 << ","
        << fit.ndf << ","
        << chi2Ndf << ","
        << (fit.ok ? 1 : 0) << ","
        << pngPath.Data()
        << "\n";

      delete h;
    }
  }

  std::cout << "[INFO] Wrote " << csvPath << "\n";
  std::cout << "[INFO] Wrote PNGs under " << pngDir << "\n";
}
