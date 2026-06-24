// MeasuredOffsetsBySetting.C
//
// Purpose:
//   For each hardcoded setting, chain the corresponding Data and Sim ROOT files,
//   build SHMS dp-bin histograms, fit Gaussian means for
//     - W
//     - Em
//     - Pmz
//     - Pmy
//   and measure the setting-wise offset
//       offset = mu_data - mu_sim
//
// Settings:
//   Setting A: runs 23839-23848 except 23843
//   Setting B: runs 23849-23851
//
// Per setting, the macro builds 6 dp-bin canvases:
//   - 1 full bin over SHMS dp in [-10, 20]
//   - 5 narrow bins spanning [dpFocusLo, dpFocusHi]
//
// Each canvas is a 2x2 layout with Data vs Sim overlays for:
//   upper left:  W
//   upper right: Em
//   lower left:  Pmz
//   lower right: Pmy
//
// Analysis notes:
//   - Data ROOT files are read from the Pass0p1 replay symlink:
//       RP_Scripts/heep_check_v2/Pass0p1_ROOTfiles/
//   - Sim is filled from the h10 tree using event-by-event Weight
//   - Sim is normalized to Data inside an adaptive peak window
//   - Data and Sim are fit with Gaussians when possible
//   - The measured offset is:
//       offset = mu_data - mu_sim
//   - All physics-valued CSV columns are written in MeV
//   - The offset uncertainty is defined from the propagated Gaussian mean-fit
//     errors for Data and Sim:
//       offset_err_MeV = sqrt(muD_err_MeV^2 + muS_err_MeV^2)
//   - The fitted Gaussian widths and their errors are still retained for
//     diagnostics, in MeV
//
// Output PNGs:
//   results/PNGs/MeasuredOffsetsBySetting/offsets_Setting<AB>_<binlabel>bin.png
//
// Output CSV:
//   results/tables/MeasuredOffsetsBySetting.csv
//
// Run from:
//   RP_Scripts/heep_check_v2/
//
// Example:
//   root -l -b -q 'macros/MeasuredOffsetsBySetting.C+( -2, 3, 5 )'

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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
#include "TROOT.h"

// -----------------------------------------------------------------------------
// Directories
// -----------------------------------------------------------------------------
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check_v2/Pass0p1_ROOTfiles/";
static const TString kSimDir  = "/home/cdaq/users/rparvez/simc_gfortran/worksim/heep_check_v2/";

static TString GetDataFileName(int run) {
  return kDataDir + Form("coin_replay_production_%d_-1.root", run);
}
static TString GetSimFileName(int run) {
  return kSimDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
}

// Trees
static const TString kDataTreeName = "T";
static const TString kSimTreeName  = "h10";

// Baseline cuts
static const char* kDataCutBase = "(H.gtr.dp>-8) && (H.gtr.dp<8)";
static const char* kSimCutBase  = "(hsdelta>-8) && (hsdelta<8)";

// Sim weight branch
static const char* kSimWeightBranch = "Weight";

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

struct DpBin {
  int idx;
  TString label;
  double lo;
  double hi;
};

static std::vector<DpBin> BuildDpBins(double dpFocusLo, double dpFocusHi, int nFocusBins) {
  std::vector<DpBin> bins;
  bins.push_back({0, "full", -10.0, 20.0});

  if (nFocusBins < 1) nFocusBins = 1;
  const double step = (dpFocusHi - dpFocusLo) / nFocusBins;
  for (int i = 0; i < nFocusBins; ++i) {
    bins.push_back({i + 1,
                    Form("b%d", i + 1),
                    dpFocusLo + i * step,
                    dpFocusLo + (i + 1) * step});
  }
  return bins;
}

struct SettingInfo {
  TString label;
  std::vector<int> runs;
};

static std::vector<SettingInfo> GetSettings() {
  return {
    {"A", {23839, 23840, 23841, 23842, 23844, 23845, 23846, 23847, 23848}},
    {"B", {23849, 23850, 23851}}
  };
}

static void StyleHists(TH1D* hData, TH1D* hSim) {
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);

  hSim->SetMarkerStyle(24);
  hSim->SetMarkerSize(0.7);
  hSim->SetLineColor(kRed + 1);

  hData->GetYaxis()->SetTitle("Counts");
  hData->GetYaxis()->SetTitleOffset(1.15);
  hData->GetXaxis()->SetTitleOffset(1.05);

  hData->GetXaxis()->SetTitleSize(0.058);
  hData->GetXaxis()->SetLabelSize(0.048);
  hData->GetYaxis()->SetTitleSize(0.058);
  hData->GetYaxis()->SetLabelSize(0.048);
}

struct PeakWindow {
  double xLo = 0;
  double xHi = 0;
  double muSeed = 0;
  double sigSeed = 0;
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

  const int bMax = h->GetMaximumBin();
  if (bMax < 1 || bMax > nb) return pw;

  const int bLo = std::max(1, bMax - seedHalfWindowBins);
  const int bHi = std::min(nb, bMax + seedHalfWindowBins);

  double sumW = 0.0, sumWX = 0.0, sumWXX = 0.0;
  for (int b = bLo; b <= bHi; ++b) {
    const double w = h->GetBinContent(b);
    const double x = h->GetXaxis()->GetBinCenter(b);
    if (w <= 0) continue;
    sumW += w;
    sumWX += w * x;
    sumWXX += w * x * x;
  }
  if (sumW <= 0) return pw;

  const double mu = sumWX / sumW;
  const double var = sumWXX / sumW - mu * mu;
  double sig = (var > 0 ? std::sqrt(var) : 0.0);

  const double binW = h->GetXaxis()->GetBinWidth(1);
  if (!(sig > 0)) sig = seedHalfWindowBins * binW * 0.35;

  double halfW = kSigma * sig;
  halfW = std::max(halfW, minHalfWidthBins * binW);
  halfW = std::min(halfW, maxHalfWidthBins * binW);

  double xLo = mu - halfW;
  double xHi = mu + halfW;

  xLo = std::max(xLo, h->GetXaxis()->GetXmin());
  xHi = std::min(xHi, h->GetXaxis()->GetXmax());

  pw.xLo = xLo;
  pw.xHi = xHi;
  pw.muSeed = mu;
  pw.sigSeed = sig;
  pw.ok = (xHi > xLo);
  return pw;
}

static double IntegralInXRange(const TH1D* h, double xLo, double xHi) {
  if (!h || xHi <= xLo) return 0.0;
  int b1 = h->GetXaxis()->FindBin(xLo);
  int b2 = h->GetXaxis()->FindBin(xHi);
  b1 = std::max(1, std::min(b1, h->GetNbinsX()));
  b2 = std::max(1, std::min(b2, h->GetNbinsX()));
  if (b2 < b1) std::swap(b1, b2);
  return h->Integral(b1, b2);
}

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
  if (!h || !(xHi > xLo) || h->GetEntries() < 30) return out;

  TF1* f = new TF1(fname, "gaus", xLo, xHi);
  f->SetLineColor(lineColor);
  f->SetLineWidth(2);

  TFitResultPtr r = h->Fit(f, "QRSN");
  out.f = f;
  out.status = r.Get() ? r->Status() : -999;

  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf = r->Ndf();
  }

  out.mu = f->GetParameter(1);
  out.muErr = f->GetParError(1);
  out.sig = f->GetParameter(2);
  out.sigErr = f->GetParError(2);
  out.ok = (out.status == 0);
  return out;
}

static void EnsureCsvHeader(const TString& csvPath) {
  if (!gSystem->AccessPathName(csvPath)) return;

  TString parent = gSystem->DirName(csvPath);
  if (!parent.IsNull() && gSystem->AccessPathName(parent)) {
    gSystem->mkdir(parent, kTRUE);
  }

  std::ofstream csv(csvPath.Data(), std::ios::out);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot create CSV: " << csvPath << "\n";
    return;
  }

  csv
    << "setting,dp_idx,dp_label,dp_lo,dp_hi,"
    << "var,"
    << "norm_win_lo,norm_win_hi,"
    << "intD_win,intS_win,sim_scale,"
    << "entriesD,entriesS,"
    << "muD_MeV,muD_err_MeV,sigD_MeV,sigD_err_MeV,statusD,chi2D,ndfD,"
    << "muS_MeV,muS_err_MeV,sigS_MeV,sigS_err_MeV,statusS,chi2S,ndfS,"
    << "offset_MeV,offset_err_MeV,fit_valid"
    << "\n";
}

static void AppendCsvRow(const TString& csvPath,
                         const TString& setting, const DpBin& dp, const VarInfo& V,
                         const PeakWindow& pw,
                         double intDwin, double intSwin, double simScale,
                         double entriesD, double entriesS,
                         const FitOut& fD, const FitOut& fS,
                         double muDMeV, double muDErrMeV,
                         double sigDMeV, double sigDErrMeV,
                         double muSMeV, double muSErrMeV,
                         double sigSMeV, double sigSErrMeV,
                         double offsetMeV, double offsetErrMeV,
                         bool fitValid) {
  std::ofstream csv(csvPath.Data(), std::ios::out | std::ios::app);
  if (!csv.is_open()) {
    std::cerr << "[ERROR] Cannot append CSV: " << csvPath << "\n";
    return;
  }

  csv.setf(std::ios::fixed);
  csv << std::setprecision(6);

  csv
    << setting.Data() << ","
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
    << muDMeV << ","
    << muDErrMeV << ","
    << sigDMeV << ","
    << sigDErrMeV << ","
    << fD.status << ","
    << fD.chi2 << ","
    << fD.ndf << ","
    << muSMeV << ","
    << muSErrMeV << ","
    << sigSMeV << ","
    << sigSErrMeV << ","
    << fS.status << ","
    << fS.chi2 << ","
    << fS.ndf << ","
    << offsetMeV << ","
    << offsetErrMeV << ","
    << (fitValid ? 1 : 0)
    << "\n";
}

static bool BuildChainsForSetting(const SettingInfo& setting, TChain& dataChain, TChain& simChain) {
  bool addedAny = false;
  for (int run : setting.runs) {
    const TString dataFile = GetDataFileName(run);
    const TString simFile  = GetSimFileName(run);

    if (gSystem->AccessPathName(dataFile)) {
      std::cerr << "[WARN] Missing data file for run " << run << ": " << dataFile << "\n";
      continue;
    }
    if (gSystem->AccessPathName(simFile)) {
      std::cerr << "[WARN] Missing sim file for run " << run << ": " << simFile << "\n";
      continue;
    }

    dataChain.Add(dataFile);
    simChain.Add(simFile);
    addedAny = true;
  }
  return addedAny;
}

static void ProcessOneSetting(const SettingInfo& setting,
                              double dpFocusLo, double dpFocusHi, int nFocusBins,
                              const TString& outDir,
                              const TString& csvPath) {
  TChain tData(kDataTreeName);
  TChain tSim(kSimTreeName);
  if (!BuildChainsForSetting(setting, tData, tSim)) {
    std::cerr << "[ERROR] No valid data/sim file pairs found for Setting " << setting.label << "\n";
    return;
  }

  std::vector<VarInfo> vars = GetVariables();
  std::vector<DpBin> dpBins = BuildDpBins(dpFocusLo, dpFocusHi, nFocusBins);

  for (const auto& dp : dpBins) {
    TCanvas* c = new TCanvas(Form("c_heep_setting%s_%s", setting.label.Data(), dp.label.Data()),
                             Form("HEEP offsets Setting %s dp %s", setting.label.Data(), dp.label.Data()),
                             1600, 1000);
    c->Divide(2, 2);

    const TString dataDpCut = Form("(P.gtr.dp>%g) && (P.gtr.dp<%g)", dp.lo, dp.hi);
    const TString simDpCut  = Form("(ssdelta>%g) && (ssdelta<%g)", dp.lo, dp.hi);
    const TString dataCut = Form("(%s) && (%s)", kDataCutBase, dataDpCut.Data());
    const TString simCut  = Form("(%s) && (%s)", kSimCutBase,  simDpCut.Data());

    for (int i = 0; i < (int)vars.size(); ++i) {
      const VarInfo& V = vars[i];
      c->cd(i + 1);
      gPad->SetLeftMargin(0.12);
      gPad->SetRightMargin(0.05);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.08);

      const TString hDname = Form("hD_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data());
      const TString hSname = Form("hS_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data());

      TH1D* hD = new TH1D(hDname, "", V.nbins, V.xmin, V.xmax);
      TH1D* hS = new TH1D(hSname, "", V.nbins, V.xmin, V.xmax);
      hD->Sumw2();
      hS->Sumw2();

      tData.Draw(Form("%s>>%s", V.data_expr.Data(), hDname.Data()), dataCut.Data(), "goff");
      const TString simSel = Form("%s*(%s)", kSimWeightBranch, simCut.Data());
      tSim.Draw(Form("%s>>%s", V.sim_expr.Data(), hSname.Data()), simSel.Data(), "goff");

      const PeakWindow pw = EstimatePeakWindow(hD, 10, 2.0, 4, 25);

      double intDwin = 0.0, intSwin = 0.0, simScale = 1.0;
      if (pw.ok) {
        intDwin = IntegralInXRange(hD, pw.xLo, pw.xHi);
        intSwin = IntegralInXRange(hS, pw.xLo, pw.xHi);
        if (intSwin > 0 && intDwin > 0) {
          simScale = intDwin / intSwin;
          hS->Scale(simScale);
        }
      }

      const double fitLoD = pw.ok ? pw.xLo : V.xmin;
      const double fitHiD = pw.ok ? pw.xHi : V.xmax;

      const PeakWindow pwS = EstimatePeakWindow(hS, 10, 2.0, 4, 25);
      const double fitLoS = pwS.ok ? pwS.xLo : fitLoD;
      const double fitHiS = pwS.ok ? pwS.xHi : fitHiD;

      const FitOut fD = FitGaus(hD,
                                Form("fD_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data()),
                                fitLoD, fitHiD, kBlack);
      const FitOut fS = FitGaus(hS,
                                Form("fS_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data()),
                                fitLoS, fitHiS, kRed + 1);

      double muDMeV = -999.0;
      double muDErrMeV = -999.0;
      double sigDMeV = -999.0;
      double sigDErrMeV = -999.0;
      double muSMeV = -999.0;
      double muSErrMeV = -999.0;
      double sigSMeV = -999.0;
      double sigSErrMeV = -999.0;
      double offsetMeV = -999.0;
      double offsetErrMeV = -999.0;
      if (fD.ok && fS.ok) {
        muDMeV = fD.mu * 1000.0;
        muDErrMeV = fD.muErr * 1000.0;
        sigDMeV = fD.sig * 1000.0;
        sigDErrMeV = fD.sigErr * 1000.0;
        muSMeV = fS.mu * 1000.0;
        muSErrMeV = fS.muErr * 1000.0;
        sigSMeV = fS.sig * 1000.0;
        sigSErrMeV = fS.sigErr * 1000.0;
        offsetMeV = muDMeV - muSMeV;
        // Propagate the uncertainty on the offset from the fitted mean errors.
        offsetErrMeV = std::sqrt(muDErrMeV * muDErrMeV + muSErrMeV * muSErrMeV);
      }

      StyleHists(hD, hS);
      hD->GetXaxis()->SetTitle(V.x_title);
      const double maxY = std::max(hD->GetMaximum(), hS->GetMaximum());
      hD->SetMaximum(1.25 * maxY);

      hD->Draw("E1");
      hS->Draw("E1 SAME");
      if (fD.f) fD.f->Draw("SAME");
      if (fS.f) fS.f->Draw("SAME");

      TLegend* L = new TLegend(0.58, 0.72, 0.90, 0.88);
      L->SetBorderSize(0);
      L->SetFillStyle(0);
      L->SetTextSize(0.042);
      L->AddEntry(hD, "Data", "lp");
      L->AddEntry(hS, "Sim (weighted, scaled)", "lp");
      L->Draw();

      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.058);
      lat.SetTextColor(kBlack);
      lat.DrawLatex(0.43, 0.92, Form("Setting %s, dp[%g,%g] (%s)",
                                     setting.label.Data(), dp.lo, dp.hi, dp.label.Data()));

      lat.SetTextSize(0.050);
      lat.DrawLatex(0.15, 0.84, Form("Data: #mu=%.5f #pm %.5f", fD.mu, fD.muErr));
      lat.DrawLatex(0.15, 0.79, Form("Sim:  #mu=%.5f #pm %.5f", fS.mu, fS.muErr));

      lat.SetTextColor(kRed + 1);
      lat.SetTextSize(0.062);
      if (fD.ok && fS.ok) {
        lat.DrawLatex(0.15, 0.72, Form("Offset = %.2f #pm %.2f MeV", offsetMeV, offsetErrMeV));
      } else {
        lat.DrawLatex(0.15, 0.72, Form("Offset invalid  (statusD=%d, statusS=%d)", fD.status, fS.status));
      }

      lat.SetTextColor(kGray + 2);
      lat.SetTextSize(0.042);
      if (pw.ok) lat.DrawLatex(0.15, 0.66, Form("NormWin: [%.5f, %.5f], SimScale=%.4g", pw.xLo, pw.xHi, simScale));

      AppendCsvRow(csvPath, setting.label, dp, V, pw,
                   intDwin, intSwin, simScale,
                   hD->GetEntries(), hS->GetEntries(),
                   fD, fS,
                   muDMeV, muDErrMeV, sigDMeV, sigDErrMeV,
                   muSMeV, muSErrMeV, sigSMeV, sigSErrMeV,
                   offsetMeV, offsetErrMeV,
                   (fD.ok && fS.ok));
    }

    const TString outName = outDir + Form("/offsets_Setting%s_%sbin.png",
                                          setting.label.Data(), dp.label.Data());
    c->SaveAs(outName);
    std::cout << "[INFO] Saved: " << outName << "\n";
    delete c;
  }
}

void MeasuredOffsetsBySetting(double dpFocusLo = -2.0,
                              double dpFocusHi =  3.0,
                              int nFocusBins = 5,
                              const char* outDirC = "results/PNGs/MeasuredOffsetsBySetting",
                              const char* csvPathC = "results/tables/MeasuredOffsetsBySetting.csv") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const TString outDir(outDirC);
  const TString csvPath(csvPathC);

  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);
  EnsureCsvHeader(csvPath);

  const std::vector<SettingInfo> settings = GetSettings();
  for (const auto& setting : settings) {
    ProcessOneSetting(setting, dpFocusLo, dpFocusHi, nFocusBins, outDir, csvPath);
  }
}
