// MeasuredOffsetsBySetting.C
//
// Purpose:
//   For each hardcoded setting, chain the corresponding Data and Sim ROOT files,
//   build SHMS dp-bin histograms, fit Gaussian means for
//     - W
//     - Em
//     - Pmz
//     - Pmx
//   and measure the setting-wise offset
//       offset = mu_data - mu_sim
//
// Settings:
//   4-pass:
//     E_beam = 8.5831 GeV
//     Setting A: runs 23839-23848 except 23843
//       HMS P = 2.28 GeV, HMS angle = 42.0 deg
//       SHMS P = -7.07 GeV, SHMS angle = 12.465 deg
//     Setting B: runs 23849-23851
//       HMS P = 3.75 GeV, HMS angle = 29.99 deg
//       SHMS P = -5.67 GeV, SHMS angle = 19.34 deg
//     Cuts:
//       H.gtr.dp / hsdelta in [-1, 1]
//       one SHMS dp bin:
//         b1 [-2,1]
//   5-pass:
//     E_beam = 10.6716 GeV
//     Current study uses only settings C, D, and E:
//     Setting C: run 25476, SHMS P = -6.7200 GeV
//     Setting D: run 25477, SHMS P = -6.3840 GeV
//     Setting E: run 25478, SHMS P = -6.6048 GeV
//     For C/D/E, HMS P = 4.72 GeV, HMS angle = 26.395 deg,
//     and SHMS angle = 18.56 deg.
//     Cuts:
//       H.gtr.dp / hsdelta in [-1, 1]
//       C b1: [-0.5, 2.5]
//       D b1: [ 4.5, 7.5]
//       E b1: [ 1.5, 4.5]
//
// Per setting, the macro builds dp-bin canvases:
//   - For 4-pass: one narrow bin [-2,1]
//   - For 5-pass: exactly one narrow bin per setting using the hardcoded
//     setting-specific SHMS dp range. This avoids low-statistics slicing
//     when using the 5-pass dataset.
//
// Each canvas is a 2x2 layout with Data vs Sim overlays for:
//   upper left:  W
//   upper right: Em
//   lower left:  Pmz
//   lower right: Pmx
//
// Analysis notes:
//   - Data and Sim ROOT directories are passed to MeasuredOffsetsBySetting().
//     By default these are local symlink directories:
//       Pass0p1_DataROOTfiles/
//       Pass0p1_SimROOTfiles/
//   - Sim is filled from the h10 tree using event-by-event Weight
//   - Observable histogram ranges are intentionally narrow:
//       W   [0.85, 1.05]
//       Em  [-0.025, 0.025]
//       Pmz [-0.025, 0.025]
//       Pmx [-0.030, 0.030]
//   - The number of histogram bins is selected from the Data entries inside
//     the observable range:
//       nBins = int(2*sqrt(nEntriesDataInRange)), constrained to [10,200]
//     The same nBins is then used for Data and Sim for that observable/bin.
//   - Data and Sim each use a two-stage Gaussian fit:
//       first fit: full observable histogram range
//       final fit: mu1 +/- 1.5*sigma1 from the first fit
//     There is no physics/window clamping in this diagnostic version.
//   - Sim is normalized to Data inside the Data final fit window
//   - If the first approximate Gaussian fit fails, the row is marked invalid;
//     no weighted-moment fallback is used
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
//   results/PNGs/MeasuredOffsetsBySetting_<pass-or-study>/offsets_Setting<label>_<binlabel>bin.png
//
// Output CSV:
//   results/tables/MeasuredOffsetsBySetting_<pass-or-study>.csv
//
// Run from:
//   RSIDIS/heep_check_v3/
//
// Examples:
//   4-pass:
//   root -l -b -q 'macros/MeasuredOffsetsBySetting.C(1, "", "", "Pass0p1_DataROOTfiles", "Pass0p1_SimROOTfiles", "4pass" )'
//
//   5-pass:
//   root -l -b -q 'macros/MeasuredOffsetsBySetting.C(1, "", "", "Pass0p1_DataROOTfiles", "Pass0p1_SimROOTfiles", "5pass" )'

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
static TString gDataDir = "Pass0p1_DataROOTfiles/";
static TString gSimDir  = "Pass0p1_SimROOTfiles/";

static TString NormalizeDir(TString dir) {
  if (!dir.EndsWith("/")) dir += "/";
  return dir;
}

static TString GetDataFileName(int run) {
  return gDataDir + Form("coin_replay_production_%d_-1.root", run);
}
static TString GetSimFileName(int run) {
  return gSimDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
}

// Trees
static const TString kDataTreeName = "T";
static const TString kSimTreeName  = "h10";

// Baseline HMS cuts. These are set by MeasuredOffsetsBySetting() from the
// requested setting set; 5-pass uses the tighter HMS dp cut from the partner
// comparison study.
static TString gDataCutBase = "(H.gtr.dp>-1) && (H.gtr.dp<1)";
static TString gSimCutBase  = "(hsdelta>-1) && (hsdelta<1)";
static bool gSkipShmsDpCut = false;

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
               100, 0.85, 1.05});
  v.push_back({"Em",
               "H.kin.secondary.emiss",
               "Em_recon",
               "Em (GeV)",
               100, -0.025, 0.025});
  v.push_back({"Pmz",
               "H.kin.secondary.pmiss_z",
               "Pmz_recon",
               "Pmz (GeV/c)",
               100, -0.025, 0.025});
  v.push_back({"Pmx",
               "H.kin.secondary.pmiss_x",
               "Pmx_recon",
               "Pmx (GeV/c)",
               100, -0.030, 0.030});
  return v;
}

struct DpBin {
  int idx;
  TString label;
  double lo;
  double hi;
};

struct SettingInfo {
  TString label;
  std::vector<int> runs;
  double eBeam = NAN;
  double hmsP = NAN;
  double hmsAngleDeg = NAN;
  double shmsP = NAN;
  double shmsAngleDeg = NAN;
  double shmsDpLo = -2.0;  // Default SHMS dp low cut
  double shmsDpHi = 3.0;   // Default SHMS dp high cut
};

static std::vector<DpBin> BuildDpBins(const SettingInfo& setting) {
  std::vector<DpBin> bins;

  if (setting.label == "A" || setting.label == "B") {
    if (gSkipShmsDpCut) {
      bins.push_back({1, "b1", -999.0, 999.0});
      return bins;
    }
    bins.push_back({1, "b1", -2.0, 1.0});
    return bins;
  }

  bins.push_back({1, "b1", setting.shmsDpLo, setting.shmsDpHi});
  return bins;
}

static std::vector<SettingInfo> Get4PassSettings() {
  const double eBeam = 8.5831;
  return {
    {"A", {23839, 23840, 23841, 23842, 23844, 23845, 23846, 23847, 23848},
     eBeam, 2.28, 42.0, -7.07, 12.465},
    {"B", {23849, 23850, 23851},
     eBeam, 3.75, 29.99, -5.67, 19.34}
  };
}

static std::vector<SettingInfo> Get5PassSettings() {
  const double eBeam = 10.6716;
  const double hmsP = 4.72;
  const double hmsAngleDeg = 26.395;
  const double shmsAngleDeg = 18.56;
  return {
    // 5-pass comparison study: C/D/E only, with partner delta windows.
    {"C", {25476},        eBeam, hmsP, hmsAngleDeg, -6.7200, shmsAngleDeg, -0.5, 2.5},
    {"D", {25477},        eBeam, hmsP, hmsAngleDeg, -6.3840, shmsAngleDeg,  4.5, 7.5},
    {"E", {25478},        eBeam, hmsP, hmsAngleDeg, -6.6048, shmsAngleDeg,  1.5, 4.5}
  };
}

static std::vector<SettingInfo> GetSettings(const char* passName) {
  TString pass(passName ? passName : "");
  pass.ToLower();

  if (pass == "4pass" || pass == "4" || pass == "pass4" ||
      pass == "4passnoshmsdpcut" || pass == "4pass_noshmsdpcut" ||
      pass == "4pass_noshms") return Get4PassSettings();
  if (pass == "5pass" || pass == "5" || pass == "pass5") return Get5PassSettings();

  std::cerr << "[ERROR] Unknown setting set \"" << passName
            << "\". Use \"4pass\", \"4passNoShmsDpCut\", or \"5pass\".\n";
  return {};
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

struct FitWindow {
  double xLo = 0;
  double xHi = 0;
  double muSeed = 0;
  double sigSeed = 0;
  double muSeedErr = 0;
  double sigSeedErr = 0;
  double chi2Seed = 0;
  int ndfSeed = 0;
  int statusSeed = -999;
  bool ok = false;
};

static const int kMinEntriesForFit = 20;

static bool IsFinitePositive(double x) {
  return std::isfinite(x) && x > 0.0;
}

static FitWindow EstimateFitWindowFromFullRange(TH1D* h,
                                                const TString& fname,
                                                double finalKSigma = 1.5,
                                                Color_t lineColor = kGray + 2) {
  FitWindow pw;
  if (!h) return pw;

  if (h->GetEntries() < kMinEntriesForFit) return pw;

  const double xMin = h->GetXaxis()->GetXmin();
  const double xMax = h->GetXaxis()->GetXmax();
  if (!(xMax > xMin)) return pw;

  const int bMax = h->GetMaximumBin();
  const double binW = h->GetXaxis()->GetBinWidth(1);
  TF1 fSeed(fname, "gaus", xMin, xMax);
  fSeed.SetLineColor(lineColor);
  fSeed.SetLineStyle(2);
  fSeed.SetLineWidth(1);
  fSeed.SetParameters(h->GetBinContent(bMax), h->GetXaxis()->GetBinCenter(bMax),
                      std::max(2.0 * binW, 0.2 * (xMax - xMin)));

  TFitResultPtr r = h->Fit(&fSeed, "QRSN");
  pw.statusSeed = r.Get() ? r->Status() : -999;
  if (r.Get()) {
    pw.chi2Seed = r->Chi2();
    pw.ndfSeed = r->Ndf();
  }
  if (pw.statusSeed != 0) return pw;

  const double mu = fSeed.GetParameter(1);
  const double sig = std::fabs(fSeed.GetParameter(2));
  const double muErr = fSeed.GetParError(1);
  const double sigErr = fSeed.GetParError(2);
  if (!std::isfinite(mu) || !IsFinitePositive(sig)) return pw;

  const double xLo = mu - finalKSigma * sig;
  const double xHi = mu + finalKSigma * sig;
  if (!std::isfinite(xLo) || !std::isfinite(xHi) || !(xHi > xLo)) return pw;
  if (xHi <= xMin || xLo >= xMax) return pw;

  pw.xLo = xLo;
  pw.xHi = xHi;
  pw.muSeed = mu;
  pw.sigSeed = sig;
  pw.muSeedErr = muErr;
  pw.sigSeedErr = sigErr;
  pw.ok = true;
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
  if (!h || !(xHi > xLo) || h->GetEntries() < kMinEntriesForFit) return out;

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
    << "setting,e_beam_GeV,hms_p_GeV,hms_angle_deg,shms_p_GeV,shms_angle_deg,"
    << "dp_idx,dp_label,dp_lo,dp_hi,bin_center,"
    << "var,"
    << "nBins,"
    << "norm_win_lo,norm_win_hi,"
    << "fitD_lo,fitD_hi,fitS_lo,fitS_hi,"
    << "intD_win,intS_win,sim_scale,"
    << "entriesD,entriesS,"
    << "statusD1,muD1_MeV,muD1_err_MeV,sigD1_MeV,sigD1_err_MeV,chi2D1,ndfD1,"
    << "statusS1,muS1_MeV,muS1_err_MeV,sigS1_MeV,sigS1_err_MeV,chi2S1,ndfS1,"
    << "muD_MeV,muD_err_MeV,sigD_MeV,sigD_err_MeV,statusD,chi2D,ndfD,"
    << "muS_MeV,muS_err_MeV,sigS_MeV,sigS_err_MeV,statusS,chi2S,ndfS,"
    << "offset_MeV,offset_err_MeV,fit_valid"
    << "\n";
}

static void AppendCsvRow(const TString& csvPath,
                         const SettingInfo& setting, const DpBin& dp, const VarInfo& V,
                         int nBins,
                         const FitWindow& pwD, const FitWindow& pwS,
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
    << setting.label.Data() << ","
    << setting.eBeam << ","
    << setting.hmsP << ","
    << setting.hmsAngleDeg << ","
    << setting.shmsP << ","
    << setting.shmsAngleDeg << ","
    << dp.idx << ","
    << dp.label.Data() << ","
    << dp.lo << ","
    << dp.hi << ","
    << 0.5 * (dp.lo + dp.hi) << ","
    << V.name.Data() << ","
    << nBins << ","
    << pwD.xLo << ","
    << pwD.xHi << ","
    << pwD.xLo << ","
    << pwD.xHi << ","
    << pwS.xLo << ","
    << pwS.xHi << ","
    << intDwin << ","
    << intSwin << ","
    << simScale << ","
    << entriesD << ","
    << entriesS << ","
    << pwD.statusSeed << ","
    << pwD.muSeed * 1000.0 << ","
    << pwD.muSeedErr * 1000.0 << ","
    << pwD.sigSeed * 1000.0 << ","
    << pwD.sigSeedErr * 1000.0 << ","
    << pwD.chi2Seed << ","
    << pwD.ndfSeed << ","
    << pwS.statusSeed << ","
    << pwS.muSeed * 1000.0 << ","
    << pwS.muSeedErr * 1000.0 << ","
    << pwS.sigSeed * 1000.0 << ","
    << pwS.sigSeedErr * 1000.0 << ","
    << pwS.chi2Seed << ","
    << pwS.ndfSeed << ","
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
                              const TString& outDir,
                              const TString& csvPath) {
  TChain tData(kDataTreeName);
  TChain tSim(kSimTreeName);
  if (!BuildChainsForSetting(setting, tData, tSim)) {
    std::cerr << "[ERROR] No valid data/sim file pairs found for Setting " << setting.label << "\n";
    return;
  }

  std::vector<VarInfo> vars = GetVariables();
  std::vector<DpBin> dpBins = BuildDpBins(setting);

  for (const auto& dp : dpBins) {
    TCanvas* c = new TCanvas(Form("c_heep_setting%s_%s", setting.label.Data(), dp.label.Data()),
                             Form("HEEP offsets Setting %s dp %s", setting.label.Data(), dp.label.Data()),
                             1600, 1000);
    c->Divide(2, 2);

    const TString dataDpCut = gSkipShmsDpCut ? "1" : Form("(P.gtr.dp>%g) && (P.gtr.dp<%g)", dp.lo, dp.hi);
    const TString simDpCut  = gSkipShmsDpCut ? "1" : Form("(ssdelta>%g) && (ssdelta<%g)", dp.lo, dp.hi);
    const TString dataCut = Form("(%s) && (%s)", gDataCutBase.Data(), dataDpCut.Data());
    const TString simCut  = Form("(%s) && (%s)", gSimCutBase.Data(),  simDpCut.Data());

    for (int i = 0; i < (int)vars.size(); ++i) {
      const VarInfo& V = vars[i];
      c->cd(i + 1);
      gPad->SetLeftMargin(0.12);
      // Keep enough space for ROOT's x-axis exponent, e.g. #times10^{-3}.
      gPad->SetRightMargin(0.12);
      gPad->SetBottomMargin(0.12);
      gPad->SetTopMargin(0.08);

      const TString hDname = Form("hD_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data());
      const TString hSname = Form("hS_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data());

      const TString dataObsCut = Form("(%s>%.12g) && (%s<%.12g)",
                                      V.data_expr.Data(), V.xmin,
                                      V.data_expr.Data(), V.xmax);
      const TString dataCountCut = Form("(%s) && (%s)", dataCut.Data(), dataObsCut.Data());
      const Long64_t nEventsDInRange = tData.GetEntries(dataCountCut);
      int nBins = (int)(2.0 * std::sqrt((double)std::max((Long64_t)0, nEventsDInRange)));
      nBins = std::max(10, std::min(nBins, 200));

      TH1D* hD = new TH1D(hDname, "", nBins, V.xmin, V.xmax);
      TH1D* hS = new TH1D(hSname, "", nBins, V.xmin, V.xmax);
      hD->Sumw2();
      hS->Sumw2();

      tData.Draw(Form("%s>>%s", V.data_expr.Data(), hDname.Data()), dataCut.Data(), "goff");
      const TString simSel = Form("%s*(%s)", kSimWeightBranch, simCut.Data());
      tSim.Draw(Form("%s>>%s", V.sim_expr.Data(), hSname.Data()), simSel.Data(), "goff");

      const FitWindow pwD = EstimateFitWindowFromFullRange(
          hD,
          Form("fD1_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data()),
          1.5, kGray + 2);

      const FitWindow pwS1 = EstimateFitWindowFromFullRange(
          hS,
          Form("fS1_%s_setting%s_%s", V.name.Data(), setting.label.Data(), dp.label.Data()),
          1.5, kMagenta + 2);

      double intDwin = 0.0, intSwin = 0.0, simScale = 1.0;
      if (pwD.ok) {
        intDwin = IntegralInXRange(hD, pwD.xLo, pwD.xHi);
        intSwin = IntegralInXRange(hS, pwD.xLo, pwD.xHi);
        if (intSwin > 0 && intDwin > 0) {
          simScale = intDwin / intSwin;
          hS->Scale(simScale);
        }
      }

      const double fitLoD = pwD.ok ? pwD.xLo : 0.0;
      const double fitHiD = pwD.ok ? pwD.xHi : 0.0;

      const FitWindow pwS = pwS1;
      const double fitLoS = pwS.ok ? pwS.xLo : 0.0;
      const double fitHiS = pwS.ok ? pwS.xHi : 0.0;

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
      const bool fitValid = (nEventsDInRange >= kMinEntriesForFit &&
                             hS->GetEntries() >= kMinEntriesForFit &&
                             pwD.ok && pwS.ok &&
                             fD.ok && fS.ok &&
                             std::isfinite(fD.mu) && std::isfinite(fD.muErr) &&
                             std::isfinite(fD.sig) && std::isfinite(fD.sigErr) &&
                             std::isfinite(fS.mu) && std::isfinite(fS.muErr) &&
                             std::isfinite(fS.sig) && std::isfinite(fS.sigErr));

      if (fitValid) {
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
      const TString dpText = gSkipShmsDpCut
        ? Form("Setting %s, no SHMS dp cut (%s)", setting.label.Data(), dp.label.Data())
        : Form("Setting %s, dp[%g,%g] (%s)", setting.label.Data(), dp.lo, dp.hi, dp.label.Data());
      lat.DrawLatex(0.43, 0.92, dpText);

      lat.SetTextSize(0.050);
      lat.DrawLatex(0.15, 0.84, Form("Data: #mu=%.5f #pm %.5f", fD.mu, fD.muErr));
      lat.DrawLatex(0.15, 0.79, Form("Sim:  #mu=%.5f #pm %.5f", fS.mu, fS.muErr));

      lat.SetTextColor(kRed + 1);
      lat.SetTextSize(0.062);
      if (fitValid) {
        lat.DrawLatex(0.15, 0.72, Form("Offset = %.2f #pm %.2f MeV", offsetMeV, offsetErrMeV));
      } else {
        lat.DrawLatex(0.15, 0.72, Form("Offset invalid  (D1=%d,S1=%d,D=%d,S=%d)",
                                       pwD.statusSeed, pwS.statusSeed, fD.status, fS.status));
      }

      lat.SetTextColor(kGray + 2);
      lat.SetTextSize(0.042);
	      if (pwD.ok) lat.DrawLatex(0.15, 0.66, Form("NormWin: [%.5f, %.5f], SimScale=%.4g", pwD.xLo, pwD.xHi, simScale));
	      lat.DrawLatex(0.15, 0.56, Form("nD=%lld, nBins=%d", nEventsDInRange, nBins));
	      if (std::isfinite(setting.shmsP)) {
	        lat.DrawLatex(0.15, 0.51, Form("HMS P=%.3f, SHMS P=%.4f GeV", setting.hmsP, setting.shmsP));
	      }
	
	      AppendCsvRow(csvPath, setting, dp, V, nBins, pwD, pwS,
                   intDwin, intSwin, simScale,
                   hD->GetEntries(), hS->GetEntries(),
                   fD, fS,
                   muDMeV, muDErrMeV, sigDMeV, sigDErrMeV,
                   muSMeV, muSErrMeV, sigSMeV, sigSErrMeV,
                   offsetMeV, offsetErrMeV,
                   fitValid);
    }

    const TString outName = outDir + Form("/offsets_Setting%s_%sbin.png",
                                          setting.label.Data(), dp.label.Data());
    c->SaveAs(outName);
    std::cout << "[INFO] Saved: " << outName << "\n";
    delete c;
  }
}

void MeasuredOffsetsBySetting(int nFocusBins = 1,
	                              const char* outDirC = "",
	                              const char* csvPathC = "",
	                              const char* dataDirC = "Pass0p1_DataROOTfiles",
	                              const char* simDirC = "Pass0p1_SimROOTfiles",
	                              const char* settingSetC = "5pass") {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TString settingSet(settingSetC ? settingSetC : "");
  settingSet.ToLower();

  TString passTag = "5pass";
  const bool is4PassNoShmsDpCut =
    (settingSet == "4passnoshmsdpcut" || settingSet == "4pass_noshmsdpcut" ||
     settingSet == "4pass_noshms");

  if (settingSet == "4pass" || settingSet == "4" || settingSet == "pass4" || is4PassNoShmsDpCut) {
    passTag = "4pass";
    if (is4PassNoShmsDpCut) passTag = "4passNoShmsDpCut";
  } else if (settingSet == "5pass" || settingSet == "5" || settingSet == "pass5" || settingSet.IsNull()) {
    passTag = "5pass";
  }

  TString outDir(outDirC ? outDirC : "");
  TString csvPath(csvPathC ? csvPathC : "");
  if (outDir.IsNull()) outDir = Form("results/PNGs/MeasuredOffsetsBySetting_%s", passTag.Data());
  if (csvPath.IsNull()) csvPath = Form("results/tables/MeasuredOffsetsBySetting_%s.csv", passTag.Data());

  gDataDir = NormalizeDir(TString(dataDirC));
  gSimDir  = NormalizeDir(TString(simDirC));

	  std::cout << "[INFO] Data dir: " << gDataDir << "\n";
	  std::cout << "[INFO] Sim dir:  " << gSimDir << "\n";
	  std::cout << "[INFO] Settings: " << settingSetC << "\n";

  if (gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir, kTRUE);
  EnsureCsvHeader(csvPath);

	  const std::vector<SettingInfo> settings = GetSettings(settingSetC);
	  if (settings.empty()) return;
  gSkipShmsDpCut = false;
	  if (settingSet == "5pass" || settingSet == "5" || settingSet == "pass5") {
	    gDataCutBase = "(H.gtr.dp>-1) && (H.gtr.dp<1)";
	    gSimCutBase  = "(hsdelta>-1) && (hsdelta<1)";
	    std::cout << "[INFO] 5-pass HMS dp cut: [-1, 1]\n";
    if (nFocusBins != 1) {
      std::cout << "[INFO] For 5-pass, forcing nFocusBins=1 to use one SHMS dp bin per setting.\n";
      nFocusBins = 1;
    }
	  } else {
	    gDataCutBase = "(H.gtr.dp>-1) && (H.gtr.dp<1)";
	    gSimCutBase  = "(hsdelta>-1) && (hsdelta<1)";
	    gSkipShmsDpCut = is4PassNoShmsDpCut;
	    std::cout << "[INFO] 4-pass HMS dp cut: [-1, 1]\n";
	    if (gSkipShmsDpCut) {
	      std::cout << "[INFO] 4-pass diagnostic: no SHMS dp cut applied\n";
	    }
    if (nFocusBins != 1) {
      std::cout << "[INFO] For 4-pass, forcing nFocusBins=1.\n";
      nFocusBins = 1;
    }
	  }
	  for (const auto& setting : settings) {
    ProcessOneSetting(setting, outDir, csvPath);
  }
}
