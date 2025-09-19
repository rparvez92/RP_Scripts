// CoincidenceRandomSubtraction.h
#pragma once
#include <memory>
#include <utility>
#include <vector>
#include <cmath>
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TAxis.h"

struct CoincidenceConfig {
  // Branch name for coincidence-time
  TString CtBranchName = "CTime.ePiCoinTime_ROC1";

  // Step 1: wide time gate you want to enforce, e.g. (30,70) ns
  double WideWindowMinNs = 20.0;
  double WideWindowMaxNs = 80.0;

  // Internal histogram for finding the peak
  int    CtHistogramNBins = 400;

  // Beam RF period (ns). 250 MHz -> 4 ns spacing of random peaks.
  double RfPeriodNs = 4.0;

  // Half-width of the coin window around the peak center (ns). (lo,hi)=(t0±PeakHalfWidthNs)
  double PeakHalfWidthNs = 1.0;

  // Try up to this many peaks on each side for random windows (only those inside wide window are used).
  int    MaxSidePeaks = 6;
};

struct CoincidenceResult {
  // Peak center we detected (ns)
  double PeakCenterNs = 0.0;

  // Coin (main) window edges (ns)
  std::pair<double,double> CoinWindowNs = {0.0, 0.0};

  // Yields and errors
  double CoinYield = 0.0;
  double CoinYieldErr = 0.0;
  double RandomMeanYield = 0.0;
  double RandomMeanYieldErr = 0.0;
  double RandomSubtractedYield = 0.0;
  double RandomSubtractedYieldErr = 0.0;

  // Random windows actually used (ns)
  std::vector<std::pair<double,double>> RandomWindowListNs;
};

// Helper: logical-AND two cut strings.
inline TString CombineCutsAND(const TString& A, const TString& B) {
  if (A.IsNull() || A.Length()==0) return B;
  if (B.IsNull() || B.Length()==0) return A;
  return Form("(%s)&&(%s)", A.Data(), B.Data());
}

// Helper: build a [lo,hi] cut on a variable
inline TString BuildRangeCut(const char* Var, double Lo, double Hi) {
  return Form("(%s > %.2f && %s < %.2f)", Var, Lo, Var, Hi);
}

// Helper: guard function to ensure if histogram calls Sumw2() already or not
inline void EnsureSumw2(TH1* h) {
  if (h && h->GetSumw2N() == 0) h->Sumw2();
}

// Compute CT peak, coin window, random windows and all three yields.
inline CoincidenceResult ComputeCoincidenceRandomSubtraction(
    TTree* Tree,
    const TString& BaseCuts,          // your existing d&d cuts (unchanged)
    const CoincidenceConfig& Config)
{
  CoincidenceResult R;

  // 1) Wide gate + base cuts
  TString WideGate = BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs);
  TString CutsWide = CombineCutsAND(BaseCuts, WideGate);

  // 2) Fill CT hist to find peak
  std::unique_ptr<TH1D> Hct(new TH1D("Hct",";Coincidence time (ns);Counts",
                                     Config.CtHistogramNBins,
                                     Config.WideWindowMinNs, Config.WideWindowMaxNs));
  //Hct->Sumw2();
  EnsureSumw2(Hct.get());
  Tree->Project("Hct", Config.CtBranchName, CutsWide);
  if (Hct->GetEntries()==0) return R;

  int    MaxBin     = Hct->GetMaximumBin();
  double PeakCenter = Hct->GetBinCenter(MaxBin);
  R.PeakCenterNs    = PeakCenter;

  // 3) Coin window integral
  double CoinLo = PeakCenter - Config.PeakHalfWidthNs;
  double CoinHi = PeakCenter + Config.PeakHalfWidthNs;
  R.CoinWindowNs   = {CoinLo, CoinHi};

  int CoinBinLo = Hct->GetXaxis()->FindFixBin(CoinLo);
  int CoinBinHi = Hct->GetXaxis()->FindFixBin(CoinHi) - 1; // make [CoinLo, CoinHi) in bin indices
  double CoinErr = 0.0;
  double CoinVal = Hct->IntegralAndError(CoinBinLo, CoinBinHi, CoinErr, "width");
  R.CoinYield     = CoinVal;
  R.CoinYieldErr  = CoinErr;

  // 4) Random windows at ±k*RF with same width; keep only those inside the wide gate
  double SumRand = 0.0, SumRandVar = 0.0;
  int    Nused   = 0;
  for (int k=1; k<=Config.MaxSidePeaks; ++k) {
    if (k == 1) continue;  // skip the first sideband (±1) from both side
    for (int sgn : {-1, +1}) {
      double Center = PeakCenter + sgn * k * Config.RfPeriodNs;
      double Lo = Center - Config.PeakHalfWidthNs;
      double Hi = Center + Config.PeakHalfWidthNs;
      if (Lo < Config.WideWindowMinNs || Hi > Config.WideWindowMaxNs) continue;

      int BinLo = Hct->GetXaxis()->FindFixBin(Lo);
      int BinHi = Hct->GetXaxis()->FindFixBin(Hi) - 1;

      double RandErr = 0.0;
      double RandVal = Hct->IntegralAndError(BinLo, BinHi, RandErr, "width");
      SumRand    += RandVal;
      SumRandVar += RandErr*RandErr;
      ++Nused;
      R.RandomWindowListNs.emplace_back(Lo,Hi);
    }
  }

  if (Nused > 0) {
    //cout << "Nused: " << Nused << endl;
    R.RandomMeanYield     = SumRand / Nused;
    R.RandomMeanYieldErr  = std::sqrt(SumRandVar) / Nused; // error of the mean (assuming independence)
    R.RandomSubtractedYield    = R.CoinYield - R.RandomMeanYield;
    R.RandomSubtractedYieldErr = std::hypot(R.CoinYieldErr, R.RandomMeanYieldErr);
  }
  return R;
}

// Make a random-subtracted histogram of some variable (e.g., "H.gtr.dp").
// The output histogram must exist with desired binning; it will be reset and filled.
inline CoincidenceResult FillRandomSubtractedHistogram(
    TTree* Tree,
    const TString& BaseCuts,              // your existing d&d cuts
    const char* VarExpression,            // e.g. "H.gtr.dp"
    TH1* OutputHist,                      // pre-booked with your binning
    const CoincidenceConfig& Config)
{
  OutputHist->Reset();
  //OutputHist->Sumw2();
  EnsureSumw2(OutputHist);

  // First compute windows & yields (also gives us t0)
  CoincidenceResult R = ComputeCoincidenceRandomSubtraction(Tree, BaseCuts, Config);

  // If we didn't find entries, return as-is
  if (R.CoinWindowNs.first >= R.CoinWindowNs.second) return R;

  // Coin histogram
  std::unique_ptr<TH1D> Hcoin(static_cast<TH1D*>(OutputHist->Clone("Hcoin")));
  Hcoin->Reset(); EnsureSumw2(Hcoin.get()); //Hcoin->Sumw2();
  TString CoinCut = CombineCutsAND(BaseCuts,
                       CombineCutsAND(
                         BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                         BuildRangeCut(Config.CtBranchName, R.CoinWindowNs.first, R.CoinWindowNs.second)));
  Tree->Project("Hcoin", VarExpression, CoinCut);

  // Random histograms (sum then average)
  std::unique_ptr<TH1D> HrandSum(static_cast<TH1D*>(OutputHist->Clone("HrandSum")));
  HrandSum->Reset(); EnsureSumw2(HrandSum.get()); //HrandSum->Sumw2();
  int M = 0;
  for (const auto& win : R.RandomWindowListNs) {
    TString RandCut = CombineCutsAND(BaseCuts,
                         CombineCutsAND(
                           BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                           BuildRangeCut(Config.CtBranchName, win.first, win.second)));
    std::unique_ptr<TH1D> Htmp(static_cast<TH1D*>(OutputHist->Clone("HtmpRand")));
    Htmp->Reset(); EnsureSumw2(Htmp.get()); //Htmp->Sumw2();
    Tree->Project("HtmpRand", VarExpression, RandCut);
    HrandSum->Add(Htmp.get());
    ++M;
  }

  if (M > 0) HrandSum->Scale(1.0 / M);   // average of random windows

  // Random-subtracted
  OutputHist->Add(Hcoin.get());
  OutputHist->Add(HrandSum.get(), -1.0);

  return R;
}
