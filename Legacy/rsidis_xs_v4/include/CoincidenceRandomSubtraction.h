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
  TString CtBranchName = "CTime_ePiCoinTime_ROC2";

  double WideWindowMinNs = 20.0;
  double WideWindowMaxNs = 80.0;

  int    CtHistogramNBins = 400;

  double RfPeriodNs = 4.0;

  double PeakHalfWidthNs = 1.0;

  int    MaxSidePeaks = 3;
};

struct CoincidenceResult {
  double PeakCenterNs = 0.0;

  std::pair<double,double> CoinWindowNs = {0.0, 0.0};

  double CoinYield = 0.0;
  double CoinYieldErr = 0.0;
  double RandomMeanYield = 0.0;
  double RandomMeanYieldErr = 0.0;
  double RandomSubtractedYield = 0.0;
  double RandomSubtractedYieldErr = 0.0;

  std::vector<std::pair<double,double>> RandomWindowListNs;
};

inline TString CombineCutsAND(const TString& A, const TString& B) {
  if (A.IsNull() || A.Length()==0) return B;
  if (B.IsNull() || B.Length()==0) return A;
  return Form("(%s)&&(%s)", A.Data(), B.Data());
}

inline TString BuildRangeCut(const char* Var, double Lo, double Hi) {
  return Form("(%s > %.2f && %s < %.2f)", Var, Lo, Var, Hi);
}

inline void EnsureSumw2(TH1* h) {
  if (h && h->GetSumw2N() == 0) h->Sumw2();
}



// Helper: unique histogram names to avoid ROOT directory name collisions.
inline TString UniqueHistName(const TString& Prefix) {
  static unsigned long long counter = 0;
  ++counter;
  return Form("%s_%llu", Prefix.Data(), counter);
}

// Build a selection expression suitable for TTree::Project where WeightExpr is multiplied
// by a boolean CutExpr. If WeightExpr is empty, return CutExpr.
inline TString BuildWeightedSelection(const TString& WeightExpr, const TString& CutExpr) {
  if (WeightExpr.IsNull() || WeightExpr.Length()==0) return CutExpr;
  return Form("(%s) * (%s)", WeightExpr.Data(), CutExpr.Data());
}

inline CoincidenceResult ComputeCoincidenceRandomSubtraction(
    TTree* Tree,
    const TString& BaseCuts,
    const CoincidenceConfig& Config)
{
  CoincidenceResult R;

  TString WideGate = BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs);
  TString CutsWide = CombineCutsAND(BaseCuts, WideGate);

  std::unique_ptr<TH1D> Hct(new TH1D(UniqueHistName("Hct").Data(),";Coincidence time (ns);Counts",
                                     Config.CtHistogramNBins,
                                     Config.WideWindowMinNs, Config.WideWindowMaxNs));
  EnsureSumw2(Hct.get());
  Tree->Project(Hct->GetName(), Config.CtBranchName, CutsWide);
  if (Hct->GetEntries()==0) return R;

  int    MaxBin     = Hct->GetMaximumBin();
  double PeakCenter = Hct->GetBinCenter(MaxBin);
  R.PeakCenterNs    = PeakCenter;

  double CoinLo = PeakCenter - Config.PeakHalfWidthNs;
  double CoinHi = PeakCenter + Config.PeakHalfWidthNs;
  R.CoinWindowNs   = {CoinLo, CoinHi};

  int CoinBinLo = Hct->GetXaxis()->FindFixBin(CoinLo);
  int CoinBinHi = Hct->GetXaxis()->FindFixBin(CoinHi) - 1;
  double CoinErr = 0.0;
  double CoinVal = Hct->IntegralAndError(CoinBinLo, CoinBinHi, CoinErr, "width");
  R.CoinYield     = CoinVal;
  R.CoinYieldErr  = CoinErr;

  double SumRand = 0.0, SumRandVar = 0.0;
  int    Nused   = 0;
  for (int k=1; k<=Config.MaxSidePeaks; ++k) {
    if (k == 1) continue;
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
    R.RandomMeanYield     = SumRand / Nused;
    R.RandomMeanYieldErr  = std::sqrt(SumRandVar) / Nused;
    R.RandomSubtractedYield    = R.CoinYield - R.RandomMeanYield;
    R.RandomSubtractedYieldErr = std::hypot(R.CoinYieldErr, R.RandomMeanYieldErr);
  }
  return R;
}


// Fill histogram using a precomputed CoincidenceResult (windows already decided).
// This avoids re-finding the peak for every variable.
// BaseCuts must be the same boolean cut string that was used to compute R.
inline void FillRandomSubtractedHistogramWithWindows(
    TTree* Tree,
    const TString& BaseCuts,
    const TString& WeightExpr,
    const char* VarExpression,
    TH1* OutputHist,
    const CoincidenceConfig& Config,
    const CoincidenceResult& R)
{
  OutputHist->Reset();
  EnsureSumw2(OutputHist);

  if (R.CoinWindowNs.first >= R.CoinWindowNs.second) return;

  std::unique_ptr<TH1D> Hcoin(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("Hcoin").Data())));
  Hcoin->Reset();
  EnsureSumw2(Hcoin.get());

  TString CoinCutBool = CombineCutsAND(BaseCuts,
                       CombineCutsAND(
                         BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                         BuildRangeCut(Config.CtBranchName, R.CoinWindowNs.first, R.CoinWindowNs.second)));

  TString CoinSelect = BuildWeightedSelection(WeightExpr, CoinCutBool);
  Tree->Project(Hcoin->GetName(), VarExpression, CoinSelect);

  std::unique_ptr<TH1D> HrandSum(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("HrandSum").Data())));
  HrandSum->Reset();
  EnsureSumw2(HrandSum.get());

  int M = 0;
  for (const auto& win : R.RandomWindowListNs) {
    TString RandCutBool = CombineCutsAND(BaseCuts,
                         CombineCutsAND(
                           BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                           BuildRangeCut(Config.CtBranchName, win.first, win.second)));

    std::unique_ptr<TH1D> Htmp(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("HtmpRand").Data())));
    Htmp->Reset();
    EnsureSumw2(Htmp.get());

    TString RandSelect = BuildWeightedSelection(WeightExpr, RandCutBool);
    Tree->Project(Htmp->GetName(), VarExpression, RandSelect);
    HrandSum->Add(Htmp.get());
    ++M;
  }

  if (M > 0) HrandSum->Scale(1.0 / M);

  OutputHist->Add(Hcoin.get());
  OutputHist->Add(HrandSum.get(), -1.0);
}

// Weighted version.
// BaseCuts must be boolean cuts only.
// WeightExpr is multiplied into the selection used by TTree::Project.
inline CoincidenceResult FillRandomSubtractedHistogram(
    TTree* Tree,
    const TString& BaseCuts,
    const TString& WeightExpr,
    const char* VarExpression,
    TH1* OutputHist,
    const CoincidenceConfig& Config)
{
  OutputHist->Reset();
  EnsureSumw2(OutputHist);

  // Window finding is done with BaseCuts only
  CoincidenceResult R = ComputeCoincidenceRandomSubtraction(Tree, BaseCuts, Config);

  if (R.CoinWindowNs.first >= R.CoinWindowNs.second) return R;

  std::unique_ptr<TH1D> Hcoin(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("Hcoin").Data())));
  Hcoin->Reset();
  EnsureSumw2(Hcoin.get());

  TString CoinCutBool = CombineCutsAND(BaseCuts,
                       CombineCutsAND(
                         BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                         BuildRangeCut(Config.CtBranchName, R.CoinWindowNs.first, R.CoinWindowNs.second)));

  TString CoinSelect = BuildWeightedSelection(WeightExpr, CoinCutBool);
  Tree->Project(Hcoin->GetName(), VarExpression, CoinSelect);

  std::unique_ptr<TH1D> HrandSum(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("HrandSum").Data())));
  HrandSum->Reset();
  EnsureSumw2(HrandSum.get());

  int M = 0;
  for (const auto& win : R.RandomWindowListNs) {
    TString RandCutBool = CombineCutsAND(BaseCuts,
                         CombineCutsAND(
                           BuildRangeCut(Config.CtBranchName, Config.WideWindowMinNs, Config.WideWindowMaxNs),
                           BuildRangeCut(Config.CtBranchName, win.first, win.second)));

    std::unique_ptr<TH1D> Htmp(static_cast<TH1D*>(OutputHist->Clone(UniqueHistName("HtmpRand").Data())));
    Htmp->Reset();
    EnsureSumw2(Htmp.get());

    TString RandSelect = BuildWeightedSelection(WeightExpr, RandCutBool);
    Tree->Project(Htmp->GetName(), VarExpression, RandSelect);
    HrandSum->Add(Htmp.get());
    ++M;
  }

  if (M > 0) HrandSum->Scale(1.0 / M);

  OutputHist->Add(Hcoin.get());
  OutputHist->Add(HrandSum.get(), -1.0);

  return R;
}

// Backward compatible version (no weight)
inline CoincidenceResult FillRandomSubtractedHistogram(
    TTree* Tree,
    const TString& BaseCuts,
    const char* VarExpression,
    TH1* OutputHist,
    const CoincidenceConfig& Config)
{
  return FillRandomSubtractedHistogram(Tree, BaseCuts, "", VarExpression, OutputHist, Config);
}
