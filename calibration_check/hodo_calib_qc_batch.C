// hodo_calib_qc_batch.C
//
// Plots per Spec:
//   HMS  : (1) H_gtr_beta vs H_dc_x_fp; (2) beta_mean & beta_sigma vs run
//   SHMS : (1) P_gtr_beta vs P_dc_x_fp; (2) beta_mean & beta_sigma vs run
//   COIN : (1) H_gtr_beta vs H_dc_x_fp; (2) P_gtr_beta vs P_dc_x_fp; (3) 1D CTime_ePiCoinTime_ROC2;
//          (4) coin_mean & coin_sigma vs run (from ROC2 peak)
//
// Notes:
//  * All physics cuts live in BuildCuts(Spec). No CTime gates are used.
//  * Uses branch-status pruning to disable unused branches and enable only what is needed.
//  * Batch mode; outputs PNGs under ./%specPNGs/ .
//
// Usage examples:
//   root -l -b -q 'hodo_calib_qc_batch.C+("hms","./ROOTfiles","6126,6128-6130")'
//   root -l -b -q 'hodo_calib_qc_batch.C+("shms","./ROOTfiles","24017-24025")'
//   root -l -b -q 'hodo_calib_qc_batch.C+("coin","./ROOTfiles","6126,6127")'

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TCut.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TString.h>
#include <TLine.h>
#include <TBox.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TF1.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include "runs_vec_hms.hh" // HMS run numbers listed in a vector
#include "runs_vec_coin.hh" // COIN run numbers listed in a vector

//-------------------------------------------------
// MakeFileName: build file name from Spec and Run
//-------------------------------------------------
static TString MakeFileName(const char* Spec, int RunNumber) {
  TString s(Spec ? Spec : "");
  //if (s == "hms")  return TString::Format("hms_coin_replay_production_%d_-1.root",  RunNumber);
  //if (s == "shms") return TString::Format("shms_coin_replay_production_%d_-1.root", RunNumber);
  //if (s == "coin") return TString::Format("coin_replay_production_%d_-1.root",      RunNumber);

  if (s == "hms")  return TString::Format("skimmed_hms_coin_replay_production_%d_-1.root",  RunNumber);
  if (s == "shms") return TString::Format("skimmed_shms_coin_replay_production_%d_-1.root", RunNumber);
  if (s == "coin") return TString::Format("skimmed_coin_replay_production_%d_-1.root",	RunNumber);

  return "";
}

//----------------------------------------------------------------
// ParseRunsList: "24329,24332-24334" → {24329,24332,24333,24334}
//----------------------------------------------------------------
static std::vector<int> ParseRunsList(const std::string &CsvLike) {
  std::vector<int> Runs; Runs.reserve(64);
  if (CsvLike.empty()) return Runs;
  std::stringstream SS(CsvLike);
  std::string Tok;
  while (std::getline(SS, Tok, ',')) {
    if (Tok.empty()) continue;
    auto Dash = Tok.find('-');
    if (Dash != std::string::npos) {
      int A = std::stoi(Tok.substr(0, Dash));
      int B = std::stoi(Tok.substr(Dash + 1));
      if (B < A) std::swap(A, B);
      for (int R = A; R <= B; ++R) Runs.push_back(R);
    } else {
      Runs.push_back(std::stoi(Tok));
    }
  }
  return Runs;
}

//--------------------------------------------------
// BuildCuts: centralized PID/track-quality cuts;
//   Spec = "hms"  → HMS electron PID
//   Spec = "shms" → SHMS pion PID
//   Spec = "coin" → HMS && SHMS PID
//--------------------------------------------------
static TCut BuildCuts(const TString &Spec) {
  TCut HmsPid = "(H_gtr_dp>-8) && (H_gtr_dp<8) && (H_gtr_beta>0) && (H_gtr_beta<1.2) && (H_cal_etottracknorm>0.7) && (H_cer_npeSum>2.0)";

  TCut ShmsBase = "(P_gtr_dp>-10) && (P_gtr_dp<22) && (P_gtr_beta>0) && (P_gtr_beta<1.2) && (P_cal_etottracknorm<0.8)";

  //TCut ShmsNGC = "(P_gtr_p>3.5) && (P_gtr_p<9.5) && (P_ngcer_npeSum>2)";
  TCut ShmsAero = "(P_gtr_p<2.7) && (P_aero_npeSum>2)";
  TCut ShmsHGC = "(P_gtr_p>=2.7) && (P_hgcer_npeSum>1) && (P_aero_npeSum>2)";

  //TCut ShmsPidMomentumLogic = "((P_gtr_p<2.84 && P_aero_npeSum>2) || (P_gtr_p>2.7 && P_gtr_p<9.5 && P_hgcer_npeSum>1))";//HGCER and Aerogel is momentum dependent

  TCut ShmsPidMomentumLogic = (ShmsAero || ShmsHGC);
  TCut ShmsPid = ShmsBase && ShmsPidMomentumLogic;
  //TCut ShmsPid = ShmsBase && ShmsNGC;

  if (Spec == "hms")  return HmsPid;
  if (Spec == "shms") return ShmsPid;
  if (Spec == "coin") return HmsPid && ShmsPid;
  return TCut("");
}

//------------------------------------------------------------------------------
// SetBranchStatusesForSpec: turn off all branches; enable only the needed ones.
//------------------------------------------------------------------------------
static void SetBranchStatusesForSpec(TTree *T, const TString &Spec) {
  if (!T) return;
  T->SetBranchStatus("*", 0); // disable everything
  if (Spec == "hms") {
    const char* Vars[] = {
      "H_gtr_dp","H_gtr_beta","H_cal_etottracknorm","H_cer_npeSum",
      "H_dc_x_fp"
    };
    for (auto v : Vars) T->SetBranchStatus(v, 1);
  } else if (Spec == "shms") {
    const char* Vars[] = {
      "P_gtr_dp","P_gtr_beta","P_cal_etottracknorm","P_ngcer_npeSum","P_hgcer_npeSum","P_aero_npeSum","P_gtr_p","P_dc_x_fp"
    };
    for (auto v : Vars) T->SetBranchStatus(v, 1);
  } else if (Spec == "coin") {
    const char* Vars[] = {
      // HMS
      "H_gtr_dp","H_gtr_beta","H_cal_etottracknorm","H_cer_npeSum","H_dc_x_fp",
      // SHMS
      "P_gtr_dp","P_gtr_beta","P_cal_etottracknorm","P_ngcer_npeSum","P_hgcer_npeSum","P_aero_npeSum","P_gtr_p","P_dc_x_fp",
      // Coin time
      "CTime_ePiCoinTime_ROC2"
    };
    for (auto v : Vars) T->SetBranchStatus(v, 1);
  }
}

//------------------------------------------------------------------
// DrawBetaVsXfp: 2D heatmap of beta vs x_fp with dashed beta bands
//------------------------------------------------------------------
static void DrawBetaVsXfp(TTree *T, const TString &Spec, int Run) {
  const char *Expr = nullptr;
  if      (Spec == "hms")  Expr = "H_gtr_beta:H_dc_x_fp";
  else if (Spec == "shms") Expr = "P_gtr_beta:P_dc_x_fp";
  else if (Spec == "coin") Expr = "H_gtr_beta:H_dc_x_fp";
  // HMS view; here coin run is getting the HMS plot. For SHMS plot, we added necessary lines at function call site.
  else { std::cerr << "[WARN] Unknown Spec in DrawBetaVsXfp: " << Spec << ""; return; }

  TCut Cuts = BuildCuts(Spec);

  TH2D *H2 = new TH2D("H2_BetaVsXfp","#beta vs x_{fp};x_{fp} (cm);#beta",80,-45,45,120,0.2,1.2);
  H2->Sumw2();
  T->Project("H2_BetaVsXfp", Expr, Cuts);
  H2->SetDirectory(nullptr);

  TCanvas *C = new TCanvas("C_BetaVsXfp","C_BetaVsXfp",900,700);
  C->SetRightMargin(0.12);
  gStyle->SetOptStat(0);
  H2->Draw("COLZ");

  double Xmin = H2->GetXaxis()->GetXmin();
  double Xmax = H2->GetXaxis()->GetXmax();
  TLine L1(Xmin,0.95,Xmax,0.95), L2(Xmin,1.05,Xmax,1.05);
  L1.SetLineStyle(2); L2.SetLineStyle(2);
  L1.SetLineColor(kBlack); L2.SetLineColor(kBlack);
  L1.SetLineWidth(5); L2.SetLineWidth(5);
  L1.Draw("SAME"); L2.Draw("SAME");

  // Build output png name
  TString Out;
  if      (Spec=="hms")  Out = TString::Format("hmsPNGs/hms_run%d_beta_vs_xfp.png", Run);
  else if (Spec=="shms") Out = TString::Format("shmsPNGs/shms_run%d_beta_vs_xfp.png", Run);
  else if (Spec=="coin") Out = TString::Format("coinPNGs/coin_hms_run%d_beta_vs_xfp.png", Run);
  C->SaveAs(Out);
  delete C; delete H2;
}

//----------------------------------
// DrawCoinTime1D: 1D CoinTime ROC2
//----------------------------------
static void DrawCoinTime1D(TTree *T, int Run){
  TCut Cuts = BuildCuts("coin");
  const char *Ct = "CTime_ePiCoinTime_ROC2";
  TH1D *H1 = new TH1D("H1_CTime","Coincidence Time (ROC2);CTime_ePiCoinTime_ROC2 (ns);Counts",400,0,100);
  H1->Sumw2();
  T->Project("H1_CTime", Ct, Cuts);

  TCanvas *C = new TCanvas("C_CTime1D","C_CTime1D",800,600);
  C->SetLeftMargin(0.12);
  gStyle->SetOptStat(0);

  int PeakBin = H1->GetMaximumBin();
  double PeakX = H1->GetBinCenter(PeakBin);
  double FitLo = std::max(0.0, PeakX - 2.0);
  double FitHi = std::min(100.0, PeakX + 2.0);
  TF1 *G = new TF1("G","gaus", FitLo, FitHi);
  H1->Fit(G, "QNR");  // fit quietly, no auto draw
  H1->Draw("HIST");
  G->SetLineColor(kRed);
  G->SetLineWidth(3);
  G->Draw("SAME");
  H1->SetDirectory(nullptr);

  // Optional green visual band (no cut is applied)
  //TBox Band(48.0, 0.0, 53.0, H1->GetMaximum()); Band.SetFillColor(kGreen+1); Band.SetFillStyle(3001); Band.Draw("SAME");
  C->SaveAs(TString::Format("coinPNGs/coin_run%d_ctime1D_ROC2.png",Run));
  delete C; delete H1;
}

//------------------------------------------------------------------------------
// ComputeBetaMetrics: fit 1D beta to get mean/sigma (robust window around peak)
//------------------------------------------------------------------------------
static bool ComputeBetaMetrics(TTree *T, const TString &Spec, double &Mean, double &Sigma, double &NEntries){
  Mean = Sigma = NEntries = std::nan("");
  const char *Var = (Spec=="shms")?"P_gtr_beta":"H_gtr_beta";
  TCut Cuts = BuildCuts(Spec);

  TH1D *H = new TH1D("H1_Beta",";#beta;Counts",200,0.2,1.2);
  H->Sumw2();
  T->Project("H1_Beta", Var, Cuts);
  H->SetDirectory(nullptr);

  NEntries = H->GetEntries();
  if (NEntries < 50) { delete H; return false; } // not enough stats to fit

  int PeakBin = H->GetMaximumBin();
  double PeakX = H->GetBinCenter(PeakBin);
  double FitLo = std::max(0.9, PeakX - 0.03);
  double FitHi = std::min(1.1, PeakX + 0.03);

  TF1 G("G","gaus", FitLo, FitHi);
  if (H->Fit(&G, "QNR") != 0) { delete H; return false; }
  Mean  = G.GetParameter(1);
  Sigma = G.GetParameter(2);
  delete H; return true;
}

//-------------------------------------------------------------
// ComputeCoinTimeMetrics: fit 1D ROC2 CTime to get mean/sigma
//-------------------------------------------------------------
static bool ComputeCoinTimeMetrics(TTree *T, double &Mean, double &Sigma, double &NEntries){
  Mean = Sigma = NEntries = std::nan("");
  TCut Cuts = BuildCuts("coin");
  const char *Ct = "CTime_ePiCoinTime_ROC2";
  TH1D *H = new TH1D("H1_CtFit",";CTime_ePiCoinTime_ROC2 (ns);Counts",400,0,100);
  H->Sumw2();
  T->Project("H1_CtFit", Ct, Cuts);
  H->SetDirectory(nullptr);

  NEntries = H->GetEntries();
  if (NEntries < 50) { delete H; return false; }

  int PeakBin = H->GetMaximumBin();
  double PeakX = H->GetBinCenter(PeakBin);
  double FitLo = std::max(0.0, PeakX - 2.0);
  double FitHi = std::min(100.0, PeakX + 2.0);

  TF1 G("G","gaus", FitLo, FitHi);
  if (H->Fit(&G, "QNR") != 0) { delete H; return false; }
  Mean  = G.GetParameter(1);
  Sigma = G.GetParameter(2);
  delete H; return true;
}

//------------------------------------------------------------------------------
// DrawBetaTrends: #beta mean (left) and #beta sigma (right) vs run (markers-only)
// Fixed QA ranges: mean in [0.5, 1.5], sigma in [0.0, 0.1]
// Uses categorical x-axis labels (each bin = one run).
//------------------------------------------------------------------------------
static void DrawBetaTrends(const std::vector<int> &Runs,
                           const std::vector<double> &Means,
                           const std::vector<double> &Sigmas,
                           const TString &Spec) {
  if (Runs.empty()) return;
  const int N = (int)Runs.size();

  // Fixed bands for QA
  const double y1Min = 0.5;   // beta mean
  const double y1Max = 1.5;
  const double sMin  = 0.0;   // beta sigma
  const double sMax  = 0.1;

  // Categorical frame with run-number labels
  TCanvas *C = new TCanvas("C_BetaTrends","C_BetaTrends",1000,600);
  TH1F *frame = new TH1F("frame_beta_trend",
                         TString::Format("%s: #beta mean / sigma vs run;Run;#beta mean", Spec.Data()),
                         N, 0.0, (double)N);

  // Show only every 'step'-th run label:
  int step = 5;  // try 5, 10, 20 depending on N
  for (int i = 0; i < N; ++i) {
    if ((i % step) == 0) frame->GetXaxis()->SetBinLabel(i+1, Form("%d", Runs[i]));
    else                 frame->GetXaxis()->SetBinLabel(i+1, "");
  }

  frame->SetMinimum(y1Min);
  frame->SetMaximum(y1Max);
  //for (int i = 1; i <= N; ++i) frame->GetXaxis()->SetBinLabel(i, Form("%d", Runs[i-1]));
  frame->GetXaxis()->LabelsOption("v"); // v for veritical, h for horizontal
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->Draw("HIST");

  // Points at bin centers; scale sigma to left axis for drawing
  std::vector<double> X(N), Ymean(N), YsigScaled(N);
  for (int i = 0; i < N; ++i) {
    X[i] = i + 0.5;
    Ymean[i] = Means[i];
    YsigScaled[i] = y1Min + (Sigmas[i] - sMin) * (y1Max - y1Min) / (sMax - sMin);
  }

  TGraph *gMean  = new TGraph(N, X.data(), Ymean.data());
  TGraph *gSigma = new TGraph(N, X.data(), YsigScaled.data());

  // Markers
  gMean->SetMarkerStyle(20); gMean->SetMarkerSize(1.1); gMean->SetLineStyle(0); gMean->SetLineWidth(0);
  gSigma->SetMarkerStyle(3); gSigma->SetMarkerSize(1.1); gSigma->SetLineStyle(0); gSigma->SetLineWidth(0);

  gMean->Draw("P SAME");
  gSigma->Draw("P SAME");

  // Right axis for true sigma values
  TGaxis *axisR = new TGaxis((double)N, y1Min, (double)N, y1Max, sMin, sMax, 510, "+L");
  axisR->SetTitle("#beta sigma");
  axisR->SetTitleOffset(1.2);
  axisR->Draw();

  TLegend *leg = new TLegend(0.12, 0.84, 0.24, 0.92);
  leg->AddEntry(gMean,  "mean",  "p");
  leg->AddEntry(gSigma, "sigma", "p");
  leg->Draw();

  // Save into per-spec folder (e.g., hmsPNGs/hms_beta_trends.png)
  C->SaveAs(TString::Format("%sPNGs/%s_beta_trends.png", Spec.Data(), Spec.Data()));

  delete leg; delete axisR; delete frame; delete gMean; delete gSigma; delete C;
}

//------------------------------------------------------------------------------
// DrawCoinTimeTrends: CTime mean (left) and CTime sigma (right) vs run (ROC2)
// Fixed QA ranges: mean in [45, 55] ns, sigma in [0, 1] ns
// Markers-only; categorical x-axis labels.
//------------------------------------------------------------------------------
static void DrawCoinTimeTrends(const std::vector<int> &Runs,
                               const std::vector<double> &Means,
                               const std::vector<double> &Sigmas) {
  if (Runs.empty()) return;
  const int N = (int)Runs.size();

  // Fixed bands for QA
  const double y1Min = 45.0;   // CTime mean (ns)
  const double y1Max = 55.0;
  const double sMin  = 0.0;    // CTime sigma (ns)
  const double sMax  = 1.0;

  // Categorical frame with run-number labels
  TCanvas *C = new TCanvas("C_CtTrends","C_CtTrends",1000,600);
  TH1F *frame = new TH1F("frame_ct_trend",
                         "COIN: CTime (ROC2) mean / sigma vs run;Run;CTime mean (ns)",
                         N, 0.0, (double)N);

  // Show only every 'step'-th run label:
  int step = 10;  // try 5, 10, 20 depending on N
  for (int i = 0; i < N; ++i) {
    if ((i % step) == 0) frame->GetXaxis()->SetBinLabel(i+1, Form("%d", Runs[i]));
    else                 frame->GetXaxis()->SetBinLabel(i+1, "");
  }

  frame->SetMinimum(y1Min);
  frame->SetMaximum(y1Max);
  //for (int i = 1; i <= N; ++i) frame->GetXaxis()->SetBinLabel(i, Form("%d", Runs[i-1]));
  frame->GetXaxis()->LabelsOption("v");
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->Draw("HIST");

  // Points at bin centers; scale sigma to left axis for drawing
  std::vector<double> X(N), Ymean(N), YsigScaled(N);
  for (int i = 0; i < N; ++i) {
    X[i] = i + 0.5;
    Ymean[i] = Means[i];
    YsigScaled[i] = y1Min + (Sigmas[i] - sMin) * (y1Max - y1Min) / (sMax - sMin);
  }

  TGraph *gMean  = new TGraph(N, X.data(), Ymean.data());
  TGraph *gSigma = new TGraph(N, X.data(), YsigScaled.data());

  // Markers
  gMean->SetMarkerStyle(20);  gMean->SetMarkerSize(1.1);  gMean->SetLineStyle(0);  gMean->SetLineWidth(0);
  gSigma->SetMarkerStyle(3);  gSigma->SetMarkerSize(1.1); gSigma->SetLineStyle(0); gSigma->SetLineWidth(0);

  gMean->Draw("P SAME");
  gSigma->Draw("P SAME");

  TGaxis *axisR = new TGaxis((double)N, y1Min, (double)N, y1Max, sMin, sMax, 510, "+L");
  axisR->SetTitle("CTime sigma (ns)");
  axisR->SetTitleOffset(1.2);
  axisR->Draw();

  TLegend *leg = new TLegend(0.12, 0.84, 0.24, 0.92);
  leg->AddEntry(gMean,  "mean",  "p");
  leg->AddEntry(gSigma, "sigma", "p");
  leg->Draw();

  // Save into coin folder
  C->SaveAs("coinPNGs/coin_ctime_trends_ROC2.png");

  delete leg; delete axisR; delete frame; delete gMean; delete gSigma; delete C;
}

//------------------------------------------------------------------------------
// ProcessOneRun: open file → prune branches → make plots → collect metrics
//------------------------------------------------------------------------------
static void ProcessOneRun(const TString &Spec, const TString &RootDir, int Run,
                          std::vector<int> &RunVec,
                          std::vector<double> &Means,
                          std::vector<double> &Sigmas,
                          bool ForCoinTime=false,
                          std::vector<double> *CoinMeans=nullptr,
                          std::vector<double> *CoinSigmas=nullptr){
  TString FileName = MakeFileName(Spec, Run);
  if (FileName.IsNull()) { std::cerr << "[WARN] Invalid Spec for run " << Run << ""; return; }
  TString Full = TString::Format("%s/%s", RootDir.Data(), FileName.Data());
  TFile *F = TFile::Open(Full, "READ");
  if (!F || F->IsZombie()) { std::cerr << "[WARN] Could not open " << Full << ""; return; }
  TTree *T = (TTree*) F->Get("T");
  if (!T) { std::cerr << "[WARN] Tree 'T' missing in " << Full << ""; F->Close(); delete F; return; }

  // Prune branches for this Spec
  SetBranchStatusesForSpec(T, Spec);

  // Plots
  if (Spec == "hms") {
    DrawBetaVsXfp(T, "hms", Run);
    double m,s,n; if (ComputeBetaMetrics(T, "hms", m,s,n)) { RunVec.push_back(Run); Means.push_back(m); Sigmas.push_back(s); }
  } else if (Spec == "shms") {
    DrawBetaVsXfp(T, "shms", Run);
    double m,s,n; if (ComputeBetaMetrics(T, "shms", m,s,n)) { RunVec.push_back(Run); Means.push_back(m); Sigmas.push_back(s); }
  } else if (Spec == "coin") {
    // HMS view
    DrawBetaVsXfp(T, "coin", Run);
    // SHMS view under coin selection: temporarily switch expr by calling DrawBetaVsXfp with shms
    // (Expr is set by Spec value; we call a dedicated SHMS draw under coin cuts)
    // Enable SHMS view explicitly by re-projecting with SHMS expression and coin cuts
    // Quick way: temporarily enable SHMS branches already on in coin
    {
      TH2D *H2 = new TH2D("H2_BetaVsXfp_SHMS","#beta vs x_{fp};x_{fp} (cm);#beta",80,-45,45,120,0.2,1.2);
      H2->Sumw2();
      T->Project("H2_BetaVsXfp_SHMS", "P_gtr_beta:P_dc_x_fp", BuildCuts("coin"));
      H2->SetDirectory(nullptr);

      TCanvas *C = new TCanvas("C_BetaVsXfp_SHMS","C_BetaVsXfp_SHMS",900,700);
      C->SetRightMargin(0.12);
      gStyle->SetOptStat(0);
      H2->Draw("COLZ");

      double Xmin = H2->GetXaxis()->GetXmin();
      double Xmax = H2->GetXaxis()->GetXmax();
      TLine L1(Xmin,0.95,Xmax,0.95), L2(Xmin,1.05,Xmax,1.05);
      L1.SetLineStyle(2); L2.SetLineStyle(2); L1.SetLineColor(kBlack); L2.SetLineColor(kBlack); L1.SetLineWidth(5); L2.SetLineWidth(5);
      L1.Draw("SAME"); L2.Draw("SAME");

      C->SaveAs(TString::Format("coinPNGs/coin_shms_run%d_beta_vs_xfp.png",Run));
      delete C; delete H2;
    }
    // Coin time 1D + metrics
    DrawCoinTime1D(T, Run);
    if (ForCoinTime && CoinMeans && CoinSigmas) {
      double m,s,n; if (ComputeCoinTimeMetrics(T, m,s,n)) { RunVec.push_back(Run); CoinMeans->push_back(m); CoinSigmas->push_back(s); }
    }
  }

  F->Close(); delete F;
}

//--------------------------------------------------------------
// Entry point: orchestrates per-run work and draws trend plots
//--------------------------------------------------------------
void hodo_calib_qc_batch(const char *Spec="", const char *RootDir="", const char *RunsList=""){
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("hmsPNGs",  true);
  gSystem->mkdir("shmsPNGs", true);
  gSystem->mkdir("coinPNGs", true);

  TString S(Spec?Spec:"");
  if (!(S=="hms" || S=="shms" || S=="coin")) { std::cerr << "[ERROR] Spec must be 'hms', 'shms', or 'coin'" << std::endl; return; }

  //std::vector<int> Runs = ParseRunsList(RunsList?RunsList:"");
  //std::vector<int> Runs = HMSRuns;
  //std::vector<int> Runs = SHMSRuns;
  std::vector<int> Runs = COINRuns;
  if (Runs.empty()) { std::cerr << "[INFO] No runs provided. Exiting."; return; }

  std::vector<int> TrendRuns; TrendRuns.reserve(Runs.size());
  std::vector<double> Means, Sigmas; Means.reserve(Runs.size()); Sigmas.reserve(Runs.size());
  std::vector<double> CoinMeans, CoinSigmas; CoinMeans.reserve(Runs.size()); CoinSigmas.reserve(Runs.size());

  for (int R : Runs) {
    if (S == "coin") {
      ProcessOneRun(S, RootDir, R, TrendRuns, Means, Sigmas, /*ForCoinTime=*/true, &CoinMeans, &CoinSigmas);
    } else {
      ProcessOneRun(S, RootDir, R, TrendRuns, Means, Sigmas);
    }
  }

  if (S == "hms" || S == "shms") {
    if (!TrendRuns.empty()) DrawBetaTrends(TrendRuns, Means, Sigmas, S);
  } else if (S == "coin") {
    if (!TrendRuns.empty()) DrawCoinTimeTrends(TrendRuns, CoinMeans, CoinSigmas);
  }
}
