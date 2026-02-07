// heep_check_WEmVsDpWindow.C
// ------------------------------------------------------------
// Make a 2x2 comparison plot for a single HEEP run:
//
//   (Top row)    W
//   (Bottom row) Emiss
//
//   (Left col)  Base cuts only  (no extra P.gtr.dp window)
//   (Right col) Base cuts + a chosen P.gtr.dp (ssdelta) window
//
// Normalization choice: normalize BOTH Data and Sim to unity (shape-only).
// Fits: Gaussian fit around peak bin ± 5 bins (for both Data and Sim).
// Empty-hist protection: annotate and skip fits if empty.
//
// Run example:
//   root -l -b -q 'heep_check_WEmVsDpWindow.C+(23839,0.5,2.0)'
// or pick from a predefined dp-window list by index:
//   root -l -b -q 'heep_check_WEmVsDpWindow.C+(23839,5)'
// ------------------------------------------------------------

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"

// Directories (mirroring your existing macro)
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check/HEEP_ROOTfiles/";
static const TString kSimDir  = "/home/cdaq/users/rparvez/simc_gfortran/worksim/heep_check/";

// File patterns (mirroring your existing macro)
static TString GetDataFileName(int run) {
  return kDataDir + Form("coin_replay_production_%d_-1.root", run);
}
static TString GetSimFileName(int run) {
  return kSimDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
}

// Tree names
static const TString kDataTreeName = "T";
static const TString kSimTreeName  = "h10";

// Base cuts (as you specified)
static const char* kDataCutBase =
  "(H.gtr.dp>-8) && (H.gtr.dp<8) && (P.gtr.dp>-10) && (P.gtr.dp<20)";
static const char* kSimCutBase  =
  "(hsdelta>-8) && (hsdelta<8) && (ssdelta>-10) && (ssdelta<20)";

// Variables (only W and Em, as requested)
struct VarInfo {
  TString name;      // internal name
  TString data_expr; // Data tree expression
  TString sim_expr;  // Sim tree expression
  TString xTitle;    // x-axis title
  int nbins;
  double xmin;
  double xmax;
  TString deltaLabel; // label for delta annotation (e.g., "#DeltaW")
};

static std::vector<VarInfo> GetVariables_W_Em() {
  std::vector<VarInfo> v;
  v.push_back({"W",
               "P.kin.primary.W",
               "W_recon",
               "W (GeV)",
               120, 0.80, 1.10,
               "#DeltaW"});
  v.push_back({"Em",
               "H.kin.secondary.emiss",
               "Em_recon",
               "E_{miss} (GeV)",
               120, -0.06, 0.06,
               "#DeltaE_{miss}"});
  return v;
}

// A clean set of dp windows to scan later (implemented here as a selectable list).
// These are non-overlapping (mostly) and include a few “interesting” tight windows near 0–2.
static std::vector<std::pair<double,double>> GetDpWindows() {
  std::vector<std::pair<double,double>> w;

  // Coarse coverage across full range [-10, 20]
  w.push_back({-10.0, -6.0});
  w.push_back({ -6.0, -3.0});
  w.push_back({ -3.0,  0.0});
  w.push_back({  0.0,  0.5});
  w.push_back({  0.5,  2.0});   // the “example” window you showed
  w.push_back({  2.0,  5.0});
  w.push_back({  5.0, 10.0});
  w.push_back({ 10.0, 20.0});

  // Optional “fine” windows near the elastic peak region (useful if you want more granularity later)
  // (Kept at the end so indices for the coarse set remain stable.)
  w.push_back({ 0.0, 1.0});
  w.push_back({ 1.0, 2.0});
  w.push_back({-1.0, 1.0});

  return w;
}

static void PrintDpWindows() {
  auto w = GetDpWindows();
  std::cout << "\n[INFO] Predefined dp windows (index : lo hi)\n";
  for (size_t i=0; i<w.size(); ++i) {
    std::cout << "  " << i << " : " << w[i].first << "  " << w[i].second << "\n";
  }
  std::cout << std::endl;
}

static void SetHistStyle(TH1D* hData, TH1D* hSim) {
  // Data
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.85);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);

  // Sim
  hSim->SetMarkerStyle(24);
  hSim->SetMarkerSize(0.75);
  hSim->SetMarkerColor(kRed+1);
  hSim->SetLineColor(kRed+1);

  // Axes
  hData->GetYaxis()->SetTitle("Normalized Yield");
  hData->GetYaxis()->SetTitleOffset(1.35);
  hData->GetXaxis()->SetTitleOffset(1.10);
}

struct FitResult {
  bool ok = false;
  double mean = 0.0;
  double sigma = 0.0;
  TF1* f = nullptr;
};

// Normalize histogram to unity integral (shape-only). Returns false if empty.
static bool NormalizeToUnity(TH1D* h) {
  double integral = h->Integral();
  if (integral <= 0) return false;
  h->Scale(1.0 / integral);
  return true;
}

// Gaussian fit around peak bin ± nBins. Returns FitResult with ok=false if cannot fit.
static FitResult FitPeakGaus(TH1D* h, const TString& fname, int nBinsAroundPeak, int lineColor) {
  FitResult r;
  if (!h) return r;
  if (h->Integral() <= 0 || h->GetMaximum() <= 0) return r;

  const int nb = h->GetNbinsX();
  int binMax = h->GetMaximumBin();
  int binLo  = std::max(1, binMax - nBinsAroundPeak);
  int binHi  = std::min(nb, binMax + nBinsAroundPeak);

  double xMin = h->GetXaxis()->GetBinLowEdge(binLo);
  double xMax = h->GetXaxis()->GetBinUpEdge(binHi);

  TF1* f = new TF1(fname, "gaus", xMin, xMax);
  f->SetLineColor(lineColor);
  f->SetLineWidth(2);

  // Quiet fit, use range, do not auto-draw
  int status = h->Fit(f, "QRN");
  if (status != 0) {
    return r;
  }

  r.ok = true;
  r.mean  = f->GetParameter(1);
  r.sigma = f->GetParameter(2);
  r.f = f;
  return r;
}

static void DrawEmptyNotice(const TString& msg) {
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.07);
  lat.SetTextColor(kRed+1);
  lat.DrawLatex(0.15, 0.55, msg);
}

// ------------------------------------------------------------
// Core worker (dpLo, dpHi specified explicitly)
// ------------------------------------------------------------
void heep_check_WEmVsDpWindow(int run, double dpLo, double dpHi) {

  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  // Ensure PNGs/ exists (mirror your macro’s behavior)
  if (gSystem->AccessPathName("PNGs")) {
    std::cerr << "[WARNING] PNGs/ does not exist. Please create it: mkdir PNGs\n";
  }

  TString dataFile = GetDataFileName(run);
  TString simFile  = GetSimFileName(run);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[ERROR] Cannot open data file: " << dataFile << "\n";
    return;
  }
  TFile* fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) {
    std::cerr << "[ERROR] Cannot open sim file: " << simFile << "\n";
    return;
  }

  TTree* tData = (TTree*)fData->Get(kDataTreeName);
  TTree* tSim  = (TTree*)fSim->Get(kSimTreeName);
  if (!tData) { std::cerr << "[ERROR] Missing data tree: " << kDataTreeName << "\n"; return; }
  if (!tSim)  { std::cerr << "[ERROR] Missing sim tree: "  << kSimTreeName  << "\n"; return; }

  // Cuts:
  // Left column: base cuts only
  TString dataCutA = kDataCutBase;
  TString simCutA  = kSimCutBase;

  // Right column: base + dp window
  TString dataDpWin = Form("(P.gtr.dp>%.4f) && (P.gtr.dp<%.4f)", dpLo, dpHi);
  TString simDpWin  = Form("(ssdelta>%.4f) && (ssdelta<%.4f)", dpLo, dpHi);

  TString dataCutB = Form("(%s) && (%s)", kDataCutBase, dataDpWin.Data());
  TString simCutB  = Form("(%s) && (%s)", kSimCutBase,  simDpWin.Data());

  // Canvas: 2x2
  TCanvas* c = new TCanvas(Form("c_WEm_dp_%d", run),
                           Form("HEEP W/Em vs dp-window, Run %d", run),
                           1600, 1000);
  c->Divide(2, 2);

  // Variables
  auto vars = GetVariables_W_Em();

  // We will fill all hists first so we can harmonize maxima per variable (left vs right)
  struct HistPack {
    TH1D* hDA = nullptr; TH1D* hSA = nullptr;
    TH1D* hDB = nullptr; TH1D* hSB = nullptr;
    bool okDA=false, okSA=false, okDB=false, okSB=false;
    double maxY = 0.0;
  };
  std::vector<HistPack> packs(vars.size());

  for (size_t iv=0; iv<vars.size(); ++iv) {
    const auto& V = vars[iv];

    TString nDA = Form("hDA_%s_%d", V.name.Data(), run);
    TString nSA = Form("hSA_%s_%d", V.name.Data(), run);
    TString nDB = Form("hDB_%s_%d", V.name.Data(), run);
    TString nSB = Form("hSB_%s_%d", V.name.Data(), run);

    packs[iv].hDA = new TH1D(nDA, "", V.nbins, V.xmin, V.xmax);
    packs[iv].hSA = new TH1D(nSA, "", V.nbins, V.xmin, V.xmax);
    packs[iv].hDB = new TH1D(nDB, "", V.nbins, V.xmin, V.xmax);
    packs[iv].hSB = new TH1D(nSB, "", V.nbins, V.xmin, V.xmax);

    packs[iv].hDA->Sumw2(); packs[iv].hSA->Sumw2();
    packs[iv].hDB->Sumw2(); packs[iv].hSB->Sumw2();

    // Fill
    tData->Draw(Form("%s>>%s", V.data_expr.Data(), nDA.Data()), dataCutA, "goff");
    tSim ->Draw(Form("%s>>%s", V.sim_expr.Data(),  nSA.Data()), simCutA,  "goff");
    tData->Draw(Form("%s>>%s", V.data_expr.Data(), nDB.Data()), dataCutB, "goff");
    tSim ->Draw(Form("%s>>%s", V.sim_expr.Data(),  nSB.Data()), simCutB,  "goff");

    // Normalize to unity (shape-only)
    packs[iv].okDA = NormalizeToUnity(packs[iv].hDA);
    packs[iv].okSA = NormalizeToUnity(packs[iv].hSA);
    packs[iv].okDB = NormalizeToUnity(packs[iv].hDB);
    packs[iv].okSB = NormalizeToUnity(packs[iv].hSB);

    // Determine max for consistent y-scale across left/right for this variable
    double m = 0.0;
    if (packs[iv].hDA) m = std::max(m, packs[iv].hDA->GetMaximum());
    if (packs[iv].hSA) m = std::max(m, packs[iv].hSA->GetMaximum());
    if (packs[iv].hDB) m = std::max(m, packs[iv].hDB->GetMaximum());
    if (packs[iv].hSB) m = std::max(m, packs[iv].hSB->GetMaximum());
    packs[iv].maxY = (m > 0 ? 1.25*m : 1.0);
  }

  // Draw pads:
  // Pad indices: (col,row): (1,1)=top-left, (2,1)=top-right, (1,2)=bottom-left, (2,2)=bottom-right
  for (size_t iv=0; iv<vars.size(); ++iv) {
    const auto& V = vars[iv];
    auto& P = packs[iv];

    // Determine row: iv=0 -> top row (W), iv=1 -> bottom row (Em)
    int row = (int)iv + 1;

    // --------------------------
    // Left pad (base cuts)
    // --------------------------
    int padLeft = (row-1)*2 + 1;
    c->cd(padLeft);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);

    SetHistStyle(P.hDA, P.hSA);
    P.hDA->GetXaxis()->SetTitle(V.xTitle);
    P.hDA->SetMaximum(P.maxY);

    // If both empty, draw frame + message
    if (!P.okDA && !P.okSA) {
      P.hDA->Draw(); // draws axes
      DrawEmptyNotice("Empty after cuts");
    } else {
      P.hDA->Draw("E1");
      P.hSA->Draw("E1 SAME");

      // Fits (peak ±5 bins, as you requested)
      FitResult fD = FitPeakGaus(P.hDA, Form("fDA_%s_%d", V.name.Data(), run), 5, kBlack);
      FitResult fS = FitPeakGaus(P.hSA, Form("fSA_%s_%d", V.name.Data(), run), 5, kRed+1);

      if (fD.ok) fD.f->Draw("SAME");
      if (fS.ok) fS.f->Draw("SAME");

      // Legend
      TLegend* L = new TLegend(0.62, 0.72, 0.92, 0.90);
      L->SetBorderSize(0);
      L->SetFillStyle(0);
      L->AddEntry(P.hDA, "Data", "lp");
      L->AddEntry(P.hSA, "Sim",  "lp");
      L->Draw();

      // Delta mean annotation (only if both fits ok)
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.07);
      lat.SetTextColor(kRed+1);

      if (fD.ok && fS.ok) {
        double dmu = fD.mean - fS.mean;
        lat.DrawLatex(0.15, 0.82, Form("%s = %.5f", V.deltaLabel.Data(), dmu));
      } else {
        lat.DrawLatex(0.15, 0.82, Form("%s = (fit failed)", V.deltaLabel.Data()));
      }
    }

    // Pad header (base)
    {
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.05);
      lat.SetTextColor(kBlack);
      lat.DrawLatex(0.38, 0.94, "No extra P.gtr.dp window");
    }

    // --------------------------
    // Right pad (dp window)
    // --------------------------
    int padRight = (row-1)*2 + 2;
    c->cd(padRight);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);

    SetHistStyle(P.hDB, P.hSB);
    P.hDB->GetXaxis()->SetTitle(V.xTitle);
    P.hDB->SetMaximum(P.maxY);

    if (!P.okDB && !P.okSB) {
      P.hDB->Draw();
      DrawEmptyNotice("Empty after dp window");
    } else {
      P.hDB->Draw("E1");
      P.hSB->Draw("E1 SAME");

      FitResult fD = FitPeakGaus(P.hDB, Form("fDB_%s_%d", V.name.Data(), run), 5, kBlack);
      FitResult fS = FitPeakGaus(P.hSB, Form("fSB_%s_%d", V.name.Data(), run), 5, kRed+1);

      if (fD.ok) fD.f->Draw("SAME");
      if (fS.ok) fS.f->Draw("SAME");

      TLegend* L = new TLegend(0.62, 0.72, 0.92, 0.90);
      L->SetBorderSize(0);
      L->SetFillStyle(0);
      L->AddEntry(P.hDB, "Data", "lp");
      L->AddEntry(P.hSB, "Sim",  "lp");
      L->Draw();

      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.07);
      lat.SetTextColor(kRed+1);

      if (fD.ok && fS.ok) {
        double dmu = fD.mean - fS.mean;
        lat.DrawLatex(0.15, 0.82, Form("%s = %.5f", V.deltaLabel.Data(), dmu));
      } else {
        lat.DrawLatex(0.15, 0.82, Form("%s = (fit failed)", V.deltaLabel.Data()));
      }
    }

    // Pad header (window)
    {
      TLatex lat;
      lat.SetNDC();
      lat.SetTextSize(0.05);
      lat.SetTextColor(kBlack);
      lat.DrawLatex(0.20, 0.94, Form("%.2f < P.gtr.dp < %.2f", dpLo, dpHi));
    }
  }

  // Global run/cut note (small, bottom-left of canvas area)
  c->cd(1);
  TLatex gl;
  gl.SetNDC();
  gl.SetTextSize(0.04);
  gl.SetTextColor(kBlack);
  gl.DrawLatex(0.13, 0.06, Form("Run %d | Base: %s", run, kDataCutBase));

  TString outName = Form("PNGs/heep_WEmVsDpWindow_run%d_dp%.2f_%.2f.png", run, dpLo, dpHi);
  //outName.ReplaceAll("-", "m").ReplaceAll(".", "p"); // safer filename
  c->SaveAs(outName);

  std::cout << "[INFO] Saved: " << outName << "\n";
}

// ------------------------------------------------------------
// Convenience overload: choose dp window by index from the list
// ------------------------------------------------------------
void heep_check_WEmVsDpWindow(int run, int winIndex) {
  PrintDpWindows();
  auto wins = GetDpWindows();
  if (winIndex < 0 || winIndex >= (int)wins.size()) {
    std::cerr << "[ERROR] winIndex out of range. Use one of the indices printed above.\n";
    return;
  }
  heep_check_WEmVsDpWindow(run, wins[winIndex].first, wins[winIndex].second);
}
