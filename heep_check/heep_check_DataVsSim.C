// heep_check_DataVsSim.C
// Compare Data vs Sim for W, Em, Pmx, Pmy, Pmz for a single HEEP run.
//
// Data tree: "T"
// Sim tree:  "h10"
//
// Data variables (T):
//   W     -> P.kin.primary.W
//   Em    -> H.kin.secondary.emiss
//   Pmx   -> H.kin.secondary.pmiss_x
//   Pmy   -> H.kin.secondary.pmiss_y
//   Pmz   -> H.kin.secondary.pmiss_z
//
// Sim variables (h10):
//   W_recon
//   Em_recon
//   Pmx_recon
//   Pmy_recon
//   Pmz_recon
//
// Cuts:
//   DataCut = "(H.gtr.dp>-8) && (H.gtr.dp<8) && (P.gtr.dp>-10) && (P.gtr.dp<20)"
//   SimCut  = "(hsdelta>-8) && (hsdelta<8) && (ssdelta>-10) && (ssdelta<20)"
//
// File patterns:
//   Data: coin_replay_production_<run>_-1.root
//   Sim:  recon_hcana_coin_replay_production_<run>_-1.root
//
// Run from: RP_Scripts/heep_check/
// PNG output: PNGs/heep_DataVsSim_run<run>.png
//
// Run Command (single run): root -l -b -q 'heep_check_DataVsSim.C+(23839)'


#include <iostream>
#include <vector>
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

// Directories
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check/HEEP_ROOTfiles/";
static const TString kSimDir  = "/home/cdaq/users/rparvez/simc_gfortran/worksim/heep_check/";

// Patterns
static TString GetDataFileName(int run) {
  return kDataDir + Form("coin_replay_production_%d_-1.root", run);
}

static TString GetSimFileName(int run) {
  return kSimDir + Form("recon_hcana_coin_replay_production_%d_-1.root", run);
}

// Tree names
static const TString kDataTreeName = "T";
static const TString kSimTreeName  = "h10";

// Cuts
static const char* kDataCut = "(H.gtr.dp>-8) && (H.gtr.dp<8) && (P.gtr.dp>-10) && (P.gtr.dp<20)";
static const char* kSimCut  = "(hsdelta>-8) && (hsdelta<8) && (ssdelta>-10) && (ssdelta<20)";

// Variable definitions
struct VarInfo {
  TString name;
  TString data_expr;
  TString sim_expr;
  TString title;
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

  v.push_back({"Pmx",
               "H.kin.secondary.pmiss_x",
               "Pmx_recon",
               "Pmx (GeV/c)",
               100, -0.3, 0.3});

  v.push_back({"Pmy",
               "H.kin.secondary.pmiss_y",
               "Pmy_recon",
               "Pmy (GeV/c)",
               100, -0.3, 0.3});

  v.push_back({"Pmz",
               "H.kin.secondary.pmiss_z",
               "Pmz_recon",
               "Pmz (GeV/c)",
               100, -0.3, 0.3});

  return v;
}

void SetHistStyle(TH1D* hData, TH1D* hSim) {
  hData->SetMarkerStyle(20);
  hData->SetMarkerSize(0.8);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);

  hData->GetYaxis()->SetTitle("Counts");
  hData->GetYaxis()->SetTitleOffset(1.3);
  hData->GetXaxis()->SetTitleOffset(1.1);

  hSim->SetMarkerStyle(24);
  hSim->SetMarkerSize(0.7);
  hSim->SetLineColor(kRed+1);
  hSim->SetLineColor(kRed+1);

}

// =======================================================
// Main function
// =======================================================

void heep_check_DataVsSim(int run) {

  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  // Check PNG directory
  if (gSystem->AccessPathName("PNGs")) {
    std::cerr << "[WARNING] PNGs/ does not exist. Please create it: mkdir PNGs\n";
  }

  TString dataFile = GetDataFileName(run);
  TString simFile  = GetSimFileName(run);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) { std::cerr << "Cannot open data file\n"; return; }

  TFile* fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) { std::cerr << "Cannot open sim file\n"; return; }

  TTree* tData = (TTree*) fData->Get(kDataTreeName);
  TTree* tSim  = (TTree*) fSim->Get(kSimTreeName);

  std::vector<VarInfo> vars = GetVariables();

  TCanvas* c = new TCanvas(Form("c_heep_%d", run),
                           Form("HEEP Data vs Sim, Run %d", run),
                           1600, 1000);
  c->Divide(3, 2);
  gPad->SetLeftMargin(0.35);
  gPad->SetBottomMargin(0.16);

  for (int i = 0; i < (int)vars.size(); i++) {

    const VarInfo& V = vars[i];
    c->cd(i+1);

    TString hD = Form("hD_%s_%d", V.name.Data(), run);
    TString hS = Form("hS_%s_%d", V.name.Data(), run);

    TH1D* hData = new TH1D(hD, "", V.nbins, V.xmin, V.xmax);
    TH1D* hSim  = new TH1D(hS, "", V.nbins, V.xmin, V.xmax);

    hData->Sumw2();
    hSim->Sumw2();

    tData->Draw(Form("%s>>%s", V.data_expr.Data(), hD.Data()), kDataCut, "goff");
    tSim ->Draw(Form("%s>>%s", V.sim_expr.Data(),  hS.Data()), kSimCut,  "goff");

    double intData = hData->Integral();
    double intSim  = hSim->Integral();
    if (intSim > 0 && intData > 0) hSim->Scale(intData / intSim);

    SetHistStyle(hData, hSim);
    hData->GetXaxis()->SetTitle(V.title);

    double maxY = std::max(hData->GetMaximum(), hSim->GetMaximum());
    hData->SetMaximum(1.2 * maxY);

    // Fitting
    // Data Fit
    int    binMaxD  = hData->GetMaximumBin();
    int    binLowD  = std::max(1, binMaxD - 5);
    int    binHighD = std::min(hData->GetNbinsX(), binMaxD + 5);

    double xMinD = hData->GetXaxis()->GetBinLowEdge(binLowD);
    double xMaxD = hData->GetXaxis()->GetBinUpEdge(binHighD);

    TF1* fitData = new TF1(Form("fD_%s_%d", V.name.Data(), run),
                         "gaus", xMinD, xMaxD);
    fitData->SetLineColor(kBlack);
    fitData->SetLineWidth(1);
    hData->Fit(fitData, "QRN");  // Q=quiet, R=use range, N=don't draw automatically

    // Sim fit
    int    binMaxS  = hSim->GetMaximumBin();
    int    binLowS  = std::max(1, binMaxS - 2);
    int    binHighS = std::min(hSim->GetNbinsX(), binMaxS + 2);

    double xMinS = hSim->GetXaxis()->GetBinLowEdge(binLowS);
    double xMaxS = hSim->GetXaxis()->GetBinUpEdge(binHighS);

    TF1* fitSim = new TF1(Form("fS_%s_%d", V.name.Data(), run),
                        "gaus", xMinS, xMaxS);
    fitSim->SetLineColor(kRed+1);
    fitSim->SetLineWidth(1);
    hSim->Fit(fitSim, "QRN");

    // Extract mean and sigma from the fits
    double meanD  = fitData->GetParameter(1);
    double sigmaD = fitData->GetParameter(2);
    double meanS  = fitSim ->GetParameter(1);
    double sigmaS = fitSim ->GetParameter(2);

    hData->Draw("E1");
    hSim ->Draw("E1 SAME");

    // Draw fitted Gaussians on top
    fitData->SetLineColor(kBlack);
    fitData->SetLineWidth(2);
    fitSim->SetLineColor(kRed+1);
    fitSim->SetLineWidth(2);

    fitData->Draw("SAME");
    fitSim ->Draw("SAME");

    TLegend* L = new TLegend(0.55, 0.70, 0.88, 0.88);
    L->SetBorderSize(0);
    L->SetFillStyle(0);
    L->AddEntry(hData, "Data", "l");
    L->AddEntry(hSim,  "Sim",  "l");
    L->Draw();

    // mean info
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.15, 0.86, Form("Data: #mu=%.4f, #sigma=%.4f", meanD, sigmaD));
    lat.DrawLatex(0.15, 0.80, Form("Sim:  #mu=%.4f, #sigma=%.4f", meanS, sigmaS));
  }

  // Last pad
  c->cd(6);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.10);
  gPad->SetTopMargin(0.20);
  gPad->SetBottomMargin(0.20);

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.05);
  lat.DrawLatex(0.10, 0.82, "HEEP Data vs Sim");
  lat.DrawLatex(0.10, 0.70, Form("Run %d", run));
  lat.SetTextSize(0.03);
  lat.DrawLatex(0.10, 0.55, Form("DataCut: %s", kDataCut));
  lat.DrawLatex(0.10, 0.45, Form("SimCut : %s", kSimCut));

  TString outName = Form("PNGs/heep_DataVsSim_run%d.png", run);
  c->SaveAs(outName);

  std::cout << "[INFO] Saved: " << outName << std::endl;
}
