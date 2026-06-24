// heep_EmVsDelta.C
// Plot missing energy Em vs deltas for a single HEEP run.
//
// Data (tree "T"):
//   Em   = H.kin.secondary.emiss
//   H dp = H.gtr.dp
//   P dp = P.gtr.dp
//
// Sim (tree "h10"):
//   Em   = Em_recon
//   hsdp = hsdelta
//   ssdp = ssdelta
//
// Cuts (same as your previous script):
//   DataCut = "(H.gtr.dp>-8) && (H.gtr.dp<8) && (P.gtr.dp>-10) && (P.gtr.dp<20)"
//   SimCut  = "(hsdelta>-8) && (hsdelta<8) && (ssdelta>-10) && (ssdelta<20)"
//
// Output PNG: PNGs/heep_EmVsDelta_run<run>.png
// Run this from:  RP_Scripts/heep_check/
//
// Run Command (single run): root -l -b -q 'heep_EmVsDelta.C+(23839)'

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

// Directories (same as before)
static const TString kDataDir = "/home/cdaq/users/rparvez/RP_Scripts/heep_check/HEEP_ROOTfiles/";
static const TString kSimDir  = "/home/cdaq/users/rparvez/simc_gfortran/worksim/heep_check/";

// File name patterns
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

void heep_EmVsDelta(int run)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  // Check PNG directory
  if (gSystem->AccessPathName("PNGs")) {
    std::cerr << "[WARNING] PNGs/ does not exist. Please create it:  mkdir PNGs\n";
  }

  // Open files
  TString dataFile = GetDataFileName(run);
  TString simFile  = GetSimFileName(run);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    std::cerr << "[ERROR] Cannot open data file: " << dataFile << std::endl;
    return;
  }
  TFile* fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) {
    std::cerr << "[ERROR] Cannot open sim file: " << simFile << std::endl;
    return;
  }

  TTree* tData = (TTree*) fData->Get(kDataTreeName);
  TTree* tSim  = (TTree*) fSim->Get(kSimTreeName);
  if (!tData || !tSim) {
    std::cerr << "[ERROR] Could not get trees T / h10\n";
    return;
  }

  // Binning & ranges
  const int nbins_dp  = 100;
  const int nbins_Em  = 100;
  const double Hdp_min = -8.0, Hdp_max = 8.0;
  const double Pdp_min = -10.0, Pdp_max = 20.0;
  const double Em_min  = -0.05, Em_max  = 0.15;

  // Histograms
  TH2D* h_Em_vs_Hdp_data = new TH2D("h_Em_vs_Hdp_data", "",
                                    nbins_dp, Hdp_min, Hdp_max,
                                    nbins_Em, Em_min, Em_max);
  TH2D* h_Em_vs_Pdp_data = new TH2D("h_Em_vs_Pdp_data", "",
                                    nbins_dp, Pdp_min, Pdp_max,
                                    nbins_Em, Em_min, Em_max);

  TH2D* h_Em_vs_hsdp_sim = new TH2D("h_Em_vs_hsdp_sim", "",
                                    nbins_dp, Hdp_min, Hdp_max,
                                    nbins_Em, Em_min, Em_max);
  TH2D* h_Em_vs_ssdp_sim = new TH2D("h_Em_vs_ssdp_sim", "",
                                    nbins_dp, Pdp_min, Pdp_max,
                                    nbins_Em, Em_min, Em_max);

  // Fill: Em vs delta
  tData->Draw("H.kin.secondary.emiss:H.gtr.dp>>h_Em_vs_Hdp_data", kDataCut, "goff");
  tData->Draw("H.kin.secondary.emiss:P.gtr.dp>>h_Em_vs_Pdp_data", kDataCut, "goff");

  tSim->Draw("Em_recon:hsdelta>>h_Em_vs_hsdp_sim", kSimCut, "goff");
  tSim->Draw("Em_recon:ssdelta>>h_Em_vs_ssdp_sim", kSimCut, "goff");

  // Set axis titles
  h_Em_vs_Hdp_data->GetXaxis()->SetTitle("H.gtr.dp");
  h_Em_vs_Hdp_data->GetYaxis()->SetTitle("Em (GeV)");

  h_Em_vs_Pdp_data->GetXaxis()->SetTitle("P.gtr.dp");
  h_Em_vs_Pdp_data->GetYaxis()->SetTitle("Em (GeV)");

  h_Em_vs_hsdp_sim->GetXaxis()->SetTitle("hsdelta");
  h_Em_vs_hsdp_sim->GetYaxis()->SetTitle("Em (GeV)");

  h_Em_vs_ssdp_sim->GetXaxis()->SetTitle("ssdelta");
  h_Em_vs_ssdp_sim->GetYaxis()->SetTitle("Em (GeV)");

  // Canvas
  TCanvas* c = new TCanvas(Form("c_EmVsDelta_%d", run),
                           Form("Em vs delta, Run %d", run),
                           1200, 900);
  c->Divide(2,2);

  // Pad 1: Data Em vs H.gtr.dp
  c->cd(1);
  gPad->SetRightMargin(0.15);
  h_Em_vs_Hdp_data->Draw("COLZ");
  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.12, 0.93, "Data: Em vs H.gtr.dp");

  // Pad 2: Data Em vs P.gtr.dp
  c->cd(2);
  gPad->SetRightMargin(0.15);
  h_Em_vs_Pdp_data->Draw("COLZ");
  lat.DrawLatex(0.12, 0.93, "Data: Em vs P.gtr.dp");

  // Pad 3: Sim Em vs hsdelta
  c->cd(3);
  gPad->SetRightMargin(0.15);
  h_Em_vs_hsdp_sim->Draw("COLZ");
  lat.DrawLatex(0.12, 0.93, "Sim: Em vs hsdelta");

  // Pad 4: Sim Em vs ssdelta
  c->cd(4);
  gPad->SetRightMargin(0.15);
  h_Em_vs_ssdp_sim->Draw("COLZ");
  lat.DrawLatex(0.12, 0.93, "Sim: Em vs ssdelta");

  c->Update();

  // Save PNG
  TString outName = Form("PNGs/heep_EmVsDelta_run%d.png", run);
  c->SaveAs(outName);
  std::cout << "[INFO] Saved " << outName << std::endl;
}
