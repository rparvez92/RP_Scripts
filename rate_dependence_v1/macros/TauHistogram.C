// File: rate_dependence_v1/macros/TauHistogram.C
//
// Run from rate_dependence_v1/:
//   root -l -b -q 'macros/TauHistogram.C()'
// or specify a different results directory:
//   root -l -b -q 'macros/TauHistogram.C("results")'
//
// Produces:
//   <resultsDir>/tau/tau_hist_all.png
//
// It scans recursively for files named:
//   yield_ratio_vs_trigger_shms34.csv
// and extracts tau_ns from the CSV header lines that start with '#'.

#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TCollection.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

namespace {
  // --------------------------------------------------------------------------
  // Recursively collect CSV paths matching the target filename
  // --------------------------------------------------------------------------
  void CollectFilesRec(const TString& dirPath,
                       const TString& targetFileName,
                       std::vector<TString>& outFiles)
  {
    TSystemDirectory dir(dirPath, dirPath);
    TList* files = dir.GetListOfFiles();
    if (!files) return;

    TIter next(files);
    while (TSystemFile* f = (TSystemFile*) next()) {
      TString name = f->GetName();
      if (name == "." || name == "..") continue;

      TString fullPath = dirPath + "/" + name;

      if (f->IsDirectory()) {
        CollectFilesRec(fullPath, targetFileName, outFiles);
      } else {
        if (name == targetFileName) {
          outFiles.push_back(fullPath);
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // Parse tau_ns from CSV header. Returns true if found & valid.
  // Guardrails: finite, >0, and <1e6 ns (to avoid parse garbage).
  // --------------------------------------------------------------------------
  bool ParseTauNsFromHeader(const TString& csvPath, double& tauNs)
  {
    std::ifstream in(csvPath.Data());
    if (!in.is_open()) return false;

    std::string line;
    bool found = false;
    double tau = std::nan("");

    while (std::getline(in, line)) {
      if (line.size() == 0) continue;
      if (line[0] != '#') break; // header ended

      // Look for: "# tau_ns: <value>"
      auto pos = line.find("tau_ns:");
      if (pos != std::string::npos) {
        std::string rhs = line.substr(pos + std::string("tau_ns:").size());
        try {
          tau = std::stod(rhs);
          found = true;
        } catch (...) {
          found = false;
        }
        break;
      }
    }

    if (!found) return false;
    if (!std::isfinite(tau)) return false;
    if (tau <= 0.0) return false;
    if (tau > 1e6) return false;

    tauNs = tau;
    return true;
  }

  // --------------------------------------------------------------------------
  // Extract pass and target from path using directory tokens.
  // --------------------------------------------------------------------------
  bool ExtractPassAndTarget(const TString& path, int& pass, TString& target)
  {
    pass = 0;
    target = "";

    if (path.Contains("/pass4/")) pass = 4;
    else if (path.Contains("/pass5/")) pass = 5;

    if (path.Contains("/LH2/")) target = "LH2";
    else if (path.Contains("/LD2/")) target = "LD2";

    return (pass != 0 && (target == "LH2" || target == "LD2"));
  }

  // --------------------------------------------------------------------------
  // Helper to add small stats box (N, mean, RMS, under/overflow)
  // --------------------------------------------------------------------------
  void DrawStatsBox(TH1D* h)
  {
    const int nb = h->GetNbinsX();
    const double under = h->GetBinContent(0);
    const double over  = h->GetBinContent(nb + 1);
    const double N     = h->GetEntries();
    const double mean  = (N > 0) ? h->GetMean() : 0.0;
    const double rms   = (N > 0) ? h->GetRMS()  : 0.0;

    TPaveText* pt = new TPaveText(0.62, 0.70, 0.95, 0.93, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.035);

    pt->AddText(Form("N = %.0f", N));
    pt->AddText(Form("Mean = %.2f ns", mean));
    pt->AddText(Form("RMS = %.2f ns", rms));
    pt->AddText(Form("Under = %.0f", under));
    pt->AddText(Form("Over  = %.0f", over));
    pt->Draw("same");
  }
} // end anon namespace

// ============================================================================
// Main macro
// ============================================================================
void TauHistogram(const char* resultsDir = "results")
{
  gStyle->SetOptStat(0);

  const TString kResultsDir = resultsDir;
  const TString kTargetName = "yield_ratio_vs_trigger_shms34.csv";

  // Collect all matching CSV files under resultsDir
  std::vector<TString> csvFiles;
  CollectFilesRec(kResultsDir, kTargetName, csvFiles);

  std::cout << "[TauHistogram] Found " << csvFiles.size()
            << " candidate CSV files named '" << kTargetName << "' under '"
            << kResultsDir << "'.\n";

  // Buckets: (pass, target)
  // Index mapping:
  // 0: pass4 LH2, 1: pass4 LD2, 2: pass5 LH2, 3: pass5 LD2
  auto idxOf = [](int pass, const TString& target) -> int {
    if (pass == 4 && target == "LH2") return 0;
    if (pass == 4 && target == "LD2") return 1;
    if (pass == 5 && target == "LH2") return 2;
    if (pass == 5 && target == "LD2") return 3;
    return -1;
  };

  std::vector<double> taus[4];

  int nRead = 0, nSkippedNoTau = 0, nSkippedBadCat = 0, nSkippedIO = 0;

  for (const auto& p : csvFiles) {
    int pass = 0;
    TString target;

    if (!ExtractPassAndTarget(p, pass, target)) {
      std::cerr << "[WARN] Cannot categorize by pass/target from path, skipping: "
                << p << "\n";
      nSkippedBadCat++;
      continue;
    }

    double tauNs = std::nan("");
    bool ok = ParseTauNsFromHeader(p, tauNs);
    if (!ok) {
      // Distinguish file open failure vs missing tau
      std::ifstream test(p.Data());
      if (!test.is_open()) {
        std::cerr << "[WARN] Cannot open file, skipping: " << p << "\n";
        nSkippedIO++;
      } else {
        std::cerr << "[WARN] Missing/invalid tau_ns in header, skipping: " << p << "\n";
        nSkippedNoTau++;
      }
      continue;
    }

    const int idx = idxOf(pass, target);
    if (idx < 0) {
      std::cerr << "[WARN] Unexpected pass/target combination, skipping: " << p << "\n";
      nSkippedBadCat++;
      continue;
    }

    taus[idx].push_back(tauNs);
    nRead++;
  }

  std::cout << "[TauHistogram] Read tau from " << nRead << " files.\n"
            << "  Skipped (bad category): " << nSkippedBadCat << "\n"
            << "  Skipped (missing/invalid tau): " << nSkippedNoTau << "\n"
            << "  Skipped (I/O open failure): " << nSkippedIO << "\n";

  // Create output directory: <resultsDir>/tau/
  const TString outDir = kResultsDir + "/tau";
  if (gSystem->AccessPathName(outDir)) {
    int mk = gSystem->mkdir(outDir, kTRUE);
    if (mk != 0) {
      std::cerr << "[ERROR] Could not create output directory: " << outDir << "\n";
      return;
    }
  }

  const TString outPng = outDir + "/tau_hist_all.png";

  // Histogram settings: 0-400 ns, 10 ns bins => 40 bins
  const int nbins = 8;
  const double xmin = 0.0, xmax = 400.0;

  TH1D* h[4];
  h[0] = new TH1D("h_p4_LH2", "Pass4, LH2;#tau (ns);Counts", nbins, xmin, xmax);
  h[1] = new TH1D("h_p4_LD2", "Pass4, LD2;#tau (ns);Counts", nbins, xmin, xmax);
  h[2] = new TH1D("h_p5_LH2", "Pass5, LH2;#tau (ns);Counts", nbins, xmin, xmax);
  h[3] = new TH1D("h_p5_LD2", "Pass5, LD2;#tau (ns);Counts", nbins, xmin, xmax);

  for (int i = 0; i < 4; i++) {
    h[i]->Sumw2(false);
    for (double v : taus[i]) h[i]->Fill(v);
  }

  // 2x2 canvas (Option A)
  TCanvas* c = new TCanvas("c_tau_hist", "Tau histogram", 1200, 900);
  c->Divide(2, 2);

  // Order: (Pass4 LH2) (Pass4 LD2)
  //        (Pass5 LH2) (Pass5 LD2)
  for (int pad = 1; pad <= 4; pad++) {
    c->cd(pad);
    gPad->SetGrid(0, 0);
    gPad->SetTicks(1, 1);

    TH1D* hh = h[pad - 1];
    hh->SetLineWidth(2);
    hh->Draw("hist");

    DrawStatsBox(hh);

    // Small note about binning/range
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.12, 0.92, "Range: 0-400 ns, bin: 50 ns");
  }

  c->SaveAs(outPng);

  std::cout << "[TauHistogram] Saved: " << outPng << "\n";
}
