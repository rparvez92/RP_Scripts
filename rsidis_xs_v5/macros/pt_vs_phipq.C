// pt_vs_phipq.C
//
// Coincidence-data acceptance maps in pT vs phi_pq for nominal z settings.
//
// This diagnostic combines all currently defined signal coin-data runs for a
// given pass and nominal z, independent of target and pion charge. Coincidence
// windows are determined run-by-run with CoincidenceRandomSubtraction.h, using
// the same peak/window logic as TableCoinXsec.C. The plotted map is
// random-subtracted bin-by-bin:
//   Hsub(pT,phi) = Hcoin(pT,phi) - average_i Hrand_i(pT,phi).
//
// Output:
//   results/PNGs/pt_vs_phipq_z0p36.png
//   results/PNGs/pt_vs_phipq_z0p5.png
//   results/PNGs/pt_vs_phipq_z0p67.png
//   results/PNGs/pt_vs_phipq_z0p9.png
//   results/tables/pt_vs_phipq.csv
//
// Run from rsidis_xs_v5:
//   root -l -b -q 'macros/pt_vs_phipq.C()'

#include "../include/CoincidenceRandomSubtraction.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace {

struct Bucket {
  std::set<int> runs;
  std::set<std::string> manifests;
  std::set<std::string> groups;
};

static bool EndsWith(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string NormalizeSlashes(std::string s) {
  std::replace(s.begin(), s.end(), '\\', '/');
  while (s.size() > 1 && s.find("//") != std::string::npos) {
    size_t p = s.find("//");
    s.erase(p, 1);
  }
  while (s.size() > 1 && s.back() == '/') s.pop_back();
  return s;
}

static std::string Dirname(const std::string& p) {
  std::string s = NormalizeSlashes(p);
  auto pos = s.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return s.substr(0, pos);
}

static std::string Basename(const std::string& p) {
  std::string s = NormalizeSlashes(p);
  auto pos = s.find_last_of('/');
  if (pos == std::string::npos) return s;
  return s.substr(pos + 1);
}

static void MkdirP(const std::string& path) {
  if (!path.empty()) gSystem->mkdir(path.c_str(), true);
}

static void CollectFilesRecursive(const std::string& dir,
                                  const std::string& suffix,
                                  std::vector<std::string>& out) {
  DIR* dp = opendir(dir.c_str());
  if (!dp) return;
  while (dirent* ent = readdir(dp)) {
    std::string name = ent->d_name;
    if (name == "." || name == "..") continue;
    std::string path = NormalizeSlashes(dir + "/" + name);
    if (ent->d_type == DT_DIR) {
      CollectFilesRecursive(path, suffix, out);
    } else if (EndsWith(name, suffix)) {
      out.push_back(path);
    } else if (ent->d_type == DT_UNKNOWN) {
      // Some filesystems do not populate d_type reliably.
      if (gSystem->AccessPathName(path.c_str(), kReadPermission) == false &&
          EndsWith(name, suffix)) {
        out.push_back(path);
      }
    }
  }
  closedir(dp);
}

static std::vector<int> ReadRunList(const std::string& path) {
  std::vector<int> runs;
  std::ifstream f(path);
  if (!f.is_open()) return runs;
  std::string line;
  while (std::getline(f, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    line.erase(std::remove_if(line.begin(), line.end(),
                              [](unsigned char c){ return std::isspace(c); }),
               line.end());
    if (line.empty()) continue;
    runs.push_back(std::atoi(line.c_str()));
  }
  return runs;
}

static std::vector<std::string> SplitPath(const std::string& p) {
  std::vector<std::string> parts;
  std::string cur;
  for (char c : NormalizeSlashes(p)) {
    if (c == '/') {
      if (!cur.empty()) parts.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) parts.push_back(cur);
  return parts;
}

static bool ExtractPassAndZ(const std::string& manifestPath,
                            std::string& pass,
                            std::string& z) {
  auto parts = SplitPath(manifestPath);
  for (size_t i = 0; i < parts.size(); ++i) {
    if (parts[i] == "settings" && i + 4 < parts.size()) {
      pass = parts[i + 1];
      z = parts[i + 4];
      return pass.rfind("pass", 0) == 0 && z.rfind("z", 0) == 0;
    }
  }
  return false;
}

static std::string DataRootPath(const std::string& dataRoot, int run) {
  char buf[1024];
  std::snprintf(buf, sizeof(buf), "%s/skimmed_coin_replay_production_%d_-1.root",
                dataRoot.c_str(), run);
  return NormalizeSlashes(buf);
}

static void ReadGroupList(const std::string& groupPath,
                          std::map<std::string, Bucket>& buckets) {
  std::ifstream f(groupPath);
  if (!f.is_open()) {
    std::cerr << "WARNING: cannot open group list: " << groupPath << "\n";
    return;
  }
  const std::string baseDir = Dirname(groupPath);
  std::string line;
  while (std::getline(f, line)) {
    auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    std::istringstream iss(line);
    std::string label, manifestRel;
    if (!(iss >> label >> manifestRel)) continue;

    std::string manifest = manifestRel;
    if (!manifest.empty() && manifest[0] != '/') {
      manifest = NormalizeSlashes(baseDir + "/" + manifest);
    }
    std::string pass, z;
    if (!ExtractPassAndZ(manifest, pass, z)) {
      std::cerr << "WARNING: cannot infer pass/z from manifest: " << manifest << "\n";
      continue;
    }
    std::string leafDir = Dirname(manifest);
    auto runs = ReadRunList(leafDir + "/runs_signal.txt");
    std::string key = pass + "|" + z;
    Bucket& b = buckets[key];
    b.groups.insert(groupPath);
    b.manifests.insert(manifest);
    for (int run : runs) b.runs.insert(run);
  }
}

static std::string ZPretty(const std::string& z) {
  if (z == "z0p36") return "z = 0.36";
  if (z == "z0p5")  return "z = 0.50";
  if (z == "z0p67") return "z = 0.67";
  if (z == "z0p9")  return "z = 0.90";
  return z;
}

static void SetHistStyle(TH2D* h) {
  if (!h) return;
  h->SetStats(false);
  h->GetXaxis()->SetLabelSize(0.0);
  h->GetYaxis()->SetLabelSize(0.0);
  h->GetXaxis()->SetTickLength(0.0);
  h->GetYaxis()->SetTickLength(0.0);
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("");
}

static void DrawPolarGrid(double ptMax) {
  const int gridColor = kGray + 1;
  TLatex lat;
  lat.SetTextSize(0.028);
  lat.SetTextAlign(22);
  lat.SetTextColor(kBlack);

  for (double r = 0.1; r <= ptMax + 1e-9; r += 0.1) {
    TEllipse* e = new TEllipse(0.0, 0.0, r, r);
    e->SetFillStyle(0);
    e->SetLineStyle(3);
    e->SetLineColor(gridColor);
    e->SetLineWidth(1);
    e->Draw("same");
    lat.SetTextAlign(12);
    lat.DrawLatex(r / std::sqrt(2.0) + 0.01, r / std::sqrt(2.0),
                  Form("p_{T}=%.1f", r));
  }

  TLine line;
  line.SetLineStyle(3);
  line.SetLineColor(gridColor);
  line.SetLineWidth(1);
  line.DrawLine(-ptMax, 0.0, ptMax, 0.0);
  line.DrawLine(0.0, -ptMax, 0.0, ptMax);

  lat.SetTextAlign(22);
  lat.DrawLatex(ptMax + 0.035, 0.0, "0#circ");
  lat.DrawLatex(0.0, ptMax + 0.035, "90#circ");
  lat.DrawLatex(-ptMax - 0.040, 0.0, "180#circ");
  lat.DrawLatex(0.0, -ptMax - 0.035, "270#circ");
}

static void SavePolarPlot(const std::string& pngPath,
                          const std::string& z,
                          TH2D* h4,
                          TH2D* h5,
                          Long64_t n4,
                          Long64_t n5,
                          double ptMax) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kGreyScale);
  TGaxis::SetMaxDigits(3);

  TCanvas c("c_pt_vs_phipq", "pT vs phi_pq", 1500, 850);
  c.SetFillColor(kWhite);
  c.Divide(2, 1, 0.045, 0.02);

  struct PadSpec { TH2D* h; const char* title; Long64_t n; };
  PadSpec pads[2] = {{h4, "pass4", n4}, {h5, "pass5", n5}};
  for (int i = 0; i < 2; ++i) {
    c.cd(i + 1);
    gPad->SetLeftMargin(0.06);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.08);
    gPad->SetTicks(1, 1);
    if (pads[i].h) {
      SetHistStyle(pads[i].h);
      pads[i].h->SetMinimum(0.0);
      pads[i].h->Draw("COL");
    } else {
      TH2D frame("frame", "", 10, -ptMax, ptMax, 10, -ptMax, ptMax);
      SetHistStyle(&frame);
      frame.Draw("AXIS");
    }
    DrawPolarGrid(ptMax);

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextAlign(22);
    lat.SetTextSize(0.048);
    lat.DrawLatex(0.50, 0.965, Form("p_{T} vs #phi_{pq} (%s)", pads[i].title));
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.50, 0.035, Form("coin-window entries: %lld", pads[i].n));
  }

  c.cd();
  TLatex title;
  title.SetNDC(true);
  title.SetTextSize(0.045);
  title.SetTextAlign(13);
  title.DrawLatex(0.075, 0.965, ZPretty(z).c_str());

  c.SaveAs(pngPath.c_str());
}

} // namespace

void pt_vs_phipq(const char* groupsRootIn = "groups",
                 const char* dataRootIn = "Pass0_SkimmedDataROOTfiles",
                 const char* resultsRootIn = "results")
{
  const std::string groupsRoot = NormalizeSlashes(groupsRootIn ? groupsRootIn : "groups");
  const std::string dataRoot = NormalizeSlashes(dataRootIn ? dataRootIn : "Pass0_SkimmedDataROOTfiles");
  const std::string resultsRoot = NormalizeSlashes(resultsRootIn ? resultsRootIn : "results");
  const std::string pngDir = resultsRoot + "/PNGs";
  const std::string tableDir = resultsRoot + "/tables";
  MkdirP(pngDir);
  MkdirP(tableDir);

  std::vector<std::string> groupFiles;
  CollectFilesRecursive(groupsRoot, ".list", groupFiles);
  std::sort(groupFiles.begin(), groupFiles.end());
  if (groupFiles.empty()) {
    std::cerr << "ERROR: no group .list files found under " << groupsRoot << "\n";
    return;
  }

  std::map<std::string, Bucket> buckets;
  for (const auto& gf : groupFiles) ReadGroupList(gf, buckets);

  const std::vector<std::string> zOrder = {"z0p36", "z0p5", "z0p67", "z0p9"};
  const std::vector<std::string> passOrder = {"pass4", "pass5"};
  const double ptMax = 0.6;
  const int nXY = 180;
  const int nPt = 60;
  const int nPhi = 72;
  const TString phiExpr =
    "((P_kin_secondary_ph_xq>=0)*P_kin_secondary_ph_xq + "
    "(P_kin_secondary_ph_xq<0)*(P_kin_secondary_ph_xq + 2*TMath::Pi()))";
  const TString dataCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_cal_etottracknorm<0.8)";
  const TString ptCut = Form("(pt>0) && (pt<%.6f)", ptMax);
  const TString baseCuts = CombineCutsAND(dataCuts, ptCut);

  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  std::ofstream csv(tableDir + "/pt_vs_phipq.csv");
  csv << "z_label,pass,pt_bin,phi_bin,pt_lo,pt_hi,pt_center,"
      << "phi_lo_rad,phi_hi_rad,phi_center_rad,"
      << "random_subtracted_counts,coin_counts_raw,random_avg_counts,"
      << "n_runs,n_manifests,n_groups\n";

  for (const auto& z : zOrder) {
    std::unique_ptr<TH2D> hxy[2];
    std::unique_ptr<TH2D> hpf[2];
    std::unique_ptr<TH2D> hpfCoin[2];
    std::unique_ptr<TH2D> hpfRand[2];
    Long64_t nCoinTotal[2] = {0, 0};

    for (size_t ip = 0; ip < passOrder.size(); ++ip) {
      const std::string pass = passOrder[ip];
      const std::string key = pass + "|" + z;
      auto bit = buckets.find(key);
      hxy[ip].reset(new TH2D(Form("hxy_%s_%s", pass.c_str(), z.c_str()),
                             "", nXY, -ptMax, ptMax, nXY, -ptMax, ptMax));
      hpf[ip].reset(new TH2D(Form("hpf_%s_%s", pass.c_str(), z.c_str()),
                             "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax));
      hpfCoin[ip].reset(new TH2D(Form("hpf_coin_%s_%s", pass.c_str(), z.c_str()),
                                 "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax));
      hpfRand[ip].reset(new TH2D(Form("hpf_rand_%s_%s", pass.c_str(), z.c_str()),
                                 "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax));

      if (bit == buckets.end()) {
        std::cerr << "WARNING: no bucket for " << pass << " " << z << "\n";
        continue;
      }

      for (int run : bit->second.runs) {
        const std::string rootPath = DataRootPath(dataRoot, run);
        if (gSystem->AccessPathName(rootPath.c_str(), kReadPermission)) {
          std::cerr << "WARNING: missing data file for run " << run << ": " << rootPath << "\n";
          continue;
        }
        TFile f(rootPath.c_str(), "READ");
        if (f.IsZombie()) {
          std::cerr << "WARNING: cannot open " << rootPath << "\n";
          continue;
        }
        TTree* T = dynamic_cast<TTree*>(f.Get("T"));
        if (!T) {
          std::cerr << "WARNING: missing tree T in " << rootPath << "\n";
          continue;
        }

        CoincidenceResult win = ComputeCoincidenceRandomSubtraction(T, baseCuts, cfg);
        if (!(win.CoinWindowNs.second > win.CoinWindowNs.first)) {
          std::cerr << "WARNING: invalid coincidence window for run " << run << "\n";
          continue;
        }
        TString xyExpr = Form("pt*sin(%s):pt*cos(%s)", phiExpr.Data(), phiExpr.Data());
        TString pfExpr = Form("pt:%s", phiExpr.Data());
        TDirectory* oldDir = gDirectory;
        TDirectory* histDir = hxy[ip]->GetDirectory();
        if (histDir) histDir->cd();

        TString wideGate = BuildRangeCut(cfg.CtBranchName, cfg.WideWindowMinNs, cfg.WideWindowMaxNs);
        TString coinGate = BuildRangeCut(cfg.CtBranchName, win.CoinWindowNs.first, win.CoinWindowNs.second);
        TString cutsCoin = CombineCutsAND(baseCuts, CombineCutsAND(wideGate, coinGate));
        nCoinTotal[ip] += T->GetEntries(cutsCoin);

        TString srun = Form("%s_%s_%d", pass.c_str(), z.c_str(), run);
        TH2D xyCoin(Form("xy_coin_%s", srun.Data()), "", nXY, -ptMax, ptMax, nXY, -ptMax, ptMax);
        TH2D xyRandSum(Form("xy_rand_sum_%s", srun.Data()), "", nXY, -ptMax, ptMax, nXY, -ptMax, ptMax);
        TH2D pfCoin(Form("pf_coin_%s", srun.Data()), "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax);
        TH2D pfRandSum(Form("pf_rand_sum_%s", srun.Data()), "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax);
        xyCoin.SetDirectory(histDir);
        xyRandSum.SetDirectory(histDir);
        pfCoin.SetDirectory(histDir);
        pfRandSum.SetDirectory(histDir);
        xyCoin.Sumw2();
        xyRandSum.Sumw2();
        pfCoin.Sumw2();
        pfRandSum.Sumw2();

        T->Draw(xyExpr + ">>" + xyCoin.GetName(), cutsCoin, "goff");
        T->Draw(pfExpr + ">>" + pfCoin.GetName(), cutsCoin, "goff");

        int nRand = 0;
        for (const auto& rwin : win.RandomWindowListNs) {
          TString randGate = BuildRangeCut(cfg.CtBranchName, rwin.first, rwin.second);
          TString cutsRand = CombineCutsAND(baseCuts, CombineCutsAND(wideGate, randGate));
          TH2D xyTmp(Form("xy_rand_tmp_%s_%d", srun.Data(), nRand),
                     "", nXY, -ptMax, ptMax, nXY, -ptMax, ptMax);
          TH2D pfTmp(Form("pf_rand_tmp_%s_%d", srun.Data(), nRand),
                     "", nPhi, 0.0, 2.0 * TMath::Pi(), nPt, 0.0, ptMax);
          xyTmp.SetDirectory(histDir);
          pfTmp.SetDirectory(histDir);
          xyTmp.Sumw2();
          pfTmp.Sumw2();
          T->Draw(xyExpr + ">>" + xyTmp.GetName(), cutsRand, "goff");
          T->Draw(pfExpr + ">>" + pfTmp.GetName(), cutsRand, "goff");
          xyRandSum.Add(&xyTmp);
          pfRandSum.Add(&pfTmp);
          ++nRand;
        }

        if (nRand > 0) {
          xyRandSum.Scale(1.0 / nRand);
          pfRandSum.Scale(1.0 / nRand);
        }

        hxy[ip]->Add(&xyCoin);
        hxy[ip]->Add(&xyRandSum, -1.0);
        hpf[ip]->Add(&pfCoin);
        hpf[ip]->Add(&pfRandSum, -1.0);
        hpfCoin[ip]->Add(&pfCoin);
        hpfRand[ip]->Add(&pfRandSum);

        if (oldDir) oldDir->cd();
      }

      hxy[ip]->SetDirectory(nullptr);
      hpf[ip]->SetDirectory(nullptr);
      hpfCoin[ip]->SetDirectory(nullptr);
      hpfRand[ip]->SetDirectory(nullptr);

      for (int ix = 1; ix <= hpf[ip]->GetNbinsX(); ++ix) {
        const double phiLo = hpf[ip]->GetXaxis()->GetBinLowEdge(ix);
        const double phiHi = hpf[ip]->GetXaxis()->GetBinUpEdge(ix);
        const double phiCtr = hpf[ip]->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hpf[ip]->GetNbinsY(); ++iy) {
          const double ptLo = hpf[ip]->GetYaxis()->GetBinLowEdge(iy);
          const double ptHi = hpf[ip]->GetYaxis()->GetBinUpEdge(iy);
          const double ptCtr = hpf[ip]->GetYaxis()->GetBinCenter(iy);
          const double counts = hpf[ip]->GetBinContent(ix, iy);
          const double coinCounts = hpfCoin[ip]->GetBinContent(ix, iy);
          const double randCounts = hpfRand[ip]->GetBinContent(ix, iy);
          csv << z << "," << pass << ","
              << (iy - 1) << "," << (ix - 1) << ","
              << std::setprecision(8) << ptLo << "," << ptHi << "," << ptCtr << ","
              << phiLo << "," << phiHi << "," << phiCtr << ","
              << counts << "," << coinCounts << "," << randCounts << ","
              << (bit == buckets.end() ? 0 : bit->second.runs.size()) << ","
              << (bit == buckets.end() ? 0 : bit->second.manifests.size()) << ","
              << (bit == buckets.end() ? 0 : bit->second.groups.size()) << "\n";
        }
      }
    }

    const std::string pngPath = pngDir + "/pt_vs_phipq_" + z + ".png";
    SavePolarPlot(pngPath, z, hxy[0].get(), hxy[1].get(),
                  nCoinTotal[0], nCoinTotal[1], ptMax);
  }

  std::cout << "Wrote " << tableDir << "/pt_vs_phipq.csv\n";
  std::cout << "Wrote " << pngDir << "/pt_vs_phipq_<z>.png\n";
}
