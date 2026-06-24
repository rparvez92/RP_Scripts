// PlotCoinXsec.C (rsidis_xs_v4)
// v4 goal: xsec vs phipq panels (pads = pT bins, curves = z bins)
// Applies v3-style subtractions on DATA:
//   - Random subtraction via CoincidenceRandomSubtraction.h
//   - Category subtraction: (signal - positron) - (dummy - positron_dummy)
// SIMC:
//   - Uses h10 tree
//   - Uses pt = sqrt(pt2) for SIM (since pt branch not present)
//   - Sim yield is computed by summing (Weight*Norm) over selected rows (no histogram lookup!)
// Notes:
//   - No beta cuts in DATA (per Radwan request)
//   - phipq in DATA is P_kin_secondary_ph_xq wrapped to [0,2pi)
//   - phipq in SIM uses phipq wrapped to [0,2pi)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TF1.h"

#include "../include/CoincidenceRandomSubtraction.h"

namespace {
// Color palette
static const Color_t zColors[] = { kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7 };
static const int nZColors = sizeof(zColors)/sizeof(zColors[0]);


static std::string NormalizeSlashes(std::string s) {
  std::string out;
  bool prev = false;
  for (char c : s) {
    if (c == '/') {
      if (!prev) out.push_back(c);
      prev = true;
    } else {
      out.push_back(c);
      prev = false;
    }
  }
  if (out.size() > 1 && out.back() == '/') out.pop_back();
  return out;
}

static std::string Dirname(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return p.substr(0, pos);
}

static std::string Basename(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return p;
  return p.substr(pos+1);
}

// settingsRoot can be relative ("settings") or absolute ("/.../settings")
static std::string MakeSettingIdFromManifestPath(const std::string& manifestPath,
                                                const std::string& settingsRoot) {
  std::string mdir  = NormalizeSlashes(Dirname(manifestPath));
  std::string sroot = NormalizeSlashes(settingsRoot);

  if (mdir.rfind(sroot + "/", 0) == 0) return mdir.substr(sroot.size() + 1);

  std::string sbase = Basename(sroot); // usually "settings"
  auto pos = mdir.find("/" + sbase + "/");
  if (pos != std::string::npos) return mdir.substr(pos + sbase.size() + 2);
  if (mdir.rfind(sbase + "/", 0) == 0) return mdir.substr(sbase.size() + 1);

  return Basename(mdir);
}

static void MkdirP(const std::string& path) {
  if (path.empty()) return;
  gSystem->mkdir(path.c_str(), true);
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

static std::vector<std::string> SplitCSVLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i=0;i<line.size();i++) {
    char c = line[i];
    if (c == '"') { inq = !inq; continue; }
    if (c == ',' && !inq) { out.push_back(cur); cur.clear(); continue; }
    cur.push_back(c);
  }
  out.push_back(cur);
  for (auto& s : out) {
    auto issp = [](unsigned char c){ return std::isspace(c)!=0; };
    while (!s.empty() && issp((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && issp((unsigned char)s.back())) s.pop_back();
  }
  return out;
}

struct RunNormInfo {
  double charge_mC = 0.0;   // BCM4A_Q
  double hms_eff   = 1.0;   // h_esing_Eff
  double ps_factor = 1.0;   // ps5
  double comp_lt   = 1.0;   // optional
};

static bool LoadRunMetadataCSV(const std::string& csvPath,
                               std::unordered_map<int, RunNormInfo>& out)
{
  std::ifstream f(csvPath);
  if (!f.is_open()) return false;

  std::string header;
  if (!std::getline(f, header)) return false;
  auto cols = SplitCSVLine(header);

  auto idx = [&](const std::string& name)->int {
    for (size_t i=0;i<cols.size();i++) if (cols[i]==name) return (int)i;
    return -1;
  };

  int i_run   = idx("run");
  int i_q     = idx("BCM4A_Q");
  int i_eff   = idx("h_esing_Eff");
  int i_ps    = idx("ps5");
  int i_lt    = idx("comp_livetime");

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto v = SplitCSVLine(line);
    if (i_run < 0 || i_run >= (int)v.size()) continue;

    int run = std::atoi(v[i_run].c_str());
    RunNormInfo info;

    if (i_q >= 0 && i_q < (int)v.size()) info.charge_mC = std::atof(v[i_q].c_str());
    if (i_eff >= 0 && i_eff < (int)v.size()) info.hms_eff = std::atof(v[i_eff].c_str());
    if (i_ps >= 0 && i_ps < (int)v.size()) info.ps_factor = std::atof(v[i_ps].c_str());
    if (i_lt >= 0 && i_lt < (int)v.size()) info.comp_lt = std::atof(v[i_lt].c_str());

    if (!std::isfinite(info.hms_eff) || info.hms_eff <= 0) info.hms_eff = 1.0;
    if (!std::isfinite(info.ps_factor) || info.ps_factor <= 0) info.ps_factor = 1.0;
    if (!std::isfinite(info.comp_lt) || info.comp_lt <= 0) info.comp_lt = 1.0;

    out[run] = info;
  }
  return true;
}

static double ReadSigmaModel(const std::string& leafDir, std::ofstream& log) {
  std::ifstream f(leafDir + "/sigma_model.txt");
  if (!f.is_open()) {
    log << "WARNING: cannot open sigma_model.txt; using 1.\n";
    return 1.0;
  }
  double s=1.0; f>>s;
  if (!std::isfinite(s) || s<=0) {
    log << "WARNING: bad sigma_model value; using 1.\n";
    return 1.0;
  }
  log << "sigma_model=" << s << "\n";
  return s;
}

static std::string DataRootPath(int run) {
  char buf[512];
  std::snprintf(buf, sizeof(buf), "./Skimmed_ROOTfiles/skimmed_coin_replay_production_%d_-1.root", run);
  return std::string(buf);
}

struct Yield {
  double sumw  = 0.0;
  double sumw2 = 0.0;
};

static Yield AddYieldFromRun(TTree* T,
                             const TString& cutsExpr,
                             const TString& weightExpr,
                             const CoincidenceConfig& cfg,
                             const CoincidenceResult& win)
{
  // Unique histogram name but **not** registered in any directory.
  static long long uid = 0;
  TString hname = Form("h_yield_%lld", uid++);

  TH1D h(hname, hname, 1, 0, 1);
  h.SetDirectory(nullptr);
  h.Sumw2(true);


  // --- DEBUG: does the selection itself work?
  Long64_t nPlain = T->GetEntries(cutsExpr);

  // Build the same boolean cuts used inside FillRandomSubtractedHistogramWithWindows
  TString wideGate = BuildRangeCut(cfg.CtBranchName, cfg.WideWindowMinNs, cfg.WideWindowMaxNs);
  TString cutsWide = CombineCutsAND(cutsExpr, wideGate);

  Long64_t nWide = T->GetEntries(cutsWide);

  // Coin window count (using the cached window)
  TString coinGate = BuildRangeCut(cfg.CtBranchName, win.CoinWindowNs.first, win.CoinWindowNs.second);
  TString cutsCoin = CombineCutsAND(cutsWide, coinGate);

  Long64_t nCoin = T->GetEntries(cutsCoin);

/*
  std::cerr << "DBG Ct: PeakCenter=" << win.PeakCenterNs
            << " Coin=[" << win.CoinWindowNs.first << "," << win.CoinWindowNs.second << "]"
            << " nRand=" << win.RandomWindowListNs.size()
            << " nWide=" << nWide
            << " nCoin=" << nCoin
            << "\n";
*/

  // Constant expression; only selection + weight matter
  Bool_t oldAdd = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kTRUE);
  FillRandomSubtractedHistogramWithWindows(T, cutsExpr, weightExpr, "0.5", &h, cfg, win);
  TH1::AddDirectory(oldAdd);


  Yield y;
  y.sumw  = h.GetBinContent(1);
  double e = h.GetBinError(1);
  y.sumw2 = e*e;
  return y;
}

// SIM yield: sum Weight*norm over selected rows (no histogram creation)
static Yield SimYieldSum(TChain& sim,
                         const TString& exprWeight,
                         const TString& cutBool)
{
  Yield out;
  Long64_t n = sim.Draw(exprWeight, cutBool, "goff");
  if (n <= 0) return out;

  const double* v = sim.GetV1();
  if (!v) return out;

  double sum = 0.0, sum2 = 0.0;
  for (Long64_t i=0; i<n; i++) {
    double w = v[i];
    if (!std::isfinite(w)) continue;
    sum  += w;
    sum2 += w*w;
  }
  out.sumw  = sum;
  out.sumw2 = sum2; // MC stat: var = sum(w^2)
  return out;
}

static double Center(double lo, double hi) { return 0.5*(lo+hi); }

static std::vector<double> PhiEdgesN(int n) {
  const double twopi = 2.0*M_PI;
  std::vector<double> e;
  for (int i=0;i<=n;i++) e.push_back(twopi*i/n);
  return e;
}


struct Binning {
  std::vector<double> ptEdges {0.0, 0.10, 0.20, 0.30, 0.40};
  std::vector<double> zEdges  {0.30, 0.40, 0.50, 0.60, 0.70};
  std::vector<double> phiEdges = PhiEdgesN(8);
};

/*
struct Binning {
  std::vector<double> ptEdges {0.0, 0.40};   // one wide bin
  std::vector<double> zEdges  {0.30, 0.70};  // one wide bin
  std::vector<double> phiEdges = PhiEdgesN(8);
};
*/

} // namespace

void PlotCoinXsec(const char* manifestPath,
                  const char* resultsRoot = "results",
                  const char* settingsRoot = "settings",
                  bool savePNGs = true,
                  bool saveTables = true)
{
  // Keep fits alive until canvas is written
  std::vector<TF1*> fitFuncs;
  fitFuncs.reserve(1000);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);
  TH1::AddDirectory(kFALSE); // keep ROOT from owning lots of histos

  if (!manifestPath || std::string(manifestPath).empty()) {
    std::cerr << "ERROR: manifestPath is required\n";
    return;
  }

  const std::string manifest = manifestPath;
  const std::string leafDir  = NormalizeSlashes(Dirname(manifest));
  const std::string settingId = MakeSettingIdFromManifestPath(manifest, settingsRoot);

  const std::string resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + settingId);
  MkdirP(resultsDir);
  MkdirP(resultsDir + "/PNGs");
  MkdirP(resultsDir + "/tables");

  std::ofstream log(resultsDir + "/log.txt", std::ios::app);
  log << "=== PlotCoinXsec v4 (phipq panels) ===\n";
  log << "manifest: " << manifest << "\n";
  log << "leafDir: " << leafDir << "\n";
  log << "setting_id: " << settingId << "\n";
  log << "resultsDir: " << resultsDir << "\n";

  auto runs_sig    = ReadRunList(leafDir + "/runs_signal.txt");
  auto runs_dum    = ReadRunList(leafDir + "/runs_dummy.txt");
  auto runs_posSig = ReadRunList(leafDir + "/runs_positron.txt");
  auto runs_posDum = ReadRunList(leafDir + "/runs_positron_dummy.txt");

  log << "runs_signal: " << runs_sig.size() << "\n";
  log << "runs_dummy: " << runs_dum.size() << "\n";
  log << "runs_positron: " << runs_posSig.size() << "\n";
  log << "runs_positron_dummy: " << runs_posDum.size() << "\n";

  std::unordered_map<int, RunNormInfo> runInfo;
  if (!LoadRunMetadataCSV(leafDir + "/run_metadata.csv", runInfo)) {
    log << "ERROR: cannot load run_metadata.csv\n";
    std::cerr << "ERROR: cannot load run_metadata.csv\n";
    return;
  }

  const double sigmaModel = ReadSigmaModel(leafDir, log);

  // DATA base cuts (underscore), NO beta cuts
  //const TCut dataCuts =
    //"(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    //"(H_cal_etottracknorm>0.7) && "
    //"(H_cer_npeSum>2.0) && "
    //"(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    //"(P_cal_etottracknorm<0.8) && "
    //"(((P_gtr_p<2.7) && (P_aero_npeSum>2)) || "
    //" ((P_gtr_p>=2.7) && (P_hgcer_npeSum>1) && (P_aero_npeSum>2)))";

  const TCut dataCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_gtr_th>-0.060) && (H_gtr_th<0.060) && "
    "(H_gtr_ph>-0.022) && (H_gtr_ph<0.022) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_gtr_th>-0.045) && (P_gtr_th<0.045) && "
    "(P_gtr_ph>-0.024) && (P_gtr_ph<0.024) && "
    "(P_cal_etottracknorm<0.8)";

  // DATA vars (skimmed)
  const std::string data_pt  = "pt";
  const std::string data_z   = "z";
  //const std::string data_phi = "(P_kin_secondary_ph_xq<0 ? P_kin_secondary_ph_xq+2*TMath::Pi() : P_kin_secondary_ph_xq)";
  // Wrap phipq into [0,2pi) WITHOUT ternary (ROOT-safe)
  const std::string data_phi = "((P_kin_secondary_ph_xq>=0)*P_kin_secondary_ph_xq + (P_kin_secondary_ph_xq<0)*(P_kin_secondary_ph_xq + 2*TMath::Pi()))";

  const double wall_thickness_ratio = 3.82;
  const double dumScale = 1.0 / wall_thickness_ratio;


  // SIM vars (SIMC h10)
  const std::string sim_pt   = "sqrt(pt2)";
  const std::string sim_z    = "z";
  //const std::string sim_phi  = "phipq";
  // Wrap SIM phipq too (in case SIMC stores it in [-pi,pi])
  const std::string sim_phi  = "((phipq>=0)*phipq + (phipq<0)*(phipq + 2*TMath::Pi()))";

  // SIM chain
  TChain sim("h10");
  {
    std::string simGlob = "./simc_worksim/" + settingId + "/*.root";
    int nAdded = sim.Add(simGlob.c_str());
    log << "SIM glob: " << simGlob << " added=" << nAdded << "\n";
    if (sim.GetEntries() <= 0) {
      std::cerr << "ERROR: SIM chain empty for: " << simGlob << "\n";
      log << "ERROR: SIM chain empty\n";
      return;
    }
    std::cout << "SIM phipq min=" << sim.GetMinimum("phipq")
              << " max=" << sim.GetMaximum("phipq") << "\n";
  }

  const Long64_t nGenSim = sim.GetEntries();
  std::cout << "SIM: total entries in h10 = " << nGenSim << std::endl;
  if (nGenSim <= 0) {
    std::cerr << "ERROR: SIM chain has zero entries. Check SIM file path/tree name.\n";
  }

  // SIM cuts
  const TCut simCutsBase = 
    "(hsdelta>-8.0) && (hsdelta<8.0) && "
    "(hsxptar>-0.060) && (hsxptar<0.060) && "
    "(hsyptar>-0.022) && (hsyptar<0.022) && "
    "(ssdelta>-10.0) && (ssdelta<22.0) && "
    "(ssxptar>-0.045) && (ssxptar<0.045) && "
    "(ssyptar>-0.024) && (ssyptar<0.024)";

  const double simNormalizationFactor = 0.842205e11;
  const double simNtried = 53411876.0;                // Ntried from .hist
  //const double simTryScale = 1000.0 / simNtried;      // matches "MeV:" convention in .hist
  const double simTryScale = (nGenSim > 0) ? (1.0 / double(nGenSim)) : 0.0;
  const TString simWeightExpr = Form("Weight*%e*%e", simNormalizationFactor, simTryScale);


  std::cout << "SIM: 1000/Ntried = " << (1000.0/simNtried)
            << " , 1/Nentries = " << (nGenSim>0 ? 1.0/double(nGenSim) : 0.0)
            << " , ratio (1/Nentries)/(1000/Ntried) = "
            << ((nGenSim>0) ? ( (1.0/double(nGenSim)) / (1000.0/simNtried) ) : 0.0)
            << std::endl;



  // Random subtraction config (skimmed branch name)
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  // Binning
  Binning binning;
  const int nPt  = (int)binning.ptEdges.size()  - 1;
  const int nZ   = (int)binning.zEdges.size()   - 1;
  const int nPhi = (int)binning.phiEdges.size() - 1;

  // CSV
  std::ofstream csv;
  if (saveTables) {
    std::string csvPath = resultsDir + "/tables/xsec_phipq_pt_z.csv";
    csv.open(csvPath);
    csv
      << "setting_id,pt_bin,z_bin,phibin,"
      << "pt_lo,pt_hi,pt_center,z_lo,z_hi,z_center,phi_lo,phi_hi,phi_center,"
      << "yield_sig,yield_pos,yield_dum,yield_posdum,yield_net,yield_net_err,"
      << "yield_sim,yield_sim_err,"
      << "xsec,xsec_stat\n";
    log << "CSV: " << csvPath << "\n";
  }

  // Coin window cache per run
  std::unordered_map<int, CoincidenceResult> coinCache;

  auto CategoryYield = [&](const std::vector<int>& runs,
                           const TCut& baseCutsBool,
                           const TCut& kinCutsBool,
                           double& totalCharge_mC) -> Yield
  {
    Yield out;
    totalCharge_mC = 0.0;

    int nOpened = 0;
    bool didSanity = false;

    const TString baseCutsStr = TString((baseCutsBool).GetTitle());
    const TString fullCutsStr = TString((baseCutsBool && kinCutsBool).GetTitle());

    for (int run : runs) {
      auto it = runInfo.find(run);
      if (it == runInfo.end()) {
        log << "WARNING: missing metadata for run " << run << "\n";
        continue;
      }
      const RunNormInfo& info = it->second;
      totalCharge_mC += info.charge_mC;

      double w = 1.0;
      w *= info.ps_factor;
      if (info.hms_eff > 0) w /= info.hms_eff;
      if (info.comp_lt > 0) w /= info.comp_lt;

      std::string path = DataRootPath(run);
      TFile* f = TFile::Open(path.c_str(), "READ");
      if (!f || f->IsZombie()) {
        log << "WARNING: cannot open " << path << "\n";
        std::cerr << "WARNING: cannot open " << path << "\n";
        if (f) f->Close();
        continue;
      }
      TTree* T = (TTree*) f->Get("T");
      if (!T) {
        log << "WARNING: missing tree T in " << path << "\n";
        std::cerr << "WARNING: missing tree T in " << path << "\n";
        f->Close();
        continue;
      }

      ++nOpened;
      if (!didSanity) {
        Long64_t nAll  = T->GetEntries();
        Long64_t nBase = T->GetEntries(baseCutsStr);
        Long64_t nFull = T->GetEntries(fullCutsStr);
        log << "SANITY (run " << run << "): entries(all)=" << nAll
            << " entries(baseCuts)=" << nBase
            << " entries(fullCuts)=" << nFull << "\n";
        //std::cerr << "SANITY (run " << run << "): entries(all)=" << nAll
                  //<< " entries(baseCuts)=" << nBase
                  //<< " entries(fullCuts)=" << nFull << "\n";
        didSanity = true;
      }

      auto itC = coinCache.find(run);
      if (itC == coinCache.end()) {
        Bool_t oldAdd = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kTRUE);
        CoincidenceResult win = ComputeCoincidenceRandomSubtraction(T, baseCutsStr, cfg);
        TH1::AddDirectory(oldAdd);

      //std::cout << "DBG Win: CoinYield=" << win.CoinYield
                //<< " RandMean=" << win.RandomMeanYield
                //<< " Sub=" << win.RandomSubtractedYield
                //<< "  (err " << win.RandomSubtractedYieldErr << ")\n";


        coinCache[run] = win;
        itC = coinCache.find(run);
      }

      TString wexpr = Form("%g", w);
      Yield y = AddYieldFromRun(T, fullCutsStr, wexpr, cfg, itC->second);
      out.sumw  += y.sumw;
      out.sumw2 += y.sumw2;

      f->Close();
    }

    if (nOpened == 0) {
      log << "ERROR: none of the data files for this category could be opened."
          << " Check your working directory (relative paths), symlinks, and file naming.\n";
      std::cerr << "ERROR: none of the data files for this category could be opened."
                << " Check working directory (relative paths), symlinks, and file naming.\n";
    }

    if (totalCharge_mC > 0) {
      out.sumw  /= totalCharge_mC;
      out.sumw2 /= (totalCharge_mC * totalCharge_mC);
    }
    if(nOpened == 0) {
      log << "ERROR: opened 0 ROOT files for this category (all yields will be zero).\n";
      std::cerr << "ERROR: opened 0 ROOT files for this category (all yields will be zero).\n";
    }

    return out;
  };

  // Canvas
  int nCols = (nPt <= 2) ? nPt : 2;
  int nRows = (int)std::ceil((double)nPt / nCols);
  TCanvas* c = new TCanvas("c_phipq_panels", "xsec vs phipq panels", 1200, 450*nRows);
  c->Divide(nCols, nRows, 0.001, 0.001);

  std::vector<double> yMaxByPt(nPt, 0.0);

  struct GraphPack { int ptBin; int zBin; TGraphErrors* g; };
  std::vector<GraphPack> graphs;

  for (int iPt=0; iPt<nPt; iPt++) {
    double ptLo = binning.ptEdges[iPt];
    double ptHi = binning.ptEdges[iPt+1];
    double ptCtr= Center(ptLo, ptHi);

    for (int iZ=0; iZ<nZ; iZ++) {
      double zLo = binning.zEdges[iZ];
      double zHi = binning.zEdges[iZ+1];
      double zCtr= Center(zLo, zHi);

      //std::vector<double> x(nPhi), y(nPhi), ex(nPhi), ey(nPhi);
      std::vector<double> x, y, ex, ey;
      x.reserve(nPhi);
      y.reserve(nPhi);
      ex.reserve(nPhi);
      ey.reserve(nPhi);

      for (int iPhi=0; iPhi<nPhi; iPhi++) {
        double phiLo = binning.phiEdges[iPhi];
        double phiHi = binning.phiEdges[iPhi+1];
        double phiCtr= Center(phiLo, phiHi);

        TCut kinData = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                            data_pt.c_str(), ptLo, data_pt.c_str(), ptHi,
                            data_z.c_str(),  zLo,  data_z.c_str(),  zHi,
                            data_phi.c_str(),phiLo, data_phi.c_str(),phiHi);

        static bool printedOnce = false;
        if (!printedOnce) {
          // check only for the signal category’s first run (or whichever you like)
          std::cout << "DEBUG kinData = " << kinData.GetTitle() << "\n";
          printedOnce = true;
        }

        double qSig=0, qPos=0, qDum=0, qPosDum=0;
        Yield ySig    = CategoryYield(runs_sig,    dataCuts, kinData, qSig);
        Yield yPos    = CategoryYield(runs_posSig, dataCuts, kinData, qPos);
        Yield yDum    = CategoryYield(runs_dum,    dataCuts, kinData, qDum);
        Yield yPosDum = CategoryYield(runs_posDum, dataCuts, kinData, qPosDum);

        //double yNet  = (ySig.sumw - yPos.sumw) - (yDum.sumw - yPosDum.sumw);
        //double eNet2 = ySig.sumw2 + yPos.sumw2 + yDum.sumw2 + yPosDum.sumw2;
        //double eNet  = (eNet2 > 0) ? std::sqrt(eNet2) : 0.0;

        // (Data − PosData) − (Dummy − PosDummy)/wall_thickness_ratio
        double yNet  = (ySig.sumw - yPos.sumw) - dumScale * (yDum.sumw - yPosDum.sumw);
        // Error propagation: dummy terms get scaled by dumScale
        double eNet2 = ySig.sumw2 + yPos.sumw2
                     + dumScale*dumScale * (yDum.sumw2 + yPosDum.sumw2);
        double eNet  = (eNet2 > 0) ? std::sqrt(eNet2) : 0.0;



        // SIM selection for this bin (boolean only)
        TCut kinSim = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                           sim_pt.c_str(), ptLo, sim_pt.c_str(), ptHi,
                           sim_z.c_str(),  zLo,  sim_z.c_str(),  zHi,
                           sim_phi.c_str(),phiLo, sim_phi.c_str(),phiHi);

        TString simBool = TString((simCutsBase && kinSim).GetTitle());
        Yield ySimY = SimYieldSum(sim, simWeightExpr, simBool);
        double ySim = ySimY.sumw;
        double eSim = (ySimY.sumw2 > 0) ? std::sqrt(ySimY.sumw2) : 0.0;


        bool mc_ok = (ySim > 0) && std::isfinite(ySim) && std::isfinite(eSim);
        double relSim = (mc_ok && ySim>0) ? (eSim / ySim) : 1e9;

        // Pick a threshold (start with 0.5, tighten later)
        bool mc_starved = (!mc_ok) || (relSim > 0.5);

        auto yerr = [](const Yield& Y){ return (Y.sumw2>0) ? std::sqrt(Y.sumw2) : 0.0; };

        //std::cout
          //<< "  ySig=" << ySig.sumw << " +/- " << yerr(ySig)
          //<< "  yPos=" << yPos.sumw << " +/- " << yerr(yPos)
          //<< "  yDum=" << yDum.sumw << " +/- " << yerr(yDum)
          //<< "  yPosDum=" << yPosDum.sumw << " +/- " << yerr(yPosDum)
          //<< "  yNet=" << yNet << " +/- " << eNet
          //<< "  ySim=" << ySim << " +/- " << eSim
          //<< "  sigmaModel=" << sigmaModel
          //<< "\n";

        // Cross section
        double xsec = 0.0, xsecErr = 0.0;

        // --- Compute xsec/xsecErr (keep your logic) ---
        if (!mc_starved) {
          xsec = (yNet / ySim) * sigmaModel;

          double rel2 = 0.0;
          if (std::abs(yNet) > 0) rel2 += (eNet / yNet) * (eNet / yNet);
          rel2 += (eSim / ySim) * (eSim / ySim);
          xsecErr = std::abs(xsec) * std::sqrt(rel2);
        } else {
          // For starved bins, we will NOT add a point to the graph
          xsec    = 0.0;
          xsecErr = 0.0;
        }

        // --- IMPORTANT CHANGE: only keep non-starved points for plotting/fitting ---
        if (!mc_starved) {
          x.push_back(phiCtr);
          ex.push_back(0.5 * (phiHi - phiLo));
          y.push_back(xsec);
          ey.push_back(xsecErr);

          // Autoscale should only consider the points we actually keep
          if (std::isfinite(xsec) && std::isfinite(xsecErr) && (xsecErr >= 0)) {
            yMaxByPt[iPt] = std::max(yMaxByPt[iPt], xsec + xsecErr);
          }
        }

        // (CSV writing can remain unchanged below; it will still record all bins,
        // including starved bins, which is usually useful for debugging.)

        if (saveTables && csv.is_open()) {
          csv
            << settingId << ","
            << iPt << "," << iZ << "," << iPhi << ","
            << ptLo << "," << ptHi << "," << ptCtr << ","
            << zLo  << "," << zHi  << "," << zCtr  << ","
            << phiLo << "," << phiHi << "," << phiCtr << ","
            << ySig.sumw << "," << yPos.sumw << "," << yDum.sumw << "," << yPosDum.sumw << ","
            << yNet << "," << eNet << ","
            << ySim << "," << eSim << ","
            << xsec << "," << xsecErr
            << "\n";
        }
      }

      //auto g = new TGraphErrors(nPhi, x.data(), y.data(), ex.data(), ey.data());
      auto g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
      g->SetName(Form("g_pt%d_z%d", iPt, iZ));
      g->SetMarkerStyle(20 + (iZ % 10));
      g->SetMarkerSize(0.9);
      g->SetLineWidth(2);
      graphs.push_back({iPt, iZ, g});
    }
  }

  // Draw panels
  for (int iPt=0; iPt<nPt; iPt++) {
    c->cd(iPt+1);
    gPad->SetTicks(1,1);

    double yMax = (yMaxByPt[iPt] > 0) ? yMaxByPt[iPt] * 1.25 : 1.0;
    TH1F* frame = (TH1F*) gPad->DrawFrame(0.0, 0.0, 2.0*TMath::Pi(), yMax);

    double ptLo = binning.ptEdges[iPt];
    double ptHi = binning.ptEdges[iPt+1];

    frame->GetXaxis()->SetTitle("#phi_{pq} (rad)");
    frame->GetYaxis()->SetTitle("d#sigma");
    frame->GetYaxis()->SetTitleOffset(1.2);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22);     // center horizontally and vertically around the point
    latex.DrawLatex(0.50, 0.92, Form("p_{T} #in [%.2f, %.2f] GeV", ptLo, ptHi));

    TLegend* leg = new TLegend(0.55, 0.62, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    int drawn = 0;
    for (const auto& gp : graphs) {
      if (gp.ptBin != iPt) continue;

      double zLo = binning.zEdges[gp.zBin];
      double zHi = binning.zEdges[gp.zBin + 1];

      // --- style per z ---
      Color_t col = zColors[gp.zBin % nZColors];
      gp.g->SetMarkerColor(col);
      gp.g->SetLineColor(col);
      gp.g->SetLineWidth(2);
      gp.g->SetMarkerStyle(20 + (gp.zBin % 6));
      gp.g->SetMarkerSize(0.9);

      // --- draw points with error bars ---
      gp.g->Draw("E1P SAME");   // E1 = error bars, P = points

      // --- fit: sigma(phi) = sigma0 * (1 + A cos(phi) + B cos(2phi)) ---
      // Only fit if enough points
      if (gp.g->GetN() >= 4) {
        TF1* f = new TF1(Form("f_pt%d_z%d", gp.ptBin, gp.zBin),
                         "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))",
                         0.0, 2.0*TMath::Pi());

        // start guesses (stable defaults)
        double ymean = 0.0;
        for (int ip = 0; ip < gp.g->GetN(); ip++) {
          double xx, yy;
          gp.g->GetPoint(ip, xx, yy);
          ymean += yy;
        }
        ymean /= gp.g->GetN();


        // Set Parameter limit
        f->SetParLimits(0, 0.0, 1e9);   // sigma0 positive
        f->SetParLimits(1, -1.0, 1.0);  // A
        f->SetParLimits(2, -1.0, 1.0);  // B

        f->SetParameters(ymean, 0.0, 0.0);
        f->SetParNames("sigma0", "Acos", "Bcos2");
        f->SetLineColor(col);
        f->SetLineWidth(2);

        // Fit
        //gp.g->Fit(f, "RQ0");
        if (gp.g->GetN() >= 4) {
          gp.g->Fit(f, "RQ0");
          f->Draw("L SAME");
        } else {
          std::cout << "Skip fit: pt[" << gp.ptBin << "] z[" << gp.zBin
                    << "] has only N=" << gp.g->GetN() << " points\n";
        }

        fitFuncs.push_back(f);
        std::cout << Form("pt[%d] z[%d]: sigma0=%.3g A=%.3g B=%.3g chi2/ndf=%.2f\n",
                          gp.ptBin, gp.zBin,
                          f->GetParameter(0), f->GetParameter(1), f->GetParameter(2),
                          f->GetChisquare() / std::max(1.0, (double)f->GetNDF()));
      }



      // Legend: keep it simple (one entry per z)
      leg->AddEntry(gp.g, Form("z #in [%.2f, %.2f]", zLo, zHi), "pe");

      drawn++;
    }

    leg->Draw();
  }

  if (savePNGs) {
    std::string outPng = resultsDir + "/PNGs/xsec_vs_phipq_panels.png";
    c->SaveAs(outPng.c_str());
    log << "Saved: " << outPng << "\n";
  }

  log << "Done.\n";
  log.close();
  if (csv.is_open()) csv.close();
}
