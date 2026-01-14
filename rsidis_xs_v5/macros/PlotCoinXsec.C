// PlotCoinXsec.C (rsidis_xs_v5)
// v5 goal:
//   - Single-setting mode (same as v4): xsec vs phipq panels (pads = pT bins, curves = z bins)
//   - Group overlay mode (.list): xsec vs phipq panels (pads = pT bins, curves = settings in the group file)
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
// Run Command Example: (group mode)
// root -l -b -q 'macros/PlotCoinXsec.C("/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/groups/pass4/pi+sidis/LH2/x0p25/q23p3/tpq2p0/grp_LH2_zOv_x0p25_Q23p3_tpq2p0.list", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/results", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/settings")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <cctype>
#include <cstdlib>
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

static bool EndsWith(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

static std::string StripExtension(const std::string& path) {
  std::string base = Basename(path);
  auto pos = base.find_last_of('.');
  if (pos == std::string::npos) return base;
  return base.substr(0, pos);
}

// Ensure a CSV file exists and has a header line.
// If the file does not exist (or is empty), write the header.
static void EnsureCSVHeader(const std::string& path, const std::string& header) {
  bool need = false;
  if (gSystem->AccessPathName(path.c_str())) {
    need = true;
  } else {
    std::ifstream in(path);
    if (!in.good()) {
      need = true;
    } else {
      in.seekg(0, std::ios::end);
      if (in.tellg() == 0) need = true;
    }
  }
  if (need) {
    std::ofstream out(path);
    out << header << "\n";
  }
}


// Parse encoded number tokens like "0p25" -> 0.25, "neg1p531" -> -1.531
static bool ParseEncodedNumber(const std::string& enc, double& out) {
  if (enc.empty()) return false;
  std::string s = enc;
  bool neg = false;
  if (s.rfind("neg", 0) == 0) { neg = true; s = s.substr(3); }
  // replace 'p' with '.'
  for (auto& ch : s) if (ch=='p') ch='.';
  char* endp = nullptr;
  double v = std::strtod(s.c_str(), &endp);
  if (endp == s.c_str()) return false;
  out = neg ? -v : v;
  return true;
}

// Extract nominal xB, Q2, z from the settingId path tokens like ".../z0p5/x0p25/q23p3/..."
static void ExtractNominalKinFromSettingId(const std::string& settingId, double& xB, double& Q2, double& z) {
  xB = NAN; Q2 = NAN; z = NAN;
  std::vector<std::string> parts;
  {
    std::string tmp; tmp.reserve(settingId.size());
    for (char c : settingId) {
      if (c=='/') { if(!tmp.empty()) { parts.push_back(tmp); tmp.clear(); } }
      else tmp.push_back(c);
    }
    if (!tmp.empty()) parts.push_back(tmp);
  }
  for (const auto& p : parts) {
    if (p.size() >= 2 && p[0]=='x') {
      double v; if (ParseEncodedNumber(p.substr(1), v)) xB = v;
    } else if (p.rfind("q2",0)==0) {
      double v; if (ParseEncodedNumber(p.substr(2), v)) Q2 = v;
    } else if (p.size() >= 2 && p[0]=='z') {
      double v; if (ParseEncodedNumber(p.substr(1), v)) z = v;
    }
  }
}


struct RunNormInfo {
  double charge_mC = 0.0;   // BCM4A_Q
  double xB = NAN;           // optional: x_Bj (per-run)
  double Q2 = NAN;           // optional: Q2 (per-run)
  double hms_eff   = 1.0;   // h_esing_Eff
  double ps_factor = 1.0;   // ps5
  double comp_lt   = 1.0;   // optional
};

static double ChargeWeightedMean(const std::vector<int>& runs,
                                 const std::unordered_map<int, RunNormInfo>& runInfo,
                                 bool wantXB)
{
  double sumQ = 0.0;
  double sumQx = 0.0;
  for (int run : runs) {
    auto it = runInfo.find(run);
    if (it == runInfo.end()) continue;
    const RunNormInfo& info = it->second;
    const double Q = info.charge_mC;
    const double v = wantXB ? info.xB : info.Q2;
    if (!std::isfinite(Q) || Q<=0) continue;
    if (!std::isfinite(v)) continue;
    sumQ  += Q;
    sumQx += Q*v;
  }
  if (sumQ <= 0) return NAN;
  return sumQx / sumQ;
}

struct GroupEntry {
  std::string label;
  std::string manifestAbs;
  std::string leafDir;
  std::string settingId;
};

static std::vector<GroupEntry> ReadGroupListFile(const std::string& groupPath,
                                                 const std::string& settingsRoot,
                                                 std::ostream& log)
{
  std::vector<GroupEntry> out;
  std::ifstream f(groupPath);
  if (!f.is_open()) {
    log << "ERROR: cannot open group file: " << groupPath << "\n";
    return out;
  }
  const std::string baseDir = NormalizeSlashes(Dirname(groupPath));
  std::string line;
  int ln = 0;
  while (std::getline(f, line)) {
    ++ln;
    // trim
    auto trim = [](std::string& s){
      while(!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
      while(!s.empty() && std::isspace((unsigned char)s.back()))  s.pop_back();
    };
    trim(line);
    if (line.empty() || line[0]=='#') continue;

    std::istringstream iss(line);
    std::string label, mpath;
    if (!(iss >> label >> mpath)) {
      log << "WARNING: malformed group line " << ln << ": " << line << "\n";
      continue;
    }
    std::string abs = mpath;
    if (!abs.empty() && abs[0] != '/') {
      abs = NormalizeSlashes(baseDir + "/" + abs);
    } else {
      abs = NormalizeSlashes(abs);
    }

    GroupEntry e;
    e.label = label;
    e.manifestAbs = abs;
    e.leafDir = NormalizeSlashes(Dirname(abs));
    e.settingId = MakeSettingIdFromManifestPath(abs, settingsRoot);
    out.push_back(e);
  }
  return out;
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

// Overload that logs and delegates to the single-argument version.
static std::vector<int> ReadRunList(const std::string& path, std::ostream& log) {
  log << "Reading runlist: " << path << "\n";
  return ReadRunList(path);
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

static void PlotCoinXsec_Single(const char* manifestPath,
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


  std::ofstream fitcsv;
  if (saveTables) {
    const std::string fitPath = resultsDir + "/tables/fit_parameters.csv";
    EnsureCSVHeader(fitPath,
      "mode,group_id,curve_label,setting_id,pt_bin,z_bin,pt_lo,pt_hi,pt_center,z_lo,z_hi,z_center,n_points,"
      "M0,M0_err,A,A_err,B,B_err,chi2,ndf,prob");
    fitcsv.open(fitPath, std::ios::app);
    log << "Fit CSV: " << fitPath << "\n";
  }
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
          ex.push_back(0.0);
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

        // Log + persist fit parameters
        double ptLoF = binning.ptEdges[gp.ptBin];
        double ptHiF = binning.ptEdges[gp.ptBin + 1];
        double ptCF  = Center(ptLoF, ptHiF);
        double zLoF  = binning.zEdges[gp.zBin];
        double zHiF  = binning.zEdges[gp.zBin + 1];
        double zCF   = Center(zLoF, zHiF);

        log << "FIT single setting_id=" << settingId
            << " pt_bin=" << gp.ptBin << " z_bin=" << gp.zBin
            << " N=" << gp.g->GetN()
            << " M0=" << f->GetParameter(0) << " +/- " << f->GetParError(0)
            << " A="  << f->GetParameter(1) << " +/- " << f->GetParError(1)
            << " B="  << f->GetParameter(2) << " +/- " << f->GetParError(2)
            << " chi2=" << f->GetChisquare() << " ndf=" << f->GetNDF()
            << " prob=" << f->GetProb() << "\n";

        if (saveTables && fitcsv.is_open()) {
          fitcsv
            << "single" << ","            // mode
            << "" << ","                  // group_id
            << Form("zbin%d", gp.zBin) << "," // curve_label
            << settingId << ","
            << gp.ptBin << "," << gp.zBin << ","
            << ptLoF << "," << ptHiF << "," << ptCF << ","
            << zLoF  << "," << zHiF  << "," << zCF << ","
            << gp.g->GetN() << ","
            << f->GetParameter(0) << "," << f->GetParError(0) << ","
            << f->GetParameter(1) << "," << f->GetParError(1) << ","
            << f->GetParameter(2) << "," << f->GetParError(2) << ","
            << f->GetChisquare() << "," << f->GetNDF() << "," << f->GetProb()
            << "\n";
        }
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



static void PlotCoinXsec_Group(const char* groupPath,
                               const char* resultsRoot = "results",
                               const char* settingsRoot = "settings",
                               bool savePNGs = true,
                               bool saveTables = true)
{
  const std::string groupFile = groupPath ? std::string(groupPath) : "";
  if (groupFile.empty()) {
    std::cerr << "ERROR: groupPath is empty\n";
    return;
  }

  const std::string groupId  = StripExtension(groupFile);
  const std::string resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + groupId);
  MkdirP(resultsDir);
  MkdirP(resultsDir + "/PNGs");
  MkdirP(resultsDir + "/tables");

  std::ofstream log(resultsDir + "/log.txt", std::ios::app);
  log << "=== PlotCoinXsec v5 (GROUP overlay) ===\n";
  log << "groupFile: " << groupFile << "\n";
  log << "resultsDir: " << resultsDir << "\n";
  log << "settingsRoot: " << settingsRoot << "\n";


  // Keep group-mode fit functions alive until canvas is written
  std::vector<TF1*> fitFuncs;
  fitFuncs.reserve(200);

  std::ofstream fitcsv;
  if (saveTables) {
    const std::string fitPath = resultsDir + "/tables/fit_parameters.csv";
    EnsureCSVHeader(fitPath,
      "mode,group_id,curve_label,setting_id,pt_bin,z_bin,pt_lo,pt_hi,pt_center,z_lo,z_hi,z_center,n_points,"
      "M0,M0_err,A,A_err,B,B_err,chi2,ndf,prob");
    fitcsv.open(fitPath, std::ios::app);
    log << "Fit CSV: " << fitPath << "\n";
  }

  // Read group entries
  auto entries = ReadGroupListFile(groupFile, settingsRoot, log);
  if (entries.empty()) {
    log << "ERROR: group file produced 0 entries.\n";
    std::cerr << "ERROR: group file produced 0 entries.\n";
    return;
  }
  log << "Group entries: " << entries.size() << "\n";
  for (auto& e : entries) {
    log << "  " << e.label << "  settingId=" << e.settingId << "  leafDir=" << e.leafDir << "\n";
  }

  // Random subtraction config
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  // Binning defaults:
  // - In group overlay mode, start with ONE wide pT bin and ONE wide z bin (z dimension comes from curves)
  Binning binning;
  //const bool overlayWidePt = true;   // set false later when you want 4 pads (0-0.1-0.2-0.3-0.4)
  const bool overlayWidePt = false;   // set false later when you want 4 pads (0-0.1-0.2-0.3-0.4)
  if (overlayWidePt) binning.ptEdges = {0.0, 0.40};
  binning.zEdges = {0.0, 1.0};       // integrate over z (curves are different z SETTINGS)

  const int nCurves = (int)entries.size();
  const int nPt  = (int)binning.ptEdges.size()  - 1;
  const int nPhi = (int)binning.phiEdges.size() - 1;

  // DATA base cuts (same as v4)
  const TCut dataCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_cal_etottracknorm<0.8)";

  // DATA vars (skimmed)
  const std::string data_pt  = "pt";
  const std::string data_z   = "z";
  // Wrap phipq into [0,2pi) WITHOUT ternary (ROOT-safe)
  const std::string data_phi =
    "((P_kin_secondary_ph_xq>=0)*P_kin_secondary_ph_xq + (P_kin_secondary_ph_xq<0)*(P_kin_secondary_ph_xq + 2*TMath::Pi()))";

  const double wall_thickness_ratio = 3.82;
  const double dumScale = 1.0 / wall_thickness_ratio;

  // SIM vars (SIMC h10)
  const std::string sim_pt   = "sqrt(pt2)";
  const std::string sim_z    = "z";
  const std::string sim_phi  =
    "((phipq>=0)*phipq + (phipq<0)*(phipq + 2*TMath::Pi()))";

  const TCut simCutsBase =
    "(hsdelta>-8.0) && (hsdelta<8.0) && "
    "(ssdelta>-10.0) && (ssdelta<22.0)";

  // SIM normalization (same as v4)
  const double simNormalizationFactor = 0.842205e11;

  // ---- Overlay validity guards (two layers) ----
  // Layer A (MC health): require non-empty SIMC yield and reasonable MC relative uncertainty
  const double kMinSimYield     = 0.05;  // in the same units as yield_sim in the CSV
  const double kMaxRelSimErr    = 0.35;  // eSim/ySim
  // Layer B (xsec quality): require reasonable relative uncertainty on the extracted xsec
  const double kMaxRelXsecErr   = 0.50;  // xsec_err/|xsec|


  // Per-curve loaded inputs
  struct CurveInputs {
    std::string label;
    std::string leafDir;
    std::string settingId;
    std::vector<int> runs_sig, runs_dum, runs_pos, runs_posdum;
    std::unordered_map<int, RunNormInfo> runInfo;
    double sigmaModel = 0.0;
    double meanXB = NAN;
    double meanQ2 = NAN;
    double nominalZ = NAN;

    std::unique_ptr<TChain> sim;
    Long64_t nGenSim = 0;

    // Cache coincidence random subtraction windows per run
    std::unordered_map<int, CoincidenceResult> coinCache;
  };

  std::vector<CurveInputs> curves;
  curves.reserve(nCurves);

  for (const auto& e : entries) {
    CurveInputs c;
    c.label = e.label;
    c.leafDir = e.leafDir;
    c.settingId = e.settingId;

    // Run lists
    c.runs_sig    = ReadRunList(c.leafDir + "/runs_signal.txt", log);
    c.runs_dum    = ReadRunList(c.leafDir + "/runs_dummy.txt", log);
    c.runs_pos    = ReadRunList(c.leafDir + "/runs_positron.txt", log);
    c.runs_posdum = ReadRunList(c.leafDir + "/runs_positron_dummy.txt", log);

    if (!LoadRunMetadataCSV(c.leafDir + "/run_metadata.csv", c.runInfo)) {
      log << "ERROR: cannot load run_metadata.csv for " << c.settingId << "\n";
      std::cerr << "ERROR: cannot load run_metadata.csv for " << c.settingId << "\n";
      return;
    }

    // sigma model
    c.sigmaModel = ReadSigmaModel(c.leafDir, log);

    // nominal kin from settingId
    double xNom, qNom, zNom;
    ExtractNominalKinFromSettingId(c.settingId, xNom, qNom, zNom);
    c.nominalZ = zNom;

    // charge-weighted means if columns exist
    c.meanXB = ChargeWeightedMean(c.runs_sig, c.runInfo, true);
    c.meanQ2 = ChargeWeightedMean(c.runs_sig, c.runInfo, false);
    if (!std::isfinite(c.meanXB)) c.meanXB = xNom;
    if (!std::isfinite(c.meanQ2)) c.meanQ2 = qNom;

    // SIM chain
    c.sim.reset(new TChain("h10"));
    std::string simGlob = "./simc_worksim/" + c.settingId + "/*.root";
    int nAdded = c.sim->Add(simGlob.c_str());
    log << "SIM glob (" << c.label << "): " << simGlob << " added=" << nAdded << "\n";
    if (c.sim->GetEntries() <= 0) {
      log << "ERROR: SIM chain empty for: " << simGlob << "\n";
      std::cerr << "ERROR: SIM chain empty for: " << simGlob << "\n";
      return;
    }
    c.nGenSim = c.sim->GetEntries();
    curves.push_back(std::move(c));
  }

  // CSV (group overlay)
  std::ofstream csv;
  if (saveTables) {
    std::string csvPath = resultsDir + "/tables/xsec_phipq_overlay.csv";
    csv.open(csvPath);
    csv
      << "group_id,curve_label,setting_id,pt_bin,phibin,"
      << "pt_lo,pt_hi,pt_center,phi_lo,phi_hi,phi_center,"
      << "yield_sig,yield_pos,yield_dum,yield_posdum,yield_net,yield_net_err,"
      << "yield_sim,yield_sim_err,"
      << "xsec,xsec_stat,mean_xB,mean_Q2,valid\n";
    log << "CSV: " << csvPath << "\n";
  }

  // Canvas (pads = pT bins)
  int nCols = (nPt <= 2) ? nPt : 2;
  int nRows = (int)std::ceil((double)nPt / nCols);
  TCanvas* c = new TCanvas("c_phipq_overlay", "xsec vs phipq (z-setting overlay)", 1200, 450*nRows);
  c->Divide(nCols, nRows, 0.001, 0.001);

  std::vector<double> yMaxByPt(nPt, 0.0);

  // Color palette
  static const Color_t curveColors[] = { kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7 };
  static const int nCurveColors = sizeof(curveColors)/sizeof(curveColors[0]);

  // Loop over pT pads
  for (int iPt=0; iPt<nPt; iPt++) {
    double ptLo = binning.ptEdges[iPt];
    double ptHi = binning.ptEdges[iPt+1];
    double ptCtr= Center(ptLo, ptHi);

    c->cd(iPt+1);
    gPad->SetGrid(1,1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.12);

    TH1F* frame = (TH1F*) gPad->DrawFrame(0.0, 0.0, 2.0*M_PI, 1.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("#phi_{pq} (rad)");
    frame->GetYaxis()->SetTitle("d#sigma");
    frame->GetYaxis()->SetTitleOffset(1.2);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22);
    latex.DrawLatex(0.50, 0.92, Form("p_{T} #in [%.2f, %.2f] GeV", ptLo, ptHi));

    TLegend* leg = new TLegend(0.52, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);

    // Loop over curves (different z SETTINGS)
    for (int ic=0; ic<nCurves; ic++) {
      auto& cv = curves[ic];

      // Per-category yields for a given kin cut
      double qSig=0, qPos=0, qDum=0, qPosDum=0;

      auto CategoryYield = [&](const std::vector<int>& runs,
                               const TCut& baseCutsBool,
                               const TCut& kinCutsBool,
                               double& totalCharge_mC) -> Yield
      {
        Yield out;
        totalCharge_mC = 0.0;
        int nOpened = 0;
        bool didSanity = false;

        for (int run : runs) {
          auto it = cv.runInfo.find(run);
          if (it == cv.runInfo.end()) continue;
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
            if (f) f->Close();
            continue;
          }
          TTree* T = (TTree*) f->Get("T");
          if (!T) {
            log << "WARNING: missing tree T in " << path << "\n";
            f->Close();
            continue;
          }

          ++nOpened;
          const TString baseCutsStr = TString(baseCutsBool.GetTitle());
          const TString fullCutsStr = TString((baseCutsBool && kinCutsBool).GetTitle());

          if (!didSanity) {
            Long64_t nAll  = T->GetEntries();
            Long64_t nBase = T->GetEntries(baseCutsStr);
            Long64_t nFull = T->GetEntries(fullCutsStr);
            log << "SANITY (" << cv.label << ", run " << run << "): entries(all)=" << nAll
                << " entries(baseCuts)=" << nBase
                << " entries(fullCuts)=" << nFull << "\n";
            didSanity = true;
          }

          // Coin window cached per run
          auto itC = cv.coinCache.find(run);
          if (itC == cv.coinCache.end()) {
            Bool_t oldAdd = TH1::AddDirectoryStatus();
            TH1::AddDirectory(kTRUE);
            CoincidenceResult win = ComputeCoincidenceRandomSubtraction(T, baseCutsStr, cfg);
            TH1::AddDirectory(oldAdd);
            cv.coinCache[run] = win;
            itC = cv.coinCache.find(run);
          }

          TString wexpr = Form("%g", w);
          Yield y = AddYieldFromRun(T, fullCutsStr, wexpr, cfg, itC->second);
          out.sumw  += y.sumw;
          out.sumw2 += y.sumw2;

          f->Close();
        }

        if (nOpened == 0) {
          log << "ERROR: opened 0 ROOT files for category in curve " << cv.label << "\n";
        }

        if (totalCharge_mC > 0) {
          out.sumw  /= totalCharge_mC;
          out.sumw2 /= (totalCharge_mC * totalCharge_mC);
        }
        return out;
      };

      // Build graph points for this curve in this pT bin
      std::vector<double> x, y, ex, ey;
      x.reserve(nPhi); y.reserve(nPhi); ex.reserve(nPhi); ey.reserve(nPhi);

      for (int iPhi=0; iPhi<nPhi; iPhi++) {
        double phiLo = binning.phiEdges[iPhi];
        double phiHi = binning.phiEdges[iPhi+1];
        double phiCtr= Center(phiLo, phiHi);

        // z integrated (wide)
        double zLo = binning.zEdges[0];
        double zHi = binning.zEdges[1];

        TCut kinData = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                            data_pt.c_str(), ptLo, data_pt.c_str(), ptHi,
                            data_z.c_str(),  zLo,  data_z.c_str(),  zHi,
                            data_phi.c_str(),phiLo, data_phi.c_str(),phiHi);

        Yield ySig    = CategoryYield(cv.runs_sig,    dataCuts, kinData, qSig);
        Yield yPos    = CategoryYield(cv.runs_pos,    dataCuts, kinData, qPos);
        Yield yDum    = CategoryYield(cv.runs_dum,    dataCuts, kinData, qDum);
        Yield yPosDum = CategoryYield(cv.runs_posdum, dataCuts, kinData, qPosDum);

        // (Data − PosData) − (Dummy − PosDummy)/wall_thickness_ratio
        double yNet  = (ySig.sumw - yPos.sumw) - dumScale * (yDum.sumw - yPosDum.sumw);
        double eNet2 = ySig.sumw2 + yPos.sumw2
                     + dumScale*dumScale * (yDum.sumw2 + yPosDum.sumw2);
        double eNet  = (eNet2 > 0) ? std::sqrt(eNet2) : 0.0;

        // SIM yield
        TCut kinSim = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                           sim_pt.c_str(), ptLo, sim_pt.c_str(), ptHi,
                           sim_z.c_str(),  zLo,  sim_z.c_str(),  zHi,
                           sim_phi.c_str(),phiLo, sim_phi.c_str(),phiHi);

        const double simTryScale = (cv.nGenSim > 0) ? (1.0 / double(cv.nGenSim)) : 0.0;
        const TString simWeightExpr = Form("Weight*%e*%e", simNormalizationFactor, simTryScale);

        Yield ySim = SimYieldSum(*cv.sim, simWeightExpr, TString((simCutsBase && kinSim).GetTitle()));
        double eSim = (ySim.sumw2 > 0) ? std::sqrt(ySim.sumw2) : 0.0;

        // ----------------------------
        // Two-layer validity guard (GROUP overlay mode)
        // ----------------------------

        // Layer A: MC health check
        const bool mc_zero = (std::abs(ySim.sumw) < 1e-15);
        const double relSim = (ySim.sumw != 0.0) ? (std::abs(eSim / ySim.sumw)) : 1e9;

        bool valid = true;
        if (mc_zero) valid = false;
        if (!std::isfinite(ySim.sumw) || !std::isfinite(eSim)) valid = false;
        if (ySim.sumw <= kMinSimYield) valid = false;
        if (relSim > kMaxRelSimErr) valid = false;

        // Compute xsec (even if invalid) for debugging/CSV, but plot only if valid.
        double xsec = 0.0, xsecErr = 0.0;
        if (!mc_zero && ySim.sumw != 0.0) {
          xsec = (yNet / ySim.sumw) * cv.sigmaModel;

          double rel2 = 0.0;
          if (std::abs(yNet) > 0) rel2 += (eNet / yNet) * (eNet / yNet);
          if (std::abs(ySim.sumw) > 0) rel2 += (eSim / ySim.sumw) * (eSim / ySim.sumw);
          xsecErr = std::abs(xsec) * std::sqrt(rel2);
        }

        // Layer B: xsec quality check
        if (!std::isfinite(xsec) || !std::isfinite(xsecErr)) valid = false;
        if (xsec <= 0.0) valid = false; // avoid plotting negative/subtracted bins in overlay view
        const double relX = (std::abs(xsec) > 0.0) ? (xsecErr / std::abs(xsec)) : 1e9;
        if (relX > kMaxRelXsecErr) valid = false;

        // Keep only valid points for plotting/fitting and for y-axis autoscale
        if (valid) {
          x.push_back(phiCtr);
          ex.push_back(0.0);
          y.push_back(xsec);
          ey.push_back(xsecErr);
          yMaxByPt[iPt] = std::max(yMaxByPt[iPt], xsec + xsecErr);
        }

        if (saveTables && csv.is_open()) {
          csv
            << groupId << ","
            << cv.label << ","
            << cv.settingId << ","
            << iPt << "," << iPhi << ","
            << ptLo << "," << ptHi << "," << ptCtr << ","
            << phiLo << "," << phiHi << "," << phiCtr << ","
            << ySig.sumw << "," << yPos.sumw << "," << yDum.sumw << "," << yPosDum.sumw << ","
            << yNet << "," << eNet << ","
            << ySim.sumw << "," << eSim << ","
            << xsec << "," << xsecErr << ","
            << cv.meanXB << "," << cv.meanQ2 << ","
            << (valid ? 1 : 0)
            << "\n";
        }
      } // phi loop

      // draw graph
      TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.9);
      g->SetLineWidth(2);
      Color_t col = curveColors[ic % nCurveColors];
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->Draw("E1P SAME");


      // Fit: sigma(phi) = M0 * (1 + A cos(phi) + B cos(2phi))
      if (g->GetN() >= 4) {
        TF1* f = new TF1(Form("f_grp_pt%d_c%d", iPt, ic),
                         "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))",
                         0.0, 2.0*TMath::Pi());

        // initial guesses
        double ymean = 0.0;
        for (int ip = 0; ip < g->GetN(); ip++) {
          double xx, yy;
          g->GetPoint(ip, xx, yy);
          ymean += yy;
        }
        ymean /= g->GetN();

        f->SetParLimits(0, 0.0, 1e9);
        f->SetParLimits(1, -1.0, 1.0);
        f->SetParLimits(2, -1.0, 1.0);
        f->SetParameters(ymean, 0.0, 0.0);
        f->SetParNames("M0", "Acos", "Bcos2");
        f->SetLineColor(col);
        f->SetLineWidth(2);

        g->Fit(f, "RQ0");
        f->Draw("L SAME");
        fitFuncs.push_back(f);

        log << "FIT group group_id=" << groupId
            << " curve=" << cv.label
            << " setting_id=" << cv.settingId
            << " pt_bin=" << iPt
            << " N=" << g->GetN()
            << " M0=" << f->GetParameter(0) << " +/- " << f->GetParError(0)
            << " A="  << f->GetParameter(1) << " +/- " << f->GetParError(1)
            << " B="  << f->GetParameter(2) << " +/- " << f->GetParError(2)
            << " chi2=" << f->GetChisquare() << " ndf=" << f->GetNDF()
            << " prob=" << f->GetProb() << "\n";

        if (saveTables && fitcsv.is_open()) {
          fitcsv
            << "group" << ","
            << groupId << ","
            << cv.label << ","
            << cv.settingId << ","
            << iPt << "," << -1 << ","
            << ptLo << "," << ptHi << "," << ptCtr << ","
            << cv.nominalZ << "," << cv.nominalZ << "," << cv.nominalZ << ","
            << g->GetN() << ","
            << f->GetParameter(0) << "," << f->GetParError(0) << ","
            << f->GetParameter(1) << "," << f->GetParError(1) << ","
            << f->GetParameter(2) << "," << f->GetParError(2) << ","
            << f->GetChisquare() << "," << f->GetNDF() << "," << f->GetProb()
            << "\n";
        }
      } else {
        log << "Skip fit (group): curve=" << cv.label << " pt_bin=" << iPt << " has N=" << g->GetN() << "\n";
      }

      // legend label
      TString leglab = Form("%s  #LT x #GT=%.3f  #LTQ^{2}#GT=%.2f", cv.label.c_str(), cv.meanXB, cv.meanQ2);
      leg->AddEntry(g, leglab, "pe");
    } // curve loop

    // scale y
    double ymax = (yMaxByPt[iPt] > 0) ? (1.20 * yMaxByPt[iPt]) : 1.0;
    frame->GetYaxis()->SetRangeUser(0.0, ymax);

    leg->Draw();
  } // pt loop

  if (savePNGs) {
    std::string outPng = resultsDir + "/PNGs/xsec_vs_phipq_overlay.png";
    c->SaveAs(outPng.c_str());
    log << "Saved: " << outPng << "\n";
  }

  log << "Done (GROUP overlay).\n";
  log.close();
  if (csv.is_open()) csv.close();
}

// Wrapper: same call signature as v4, but supports group .list files in v5
void PlotCoinXsec(const char* manifestPath,
                  const char* resultsRoot = "results",
                  const char* settingsRoot = "settings",
                  bool savePNGs = true,
                  bool saveTables = true)
{
  if (!manifestPath || std::string(manifestPath).empty()) {
    std::cerr << "ERROR: manifestPath is required\n";
    return;
  }
  std::string p = manifestPath;
  if (EndsWith(p, ".list")) {
    PlotCoinXsec_Group(manifestPath, resultsRoot, settingsRoot, savePNGs, saveTables);
  } else {
    PlotCoinXsec_Single(manifestPath, resultsRoot, settingsRoot, savePNGs, saveTables);
  }
}

