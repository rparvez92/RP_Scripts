// TableCoinXsec.C (rsidis_xs_v5)
//
// Purpose (v5):
//   Produce ONLY the cross-section table(s) used by downstream plot/fit macros.
//   This is a compute-heavy step; visualization should read these CSVs and apply
//   filtering/guards/outlier handling separately.
//
// Outputs:
//   Single-setting mode (manifestPath is a setting manifest):
//     rsidis_xs_v5/results/<setting_id>/tables/xsec_phipq_z_pt_overlayed_single.csv
//
//   Group overlay mode (manifestPath ends with .list):
//     rsidis_xs_v5/results/<group_id>/tables/xsec_phipq_z_pt_overlayed_group.csv
//
// Notes:
//   - Keeps ALL bins in the CSV (including MC-starved, negative, large-err, etc.).
//     Plotting macros decide what to draw.
//   - Uses a flag bitmask (flag_bits) plus a small valid_default field.
//   - Uses a sentinel (kMissing = -999) ONLY for values that are undefined/uncomputed
//     (e.g., xsec when ySim<=0, relative errors when denominator is 0, etc.).
//     "Bad quality" but finite values are preserved.
//
// Run examples:
//   root -l -b -q 'macros/TableCoinXsec.C("settings/.../manifest.txt", "results", "settings")'
//   root -l -b -q 'macros/TableCoinXsec.C("groups/.../grp_*.list", "results", "settings")'
//   root -l -b -q 'macros/TableCoinXsec.C("/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/groups/pass4/pi+sidis/LH2/x0p25/q23p3/thpq2/grp_pass4_piplus_LH2_zOv_x0p25_q23p3_thpq2.list", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/results", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/settings")'

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
#include "TMath.h"

#include "../include/CoincidenceRandomSubtraction.h"

namespace {

// ---------- small utilities copied/adapted from PlotCoinXsec.C ----------

static const double kMissing = -999.0;

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

static bool EndsWith(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

static std::string StripExtension(const std::string& path) {
  std::string base = Basename(path);
  auto pos = base.find_last_of('.');
  if (pos == std::string::npos) return base;
  return base.substr(0, pos);
}

static std::string SwapExtension(const std::string& path, const std::string& newExtNoDot) {
  auto pos = path.find_last_of('.');
  if (pos == std::string::npos) return path + "." + newExtNoDot;
  return path.substr(0, pos + 1) + newExtNoDot;
}

static std::string GetFirstRootFileFromChain(TChain& ch, std::ofstream& log, const std::string& simGlobForMsg) {
  // Ensure the first file is actually opened.
  if (ch.GetEntries() <= 0) {
    log << "ERROR: SIM chain has 0 entries for: " << simGlobForMsg << "\n";
    std::cerr << "ERROR: SIM chain has 0 entries for: " << simGlobForMsg << "\n";
    gSystem->Exit(1);
  }
  ch.GetEntry(0);
  TFile* f = ch.GetCurrentFile();
  if (!f) {
    log << "ERROR: could not open current SIM ROOT file for: " << simGlobForMsg << "\n";
    std::cerr << "ERROR: could not open current SIM ROOT file for: " << simGlobForMsg << "\n";
    gSystem->Exit(1);
  }
  return NormalizeSlashes(std::string(f->GetName()));
}

static double GrabNormfacFromHist(const std::string& simRootPath,
                                  std::string* histPathOut,
                                  std::ofstream& log)
{
  std::string histPath = simRootPath;
  if (EndsWith(histPath, ".root")) {
    histPath = SwapExtension(histPath, "hist");
  } else {
    // Defensive: still try to read a .hist next to whatever we were given.
    histPath = simRootPath + ".hist";
  }
  if (histPathOut) *histPathOut = histPath;

  std::ifstream f(histPath);
  if (!f) {
    log << "ERROR: cannot open SIMC .hist file: " << histPath << "\n";
    std::cerr << "ERROR: cannot open SIMC .hist file: " << histPath << "\n";
    gSystem->Exit(1);
  }

  bool inMisc = false;
  std::string line;
  while (std::getline(f, line)) {
    if (!inMisc) {
      if (line.find("MISCELLANEOUS") != std::string::npos) inMisc = true;
      continue;
    }

    // normfac is always within MISCELLANEOUS, before RECON SUMMARY in SIMC output.
    if (line.find("RECON SUMMARY") != std::string::npos) break;

    auto p = line.find("normfac");
    if (p == std::string::npos) continue;
    auto eq = line.find('=', p);
    if (eq == std::string::npos) continue;

    std::string rhs = line.substr(eq + 1);
    std::istringstream iss(rhs);
    double v = 0.0;
    if (!(iss >> v)) continue;
    if (!std::isfinite(v) || v <= 0.0) {
      log << "ERROR: parsed non-finite or non-positive normfac=" << v << " from: " << histPath << "\n";
      std::cerr << "ERROR: parsed non-finite or non-positive normfac=" << v << " from: " << histPath << "\n";
      gSystem->Exit(1);
    }
    return v;
  }

  log << "ERROR: could not find 'normfac' under MISCELLANEOUS in: " << histPath << "\n";
  std::cerr << "ERROR: could not find 'normfac' under MISCELLANEOUS in: " << histPath << "\n";
  gSystem->Exit(1);
  return 0.0; // unreachable
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

// Extract nominal xB, Q2, z from settingId tokens like ".../z0p5/x0p25/q23p3/..."
static bool ParseEncodedNumber(const std::string& enc, double& out) {
  if (enc.empty()) return false;
  std::string s = enc;
  bool neg = false;
  if (s.rfind("neg", 0) == 0) { neg = true; s = s.substr(3); }
  for (auto& ch : s) if (ch=='p') ch='.';
  char* endp = nullptr;
  double v = std::strtod(s.c_str(), &endp);
  if (endp == s.c_str()) return false;
  out = neg ? -v : v;
  return true;
}

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
  double xB = NAN;          // optional
  double Q2 = NAN;          // optional
  double hms_eff   = 1.0;   // h_esing_Eff
  double ps_factor = 1.0;   // ps5
  double comp_lt   = 1.0;   // comp_livetime
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
  // Optional (if present in future tables)
  int i_xb    = idx("x_Bj");
  int i_q2    = idx("Q2");

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto v = SplitCSVLine(line);
    if (i_run < 0 || i_run >= (int)v.size()) continue;

    int run = std::atoi(v[i_run].c_str());
    RunNormInfo info;

    if (i_q  >= 0 && i_q  < (int)v.size()) info.charge_mC = std::atof(v[i_q].c_str());
    if (i_eff>= 0 && i_eff< (int)v.size()) info.hms_eff   = std::atof(v[i_eff].c_str());
    if (i_ps >= 0 && i_ps < (int)v.size()) info.ps_factor = std::atof(v[i_ps].c_str());
    if (i_lt >= 0 && i_lt < (int)v.size()) info.comp_lt   = std::atof(v[i_lt].c_str());
    if (i_xb >= 0 && i_xb < (int)v.size()) info.xB        = std::atof(v[i_xb].c_str());
    if (i_q2 >= 0 && i_q2 < (int)v.size()) info.Q2        = std::atof(v[i_q2].c_str());

    if (!std::isfinite(info.hms_eff)   || info.hms_eff <= 0)   info.hms_eff   = 1.0;
    if (!std::isfinite(info.ps_factor) || info.ps_factor <= 0) info.ps_factor = 1.0;
    if (!std::isfinite(info.comp_lt)   || info.comp_lt <= 0)   info.comp_lt   = 1.0;

    out[run] = info;
  }
  return true;
}

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

static double ReadSigmaModel(const std::string& leafDir, std::ostream& log) {
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
  std::snprintf(buf, sizeof(buf), "./Pass0_SkimmedDataROOTfiles/skimmed_coin_replay_production_%d_-1.root", run);
  return std::string(buf);
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
  std::vector<double> phiEdges = PhiEdgesN(16);
};

// ---------- yields ----------

struct YieldEx {
  double sumw  = 0.0;
  double sumw2 = 0.0;
  Long64_t nAll  = 0; // T->GetEntries(all)
  Long64_t nBase = 0; // T->GetEntries(baseCuts)
  Long64_t nFull = 0; // T->GetEntries(fullCuts)
  Long64_t nWide = 0; // within wide coincidence gate
  Long64_t nCoin = 0; // within coin window (subset of wide)
};

static YieldEx AddYieldFromRun(TTree* T,
                              const TString& cutsExpr,
                              const TString& weightExpr,
                              const CoincidenceConfig& cfg,
                              const CoincidenceResult& win,
                              const TString& baseCutsExpr)
{
  static long long uid = 0;
  TString hname = Form("h_yield_%lld", uid++);

  TH1D h(hname, hname, 1, 0, 1);
  h.SetDirectory(nullptr);
  h.Sumw2(true);

  // Counts (useful for later filters / diagnostics)
  YieldEx y;
  y.nAll  = T->GetEntries();
  y.nBase = T->GetEntries(baseCutsExpr);
  y.nFull = T->GetEntries(cutsExpr);

  TString wideGate = BuildRangeCut(cfg.CtBranchName, cfg.WideWindowMinNs, cfg.WideWindowMaxNs);
  TString cutsWide = CombineCutsAND(cutsExpr, wideGate);
  y.nWide = T->GetEntries(cutsWide);

  TString coinGate = BuildRangeCut(cfg.CtBranchName, win.CoinWindowNs.first, win.CoinWindowNs.second);
  TString cutsCoin = CombineCutsAND(cutsWide, coinGate);
  y.nCoin = T->GetEntries(cutsCoin);

  Bool_t oldAdd = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kTRUE);
  FillRandomSubtractedHistogramWithWindows(T, cutsExpr, weightExpr, "0.5", &h, cfg, win);
  TH1::AddDirectory(oldAdd);

  y.sumw  = h.GetBinContent(1);
  const double e = h.GetBinError(1);
  y.sumw2 = e*e;
  return y;
}

struct SimYieldEx {
  double sumw  = 0.0;
  double sumw2 = 0.0;
  Long64_t nSel = 0;
};

static SimYieldEx SimYieldSum(TChain& sim,
                              const TString& exprWeight,
                              const TString& cutBool)
{
  SimYieldEx out;
  Long64_t n = sim.Draw(exprWeight, cutBool, "goff");
  out.nSel = n;
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
  out.sumw2 = sum2;
  return out;
}

// ---------- flags + note ----------

enum FlagBits : unsigned int {
  MC_EMPTY        = 1u << 0,
  MC_STARVED      = 1u << 1,
  MC_BAD_RELERR   = 1u << 2,
  NET_NEGATIVE    = 1u << 3,
  NET_NANINF      = 1u << 4,
  SIM_NANINF      = 1u << 5,
  XSEC_NANINF     = 1u << 6,
  XSEC_NEGATIVE   = 1u << 7,
  XSEC_BAD_RELERR = 1u << 8,
  XSEC_ZERO_ERR   = 1u << 9,
  SIM_ZERO        = 1u << 10,
  NET_ZERO        = 1u << 11
};

static void AppendNote(std::string& note, const char* token) {
  if (!note.empty()) note.push_back(';');
  note += token;
}

static double SafeRelErr(double v, double e) {
  if (!std::isfinite(v) || !std::isfinite(e)) return kMissing;
  if (v == 0.0) return kMissing;
  return std::abs(e / v);
}

// ---------- group file reader ----------

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

// ---------- table header ----------

static std::string TableHeader() {
  // Keep stable ordering: plotters can rely on it.
  return
    "schema_version,mode,group_id,curve_label,setting_id,"
    "pt_bin,z_bin,phi_bin,"
    "pt_lo,pt_hi,pt_center,z_lo,z_hi,z_center,phi_lo,phi_hi,phi_center,"
    "n_sig,n_pos,n_dum,n_posdum,n_net,n_sim,"
    "y_sig,y_pos,y_dum,y_posdum,y_net,y_net_err,y_sim,y_sim_err,"
    "sim_scale,sigma_model,mean_xB,mean_Q2,nominal_z,"
    "xsec,xsec_err,rel_net_err,rel_sim_err,rel_xsec_err,"
    "flag_bits,valid_default,note";
}

// ---------- core computation for one bin ----------

struct BinOut {
  // data
  YieldEx sig, pos, dum, posdum;
  double yNet = 0.0;
  double eNet = 0.0;
  Long64_t nNet = 0;

  // sim
  SimYieldEx sim;
  double ySim = 0.0;
  double eSim = 0.0;

  // xsec
  double xsec = kMissing;
  double xsecErr = kMissing;

  // derived
  double relNet = kMissing;
  double relSim = kMissing;
  double relXsec = kMissing;

  unsigned int flags = 0;
  int valid_default = 0;
  std::string note;
};

// Guard thresholds for "valid_default" (table keeps everything regardless)
struct GuardCfg {
  double minSimYield   = 0.05; // same default as group macro
  double maxRelSimErr  = 0.35;
  double maxRelXsecErr = 0.50;
};

static void FinalizeFlagsAndValidity(BinOut& o, const GuardCfg& gc) {
  // Notes reflect flags
  if (o.flags & MC_EMPTY)        AppendNote(o.note, "MC_EMPTY");
  if (o.flags & MC_STARVED)      AppendNote(o.note, "MC_STARVED");
  if (o.flags & MC_BAD_RELERR)   AppendNote(o.note, "MC_BAD_RELERR");
  if (o.flags & SIM_ZERO)        AppendNote(o.note, "SIM_ZERO");
  if (o.flags & SIM_NANINF)      AppendNote(o.note, "SIM_NANINF");
  if (o.flags & NET_ZERO)        AppendNote(o.note, "NET_ZERO");
  if (o.flags & NET_NEGATIVE)    AppendNote(o.note, "NET_NEGATIVE");
  if (o.flags & NET_NANINF)      AppendNote(o.note, "NET_NANINF");
  if (o.flags & XSEC_NEGATIVE)   AppendNote(o.note, "XSEC_NEGATIVE");
  if (o.flags & XSEC_ZERO_ERR)   AppendNote(o.note, "XSEC_ZERO_ERR");
  if (o.flags & XSEC_BAD_RELERR) AppendNote(o.note, "XSEC_BAD_RELERR");
  if (o.flags & XSEC_NANINF)     AppendNote(o.note, "XSEC_NANINF");

  const unsigned int catastrophic = (MC_EMPTY | NET_NANINF | SIM_NANINF | XSEC_NANINF);
  const unsigned int defaultReject = catastrophic | MC_STARVED | MC_BAD_RELERR | XSEC_BAD_RELERR | SIM_ZERO;
  o.valid_default = ((o.flags & defaultReject) == 0u) ? 1 : 0;
}

// Compute one bin (DATA yields + SIM yield + xsec), but do not decide plotting.
static BinOut ComputeOneBin(
  // data context
  const std::vector<int>& runs_sig,
  const std::vector<int>& runs_pos,
  const std::vector<int>& runs_dum,
  const std::vector<int>& runs_posdum,
  const std::unordered_map<int, RunNormInfo>& runInfo,
  std::unordered_map<int, CoincidenceResult>& coinCache,
  const CoincidenceConfig& cfg,
  const TCut& dataCuts,
  const TCut& kinData,
  double dumScale,
  // sim context
  TChain& sim,
  const TCut& simCutsBase,
  const TCut& kinSim,
  const TString& simWeightExpr,
  // physics
  double sigmaModel,
  const GuardCfg& guards,
  std::ostream& log,
  bool doSanityOnce,
  bool& sanityPrinted
)
{
  BinOut o;

  auto CategoryYield = [&](const std::vector<int>& runs,
                           const TCut& baseCutsBool,
                           const TCut& kinCutsBool,
                           YieldEx& out,
                           double& totalCharge_mC) -> void
  {
    out = YieldEx{};
    totalCharge_mC = 0.0;
    int nOpened = 0;

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

      // Compute coincidence windows once per run
      auto itC = coinCache.find(run);
      if (itC == coinCache.end()) {
        Bool_t oldAdd = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kTRUE);
        CoincidenceResult win = ComputeCoincidenceRandomSubtraction(T, baseCutsStr, cfg);
        TH1::AddDirectory(oldAdd);
        coinCache[run] = win;
        itC = coinCache.find(run);
      }

      // Optional sanity print (first opened file only)
      if (doSanityOnce && !sanityPrinted) {
        Long64_t nAll  = T->GetEntries();
        Long64_t nBase = T->GetEntries(baseCutsStr);
        Long64_t nFull = T->GetEntries(fullCutsStr);
        log << "SANITY (run " << run << "): entries(all)=" << nAll
            << " entries(baseCuts)=" << nBase
            << " entries(fullCuts)=" << nFull << "\n";
        sanityPrinted = true;
      }

      TString wexpr = Form("%g", w);
      YieldEx y = AddYieldFromRun(T, fullCutsStr, wexpr, cfg, itC->second, baseCutsStr);

      out.sumw  += y.sumw;
      out.sumw2 += y.sumw2;
      out.nAll  += y.nAll;
      out.nBase += y.nBase;
      out.nFull += y.nFull;
      out.nWide += y.nWide;
      out.nCoin += y.nCoin;

      f->Close();
    }

    if (nOpened == 0) {
      log << "ERROR: opened 0 ROOT files for this category (yields will be zero).\n";
    }

    if (totalCharge_mC > 0) {
      out.sumw  /= totalCharge_mC;
      out.sumw2 /= (totalCharge_mC * totalCharge_mC);
    }
  };

  double qSig=0, qPos=0, qDum=0, qPosDum=0;
  bool localSanity = doSanityOnce;

  CategoryYield(runs_sig,    dataCuts, kinData, o.sig,    qSig);
  CategoryYield(runs_pos,    dataCuts, kinData, o.pos,    qPos);
  CategoryYield(runs_dum,    dataCuts, kinData, o.dum,    qDum);
  CategoryYield(runs_posdum, dataCuts, kinData, o.posdum, qPosDum);

  // net yield: (Data − PosData) − dumScale*(Dummy − PosDummy)
  o.yNet = (o.sig.sumw - o.pos.sumw) - dumScale * (o.dum.sumw - o.posdum.sumw);
  const double eNet2 = o.sig.sumw2 + o.pos.sumw2 + dumScale*dumScale * (o.dum.sumw2 + o.posdum.sumw2);
  o.eNet = (eNet2 > 0) ? std::sqrt(eNet2) : 0.0;
  o.nNet = o.sig.nCoin; // heuristic: signal coin count as a proxy

  if (!std::isfinite(o.yNet) || !std::isfinite(o.eNet)) o.flags |= NET_NANINF;
  if (o.yNet < 0) o.flags |= NET_NEGATIVE;
  if (o.yNet == 0.0) o.flags |= NET_ZERO;
  o.relNet = SafeRelErr(o.yNet, o.eNet);

  // SIM
  TString simBool = TString((simCutsBase && kinSim).GetTitle());
  o.sim = SimYieldSum(sim, simWeightExpr, simBool);
  o.ySim = o.sim.sumw;
  o.eSim = (o.sim.sumw2 > 0) ? std::sqrt(o.sim.sumw2) : 0.0;
  if (!std::isfinite(o.ySim) || !std::isfinite(o.eSim)) o.flags |= SIM_NANINF;
  if (o.sim.nSel <= 0) o.flags |= MC_EMPTY;
  if (!(o.ySim > 0.0)) o.flags |= SIM_ZERO;

  o.relSim = SafeRelErr(o.ySim, o.eSim);

  // Health/starvation checks (do NOT delete points; only flag)
  bool mc_ok = (o.ySim > 0.0) && std::isfinite(o.ySim) && std::isfinite(o.eSim);
  bool mc_starved = true;
  if (mc_ok) {
    // starved if below min yield or too-large relative uncertainty
    const double rel = (o.relSim != kMissing) ? o.relSim : 1e9;
    mc_starved = (o.ySim < guards.minSimYield) || (rel > guards.maxRelSimErr);
    if (rel > guards.maxRelSimErr) o.flags |= MC_BAD_RELERR;
  }
  if (mc_starved) o.flags |= MC_STARVED;

  // xsec
  if (mc_ok) {
    o.xsec = (o.yNet / o.ySim) * sigmaModel;
    double rel2 = 0.0;
    if (o.yNet != 0.0) rel2 += (o.eNet / o.yNet) * (o.eNet / o.yNet);
    rel2 += (o.eSim / o.ySim) * (o.eSim / o.ySim);
    o.xsecErr = std::abs(o.xsec) * std::sqrt(std::abs(rel2));
  } else {
    o.xsec = kMissing;
    o.xsecErr = kMissing;
  }

  if (!std::isfinite(o.xsec) || !std::isfinite(o.xsecErr)) {
    o.flags |= XSEC_NANINF;
    o.xsec = kMissing;
    o.xsecErr = kMissing;
  }
  if (o.xsec != kMissing && o.xsec < 0.0) o.flags |= XSEC_NEGATIVE;
  if (o.xsec != kMissing && o.xsecErr == 0.0 && o.xsec != 0.0) o.flags |= XSEC_ZERO_ERR;

  o.relXsec = (o.xsec == kMissing) ? kMissing : SafeRelErr(o.xsec, o.xsecErr);
  if (o.relXsec != kMissing && o.relXsec > guards.maxRelXsecErr) o.flags |= XSEC_BAD_RELERR;

  FinalizeFlagsAndValidity(o, guards);
  return o;
}

// ---------- SINGLE mode table ----------

static void TableCoinXsec_Single(const char* manifestPath,
                                 const char* resultsRoot,
                                 const char* settingsRoot)
{
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  const std::string manifest = manifestPath ? std::string(manifestPath) : "";
  const std::string leafDir  = NormalizeSlashes(Dirname(manifest));
  const std::string settingId = MakeSettingIdFromManifestPath(manifest, settingsRoot);
  const std::string resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + settingId);

  MkdirP(resultsDir);
  MkdirP(resultsDir + "/tables");

  std::ofstream log(resultsDir + "/log_table.txt", std::ios::out);
  log << "=== TableCoinXsec v5 (SINGLE) ===\n";
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

  // Nominal kin from settingId
  double xNom, qNom, zNom;
  ExtractNominalKinFromSettingId(settingId, xNom, qNom, zNom);
  double meanXB = ChargeWeightedMean(runs_sig, runInfo, true);
  double meanQ2 = ChargeWeightedMean(runs_sig, runInfo, false);
  if (!std::isfinite(meanXB)) meanXB = xNom;
  if (!std::isfinite(meanQ2)) meanQ2 = qNom;

  // DATA base cuts (same as PlotCoinXsec_Single)
  const TCut dataCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_gtr_th>-0.15) && (H_gtr_th<0.15) && "
    "(H_gtr_ph>-0.10) && (H_gtr_ph<0.10) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_gtr_th>-0.15) && (P_gtr_th<0.15) && "
    "(P_gtr_ph>-0.10) && (P_gtr_ph<0.10) && "
    "(P_cal_etottracknorm<0.8)";

  const std::string data_pt  = "pt";
  const std::string data_z   = "z";
  const std::string data_phi =
    "((P_kin_secondary_ph_xq>=0)*P_kin_secondary_ph_xq + (P_kin_secondary_ph_xq<0)*(P_kin_secondary_ph_xq + 2*TMath::Pi()))";

  const double wall_thickness_ratio = 3.82;
  const double dumScale = 1.0 / wall_thickness_ratio;

  // SIM vars
  const std::string sim_pt  = "sqrt(pt2)";
  const std::string sim_z   = "z";
  const std::string sim_phi = "((phipq>=0)*phipq + (phipq<0)*(phipq + 2*TMath::Pi()))";

  // SIM chain
  TChain sim("h10");
  const std::string projectRoot = NormalizeSlashes(Dirname(NormalizeSlashes(settingsRoot)));
  std::string simGlob = projectRoot + "/Pass0_SimROOTfiles/" + settingId + "/*.root";
  int nAdded = sim.Add(simGlob.c_str());
  log << "SIM glob: " << simGlob << " added=" << nAdded << "\n";
  const Long64_t nGenSim = sim.GetEntries();
  if (nGenSim <= 0) {
    log << "ERROR: SIM chain empty\n";
    std::cerr << "ERROR: SIM chain empty for: " << simGlob << "\n";
    return;
  }

  const TCut simCutsBase =
    "(hsdelta>-8.0) && (hsdelta<8.0) && "
    "(hsxptar>-0.15) && (hsxptar<0.15) && "
    "(hsyptar>-0.10) && (hsyptar<0.10) && "
    "(ssdelta>-10.0) && (ssdelta<22.0) && "
    "(ssxptar>-0.15) && (ssxptar<0.15) && "
    "(ssyptar>-0.10) && (ssyptar<0.10)";

  // --- SIM normalization: MUST come from the corresponding SIMC .hist file ---
  // Directory layout assumption (your case):
  //   <projectRoot>/Pass0_SimROOTfiles/<settingId>/<stem>.root
  //   <projectRoot>/Pass0_SimROOTfiles/<settingId>/<stem>.hist
  const std::string simRootPath = GetFirstRootFileFromChain(sim, log, simGlob);
  std::string simHistPath;
  const double simNormalizationFactor = GrabNormfacFromHist(simRootPath, &simHistPath, log);
  log << "SIM root file: " << simRootPath << "\n";
  log << "SIM hist file: " << simHistPath << "\n";
  log << "SIM normfac:  " << simNormalizationFactor << "\n";

  const double simTryScale = (nGenSim > 0) ? (1.0 / double(nGenSim)) : 0.0;
  const double simEffectiveScale = simNormalizationFactor * simTryScale;
  log << "SIM nGenSim:  " << nGenSim << "\n";
  log << "SIM scale:   " << simTryScale << "  (1/nGenSim)\n";
  log << "SIM eff:     " << simEffectiveScale << "  (normfac/nGenSim)\n";

  std::cout << "SIM root file: " << simRootPath << "\n";
  std::cout << "SIM hist file: " << simHistPath << "\n";
  std::cout << "SIM normfac:  " << simNormalizationFactor << "\n";
  std::cout << "SIM nGenSim:  " << nGenSim << "\n";
  std::cout << "SIM eff:     " << simEffectiveScale << "  (normfac/nGenSim)\n";

  const TString simWeightExpr = Form("Weight*%e*%e", simNormalizationFactor, simTryScale);

  // Random subtraction config
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  // Binning
  Binning binning;
  const int nPt  = (int)binning.ptEdges.size()  - 1;
  const int nZ   = (int)binning.zEdges.size()   - 1;
  const int nPhi = (int)binning.phiEdges.size() - 1;

  // Guards
  GuardCfg guards;

  // Coin window cache per run
  std::unordered_map<int, CoincidenceResult> coinCache;

  const std::string outCsvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_single.csv";
  std::ofstream csv(outCsvPath, std::ios::out);
  csv << TableHeader() << "\n";
  log << "CSV: " << outCsvPath << "\n";

  bool sanityPrinted = false;
  const bool doSanityOnce = true;

  const int schema_version = 1;
  const std::string mode = "single";
  const std::string group_id = settingId;
  const std::string curve_label = "";

  for (int iPt=0; iPt<nPt; iPt++) {
    const double ptLo = binning.ptEdges[iPt];
    const double ptHi = binning.ptEdges[iPt+1];
    const double ptCtr= Center(ptLo, ptHi);

    for (int iZ=0; iZ<nZ; iZ++) {
      const double zLo = binning.zEdges[iZ];
      const double zHi = binning.zEdges[iZ+1];
      const double zCtr= Center(zLo, zHi);

      for (int iPhi=0; iPhi<nPhi; iPhi++) {
        const double phiLo = binning.phiEdges[iPhi];
        const double phiHi = binning.phiEdges[iPhi+1];
        const double phiCtr= Center(phiLo, phiHi);

        const TCut kinData = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                                  data_pt.c_str(), ptLo, data_pt.c_str(), ptHi,
                                  data_z.c_str(),  zLo,  data_z.c_str(),  zHi,
                                  data_phi.c_str(),phiLo, data_phi.c_str(),phiHi);

        const TCut kinSim = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                                 sim_pt.c_str(), ptLo, sim_pt.c_str(), ptHi,
                                 sim_z.c_str(),  zLo,  sim_z.c_str(),  zHi,
                                 sim_phi.c_str(),phiLo, sim_phi.c_str(),phiHi);

        BinOut o = ComputeOneBin(
          runs_sig, runs_posSig, runs_dum, runs_posDum,
          runInfo, coinCache, cfg,
          dataCuts, kinData, dumScale,
          sim, simCutsBase, kinSim, simWeightExpr,
          sigmaModel, guards, log,
          doSanityOnce, sanityPrinted);

        // Counts per category (coin counts as a useful proxy)
        const Long64_t nSig = o.sig.nCoin;
        const Long64_t nPos = o.pos.nCoin;
        const Long64_t nDum = o.dum.nCoin;
        const Long64_t nPosD= o.posdum.nCoin;
        const Long64_t nNet = o.nNet;
        const Long64_t nSim = o.sim.nSel;

        // sim_scale recorded for reproducibility
        const double sim_scale = simTryScale;

        csv
          << schema_version << "," << mode << "," << group_id << "," << curve_label << "," << settingId << ","
          << iPt << "," << iZ << "," << iPhi << ","
          << ptLo << "," << ptHi << "," << ptCtr << ","
          << zLo  << "," << zHi  << "," << zCtr  << ","
          << phiLo << "," << phiHi << "," << phiCtr << ","
          << nSig << "," << nPos << "," << nDum << "," << nPosD << "," << nNet << "," << nSim << ","
          << o.sig.sumw << "," << o.pos.sumw << "," << o.dum.sumw << "," << o.posdum.sumw << ","
          << o.yNet << "," << o.eNet << "," << o.ySim << "," << o.eSim << ","
          << sim_scale << "," << sigmaModel << "," << meanXB << "," << meanQ2 << "," << zNom << ","
          << o.xsec << "," << o.xsecErr << "," << o.relNet << "," << o.relSim << "," << o.relXsec << ","
          << o.flags << "," << o.valid_default << "," << o.note
          << "\n";
      }
    }
  }

  log << "Done (SINGLE table).\n";
  log.close();
  csv.close();
}

// ---------- GROUP mode table ----------

static void TableCoinXsec_Group(const char* groupPath,
                                const char* resultsRoot,
                                const char* settingsRoot)
{
  gROOT->SetBatch(kTRUE);
  TH1::AddDirectory(kFALSE);

  const std::string groupFile = groupPath ? std::string(groupPath) : "";
  if (groupFile.empty()) {
    std::cerr << "ERROR: groupPath is empty\n";
    return;
  }

  const std::string groupId = StripExtension(groupFile);
  const std::string resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + groupId);
  MkdirP(resultsDir);
  MkdirP(resultsDir + "/tables");

  std::ofstream log(resultsDir + "/log_table.txt", std::ios::out);
  log << "=== TableCoinXsec v5 (GROUP overlay) ===\n";
  log << "groupFile: " << groupFile << "\n";
  log << "resultsDir: " << resultsDir << "\n";
  log << "settingsRoot: " << settingsRoot << "\n";

  auto entries = ReadGroupListFile(groupFile, settingsRoot, log);
  if (entries.empty()) {
    log << "ERROR: group file produced 0 entries.\n";
    std::cerr << "ERROR: group file produced 0 entries.\n";
    return;
  }
  log << "Group entries: " << entries.size() << "\n";

  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  // Binning (match PlotCoinXsec_Group defaults)
  Binning binning;
  const bool overlayWidePt = false;
  if (overlayWidePt) binning.ptEdges = {0.0, 0.40};
  binning.zEdges = {0.0, 1.0};

  const int nPt  = (int)binning.ptEdges.size()  - 1;
  const int nZ   = (int)binning.zEdges.size()   - 1; // should be 1
  const int nPhi = (int)binning.phiEdges.size() - 1;

  // DATA base cuts (match PlotCoinXsec_Group)
  const TCut dataCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_cal_etottracknorm<0.8)";

  const std::string data_pt  = "pt";
  const std::string data_z   = "z";
  const std::string data_phi =
    "((P_kin_secondary_ph_xq>=0)*P_kin_secondary_ph_xq + (P_kin_secondary_ph_xq<0)*(P_kin_secondary_ph_xq + 2*TMath::Pi()))";

  const double wall_thickness_ratio = 3.82;
  const double dumScale = 1.0 / wall_thickness_ratio;

  const std::string sim_pt  = "sqrt(pt2)";
  const std::string sim_z   = "z";
  const std::string sim_phi = "((phipq>=0)*phipq + (phipq<0)*(phipq + 2*TMath::Pi()))";

  const TCut simCutsBase =
    "(hsdelta>-8.0) && (hsdelta<8.0) && "
    "(ssdelta>-10.0) && (ssdelta<22.0)";

  GuardCfg guards;
  // Use the same defaults as the plotting macro (can be tuned later in plotter)
  guards.minSimYield   = 0.05;
  guards.maxRelSimErr  = 0.35;
  guards.maxRelXsecErr = 0.50;

  const std::string outCsvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_group.csv";
  std::ofstream csv(outCsvPath, std::ios::out);
  csv << TableHeader() << "\n";
  log << "CSV: " << outCsvPath << "\n";

  const int schema_version = 1;
  const std::string mode = "group";

  // For each curve (setting), load its inputs and write all bins
  for (const auto& e : entries) {
    const std::string settingId = e.settingId;
    const std::string leafDir   = e.leafDir;
    const std::string curve_label = e.label;

    log << "--- Curve: " << curve_label << " settingId=" << settingId << "\n";

    // Run lists
    auto runs_sig    = ReadRunList(leafDir + "/runs_signal.txt", log);
    auto runs_dum    = ReadRunList(leafDir + "/runs_dummy.txt", log);
    auto runs_pos    = ReadRunList(leafDir + "/runs_positron.txt", log);
    auto runs_posdum = ReadRunList(leafDir + "/runs_positron_dummy.txt", log);

    std::unordered_map<int, RunNormInfo> runInfo;
    if (!LoadRunMetadataCSV(leafDir + "/run_metadata.csv", runInfo)) {
      log << "ERROR: cannot load run_metadata.csv for " << settingId << "\n";
      std::cerr << "ERROR: cannot load run_metadata.csv for " << settingId << "\n";
      return;
    }

    const double sigmaModel = ReadSigmaModel(leafDir, log);

    // Nominal/mean kin
    double xNom, qNom, zNom;
    ExtractNominalKinFromSettingId(settingId, xNom, qNom, zNom);
    double meanXB = ChargeWeightedMean(runs_sig, runInfo, true);
    double meanQ2 = ChargeWeightedMean(runs_sig, runInfo, false);
    if (!std::isfinite(meanXB)) meanXB = xNom;
    if (!std::isfinite(meanQ2)) meanQ2 = qNom;

    // SIM chain
    std::unique_ptr<TChain> sim(new TChain("h10"));
    const std::string projectRoot = NormalizeSlashes(Dirname(NormalizeSlashes(settingsRoot)));
    std::string simGlob = projectRoot + "/Pass0_SimROOTfiles/" + settingId + "/*.root";
    int nAdded = sim->Add(simGlob.c_str());
    log << "SIM glob: " << simGlob << " added=" << nAdded << "\n";
    const Long64_t nGenSim = sim->GetEntries();
    if (nGenSim <= 0) {
      log << "ERROR: SIM chain empty for: " << simGlob << "\n";
      std::cerr << "ERROR: SIM chain empty for: " << simGlob << "\n";
      return;
    }

    // --- SIM normalization: MUST come from the corresponding SIMC .hist file ---
    const std::string simRootPath = GetFirstRootFileFromChain(*sim, log, simGlob);
    std::string simHistPath;
    const double simNormalizationFactor = GrabNormfacFromHist(simRootPath, &simHistPath, log);
    log << "SIM root file: " << simRootPath << "\n";
    log << "SIM hist file: " << simHistPath << "\n";
    log << "SIM normfac:  " << simNormalizationFactor << "\n";

    const double simTryScale = (nGenSim > 0) ? (1.0 / double(nGenSim)) : 0.0;
    const double simEffectiveScale = simNormalizationFactor * simTryScale;
    log << "SIM nGenSim:  " << nGenSim << "\n";
    log << "SIM scale:   " << simTryScale << "  (1/nGenSim)\n";
    log << "SIM eff:     " << simEffectiveScale << "  (normfac/nGenSim)\n";

    std::cout << "[GROUP] " << curve_label << " SIM root file: " << simRootPath << "\n";
    std::cout << "[GROUP] " << curve_label << " SIM hist file: " << simHistPath << "\n";
    std::cout << "[GROUP] " << curve_label << " SIM normfac:  " << simNormalizationFactor << "\n";
    std::cout << "[GROUP] " << curve_label << " SIM nGenSim:  " << nGenSim << "\n";
    std::cout << "[GROUP] " << curve_label << " SIM eff:     " << simEffectiveScale << "  (normfac/nGenSim)\n";

    const TString simWeightExpr = Form("Weight*%e*%e", simNormalizationFactor, simTryScale);

    // Coin cache for this curve
    std::unordered_map<int, CoincidenceResult> coinCache;
    bool sanityPrinted = false;
    const bool doSanityOnce = true;

    for (int iPt=0; iPt<nPt; iPt++) {
      const double ptLo = binning.ptEdges[iPt];
      const double ptHi = binning.ptEdges[iPt+1];
      const double ptCtr= Center(ptLo, ptHi);

      for (int iZ=0; iZ<nZ; iZ++) {
        const double zLo = binning.zEdges[iZ];
        const double zHi = binning.zEdges[iZ+1];
        const double zCtr= Center(zLo, zHi);

        for (int iPhi=0; iPhi<nPhi; iPhi++) {
          const double phiLo = binning.phiEdges[iPhi];
          const double phiHi = binning.phiEdges[iPhi+1];
          const double phiCtr= Center(phiLo, phiHi);

          // Group mode integrates over z (zEdges = {0,1}), but keep the formal dimension.
          const TCut kinData = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                                    data_pt.c_str(), ptLo, data_pt.c_str(), ptHi,
                                    data_z.c_str(),  zLo,  data_z.c_str(),  zHi,
                                    data_phi.c_str(),phiLo, data_phi.c_str(),phiHi);

          const TCut kinSim = Form("(%s>=%g && %s<%g) && (%s>=%g && %s<%g) && (%s>=%g && %s<%g)",
                                   sim_pt.c_str(), ptLo, sim_pt.c_str(), ptHi,
                                   sim_z.c_str(),  zLo,  sim_z.c_str(),  zHi,
                                   sim_phi.c_str(),phiLo, sim_phi.c_str(),phiHi);

          BinOut o = ComputeOneBin(
            runs_sig, runs_pos, runs_dum, runs_posdum,
            runInfo, coinCache, cfg,
            dataCuts, kinData, dumScale,
            *sim, simCutsBase, kinSim, simWeightExpr,
            sigmaModel, guards, log,
            doSanityOnce, sanityPrinted);

          const Long64_t nSig = o.sig.nCoin;
          const Long64_t nPos = o.pos.nCoin;
          const Long64_t nDum = o.dum.nCoin;
          const Long64_t nPosD= o.posdum.nCoin;
          const Long64_t nNet = o.nNet;
          const Long64_t nSim = o.sim.nSel;

          const double sim_scale = simTryScale;

          csv
            << schema_version << "," << mode << "," << groupId << "," << curve_label << "," << settingId << ","
            << iPt << "," << iZ << "," << iPhi << ","
            << ptLo << "," << ptHi << "," << ptCtr << ","
            << zLo  << "," << zHi  << "," << zCtr  << ","
            << phiLo << "," << phiHi << "," << phiCtr << ","
            << nSig << "," << nPos << "," << nDum << "," << nPosD << "," << nNet << "," << nSim << ","
            << o.sig.sumw << "," << o.pos.sumw << "," << o.dum.sumw << "," << o.posdum.sumw << ","
            << o.yNet << "," << o.eNet << "," << o.ySim << "," << o.eSim << ","
            << sim_scale << "," << sigmaModel << "," << meanXB << "," << meanQ2 << "," << zNom << ","
            << o.xsec << "," << o.xsecErr << "," << o.relNet << "," << o.relSim << "," << o.relXsec << ","
            << o.flags << "," << o.valid_default << "," << o.note
            << "\n";
        }
      }
    }
  }

  log << "Done (GROUP table).\n";
  log.close();
  csv.close();
}

} // namespace

// Wrapper: mirrors PlotCoinXsec call signature (sans PNG/table toggles)
void TableCoinXsec(const char* manifestPath,
                   const char* resultsRoot = "results",
                   const char* settingsRoot = "settings")
{
  if (!manifestPath || std::string(manifestPath).empty()) {
    std::cerr << "ERROR: manifestPath is required\n";
    return;
  }
  std::string p = manifestPath;
  if (EndsWith(p, ".list")) {
    TableCoinXsec_Group(manifestPath, resultsRoot, settingsRoot);
  } else {
    TableCoinXsec_Single(manifestPath, resultsRoot, settingsRoot);
  }
}
