// macros/YieldRatioVsTrigger.C
//
// Usage (from rate_dependence_v1/):
//   root -l -b -q 'macros/YieldRatioVsTrigger.C("settings/.../manifest.json","results")'
//
// Reads (from results/<rel_under_settings>/tables/):
//   yield_vs_trigger_shms34.csv
//   hms_elclean_chargeNorm_vs_current.csv
//
// Writes (to results/<rel_under_settings>/):
//   tables/yield_ratio_vs_trigger_shms34.csv
//   PNGs/yield_ratio_vs_trigger_shms34.png
//   canvases/yield_ratio_vs_trigger_shms34.root
//   logs/YieldRatioVsTrigger.log
//
// Method (your requested "fit intercept free, then normalize"):
//   1) Build raw_ratio r = yield / elclean_per_mC
//   2) Fit raw_ratio vs rate with free intercept: r(R) = a + b R   (weighted least squares)
//   3) Define normalized ratio: y_norm(R) = r(R) / a   so y_norm(0)=1
//      Propagate error including a uncertainty:
//        (σ_y/y)^2 = (σ_r/r)^2 + (σ_a/a)^2
//   4) Fit normalized ratio with intercept fixed: y_norm(R) = 1 + m R  (weighted)
//   5) Convert m (per kHz) -> tau:
//        tau [s]  = -m/1000
//        tau [ns] = tau*1e9
//
// Notes:
// - signal only, status==OK only
// - no outlier rejection (C1)
// - full valid range (A1)
// - warnings:
//     warn if prob < 0.01
//     strong warn if prob < 1e-4

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TGaxis.h"

// ---------- helpers ----------
static std::string Trim(const std::string& s) {
  size_t b=0,e=s.size();
  while (b<e && std::isspace((unsigned char)s[b])) b++;
  while (e>b && std::isspace((unsigned char)s[e-1])) e--;
  return s.substr(b,e-b);
}

// CSV split (quotes supported, "" escapes)
static std::vector<std::string> SplitCSVLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq=false;
  for (size_t i=0;i<line.size();++i) {
    char c=line[i];
    if (c=='"') {
      if (inq && i+1<line.size() && line[i+1]=='"') { cur.push_back('"'); i++; }
      else inq=!inq;
    } else if (c==',' && !inq) {
      out.push_back(Trim(cur));
      cur.clear();
    } else cur.push_back(c);
  }
  out.push_back(Trim(cur));
  return out;
}

static double SafeToD(const std::string& s) {
  try {
    size_t pos=0;
    std::string t=Trim(s);
    double v=std::stod(t,&pos);
    if (pos!=t.size()) return NAN;
    return v;
  } catch (...) { return NAN; }
}
static int SafeToI(const std::string& s, bool& ok) {
  ok=false;
  try {
    size_t pos=0;
    std::string t=Trim(s);
    int v=std::stoi(t,&pos);
    if (pos!=t.size()) return -1;
    ok=true;
    return v;
  } catch (...) { return -1; }
}

static void EnsureDir(const std::string& p) { gSystem->mkdir(p.c_str(), kTRUE); }
static std::string DirName(const std::string& p) { return std::string(gSystem->DirName(p.c_str())); }
static std::string BaseName(const std::string& p) { return std::string(gSystem->BaseName(p.c_str())); }

// Robust rel path under settings: works for absolute "/.../settings/..." and relative "settings/..."
static std::string RelUnderSettings(const std::string& settingDir) {
  // try absolute form first
  const std::string needleAbs = "/settings/";
  auto pos = settingDir.find(needleAbs);
  if (pos != std::string::npos) {
    std::string rel = settingDir.substr(pos + needleAbs.size());
    while (!rel.empty() && rel.front()=='/') rel.erase(0,1);
    return rel.empty() ? BaseName(settingDir) : rel;
  }
  // try relative form
  const std::string needleRel = "settings/";
  pos = settingDir.find(needleRel);
  if (pos != std::string::npos) {
    std::string rel = settingDir.substr(pos + needleRel.size());
    while (!rel.empty() && rel.front()=='/') rel.erase(0,1);
    return rel.empty() ? BaseName(settingDir) : rel;
  }
  return BaseName(settingDir);
}

static std::vector<std::string> SplitPath(const std::string& s) {
  std::vector<std::string> out; std::string cur;
  for (char c : s) {
    if (c=='/') { if (!cur.empty()) out.push_back(cur); cur.clear(); }
    else cur.push_back(c);
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

static void BuildSettingLabel(const std::string& relUnderSettings,
                              std::string& line1, std::string& line2)
{
  auto tok = SplitPath(relUnderSettings);
  if (tok.size() >= 7) {
    line1 = tok[0] + " / " + tok[1] + " / " + tok[2] + " / " +
            tok[3] + " / " + tok[4] + " / " + tok[5];
    line2 = tok[6];
  } else {
    line1 = relUnderSettings;
    line2 = "";
  }
}

static void SortByX(std::vector<double>& x, std::vector<double>& y,
                    std::vector<double>& ex, std::vector<double>& ey) {
  if (x.size() <= 1) return;
  std::vector<size_t> idx(x.size());
  for (size_t i=0;i<idx.size();++i) idx[i]=i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return x[a] < x[b]; });
  auto reorder = [&](std::vector<double>& v){
    std::vector<double> tmp; tmp.reserve(v.size());
    for (auto i: idx) tmp.push_back(v[i]);
    v.swap(tmp);
  };
  reorder(x); reorder(y); reorder(ex); reorder(ey);
}

// ---------- structs ----------
struct YieldRow { bool present=false; std::string status; double yield=NAN, yield_err=NAN, rate_kHz=NAN; };
struct ElRow    { bool present=false; std::string status; double el=NAN, el_err=NAN; };

struct OutRow {
  int run=-1;
  double rate_kHz=NAN;
  double yield=NAN, yield_err=NAN;
  double el=NAN, el_err=NAN;

  // raw ratio
  double raw_ratio=NAN, raw_ratio_err=NAN;

  // normalized ratio (raw/a)
  double ynorm=NAN, ynorm_err=NAN;

  std::string status_yield="MISSING";
  std::string status_el="MISSING";
  std::string join_flag="";
  int use_in_fit=0;
  std::string reason="";
};

static bool ReadYieldCSV(const std::string& path, std::map<int,YieldRow>& out, std::ostream& log) {
  std::ifstream f(path);
  if (!f.is_open()) { log << "ERROR: cannot open " << path << "\n"; return false; }
  std::string header; if (!std::getline(f, header)) { log << "ERROR: empty " << path << "\n"; return false; }

  auto cols = SplitCSVLine(header);
  std::unordered_map<std::string,int> idx;
  for (int i=0;i<(int)cols.size();++i) idx[cols[i]] = i;

  const char* req[] = {"category","run","yield","yield_err","shms34_rate_kHz","status"};
  for (auto* name : req) {
    if (idx.find(name) == idx.end()) {
      log << "ERROR: missing column '" << name << "' in " << path << "\n";
      return false;
    }
  }

  auto get = [&](const std::vector<std::string>& r, const std::string& name)->std::string{
    auto it = idx.find(name);
    if (it==idx.end()) return "";
    int j = it->second;
    if (j<0 || j>=(int)r.size()) return "";
    return r[j];
  };

  std::string line;
  int nAll=0, nSig=0;
  while (std::getline(f, line)) {
    if (Trim(line).empty()) continue;
    nAll++;
    auto r = SplitCSVLine(line);
    if (Trim(get(r,"category")) != "signal") continue;
    nSig++;

    bool ok=false;
    int run = SafeToI(get(r,"run"), ok);
    if (!ok) continue;

    YieldRow y;
    y.present=true;
    y.status = Trim(get(r,"status"));
    y.yield = SafeToD(get(r,"yield"));
    y.yield_err = SafeToD(get(r,"yield_err"));
    y.rate_kHz = SafeToD(get(r,"shms34_rate_kHz"));
    out[run]=y;
  }

  log << "Loaded " << path << ": rows=" << nAll << ", signal=" << nSig << ", runs=" << out.size() << "\n";
  return true;
}

static bool ReadElCSV(const std::string& path, std::map<int,ElRow>& out, std::ostream& log) {
  std::ifstream f(path);
  if (!f.is_open()) { log << "ERROR: cannot open " << path << "\n"; return false; }
  std::string header; if (!std::getline(f, header)) { log << "ERROR: empty " << path << "\n"; return false; }

  auto cols = SplitCSVLine(header);
  std::unordered_map<std::string,int> idx;
  for (int i=0;i<(int)cols.size();++i) idx[cols[i]] = i;

  const char* req[] = {"category","run","hms_elclean_counts_per_mC","hms_elclean_counts_per_mC_err","status"};
  for (auto* name : req) {
    if (idx.find(name) == idx.end()) {
      log << "ERROR: missing column '" << name << "' in " << path << "\n";
      return false;
    }
  }

  auto get = [&](const std::vector<std::string>& r, const std::string& name)->std::string{
    auto it = idx.find(name);
    if (it==idx.end()) return "";
    int j = it->second;
    if (j<0 || j>=(int)r.size()) return "";
    return r[j];
  };

  std::string line;
  int nAll=0, nSig=0;
  while (std::getline(f, line)) {
    if (Trim(line).empty()) continue;
    nAll++;
    auto r = SplitCSVLine(line);
    if (Trim(get(r,"category")) != "signal") continue;
    nSig++;

    bool ok=false;
    int run = SafeToI(get(r,"run"), ok);
    if (!ok) continue;

    ElRow e;
    e.present=true;
    e.status = Trim(get(r,"status"));
    e.el = SafeToD(get(r,"hms_elclean_counts_per_mC"));
    e.el_err = SafeToD(get(r,"hms_elclean_counts_per_mC_err"));
    out[run]=e;
  }

  log << "Loaded " << path << ": rows=" << nAll << ", signal=" << nSig << ", runs=" << out.size() << "\n";
  return true;
}

// ---------- main ----------
void YieldRatioVsTrigger(const char* manifestPath, const char* resultsDir)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsDir   ? resultsDir   : "";
  if (manifestP.empty() || outRoot.empty()) {
    std::cerr << "ERROR: Usage: YieldRatioVsTrigger(\".../manifest.json\",\"results\")\n";
    return;
  }

  const std::string settingDir = DirName(manifestP);
  const std::string rel = RelUnderSettings(settingDir);
  const std::string outBase = outRoot + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTabs = outBase + "/tables";
  const std::string outCanv = outBase + "/canvases";
  const std::string outLogs = outBase + "/logs";
  EnsureDir(outPNGs);
  EnsureDir(outTabs);
  EnsureDir(outCanv);
  EnsureDir(outLogs);

  const std::string logPath = outLogs + "/YieldRatioVsTrigger.log";
  std::ofstream log(logPath.c_str());

  auto warn = [&](const std::string& msg){
    std::cerr << "WARNING: " << msg << "\n";
    log       << "WARNING: " << msg << "\n";
  };
  auto strong_warn = [&](const std::string& msg){
    std::cerr << "WARNING (STRONG): " << msg << "\n";
    log       << "WARNING (STRONG): " << msg << "\n";
  };

  log << "YieldRatioVsTrigger (v1, signal-only)\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "results_base: " << outBase << "\n";
  log << "log: " << logPath << "\n\n";

  const std::string inYieldCSV = outTabs + "/yield_vs_trigger_shms34.csv";
  const std::string inElCSV    = outTabs + "/hms_elclean_chargeNorm_vs_current.csv";

  std::map<int,YieldRow> ymap;
  std::map<int,ElRow> emap;

  if (!ReadYieldCSV(inYieldCSV, ymap, log)) return;
  if (!ReadElCSV(inElCSV, emap, log)) return;

  // Collect union of runs
  std::set<int> allRuns;
  for (const auto& kv : ymap) allRuns.insert(kv.first);
  for (const auto& kv : emap) allRuns.insert(kv.first);

  std::vector<OutRow> rows;
  rows.reserve(allRuns.size());

  // First pass: build raw ratio candidates for fit1
  std::vector<double> x_raw, r_raw, ex_raw, er_raw;

  int nTotal=0;
  for (int run : allRuns) {
    nTotal++;
    OutRow o; o.run = run;

    auto itY = ymap.find(run);
    auto itE = emap.find(run);

    bool hasY = (itY!=ymap.end() && itY->second.present);
    bool hasE = (itE!=emap.end() && itE->second.present);

    if (!hasY) {
      o.join_flag="MISSING_YIELD";
      o.reason="missing in yield_vs_trigger_shms34.csv";
      rows.push_back(o);
      continue;
    }
    if (!hasE) {
      o.join_flag="MISSING_ELCLEAN";
      o.reason="missing in hms_elclean_chargeNorm_vs_current.csv";
      o.status_yield = itY->second.status;
      o.rate_kHz = itY->second.rate_kHz;
      o.yield = itY->second.yield;
      o.yield_err = itY->second.yield_err;
      rows.push_back(o);
      continue;
    }

    // Fill basics
    o.status_yield = itY->second.status;
    o.status_el    = itE->second.status;
    o.rate_kHz     = itY->second.rate_kHz;
    o.yield        = itY->second.yield;
    o.yield_err    = itY->second.yield_err;
    o.el           = itE->second.el;
    o.el_err       = itE->second.el_err;

    if (Trim(o.status_yield)!="OK" || Trim(o.status_el)!="OK") {
      o.join_flag="BAD_STATUS";
      std::ostringstream ss; ss << "status_yield=" << o.status_yield << ", status_el=" << o.status_el;
      o.reason = ss.str();
      rows.push_back(o);
      continue;
    }

    // Validate
    auto finitePos = [&](double v){ return std::isfinite(v) && v>0; };
    if (!std::isfinite(o.rate_kHz)) { o.join_flag="INVALID_VALUES"; o.reason="non-finite rate"; rows.push_back(o); continue; }
    if (!finitePos(o.yield) || !finitePos(o.el) || !finitePos(o.yield_err) || !finitePos(o.el_err)) {
      o.join_flag="INVALID_VALUES";
      o.reason="non-finite or non-positive yield/el/errs";
      rows.push_back(o);
      continue;
    }

    // Raw ratio and error
    o.raw_ratio = o.yield / o.el;
    double relY = o.yield_err / o.yield;
    double relE = o.el_err / o.el;
    o.raw_ratio_err = std::fabs(o.raw_ratio) * std::sqrt(relY*relY + relE*relE);

    if (!std::isfinite(o.raw_ratio) || !std::isfinite(o.raw_ratio_err) || o.raw_ratio_err<=0) {
      o.join_flag="INVALID_VALUES";
      o.reason="invalid raw_ratio or raw_ratio_err";
      rows.push_back(o);
      continue;
    }

    o.join_flag="OK";
    rows.push_back(o);

    x_raw.push_back(o.rate_kHz);
    r_raw.push_back(o.raw_ratio);
    ex_raw.push_back(0.0);
    er_raw.push_back(o.raw_ratio_err);
  }

  // Sort by x for fit stability / plotting
  SortByX(x_raw, r_raw, ex_raw, er_raw);

  // Output paths
  const std::string outCsv      = outTabs + "/yield_ratio_vs_trigger_shms34.csv";
  const std::string outPng      = outPNGs + "/yield_ratio_vs_trigger_shms34.png";
  const std::string outRootFile = outCanv + "/yield_ratio_vs_trigger_shms34.root";

  // Fit results
  // Fit1: raw ratio r = a + bR
  double a=NAN, aerr=NAN, b=NAN, berr=NAN, chi2_1=NAN, prob_1=NAN;
  int ndf_1=-1;

  // Fit2: normalized y = 1 + mR
  double m=NAN, merr=NAN, chi2_2=NAN, prob_2=NAN;
  int ndf_2=-1;
  double tau_ns=NAN, tau_ns_err=NAN;

  // Need at least 2 points for fit1
  if ((int)x_raw.size() < 2) {
    warn("Not enough valid points for fit1 (need >=2). Will write CSV and produce an empty plot.");
  } else {
    TGraphErrors gRaw((int)x_raw.size(), x_raw.data(), r_raw.data(), ex_raw.data(), er_raw.data());
    const double fxmin = x_raw.front();
    const double fxmax = x_raw.back();

    TF1 f1("f_raw_pol1", "pol1", fxmin, fxmax);
    // reasonable starting values
    f1.SetParameter(0, r_raw.front());
    f1.SetParameter(1, 0.0);

    // weighted fit (uses y-errors)
    gRaw.Fit(&f1, "Q"); // quiet

    a    = f1.GetParameter(0);
    aerr = f1.GetParError(0);
    b    = f1.GetParameter(1);
    berr = f1.GetParError(1);
    chi2_1 = f1.GetChisquare();
    ndf_1  = f1.GetNDF();
    prob_1 = (ndf_1 > 0) ? TMath::Prob(chi2_1, ndf_1) : NAN;

    log << "Fit1 (raw): r = a + b*R\n";
    log << "  range_kHz: [" << fxmin << ", " << fxmax << "]\n";
    log << std::setprecision(12);
    log << "  a = " << a << " +/- " << aerr << "\n";
    log << "  b = " << b << " +/- " << berr << " (per kHz)\n";
    if (ndf_1 > 0) log << "  chi2/ndf = " << (chi2_1/ndf_1) << " (" << ndf_1 << "), prob=" << prob_1 << "\n\n";

    if (std::isfinite(prob_1)) {
      if (prob_1 < 1e-4) strong_warn("Fit1 probability is very small (<1e-4). Linear model may be insufficient or data inconsistent.");
      else if (prob_1 < 0.01) warn("Fit1 probability is small (<0.01). Interpret a (intercept) with caution.");
    }

    if (!std::isfinite(a) || a == 0.0 || !std::isfinite(aerr) || aerr <= 0.0) {
      strong_warn("Fit1 returned invalid a or aerr; cannot normalize. Will write CSV/empty plot.");
    }
  }

  // Second pass: compute ynorm, ynorm_err and build fit2 vectors (only if a is valid)
  std::vector<double> x, y, ex, ey;
  int nUsed2=0;

  const bool canNorm = (std::isfinite(a) && std::isfinite(aerr) && a != 0.0 && aerr > 0.0);
  for (auto& o : rows) {
    if (o.join_flag != "OK") { o.use_in_fit=0; continue; }
    if (!canNorm) {
      o.use_in_fit=0;
      o.reason = "cannot normalize (invalid a/aerr)";
      continue;
    }
    // y_norm = raw/a
    o.ynorm = o.raw_ratio / a;

    // include a uncertainty:
    // (σy/y)^2 = (σr/r)^2 + (σa/a)^2
    const double term1 = o.raw_ratio_err / o.raw_ratio;
    const double term2 = aerr / a;
    o.ynorm_err = std::fabs(o.ynorm) * std::sqrt(term1*term1 + term2*term2);

    if (!std::isfinite(o.ynorm) || !std::isfinite(o.ynorm_err) || o.ynorm_err <= 0) {
      o.use_in_fit=0;
      o.join_flag="INVALID_VALUES";
      o.reason="invalid ynorm or ynorm_err";
      continue;
    }

    o.use_in_fit=1;
    x.push_back(o.rate_kHz);
    y.push_back(o.ynorm);
    ex.push_back(0.0);
    ey.push_back(o.ynorm_err);
    nUsed2++;
  }

  SortByX(x, y, ex, ey);

  // ---------- Fit2 + plot (style matched to YieldVsCurrent.C) ----------
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(3);

  // Compute ranges for normalized plot
  bool any = (!x.empty());
  double minX=0,maxX=0,minY=0,maxY=0;
  if (any) {
    minX=maxX=x[0];
    minY=y[0]-ey[0];
    maxY=y[0]+ey[0];
    for (size_t i=0;i<x.size();++i) {
      minX = std::min(minX, x[i]);
      maxX = std::max(maxX, x[i]);
      minY = std::min(minY, y[i]-ey[i]);
      maxY = std::max(maxY, y[i]+ey[i]);
    }
  }

  const double dx = (any ? (maxX-minX) : 1.0);
  const double dy = (any ? (maxY-minY) : 1.0);
  const double xpad = 0.06 * (dx>0?dx:1.0);
  const double ypad = 0.10 * (dy>0?dy:1.0);

  const double xmin = minX - xpad;
  const double xmax = maxX + xpad;
  double ymin = minY - ypad;
  double ymax = maxY + ypad;
  if (ymax < 1.02) ymax = 1.02;

  TCanvas c("c_yield_ratio_vs_trigger", "Yield Ratio vs SHMS34 (signal)", 1200, 850);
  c.SetTopMargin(0.22);
  c.SetRightMargin(0.06);
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);

  TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("SHMS 3/4 Trigger Rate (kHz)");
  frame->GetYaxis()->SetTitle("Yield Ratio");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);

  // Normalized graph (red)
  TGraphErrors g((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.55);
  g.SetMarkerColor(kRed+1);
  g.SetLineColor(kRed+1);
  g.Draw("P SAME");

  // y=1 reference line
  const double xlo = frame->GetXaxis()->GetXmin();
  const double xhi = frame->GetXaxis()->GetXmax();
  TLine y1(xlo, 1.0, xhi, 1.0);
  y1.SetLineColor(kBlack);
  y1.SetLineStyle(2);
  y1.SetLineWidth(2);
  y1.Draw("SAME");

  // Fit2 if enough points
  TF1* f2 = nullptr;
  if ((int)x.size() >= 2) {
    const double fxmin = x.front();
    const double fxmax = x.back();

    f2 = new TF1("f_norm", "1 + [0]*x", fxmin, fxmax);
    f2->SetParameter(0, -1e-4);
    f2->SetLineColor(kBlack);
    f2->SetLineStyle(1);
    f2->SetLineWidth(2);

    g.Fit(f2, "Q");

    m    = f2->GetParameter(0);
    merr = f2->GetParError(0);
    chi2_2 = f2->GetChisquare();
    ndf_2  = f2->GetNDF();
    prob_2 = (ndf_2 > 0) ? TMath::Prob(chi2_2, ndf_2) : NAN;

    // tau conversion
    double tau_s     = -m / 1000.0;
    double tau_s_err =  merr / 1000.0;
    tau_ns     = tau_s * 1e9;
    tau_ns_err = tau_s_err * 1e9;

    f2->Draw("SAME");

    log << "Fit2 (norm): y = 1 + m*R\n";
    log << "  range_kHz: [" << fxmin << ", " << fxmax << "]\n";
    log << std::setprecision(12);
    log << "  m = " << m << " +/- " << merr << " (per kHz)\n";
    log << "  tau_ns = " << tau_ns << " +/- " << tau_ns_err << "\n";
    if (ndf_2 > 0) log << "  chi2/ndf = " << (chi2_2/ndf_2) << " (" << ndf_2 << "), prob=" << prob_2 << "\n\n";

    if (std::isfinite(prob_2)) {
      if (prob_2 < 1e-4) strong_warn("Fit2 probability is very small (<1e-4). Linear deadtime approx may be insufficient or data inconsistent.");
      else if (prob_2 < 0.01) warn("Fit2 probability is small (<0.01). Interpret slope/tau with caution.");
    }

  } else {
    warn("Not enough points to do Fit2 (need >=2).");
  }

  // Legend
  TLegend leg(0.10, 0.82, 0.50, 0.98);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.035);
  leg.AddEntry(&g, "signal", "p");
  leg.AddEntry(&y1, "y = 1", "l");
  if (f2) leg.AddEntry(f2, "fit: y = 1 + m x", "l");
  leg.Draw();

  // Setting label
  std::string lbl1, lbl2;
  BuildSettingLabel(rel, lbl1, lbl2);

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(13);
  t.SetTextSize(0.028);
  t.DrawLatex(0.39, 0.94, lbl1.c_str());
  t.DrawLatex(0.39, 0.89, lbl2.c_str());

  // Middle text (blue) — show Fit2 plus show 'a' from Fit1 (so you see how far off intercept was)
  if (f2 && std::isfinite(chi2_2) && ndf_2 > 0) {
    TLatex tf;
    tf.SetNDC();
    tf.SetTextAlign(22);
    tf.SetTextColor(kBlue+2);
    tf.SetTextSize(0.040);
    tf.DrawLatex(0.60, 0.52, Form("#chi^{2}/ndf = %.2f (%d), Prob = %.3g", chi2_2/ndf_2, ndf_2, prob_2));
    tf.SetTextSize(0.033);
    tf.DrawLatex(0.60, 0.47, Form("m = %.4g #pm %.2g (per kHz)    #tau = %.3g #pm %.3g ns",
                                  m, merr, tau_ns, tau_ns_err));
    // smaller line for fit1 intercept
    if (std::isfinite(a) && std::isfinite(aerr)) {
      tf.SetTextSize(0.028);
      tf.DrawLatex(0.60, 0.42, Form("Fit1 raw intercept a = %.4g #pm %.2g", a, aerr));
    }
  }

  // Save outputs
  c.SaveAs(outPng.c_str());
  c.SaveAs(outRootFile.c_str());

  // ---------- Write CSV ----------
  {
    std::ofstream csv(outCsv.c_str());
    if (!csv.is_open()) {
      std::cerr << "ERROR: cannot write " << outCsv << "\n";
      if (f2) delete f2;
      return;
    }

    csv << "# YieldRatioVsTrigger\n";
    csv << "# manifest: " << manifestP << "\n";
    csv << "# rel_under_settings: " << rel << "\n";
    csv << "# inputs:\n";
    csv << "#   " << inYieldCSV << "\n";
    csv << "#   " << inElCSV << "\n";
    csv << "# N_total: " << nTotal << "\n";
    csv << "# N_valid_for_fit1: " << x_raw.size() << "\n";
    csv << "# N_valid_for_fit2: " << x.size() << "\n";
    csv << std::setprecision(12);

    // Fit1 summary
    if (std::isfinite(a) && std::isfinite(b) && ndf_1 > 0) {
      csv << "# fit1_model: raw_ratio = a + b*x   (x in kHz)\n";
      csv << "# fit1_a: " << a << "\n";
      csv << "# fit1_a_err: " << aerr << "\n";
      csv << "# fit1_b_per_kHz: " << b << "\n";
      csv << "# fit1_b_err_per_kHz: " << berr << "\n";
      csv << "# fit1_chi2: " << chi2_1 << "\n";
      csv << "# fit1_ndf: " << ndf_1 << "\n";
      csv << "# fit1_prob: " << prob_1 << "\n";
    } else {
      csv << "# fit1_model: raw_ratio = a + b*x   (x in kHz)\n";
      csv << "# fit1_status: FAILED_OR_INSUFFICIENT_POINTS\n";
    }

    // Fit2 summary
    if (f2 && ndf_2 > 0) {
      csv << "# fit2_model: y_norm = 1 + m*x   (x in kHz)\n";
      csv << "# fit2_m_per_kHz: " << m << "\n";
      csv << "# fit2_m_err_per_kHz: " << merr << "\n";
      csv << "# tau_ns: " << tau_ns << "\n";
      csv << "# tau_ns_err: " << tau_ns_err << "\n";
      csv << "# fit2_chi2: " << chi2_2 << "\n";
      csv << "# fit2_ndf: " << ndf_2 << "\n";
      csv << "# fit2_prob: " << prob_2 << "\n";
    } else {
      csv << "# fit2_model: y_norm = 1 + m*x   (x in kHz)\n";
      csv << "# fit2_status: FAILED_OR_INSUFFICIENT_POINTS\n";
    }

    csv << "\n";

    // Columns (include raw columns + ynorm columns)
    csv << "run,shms34_rate_kHz,"
        << "yield,yield_err,"
        << "hms_elclean_counts_per_mC,hms_elclean_counts_per_mC_err,"
        << "raw_ratio,raw_ratio_err,"
        << "ynorm,ynorm_err,"
        << "status_yield_csv,status_elclean_csv,join_flag,use_in_fit,reason\n";

    csv << std::setprecision(10);
    for (const auto& r : rows) {
      csv << r.run << ","
          << r.rate_kHz << ","
          << r.yield << ","
          << r.yield_err << ","
          << r.el << ","
          << r.el_err << ","
          << r.raw_ratio << ","
          << r.raw_ratio_err << ","
          << r.ynorm << ","
          << r.ynorm_err << ","
          << r.status_yield << ","
          << r.status_el << ","
          << r.join_flag << ","
          << r.use_in_fit << ","
          << "\"" << r.reason << "\"\n";
    }
  }

  log << "Wrote CSV: " << outCsv << "\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT: " << outRootFile << "\n";
  log.close();

  if (f2) delete f2;
}
