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
// Method:
//   raw_ratio r = yield / elclean_per_mC
//   r_ref = weighted mean of lowest N=2 rate points
//   scaled_ratio y = r / r_ref  (so y ~ 1 at low rate; "y-intercept = 1" normalization)
//   Fit scaled_ratio with intercept fixed: y(x)=1+m x, x in kHz
//   tau[s] = -m/1000 ; tau[ns] = tau*1e9

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

static std::string RelUnderSettings(const std::string& settingDirAbs) {
  const std::string needle = "/settings/";
  auto pos = settingDirAbs.find(needle);
  if (pos == std::string::npos) return BaseName(settingDirAbs);
  std::string rel = settingDirAbs.substr(pos + needle.size());
  while (!rel.empty() && rel.front()=='/') rel.erase(0,1);
  return rel.empty() ? BaseName(settingDirAbs) : rel;
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
  double raw_ratio=NAN, raw_ratio_err=NAN;
  double scaled_ratio=NAN, scaled_ratio_err=NAN;
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
  const int N_LOW = 2; // user-selected

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

  auto warn = [&](int run, const std::string& msg){
    std::cerr << "WARNING [run " << run << "]: " << msg << "\n";
    log       << "WARNING [run " << run << "]: " << msg << "\n";
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

  // Candidates for reference selection
  struct Cand { int run; double rate; double r; double rerr; };
  std::vector<Cand> cands;

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

    cands.push_back({run, o.rate_kHz, o.raw_ratio, o.raw_ratio_err});
  }

  // Select lowest N_LOW by rate
  std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b){ return a.rate < b.rate; });

  if ((int)cands.size() < N_LOW) {
    warn(-1, "Not enough valid runs to define reference (need N_LOW=2). Will write CSV without fit.");
  }

  std::vector<Cand> refPts;
  for (int i=0;i<std::min(N_LOW,(int)cands.size());++i) refPts.push_back(cands[i]);

  // Weighted mean r_ref and its error
  double sumW=0, sumWR=0;
  for (const auto& c : refPts) {
    double w = 1.0/(c.rerr*c.rerr);
    sumW  += w;
    sumWR += w*c.r;
  }
  double r_ref = (sumW>0) ? (sumWR/sumW) : NAN;
  double r_ref_err = (sumW>0) ? std::sqrt(1.0/sumW) : NAN;

  log << "Reference (lowest N=" << N_LOW << " by rate):\n";
  for (size_t i=0;i<refPts.size();++i) {
    log << "  REF" << i+1 << " run=" << refPts[i].run
        << " rate_kHz=" << refPts[i].rate
        << " raw_ratio=" << refPts[i].r << " +/- " << refPts[i].rerr << "\n";
  }
  log << "  r_ref=" << std::setprecision(12) << r_ref << " +/- " << r_ref_err << "\n\n";

  // Build plot/fit vectors using scaled ratio
  std::vector<double> x, y, ex, ey;
  int nUsed=0;

  for (auto& o : rows) {
    if (o.join_flag != "OK") { o.use_in_fit=0; continue; }
    if (!std::isfinite(r_ref) || r_ref<=0 || !std::isfinite(r_ref_err)) {
      o.use_in_fit=0;
      o.join_flag="INVALID_VALUES";
      o.reason="invalid r_ref";
      continue;
    }

    o.scaled_ratio = o.raw_ratio / r_ref;

    // include ref uncertainty
    // (σy/y)^2 = (σr/r)^2 + (σref/ref)^2
    double term1 = o.raw_ratio_err / o.raw_ratio;
    double term2 = r_ref_err / r_ref;
    o.scaled_ratio_err = std::fabs(o.scaled_ratio) * std::sqrt(term1*term1 + term2*term2);

    if (!std::isfinite(o.scaled_ratio) || !std::isfinite(o.scaled_ratio_err) || o.scaled_ratio_err<=0) {
      o.use_in_fit=0;
      o.join_flag="INVALID_VALUES";
      o.reason="invalid scaled_ratio or scaled_ratio_err";
      continue;
    }

    o.use_in_fit=1;
    x.push_back(o.rate_kHz);
    y.push_back(o.scaled_ratio);
    ex.push_back(0.0);
    ey.push_back(o.scaled_ratio_err);
    nUsed++;
  }

  // Sort by x for plotting
  SortByX(x, y, ex, ey);

  // Output paths
  const std::string outCsv      = outTabs + "/yield_ratio_vs_trigger_shms34.csv";
  const std::string outPng      = outPNGs + "/yield_ratio_vs_trigger_shms34.png";
  const std::string outRootFile = outCanv + "/yield_ratio_vs_trigger_shms34.root";

  // ---------- Fit + plot (style matched to YieldVsCurrent.C) ----------
  double fit_m=NAN, fit_merr=NAN, chi2=NAN, prob=NAN;
  int ndf=-1;
  double tau_ns=NAN, tau_ns_err=NAN;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TGaxis::SetMaxDigits(3);

  // Compute ranges
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
  const double ymin = minY - ypad;
  const double ymax = maxY + ypad;

  TCanvas c("c_yield_ratio_vs_trigger", "Yield Ratio vs SHMS34 (signal)", 1200, 850);
  c.SetTopMargin(0.22);     // match YieldVsCurrent.C
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

  // Graph (red, like other v1 plots)
  TGraphErrors g((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.55);
  g.SetMarkerColor(kRed+1);
  g.SetLineColor(kRed+1);
  g.Draw("P SAME");

  // y=1 reference line (dashed)
  TLine y1(xmin, 1.0, xmax, 1.0);
  y1.SetLineColor(kBlack);
  y1.SetLineStyle(2);
  y1.SetLineWidth(2);
  y1.Draw("SAME");

  // Fit (if enough points)
  TF1* f = nullptr;
  if ((int)x.size() >= 2) {
    // Fit only on the data x-range (not padded)
    const double fxmin = x.front();
    const double fxmax = x.back();

    f = new TF1("f_ratio", "1 + [0]*x", fxmin, fxmax);
    f->SetParameter(0, -1e-4);
    f->SetLineColor(kBlack);
    f->SetLineStyle(1);
    f->SetLineWidth(2);

    g.Fit(f, "Q"); // quiet

    fit_m    = f->GetParameter(0);
    fit_merr = f->GetParError(0);
    chi2 = f->GetChisquare();
    ndf  = f->GetNDF();
    prob = (ndf>0) ? TMath::Prob(chi2, ndf) : NAN;

    // tau conversion
    double tau_s     = -fit_m / 1000.0;
    double tau_s_err =  fit_merr / 1000.0;
    tau_ns     = tau_s * 1e9;
    tau_ns_err = tau_s_err * 1e9;

    f->Draw("SAME");
  } else {
    warn(-1, "Not enough points to fit (need >=2).");
  }

  // Legend (same placement/style)
  TLegend leg(0.10, 0.82, 0.50, 0.98);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.035);
  leg.AddEntry(&g, "signal", "p");
  leg.AddEntry(&y1, "y = 1", "l");
  if (f) leg.AddEntry(f, "fit: y = 1 + m x", "l");
  leg.Draw();

  // Setting label (same placement/style)
  std::string lbl1, lbl2;
  BuildSettingLabel(rel, lbl1, lbl2);

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(13);
  t.SetTextSize(0.028);
  t.DrawLatex(0.39, 0.94, lbl1.c_str());
  t.DrawLatex(0.39, 0.89, lbl2.c_str());

  // Middle text (blue), like your other plots
  if (f && std::isfinite(chi2) && ndf > 0) {
    TLatex tf;
    tf.SetNDC();
    tf.SetTextAlign(22);        // centered
    tf.SetTextColor(kBlue+2);
    tf.SetTextSize(0.040);
    tf.DrawLatex(0.60, 0.52, Form("#chi^{2}/ndf = %.2f (%d), Prob = %.3g", chi2/ndf, ndf, prob));
    tf.SetTextSize(0.033);
    tf.DrawLatex(0.60, 0.47, Form("m = %.4g #pm %.2g (per kHz)    #tau = %.3g #pm %.3g ns",
                                  fit_m, fit_merr, tau_ns, tau_ns_err));
  }

  // Save outputs
  c.SaveAs(outPng.c_str());
  c.SaveAs(outRootFile.c_str());

  // ---------- Write CSV (with REF diagnostics) ----------
  {
    std::ofstream csv(outCsv.c_str());
    if (!csv.is_open()) {
      std::cerr << "ERROR: cannot write " << outCsv << "\n";
      return;
    }

    csv << "# YieldRatioVsTrigger\n";
    csv << "# manifest: " << manifestP << "\n";
    csv << "# rel_under_settings: " << rel << "\n";
    csv << "# inputs:\n";
    csv << "#   " << inYieldCSV << "\n";
    csv << "#   " << inElCSV << "\n";
    csv << "# reference_lowN: " << N_LOW << "\n";
    csv << std::setprecision(12);
    csv << "# r_ref: " << r_ref << "\n";
    csv << "# r_ref_err: " << r_ref_err << "\n";
    csv << "# N_total: " << nTotal << "\n";
    csv << "# N_used_in_fit: " << nUsed << "\n";

    // --- Requested diagnostic lines: exactly which 2 runs formed the reference ---
    csv << "#REF,ref_index,run,rate_kHz,raw_ratio,raw_ratio_err\n";
    for (size_t i=0;i<refPts.size();++i) {
      csv << "#REF," << (i+1) << ","
          << refPts[i].run << ","
          << refPts[i].rate << ","
          << refPts[i].r << ","
          << refPts[i].rerr << "\n";
    }

    if (f && ndf > 0) {
      csv << "# fit_model: y = 1 + m * x   (x in kHz)\n";
      csv << "# m_per_kHz: " << fit_m << "\n";
      csv << "# m_per_kHz_err: " << fit_merr << "\n";
      csv << "# tau_ns: " << tau_ns << "\n";
      csv << "# tau_ns_err: " << tau_ns_err << "\n";
      csv << "# chi2: " << chi2 << "\n";
      csv << "# ndf: " << ndf << "\n";
      csv << "# prob: " << prob << "\n";
    }
    csv << "\n";

    csv << "run,shms34_rate_kHz,"
        << "yield,yield_err,"
        << "hms_elclean_counts_per_mC,hms_elclean_counts_per_mC_err,"
        << "raw_ratio,raw_ratio_err,"
        << "scaled_ratio,scaled_ratio_err,"
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
          << r.scaled_ratio << ","
          << r.scaled_ratio_err << ","
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

  if (f) delete f;
}
