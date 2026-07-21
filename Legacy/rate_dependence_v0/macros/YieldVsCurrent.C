#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TString.h"
#include "TPad.h"
#include "TH1.h"
#include "TLatex.h"


#include "../include/CoincidenceRandomSubtraction.h"

// ---------- helpers ----------
static double SafeToD(const std::string& s) {
  try { return std::stod(s); } catch (...) { return NAN; }
}
static bool IsBadSentinel(double v) {
  if (!std::isfinite(v)) return true;
  // treat -999 as "bad"
  return (std::fabs(v + 999.0) < 1e-9);
}

static std::vector<std::string> SplitCSVLine(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool inq=false;
  for (size_t i=0;i<s.size();++i) {
    char c=s[i];
    if (c=='"') { inq=!inq; continue; }
    if (c==',' && !inq) { out.push_back(cur); cur.clear(); }
    else cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

static void EnsureDir(const std::string& p) { gSystem->mkdir(p.c_str(), kTRUE); }
static std::string DirName(const std::string& p) { return std::string(gSystem->DirName(p.c_str())); }
static std::string BaseName(const std::string& p) { return std::string(gSystem->BaseName(p.c_str())); }

static std::string RelUnderSettings(const std::string& settingDirAbs) {
  // Expect: .../settings/<pass>/<run_type>/<target>/<z>/<x>/<Q2>/<setting_id>
  const std::string needle = "/settings/";
  auto pos = settingDirAbs.find(needle);
  if (pos == std::string::npos) return BaseName(settingDirAbs); // fallback
  std::string rel = settingDirAbs.substr(pos + needle.size());
  while (!rel.empty() && rel.front()=='/') rel.erase(0,1);
  return rel.empty() ? BaseName(settingDirAbs) : rel;
}


static std::vector<std::string> SplitPath(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c=='/') {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else cur.push_back(c);
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

static void BuildSettingLabel(const std::string& relUnderSettings,
                              std::string& line1, std::string& line2)
{
  // relUnderSettings expected:
  // pass4/pi+sidis/LH2/z0p36/x0p25/Q23p3/hmsPneg...
  auto tok = SplitPath(relUnderSettings);

  if (tok.size() >= 7) {
    // line 1: pass/run_type/target + z x Q
    line1 = tok[0] + " / " + tok[1] + " / " + tok[2] + " / " +
            tok[3] + " / " + tok[4] + " / " + tok[5];
    // line 2: setting_id
    line2 = tok[6];
  } else {
    // fallback
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

static std::string FormatWindows(const std::vector<std::pair<double,double>>& wins) {
  // "lo1-hi1;lo2-hi2;..."
  std::ostringstream ss;
  for (size_t i=0;i<wins.size();++i) {
    if (i) ss << ";";
    ss << wins[i].first << "-" << wins[i].second;
  }
  return ss.str();
}

// ---------- metadata row ----------
struct MetaRow {
  std::string category;
  int run = -1;

  double BCM2_I = NAN;
  double BCM2_Q = NAN;
  double comp_livetime = NAN;
  double h_esing_Eff   = NAN;
  double boil_corr     = NAN;

  double ps_factor     = NAN;
  std::string ps_choice;

  std::string target_raw;
  std::string run_type_raw;
};

static std::vector<MetaRow> ReadRunMetadata(const std::string& path, std::ostream& log) {
  std::vector<MetaRow> rows;
  std::ifstream f(path);
  if (!f.is_open()) {
    log << "ERROR: cannot open " << path << "\n";
    return rows;
  }

  std::string header;
  if (!std::getline(f, header)) return rows;

  auto cols = SplitCSVLine(header);
  std::unordered_map<std::string,int> idx;
  for (int i=0;i<(int)cols.size();++i) idx[cols[i]] = i;

  auto col = [&](const std::vector<std::string>& r, const std::string& name)->std::string{
    auto it = idx.find(name);
    if (it==idx.end()) return "";
    int j = it->second;
    if (j < 0 || j >= (int)r.size()) return "";
    return r[j];
  };

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto r = SplitCSVLine(line);

    MetaRow m;
    m.category = col(r, "category");
    try { m.run = std::stoi(col(r, "run")); } catch (...) { continue; }

    m.BCM2_I        = SafeToD(col(r, "BCM2_I"));
    m.BCM2_Q        = SafeToD(col(r, "BCM2_Q"));
    m.comp_livetime = SafeToD(col(r, "comp_livetime"));
    m.h_esing_Eff   = SafeToD(col(r, "h_esing_Eff"));
    m.boil_corr     = SafeToD(col(r, "boil_corr"));

    m.ps_factor     = SafeToD(col(r, "ps_factor"));
    m.ps_choice     = col(r, "ps_choice");

    m.target_raw    = col(r, "target_raw");
    m.run_type_raw  = col(r, "run_type_raw");

    rows.push_back(m);
  }

  log << "Loaded rows from run_metadata.csv: " << rows.size() << "\n";
  return rows;
}

// ---------- main ----------
void YieldVsCurrent(const char* manifestPath,
                    const char* resultsDir)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsDir   ? resultsDir   : "";

  const std::string settingDir = DirName(manifestP);
  const std::string setting_id = BaseName(settingDir);

  // Mirror everything after /settings/
  const std::string rel = RelUnderSettings(settingDir);
  const std::string outBase = outRoot + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTabs = outBase + "/tables";
  const std::string outCanv = outBase + "/canvases";
  EnsureDir(outPNGs);
  EnsureDir(outTabs);
  EnsureDir(outCanv);

  std::ofstream log((outBase + "/log.txt").c_str());
  log << "YieldVsCurrent (myCTime only)\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "results_base: " << outBase << "\n\n";

  // Read metadata
  const std::string metaPath = settingDir + "/run_metadata.csv";
  auto meta = ReadRunMetadata(metaPath, log);
  if (meta.empty()) {
    log << "ERROR: no metadata rows. Nothing to do.\n";
    return;
  }

  // Physics cuts (as given)
  const TString baseCuts =
      "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
      "(H_cal_etottracknorm>0.7) && "
      "(H_cer_npeSum>2.0) && "
      "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
      "(P_cal_etottracknorm<0.8)";

  // Coincidence random subtraction config
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";
  // keep defaults for windows unless you later want to tune:
  // cfg.WideWindowMinNs, cfg.WideWindowMaxNs, cfg.RfPeriodNs, cfg.PeakHalfWidthNs, cfg.MaxSidePeaks ...

  // Output CSV
  const std::string outCsv = outTabs + "/yield_vs_current_overlay.csv";
  std::ofstream csv(outCsv.c_str());
  csv << "category,run,BCM2_I,"
      << "Nnet,Nnet_err,"
      << "norm_factor,yield,yield_err,"
      << "BCM2_Q_mC,comp_livetime,h_esing_Eff,boil_corr,ps_choice,ps_factor,"
      << "PeakCenterNs,CoinLoNs,CoinHiNs,RandomWindowListNs,rootfile,status\n";

  // For plotting: per category series
  struct Series { std::vector<double> x,y,ex,ey; };
  std::map<std::string, Series> series;

  auto warn = [&](int run, const std::string& msg){
    std::cerr << "WARNING [run " << run << "]: " << msg << "\n";
    log       << "WARNING [run " << run << "]: " << msg << "\n";
  };

  // ROOT file location (as you confirmed)
  const std::string skimDir = "./Skimmed_ROOTfiles";

  // Loop runs
  for (const auto& m : meta) {
    const int run = m.run;
    const std::string rootfile = Form("%s/skimmed_coin_replay_production_%d_-1.root", skimDir.c_str(), run);

    std::string status = "OK";

    // Validate normalization scalars (flag -999 etc.)
    bool badNorm = false;

    if (IsBadSentinel(m.BCM2_I)) { warn(run, "BCM2_I is invalid (-999/NaN)."); badNorm = true; }
    if (IsBadSentinel(m.BCM2_Q) || m.BCM2_Q <= 0) {
      warn(run, Form("BCM2_Q invalid (%.6g).", m.BCM2_Q)); badNorm = true;
    }
    if (IsBadSentinel(m.comp_livetime) || m.comp_livetime <= 0) {
      warn(run, Form("comp_livetime invalid (%.6g).", m.comp_livetime)); badNorm = true;
    }
    if (IsBadSentinel(m.h_esing_Eff) || m.h_esing_Eff <= 0) {
      warn(run, Form("h_esing_Eff invalid (%.6g).", m.h_esing_Eff)); badNorm = true;
    }
    if (IsBadSentinel(m.boil_corr) || m.boil_corr <= 0) {
      warn(run, Form("boil_corr invalid (%.6g).", m.boil_corr)); badNorm = true;
    }
    if (IsBadSentinel(m.ps_factor) || m.ps_factor <= 0) {
      warn(run, Form("ps_factor invalid (%.6g).", m.ps_factor)); badNorm = true;
    }

    // Open ROOT file
    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      warn(run, "Cannot open ROOT file: " + rootfile);
      status = "NOFILE";
      // still write CSV row (with NaNs)
      csv << m.category << "," << run << "," << m.BCM2_I << ","
          << "nan,nan,"
          << "nan,nan,nan,"
          << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
          << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
          << "nan,nan,nan,"
          << "\"\","  // random list
          << "\"" << rootfile << "\","
          << status << "\n";
      continue;
    }

    TTree* T = (TTree*)f->Get("T");
    if (!T) {
      warn(run, "Tree 'T' not found in file: " + rootfile);
      status = "NOTREE";
      csv << m.category << "," << run << "," << m.BCM2_I << ","
          << "nan,nan,"
          << "nan,nan,nan,"
          << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
          << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
          << "nan,nan,nan,"
          << "\"\","  // random list
          << "\"" << rootfile << "\","
          << status << "\n";
      continue;
    }

    // Compute coincidence random subtraction (net yield in coin peak minus mean random)
    CoincidenceResult R = ComputeCoincidenceRandomSubtraction(T, baseCuts, cfg);

    const double Nnet    = R.RandomSubtractedYield;
    const double NnetErr = R.RandomSubtractedYieldErr;

    // Build normalization factor
    double norm_factor = NAN;
    double y = NAN, yerr = NAN;

    if (!badNorm) {
      norm_factor = (m.boil_corr * m.ps_factor) / (m.comp_livetime * m.h_esing_Eff * m.BCM2_Q);
      y    = Nnet    * norm_factor;
      yerr = NnetErr * norm_factor;
    } else {
      status = "BADMETA";
    }

    // Random window list as a string
    const std::string randList = FormatWindows(R.RandomWindowListNs);

    // Write CSV row (always)
    csv << m.category << "," << run << "," << m.BCM2_I << ","
        << Nnet << "," << NnetErr << ","
        << norm_factor << "," << y << "," << yerr << ","
        << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
        << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
        << R.PeakCenterNs << "," << R.CoinWindowNs.first << "," << R.CoinWindowNs.second << ","
        << "\"" << randList << "\"" << ","
        << "\"" << rootfile << "\","
        << status << "\n";

    // Add to plot series only if finite
    if (std::isfinite(m.BCM2_I) && std::isfinite(y) && std::isfinite(yerr)) {
      series[m.category].x.push_back(m.BCM2_I);
      series[m.category].y.push_back(y);
      series[m.category].ex.push_back(0.0);
      series[m.category].ey.push_back(yerr);
    }
  }

  csv.close();
  log << "\nWrote CSV: " << outCsv << "\n";

  // Sort series by current
  for (auto& kv : series) {
    SortByX(kv.second.x, kv.second.y, kv.second.ex, kv.second.ey);
    log << "Plot points: category=" << kv.first << " n=" << kv.second.x.size() << "\n";
  }

  // Compute global plot ranges (Option A)
  bool any = false;
  double minX=0,maxX=0,minY=0,maxY=0;
  for (const auto& kv : series) {
    const auto& S = kv.second;
    for (size_t i=0;i<S.x.size();++i) {
      double x = S.x[i];
      double y = S.y[i];
      double ey= (i < S.ey.size() ? S.ey[i] : 0.0);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(ey)) continue;
      double ylo = y - ey;
      double yhi = y + ey;
      if (!any) { minX=maxX=x; minY=ylo; maxY=yhi; any=true; }
      else {
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
        minY = std::min(minY, ylo);
        maxY = std::max(maxY, yhi);
      }
    }
  }

  const std::string outPng      = outPNGs + "/yield_vs_current_overlay.png";
  const std::string outRootFile = outCanv + "/yield_vs_current_overlay.root";
  if (!any) {
    warn(-1, "No finite points to plot. Saving empty canvas.");
    TCanvas c("c_yield_vs_current_overlay", "Yield vs Current (myCTime)", 1200, 850);

    c.SaveAs(outPng.c_str());
    c.SaveAs(outRootFile.c_str());

    log << "Wrote PNG (empty): " << outPng << "\n";
    log << "Wrote canvas ROOT (empty): " << outRootFile << "\n";

    log.close();
    return;
  }



  double dx = maxX - minX; if (dx <= 0) dx = 1.0;
  double dy = maxY - minY; if (dy <= 0) dy = 1.0;
  double xpad = 0.06 * dx;
  double ypad = 0.10 * dy;

  double xmin = minX - xpad;
  double xmax = maxX + xpad;
  double ymin = minY - ypad;
  double ymax = maxY + ypad;

  // Styles: colored markers + distinct shapes + larger size
  struct Style { int mstyle; int color; double msize; };
  std::unordered_map<std::string, Style> styles = {
    {"signal",         {20,  2, 1.55}}, // red
    {"positron",       {21,  4, 1.55}}, // blue
    {"dummy",          {22,  8, 1.55}}, // green
    {"positron_dummy", {23,  6, 1.55}}  // magenta
  };

  TCanvas c("c_yield_vs_current_overlay", "Yield vs Current (myCTime)", 1200, 850);
  gStyle->SetOptStat(0);

  c.SetTopMargin(0.22);
  c.SetRightMargin(0.06);
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);

  TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("BCM2_I (mean current)");
  frame->GetYaxis()->SetTitle("Normalized yield (myCTime)");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);

  // Legend in top margin band (won't cover points)
  TLegend leg(0.10, 0.82, 0.55, 0.98);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.035);

  std::vector<std::string> drawOrder = {"signal","positron","dummy","positron_dummy"};
  std::vector<TGraphErrors*> keepAlive;
  keepAlive.reserve(8);

  for (const auto& cat : drawOrder) {
    auto it = series.find(cat);
    if (it == series.end()) continue;
    const auto& S = it->second;
    if (S.x.empty()) continue;

    Style st = styles.count(cat) ? styles[cat] : Style{20,1,1.55};

    TGraphErrors* g = new TGraphErrors((int)S.x.size(),
                                       S.x.data(), S.y.data(),
                                       S.ex.data(), S.ey.data());
    keepAlive.push_back(g);

    g->SetMarkerStyle(st.mstyle);
    g->SetMarkerSize(st.msize);
    g->SetMarkerColor(st.color);
    g->SetLineColor(st.color);

    g->Draw("P SAME");
    leg.AddEntry(g, cat.c_str(), "p");
  }

  // Setting label box in top margin (beside legend)
  std::string lbl1, lbl2;
  BuildSettingLabel(rel, lbl1, lbl2);

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(13);     // left aligned
  t.SetTextSize(0.028);   // tune as needed

  // two lines in the top margin; adjust the y values to control gap
  t.DrawLatex(0.39, 0.94, lbl1.c_str());
  t.DrawLatex(0.39, 0.89, lbl2.c_str());

  // Only draw legend entries for non-empty categories (already handled by AddEntry)
  leg.Draw();

  c.SaveAs(outPng.c_str());
  c.SaveAs(outRootFile.c_str());

  log << "Global ranges: X=[" << xmin << "," << xmax << "]  Y=[" << ymin << "," << ymax << "]\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote Root File: " << outRootFile << "\n";

  for (auto* g : keepAlive) delete g;
  log.close();
}
