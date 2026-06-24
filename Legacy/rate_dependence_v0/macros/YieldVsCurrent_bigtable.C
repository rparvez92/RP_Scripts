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
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TString.h"
#include "TPad.h"
#include "TH1.h"

// ---------- helpers ----------
static double SafeToD(const std::string& s) {
  try { return std::stod(s); } catch (...) { return NAN; }
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

struct Row {
  std::string category;
  int run = -1;
  double I  = NAN;
  double y  = NAN;
  double ey = NAN;

  // extras for csv/provenance (not used in plot)
  std::string target_raw;
  std::string run_type_raw;
  double Q_mC = NAN;
  double lt = NAN;
  double trk = NAN;
  double boil = NAN;
  std::string ps_choice;
  double ps_factor = NAN;
};

static std::vector<Row> ReadRunMetadata(const std::string& path, std::ostream& log) {
  std::vector<Row> rows;
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

    Row row;
    row.category = col(r, "category");
    try { row.run = std::stoi(col(r, "run")); } catch (...) { continue; }

    row.I  = SafeToD(col(r, "BCM2_I"));
    row.y  = SafeToD(col(r, "normyield"));
    row.ey = SafeToD(col(r, "normyield_err"));

    row.target_raw   = col(r, "target_raw");
    row.run_type_raw = col(r, "run_type_raw");

    row.Q_mC      = SafeToD(col(r, "BCM2_Q"));
    row.lt        = SafeToD(col(r, "comp_livetime"));
    row.trk       = SafeToD(col(r, "h_esing_Eff"));
    row.boil      = SafeToD(col(r, "boil_corr"));
    row.ps_choice = col(r, "ps_choice");
    row.ps_factor = SafeToD(col(r, "ps_factor"));

    rows.push_back(row);
  }

  log << "Loaded rows from run_metadata.csv: " << rows.size() << "\n";
  return rows;
}

static void WriteOverlayCSV(const std::string& outCsv, const std::vector<Row>& rows, std::ostream& log) {
  std::ofstream csv(outCsv.c_str());
  csv << "category,run,BCM2_I,normyield,normyield_err,"
      << "target_raw,run_type_raw,BCM2_Q_mC,comp_livetime,h_esing_Eff,boil_corr,ps_choice,ps_factor\n";
  for (const auto& r : rows) {
    csv << r.category << ","
        << r.run << ","
        << r.I << ","
        << r.y << ","
        << r.ey << ","
        << "\"" << r.target_raw << "\"" << ","
        << "\"" << r.run_type_raw << "\"" << ","
        << r.Q_mC << ","
        << r.lt << ","
        << r.trk << ","
        << r.boil << ","
        << "\"" << r.ps_choice << "\"" << ","
        << r.ps_factor
        << "\n";
  }
  csv.close();
  log << "Wrote CSV: " << outCsv << "\n";
}

// ---------- main ----------
void YieldVsCurrent_bigtable(const char* manifestPath,
                             const char* resultsBigtableDir)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsBigtableDir ? resultsBigtableDir : "";

  // Determine setting directory and id
  const std::string settingDir = DirName(manifestP);
  const std::string setting_id = BaseName(settingDir);

  // Mirror everything after /settings/ under results_bigtable
  const std::string rel = RelUnderSettings(settingDir);
  const std::string outBase = outRoot + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTabs = outBase + "/tables";
  EnsureDir(outPNGs);
  EnsureDir(outTabs);

  std::ofstream log((outBase + "/log.txt").c_str());
  log << "YieldVsCurrent_bigtable\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "results_base: " << outBase << "\n\n";

  // Read run_metadata.csv from the setting directory
  const std::string metaPath = settingDir + "/run_metadata.csv";
  auto rows = ReadRunMetadata(metaPath, log);
  if (rows.empty()) {
    log << "ERROR: no rows found; nothing to plot.\n";
    return;
  }

  // Always write overlay CSV (no filtering)
  const std::string outCsv = outTabs + "/yield_vs_current_overlay.csv";
  WriteOverlayCSV(outCsv, rows, log);

  // Build per-category series for plotting
  struct Series { std::vector<double> x,y,ex,ey; };
  std::map<std::string, Series> series;

  for (const auto& r : rows) {
    if (!std::isfinite(r.I) || !std::isfinite(r.y)) continue;     // for plotting only
    double e = (std::isfinite(r.ey) ? r.ey : 0.0);
    series[r.category].x.push_back(r.I);
    series[r.category].y.push_back(r.y);
    series[r.category].ex.push_back(0.0);
    series[r.category].ey.push_back(e);
  }

  // Sort each category by current and log counts
  for (auto& kv : series) {
    SortByX(kv.second.x, kv.second.y, kv.second.ex, kv.second.ey);
    log << "Plot points: category=" << kv.first << " n=" << kv.second.x.size() << "\n";
  }

  // Compute global axis ranges over ALL categories (Option A)
  bool any = false;
  double minX = 0, maxX = 0, minY = 0, maxY = 0;

  for (const auto& kv : series) {
    const auto& S = kv.second;
    for (size_t i=0;i<S.x.size();++i) {
      const double x  = S.x[i];
      const double y  = S.y[i];
      const double ey = (i < S.ey.size() ? S.ey[i] : 0.0);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(ey)) continue;

      const double ylo = y - ey;
      const double yhi = y + ey;

      if (!any) {
        minX = maxX = x;
        minY = ylo;
        maxY = yhi;
        any = true;
      } else {
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
        minY = std::min(minY, ylo);
        maxY = std::max(maxY, yhi);
      }
    }
  }

  if (!any) {
    log << "WARNING: no finite (BCM2_I, normyield) points to plot.\n";
    // Still save an empty canvas
    TCanvas c("c_yield_vs_current_overlay", "Yield vs Current (bigtable)", 1200, 850);
    c.SaveAs((outPNGs + "/yield_vs_current_overlay.png").c_str());
    log.close();
    return;
  }

  // Add padding
  double dx = maxX - minX;
  double dy = maxY - minY;
  if (dx <= 0) dx = 1.0;
  if (dy <= 0) dy = 1.0;

  const double xpad = 0.06 * dx;
  const double ypad = 0.10 * dy;

  const double xmin = minX - xpad;
  const double xmax = maxX + xpad;
  const double ymin = minY - ypad;
  const double ymax = maxY + ypad;

  // Styles (colored markers + distinct shapes + large markers)
  struct Style { int mstyle; int color; double msize; };
  std::unordered_map<std::string, Style> styles = {
    {"signal",         {20,  2, 1.55}}, // red
    {"positron",       {21,  4, 1.55}}, // blue
    {"dummy",          {22,  8, 1.55}}, // green
    {"positron_dummy", {23,  6, 1.55}}  // magenta
  };

  // Draw
  TCanvas c("c_yield_vs_current_overlay", "Yield vs Current (bigtable)", 1200, 850);
  gStyle->SetOptStat(0);

  // Top margin space for legend
  c.SetTopMargin(0.22);
  c.SetRightMargin(0.06);
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);

  // Create an empty frame with GLOBAL ranges
  //gPad->DrawFrame(xmin, ymin, xmax, ymax);
  //gPad->Update();

  // Axis titles / sizes
  //gPad->GetFrame()->GetXaxis()->SetTitle("BCM2_I (mean current)");
  //gPad->GetFrame()->GetYaxis()->SetTitle("Normalized yield (normyield)");
  //gPad->GetFrame()->GetXaxis()->SetTitleSize(0.045);
  //gPad->GetFrame()->GetYaxis()->SetTitleSize(0.045);
  //gPad->GetFrame()->GetXaxis()->SetLabelSize(0.04);
  //gPad->GetFrame()->GetYaxis()->SetLabelSize(0.04);


  TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  // no need for gPad->Update() here

  frame->GetXaxis()->SetTitle("BCM2_I (mean current)");
  frame->GetYaxis()->SetTitle("Normalized yield (normyield)");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);

  // Legend in the top margin region (so it doesn't cover points)
  TLegend leg(0.15, 0.82, 0.60, 0.98);
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

    auto stIt = styles.find(cat);
    Style st = (stIt != styles.end()) ? stIt->second : Style{20, 1, 1.55};

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

  leg.Draw();

  const std::string outPng = outPNGs + "/yield_vs_current_overlay.png";
  c.SaveAs(outPng.c_str());
  log << "Global ranges: X=[" << xmin << "," << xmax << "]  Y=[" << ymin << "," << ymax << "]\n";
  log << "Wrote PNG: " << outPng << "\n";

  for (auto* g : keepAlive) delete g;
  log.close();
}
