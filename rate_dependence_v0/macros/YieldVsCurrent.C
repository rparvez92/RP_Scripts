#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h"

#include "../include/CoincidenceRandomSubtraction.h"

// -------- helpers --------
static std::vector<int> ReadRunList(const std::string& path) {
  std::ifstream f(path);
  std::vector<int> runs;
  if (!f.is_open()) return runs;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    try { runs.push_back(std::stoi(line)); } catch (...) {}
  }
  std::sort(runs.begin(), runs.end());
  runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
  return runs;
}

static std::vector<std::string> SplitCSVLine(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i=0;i<s.size();++i) {
    char c=s[i];
    if (c=='"') { inq=!inq; continue; }
    if (c==',' && !inq) { out.push_back(cur); cur.clear(); }
    else cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

static double SafeToD(const std::string& s) {
  try { return std::stod(s); } catch (...) { return NAN; }
}

static std::string ReadAll(const std::string& path) {
  std::ifstream f(path);
  if (!f.is_open()) return "";
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

static std::string JsonGetString(const std::string& json, const std::string& key, const std::string& def="") {
  std::string pat = "\"" + key + "\"";
  auto p = json.find(pat);
  if (p==std::string::npos) return def;
  p = json.find(":", p);
  if (p==std::string::npos) return def;
  p = json.find("\"", p);
  if (p==std::string::npos) return def;
  auto q = json.find("\"", p+1);
  if (q==std::string::npos) return def;
  return json.substr(p+1, q-p-1);
}

static std::string DataRootPath(int run) {
  char buf[512];
  std::snprintf(buf, sizeof(buf), "./Skimmed_ROOTfiles/skimmed_coin_replay_production_%d_-1.root", run);
  return std::string(buf);
}

static void SortByX(std::vector<double>& x, std::vector<double>& y, std::vector<double>& ex, std::vector<double>& ey) {
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

struct RunMeta {
  std::string category;
  double I_uA = NAN;
  double Q_mC = NAN;
  double lt   = NAN;
  double trk  = NAN;
  double boil = NAN;
  double ps   = 1.0;
  std::string ps_choice;
  double ny   = NAN;
  double ny_e = NAN;
};

static std::unordered_map<int, RunMeta> ReadRunMetadataCSV(const std::string& path, std::ostream& log) {
  std::unordered_map<int, RunMeta> m;
  std::ifstream f(path);
  if (!f.is_open()) {
    log << "ERROR: cannot open run_metadata.csv: " << path << "\n";
    return m;
  }
  std::string header;
  if (!std::getline(f, header)) return m;
  auto cols = SplitCSVLine(header);
  std::unordered_map<std::string,int> idx;
  for (int i=0;i<(int)cols.size();++i) idx[cols[i]] = i;

  auto col = [&](const std::vector<std::string>& row, const std::string& name)->std::string{
    auto it = idx.find(name);
    if (it==idx.end()) return "";
    int j = it->second;
    if (j<0 || j>=(int)row.size()) return "";
    return row[j];
  };

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto row = SplitCSVLine(line);
    int run=-1;
    try { run = std::stoi(col(row,"run")); } catch (...) { continue; }

    RunMeta r;
    r.category = col(row,"category");
    r.I_uA = SafeToD(col(row,"BCM2_I"));
    r.Q_mC = SafeToD(col(row,"BCM2_Q"));
    r.lt   = SafeToD(col(row,"comp_livetime"));
    r.trk  = SafeToD(col(row,"h_esing_Eff"));
    r.boil = SafeToD(col(row,"boil_corr"));
    r.ps   = SafeToD(col(row,"ps_factor"));
    r.ps_choice = col(row,"ps_choice");
    if (!std::isfinite(r.ps) || r.ps<=0) r.ps = 1.0;
    r.ny   = SafeToD(col(row,"normyield"));
    r.ny_e = SafeToD(col(row,"normyield_err"));

    m[run] = r;
  }
  log << "Loaded run metadata entries: " << m.size() << "\n";
  return m;
}

struct CatOut {
  std::vector<double> x_my,y_my,ex_my,ey_my;
  std::vector<double> x_bt,y_bt,ex_bt,ey_bt;
};

static void ProcessCategory(const std::string& cat,
                            const std::string& runlistPath,
                            const std::unordered_map<int, RunMeta>& meta,
                            const TString& baseCuts,
                            const CoincidenceConfig& cfg,
                            const std::string& outTabs,
                            const std::string& outPNGs,
                            const std::string& modeS,
                            std::ostream& log,
                            CatOut& out)
{
  auto runs = ReadRunList(runlistPath);
  log << "\n========== CATEGORY: " << cat << " ==========\n";
  log << "runlist: " << runlistPath << "  n=" << runs.size() << "\n";

  std::string csvPath = outTabs + "/yield_vs_current_" + cat + ".csv";
  std::ofstream csv(csvPath.c_str());
  csv << "category,run,BCM2_I,BCM2_Q_mC,comp_livetime,h_esing_Eff,boil_corr,ps_choice,ps_factor,"
      << "PeakCenterNs,CoinLoNs,CoinHiNs,RandomWindowsNs,"
      << "CoinYield,CoinYieldErr,RandomMeanYield,RandomMeanYieldErr,"
      << "Nnet_myCTime,NnetErr_myCTime,yield_myCTime,yieldErr_myCTime,"
      << "normyield_bigtable,normyieldErr_bigtable,ratio_myCTime_over_bigtable\n";

  for (int run : runs) {
    auto it = meta.find(run);
    if (it == meta.end()) {
      log << "WARNING: run " << run << " missing in run_metadata.csv; skipping.\n";
      continue;
    }
    const RunMeta& M = it->second;

    double yMy=NAN,yMyE=NAN,ratio=NAN;
    CoincidenceResult R;
    R.PeakCenterNs=NAN; R.CoinWindowNs={NAN,NAN};

    if (modeS != "bigtable") {
      const std::string fpath = DataRootPath(run);
      TFile* f = TFile::Open(fpath.c_str(), "READ");
      if (!f || f->IsZombie()) {
        log << "ERROR: cannot open ROOT file run " << run << ": " << fpath << "\n";
        if (f) { f->Close(); delete f; }
        continue;
      }
      TTree* T = (TTree*) f->Get("T");
      if (!T) {
        log << "ERROR: cannot find TTree 'T' in " << fpath << " (run " << run << ")\n";
        f->Close(); delete f;
        continue;
      }

      R = ComputeCoincidenceRandomSubtraction(T, baseCuts, cfg);

      const double Nnet   = R.RandomSubtractedYield;
      const double NnetEr = R.RandomSubtractedYieldErr;

      double norm=NAN;
      if (std::isfinite(M.Q_mC) && M.Q_mC>0 &&
          std::isfinite(M.lt)   && M.lt>0   &&
          std::isfinite(M.trk)  && M.trk>0  &&
          std::isfinite(M.boil) && M.boil>0) {
        norm = (M.boil * M.ps) / (M.lt * M.trk * M.Q_mC);
      }
      if (std::isfinite(norm) && std::isfinite(Nnet)) {
        yMy  = Nnet * norm;
        yMyE = NnetEr * norm;
      }
      if (std::isfinite(yMy) && std::isfinite(M.ny) && M.ny!=0) ratio = yMy / M.ny;

      f->Close(); delete f;
    }

    std::ostringstream rws;
    for (size_t i=0;i<R.RandomWindowListNs.size();++i) {
      if (i) rws << "|";
      rws << R.RandomWindowListNs[i].first << ":" << R.RandomWindowListNs[i].second;
    }

    csv << cat << ","
        << run << ","
        << M.I_uA << ","
        << M.Q_mC << ","
        << M.lt   << ","
        << M.trk  << ","
        << M.boil << ","
        << "\"" << M.ps_choice << "\"" << ","
        << M.ps   << ","
        << R.PeakCenterNs << ","
        << R.CoinWindowNs.first << ","
        << R.CoinWindowNs.second << ","
        << "\"" << rws.str() << "\"" << ","
        << R.CoinYield << ","
        << R.CoinYieldErr << ","
        << R.RandomMeanYield << ","
        << R.RandomMeanYieldErr << ","
        << R.RandomSubtractedYield << ","
        << R.RandomSubtractedYieldErr << ","
        << yMy << ","
        << yMyE << ","
        << M.ny << ","
        << M.ny_e << ","
        << ratio << "\n";

    if (std::isfinite(M.I_uA)) {
      if ((modeS=="myCTime" || modeS=="both") && std::isfinite(yMy)) {
        out.x_my.push_back(M.I_uA); out.y_my.push_back(yMy);
        out.ex_my.push_back(0.0); out.ey_my.push_back(std::isfinite(yMyE)?yMyE:0.0);
      }
      if ((modeS=="bigtable" || modeS=="both") && std::isfinite(M.ny)) {
        out.x_bt.push_back(M.I_uA); out.y_bt.push_back(M.ny);
        out.ex_bt.push_back(0.0); out.ey_bt.push_back(std::isfinite(M.ny_e)?M.ny_e:0.0);
      }
    }

    log << "run " << run
        << " I=" << M.I_uA
        << " Q=" << M.Q_mC
        << " ps=" << M.ps << " (" << M.ps_choice << ")"
        << " lt=" << M.lt << " trk=" << M.trk << " boil=" << M.boil;
    if (modeS != "bigtable") {
      log << " peak=" << R.PeakCenterNs
          << " Nnet=" << R.RandomSubtractedYield << " +/- " << R.RandomSubtractedYieldErr
          << " yMy=" << yMy << " +/- " << yMyE;
    }
    log << " yBT=" << M.ny << " +/- " << M.ny_e << "\n";
  }

  csv.close();
  log << "Wrote CSV: " << csvPath << "\n";

  if (!out.x_my.empty()) SortByX(out.x_my,out.y_my,out.ex_my,out.ey_my);
  if (!out.x_bt.empty()) SortByX(out.x_bt,out.y_bt,out.ex_bt,out.ey_bt);

  // Per-category plot
  TCanvas c(("c_"+cat).c_str(),"Yield vs Current",1100,800);
  gStyle->SetOptStat(0);
  TLegend leg(0.62,0.75,0.88,0.88);
  leg.SetBorderSize(0);
  bool first=true;

  if ((modeS=="myCTime" || modeS=="both") && !out.x_my.empty()) {
    TGraphErrors g((int)out.x_my.size(), out.x_my.data(), out.y_my.data(), out.ex_my.data(), out.ey_my.data());
    g.SetMarkerStyle(20); g.SetMarkerSize(1.1);
    g.SetTitle("");
    g.GetXaxis()->SetTitle("BCM2_I (mean current)");
    g.GetYaxis()->SetTitle("Normalized yield");
    g.Draw(first?"AP":"P SAME"); first=false;
    leg.AddEntry(&g, "myCTime", "p");
  }
  if ((modeS=="bigtable" || modeS=="both") && !out.x_bt.empty()) {
    TGraphErrors g2((int)out.x_bt.size(), out.x_bt.data(), out.y_bt.data(), out.ex_bt.data(), out.ey_bt.data());
    g2.SetMarkerStyle(24); g2.SetMarkerSize(1.1);
    g2.SetTitle("");
    if (first) {
      g2.GetXaxis()->SetTitle("BCM2_I (mean current)");
      g2.GetYaxis()->SetTitle("Normalized yield");
      g2.Draw("AP"); first=false;
    } else {
      g2.Draw("P SAME");
    }
    leg.AddEntry(&g2, "bigtable", "p");
  }
  leg.Draw();
  std::string outPng = outPNGs + "/yield_vs_current_" + cat + ".png";
  c.SaveAs(outPng.c_str());
  log << "Wrote PNG: " << outPng << "\n";
}

static void OverlayPlot(const std::string& outPNGs,
                        const std::string& whichSeries,
                        const std::unordered_map<std::string, CatOut>& cats,
                        std::ostream& log)
{
  TCanvas c(("c_overlay_"+whichSeries).c_str(),"Overlay",1200,850);
  gStyle->SetOptStat(0);
  TLegend leg(0.62,0.70,0.88,0.88);
  leg.SetBorderSize(0);

  std::vector<std::pair<std::string,int>> styles = {
    {"signal", 20},
    {"positron", 21},
    {"dummy", 22},
    {"positron_dummy", 23}
  };

  bool first=true;
  for (auto& cs : styles) {
    const auto& cat = cs.first;
    auto it = cats.find(cat);
    if (it==cats.end()) continue;
    const CatOut& co = it->second;

    const std::vector<double>* x=nullptr; const std::vector<double>* y=nullptr;
    const std::vector<double>* ex=nullptr; const std::vector<double>* ey=nullptr;

    if (whichSeries=="myCTime") { x=&co.x_my; y=&co.y_my; ex=&co.ex_my; ey=&co.ey_my; }
    else { x=&co.x_bt; y=&co.y_bt; ex=&co.ex_bt; ey=&co.ey_bt; }

    if (x->empty()) continue;

    TGraphErrors g((int)x->size(), x->data(), y->data(), ex->data(), ey->data());
    g.SetMarkerStyle(cs.second);
    g.SetMarkerSize(1.1);
    g.SetTitle("");
    if (first) {
      g.GetXaxis()->SetTitle("BCM2_I (mean current)");
      g.GetYaxis()->SetTitle("Normalized yield");
      g.Draw("AP"); first=false;
    } else {
      g.Draw("P SAME");
    }
    leg.AddEntry(&g, cat.c_str(), "p");
  }
  leg.Draw();
  std::string outPng = outPNGs + "/yield_vs_current_overlay_" + whichSeries + ".png";
  c.SaveAs(outPng.c_str());
  log << "Wrote overlay PNG: " << outPng << "\n";
}

void YieldVsCurrent(const char* manifestPath,
                    const char* settingsDir,
                    const char* resultsDir,
                    const char* mode = "both")
{
  std::string modeS(mode ? mode : "both");
  if (modeS!="myCTime" && modeS!="bigtable" && modeS!="both") modeS="both";

  std::string manifestP(manifestPath);
  std::string settingDir = gSystem->DirName(manifestP.c_str());

  std::string mtxt = ReadAll(manifestP);
  std::string setting_id = JsonGetString(mtxt, "setting_id", gSystem->BaseName(settingDir.c_str()));
  if (setting_id.empty()) setting_id = gSystem->BaseName(settingDir.c_str());

  std::string outBase = std::string(resultsDir) + "/" + setting_id;
  std::string outPNGs = outBase + "/PNGs";
  std::string outTabs = outBase + "/tables";
  gSystem->mkdir(outPNGs.c_str(), kTRUE);
  gSystem->mkdir(outTabs.c_str(), kTRUE);

  std::ofstream log(std::string(outBase + "/log.txt").c_str());
  log << "YieldVsCurrent (v5-like categories)\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "mode: " << modeS << "\n\n";

  std::string metaPath = settingDir + "/run_metadata.csv";
  auto meta = ReadRunMetadataCSV(metaPath, log);
  if (meta.empty()) {
    log << "ERROR: no run metadata loaded.\n";
    return;
  }

  const TString baseCuts =
    "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
    "(H_cal_etottracknorm>0.7) && "
    "(H_cer_npeSum>2.0) && "
    "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
    "(P_cal_etottracknorm<0.8)";
  log << "BaseCuts: " << baseCuts.Data() << "\n";

  CoincidenceConfig cfg;

  std::vector<std::string> cats = {"signal","positron","dummy","positron_dummy"};
  std::unordered_map<std::string, CatOut> catOuts;

  for (const auto& cat : cats) {
    std::string runlistPath = settingDir + "/runs_" + cat + ".txt";
    CatOut out;
    ProcessCategory(cat, runlistPath, meta, baseCuts, cfg, outTabs, outPNGs, modeS, log, out);
    catOuts[cat] = std::move(out);
  }

  if (modeS=="myCTime" || modeS=="both") OverlayPlot(outPNGs, "myCTime", catOuts, log);
  if (modeS=="bigtable" || modeS=="both") OverlayPlot(outPNGs, "bigtable", catOuts, log);

  log.close();
}
