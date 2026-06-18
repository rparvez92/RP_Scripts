// CoincidenceBlockingByRun.C
//
// Purpose:
//   For each run listed in a setting's run_metadata.csv, measure a simple
//   coincidence-blocking ratio from the pass0 replay files using
//   CTime.CoinTime_RAW_ROC2 over [-400, 400].
//
// Definitions:
//   raw window:
//     (CTime.CoinTime_RAW_ROC2 > -400) && (CTime.CoinTime_RAW_ROC2 < 400)
//
//   with-cut:
//     raw window &&
//     (P.hod.goodstarttime == 1) &&
//     (H.hod.goodstarttime == 1)
//
//   good time:
//     with-cut &&
//     (CTime.CoinTime_RAW_ROC2 > 5) &&
//     (CTime.CoinTime_RAW_ROC2 < 70)
//
//   coincidence_blocking_ratio = good time / with-cut
//
// Inputs:
//   CoincidenceBlockingByRun("<abs path>/settings/.../manifest.json", "<abs path>/results")
//
// Reads:
//   <setting_dir>/run_metadata.csv
//   ./Pass0p1_ROOTfiles/coin_replay_production_<RUN>_-1.root
//
// Writes:
//   <results>/<same rel path as manifest>/tables/coincidence_blocking_by_run.csv
//   <results>/<same rel path as manifest>/logs/CoincidenceBlockingByRun.log
//   <results>/<same rel path as manifest>/PNGs/coincidence_blocking_ratio_vs_run.png
//   <results>/<same rel path as manifest>/canvases/coincidence_blocking_ratio_vs_run.root
//
// Single-setting run example:
//   root -l -b -q 'macros/CoincidenceBlockingByRun.C("settings/pass4/pi+sidis/LH2/z0p36/x0p25/Q23p3/hmsPneg1p531_shmsP2p615_hmsTh29p045_shmsTh7p865_thpq2/manifest.json","results")'
//
// tcsh batch run over all settings:
//   set resabs = `pwd`/results
//   foreach m (`find settings -name manifest.json | sort`)
//     set mabs = `realpath $m`
//     set rel  = `echo "$m" | sed -e 's#^settings/##' -e 's#/manifest.json$##'`
//     set logdir = "$resabs/$rel/logs"
//     mkdir -p "$logdir"
//     root -l -b -q 'macros/CoincidenceBlockingByRun.C("'"$mabs"'","'"$resabs"'")' >&! "$logdir/CoincidenceBlockingByRun.batch.log"
//   end

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

static std::string SafeTrim(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace((unsigned char)s[b])) ++b;
  while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
  return s.substr(b, e - b);
}

static std::vector<std::string> SplitCSVLine_Blocking(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < s.size(); ++i) {
    char c = s[i];
    if (c == '"') {
      inq = !inq;
      continue;
    }
    if (c == ',' && !inq) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

static void EnsureDir_Blocking(const std::string& p) { gSystem->mkdir(p.c_str(), kTRUE); }
static std::string DirName_Blocking(const std::string& p) { return std::string(gSystem->DirName(p.c_str())); }
static std::string BaseName_Blocking(const std::string& p) { return std::string(gSystem->BaseName(p.c_str())); }

static std::string RelUnderSettings_Blocking(const std::string& settingDirAbs) {
  const std::string needle = "/settings/";
  auto pos = settingDirAbs.find(needle);
  if (pos == std::string::npos) return BaseName_Blocking(settingDirAbs);
  std::string rel = settingDirAbs.substr(pos + needle.size());
  while (!rel.empty() && rel.front() == '/') rel.erase(0, 1);
  return rel.empty() ? BaseName_Blocking(settingDirAbs) : rel;
}

static std::vector<std::string> SplitPath_Blocking(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == '/') {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

static void BuildSettingLabel_Blocking(const std::string& relUnderSettings,
                                       std::string& line1,
                                       std::string& line2) {
  auto tok = SplitPath_Blocking(relUnderSettings);
  if (tok.size() >= 7) {
    line1 = tok[0] + " / " + tok[1] + " / " + tok[2] + " / " +
            tok[3] + " / " + tok[4] + " / " + tok[5];
    line2 = tok[6];
  } else {
    line1 = relUnderSettings;
    line2 = "";
  }
}

struct MetaRowBlocking {
  std::string category;
  int run = -1;
};

static std::vector<MetaRowBlocking> ReadRunMetadata_Blocking(const std::string& path, std::ostream& log) {
  std::vector<MetaRowBlocking> rows;
  std::ifstream f(path);
  if (!f.is_open()) {
    log << "ERROR: cannot open " << path << "\n";
    return rows;
  }

  std::string header;
  if (!std::getline(f, header)) return rows;

  auto cols = SplitCSVLine_Blocking(header);
  std::unordered_map<std::string, int> idx;
  for (int i = 0; i < (int)cols.size(); ++i) idx[cols[i]] = i;

  auto col = [&](const std::vector<std::string>& r, const std::string& name) -> std::string {
    auto it = idx.find(name);
    if (it == idx.end()) return "";
    const int j = it->second;
    if (j < 0 || j >= (int)r.size()) return "";
    return r[j];
  };

  std::string line;
  while (std::getline(f, line)) {
    if (SafeTrim(line).empty()) continue;
    auto r = SplitCSVLine_Blocking(line);

    MetaRowBlocking m;
    m.category = col(r, "category");
    try {
      m.run = std::stoi(col(r, "run"));
    } catch (...) {
      continue;
    }
    rows.push_back(m);
  }

  log << "Loaded rows from run_metadata.csv: " << rows.size() << "\n";
  return rows;
}

void CoincidenceBlockingByRun(const char* manifestPath,
                              const char* resultsDir)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsDir   ? resultsDir   : "";

  const std::string settingDir = DirName_Blocking(manifestP);
  const std::string setting_id = BaseName_Blocking(settingDir);
  const std::string rel        = RelUnderSettings_Blocking(settingDir);
  const std::string outBase    = outRoot + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTabs = outBase + "/tables";
  const std::string outCanv = outBase + "/canvases";
  const std::string outLogs = outBase + "/logs";
  EnsureDir_Blocking(outPNGs);
  EnsureDir_Blocking(outTabs);
  EnsureDir_Blocking(outCanv);
  EnsureDir_Blocking(outLogs);

  const std::string outCsv  = outTabs + "/coincidence_blocking_by_run.csv";
  const std::string logPath = outLogs + "/CoincidenceBlockingByRun.log";
  std::ofstream log(logPath.c_str());

  log << "CoincidenceBlockingByRun\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "results_base: " << outBase << "\n";
  log << "log: " << logPath << "\n\n";

  const std::string metaPath = settingDir + "/run_metadata.csv";
  auto meta = ReadRunMetadata_Blocking(metaPath, log);
  if (meta.empty()) {
    log << "ERROR: no metadata rows. Nothing to do.\n";
    return;
  }

  std::ofstream csv(outCsv.c_str());
  csv << "category,run,ctime_raw_withstart,ctime_good_withstart,coinc_blocking_ratio,rootfile,status\n";

  const std::string pass0Dir = "./Pass0p1_ROOTfiles";
  const TString ctExpr = "CTime.CoinTime_RAW_ROC2";
  const TString withStartCut =
      "((CTime.CoinTime_RAW_ROC2>-400) && (CTime.CoinTime_RAW_ROC2<400) && "
      "(P.hod.goodstarttime==1) && "
      "(H.hod.goodstarttime==1))";
  const TString goodTimeCut =
      "((CTime.CoinTime_RAW_ROC2>5) && (CTime.CoinTime_RAW_ROC2<70) && "
      "(P.hod.goodstarttime==1) && "
      "(H.hod.goodstarttime==1))";

  std::vector<double> runPlot;
  std::vector<double> ratioPlot;

  for (const auto& m : meta) {
    const int run = m.run;
    const std::string rootfile =
        Form("%s/coin_replay_production_%d_-1.root", pass0Dir.c_str(), run);

    double withStartYield = NAN;
    double goodTimeYield = NAN;
    double ratio = NAN;
    std::string status = "OK";

    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      status = "NOFILE";
      log << "WARNING [run " << run << "]: cannot open ROOT file: " << rootfile << "\n";
    } else {
      TTree* T = (TTree*)f->Get("T");
      if (!T) {
        status = "NOTREE";
        log << "WARNING [run " << run << "]: tree T not found in " << rootfile << "\n";
      } else {
        const TString hWithStartName = Form("hCoinRawWithStart_%d", run);
        const TString hGoodName = Form("hCoinRawGood_%d", run);
        std::unique_ptr<TH1D> hWithStart(new TH1D(hWithStartName, "", 800, -400.0, 400.0));
        std::unique_ptr<TH1D> hGood(new TH1D(hGoodName, "", 800, -400.0, 400.0));

        T->Project(hWithStartName, ctExpr, withStartCut);
        T->Project(hGoodName, ctExpr, goodTimeCut);

        withStartYield = hWithStart->Integral();
        goodTimeYield = hGood->Integral();
        if (std::isfinite(withStartYield) && withStartYield > 0.0 && std::isfinite(goodTimeYield)) {
          ratio = goodTimeYield / withStartYield;
        } else {
          status = "BAD_RATIO";
          log << "WARNING [run " << run << "]: invalid with-start/good-time counts. with-start="
              << withStartYield << " good-time=" << goodTimeYield << "\n";
        }
      }
    }

    csv << m.category << ","
        << run << ","
        << withStartYield << ","
        << goodTimeYield << ","
        << ratio << ","
        << "\"" << rootfile << "\","
        << status << "\n";

    if (status == "OK" && std::isfinite(ratio)) {
      runPlot.push_back((double)run);
      ratioPlot.push_back(ratio);
    }
  }

  csv.close();
  log << "\nWrote CSV: " << outCsv << "\n";
  std::cout << "Wrote CSV: " << outCsv << "\n";

  const std::string outPng = outPNGs + "/coincidence_blocking_ratio_vs_run.png";
  const std::string outRootCanvas = outCanv + "/coincidence_blocking_ratio_vs_run.root";

  gStyle->SetOptStat(0);
  TCanvas c("c_coincidence_blocking_ratio_vs_run", "Coincidence Blocking Ratio vs Run", 1200, 850);
  c.SetGrid();

  if (runPlot.empty()) {
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);
    lat.DrawLatex(0.15, 0.55, "No valid coincidence_blocking_ratio points");
  } else {
    TGraph g((int)runPlot.size(), runPlot.data(), ratioPlot.data());
    g.SetTitle("");
    g.SetMarkerStyle(20);
    g.SetMarkerSize(1.0);
    g.SetMarkerColor(kBlue + 1);
    g.SetLineColor(kBlue + 1);
    g.GetXaxis()->SetTitle("Run number");
    g.GetYaxis()->SetTitle("coincidence_blocking_ratio");
    g.GetYaxis()->SetTitleOffset(1.2);
    g.Draw("AP");

    std::string line1, line2;
    BuildSettingLabel_Blocking(rel, line1, line2);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.13, 0.93, line1.c_str());
    if (!line2.empty()) lat.DrawLatex(0.13, 0.885, line2.c_str());
  }

  c.SaveAs(outPng.c_str());
  c.SaveAs(outRootCanvas.c_str());
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT canvas: " << outRootCanvas << "\n";
  log.close();
}
