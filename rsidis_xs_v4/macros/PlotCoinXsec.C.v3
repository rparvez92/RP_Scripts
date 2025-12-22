// PlotCoinXsec.C
//
// Example run command:
// root -l -b -q 'macros/PlotCoinXsec.C+("/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v3/settings/pass4/pi+sidis/LH2/z0p5/x0p25/q23p3/hmsPneg1p531_shmsP3p632_hmsTh29p045_shmsTh7p865_thpq2/manifest.json", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v3/results/test_one_setting")'

// Inputs expected in the setting directory:
//   manifest.json
//   runs_signal.txt
//   runs_dummy.txt
//   runs_positron.txt
//   runs_positron_dummy.txt
//   run_metadata.csv
//   sigma_model.txt
//
// Output:
//   Writes log.txt into outputDir
//   Writes plots into outputDir/PNGs

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cctype>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include "../include/Mapping.h"
#include "../include/PlotComparisonAndRatio.h"
#include "../include/CoincidenceRandomSubtraction.h"

namespace {

  std::vector<std::unique_ptr<TH1>> g_keep_hists;

  struct RunNormInfo {
    double charge_mC = 0.0;
    double hms_eff = 1.0;
    double comp_livetime = 1.0;

    double ps5 = 0.0;
    double ps6 = 0.0;
    std::string ps_choice = "";
    double ps_factor = 1.0;

    double fan_mean = 1.0;
    double boil_corr = 1.0;
    double electr_deadtime = 0.0;
  };

  static const double kTolPs = 1e-6;

  bool IsOne(double x) {
    return std::isfinite(x) && std::abs(x - 1.0) < kTolPs;
  }

  std::pair<std::string,double> ChoosePsFactor(double ps5, double ps6) {
    bool is5 = IsOne(ps5);
    bool is6 = IsOne(ps6);
    if (is5 && !is6) return {"ps5", ps5};
    if (is6 && !is5) return {"ps6", ps6};
    if (is5 && is6)  return {"ps5", ps5};
    return {"ps5", ps5};
  }

  double ToDoubleSafe(const std::string& s, double fallback) {
    try { return std::stod(s); } catch (...) { return fallback; }
  }

  std::vector<std::string> SplitCSVLine(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(item);
    return out;
  }

  std::vector<int> ReadRunListTxt(const std::string& path) {
    std::vector<int> runs;
    std::ifstream f(path);
    if (!f.is_open()) return runs;
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      if (line[0] == '#') continue;
      runs.push_back(std::stoi(line));
    }
    std::sort(runs.begin(), runs.end());
    runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
    return runs;
  }

  std::vector<std::string> ReadCorrectionsEnabledFromOutputMeta(const std::string& outDir) {
    std::vector<std::string> enabled;
    std::ifstream f(outDir + "/output_meta.json");
    if (!f.is_open()) return enabled;

    std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    auto pos = all.find("\"corrections_enabled\"");
    if (pos == std::string::npos) return enabled;
    pos = all.find('[', pos);
    if (pos == std::string::npos) return enabled;
    auto end = all.find(']', pos);
    if (end == std::string::npos) return enabled;

    std::string arr = all.substr(pos + 1, end - pos - 1);
    std::stringstream ss(arr);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      tok.erase(std::remove(tok.begin(), tok.end(), '"'), tok.end());
      tok.erase(std::remove_if(tok.begin(), tok.end(),
                               [](unsigned char c){ return std::isspace(c) != 0; }),
                tok.end());
      if (!tok.empty()) enabled.push_back(tok);
    }
    return enabled;
  }

  bool HasEnabled(const std::vector<std::string>& enabled, const std::string& key) {
    for (const auto& e : enabled) if (e == key) return true;
    return false;
  }

  double ComputePerRunWeight(const RunNormInfo& info, const std::vector<std::string>& enabled) {
    double w = 1.0;

    if (HasEnabled(enabled, "ps_factor")) w *= info.ps_factor;

    if (HasEnabled(enabled, "h_esing_Eff")) {
      if (info.hms_eff > 0.0) w /= info.hms_eff;
    }

    if (HasEnabled(enabled, "comp_livetime")) {
      if (info.comp_livetime > 0.0) w /= info.comp_livetime;
    }

    if (HasEnabled(enabled, "fan_mean"))  w *= info.fan_mean;
    if (HasEnabled(enabled, "boil_corr")) w *= info.boil_corr;

    return w;
  }

  bool LoadRunMetadataCSV(const std::string& csvPath,
                          std::unordered_map<int, RunNormInfo>& out)
  {
    std::ifstream f(csvPath);
    if (!f.is_open()) return false;

    std::string header;
    if (!std::getline(f, header)) return false;
    auto cols = SplitCSVLine(header);

    auto idx = [&](const std::string& name) -> int {
      for (size_t i = 0; i < cols.size(); ++i) if (cols[i] == name) return (int)i;
      return -1;
    };

    int i_run       = idx("run");
    int i_charge    = idx("BCM4A_Q");
    int i_eff       = idx("h_esing_Eff");
    int i_live      = idx("comp_livetime");
    int i_ps5       = idx("ps5");
    int i_ps6       = idx("ps6");
    int i_ps_choice = idx("ps_choice");
    int i_ps_factor = idx("ps_factor");
    int i_fan       = idx("fan_mean");
    int i_boil      = idx("boil_corr");
    int i_edead     = idx("electr_deadtime");

    if (i_run < 0 || i_charge < 0 || i_eff < 0) return false;

    std::string line;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      auto fields = SplitCSVLine(line);

      if (i_run >= (int)fields.size()) continue;
      int run = std::stoi(fields[i_run]);

      RunNormInfo info;
      info.charge_mC       = (i_charge >= 0 && i_charge < (int)fields.size()) ? ToDoubleSafe(fields[i_charge], 0.0) : 0.0;
      info.hms_eff         = (i_eff    >= 0 && i_eff    < (int)fields.size()) ? ToDoubleSafe(fields[i_eff], 1.0)    : 1.0;
      info.comp_livetime   = (i_live   >= 0 && i_live   < (int)fields.size()) ? ToDoubleSafe(fields[i_live], 1.0)   : 1.0;
      info.ps5             = (i_ps5    >= 0 && i_ps5    < (int)fields.size()) ? ToDoubleSafe(fields[i_ps5], 0.0)    : 0.0;
      info.ps6             = (i_ps6    >= 0 && i_ps6    < (int)fields.size()) ? ToDoubleSafe(fields[i_ps6], 0.0)    : 0.0;
      info.ps_choice       = (i_ps_choice>=0 && i_ps_choice<(int)fields.size()) ? fields[i_ps_choice] : "";
      info.ps_factor       = (i_ps_factor>=0 && i_ps_factor<(int)fields.size()) ? ToDoubleSafe(fields[i_ps_factor], 1.0) : 1.0;
      info.fan_mean        = (i_fan  >= 0 && i_fan  < (int)fields.size()) ? ToDoubleSafe(fields[i_fan], 1.0) : 1.0;
      info.boil_corr       = (i_boil >= 0 && i_boil < (int)fields.size()) ? ToDoubleSafe(fields[i_boil], 1.0) : 1.0;
      info.electr_deadtime = (i_edead>=0 && i_edead<(int)fields.size()) ? ToDoubleSafe(fields[i_edead], 0.0) : 0.0;

      if (!std::isfinite(info.ps_factor) || info.ps_factor <= 0.0) {
        auto ch = ChoosePsFactor(info.ps5, info.ps6);
        info.ps_choice = ch.first;
        info.ps_factor = ch.second;
      }

      out[run] = info;
    }

    return !out.empty();
  }

  double ReadSigmaModelUb(const std::string& settingDir, std::ofstream& log) {
    std::string path = settingDir + "/sigma_model.txt";
    std::ifstream f(path);
    if (!f.is_open()) {
      log << "ERROR cannot open sigma_model.txt at " << path << "\n";
      return 0.0;
    }
    double sigma = 0.0;
    f >> sigma;
    if (!std::isfinite(sigma) || sigma <= 0.0) {
      log << "ERROR bad sigma_model value: " << sigma << "\n";
      return 0.0;
    }
    log << "sigma_model_ub=" << sigma << "\n";
    return sigma;
  }

  std::string DataRootPath(int run) {
    char buf[512];
    std::snprintf(buf, sizeof(buf), "./Pass0_ROOTfiles/coin_replay_production_%d_-1.root", run);
    return std::string(buf);
  }

  std::unique_ptr<TH1D> ProjectRun(
      int run,
      const std::string& dataVarExpr,
      const std::string& histName,
      int nbins, double xmin, double xmax,
      const TCut& baseCutsBool,
      double& totalCharge_mC,
      const std::unordered_map<int, RunNormInfo>& infoMap,
      const std::vector<std::string>& enabled,
      std::unordered_map<int, CoincidenceResult>& coinCache,
      std::ofstream& log)
  {
    auto it = infoMap.find(run);
    if (it == infoMap.end()) {
      log << "ERROR missing metadata for run " << run << "\n";
      return nullptr;
    }
    const RunNormInfo& info = it->second;

    totalCharge_mC += info.charge_mC;

    double w = ComputePerRunWeight(info, enabled);

    log << "RUN " << run
        << " charge_mC=" << info.charge_mC
        << " ps_choice=" << info.ps_choice
        << " ps_factor=" << info.ps_factor
        << " hms_eff=" << info.hms_eff
        << " comp_livetime=" << info.comp_livetime
        << " per_run_weight=" << w
        << "\n";

    std::string path = DataRootPath(run);
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
      log << "ERROR cannot open " << path << "\n";
      if (f) f->Close();
      return nullptr;
    }

    TTree* T = dynamic_cast<TTree*>(f->Get("T"));
    if (!T) {
      log << "ERROR missing tree T in " << path << "\n";
      f->Close();
      return nullptr;
    }

    auto h = std::make_unique<TH1D>(histName.c_str(), histName.c_str(), nbins, xmin, xmax);
    h->Sumw2(true);

    CoincidenceConfig cfg;
    TString baseCutsStr = TString(baseCutsBool.GetTitle());
    TString weightExpr = Form("%g", w);

    auto itC = coinCache.find(run);
    if (itC == coinCache.end()) {
      CoincidenceResult Rwin = ComputeCoincidenceRandomSubtraction(T, baseCutsStr, cfg);
      coinCache[run] = Rwin;
      itC = coinCache.find(run);
    }

    FillRandomSubtractedHistogramWithWindows(T,
                                            baseCutsStr,
                                            weightExpr,
                                            dataVarExpr.c_str(),
                                            h.get(),
                                            cfg,
                                            itC->second);

    h->SetDirectory(nullptr);
    f->Close();
    return h;
  }

  std::unique_ptr<TH1D> BuildYield(
      const std::vector<int>& runs,
      const std::string& dataVarExpr,
      const std::string& baseName,
      int nbins, double xmin, double xmax,
      const TCut& baseCutsBool,
      double& totalCharge_mC,
      const std::unordered_map<int, RunNormInfo>& infoMap,
      const std::vector<std::string>& enabled,
      std::unordered_map<int, CoincidenceResult>& coinCache,
      std::ofstream& log)
  {
    totalCharge_mC = 0.0;

    auto sum = std::make_unique<TH1D>((baseName + "_sum").c_str(), "", nbins, xmin, xmax);
    sum->Sumw2(true);
    sum->Reset();
    sum->SetDirectory(nullptr);

    for (int run : runs) {
      auto h = ProjectRun(run, dataVarExpr, baseName + Form("_run%d", run),
                          nbins, xmin, xmax, baseCutsBool,
                          totalCharge_mC, infoMap, enabled, coinCache, log);
      if (!h) continue;
      sum->Add(h.get());
      g_keep_hists.push_back(std::move(h));
    }

    if (totalCharge_mC > 0.0) sum->Scale(1.0 / totalCharge_mC);
    return sum;
  }

  std::unique_ptr<TH1D> BuildSimHist(
      const std::string& simVar,
      TTree* simTree,
      const std::string& name,
      int nbins, double xmin, double xmax,
      const TCut& simCuts)
  {
    auto h = std::make_unique<TH1D>(name.c_str(), name.c_str(), nbins, xmin, xmax);
    h->Sumw2(true);

    TString drawCmd = TString::Format("%s>>%s(%d,%f,%f)", simVar.c_str(), name.c_str(), nbins, xmin, xmax);

    Long64_t ngen = simTree->GetEntries();
    simTree->Draw(drawCmd, simCuts, "goff");
    if (ngen > 0) h->Scale(1.0 / (double)ngen);

    h->SetDirectory(nullptr);
    return h;
  }

  struct VarConfig {
    std::string simVar;
    int nbins = 0;
    double xmin = 0.0;
    double xmax = 0.0;
  };

  struct YieldSet {
    std::unordered_map<std::string, std::unique_ptr<TH1D>> hSig;
    std::unordered_map<std::string, std::unique_ptr<TH1D>> hDum;
    std::unordered_map<std::string, std::unique_ptr<TH1D>> hPosSig;
    std::unordered_map<std::string, std::unique_ptr<TH1D>> hPosDum;
    double qSig = 0.0;
    double qDum = 0.0;
    double qPosSig = 0.0;
    double qPosDum = 0.0;
  };

  std::unordered_map<std::string, std::unique_ptr<TH1D>> BookSums(
      const std::vector<VarConfig>& vars,
      const std::string& prefix)
  {
    std::unordered_map<std::string, std::unique_ptr<TH1D>> out;
    for (const auto& v : vars) {
      std::string name = prefix + "_" + v.simVar + "_sum";
      auto h = std::make_unique<TH1D>(name.c_str(), "", v.nbins, v.xmin, v.xmax);
      h->Sumw2(true);
      h->Reset();
      h->SetDirectory(nullptr);
      out[v.simVar] = std::move(h);
    }
    return out;
  }

  void AccumulateCategory(
      const std::vector<int>& runs,
      const std::vector<VarConfig>& vars,
      std::unordered_map<std::string, std::unique_ptr<TH1D>>& outSums,
      double& totalCharge_mC,
      const TCut& baseCutsBool,
      const std::unordered_map<int, RunNormInfo>& infoMap,
      const std::vector<std::string>& enabled,
      std::unordered_map<int, CoincidenceResult>& coinCache,
      std::ofstream& log)
  {
    totalCharge_mC = 0.0;

    for (int run : runs) {
      auto it = infoMap.find(run);
      if (it == infoMap.end()) {
        log << "ERROR missing metadata for run " << run << "\n";
        continue;
      }
      const RunNormInfo& info = it->second;

      totalCharge_mC += info.charge_mC;

      double w = ComputePerRunWeight(info, enabled);

      log << "RUN " << run
          << " charge_mC=" << info.charge_mC
          << " ps_choice=" << info.ps_choice
          << " ps_factor=" << info.ps_factor
          << " hms_eff=" << info.hms_eff
          << " comp_livetime=" << info.comp_livetime
          << " per_run_weight=" << w
          << "\n";

      std::string path = DataRootPath(run);
      TFile* f = TFile::Open(path.c_str(), "READ");
      if (!f || f->IsZombie()) {
        log << "ERROR cannot open " << path << "\n";
        if (f) f->Close();
        continue;
      }

      TTree* T = dynamic_cast<TTree*>(f->Get("T"));
      if (!T) {
        log << "ERROR missing tree T in " << path << "\n";
        f->Close();
        continue;
      }

      CoincidenceConfig cfg;
      TString baseCutsStr = TString(baseCutsBool.GetTitle());
      TString weightExpr = Form("%g", w);

      auto itC = coinCache.find(run);
      if (itC == coinCache.end()) {
        CoincidenceResult Rwin = ComputeCoincidenceRandomSubtraction(T, baseCutsStr, cfg);
        coinCache[run] = Rwin;
        itC = coinCache.find(run);
      }

      for (const auto& v : vars) {
        std::string dataExpr = SimToDataMap(v.simVar);

        std::string hname = "hTmp_" + v.simVar + Form("_run%d", run);
        auto htmp = std::make_unique<TH1D>(hname.c_str(), hname.c_str(), v.nbins, v.xmin, v.xmax);
        htmp->Sumw2(true);
        htmp->SetDirectory(nullptr);

        FillRandomSubtractedHistogramWithWindows(T,
                                                baseCutsStr,
                                                weightExpr,
                                                dataExpr.c_str(),
                                                htmp.get(),
                                                cfg,
                                                itC->second);

        auto itSum = outSums.find(v.simVar);
        if (itSum != outSums.end() && itSum->second) {
          itSum->second->Add(htmp.get());
        }

        g_keep_hists.push_back(std::move(htmp));
      }

      f->Close();
    }

    if (totalCharge_mC > 0.0) {
      for (auto& kv : outSums) {
        if (kv.second) kv.second->Scale(1.0 / totalCharge_mC);
      }
    }
  }

  YieldSet BuildAllYields(
      const std::vector<int>& runsSignal,
      const std::vector<int>& runsDummy,
      const std::vector<int>& runsPosSignal,
      const std::vector<int>& runsPosDummy,
      const std::vector<VarConfig>& vars,
      const TCut& baseCutsBool,
      const std::unordered_map<int, RunNormInfo>& infoMap,
      const std::vector<std::string>& enabled,
      std::unordered_map<int, CoincidenceResult>& coinCache,
      std::ofstream& log)
  {
    YieldSet ys;
    ys.hSig    = BookSums(vars, "hSig");
    ys.hDum    = BookSums(vars, "hDum");
    ys.hPosSig = BookSums(vars, "hPosSig");
    ys.hPosDum = BookSums(vars, "hPosDum");

    log << "Accumulating signal runs\n";
    AccumulateCategory(runsSignal, vars, ys.hSig, ys.qSig, baseCutsBool, infoMap, enabled, coinCache, log);

    log << "Accumulating dummy runs\n";
    AccumulateCategory(runsDummy, vars, ys.hDum, ys.qDum, baseCutsBool, infoMap, enabled, coinCache, log);

    log << "Accumulating positron signal runs\n";
    AccumulateCategory(runsPosSignal, vars, ys.hPosSig, ys.qPosSig, baseCutsBool, infoMap, enabled, coinCache, log);

    log << "Accumulating positron dummy runs\n";
    AccumulateCategory(runsPosDummy, vars, ys.hPosDum, ys.qPosDum, baseCutsBool, infoMap, enabled, coinCache, log);

    return ys;
  }

  // Cuts
  TCut BuildCoinCut()
  {
    TCut HmsPid =
      "(H.gtr.dp>-8) && (H.gtr.dp<8) && "
      "(H.cal.etottracknorm>0.7) && "
      "(H.cer.npeSum>2.0)";

    TCut ShmsBase =
      "(P.gtr.dp>-10) && (P.gtr.dp<22) && "
      "(P.cal.etottracknorm<0.8)";

    TCut ShmsAero = "(P.gtr.p<2.7) && (P.aero.npeSum>2)";
    TCut ShmsHGC  = "(P.gtr.p>=2.7) && (P.hgcer.npeSum>1) && (P.aero.npeSum>2)";

    TCut ShmsPid = ShmsBase && (ShmsAero || ShmsHGC);

    return (HmsPid && ShmsPid);
  }

  TCut BuildSimDeltaCuts(double simNormalizationFactor)
  {
    TCut sim_delta = "((hsdelta>-8.0) && (hsdelta<8.0) && (ssdelta>-10.0) && (ssdelta<22.0))";
    TCut sim_norm  = Form("Weight * %e", simNormalizationFactor);
    return sim_delta * sim_norm; // cut * weight
  }

} // namespace


void PlotCoinXsec(const char* manifestPath,
                  const char* outputDir,
                  bool makeComparisonRatio = true,
                  bool makeXsec = true)
{
  gROOT->SetBatch(kTRUE);

  std::string manifest = manifestPath ? manifestPath : "";
  std::string outDir = outputDir ? outputDir : "";
  std::string settingDir = gSystem->DirName(manifest.c_str());

  gSystem->mkdir(outDir.c_str(), true);
  std::ofstream log(outDir + "/log.txt");

  log << "PlotCoinXsec\n";
  log << "manifest=" << manifest << "\n";
  log << "setting_dir=" << settingDir << "\n";
  log << "output_dir=" << outDir << "\n";

  std::vector<std::string> enabled = ReadCorrectionsEnabledFromOutputMeta(outDir);
  log << "corrections_enabled=";
  for (size_t i = 0; i < enabled.size(); ++i) {
    log << enabled[i];
    if (i + 1 < enabled.size()) log << ",";
  }
  log << "\n";

  std::vector<int> runsSignal    = ReadRunListTxt(settingDir + "/runs_signal.txt");
  std::vector<int> runsDummy     = ReadRunListTxt(settingDir + "/runs_dummy.txt");
  std::vector<int> runsPosSignal = ReadRunListTxt(settingDir + "/runs_positron.txt");
  std::vector<int> runsPosDummy  = ReadRunListTxt(settingDir + "/runs_positron_dummy.txt");

  log << "runs_signal_n=" << runsSignal.size() << "\n";
  log << "runs_dummy_n=" << runsDummy.size() << "\n";
  log << "runs_positron_signal_n=" << runsPosSignal.size() << "\n";
  log << "runs_positron_dummy_n=" << runsPosDummy.size() << "\n";

  if (!runsPosSignal.empty() && runsPosDummy.empty()) {
    log << "WARNING runs_positron_dummy.txt is empty. Dummy positron subtraction will be skipped.\n";
  }

  std::unordered_map<int, RunNormInfo> infoMap;
  if (!LoadRunMetadataCSV(settingDir + "/run_metadata.csv", infoMap)) {
    log << "ERROR cannot read run_metadata.csv\n";
    return;
  }

  double sigma_model_ub = ReadSigmaModelUb(settingDir, log);
  if (!(sigma_model_ub > 0.0)) {
    log << "ERROR sigma_model_ub is not valid. Stop.\n";
    return;
  }

  std::string simFilePath = "./simc_worksim/coin_7p87deg_3p632gev_hyd_rsidis.root";
  double simNormalizationFactor = 0.842205e11;

  TFile* fSim = TFile::Open(simFilePath.c_str(), "READ");
  if (!fSim || fSim->IsZombie()) {
    log << "ERROR cannot open SIM file " << simFilePath << "\n";
    if (fSim) fSim->Close();
    return;
  }

  TTree* tSim = dynamic_cast<TTree*>(fSim->Get("h10"));
  if (!tSim) {
    log << "ERROR missing SIM tree h10\n";
    fSim->Close();
    return;
  }

  //TCut sim_delta_cuts = "((hsdelta>-8.0) && (hsdelta<8.0))";
  //TCut sim_norm_cuts  = Form("Weight * %e", simNormalizationFactor);

  /*TCut data_base_cuts =
      "((H.gtr.dp>-8.0) && (H.gtr.dp<8.0) && "
      "(H.cal.etottracknorm>0.7) && (H.cer.npeSum>2.0))";*/
  TCut data_base_cuts = BuildCoinCut();
  TCut simCuts = BuildSimDeltaCuts(simNormalizationFactor);

  double wall_thickness_ratio = 3.82;
  std::vector<VarConfig> vars;
  //vars.push_back({"Q2", 30, 1.0, 6.0});
  //vars.push_back({"W",  30, 1.0, 6.0});
  //vars.push_back({"xbj",30, 0.0, 1.0});
  vars.push_back({"phipq",16, 0.0, 7.0});
  //vars.push_back({"hsdelta",30, -15.0, 15.0});


  for (const auto& v : vars) {
    log << "plot_var=" << v.simVar << " nbins=" << v.nbins << " xmin=" << v.xmin << " xmax=" << v.xmax << "\n";
  }

  std::unordered_map<int, CoincidenceResult> coinCache;

  YieldSet ys = BuildAllYields(runsSignal, runsDummy, runsPosSignal, runsPosDummy,
                              vars,
                              data_base_cuts,
                              infoMap, enabled,
                              coinCache,
                              log);

  std::unordered_map<std::string, std::unique_ptr<TH1D>> hSimMap;
  for (const auto& v : vars) {
    std::string name = "hSim_" + v.simVar;
    auto hsim = BuildSimHist(v.simVar, tSim, name, v.nbins, v.xmin, v.xmax, simCuts);
    if (hsim) {
      hsim->SetDirectory(nullptr);
    }
    hSimMap[v.simVar] = std::move(hsim);
  }

  gSystem->mkdir((outDir + "/PNGs").c_str(), true);

  for (const auto& v : vars) {
    TH1D* hSim = hSimMap[v.simVar].get();
    TH1D* hSig = ys.hSig[v.simVar].get();
    TH1D* hDum = ys.hDum[v.simVar].get();
    TH1D* hPosSig = ys.hPosSig[v.simVar].get();
    TH1D* hPosDum = ys.hPosDum[v.simVar].get();

    auto hSigMinusPos = std::make_unique<TH1D>(*hSig);
    hSigMinusPos->SetName(("hSigMinusPos_" + v.simVar).c_str());
    if (!runsPosSignal.empty()) hSigMinusPos->Add(hPosSig, -1.0);

    auto hDumMinusPos = std::make_unique<TH1D>(*hDum);
    hDumMinusPos->SetName(("hDumMinusPos_" + v.simVar).c_str());
    if (!runsPosDummy.empty()) hDumMinusPos->Add(hPosDum, -1.0);

    auto hFinal = std::make_unique<TH1D>(*hSigMinusPos);
    hFinal->SetName(("hFinal_" + v.simVar).c_str());
    hFinal->Add(hDumMinusPos.get(), -(1.0 / wall_thickness_ratio));

    if (makeComparisonRatio) {
      PlotComparisonAndRatio(hFinal.get(), hSim, v.simVar, outDir.c_str());
    }
    if (makeXsec) {
      PlotXsec(hFinal.get(), hSim, sigma_model_ub, v.simVar, outDir.c_str());
    }

    g_keep_hists.push_back(std::move(hSigMinusPos));
    g_keep_hists.push_back(std::move(hDumMinusPos));
    g_keep_hists.push_back(std::move(hFinal));
  }

  fSim->Close();
  log << "DONE\n";

}

void PlotCoinXsec()
{
  std::cout << "Loaded PlotCoinXsec. Use PlotCoinXsec(\"manifest.json\",\"outputDir\")\n";
}
