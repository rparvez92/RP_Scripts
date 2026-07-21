// get_good_coin_ev_check.C
//
// Quick diagnostic comparing the online/quick good-coin accounting from
// get_good_coin_ev.C against the RP-side CoincidenceRandomSubtraction.h
// analysis-window calculation.
//
// Run:
//   root -l -b -q 'macros/get_good_coin_ev_check.C'

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TTree.h"

#include "../include/CoincidenceRandomSubtraction.h"

namespace {

struct RunSetting {
  int run = 0;
  double beamE = 0.0;
  std::string target;
  std::string polarity;
  double x = 0.0;
  double Q2 = 0.0;
  double z = 0.0;
};

struct GGRow {
  double coin = NAN;
  double randoms = NAN;
  double goodcoin = NAN;
  double goodcoinErr = NAN;
  double normyield = NAN;
  double normyieldErr = NAN;
  double ctmean = NAN;
  double ctsigma = NAN;
};

struct NormInputs {
  double charge_mC = NAN;
  double compLT = NAN;
  double hmsTrackEff = NAN;
  double shmsTrackEff = NAN;
  double trigEff = 1.0;
};

struct CountResult {
  double coin = 0.0;
  double random = 0.0;
  double good = 0.0;
  double goodErr = 0.0;
  double peakCenter = NAN;
  int nRandomWindows = 0;
};

std::string Trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

std::vector<std::string> SplitCSVLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (char c : line) {
    if (c == '"') {
      inq = !inq;
      continue;
    }
    if (c == ',' && !inq) {
      out.push_back(Trim(cur));
      cur.clear();
      continue;
    }
    cur.push_back(c);
  }
  out.push_back(Trim(cur));
  return out;
}

std::string NormalizeSpaces(const std::string& s) {
  std::string out;
  bool prevSpace = false;
  for (char c : s) {
    const bool isSpace = (c == ' ' || c == '\t' || c == '\r' || c == '\n');
    if (isSpace) {
      if (!prevSpace) out.push_back(' ');
      prevSpace = true;
    } else {
      out.push_back(c);
      prevSpace = false;
    }
  }
  return Trim(out);
}

double ExtractValueFromReportFile(const std::string& filename,
                                  const std::string& key,
                                  const char delimiter,
                                  int skipCount = 0)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "ERROR: unable to open report file: " << filename << "\n";
    return NAN;
  }

  const std::string normalizedKey = NormalizeSpaces(key);
  std::string line;
  int matchCount = 0;
  while (std::getline(file, line)) {
    const std::string normalizedLine = NormalizeSpaces(line);
    const size_t pos = normalizedLine.find(normalizedKey);
    if (pos == std::string::npos) continue;

    bool validMatch = false;
    if (pos == 0 || std::isspace((unsigned char)normalizedLine[pos - 1])) {
      size_t searchPos = pos + normalizedKey.length();
      while (searchPos < normalizedLine.length() &&
             std::isspace((unsigned char)normalizedLine[searchPos])) {
        ++searchPos;
      }
      validMatch = (searchPos < normalizedLine.length() &&
                    normalizedLine[searchPos] == delimiter);
    }
    if (!validMatch) continue;

    if (matchCount < skipCount) {
      ++matchCount;
      continue;
    }

    std::istringstream iss(normalizedLine.substr(normalizedLine.find(delimiter) + 1));
    double value = NAN;
    iss >> value;
    return value;
  }

  std::cerr << "ERROR: key not found in report: " << key
            << " after skipping " << skipCount << " matches in " << filename << "\n";
  return NAN;
}

bool LoadGGCSV(const std::string& path, GGRow& row) {
  std::ifstream f(path);
  if (!f.is_open()) {
    std::cerr << "ERROR: cannot open GG CSV: " << path << "\n";
    return false;
  }

  std::string header, line;
  if (!std::getline(f, header) || !std::getline(f, line)) {
    std::cerr << "ERROR: malformed GG CSV: " << path << "\n";
    return false;
  }

  const auto h = SplitCSVLine(header);
  const auto v = SplitCSVLine(line);
  std::unordered_map<std::string, int> idx;
  for (int i = 0; i < (int)h.size(); ++i) idx[h[i]] = i;

  auto get = [&](const char* name) -> double {
    auto it = idx.find(name);
    if (it == idx.end() || it->second < 0 || it->second >= (int)v.size()) return NAN;
    return std::atof(v[it->second].c_str());
  };

  row.coin = get("coin");
  row.randoms = get("randoms");
  row.goodcoin = get("ransubcoin");
  row.goodcoinErr = get("ransubcoin_err");
  row.normyield = get("normyield");
  row.normyieldErr = get("normyield_err");
  row.ctmean = get("ctmean");
  row.ctsigma = get("ctsigma");
  return true;
}

NormInputs LoadNormInputs(const std::string& reportPath) {
  NormInputs n;
  n.charge_mC = ExtractValueFromReportFile(reportPath, "SHMS BCM2 Beam Cut Charge", ':', 0);
  n.compLT = ExtractValueFromReportFile(reportPath, "HMS Computer Dead Time", ':', 0) / 100.0;
  n.hmsTrackEff = ExtractValueFromReportFile(reportPath, "E SING FID TRACK EFFIC", ':', 1);
  n.shmsTrackEff = ExtractValueFromReportFile(reportPath, "HADRON SING FID TRACK EFFIC", ':', 0);
  n.trigEff = 1.0;
  return n;
}

double NormFactor(const NormInputs& n) {
  if (!(n.charge_mC > 0.0) || !(n.compLT > 0.0) ||
      !(n.hmsTrackEff > 0.0) || !(n.shmsTrackEff > 0.0) || !(n.trigEff > 0.0)) {
    return NAN;
  }
  return 1.0 / (n.charge_mC * n.compLT * n.hmsTrackEff * n.shmsTrackEff * n.trigEff);
}

CountResult CountAnalysisWindow(TTree* T, const TString& cuts) {
  CountResult r;
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime.ePiCoinTime_ROC2";

  CoincidenceResult cr = ComputeCoincidenceRandomSubtraction(T, cuts, cfg);
  r.coin = cr.CoinYield;
  r.random = cr.RandomMeanYield;
  r.good = cr.RandomSubtractedYield;
  r.goodErr = cr.RandomSubtractedYieldErr;
  r.peakCenter = cr.PeakCenterNs;
  r.nRandomWindows = (int)cr.RandomWindowListNs.size();
  return r;
}

double SafeRatio(double num, double den) {
  if (!std::isfinite(num) || !std::isfinite(den) || den == 0.0) return NAN;
  return num / den;
}

void WriteCSVHeader(std::ofstream& out) {
  out
    << "run,beam_energy,target,polarity,x,Q2,z,"
    << "charge_mC,livetime,hms_eff,shms_eff,"
    << "RP_peak,RP_n_random_window,"
    << "RP_coin,RP_random_mean,RP_goodcoin,RP_goodcoin_err,"
    << "RP_normyield,RP_normyield_err,"
    << "GG_goodcoin,GG_goodcoin_err,GG_normyield,GG_normyield_err,"
    << "ratio_RP_by_GG\n";
}

} // namespace

int get_good_coin_ev_check(
  const std::string& rootDir = "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/ROOTfiles",
  const std::string& reportDir = "/net/cdaq/cdaql3data/cdaq/hallc-online-rsidis2025/REPORT_OUTPUT/COIN/PRODUCTION",
  const std::string& outCsv = "")
{
  const std::vector<RunSetting> runs = {
    {27122, 6.449, "LH2", "-", 0.22, 2.2, 0.50},
    {27163, 6.449, "LH2", "+", 0.22, 2.2, 0.50},
    {27300, 6.449, "LD2", "+", 0.44, 4.4, 0.52},
    {27358, 6.449, "LH2", "-", 0.44, 4.4, 0.67},
  };

  const TString anaCuts =
    "(P.gtr.p<=2.9||P.hgcer.npeSum>1)&&"
    "P.aero.npeSum>2&&"
    "H.cer.npeSum>2&&"
    "H.cal.etottracknorm>0.7&&"
    "P.cal.etottracknorm<0.8&&"
    "abs(P.gtr.dp-5.)<15.&&"
    "abs(H.gtr.dp)<8.";

  std::string finalOutCsv = outCsv;
  if (finalOutCsv.empty()) {
    const TString macroDir = gSystem->DirName(__FILE__);
    const TString projectDir = gSystem->DirName(macroDir);
    finalOutCsv = std::string(projectDir.Data()) + "/results/tables/get_good_coin_ev_check.csv";
  }

  const std::string outDir = finalOutCsv.substr(0, finalOutCsv.find_last_of('/'));
  if (!outDir.empty() && outDir != finalOutCsv) gSystem->mkdir(outDir.c_str(), true);

  std::ofstream out(finalOutCsv);
  if (!out.is_open()) {
    std::cerr << "ERROR: cannot write output CSV: " << finalOutCsv << "\n";
    return 1;
  }
  WriteCSVHeader(out);

  for (const auto& rs : runs) {
    const std::string rootPath =
      rootDir + Form("/coin_replay_production_%d_-1.root", rs.run);
    const std::string reportPath =
      reportDir + Form("/replay_coin_production_%d_-1.report", rs.run);
    const std::string ggCsvPath =
      reportDir + Form("/output_get_good_coin_ev_%d_-1.csv", rs.run);

    GGRow gg;
    LoadGGCSV(ggCsvPath, gg);

    const NormInputs norm = LoadNormInputs(reportPath);
    const double nf = NormFactor(norm);

    TFile* f = TFile::Open(rootPath.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: cannot open ROOT file: " << rootPath << "\n";
      if (f) f->Close();
      continue;
    }

    TTree* T = (TTree*)f->Get("T");
    if (!T) {
      std::cerr << "ERROR: missing tree T in: " << rootPath << "\n";
      f->Close();
      continue;
    }

    const CountResult rpAna = CountAnalysisWindow(T, anaCuts);

    const double rpAnaYield = rpAna.good * nf;
    const double rpAnaYieldErr = rpAna.goodErr * nf;

    out
      << rs.run << ","
      << rs.beamE << ","
      << rs.target << ","
      << rs.polarity << ","
      << rs.x << ","
      << rs.Q2 << ","
      << rs.z << ","
      << norm.charge_mC << ","
      << norm.compLT << ","
      << norm.hmsTrackEff << ","
      << norm.shmsTrackEff << ","
      << rpAna.peakCenter << ","
      << rpAna.nRandomWindows << ","
      << rpAna.coin << ","
      << rpAna.random << ","
      << rpAna.good << ","
      << rpAna.goodErr << ","
      << rpAnaYield << ","
      << rpAnaYieldErr << ","
      << gg.goodcoin << ","
      << gg.goodcoinErr << ","
      << gg.normyield << ","
      << gg.normyieldErr << ","
      << SafeRatio(rpAnaYield, gg.normyield) << "\n";

    std::cout << "Processed run " << rs.run
              << " | RP_Ana/GG=" << SafeRatio(rpAnaYield, gg.normyield)
              << " | coin=" << rpAna.coin
              << " random_mean=" << rpAna.random
              << " good=" << rpAna.good << "\n";

    f->Close();
  }

  std::cout << "Wrote " << finalOutCsv << "\n";
  return 0;
}
