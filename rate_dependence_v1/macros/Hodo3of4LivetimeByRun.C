// Hodo3of4LivetimeByRun.C
//
// Purpose:
//   For each run listed in a setting's run_metadata.csv, estimate the SHMS
//   hodoscope 3-of-4 trigger electronic livetime using the Dave Mack
//   combinatoric/incoherent-rate approximation.
//
// Physics model:
//   The SHMS hodo 3-of-4 trigger is formed from four plane logic signals:
//     P1X, P1Y, P2X, P2Y
//
//   For each plane i:
//     D_i = R_i * DPR
//     L_i = 1 - D_i
//
//   where:
//     R_i = plane rate read directly from the replay report, converted to Hz
//     DPR = double-pulse resolution in seconds
//
//   The default DPR is 48 ns. The macro accepts dpr_ns as an argument so the
//   user can rerun the correction if the trigger history or group consensus
//   changes.
//
//   The incoherent 3-of-4 hodo livetime is:
//     LT_3of4 =
//       L1*L2*L3*L4
//       + D1*L2*L3*L4
//       + L1*D2*L3*L4
//       + L1*L2*D3*L4
//       + L1*L2*L3*D4
//
//   This macro uses the report rates directly. It does not iteratively correct
//   the plane scaler rates for this same deadtime.
//
// Inputs:
//   Hodo3of4LivetimeByRun("<abs path>/settings/.../manifest.json",
//                         "<abs path>/results",
//                         48.0)
//
// Reads:
//   <setting_dir>/run_metadata.csv
//   ./Pass0p1_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report
//
//   Report parsing is label-based, not line-number-based. The macro searches
//   for the P1X, P1Y, P2X, and P2Y lines and extracts the value inside:
//     [ <rate> kHz ]
//
// Writes:
//   <results>/<same rel path as manifest>/tables/hodo3of4_livetime_by_run.csv
//   <results>/<same rel path as manifest>/logs/Hodo3of4LivetimeByRun.log
//
// Workflow:
//   1. Run CoincidenceBlockingByRun.C to prepare coincidence-blocking ratios.
//   2. Run this macro to prepare hodo 3-of-4 livetimes.
//   3. Run YieldVsCurrent.C, which reads both CSVs by run number and writes
//      staged yields:
//        yield_no_CB_corr
//        yield_CB_corr
//        yield_no_hodo3of4_corr
//        yield_hodo3of4_corr
//
// Single-setting run example:
//   root -l -b -q 'macros/Hodo3of4LivetimeByRun.C("settings/pass4/pi+sidis/LH2/z0p36/x0p25/Q23p3/hmsPneg1p531_shmsP2p615_hmsTh29p045_shmsTh7p865_thpq2/manifest.json","results",48.0)'
//
// tcsh batch run over all settings:
//   set RESULTS = "$cwd/results"
//   foreach mf (`find settings -name manifest.json | sort`)
//     set rel = `echo "$mf" | sed 's|^settings/||'`
//     set rel_dir = `echo "$rel" | sed 's|/manifest.json$||'`
//     mkdir -p "$RESULTS/$rel_dir/logs"
//     echo "RUN Hodo3of4LivetimeByRun: $mf"
//     root -l -b -q 'macros/Hodo3of4LivetimeByRun.C("'"$cwd/$mf"'","'"$RESULTS"'",48.0)' \
//       >&! "$RESULTS/$rel_dir/logs/Hodo3of4LivetimeByRun.batch.log"
//   end

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "TString.h"
#include "TSystem.h"

namespace {

static std::string TrimHodoLT(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace((unsigned char)s[b])) ++b;
  while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
  return s.substr(b, e - b);
}

static std::vector<std::string> SplitCSVLineHodoLT(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
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

static void EnsureDirHodoLT(const std::string& p) { gSystem->mkdir(p.c_str(), kTRUE); }
static std::string DirNameHodoLT(const std::string& p) { return std::string(gSystem->DirName(p.c_str())); }
static std::string BaseNameHodoLT(const std::string& p) { return std::string(gSystem->BaseName(p.c_str())); }

static std::string RelUnderSettingsHodoLT(const std::string& settingDirAbs) {
  const std::string needle = "/settings/";
  auto pos = settingDirAbs.find(needle);
  if (pos == std::string::npos) return BaseNameHodoLT(settingDirAbs);
  std::string rel = settingDirAbs.substr(pos + needle.size());
  while (!rel.empty() && rel.front() == '/') rel.erase(0, 1);
  return rel.empty() ? BaseNameHodoLT(settingDirAbs) : rel;
}

struct MetaRowHodoLT {
  std::string category;
  int run = -1;
};

static std::vector<MetaRowHodoLT>
ReadRunMetadataHodoLT(const std::string& path, std::ostream& log) {
  std::vector<MetaRowHodoLT> rows;
  std::ifstream f(path);
  if (!f.is_open()) {
    log << "ERROR: cannot open " << path << "\n";
    return rows;
  }

  std::string header;
  if (!std::getline(f, header)) return rows;

  auto cols = SplitCSVLineHodoLT(header);
  std::unordered_map<std::string, int> idx;
  for (int i = 0; i < (int)cols.size(); ++i) idx[TrimHodoLT(cols[i])] = i;

  auto col = [&](const std::vector<std::string>& r, const std::string& name) -> std::string {
    auto it = idx.find(name);
    if (it == idx.end()) return "";
    const int j = it->second;
    if (j < 0 || j >= (int)r.size()) return "";
    return TrimHodoLT(r[j]);
  };

  std::string line;
  while (std::getline(f, line)) {
    if (TrimHodoLT(line).empty()) continue;
    auto r = SplitCSVLineHodoLT(line);

    MetaRowHodoLT m;
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

static bool StartsWithLabelHodoLT(const std::string& line, const std::string& label) {
  const std::string t = TrimHodoLT(line);
  if (t.rfind(label, 0) != 0) return false;
  const size_t p = label.size();
  return p < t.size() && std::isspace((unsigned char)t[p]);
}

static bool ParsePlaneRateFromReportHodoLT(const std::string& reportPath,
                                           const std::string& label,
                                           double& rate_kHz) {
  std::ifstream in(reportPath);
  if (!in) return false;

  std::string line;
  while (std::getline(in, line)) {
    if (!StartsWithLabelHodoLT(line, label)) continue;
    const size_t lb = line.find('[');
    const size_t rb = line.find(']', lb == std::string::npos ? 0 : lb);
    if (lb == std::string::npos || rb == std::string::npos || rb <= lb) continue;

    std::string inside = TrimHodoLT(line.substr(lb + 1, rb - lb - 1));
    if (inside.find("kHz") == std::string::npos) continue;

    std::stringstream ss(inside);
    double val = NAN;
    ss >> val;
    if (ss.fail() || !std::isfinite(val)) continue;
    rate_kHz = val;
    return true;
  }
  return false;
}

struct HodoLTRow {
  std::string category;
  int run = -1;
  double p1x_kHz = NAN;
  double p1y_kHz = NAN;
  double p2x_kHz = NAN;
  double p2y_kHz = NAN;
  double dpr_ns = NAN;
  double d1 = NAN, d2 = NAN, d3 = NAN, d4 = NAN;
  double l1 = NAN, l2 = NAN, l3 = NAN, l4 = NAN;
  double lt = NAN;
  std::string reportfile;
  std::string status;
};

static void ComputeHodoLT(HodoLTRow& r) {
  const double dpr_s = r.dpr_ns * 1e-9;
  r.d1 = r.p1x_kHz * 1000.0 * dpr_s;
  r.d2 = r.p1y_kHz * 1000.0 * dpr_s;
  r.d3 = r.p2x_kHz * 1000.0 * dpr_s;
  r.d4 = r.p2y_kHz * 1000.0 * dpr_s;

  r.l1 = 1.0 - r.d1;
  r.l2 = 1.0 - r.d2;
  r.l3 = 1.0 - r.d3;
  r.l4 = 1.0 - r.d4;

  r.lt = r.l1 * r.l2 * r.l3 * r.l4
       + r.d1 * r.l2 * r.l3 * r.l4
       + r.l1 * r.d2 * r.l3 * r.l4
       + r.l1 * r.l2 * r.d3 * r.l4
       + r.l1 * r.l2 * r.l3 * r.d4;
}

} // namespace

void Hodo3of4LivetimeByRun(const char* manifestPath,
                            const char* resultsDir,
                            double dpr_ns = 48.0)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsDir   ? resultsDir   : "";

  const std::string settingDir = DirNameHodoLT(manifestP);
  const std::string setting_id = BaseNameHodoLT(settingDir);
  const std::string rel        = RelUnderSettingsHodoLT(settingDir);
  const std::string outBase    = outRoot + "/" + rel;

  const std::string outTabs = outBase + "/tables";
  const std::string outLogs = outBase + "/logs";
  EnsureDirHodoLT(outTabs);
  EnsureDirHodoLT(outLogs);

  const std::string outCsv  = outTabs + "/hodo3of4_livetime_by_run.csv";
  const std::string logPath = outLogs + "/Hodo3of4LivetimeByRun.log";
  std::ofstream log(logPath.c_str());

  log << "Hodo3of4LivetimeByRun\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "results_base: " << outBase << "\n";
  log << "DPR_ns: " << dpr_ns << "\n";
  log << "log: " << logPath << "\n\n";

  const std::string metaPath = settingDir + "/run_metadata.csv";
  auto meta = ReadRunMetadataHodoLT(metaPath, log);
  if (meta.empty()) {
    log << "ERROR: no metadata rows. Nothing to do.\n";
    return;
  }

  const std::string reportDir = "./Pass0p1_REPORTfiles/COIN/PRODUCTION";
  log << "Report dir: " << reportDir << "\n";

  std::vector<HodoLTRow> rows;
  rows.reserve(meta.size());

  int nOK = 0, nMissingReport = 0, nMissingPlane = 0, nBadLT = 0;
  for (const auto& m : meta) {
    HodoLTRow r;
    r.category = m.category;
    r.run = m.run;
    r.dpr_ns = dpr_ns;
    r.reportfile = Form("%s/replay_coin_production_%d_-1.report", reportDir.c_str(), r.run);
    r.status = "OK";

    std::ifstream test(r.reportfile.c_str());
    if (!test) {
      r.status = "MISSING_REPORT";
      ++nMissingReport;
      log << "WARNING [run " << r.run << "]: missing report file: " << r.reportfile << "\n";
      rows.push_back(r);
      continue;
    }

    const bool gotP1X = ParsePlaneRateFromReportHodoLT(r.reportfile, "P1X", r.p1x_kHz);
    const bool gotP1Y = ParsePlaneRateFromReportHodoLT(r.reportfile, "P1Y", r.p1y_kHz);
    const bool gotP2X = ParsePlaneRateFromReportHodoLT(r.reportfile, "P2X", r.p2x_kHz);
    const bool gotP2Y = ParsePlaneRateFromReportHodoLT(r.reportfile, "P2Y", r.p2y_kHz);

    if (!(gotP1X && gotP1Y && gotP2X && gotP2Y)) {
      r.status = "MISSING";
      if (!gotP1X) r.status += "_P1X";
      if (!gotP1Y) r.status += "_P1Y";
      if (!gotP2X) r.status += "_P2X";
      if (!gotP2Y) r.status += "_P2Y";
      ++nMissingPlane;
      log << "WARNING [run " << r.run << "]: " << r.status << " in " << r.reportfile << "\n";
      rows.push_back(r);
      continue;
    }

    ComputeHodoLT(r);
    if (!std::isfinite(r.lt) || r.lt <= 0.0 || r.lt > 1.0 ||
        r.d1 < 0.0 || r.d2 < 0.0 || r.d3 < 0.0 || r.d4 < 0.0 ||
        r.d1 >= 1.0 || r.d2 >= 1.0 || r.d3 >= 1.0 || r.d4 >= 1.0) {
      r.status = "BAD_LT";
      ++nBadLT;
      log << "WARNING [run " << r.run << "]: invalid hodo3of4_LT=" << r.lt
          << " D=(" << r.d1 << "," << r.d2 << "," << r.d3 << "," << r.d4 << ")\n";
    } else {
      ++nOK;
    }

    rows.push_back(r);
  }

  log << "Rows: OK=" << nOK
      << " missingReport=" << nMissingReport
      << " missingPlane=" << nMissingPlane
      << " badLT=" << nBadLT << "\n";

  std::ofstream csv(outCsv.c_str());
  csv << std::setprecision(12);
  csv << "category,run,"
      << "P1X_rate_kHz,P1Y_rate_kHz,P2X_rate_kHz,P2Y_rate_kHz,"
      << "DPR_ns,"
      << "P1X_dead_prob,P1Y_dead_prob,P2X_dead_prob,P2Y_dead_prob,"
      << "P1X_live_prob,P1Y_live_prob,P2X_live_prob,P2Y_live_prob,"
      << "hodo3of4_LT,reportfile,status\n";

  for (const auto& r : rows) {
    csv << r.category << ","
        << r.run << ","
        << r.p1x_kHz << ","
        << r.p1y_kHz << ","
        << r.p2x_kHz << ","
        << r.p2y_kHz << ","
        << r.dpr_ns << ","
        << r.d1 << ","
        << r.d2 << ","
        << r.d3 << ","
        << r.d4 << ","
        << r.l1 << ","
        << r.l2 << ","
        << r.l3 << ","
        << r.l4 << ","
        << r.lt << ","
        << "\"" << r.reportfile << "\","
        << r.status << "\n";
  }
  csv.close();
  log << "Wrote CSV: " << outCsv << "\n";
  std::cout << "Wrote CSV: " << outCsv << "\n";
  log.close();
}
