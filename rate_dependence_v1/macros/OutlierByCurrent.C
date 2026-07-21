// OutlierByCurrent.C
//
// Zero-argument macro, intended to be run from rate_dependence_v1/:
//   root -l -b -q macros/OutlierByCurrent.C
//
// Scans:
//   results/pass4/**/tables/yield_vs_current_signal.csv
//   results/pass5/**/tables/yield_vs_current_signal.csv
//
// Conservative reporting:
//   1) outlier_by_pull
//      status == OK, fit_included == 1, yield_err > 0, and abs_pull > 3
//   2) excluded_from_fit
//      fit_included != 1, or status != OK, or invalid yield / yield_err inputs
//
// Computes when valid:
//   pull     = (yield - fit_const) / yield_err
//   abs_pull = |pull|
//
// Writes flagged rows to:
//   results/outliers/outliers_by_current.csv

#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

static std::string Trim(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace((unsigned char)s[b])) b++;
  while (e > b && std::isspace((unsigned char)s[e - 1])) e--;
  return s.substr(b, e - b);
}

static std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < line.size(); i++) {
    const char c = line[i];
    if (c == '"') {
      inq = !inq;
      continue;
    }
    if (c == ',' && !inq) {
      out.push_back(Trim(cur));
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(Trim(cur));
  return out;
}

static bool SafeToD(const std::string& s, double& out) {
  char* end = nullptr;
  const char* c = s.c_str();
  out = std::strtod(c, &end);
  if (c == end) return false;
  while (end && *end != '\0') {
    if (!std::isspace((unsigned char)*end)) return false;
    ++end;
  }
  return std::isfinite(out);
}

static bool SafeToI(const std::string& s, int& out) {
  char* end = nullptr;
  const char* c = s.c_str();
  long v = std::strtol(c, &end, 10);
  if (c == end) return false;
  while (end && *end != '\0') {
    if (!std::isspace((unsigned char)*end)) return false;
    ++end;
  }
  out = (int)v;
  return true;
}

static bool EnsureDir(const std::string& d) {
  return gSystem->mkdir(d.c_str(), true) == 0;
}

static bool StartsWith(const std::string& s, const std::string& prefix) {
  return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static void FindFilesRecursive(const std::string& dir,
                               const std::string& needle,
                               std::vector<std::string>& out) {
  DIR* dp = opendir(dir.c_str());
  if (!dp) return;

  while (dirent* ent = readdir(dp)) {
    const std::string name = ent->d_name;
    if (name == "." || name == "..") continue;

    const std::string path = dir + "/" + name;
    if (ent->d_type == DT_DIR) {
      FindFilesRecursive(path, needle, out);
    } else if (ent->d_type == DT_REG) {
      if (name == needle) out.push_back(path);
    } else {
      void* test = gSystem->OpenDirectory(path.c_str());
      if (test) {
        gSystem->FreeDirectory(test);
        FindFilesRecursive(path, needle, out);
      } else if (name == needle) {
        out.push_back(path);
      }
    }
  }
  closedir(dp);
}

struct OutlierRow {
  std::string setting;
  int run = -1;
  std::string flag_type;
  std::string reason;
  double BCM2_I = NAN;
  double yield = NAN;
  double yield_err = NAN;
  double fit_const = NAN;
  double abs_pull = NAN;
};

static bool ReadOutliersFromCSV(const std::string& csvPath,
                                const std::string& settingRel,
                                std::vector<OutlierRow>& out,
                                std::ostream& log) {
  std::ifstream in(csvPath.c_str());
  if (!in) {
    log << "WARN: cannot open " << csvPath << "\n";
    return false;
  }

  std::string line;
  std::string header;
  while (std::getline(in, line)) {
    const std::string t = Trim(line);
    if (t.empty()) continue;
    if (StartsWith(t, "#FIT")) continue;
    if (t[0] == '#') continue;
    header = t;
    break;
  }

  if (header.empty()) {
    log << "WARN: missing header in " << csvPath << "\n";
    return false;
  }

  const auto cols = SplitCSV(header);
  std::map<std::string, int> idx;
  for (int i = 0; i < (int)cols.size(); i++) idx[cols[i]] = i;

  auto need = [&](const std::string& name) -> int {
    auto it = idx.find(name);
    return (it == idx.end()) ? -1 : it->second;
  };

  const int iRun = need("run");
  const int iI = need("BCM2_I");
  const int iYield = need("yield");
  const int iYerr = need("yield_err");
  const int iStatus = need("status");
  const int iFitIncluded = need("fit_included");
  const int iFitExclReason = need("fit_excluded_reason");

  if (iRun < 0 || iI < 0 || iYield < 0 || iYerr < 0 || iStatus < 0 || iFitIncluded < 0) {
    log << "WARN: missing required columns in " << csvPath << "\n";
    return false;
  }

  std::vector<std::vector<std::string>> rows;
  double fit_const = NAN;

  while (std::getline(in, line)) {
    const std::string t = Trim(line);
    if (t.empty()) continue;

    if (StartsWith(t, "#FIT")) {
      const auto fitFields = SplitCSV(t);
      if (fitFields.size() > 1) {
        double tmp = NAN;
        if (SafeToD(fitFields[1], tmp)) fit_const = tmp;
      }
      continue;
    }

    if (t[0] == '#') continue;

    const auto fields = SplitCSV(t);
    if ((int)fields.size() < (int)cols.size()) {
      log << "WARN: short line in " << csvPath << ": " << t << "\n";
      continue;
    }
    rows.push_back(fields);
  }

  if (!std::isfinite(fit_const)) {
    log << "WARN: no numeric fit_const found in " << csvPath << "\n";
    return false;
  }

  int kept = 0;
  for (const auto& r : rows) {
    int run = -1;
    double bcm2_i = NAN, yy = NAN, yerr = NAN;
    const bool hasRun = SafeToI(r[iRun], run);
    const bool hasI = SafeToD(r[iI], bcm2_i);
    const bool hasY = SafeToD(r[iYield], yy);
    const bool hasYerr = SafeToD(r[iYerr], yerr);

    const std::string status = r[iStatus];
    const std::string fitIncluded = r[iFitIncluded];
    const std::string fitExcl = (iFitExclReason >= 0 && iFitExclReason < (int)r.size()) ? r[iFitExclReason] : "";

    const bool statusOK = (status == "OK");
    const bool fitIncludedOK = (fitIncluded == "1");
    const bool yieldValid = hasY && std::isfinite(yy);
    const bool yerrValid = hasYerr && std::isfinite(yerr) && (yerr > 0.0);

    OutlierRow o;
    o.setting = settingRel;
    o.run = hasRun ? run : -1;
    o.BCM2_I = hasI ? bcm2_i : NAN;
    o.yield = yieldValid ? yy : NAN;
    o.yield_err = hasYerr ? yerr : NAN;
    o.fit_const = fit_const;
    o.abs_pull = NAN;

    if (!statusOK || !fitIncludedOK || !yieldValid || !yerrValid) {
      o.flag_type = "excluded_from_fit";
      if (!Trim(fitExcl).empty()) o.reason = Trim(fitExcl);
      else if (!statusOK) o.reason = status;
      else if (!fitIncludedOK) o.reason = "FIT_INCLUDED_0";
      else if (!yieldValid) o.reason = "BAD_YIELD";
      else if (!yerrValid) o.reason = "BAD_YERR";
      else o.reason = "EXCLUDED";
      out.push_back(o);
      kept++;
      continue;
    }

    const double pull = (yy - fit_const) / yerr;
    const double abs_pull = std::fabs(pull);
    if (abs_pull > 3.0) {
      o.flag_type = "outlier_by_pull";
      o.reason = "abs_pull>3";
      o.abs_pull = abs_pull;
      out.push_back(o);
      kept++;
    }
  }

  log << "Scanned " << csvPath << " -> flagged rows kept: " << kept << "\n";
  return true;
}

}  // namespace

void OutlierByCurrent() {
  const std::string resultsDir = "results";
  const std::string outDir = resultsDir + "/outliers";
  const std::string outCsv = outDir + "/outliers_by_current.csv";

  EnsureDir(outDir);

  std::vector<std::string> csvFiles;
  FindFilesRecursive(resultsDir + "/pass4", "yield_vs_current_signal.csv", csvFiles);
  FindFilesRecursive(resultsDir + "/pass5", "yield_vs_current_signal.csv", csvFiles);
  std::sort(csvFiles.begin(), csvFiles.end());

  std::vector<OutlierRow> all;
  std::ostringstream log;
  log << "OutlierByCurrent\n";
  log << "Results dir: " << resultsDir << "\n";
  log << "Selection:\n";
  log << "  outlier_by_pull => status==OK && fit_included==1 && yield_err>0 && abs((yield-fit_const)/yield_err)>3\n";
  log << "  excluded_from_fit => fit_included!=1 or status!=OK or invalid yield/yield_err\n";
  log << "CSV files found: " << csvFiles.size() << "\n";

  for (const auto& csvPath : csvFiles) {
    std::string settingRel = csvPath;
    const std::string prefix = resultsDir + "/";
    if (StartsWith(settingRel, prefix)) settingRel.erase(0, prefix.size());
    const std::string suffix = "/tables/yield_vs_current_signal.csv";
    if (settingRel.size() >= suffix.size() &&
        settingRel.compare(settingRel.size() - suffix.size(), suffix.size(), suffix) == 0) {
      settingRel.erase(settingRel.size() - suffix.size());
    }
    ReadOutliersFromCSV(csvPath, settingRel, all, log);
  }

  std::sort(all.begin(), all.end(), [](const OutlierRow& a, const OutlierRow& b) {
    if (a.flag_type != b.flag_type) return a.flag_type < b.flag_type;
    if (a.flag_type == "outlier_by_pull" && b.flag_type == "outlier_by_pull" && a.abs_pull != b.abs_pull) {
      return a.abs_pull > b.abs_pull;
    }
    if (a.reason != b.reason) return a.reason < b.reason;
    if (a.setting != b.setting) return a.setting < b.setting;
    return a.run < b.run;
  });

  std::ofstream out(outCsv.c_str());
  if (!out) {
    std::cerr << "ERROR: cannot write " << outCsv << "\n";
    std::cerr << log.str();
    return;
  }

  out << "# OutlierByCurrent output\n";
  out << "# source = results/pass4/**/tables/yield_vs_current_signal.csv ; results/pass5/**/tables/yield_vs_current_signal.csv\n";
  out << "# flag_type=outlier_by_pull => status==OK && fit_included==1 && yield_err>0 && abs((yield-fit_const)/yield_err)>3\n";
  out << "# flag_type=excluded_from_fit => fit_included!=1 or status!=OK or invalid yield/yield_err\n";
  out << "# threshold = abs_pull > 3\n";
  out << "# sort = flag_type, then abs_pull descending within outlier_by_pull, then reason/setting/run\n";
  out << "setting,run,flag_type,reason,BCM2_I,yield,yield_err,fit_const,abs_pull\n";
  out << std::fixed << std::setprecision(6);
  for (const auto& r : all) {
    out << '"' << r.setting << '"' << ","
        << r.run << ","
        << '"' << r.flag_type << '"' << ","
        << '"' << r.reason << '"' << ","
        << r.BCM2_I << ","
        << r.yield << ","
        << r.yield_err << ","
        << r.fit_const << ","
        << r.abs_pull << "\n";
  }
  out.close();

  std::cout << log.str();
  std::cout << "Wrote " << all.size() << " flagged rows to " << outCsv << "\n";
}
