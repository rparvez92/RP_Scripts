// YieldVsTrigger.C
// v1: Signal-only Yield vs SHMS 3/4 Trigger Rate (kHz), constant fit (pol0)
// Inputs:
//   YieldVsTrigger("<abs path>/settings/.../manifest.json", "<abs path>/results")
//
// Reads:
//   <results>/<same rel path as manifest>/tables/yield_vs_current_signal.csv
// Report files:
//   <repo_root>/Pass0p1_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report
//
// Writes:
//   <results>/<rel>/tables/yield_vs_trigger_shms34.csv
//   <results>/<rel>/PNGs/yield_vs_trigger_shms34.png
//   <results>/<rel>/canvases/yield_vs_trigger_shms34.root
//   <results>/<rel>/logs/YieldVsTrigger.log

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>
#include <TFile.h>

#include <algorithm>
#include <cctype>
#include <cmath>
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

// Minimal CSV splitter that respects quotes
static std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < line.size(); i++) {
    char c = line[i];
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

static bool IsFinite(double x) { return std::isfinite(x); }

static bool EnsureDir(const std::string& d) {
  // ROOT's gSystem is available in ROOT macros
  return (gSystem->mkdir(d.c_str(), true) == 0);
}

static std::string Dirname(const std::string& path) {
  auto pos = path.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return path.substr(0, pos);
}

// Build results subdir by taking manifest path relative to ".../settings/"
static std::string ComputeRelFromManifest(const std::string& manifestPath) {
  // Expect: .../settings/<REL>/manifest.json
  const std::string key = "/settings/";
  auto kpos = manifestPath.find(key);
  std::string rel;
  if (kpos != std::string::npos) {
    rel = manifestPath.substr(kpos + key.size());
  } else {
    // fallback: use last ~7 directories
    rel = manifestPath;
  }
  // strip trailing "/manifest.json"
  const std::string suf = "/manifest.json";
  if (rel.size() >= suf.size() && rel.compare(rel.size() - suf.size(), suf.size(), suf) == 0) {
    rel = rel.substr(0, rel.size() - suf.size());
  } else {
    // strip filename
    rel = Dirname(rel);
  }
  return rel;
}

static std::vector<std::string> SplitPath(const std::string& rel) {
  std::vector<std::string> parts;
  std::stringstream ss(rel);
  std::string item;
  while (std::getline(ss, item, '/')) {
    if (!item.empty()) parts.push_back(item);
  }
  return parts;
}

static std::string JoinSettingHeaderLine1(const std::vector<std::string>& parts) {
  // Expect: pass/run_type/target/z/x/Q2/setting_id
  // We want a compact readable line:
  // pass4 / pi+sidis / LH2 / z0p36 / x0p25 / Q23p3
  std::string s;
  if (parts.size() >= 6) {
    s = parts[0] + " / " + parts[1] + " / " + parts[2] + " / " + parts[3] + " / " + parts[4] + " / " + parts[5];
  } else {
    s = "setting: " + (parts.empty() ? std::string("unknown") : parts[0]);
  }
  return s;
}

static std::string JoinSettingHeaderLine2(const std::vector<std::string>& parts) {
  if (parts.size() >= 7) return parts[6];
  return "";
}

static bool ParseReport_SHMS34_kHz(const std::string& reportPath, double& out) {
  std::ifstream in(reportPath);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.find("SHMS 3/4 Trigger Rate") != std::string::npos && line.find("kHz") != std::string::npos) {
      // example: "SHMS 3/4 Trigger Rate         : 864.802 kHz"
      // find ':', then parse double
      auto cpos = line.find(':');
      if (cpos == std::string::npos) continue;
      std::string rhs = line.substr(cpos + 1);
      rhs = Trim(rhs);
      // rhs begins with number
      std::stringstream ss(rhs);
      double val;
      ss >> val;
      if (ss.fail()) continue;
      out = val;
      return true;
    }
  }
  return false;
}

static bool ParseReport_SHMS_ELCLEAN_kHz(const std::string& reportPath, double& out) {
  std::ifstream in(reportPath);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    // example: "SHMS_pEL_CLEAN :    11827351    [ 6.555 kHz ]"
    if (line.find("SHMS_pEL_CLEAN") != std::string::npos && line.find('[') != std::string::npos &&
        line.find("kHz") != std::string::npos) {
      auto lb = line.find('[');
      auto rb = line.find(']', lb == std::string::npos ? 0 : lb);
      if (lb == std::string::npos || rb == std::string::npos || rb <= lb) continue;
      std::string inside = line.substr(lb + 1, rb - lb - 1); // " 6.555 kHz "
      inside = Trim(inside);
      std::stringstream ss(inside);
      double val;
      ss >> val;
      if (ss.fail()) continue;
      out = val;
      return true;
    }
  }
  return false;
}

struct Row {
  int run = -1;
  double current = NAN;      // BCM2_I
  double yield = NAN;        // myCTime normalized yield
  double yerr = NAN;         // yield_err
  double shms34 = NAN;       // kHz
  double shms_elclean = NAN; // kHz
  std::string status;        // OK / skipped reason
  std::string rootfile;
};

static bool ReadYieldVsCurrentSignalCSV(const std::string& csvPath, std::vector<Row>& rows, std::ostream& log) {
  std::ifstream in(csvPath);
  if (!in) {
    log << "ERROR: cannot open input CSV: " << csvPath << "\n";
    return false;
  }

  std::string header;
  if (!std::getline(in, header)) {
    log << "ERROR: empty CSV: " << csvPath << "\n";
    return false;
  }
  auto cols = SplitCSV(header);
  std::map<std::string, int> idx;
  for (int i = 0; i < (int)cols.size(); i++) idx[cols[i]] = i;

  auto need = [&](const std::string& name) -> int {
    auto it = idx.find(name);
    return (it == idx.end()) ? -1 : it->second;
  };

  const int iRun = need("run");
  const int iI = need("BCM2_I");
  const int iY = need("yield");
  const int iYe = need("yield_err");
  const int iStatus = need("status");
  const int iRoot = need("rootfile");

  if (iRun < 0 || iI < 0 || iY < 0 || iYe < 0) {
    log << "ERROR: CSV missing required columns. Need run, BCM2_I, yield, yield_err.\n";
    log << "Header was: " << header << "\n";
    return false;
  }

  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    if (Trim(line).empty()) continue;
    auto v = SplitCSV(line);
    if ((int)v.size() < (int)cols.size()) {
      // allow short lines, but skip
      log << "WARN: skipping short CSV line: " << line << "\n";
      continue;
    }

    Row r;
    try {
      r.run = std::stoi(v[iRun]);
    } catch (...) {
      log << "WARN: bad run value, skipping line: " << line << "\n";
      continue;
    }
    try {
      r.current = std::stod(v[iI]);
      r.yield = std::stod(v[iY]);
      r.yerr = std::stod(v[iYe]);
    } catch (...) {
      r.status = "BAD_NUMERIC";
    }
    if (iStatus >= 0) r.status = v[iStatus];
    if (iRoot >= 0) r.rootfile = v[iRoot];

    rows.push_back(r);
    n++;
  }

  log << "Read " << n << " rows from: " << csvPath << "\n";
  return true;
}

// Try to locate repo root from manifest path: go up until ".../rate_dependence_v1" directory exists in string
static std::string GuessRepoRootFromManifest(const std::string& manifestPath) {
  // Common case: .../rate_dependence_v1/settings/...
  const std::string key = "/rate_dependence_v1/";
  auto pos = manifestPath.find(key);
  if (pos != std::string::npos) {
    return manifestPath.substr(0, pos + key.size() - 1); // strip trailing '/'
  }
  // fallback: just use cwd
  return gSystem->WorkingDirectory();
}

} // namespace

void YieldVsTrigger(const char* manifest_json, const char* results_top) {
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const std::string manifestPath = manifest_json;
  const std::string resultsTop = results_top;

  const std::string rel = ComputeRelFromManifest(manifestPath);
  const std::string outBase = resultsTop + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTables = outBase + "/tables";
  const std::string outCanv = outBase + "/canvases";
  const std::string outLogs = outBase + "/logs";

  EnsureDir(outBase);
  EnsureDir(outPNGs);
  EnsureDir(outTables);
  EnsureDir(outCanv);
  EnsureDir(outLogs);

  const std::string logPath = outLogs + "/YieldVsTrigger.log";
  std::ofstream log(logPath.c_str());
  log << "YieldVsTrigger.C\n";
  log << "manifest: " << manifestPath << "\n";
  log << "resultsTop: " << resultsTop << "\n";
  log << "rel: " << rel << "\n";

  // Input CSV from YieldVsCurrent
  const std::string inCSV = outTables + "/yield_vs_current_signal.csv";
  std::vector<Row> rows;
  if (!ReadYieldVsCurrentSignalCSV(inCSV, rows, log)) {
    std::cerr << "ERROR: failed reading " << inCSV << "\n";
    log.close();
    return;
  }

  // Report directory
  const std::string repoRoot = GuessRepoRootFromManifest(manifestPath);
  const std::string reportDir = repoRoot + "/Pass0p1_REPORTfiles/COIN/PRODUCTION";

  // Parse SHMS 3/4 trigger rate for each run
  int nOK = 0, nSkip = 0, nBadErr = 0, nNoReport = 0, nNoLine = 0;
  for (auto& r : rows) {
    // Validate yield entry
    if (!IsFinite(r.current) || !IsFinite(r.yield) || !IsFinite(r.yerr)) {
      r.status = "NONFINITE";
      nSkip++;
      continue;
    }
    if (r.yerr <= 0) {
      r.status = "BAD_ERR";
      nBadErr++;
      continue;
    }

    const std::string rep = reportDir + "/replay_coin_production_" + std::to_string(r.run) + "_-1.report";
    double sh34 = NAN;
    if (!ParseReport_SHMS34_kHz(rep, sh34)) {
      // Distinguish missing file vs missing line
      std::ifstream test(rep);
      if (!test) {
        r.status = "MISSING_REPORT";
        nNoReport++;
      } else {
        r.status = "MISSING_SHMS34_LINE";
        nNoLine++;
      }
      nSkip++;
      continue;
    }
    r.shms34 = sh34;
    r.status = "OK";
    nOK++;
  }

  log << "Report dir: " << reportDir << "\n";
  log << "Points: OK=" << nOK << "  skipped=" << nSkip << "  badErr=" << nBadErr
      << "  missingReport=" << nNoReport << "  missingLine=" << nNoLine << "\n";

  // Write table (include skipped rows too, but mark status)
  const std::string outCSV = outTables + "/yield_vs_trigger_shms34.csv";
  std::ofstream csv(outCSV.c_str());
  csv << "category,run,BCM2_I,yield,yield_err,shms34_rate_kHz,rootfile,status\n";
  for (const auto& r : rows) {
    csv << "signal,"
        << r.run << ","
        << std::setprecision(10) << r.current << ","
        << std::setprecision(10) << r.yield << ","
        << std::setprecision(10) << r.yerr << ","
        << std::setprecision(10) << (IsFinite(r.shms34) ? r.shms34 : NAN) << ","
        << "\"" << r.rootfile << "\"" << ","
        << r.status << "\n";
  }

  // Build graph from OK points
  std::vector<double> xs, ys, exs, eys;
  xs.reserve(rows.size());
  ys.reserve(rows.size());
  exs.reserve(rows.size());
  eys.reserve(rows.size());

  double xmin = +1e99, xmax = -1e99, ymin = +1e99, ymax = -1e99;
  for (const auto& r : rows) {
    if (r.status != "OK") continue;
    if (!IsFinite(r.shms34)) continue;
    xs.push_back(r.shms34);
    ys.push_back(r.yield);
    exs.push_back(0.0);
    eys.push_back(r.yerr);

    xmin = std::min(xmin, r.shms34);
    xmax = std::max(xmax, r.shms34);
    ymin = std::min(ymin, r.yield - r.yerr);
    ymax = std::max(ymax, r.yield + r.yerr);
  }

  const std::string outPng = outPNGs + "/yield_vs_trigger_shms34.png";
  const std::string outRoot = outCanv + "/yield_vs_trigger_shms34.root";

  TCanvas c("c_yield_vs_trigger_shms34", "Yield vs SHMS 3/4 Trigger Rate", 1200, 850);
  c.SetTopMargin(0.22);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.06);
  c.SetBottomMargin(0.12);

  if (xs.empty()) {
    log << "WARN: No points to plot.\n";
    // Save empty
    c.SaveAs(outPng.c_str());
    TFile fout(outRoot.c_str(), "RECREATE");
    c.Write("c_yield_vs_trigger_shms34");
    fout.Close();
    csv << "# FIT,NO_POINTS\n";
    csv.close();
    std::cout << "Wrote CSV: " << outCSV << "\n";
    std::cout << "Wrote ROOT: " << outRoot << "\n";
    log << "Wrote PNG (empty): " << outPng << "\n";
    log << "Wrote ROOT: " << outRoot << "\n";
    log << "Wrote CSV: " << outCSV << "\n";
    log.close();
    return;
  }

  // Add padding to axis ranges
  const double xpad = 0.08 * (xmax - xmin);
  const double ypad = 0.12 * (ymax - ymin);
  const double xlo = xmin - xpad;
  const double xhi = xmax + xpad;
  const double ylo = ymin - ypad;
  const double yhi = ymax + ypad;

  TGraphErrors g((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
  g.SetName("g_yield_vs_trigger_shms34");
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.35);
  g.SetMarkerColor(kRed + 1);
  g.SetLineColor(kRed + 1);

  // Frame via Draw("AP") then set axis titles
  g.Draw("AP");
  g.GetXaxis()->SetTitle("SHMS 3/4 Trigger Rate (kHz)");
  g.GetYaxis()->SetTitle("Normalized Yield (myCTime)");
  g.GetXaxis()->SetTitleSize(0.045);
  g.GetYaxis()->SetTitleSize(0.045);
  g.GetXaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetTitleOffset(1.25);
  g.GetXaxis()->SetLimits(xlo, xhi);
  g.SetMinimum(ylo);
  g.SetMaximum(yhi);

  // Constant fit
  TF1 f("f_const", "pol0", xlo, xhi);
  f.SetLineColor(kBlack);
  f.SetLineStyle(2); // dashed
  f.SetLineWidth(2);

  // Fit quietly (no spam), but still compute stats
  g.Fit(&f, "Q"); // Q = quiet

  const double C = f.GetParameter(0);
  const double Cerr = f.GetParError(0);
  const double chi2 = f.GetChisquare();
  const int ndf = f.GetNDF();
  const double prob = f.GetProb();
  const double chi2ndf = (ndf > 0 ? chi2 / (double)ndf : NAN);

  // Legend in top margin
  TLegend leg(0.12, 0.80, 0.45, 0.95);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(&g, "signal", "p");
  leg.AddEntry(&f, Form("const fit: C=%.3g", C), "l");
  leg.Draw();

  // Setting info in top margin (right side)
  auto parts = SplitPath(rel);
  const std::string line1 = JoinSettingHeaderLine1(parts);
  const std::string line2 = JoinSettingHeaderLine2(parts);

  TLatex tl;
  tl.SetNDC(true);
  tl.SetTextFont(42);
  tl.SetTextSize(0.03);
  tl.SetTextAlign(13); // left-top
  tl.DrawLatex(0.35, 0.94, line1.c_str());
  if (!line2.empty()) tl.DrawLatex(0.35, 0.89, line2.c_str());

  // Fit stats text (colored) inside plot
  TLatex chiText;
  chiText.SetNDC(true);
  chiText.SetTextFont(42);
  chiText.SetTextColor(kBlue + 2);
  chiText.SetTextSize(0.040);
  chiText.SetTextAlign(22); // center
  chiText.DrawLatex(0.52, 0.30, Form("#chi^{2}/ndf=%.2f (%d), Prob=%.3g", chi2ndf, ndf, prob));

  // Write fit summary into CSV
  csv << "# FIT_MODEL,pol0\n";
  csv << "# FIT_C," << std::setprecision(12) << C << "\n";
  csv << "# FIT_Cerr," << std::setprecision(12) << Cerr << "\n";
  csv << "# FIT_chi2," << std::setprecision(12) << chi2 << "\n";
  csv << "# FIT_ndf," << ndf << "\n";
  csv << "# FIT_chi2ndf," << std::setprecision(12) << chi2ndf << "\n";
  csv << "# FIT_prob," << std::setprecision(12) << prob << "\n";
  csv.close();

  // Save outputs
  c.SaveAs(outPng.c_str());
  TFile fout(outRoot.c_str(), "RECREATE");
  c.Write("c_yield_vs_trigger_shms34");
  g.Write("g_yield_vs_trigger_shms34");
  f.Write("f_const");
  fout.Close();

  log << "Fit: C=" << C << " +/- " << Cerr << "  chi2/ndf=" << chi2ndf << " (" << ndf << ")  prob=" << prob
      << "\n";
  log << "Wrote CSV: " << outCSV << "\n";
  std::cout << "Wrote CSV: " << outCSV << "\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT: " << outRoot << "\n";
  std::cout << "Wrote ROOT: " << outRoot << "\n";
  log.close();
}
