// TriggerVsCurrent.C
// HMS_hEL_CLEAN (total counts) / BCM2_Q (mC) vs BCM2_I (current), constant fit (pol0)
// Includes Poisson statistical uncertainties: sigma_y = sqrt(N)/Q
//
// Purpose:
//   Rate-dependence diagnostic using a "signal-like" trigger (HMS_hEL_CLEAN) normalized by beam charge.
//
// Inputs:
//   TriggerVsCurrent("<abs path>/settings/.../manifest.json", "<abs path>/results")
//
// Reads:
//   <results>/<same rel path as manifest>/tables/yield_vs_current_signal.csv
// Bigtable:
//   <repo_root>/bigtable/rsidis_bigtable_pass0p1.csv  (columns: run, BCM2_Q)
//
// Report files:
//   <repo_root>/Pass0p1_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report
//
// Extracts from report line like:
//   HMS_hEL_CLEAN : 590009    [ 0.327 kHz ]
// Uses the total counts (590009) as N.
//
// Computes:
//   y = N / Q   (counts/mC)
//   sigma_y = sqrt(N) / Q    (Poisson, charge uncertainty neglected)
//
// Writes:
//   <results>/<rel>/tables/hms_elclean_chargeNorm_vs_current.csv
//   <results>/<rel>/PNGs/hms_elclean_chargeNorm_vs_current.png
//   <results>/<rel>/canvases/hms_elclean_chargeNorm_vs_current.root
//   <results>/<rel>/logs/TriggerVsCurrent.log

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TF1.h>
#include <TFile.h>
#include <TSystem.h>
#include <TAxis.h>

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

static std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i = 0; i < line.size(); i++) {
    char c = line[i];
    if (c == '"') { inq = !inq; continue; }
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

static bool EnsureDir(const std::string& d) { return (gSystem->mkdir(d.c_str(), true) == 0); }

static std::string Dirname(const std::string& path) {
  auto pos = path.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return path.substr(0, pos);
}

static std::string ComputeRelFromManifest(const std::string& manifestPath) {
  const std::string key = "/settings/";
  auto kpos = manifestPath.find(key);
  std::string rel;
  if (kpos != std::string::npos) rel = manifestPath.substr(kpos + key.size());
  else rel = manifestPath;

  const std::string suf = "/manifest.json";
  if (rel.size() >= suf.size() && rel.compare(rel.size() - suf.size(), suf.size(), suf) == 0) {
    rel = rel.substr(0, rel.size() - suf.size());
  } else {
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

// Parse HMS_hEL_CLEAN total counts from report line:
//   HMS_hEL_CLEAN : 590009    [ 0.327 kHz ]
static bool ParseReport_HMS_ELCLEAN_counts(const std::string& reportPath, double& outCounts) {
  std::ifstream in(reportPath);
  if (!in) return false;

  std::string line;
  while (std::getline(in, line)) {
    if (line.find("HMS_hEL_CLEAN") == std::string::npos) continue;

    auto cpos = line.find(':');
    auto lb   = line.find('[');
    if (cpos == std::string::npos || lb == std::string::npos || lb <= cpos) continue;

    std::string between = Trim(line.substr(cpos + 1, lb - (cpos + 1))); // "590009    "
    if (between.empty()) continue;

    std::stringstream ss(between);
    long long n = 0;
    ss >> n;
    if (ss.fail()) continue;

    outCounts = (double)n;
    return true;
  }
  return false;
}

struct Row {
  int run = -1;
  double current = NAN;     // BCM2_I
  double yield = NAN;       // for context only
  double yerr = NAN;        // for context only
  double bcm2_q = NAN;      // mC, from bigtable

  double hms_elclean_counts = NAN;       // N
  double hms_elclean_cnorm  = NAN;       // y = N/Q (counts/mC)
  double hms_elclean_cnorm_err = NAN;    // sigma_y = sqrt(N)/Q

  std::string status;
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
  const int iI   = need("BCM2_I");
  const int iY   = need("yield");
  const int iYe  = need("yield_err");
  const int iStatus = need("status");
  const int iRoot   = need("rootfile");

  if (iRun < 0 || iI < 0) {
    log << "ERROR: CSV missing required columns. Need run, BCM2_I.\n";
    log << "Header was: " << header << "\n";
    return false;
  }

  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    if (Trim(line).empty()) continue;
    auto v = SplitCSV(line);
    if ((int)v.size() < (int)cols.size()) {
      log << "WARN: skipping short CSV line: " << line << "\n";
      continue;
    }

    Row r;
    try { r.run = std::stoi(v[iRun]); }
    catch (...) { log << "WARN: bad run value, skipping line: " << line << "\n"; continue; }

    try { r.current = std::stod(v[iI]); }
    catch (...) { r.status = "BAD_CURRENT"; }

    if (iY >= 0)  { try { r.yield = std::stod(v[iY]); }  catch (...) {} }
    if (iYe >= 0) { try { r.yerr  = std::stod(v[iYe]); } catch (...) {} }

    if (iStatus >= 0) r.status = v[iStatus];
    if (iRoot   >= 0) r.rootfile = v[iRoot];

    rows.push_back(r);
    n++;
  }

  log << "Read " << n << " rows from: " << csvPath << "\n";
  return true;
}

// Read run -> BCM2_Q from bigtable CSV
static bool ReadBigtableBCM2Q(const std::string& bigtableCSV,
                             std::map<int,double>& run2q,
                             std::ostream& log)
{
  std::ifstream in(bigtableCSV);
  if (!in) {
    log << "ERROR: cannot open bigtable CSV: " << bigtableCSV << "\n";
    return false;
  }
  std::string header;
  if (!std::getline(in, header)) {
    log << "ERROR: empty bigtable CSV: " << bigtableCSV << "\n";
    return false;
  }

  auto cols = SplitCSV(header);
  std::map<std::string, int> idx;
  for (int i=0;i<(int)cols.size();++i) idx[cols[i]] = i;

  auto need = [&](const std::string& name)->int{
    auto it = idx.find(name);
    return (it==idx.end()) ? -1 : it->second;
  };

  const int iRun = need("run");
  const int iQ   = need("BCM2_Q");
  if (iRun < 0 || iQ < 0) {
    log << "ERROR: bigtable missing required columns: run, BCM2_Q\n";
    log << "Header was: " << header << "\n";
    return false;
  }

  int n=0, nBad=0;
  std::string line;
  while (std::getline(in, line)) {
    if (Trim(line).empty()) continue;
    auto v = SplitCSV(line);
    if ((int)v.size() < std::max(iRun,iQ)+1) { nBad++; continue; }
    int run=-1;
    double q=NAN;
    try { run = std::stoi(v[iRun]); } catch (...) { nBad++; continue; }
    try { q = std::stod(v[iQ]); } catch (...) { nBad++; continue; }
    if (!IsFinite(q)) { nBad++; continue; }
    run2q[run] = q;
    n++;
  }
  log << "Loaded BCM2_Q for " << n << " runs from bigtable. Bad lines: " << nBad << "\n";
  return true;
}

static std::string GuessRepoRootFromManifest(const std::string& manifestPath) {
  const std::string key = "/rate_dependence_v1/";
  auto pos = manifestPath.find(key);
  if (pos != std::string::npos) {
    return manifestPath.substr(0, pos + key.size() - 1); // path ending at .../rate_dependence_v1
  }
  return gSystem->WorkingDirectory();
}

} // namespace

void TriggerVsCurrent(const char* manifest_json, const char* results_top) {
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  const std::string manifestPath = manifest_json ? manifest_json : "";
  const std::string resultsTop   = results_top   ? results_top   : "";

  const std::string rel = ComputeRelFromManifest(manifestPath);
  const std::string outBase = resultsTop + "/" + rel;

  const std::string outPNGs  = outBase + "/PNGs";
  const std::string outTables= outBase + "/tables";
  const std::string outCanv  = outBase + "/canvases";
  const std::string outLogs  = outBase + "/logs";

  EnsureDir(outBase);
  EnsureDir(outPNGs);
  EnsureDir(outTables);
  EnsureDir(outCanv);
  EnsureDir(outLogs);

  const std::string logPath = outLogs + "/TriggerVsCurrent.log";
  std::ofstream log(logPath.c_str());
  log << "TriggerVsCurrent.C\n";
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

  // Bigtable run -> BCM2_Q
  const std::string repoRoot = GuessRepoRootFromManifest(manifestPath);
  const std::string bigtableCSV = repoRoot + "/bigtable/rsidis_bigtable_pass0p1.csv";
  std::map<int,double> run2q;
  if (!ReadBigtableBCM2Q(bigtableCSV, run2q, log)) {
    std::cerr << "ERROR: failed reading bigtable BCM2_Q: " << bigtableCSV << "\n";
    log.close();
    return;
  }

  // Report dir
  const std::string reportDir = repoRoot + "/Pass0p1_REPORTfiles/COIN/PRODUCTION";

  int nOK = 0, nSkip = 0, nNoReport = 0, nNoLine = 0, nNoQ = 0, nBadQ = 0, nBadN = 0;

  for (auto& r : rows) {
    if (!IsFinite(r.current)) {
      r.status = "NONFINITE_CURRENT";
      nSkip++;
      continue;
    }

    // charge from bigtable
    auto itq = run2q.find(r.run);
    if (itq == run2q.end()) {
      r.status = "MISSING_BCM2_Q";
      nNoQ++;
      nSkip++;
      continue;
    }
    r.bcm2_q = itq->second;
    if (!IsFinite(r.bcm2_q) || r.bcm2_q <= 0) {
      r.status = "BAD_BCM2_Q";
      nBadQ++;
      nSkip++;
      continue;
    }

    // report parse
    const std::string rep = reportDir + "/replay_coin_production_" + std::to_string(r.run) + "_-1.report";
    double counts = NAN;
    if (!ParseReport_HMS_ELCLEAN_counts(rep, counts)) {
      std::ifstream test(rep);
      if (!test) { r.status = "MISSING_REPORT"; nNoReport++; }
      else       { r.status = "MISSING_ELCLEAN_LINE"; nNoLine++; }
      nSkip++;
      continue;
    }
    if (!IsFinite(counts) || counts < 0) {
      r.status = "BAD_COUNTS";
      nBadN++;
      nSkip++;
      continue;
    }

    r.hms_elclean_counts = counts;

    // y = N/Q and sigma_y = sqrt(N)/Q (Poisson)
    r.hms_elclean_cnorm = counts / r.bcm2_q;
    r.hms_elclean_cnorm_err = (counts > 0 ? std::sqrt(counts) / r.bcm2_q : 0.0);

    r.status = "OK";
    nOK++;
  }

  log << "Repo root:   " << repoRoot << "\n";
  log << "Bigtable:    " << bigtableCSV << "\n";
  log << "Report dir:  " << reportDir << "\n";
  log << "Points: OK=" << nOK
      << "  skipped=" << nSkip
      << "  missingReport=" << nNoReport
      << "  missingLine=" << nNoLine
      << "  missingBCM2_Q=" << nNoQ
      << "  badBCM2_Q=" << nBadQ
      << "  badCounts=" << nBadN
      << "\n";

  // Output naming
  const std::string stem = "hms_elclean_chargeNorm_vs_current";
  const std::string outCSV  = outTables + "/" + stem + ".csv";
  const std::string outPng  = outPNGs   + "/" + stem + ".png";
  const std::string outRoot = outCanv   + "/" + stem + ".root";

  // Write table (include all rows)
  std::ofstream csv(outCSV.c_str());
  csv << "category,run,BCM2_I,BCM2_Q_mC,hms_elclean_counts,"
         "hms_elclean_counts_per_mC,hms_elclean_counts_per_mC_err,"
         "yield,yield_err,rootfile,status\n";
  for (const auto& r : rows) {
    csv << "signal,"
        << r.run << ","
        << std::setprecision(10) << r.current << ","
        << std::setprecision(10) << (IsFinite(r.bcm2_q) ? r.bcm2_q : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.hms_elclean_counts) ? r.hms_elclean_counts : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.hms_elclean_cnorm)  ? r.hms_elclean_cnorm  : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.hms_elclean_cnorm_err) ? r.hms_elclean_cnorm_err : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.yield) ? r.yield : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.yerr)  ? r.yerr  : NAN) << ","
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
    if (!IsFinite(r.current) || !IsFinite(r.hms_elclean_cnorm) || !IsFinite(r.hms_elclean_cnorm_err)) continue;

    xs.push_back(r.current);
    ys.push_back(r.hms_elclean_cnorm);
    exs.push_back(0.0);
    eys.push_back(r.hms_elclean_cnorm_err);

    xmin = std::min(xmin, r.current);
    xmax = std::max(xmax, r.current);
    ymin = std::min(ymin, r.hms_elclean_cnorm);
    ymax = std::max(ymax, r.hms_elclean_cnorm);
  }

  TCanvas c("c_hms_elclean_chargeNorm_vs_current",
            "HMS EL_CLEAN (counts/charge) vs Current", 1200, 850);
  c.SetTopMargin(0.22);
  c.SetLeftMargin(0.18);   // prevents y-title cutoff
  c.SetRightMargin(0.06);
  c.SetBottomMargin(0.12);

  if (xs.empty()) {
    log << "WARN: No points to plot.\n";
    c.SaveAs(outPng.c_str());
    TFile fout(outRoot.c_str(), "RECREATE");
    c.Write("c_hms_elclean_chargeNorm_vs_current");
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

  const double xpad = 0.10 * (xmax - xmin);
  const double ypad = 0.20 * (ymax - ymin);
  const double xlo = xmin - xpad;
  const double xhi = xmax + xpad;
  const double ylo = ymin - ypad;
  const double yhi = ymax + ypad;

  TGraphErrors g((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
  g.SetName("g_hms_elclean_chargeNorm_vs_current");
  g.SetMarkerStyle(21);
  g.SetMarkerSize(1.35);
  g.SetMarkerColor(kBlue + 1);
  g.SetLineColor(kBlue + 1);

  g.Draw("AP");
  g.GetXaxis()->SetTitle("BCM2_I (current)");
  g.GetYaxis()->SetTitle("HMS_hEL_CLEAN (counts/mC)");
  g.GetXaxis()->SetTitleSize(0.045);
  g.GetYaxis()->SetTitleSize(0.045);
  g.GetXaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetTitleOffset(1.6);
  g.GetXaxis()->SetLimits(xlo, xhi);
  g.SetMinimum(ylo);
  g.SetMaximum(yhi);

  // Constant fit (uses y-errors)
  TF1 f("f_const", "pol0", xlo, xhi);
  f.SetLineColor(kBlack);
  f.SetLineStyle(2);
  f.SetLineWidth(2);

  g.Fit(&f, "Q"); // quiet

  const double c0    = f.GetParameter(0);
  const double c0err = f.GetParError(0);

  const double chi2 = f.GetChisquare();
  const int ndf = f.GetNDF();
  const double prob = f.GetProb();
  const double chi2ndf = (ndf > 0 ? chi2 / (double)ndf : NAN);

  // Legend
  TLegend leg(0.12, 0.80, 0.55, 0.95);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(&g, "signal", "p");
  leg.AddEntry(&f, "const fit: y = c", "l");
  leg.Draw();

  // Setting info
  auto parts = SplitPath(rel);
  const std::string line1 = JoinSettingHeaderLine1(parts);
  const std::string line2 = JoinSettingHeaderLine2(parts);

  TLatex tl;
  tl.SetNDC(true);
  tl.SetTextFont(42);
  tl.SetTextSize(0.03);
  tl.SetTextAlign(13);
  tl.DrawLatex(0.35, 0.94, line1.c_str());
  if (!line2.empty()) tl.DrawLatex(0.35, 0.89, line2.c_str());

  // Fit stats text inside plot
  TLatex txt;
  txt.SetNDC(true);
  txt.SetTextFont(42);
  txt.SetTextColor(kBlue + 2);
  txt.SetTextSize(0.040);
  txt.SetTextAlign(22);
  txt.DrawLatex(0.52, 0.30, Form("#chi^{2}/ndf=%.2f (%d), Prob=%.3g", chi2ndf, ndf, prob));
  txt.DrawLatex(0.52, 0.24, Form("c=%.6g #pm %.3g", c0, c0err));

  // CSV fit summary
  csv << "# FIT_MODEL,pol0\n";
  csv << "# FIT_c," << std::setprecision(12) << c0 << "\n";
  csv << "# FIT_cerr," << std::setprecision(12) << c0err << "\n";
  csv << "# FIT_chi2," << std::setprecision(12) << chi2 << "\n";
  csv << "# FIT_ndf," << ndf << "\n";
  csv << "# FIT_chi2ndf," << std::setprecision(12) << chi2ndf << "\n";
  csv << "# FIT_prob," << std::setprecision(12) << prob << "\n";
  csv.close();

  c.SaveAs(outPng.c_str());
  TFile fout(outRoot.c_str(), "RECREATE");
  c.Write("c_hms_elclean_chargeNorm_vs_current");
  g.Write("g_hms_elclean_chargeNorm_vs_current");
  f.Write("f_const");
  fout.Close();

  log << "Fit: c=" << c0 << " +/- " << c0err << "\n";
  log << "chi2/ndf=" << chi2ndf << " (" << ndf << ")  prob=" << prob << "\n";
  log << "Wrote CSV: " << outCSV << "\n";
  std::cout << "Wrote CSV: " << outCSV << "\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT: " << outRoot << "\n";
  std::cout << "Wrote ROOT: " << outRoot << "\n";
  log.close();
}
