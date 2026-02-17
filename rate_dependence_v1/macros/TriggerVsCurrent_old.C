// TriggerVsCurrent.C
// v1: Signal-only SHMS_EL_CLEAN rate (kHz) vs BCM2_I (current), linear fit (pol1)
// Includes slope significance = slope / slope_err
// Inputs:
//   TriggerVsCurrent("<abs path>/settings/.../manifest.json", "<abs path>/results")
//
// Reads:
//   <results>/<same rel path as manifest>/tables/yield_vs_current_signal.csv
// Report files:
//   <repo_root>/Pass0_REPORTfiles/COIN/PRODUCTION/replay_coin_production_<RUN>_-1.report
//
// Writes:
//   <results>/<rel>/tables/shms_elclean_vs_current.csv
//   <results>/<rel>/PNGs/shms_elclean_vs_current.png
//   <results>/<rel>/canvases/shms_elclean_vs_current.root
//   <results>/<rel>/logs/TriggerVsCurrent.log

#include <TCanvas.h>
#include <TGraph.h>
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
  if (kpos != std::string::npos) {
    rel = manifestPath.substr(kpos + key.size());
  } else {
    rel = manifestPath;
  }
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

static bool ParseReport_SHMS_ELCLEAN_kHz(const std::string& reportPath, double& out) {
  std::ifstream in(reportPath);
  if (!in) return false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.find("SHMS_pEL_CLEAN") != std::string::npos && line.find('[') != std::string::npos &&
        line.find("kHz") != std::string::npos) {
      auto lb = line.find('[');
      auto rb = line.find(']', lb == std::string::npos ? 0 : lb);
      if (lb == std::string::npos || rb == std::string::npos || rb <= lb) continue;
      std::string inside = line.substr(lb + 1, rb - lb - 1);
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
  double yield = NAN;        // not used for fit, but saved in CSV for context
  double yerr = NAN;         // not used
  double shms_elclean = NAN; // kHz
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
  const int iI = need("BCM2_I");
  const int iY = need("yield");
  const int iYe = need("yield_err");
  const int iStatus = need("status");
  const int iRoot = need("rootfile");

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
    try {
      r.run = std::stoi(v[iRun]);
    } catch (...) {
      log << "WARN: bad run value, skipping line: " << line << "\n";
      continue;
    }
    try {
      r.current = std::stod(v[iI]);
    } catch (...) {
      r.status = "BAD_CURRENT";
    }
    if (iY >= 0) {
      try { r.yield = std::stod(v[iY]); } catch (...) {}
    }
    if (iYe >= 0) {
      try { r.yerr = std::stod(v[iYe]); } catch (...) {}
    }
    if (iStatus >= 0) r.status = v[iStatus];
    if (iRoot >= 0) r.rootfile = v[iRoot];

    rows.push_back(r);
    n++;
  }

  log << "Read " << n << " rows from: " << csvPath << "\n";
  return true;
}

static std::string GuessRepoRootFromManifest(const std::string& manifestPath) {
  const std::string key = "/rate_dependence_v1/";
  auto pos = manifestPath.find(key);
  if (pos != std::string::npos) {
    return manifestPath.substr(0, pos + key.size() - 1);
  }
  return gSystem->WorkingDirectory();
}

} // namespace

void TriggerVsCurrent(const char* manifest_json, const char* results_top) {
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

  const std::string repoRoot = GuessRepoRootFromManifest(manifestPath);
  const std::string reportDir = repoRoot + "/Pass0_REPORTfiles/COIN/PRODUCTION";

  int nOK = 0, nSkip = 0, nNoReport = 0, nNoLine = 0;
  for (auto& r : rows) {
    if (!IsFinite(r.current)) {
      r.status = "NONFINITE_CURRENT";
      nSkip++;
      continue;
    }

    const std::string rep = reportDir + "/replay_coin_production_" + std::to_string(r.run) + "_-1.report";
    double elc = NAN;
    if (!ParseReport_SHMS_ELCLEAN_kHz(rep, elc)) {
      std::ifstream test(rep);
      if (!test) {
        r.status = "MISSING_REPORT";
        nNoReport++;
      } else {
        r.status = "MISSING_ELCLEAN_LINE";
        nNoLine++;
      }
      nSkip++;
      continue;
    }
    r.shms_elclean = elc;
    r.status = "OK";
    nOK++;
  }

  log << "Report dir: " << reportDir << "\n";
  log << "Points: OK=" << nOK << "  skipped=" << nSkip << "  missingReport=" << nNoReport
      << "  missingLine=" << nNoLine << "\n";

  // Write table (include all rows)
  const std::string outCSV = outTables + "/shms_elclean_vs_current.csv";
  std::ofstream csv(outCSV.c_str());
  csv << "category,run,BCM2_I,shms_elclean_kHz,yield,yield_err,rootfile,status\n";
  for (const auto& r : rows) {
    csv << "signal,"
        << r.run << ","
        << std::setprecision(10) << r.current << ","
        << std::setprecision(10) << (IsFinite(r.shms_elclean) ? r.shms_elclean : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.yield) ? r.yield : NAN) << ","
        << std::setprecision(10) << (IsFinite(r.yerr) ? r.yerr : NAN) << ","
        << "\"" << r.rootfile << "\"" << ","
        << r.status << "\n";
  }

  // Build graph from OK points
  std::vector<double> xs, ys;
  xs.reserve(rows.size());
  ys.reserve(rows.size());

  double xmin = +1e99, xmax = -1e99, ymin = +1e99, ymax = -1e99;
  for (const auto& r : rows) {
    if (r.status != "OK") continue;
    if (!IsFinite(r.current) || !IsFinite(r.shms_elclean)) continue;

    xs.push_back(r.current);
    ys.push_back(r.shms_elclean);

    xmin = std::min(xmin, r.current);
    xmax = std::max(xmax, r.current);
    ymin = std::min(ymin, r.shms_elclean);
    ymax = std::max(ymax, r.shms_elclean);
  }

  const std::string outPng = outPNGs + "/shms_elclean_vs_current.png";
  const std::string outRoot = outCanv + "/shms_elclean_vs_current.root";

  TCanvas c("c_shms_elclean_vs_current", "SHMS EL_CLEAN vs Current", 1200, 850);
  c.SetTopMargin(0.22);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.06);
  c.SetBottomMargin(0.12);

  if (xs.empty()) {
    log << "WARN: No points to plot.\n";
    c.SaveAs(outPng.c_str());
    TFile fout(outRoot.c_str(), "RECREATE");
    c.Write("c_shms_elclean_vs_current");
    fout.Close();
    csv << "# FIT,NO_POINTS\n";
    csv.close();
    log << "Wrote PNG (empty): " << outPng << "\n";
    log << "Wrote ROOT: " << outRoot << "\n";
    log.close();
    return;
  }

  const double xpad = 0.10 * (xmax - xmin);
  const double ypad = 0.20 * (ymax - ymin);
  const double xlo = xmin - xpad;
  const double xhi = xmax + xpad;
  const double ylo = ymin - ypad;
  const double yhi = ymax + ypad;

  TGraph g((int)xs.size(), xs.data(), ys.data());
  g.SetName("g_shms_elclean_vs_current");
  g.SetMarkerStyle(21);
  g.SetMarkerSize(1.35);
  g.SetMarkerColor(kBlue + 1);
  g.SetLineColor(kBlue + 1);

  g.Draw("AP");
  g.GetXaxis()->SetTitle("BCM2_I (current)");
  g.GetYaxis()->SetTitle("SHMS_EL_CLEAN rate (kHz)");
  g.GetXaxis()->SetTitleSize(0.045);
  g.GetYaxis()->SetTitleSize(0.045);
  g.GetXaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetLabelSize(0.04);
  g.GetYaxis()->SetTitleOffset(1.25);
  g.GetXaxis()->SetLimits(xlo, xhi);
  g.SetMinimum(ylo);
  g.SetMaximum(yhi);

  // Linear fit
  TF1 f("f_lin", "pol1", xlo, xhi);
  f.SetLineColor(kBlack);
  f.SetLineStyle(2);
  f.SetLineWidth(2);

  g.Fit(&f, "Q"); // quiet

  const double a = f.GetParameter(0);
  const double aerr = f.GetParError(0);
  const double b = f.GetParameter(1);
  const double berr = f.GetParError(1);
  const double sig = (berr > 0 ? b / berr : NAN);

  const double chi2 = f.GetChisquare();
  const int ndf = f.GetNDF();
  const double prob = f.GetProb();
  const double chi2ndf = (ndf > 0 ? chi2 / (double)ndf : NAN);

  // Legend
  TLegend leg(0.12, 0.80, 0.50, 0.95);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(&g, "signal", "p");
  leg.AddEntry(&f, Form("lin fit: y = a + b x"), "l");
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

  // Fit stats text inside plot (colored)
  TLatex txt;
  txt.SetNDC(true);
  txt.SetTextFont(42);
  txt.SetTextColor(kBlue + 2);
  txt.SetTextSize(0.040);
  txt.SetTextAlign(22);
  txt.DrawLatex(0.52, 0.30, Form("#chi^{2}/ndf=%.2f (%d), Prob=%.3g", chi2ndf, ndf, prob));
  txt.DrawLatex(0.52, 0.24, Form("b=%.3g #pm %.3g  (b/#sigma_{b}=%.2f)", b, berr, sig));

  // CSV fit summary
  csv << "# FIT_MODEL,pol1\n";
  csv << "# FIT_a," << std::setprecision(12) << a << "\n";
  csv << "# FIT_aerr," << std::setprecision(12) << aerr << "\n";
  csv << "# FIT_b," << std::setprecision(12) << b << "\n";
  csv << "# FIT_berr," << std::setprecision(12) << berr << "\n";
  csv << "# FIT_b_over_berr," << std::setprecision(12) << sig << "\n";
  csv << "# FIT_chi2," << std::setprecision(12) << chi2 << "\n";
  csv << "# FIT_ndf," << ndf << "\n";
  csv << "# FIT_chi2ndf," << std::setprecision(12) << chi2ndf << "\n";
  csv << "# FIT_prob," << std::setprecision(12) << prob << "\n";
  csv.close();

  c.SaveAs(outPng.c_str());
  TFile fout(outRoot.c_str(), "RECREATE");
  c.Write("c_shms_elclean_vs_current");
  g.Write("g_shms_elclean_vs_current");
  f.Write("f_lin");
  fout.Close();

  log << "Fit: a=" << a << " +/- " << aerr << "  b=" << b << " +/- " << berr << "  b/berr=" << sig << "\n";
  log << "chi2/ndf=" << chi2ndf << " (" << ndf << ")  prob=" << prob << "\n";
  log << "Wrote CSV: " << outCSV << "\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT: " << outRoot << "\n";
  log.close();
}
