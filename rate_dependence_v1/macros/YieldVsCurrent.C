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
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TString.h"
#include "TPad.h"
#include "TH1.h"
#include "TLatex.h"
#include "TF1.h"
#include "TMath.h"

#include "../include/CoincidenceRandomSubtraction.h"

// ---------- helpers ----------
static double SafeToD(const std::string& s) {
  try { return std::stod(s); } catch (...) { return NAN; }
}
static bool IsBadSentinel(double v) {
  if (!std::isfinite(v)) return true;
  return (std::fabs(v + 999.0) < 1e-9);
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
  const std::string needle = "/settings/";
  auto pos = settingDirAbs.find(needle);
  if (pos == std::string::npos) return BaseName(settingDirAbs);
  std::string rel = settingDirAbs.substr(pos + needle.size());
  while (!rel.empty() && rel.front()=='/') rel.erase(0,1);
  return rel.empty() ? BaseName(settingDirAbs) : rel;
}
static std::vector<std::string> SplitPath(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c=='/') {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else cur.push_back(c);
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}
static void BuildSettingLabel(const std::string& relUnderSettings,
                              std::string& line1, std::string& line2)
{
  // pass4/pi+sidis/LH2/z0p36/x0p25/Q23p3/<setting_id>
  auto tok = SplitPath(relUnderSettings);
  if (tok.size() >= 7) {
    line1 = tok[0] + " / " + tok[1] + " / " + tok[2] + " / " +
            tok[3] + " / " + tok[4] + " / " + tok[5];
    line2 = tok[6];
  } else {
    line1 = relUnderSettings;
    line2 = "";
  }
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
static std::string FormatWindows(const std::vector<std::pair<double,double>>& wins) {
  std::ostringstream ss;
  for (size_t i=0;i<wins.size();++i) {
    if (i) ss << ";";
    ss << wins[i].first << "-" << wins[i].second;
  }
  return ss.str();
}

// ---------- metadata row ----------
struct MetaRow {
  std::string category;
  int run = -1;

  double BCM2_I = NAN;
  double BCM2_Q = NAN;
  double comp_livetime = NAN;
  double h_esing_Eff   = NAN;
  double boil_corr     = NAN;

  double ps_factor     = NAN;
  std::string ps_choice;

  std::string target_raw;
  std::string run_type_raw;
};

static std::vector<MetaRow> ReadRunMetadata(const std::string& path, std::ostream& log) {
  std::vector<MetaRow> rows;
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

    MetaRow m;
    m.category = col(r, "category");
    try { m.run = std::stoi(col(r, "run")); } catch (...) { continue; }

    m.BCM2_I        = SafeToD(col(r, "BCM2_I"));
    m.BCM2_Q        = SafeToD(col(r, "BCM2_Q"));
    m.comp_livetime = SafeToD(col(r, "comp_livetime"));
    m.h_esing_Eff   = SafeToD(col(r, "h_esing_Eff"));
    m.boil_corr     = SafeToD(col(r, "boil_corr"));

    m.ps_factor     = SafeToD(col(r, "ps_factor"));
    m.ps_choice     = col(r, "ps_choice");

    m.target_raw    = col(r, "target_raw");
    m.run_type_raw  = col(r, "run_type_raw");

    rows.push_back(m);
  }

  log << "Loaded rows from run_metadata.csv: " << rows.size() << "\n";
  return rows;
}

// ---------- main ----------
void YieldVsCurrent(const char* manifestPath,
                    const char* resultsDir)
{
  const std::string manifestP = manifestPath ? manifestPath : "";
  const std::string outRoot   = resultsDir   ? resultsDir   : "";

  const std::string settingDir = DirName(manifestP);
  const std::string setting_id = BaseName(settingDir);

  const std::string rel = RelUnderSettings(settingDir);
  const std::string outBase = outRoot + "/" + rel;

  const std::string outPNGs = outBase + "/PNGs";
  const std::string outTabs = outBase + "/tables";
  const std::string outCanv = outBase + "/canvases";
  const std::string outLogs = outBase + "/logs";
  EnsureDir(outPNGs);
  EnsureDir(outTabs);
  EnsureDir(outCanv);
  EnsureDir(outLogs);

  const std::string logPath = outLogs + "/YieldVsCurrent.log";
  std::ofstream log(logPath.c_str());

  auto warn = [&](int run, const std::string& msg){
    std::cerr << "WARNING [run " << run << "]: " << msg << "\n";
    log       << "WARNING [run " << run << "]: " << msg << "\n";
  };

  log << "YieldVsCurrent (v1, signal-only)\n";
  log << "manifest: " << manifestP << "\n";
  log << "setting_dir: " << settingDir << "\n";
  log << "setting_id: " << setting_id << "\n";
  log << "results_base: " << outBase << "\n";
  log << "log: " << logPath << "\n\n";

  // Read metadata
  const std::string metaPath = settingDir + "/run_metadata.csv";
  auto meta = ReadRunMetadata(metaPath, log);
  if (meta.empty()) {
    log << "ERROR: no metadata rows. Nothing to do.\n";
    return;
  }

  // Physics cuts
  const TString baseCuts =
      "(H_gtr_dp>-8) && (H_gtr_dp<8) && "
      "(H_cal_etottracknorm>0.7) && "
      "(H_cer_npeSum>2.0) && "
      "(P_gtr_dp>-10) && (P_gtr_dp<22) && "
      "(P_cal_etottracknorm<0.8)";

  // Coincidence random subtraction config
  CoincidenceConfig cfg;
  cfg.CtBranchName = "CTime_ePiCoinTime_ROC2";

  // Output CSV
  const std::string outCsv = outTabs + "/yield_vs_current_signal.csv";
  std::ofstream csv(outCsv.c_str());
  csv << "category,run,BCM2_I,"
      << "Nnet,Nnet_err,"
      << "norm_factor,yield,yield_err,"
      << "BCM2_Q_mC,comp_livetime,h_esing_Eff,boil_corr,ps_choice,ps_factor,"
      << "PeakCenterNs,CoinLoNs,CoinHiNs,RandomWindowListNs,rootfile,status,"
      << "fit_included,fit_excluded_reason\n";

  // For plotting + fitting (signal only)
  std::vector<double> x, y, ex, ey;
  std::vector<int>    run_used;

  // Global range calc only from signal plot points (finite)
  bool any = false;
  double minX=0,maxX=0,minY=0,maxY=0;

  // ROOT file location
  const std::string skimDir = "./Skimmed_ROOTfiles";

  // Loop runs
  for (const auto& m : meta) {
    if (m.category != "signal") continue; // v1: signal only

    const int run = m.run;
    const std::string rootfile = Form("%s/skimmed_coin_replay_production_%d_-1.root", skimDir.c_str(), run);

    std::string status = "OK";

    // Validate normalization scalars
    bool badNorm = false;
    if (IsBadSentinel(m.BCM2_I)) { warn(run, "BCM2_I is invalid (-999/NaN)."); badNorm = true; }
    if (IsBadSentinel(m.BCM2_Q) || m.BCM2_Q <= 0) {
      warn(run, Form("BCM2_Q invalid (%.6g).", m.BCM2_Q)); badNorm = true;
    }
    if (IsBadSentinel(m.comp_livetime) || m.comp_livetime <= 0) {
      warn(run, Form("comp_livetime invalid (%.6g).", m.comp_livetime)); badNorm = true;
    }
    if (IsBadSentinel(m.h_esing_Eff) || m.h_esing_Eff <= 0) {
      warn(run, Form("h_esing_Eff invalid (%.6g).", m.h_esing_Eff)); badNorm = true;
    }
    if (IsBadSentinel(m.boil_corr) || m.boil_corr <= 0) {
      warn(run, Form("boil_corr invalid (%.6g).", m.boil_corr)); badNorm = true;
    }
    if (IsBadSentinel(m.ps_factor) || m.ps_factor <= 0) {
      warn(run, Form("ps_factor invalid (%.6g).", m.ps_factor)); badNorm = true;
    }

    // Open ROOT file
    std::unique_ptr<TFile> f(TFile::Open(rootfile.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      warn(run, "Cannot open ROOT file: " + rootfile);
      status = "NOFILE";
      csv << m.category << "," << run << "," << m.BCM2_I << ","
          << "nan,nan,"
          << "nan,nan,nan,"
          << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
          << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
          << "nan,nan,nan,"
          << "\"\","
          << "\"" << rootfile << "\","
          << status << ","
          << 0 << ",\"NOFILE\"\n";
      continue;
    }

    TTree* T = (TTree*)f->Get("T");
    if (!T) {
      warn(run, "Tree 'T' not found in file: " + rootfile);
      status = "NOTREE";
      csv << m.category << "," << run << "," << m.BCM2_I << ","
          << "nan,nan,"
          << "nan,nan,nan,"
          << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
          << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
          << "nan,nan,nan,"
          << "\"\","
          << "\"" << rootfile << "\","
          << status << ","
          << 0 << ",\"NOTREE\"\n";
      continue;
    }

    // Coincidence random subtraction
    CoincidenceResult R = ComputeCoincidenceRandomSubtraction(T, baseCuts, cfg);

    const double Nnet    = R.RandomSubtractedYield;
    const double NnetErr = R.RandomSubtractedYieldErr;

    double norm_factor = NAN;
    double yy = NAN, yerr = NAN;

    if (!badNorm) {
      norm_factor = (m.boil_corr * m.ps_factor) / (m.comp_livetime * m.h_esing_Eff * m.BCM2_Q);
      yy   = Nnet    * norm_factor;
      yerr = NnetErr * norm_factor;
    } else {
      status = "BADMETA";
    }

    // Fit inclusion logic: we skip zero/nan errors for fit, but still report
    int fit_included = 0;
    std::string fit_excl = "";
    if (!std::isfinite(m.BCM2_I) || !std::isfinite(yy)) {
      fit_excl = "NONFINITE_XY";
    } else if (!std::isfinite(yerr) || yerr <= 0) {
      fit_excl = "BAD_YERR";
      warn(run, Form("yield_err invalid for fit (%.6g). Will exclude from fit.", yerr));
    } else if (status != "OK") {
      fit_excl = "STATUS_NOT_OK";
    } else {
      fit_included = 1;
    }

    // Random window list string
    const std::string randList = FormatWindows(R.RandomWindowListNs);

    // Write CSV row
    csv << m.category << "," << run << "," << m.BCM2_I << ","
        << Nnet << "," << NnetErr << ","
        << norm_factor << "," << yy << "," << yerr << ","
        << m.BCM2_Q << "," << m.comp_livetime << "," << m.h_esing_Eff << "," << m.boil_corr << ","
        << "\"" << m.ps_choice << "\"," << m.ps_factor << ","
        << R.PeakCenterNs << "," << R.CoinWindowNs.first << "," << R.CoinWindowNs.second << ","
        << "\"" << randList << "\"" << ","
        << "\"" << rootfile << "\","
        << status << ","
        << fit_included << ","
        << "\"" << fit_excl << "\"\n";

    // Add to plot vectors if finite (plot can include points even if yerr bad; but ey needs finite)
    // We keep plotting points if yy finite and BCM2_I finite and yerr finite (>=0 allowed for plotting).
    if (std::isfinite(m.BCM2_I) && std::isfinite(yy) && std::isfinite(yerr)) {
      x.push_back(m.BCM2_I);
      y.push_back(yy);
      ex.push_back(0.0);
      ey.push_back(yerr);
      run_used.push_back(run);

      // Update global plot ranges with y±err only if err finite
      double ylo = yy - yerr;
      double yhi = yy + yerr;
      if (!any) { minX=maxX=m.BCM2_I; minY=ylo; maxY=yhi; any=true; }
      else {
        minX = std::min(minX, m.BCM2_I);
        maxX = std::max(maxX, m.BCM2_I);
        minY = std::min(minY, ylo);
        maxY = std::max(maxY, yhi);
      }
    }
  }

  csv.close();
  log << "\nWrote CSV: " << outCsv << "\n";

  const std::string outPng      = outPNGs + "/yield_vs_current.png";
  const std::string outRootFile = outCanv + "/yield_vs_current.root";

  if (!any || x.empty()) {
    warn(-1, "No finite points to plot. Saving empty canvas.");
    TCanvas c("c_yield_vs_current", "Yield vs Current (signal)", 1200, 850);
    c.SaveAs(outPng.c_str());
    c.SaveAs(outRootFile.c_str());
    log << "Wrote PNG (empty): " << outPng << "\n";
    log << "Wrote ROOT (empty): " << outRootFile << "\n";
    log.close();
    return;
  }

  // Sort by current
  SortByX(x, y, ex, ey);

  // Compute padded ranges
  double dx = maxX - minX; if (dx <= 0) dx = 1.0;
  double dy = maxY - minY; if (dy <= 0) dy = 1.0;
  double xpad = 0.06 * dx;
  double ypad = 0.10 * dy;

  double xmin = minX - xpad;
  double xmax = maxX + xpad;
  double ymin = minY - ypad;
  double ymax = maxY + ypad;

  // Build fit-only graph from points with ey>0 and finite
  std::vector<double> fx, fy, fex, fey;
  fx.reserve(x.size()); fy.reserve(x.size()); fex.reserve(x.size()); fey.reserve(x.size());

  for (size_t i=0;i<x.size();++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(ey[i]) || ey[i] <= 0) continue;
    fx.push_back(x[i]);
    fy.push_back(y[i]);
    fex.push_back(0.0);
    fey.push_back(ey[i]);
  }

  double fitC = NAN, fitCerr = NAN;
  double chi2 = NAN;
  int ndf = -1;
  double prob = NAN;

  TGraphErrors* gFit = nullptr;
  TF1* fconst = nullptr;

  if (fx.size() >= 2) {
    gFit = new TGraphErrors((int)fx.size(), fx.data(), fy.data(), fex.data(), fey.data());

    // Constant fit over x-range (doesn't matter for a constant, but ROOT wants it)
    fconst = new TF1("fconst", "[0]", xmin, xmax);

    // Seed with simple average (optional)
    double avg = 0;
    for (double v : fy) avg += v;
    avg /= std::max<size_t>(1, fy.size());
    fconst->SetParameter(0, avg);

    // Fit (quiet) using errors
    gFit->Fit(fconst, "Q");

    fitC    = fconst->GetParameter(0);
    fitCerr = fconst->GetParError(0);
    chi2 = fconst->GetChisquare();
    ndf  = fconst->GetNDF();
    prob = TMath::Prob(chi2, ndf);

    log << "Fit constant C = " << fitC << " +- " << fitCerr << "\n";
    log << "chi2=" << chi2 << " ndf=" << ndf
        << " chi2/ndf=" << (ndf>0 ? chi2/ndf : NAN)
        << " prob=" << prob << "\n";
  } else {
    warn(-1, "Not enough fit-eligible points (need >=2 with yield_err>0). Will plot without fit.");
    log << "NOTE: not enough fit-eligible points for constant fit.\n";
  }

  // Plot
  TCanvas c("c_yield_vs_current", "Yield vs Current (signal)", 1200, 850);
  gStyle->SetOptStat(0);
  c.SetTopMargin(0.22);
  c.SetRightMargin(0.06);
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);

  TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  frame->GetXaxis()->SetTitle("BCM2_I (current)");
  frame->GetYaxis()->SetTitle("Normalized yield (myCTime)");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);

  // Main graph
  TGraphErrors g((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g.SetMarkerStyle(20);
  g.SetMarkerSize(1.55);
  g.SetMarkerColor(kRed+1);
  g.SetLineColor(kRed+1);
  g.Draw("P SAME");

  // Fit line (dashed)
  if (fconst && std::isfinite(fitC)) {
    fconst->SetLineColor(kBlack);
    fconst->SetLineStyle(2);  // dashed
    fconst->SetLineWidth(2);
    fconst->Draw("SAME");
  }

  // Legend
  TLegend leg(0.10, 0.82, 0.50, 0.98);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.035);
  leg.AddEntry(&g, "signal", "p");
  if (fconst && std::isfinite(fitC)) {
    leg.AddEntry(fconst, Form("const fit: C=%.3g", fitC), "l");
  }
  leg.Draw();

  // Setting label (TLatex)
  std::string lbl1, lbl2;
  BuildSettingLabel(rel, lbl1, lbl2);

  TLatex t;
  t.SetNDC();
  t.SetTextAlign(13);
  t.SetTextSize(0.028);
  t.DrawLatex(0.39, 0.94, lbl1.c_str());
  t.DrawLatex(0.39, 0.89, lbl2.c_str());

  // Fit summary text (optional small line)
  if (fconst && std::isfinite(chi2) && ndf > 0) {
    TLatex tf;
    tf.SetNDC();
    tf.SetTextFont(42);
    tf.SetTextColor(kBlue+2);
    tf.SetTextAlign(22);
    tf.SetTextSize(0.040);
    tf.DrawLatex(0.52, 0.30,
      Form("#chi^{2}/ndf=%.2f (%g), Prob=%.3g", chi2/ndf, (double)ndf, prob));
  }

  c.SaveAs(outPng.c_str());
  c.SaveAs(outRootFile.c_str());

  log << "Global ranges: X=[" << xmin << "," << xmax << "]  Y=[" << ymin << "," << ymax << "]\n";
  log << "Wrote PNG: " << outPng << "\n";
  log << "Wrote ROOT: " << outRootFile << "\n";

  // Append fit summary to CSV as comment footer
  {
    std::ofstream csv2(outCsv.c_str(), std::ios::app);
    csv2 << "#FIT,fit_const,fit_const_err,chi2,ndf,chi2_ndf,prob,n_fit_points\n";
    if (fconst && std::isfinite(fitC) && ndf >= 0) {
      csv2 << "#FIT," << fitC << "," << fitCerr << "," << chi2 << "," << ndf << ","
           << (ndf>0 ? chi2/ndf : NAN) << "," << prob << "," << fx.size() << "\n";
    } else {
      csv2 << "#FIT,nan,nan,nan,-1,nan,nan," << fx.size() << "\n";
    }
  }

  if (gFit) delete gFit;
  if (fconst) delete fconst;

  log.close();
}
