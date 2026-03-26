// PlotOffsetSummaryBySetting.C
//
// Purpose:
//   Read heep_check_offsetSummary_bySetting.csv and plot mean measured offsets
//   vs dp bin center, separately for Setting A and Setting B.
//
// Input:
//   results/tables/heep_check_offsetSummary_bySetting.csv
//
// Output PNGs:
//   results/PNGs/PlotOffsetSummaryBySetting/offset_summary_settingA.png
//   results/PNGs/PlotOffsetSummaryBySetting/offset_summary_settingB.png
//
// Example:
//   root -l -b -q 'macros/PlotOffsetSummaryBySetting.C+("results/tables/heep_check_offsetSummary_bySetting.csv")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStyle.h"

struct Row {
  TString setting = "";
  int dp_idx = -999;
  TString dp_label = "";
  double dp_lo = 0.0;
  double dp_hi = 0.0;
  TString var = "";

  int n_valid = 0;
  double mean_dmu_MeV = 0.0;
  double rms_dmu_MeV = 0.0;
  double sem_dmu_MeV = 0.0;
};

static std::vector<TString> SplitCSV(const std::string& line) {
  std::vector<TString> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(TString(item));
  }
  return out;
}

static int VarOrder(const TString& v) {
  if (v == "W")   return 0;
  if (v == "Em")  return 1;
  if (v == "Pmz") return 2;
  if (v == "Pmy") return 3;
  return 99;
}

static int VarColor(const TString& v) {
  if (v == "W")   return kBlack;
  if (v == "Em")  return kRed + 1;
  if (v == "Pmz") return kBlue + 1;
  if (v == "Pmy") return kGreen + 2;
  return kGray + 2;
}

static int VarMarker(const TString& v) {
  if (v == "W")   return 20;
  if (v == "Em")  return 21;
  if (v == "Pmz") return 22;
  if (v == "Pmy") return 23;
  return 24;
}

static TString PrettyVar(const TString& v) {
  if (v == "W")   return "W";
  if (v == "Em")  return "Em";
  if (v == "Pmz") return "Pmz";
  if (v == "Pmy") return "Pmy";
  return v;
}

static void MakeOnePlot(const std::vector<Row>& rows, const TString& setting, const TString& outPng) {
  std::vector<TString> vars = {"W", "Em", "Pmz", "Pmy"};

  // Determine x/y ranges
  bool hasAny = false;
  double xmin =  1e9, xmax = -1e9;
  double ymin =  1e9, ymax = -1e9;

  for (const auto& r : rows) {
    if (r.setting != setting) continue;
    if (r.n_valid <= 0) continue;

    double x = 0.5 * (r.dp_lo + r.dp_hi);
    double y = r.mean_dmu_MeV;
    double ey = r.sem_dmu_MeV;

    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y - ey);
    ymax = std::max(ymax, y + ey);
    hasAny = true;
  }

  if (!hasAny) {
    std::cerr << "[WARN] No valid rows found for Setting " << setting << "\n";
    return;
  }

  double xpad = 0.15 * std::max(1.0, xmax - xmin);
  double ypad = 0.20 * std::max(1.0, ymax - ymin);

  xmin -= xpad;
  xmax += xpad;
  ymin -= ypad;
  ymax += ypad;

  if (ymin > 0.0) ymin = std::min(0.0, ymin - 1.0);

  TCanvas* c = new TCanvas(Form("c_%s", setting.Data()),
                           Form("Offset summary Setting %s", setting.Data()),
                           1100, 800);
  c->SetMargin(0.12, 0.05, 0.12, 0.08);
  c->SetTicks(1, 1);

  TH1D* frame = new TH1D(Form("frame_%s", setting.Data()), "", 100, xmin, xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->SetTitle("SHMS dp bin center (%)");
  frame->GetYaxis()->SetTitle("Mean measured offset  #Delta#mu  (MeV)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->Draw();

  TLegend* leg = new TLegend(0.68, 0.68, 0.90, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);

  std::vector<TGraphErrors*> graphs;

  std::sort(vars.begin(), vars.end(),
            [](const TString& a, const TString& b) { return VarOrder(a) < VarOrder(b); });

  for (const auto& var : vars) {
    std::vector<double> xFocus, exFocus, yFocus, eyFocus;
    std::vector<double> xFull,  exFull,  yFull,  eyFull;

    for (const auto& r : rows) {
      if (r.setting != setting) continue;
      if (r.var != var) continue;
      if (r.n_valid <= 0) continue;

      double x = 0.5 * (r.dp_lo + r.dp_hi);

      if (r.dp_label == "full") {
        xFull.push_back(x);
        exFull.push_back(0.0);
        yFull.push_back(r.mean_dmu_MeV);
        eyFull.push_back(r.sem_dmu_MeV);
      } else {
        xFocus.push_back(x);
        exFocus.push_back(0.0);
        yFocus.push_back(r.mean_dmu_MeV);
        eyFocus.push_back(r.sem_dmu_MeV);
      }
    }

    // Draw focus-bin trend only (b1..b5)
    TGraphErrors* gFocus = nullptr;
    if (!xFocus.empty()) {
      gFocus = new TGraphErrors((int)xFocus.size());
      for (int i = 0; i < (int)xFocus.size(); ++i) {
        gFocus->SetPoint(i, xFocus[i], yFocus[i]);
        gFocus->SetPointError(i, exFocus[i], eyFocus[i]);
      }

      gFocus->SetLineColor(VarColor(var));
      gFocus->SetMarkerColor(VarColor(var));
      gFocus->SetMarkerStyle(VarMarker(var));
      gFocus->SetMarkerSize(1.3);
      gFocus->SetLineWidth(2);

      gFocus->Draw("PL SAME");
      leg->AddEntry(gFocus, PrettyVar(var), "pl");
      graphs.push_back(gFocus);
    }

    // Draw full-bin point separately (not connected)
    if (!xFull.empty()) {
      TGraphErrors* gFull = new TGraphErrors((int)xFull.size());
      for (int i = 0; i < (int)xFull.size(); ++i) {
        gFull->SetPoint(i, xFull[i], yFull[i]);
        gFull->SetPointError(i, exFull[i], eyFull[i]);
      }

      gFull->SetLineColor(VarColor(var));
      gFull->SetMarkerColor(VarColor(var));
      gFull->SetMarkerStyle(24);   // open marker for "full"
      gFull->SetMarkerSize(1.5);
      gFull->SetLineWidth(0);

      gFull->Draw("P SAME");
      graphs.push_back(gFull);
    }
  }

  leg->Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.14, 0.93, Form("HEEP check offset summary  -  Setting %s", setting.Data()));
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.14, 0.885, "Y error bars = SEM = RMS / #sqrt{N_{valid runs}}");

  c->SaveAs(outPng);

  delete frame;
  delete leg;
  for (auto* g : graphs) delete g;
  delete c;
}

void PlotOffsetSummaryBySetting(
    const char* inCsvC = "results/tables/heep_check_offsetSummary_bySetting.csv") {

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);

  TString inCsv(inCsvC);

  std::ifstream fin(inCsv.Data());
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Cannot open input CSV: " << inCsv << "\n";
    return;
  }

  TString outDir = "results/PNGs/PlotOffsetSummaryBySetting";
  if (gSystem->AccessPathName(outDir)) {
    gSystem->mkdir(outDir, kTRUE);
  }

  std::string line;
  if (!std::getline(fin, line)) {
    std::cerr << "[ERROR] Empty CSV: " << inCsv << "\n";
    return;
  }

  std::vector<TString> header = SplitCSV(line);
  std::map<std::string,int> col;
  for (int i = 0; i < (int)header.size(); ++i) {
    col[header[i].Data()] = i;
  }

  auto need = [&](const char* name) {
    if (!col.count(name)) {
      std::cerr << "[ERROR] Missing required column: " << name << "\n";
      return false;
    }
    return true;
  };

  if (!need("setting")       ||
      !need("dp_idx")        ||
      !need("dp_label")      ||
      !need("dp_lo")         ||
      !need("dp_hi")         ||
      !need("var")           ||
      !need("n_valid")       ||
      !need("mean_dmu_MeV")  ||
      !need("rms_dmu_MeV")   ||
      !need("sem_dmu_MeV")) {
    return;
  }

  std::vector<Row> rows;

  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    std::vector<TString> tok = SplitCSV(line);
    if ((int)tok.size() < (int)header.size()) continue;

    Row r;
    r.setting      = tok[col["setting"]];
    r.dp_idx       = tok[col["dp_idx"]].Atoi();
    r.dp_label     = tok[col["dp_label"]];
    r.dp_lo        = tok[col["dp_lo"]].Atof();
    r.dp_hi        = tok[col["dp_hi"]].Atof();
    r.var          = tok[col["var"]];
    r.n_valid      = tok[col["n_valid"]].Atoi();
    r.mean_dmu_MeV = tok[col["mean_dmu_MeV"]].Atof();
    r.rms_dmu_MeV  = tok[col["rms_dmu_MeV"]].Atof();
    r.sem_dmu_MeV  = tok[col["sem_dmu_MeV"]].Atof();

    rows.push_back(r);
  }
  fin.close();

  MakeOnePlot(rows, "A", outDir + "/offset_summary_settingA.png");
  MakeOnePlot(rows, "B", outDir + "/offset_summary_settingB.png");

  std::cout << "[INFO] Wrote PNGs to: " << outDir << "\n";
}
