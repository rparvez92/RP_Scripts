// PlotFitParameters.C
//
// Summary plots of fit parameters (M0,A,B) and key diagnostics from:
//   results/<setting_id or group_id>/tables/fit_parameters_<mode>_tierOn.csv
//
// Output:
//   results/<setting_id or group_id>/PNGs/FitParamPlots_tierOn/ (or _tierOff)/*.png
//
// Usage (from rsidis_xs_v5/):
//   root -l -q 'macros/PlotFitParameters.C("results/.../tables/fit_parameters_group_tierOn.csv")'
//
// This version:
//   - Uses only physics-usable fits: valid_fit_default==1
//   - Places legend in the top margin (right), and title in the top margin (left)
//   - Increases marker size
//   - Adds M0/A/B vs pt (in z bins) and chi2/ndf vs pt
//   - Adds small "missing bin" note when an expected z or pt bin has no values
//
// Run command:
//   root -l -b -q 'macros/PlotFitParameters.C(<csv link>)'
//   root -l -b -q 'macros/PlotFitParameters.C("results/grp_pass4_piplus_LH2_zOv_x0p25_q23p3_thpq2/tables/fit_parameters_group_tierOff.csv")'
//   root -l -b -q 'macros/PlotFitParameters.C("results/grp_pass4_piplus_LH2_zOv_x0p25_q23p3_thpq2/tables/fit_parameters_group_tierOn.csv")'

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TMath.h>
#include <TSystem.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

static constexpr double kMissing = -999.0;
static constexpr double kEps = 1e-9;

// Expected binning (used only for "missing bin" notes).
// If your binning changes, update these.
static const std::vector<double> kExpectedZCenters = {0.36, 0.50, 0.67, 0.90};
static const std::vector<std::pair<double,double>> kExpectedPtBins = {
  {0.0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}
};

static std::string Trim(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e-1]))) --e;
  return s.substr(b, e - b);
}

static std::vector<std::string> SplitCSV(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inQuotes = false;
  for (size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
    if (c == '"') {
      inQuotes = !inQuotes;
      continue;
    }
    if (c == ',' && !inQuotes) {
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

static bool NearlyEqual(double a, double b, double eps=1e-3) {
  return std::fabs(a-b) < eps;
}


static int InferBinIndex(double center, const std::vector<double>& centers, double tol=0.03) {
  if (!IsFinite(center) || center <= kMissing + 1) return -1;
  int best = -1;
  double bestd = 1e9;
  for (size_t i = 0; i < centers.size(); ++i) {
    double d = std::fabs(center - centers[i]);
    if (d < bestd) { bestd = d; best = static_cast<int>(i); }
  }
  if (best >= 0 && bestd < tol) return best;
  return -1;
}

// Use explicit bin indices when they exist in the CSV, otherwise infer from centers.
// This ensures graph keys and legend exemplar keys are consistent.
static std::string PtLabel(double lo, double hi) {
  std::ostringstream os;
  os << std::fixed << std::setprecision(2) << lo << ", " << hi;
  return os.str();
}

static std::string ZLabel(double z) {
  std::ostringstream os;
  os << std::fixed << std::setprecision(2) << z;
  return os.str();
}

struct Row {
  std::string mode;
  std::string group_id;
  std::string curve_label;
  int setting_id = -1;

  int pt_bin = -1;
  int z_bin  = -1;
  double pt_lo = kMissing, pt_hi = kMissing, pt_center = kMissing;
  double z_lo  = kMissing, z_hi  = kMissing, z_center  = kMissing;

  int npts = -1;
  int npts_fit = -1;

  double M0 = kMissing, M0_err = kMissing;
  double A  = kMissing, A_err  = kMissing;
  double B  = kMissing, B_err  = kMissing;

  double A_sig = kMissing;
  double B_sig = kMissing;

  double chi2 = kMissing;
  int ndf = -1;
  double chi2_ndf = kMissing;
  double prob = kMissing;

  int fit_tier = -1;
  int fit_status = -999;
  int cov_status = -999;
  int fit_flag_bits = 0;

  int valid_fit_default = 0;
  int valid_M0_default = 0;
  int valid_A_default  = 0;
  int valid_B_default  = 0;

  int attempts = -1;
  std::string note;
};

// Convert a row into stable integer keys for grouping. Prefer explicit bin indices if present;
// otherwise infer from the expected bin centers.
static int ZKeyFromRow(const Row& r) {
  if (r.z_bin >= 0) return r.z_bin;
  return InferBinIndex(r.z_center, kExpectedZCenters, 0.03);
}

static int PtKeyFromRow(const Row& r) {
  if (r.pt_bin >= 0) return r.pt_bin;

  // Build expected pT centers from the hard-coded expected pT bin edges.
  static const std::vector<double> ptCenters = []{
    std::vector<double> c;
    c.reserve(kExpectedPtBins.size());
    for (const auto& b : kExpectedPtBins) c.push_back(0.5*(b.first + b.second));
    return c;
  }();

  return InferBinIndex(r.pt_center, ptCenters, 0.05);
}

// Fallback labels when we only have integer keys (no exemplar row).
static std::string ZLabelFromIndex(int idx) {
  if (idx >= 0 && idx < (int)kExpectedZCenters.size())
    return Form("z = %.2f", kExpectedZCenters[idx]);
  return Form("zbin %d", idx);
}

static std::string PtBinLabelFromIndex(int idx) {
  if (idx >= 0 && idx < (int)kExpectedPtBins.size()) {
    const auto& b = kExpectedPtBins[idx];
    return Form("p_{T} [%.2f, %.2f]", b.first, b.second);
  }
  return Form("ptbin %d", idx);
}


static bool SafeAtoi(const std::string& s, int& out) {
  try {
    size_t idx = 0;
    int v = std::stoi(s, &idx);
    if (idx != s.size()) return false;
    out = v;
    return true;
  } catch (...) {
    return false;
  }
}

static bool SafeAtof(const std::string& s, double& out) {
  try {
    size_t idx = 0;
    double v = std::stod(s, &idx);
    if (idx != s.size()) return false;
    out = v;
    return true;
  } catch (...) {
    return false;
  }
}

static std::vector<Row> ReadRows(const std::string& csvPath) {
  std::ifstream in(csvPath);
  if (!in) {
    std::cerr << "ERROR: cannot open CSV: " << csvPath << "\n";
    return {};
  }

  std::string header;
  if (!std::getline(in, header)) {
    std::cerr << "ERROR: empty CSV: " << csvPath << "\n";
    return {};
  }
  auto cols = SplitCSV(header);
  std::map<std::string,int> col;
  for (size_t i = 0; i < cols.size(); ++i) col[cols[i]] = (int)i;

  auto idx = [&](const std::string& name)->int {
    auto it = col.find(name);
    return (it==col.end()) ? -1 : it->second;
  };

  auto get = [&](const std::vector<std::string>& v, const std::string& name)->std::string {
    int i = idx(name);
    if (i < 0 || i >= (int)v.size()) return "";
    return v[i];
  };

  std::vector<Row> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto v = SplitCSV(line);
    Row r;
    r.mode = get(v,"mode");
    r.group_id = get(v,"group_id");
    r.curve_label = get(v,"curve_label");

    SafeAtoi(get(v,"setting_id"), r.setting_id);
    SafeAtoi(get(v,"pt_bin"), r.pt_bin);
    SafeAtoi(get(v,"z_bin"),  r.z_bin);

    SafeAtof(get(v,"pt_lo"), r.pt_lo);
    SafeAtof(get(v,"pt_hi"), r.pt_hi);
    SafeAtof(get(v,"pt_center"), r.pt_center);

    SafeAtof(get(v,"z_lo"), r.z_lo);
    SafeAtof(get(v,"z_hi"), r.z_hi);
    SafeAtof(get(v,"z_center"), r.z_center);

    SafeAtoi(get(v,"npts"), r.npts);
    SafeAtoi(get(v,"npts_fit"), r.npts_fit);

    SafeAtof(get(v,"M0"), r.M0);
    SafeAtof(get(v,"M0_err"), r.M0_err);
    SafeAtof(get(v,"A"), r.A);
    SafeAtof(get(v,"A_err"), r.A_err);
    SafeAtof(get(v,"B"), r.B);
    SafeAtof(get(v,"B_err"), r.B_err);

    SafeAtof(get(v,"A_sig"), r.A_sig);
    SafeAtof(get(v,"B_sig"), r.B_sig);

    SafeAtof(get(v,"chi2"), r.chi2);
    SafeAtoi(get(v,"ndf"), r.ndf);
    SafeAtof(get(v,"chi2_ndf"), r.chi2_ndf);
    SafeAtof(get(v,"prob"), r.prob);

    SafeAtoi(get(v,"fit_tier"), r.fit_tier);
    SafeAtoi(get(v,"fit_status"), r.fit_status);
    SafeAtoi(get(v,"cov_status"), r.cov_status);
    SafeAtoi(get(v,"fit_flag_bits"), r.fit_flag_bits);

    SafeAtoi(get(v,"valid_fit_default"), r.valid_fit_default);
    SafeAtoi(get(v,"valid_M0_default"), r.valid_M0_default);
    SafeAtoi(get(v,"valid_A_default"), r.valid_A_default);
    SafeAtoi(get(v,"valid_B_default"), r.valid_B_default);

    SafeAtoi(get(v,"attempts"), r.attempts);
    r.note = get(v,"note");

    rows.push_back(std::move(r));
  }
  return rows;
}

static std::string InferOutDir(const std::string& csvPath) {
  // Decide tier tag from the input CSV path (keeps runtime simple)
  std::string lower = csvPath;
  for (auto& ch : lower) ch = std::tolower(static_cast<unsigned char>(ch));

  std::string subdir = "FitParamPlots";
  if (lower.find("tieron") != std::string::npos)  subdir = "FitParamPlots_tierOn";
  if (lower.find("tieroff") != std::string::npos) subdir = "FitParamPlots_tierOff";

  std::string out = csvPath;
  // Replace /tables/<file>.csv -> /PNGs/<subdir>
  auto pos = out.rfind("/tables/");
  if (pos != std::string::npos) {
    out = out.substr(0, pos) + "/PNGs/" + subdir;
  } else {
    // fallback: sibling PNGs
    out = std::string("PNGs/") + subdir;
  }
  gSystem->mkdir(out.c_str(), kTRUE);
  return out;
}

// --- plotting config ---
struct PlotConfig {
  bool require_valid_fit_default = true; // physics plots
  double marker_size = 1.8;

  // Margins (larger top margin to host legend + title)
  double top_margin = 0.22;
  double left_margin = 0.14;
  double right_margin = 0.05;
  double bottom_margin = 0.12;

  double title_text_size = 0.040;
  double note_text_size  = 0.028;
  double legend_text_size = 0.032;

  // Symmetric around 0 for A/B plots
  bool symmetric_AB = true;
  double symmetric_AB_min_abs = 0.25; // ensures negative shown even if all positive
};

static bool PassRowPhysics(const Row& r, const PlotConfig& cfg) {
  if (cfg.require_valid_fit_default && r.valid_fit_default != 1) return false;
  return true;
}

static void ApplyCanvasStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

static void DrawHeader(const PlotConfig& cfg, const std::string& title) {
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(cfg.title_text_size);
  lat.SetTextAlign(13); // left-top
  lat.DrawLatex(0.12, 1.0 - 0.02, title.c_str());
}

static void DrawMissingNote(const PlotConfig& cfg, const std::string& msg) {
  if (msg.empty()) return;
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(cfg.note_text_size);
  lat.SetTextAlign(11); // left-bottom
  // Put it inside the top margin, below the title.
  lat.DrawLatex(0.12, 1.0 - cfg.top_margin + 0.02, msg.c_str());
}

static TLegend* MakeTopLegend(const PlotConfig& cfg, double x1=0.55, double y1=0.80, double x2=0.95, double y2=0.98) {
  auto* leg = new TLegend(x1,y1,x2,y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(cfg.legend_text_size);
  return leg;
}

static std::pair<double,double> ComputeYRangeSymmetric0(const std::vector<double>& y, const std::vector<double>& ey, double minAbs) {
  double ymin = +1e99, ymax = -1e99;
  for (size_t i=0;i<y.size();++i) {
    double yi = y[i];
    double e  = (i<ey.size()) ? ey[i] : 0.0;
    if (!IsFinite(yi) || !IsFinite(e)) continue;
    ymin = std::min(ymin, yi - e);
    ymax = std::max(ymax, yi + e);
  }
  if (!(ymin < ymax)) {
    // fallback
    return {-minAbs, +minAbs};
  }
  double absMax = std::max(std::fabs(ymin), std::fabs(ymax));
  absMax = std::max(absMax, minAbs);
  absMax *= 1.15; // pad
  return {-absMax, +absMax};
}

static std::pair<double,double> ComputeYRangeAuto(const std::vector<double>& y, const std::vector<double>& ey, double padFrac=0.20) {
  double ymin = +1e99, ymax = -1e99;
  for (size_t i=0;i<y.size();++i) {
    double yi = y[i];
    double e  = (i<ey.size()) ? ey[i] : 0.0;
    if (!IsFinite(yi) || !IsFinite(e)) continue;
    ymin = std::min(ymin, yi - e);
    ymax = std::max(ymax, yi + e);
  }
  if (!(ymin < ymax)) return {0.0, 1.0};
  double span = ymax - ymin;
  if (span < 1e-12) span = std::max(std::fabs(ymax), 1.0);
  double pad = padFrac * span;
  return {ymin - pad, ymax + pad};
}

static std::string MissingZNote(const std::set<double>& present, const std::string& prefix="No data for z=") {
  std::vector<std::string> miss;
  for (double z : kExpectedZCenters) {
    bool ok=false;
    for (double p : present) if (NearlyEqual(p,z,5e-3)) { ok=true; break; }
    if (!ok) miss.push_back(ZLabel(z));
  }
  if (miss.empty()) return "";
  std::ostringstream os;
  os << prefix << "{";
  for (size_t i=0;i<miss.size();++i) {
    if (i) os << ",";
    os << miss[i];
  }
  os << "}";
  return os.str();
}

static std::string MissingPtNote(const std::set<std::pair<double,double>>& present, const std::string& prefix="No data for p_{T} bins=") {
  std::vector<std::string> miss;
  for (auto [lo,hi] : kExpectedPtBins) {
    bool ok=false;
    for (auto p : present) {
      if (NearlyEqual(p.first, lo, 5e-4) && NearlyEqual(p.second, hi, 5e-4)) { ok=true; break; }
    }
    if (!ok) {
      std::ostringstream lab;
      lab << "[" << std::fixed << std::setprecision(2) << lo << "," << hi << "]";
      miss.push_back(lab.str());
    }
  }
  if (miss.empty()) return "";
  std::ostringstream os;
  os << prefix << "{";
  for (size_t i=0;i<miss.size();++i) {
    if (i) os << ",";
    os << miss[i];
  }
  os << "}";
  return os.str();
}

// --- graph containers ---
struct GraphSet {
  // key -> graph
  // For vs_z: key is pt_bin
  // For vs_pt: key is z_bin
  std::map<int, TGraphErrors*> g;
};

static void DeleteGraphs(GraphSet& gs) {
  for (auto& kv : gs.g) delete kv.second;
  gs.g.clear();
}

// style palettes
static const int kColors[]  = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7};
static const int kMarkers[] = {20, 21, 22, 23, 33, 29};

static std::string PtLegendLabelFromAnyRow(const Row& r) {
  std::ostringstream os;
  os << "p_{T} [" << std::fixed << std::setprecision(2) << r.pt_lo << ", " << r.pt_hi << "]";
  return os.str();
}

static std::string ZLegendLabelFromAnyRow(const Row& r) {
  std::ostringstream os;
  os << "z = " << std::fixed << std::setprecision(2) << r.z_center;
  return os.str();
}

// Fill graphs for y vs z (one curve per pt bin)
static GraphSet BuildGraphsVsZ(const std::vector<Row>& rows, const PlotConfig& cfg,
                              const std::string& yField, const std::string& eyField,
                              std::set<double>* presentZ=nullptr,
                              std::set<std::pair<double,double>>* presentPt=nullptr) {
  GraphSet gs;
  for (const auto& r : rows) {
    if (!PassRowPhysics(r,cfg)) continue;

    double x = r.z_center;
    if (!IsFinite(x) || x==kMissing) continue;

    double y = kMissing, ey = 0.0;
    if (yField=="M0") y = r.M0;
    else if (yField=="A") y = r.A;
    else if (yField=="B") y = r.B;
    else if (yField=="M0_err") y = r.M0_err;
    else if (yField=="A_err")  y = r.A_err;
    else if (yField=="B_err")  y = r.B_err;
    else if (yField=="A_sig")  y = r.A_sig;
    else if (yField=="B_sig")  y = r.B_sig;
    else if (yField=="chi2_ndf") y = r.chi2_ndf;
    else if (yField=="prob") y = r.prob;
    else if (yField=="fit_tier") y = (double)r.fit_tier;
    else if (yField=="npts_fit") y = (double)r.npts_fit;

    if (eyField=="M0_err") ey = r.M0_err;
    else if (eyField=="A_err")  ey = r.A_err;
    else if (eyField=="B_err")  ey = r.B_err;
    else ey = 0.0;

    if (!IsFinite(y) || y==kMissing) continue;
    if (!IsFinite(ey) || ey==kMissing) ey = 0.0;

    int key = r.pt_bin;
    auto it = gs.g.find(key);
    if (it == gs.g.end()) {
      auto* gr = new TGraphErrors();
      gr->SetName(Form("g_%s_vs_z_ptbin%d", yField.c_str(), key));
      gs.g[key] = gr;
      it = gs.g.find(key);
    }

    auto* gr = it->second;
    int n = gr->GetN();
    gr->SetPoint(n, x, y);
    gr->SetPointError(n, 0.0, ey);

    if (presentZ) presentZ->insert(x);
    if (presentPt) presentPt->insert({r.pt_lo, r.pt_hi});
  }

  // sort points in x (per graph)
  for (auto& kv : gs.g) {
    auto* gr = kv.second;
    // Simple bubble-sort for small N
    int N = gr->GetN();
    for (int i=0;i<N;i++) for (int j=i+1;j<N;j++) {
      double xi, yi, xj, yj;
      gr->GetPoint(i, xi, yi);
      gr->GetPoint(j, xj, yj);
      if (xj < xi) {
        double exi = gr->GetErrorX(i), eyi = gr->GetErrorY(i);
        double exj = gr->GetErrorX(j), eyj = gr->GetErrorY(j);
        gr->SetPoint(i, xj, yj); gr->SetPointError(i, exj, eyj);
        gr->SetPoint(j, xi, yi); gr->SetPointError(j, exi, eyi);
      }
    }
  }

  return gs;
}

// Fill graphs for y vs pt (one curve per z bin)
static GraphSet BuildGraphsVsPt(const std::vector<Row>& rows, const PlotConfig& cfg,
                               const std::string& yField, const std::string& eyField,
                               std::set<double>* presentZ=nullptr,
                               std::set<std::pair<double,double>>* presentPt=nullptr) {
  GraphSet gs;
  for (const auto& r : rows) {
    if (!PassRowPhysics(r,cfg)) continue;

    double x = r.pt_center;
    if (!IsFinite(x) || x==kMissing) continue;

    double y = kMissing, ey = 0.0;
    if (yField=="M0") y = r.M0;
    else if (yField=="A") y = r.A;
    else if (yField=="B") y = r.B;
    else if (yField=="chi2_ndf") y = r.chi2_ndf;
    else if (yField=="prob") y = r.prob;
    else if (yField=="fit_tier") y = (double)r.fit_tier;
    else if (yField=="npts_fit") y = (double)r.npts_fit;

    if (eyField=="M0_err") ey = r.M0_err;
    else if (eyField=="A_err")  ey = r.A_err;
    else if (eyField=="B_err")  ey = r.B_err;
    else ey = 0.0;

    if (!IsFinite(y) || y==kMissing) continue;
    if (!IsFinite(ey) || ey==kMissing) ey = 0.0;

    // Key consistently by z-bin. Prefer explicit z_bin if present, otherwise infer from z_center.
    int key = ZKeyFromRow(r);
    auto it = gs.g.find(key);
    if (it == gs.g.end()) {
      auto* gr = new TGraphErrors();
      gr->SetName(Form("g_%s_vs_pt_zbin%d", yField.c_str(), key));
      gs.g[key] = gr;
      it = gs.g.find(key);
    }

    auto* gr = it->second;
    int n = gr->GetN();
    gr->SetPoint(n, x, y);
    gr->SetPointError(n, 0.0, ey);

    if (presentZ) presentZ->insert(r.z_center);
    if (presentPt) presentPt->insert({r.pt_lo, r.pt_hi});
  }

  // sort points in x (per graph)
  for (auto& kv : gs.g) {
    auto* gr = kv.second;
    int N = gr->GetN();
    for (int i=0;i<N;i++) for (int j=i+1;j<N;j++) {
      double xi, yi, xj, yj;
      gr->GetPoint(i, xi, yi);
      gr->GetPoint(j, xj, yj);
      if (xj < xi) {
        double exi = gr->GetErrorX(i), eyi = gr->GetErrorY(i);
        double exj = gr->GetErrorX(j), eyj = gr->GetErrorY(j);
        gr->SetPoint(i, xj, yj); gr->SetPointError(i, exj, eyj);
        gr->SetPoint(j, xi, yi); gr->SetPointError(j, exi, eyi);
      }
    }
  }

  return gs;
}

static void StyleGraph(TGraphErrors* gr, int idx, double msize) {
  int ci = kColors[idx % (int)(sizeof(kColors)/sizeof(kColors[0]))];
  int mi = kMarkers[idx % (int)(sizeof(kMarkers)/sizeof(kMarkers[0]))];
  gr->SetMarkerColor(ci);
  gr->SetLineColor(ci);
  gr->SetMarkerStyle(mi);
  gr->SetMarkerSize(msize);
  gr->SetLineWidth(2);
}

static void CollectAllY(const GraphSet& gs, std::vector<double>& y, std::vector<double>& ey) {
  y.clear(); ey.clear();
  for (const auto& kv : gs.g) {
    const auto* gr = kv.second;
    int N = gr->GetN();
    for (int i=0;i<N;i++) {
      double x, yi;
      gr->GetPoint(i, x, yi);
      y.push_back(yi);
      ey.push_back(gr->GetErrorY(i));
    }
  }
}

static void Plot1D_MultiCurves(const std::string& outPng,
                              const PlotConfig& cfg,
                              GraphSet& gs,
                              const std::string& title,
                              const std::string& xTitle,
                              const std::string& yTitle,
                              bool logy,
                              const std::string& missingNote,
                              bool legendIsPtBin,
                              const std::map<int, Row>& exemplarRowByKey) {
  ApplyCanvasStyle();

  auto* c = new TCanvas("c","c",1100,850);
  c->SetLeftMargin(cfg.left_margin);
  c->SetRightMargin(cfg.right_margin);
  c->SetBottomMargin(cfg.bottom_margin);
  c->SetTopMargin(cfg.top_margin);
  if (logy) c->SetLogy(true);

  // axis frame
  double xmin=0, xmax=1, ymin=0, ymax=1;
  bool first=true;
  for (auto& kv : gs.g) {
    auto* gr = kv.second;
    int N=gr->GetN();
    for (int i=0;i<N;i++) {
      double x,y;
      gr->GetPoint(i,x,y);
      xmin = first?x:std::min(xmin,x);
      xmax = first?x:std::max(xmax,x);
      ymin = first?y:std::min(ymin,y);
      ymax = first?y:std::max(ymax,y);
      first=false;
    }
  }
  if (first) {
    std::cerr << "WARN: no points for plot " << outPng << "\n";
    delete c;
    return;
  }

  // Compute y-range
  std::vector<double> allY, allEy;
  CollectAllY(gs, allY, allEy);
  std::pair<double,double> yr;
  if ((yTitle=="A") || (yTitle=="B")) {
    yr = ComputeYRangeSymmetric0(allY, allEy, cfg.symmetric_AB_min_abs);
  } else {
    yr = ComputeYRangeAuto(allY, allEy);
  }

  // Expand x a bit
  double xpad = 0.06*(xmax-xmin);
  if (xpad < 1e-3) xpad = 0.05;

  // In logy, enforce positive min
  if (logy) {
    double minPos = 1e99;
    for (double v : allY) if (IsFinite(v) && v>0) minPos = std::min(minPos, v);
    if (minPos < 1e99) {
      yr.first = std::max(yr.first, minPos/10.0);
      yr.second = std::max(yr.second, minPos*10.0);
    }
  }

  auto* frame = c->DrawFrame(xmin-xpad, yr.first, xmax+xpad, yr.second);
  frame->GetXaxis()->SetTitle(xTitle.c_str());
  frame->GetYaxis()->SetTitle(yTitle.c_str());
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->GetXaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetLabelSize(0.040);

  // Legend and header
  DrawHeader(cfg, title);
  DrawMissingNote(cfg, missingNote);
  auto* leg = MakeTopLegend(cfg, 0.55, 1.0-cfg.top_margin+0.02, 0.95, 0.98);

  // Draw graphs
  int idx=0;
  for (auto& kv : gs.g) {
    auto* gr = kv.second;
    StyleGraph(gr, idx, cfg.marker_size);
    gr->Draw("P SAME");

    std::string lab;
    auto itR = exemplarRowByKey.find(kv.first);
    if (itR != exemplarRowByKey.end()) {
      lab = legendIsPtBin ? PtLegendLabelFromAnyRow(itR->second) : ZLegendLabelFromAnyRow(itR->second);
    } else {
      if (legendIsPtBin) {
        // Fallback: show the physical pT bin edges
        lab = PtBinLabelFromIndex(kv.first);
      } else {
        // Fallback: show the physical z center
        const int k = kv.first;
        if (k >= 0 && k < (int)kExpectedZCenters.size()) lab = Form("z = %.2f", kExpectedZCenters[k]);
        else lab = Form("zbin %d", k);
      }
    }
    leg->AddEntry(gr, lab.c_str(), "p");
    ++idx;
  }
  leg->Draw();

  c->SaveAs(outPng.c_str());
  delete leg;
  delete c;
}

// Build exemplar rows for legend labels
static std::map<int, Row> ExemplarByPtBin(const std::vector<Row>& rows, const PlotConfig& cfg) {
  std::map<int, Row> out;
  for (const auto& r : rows) {
    if (!PassRowPhysics(r,cfg)) continue;
    const int k = PtKeyFromRow(r);
    if (k < 0) continue;
    if (out.find(k)==out.end()) out[k]=r;
  }
  return out;
}

static std::map<int, Row> ExemplarByZBin(const std::vector<Row>& rows, const PlotConfig& cfg) {
  std::map<int, Row> out;
  for (const auto& r : rows) {
    if (!PassRowPhysics(r,cfg)) continue;
    const int k = ZKeyFromRow(r);
    if (k < 0) continue;
    if (out.find(k)==out.end()) out[k]=r;
  }
  return out;
}

} // namespace

// --- main entry ---
void PlotFitParameters(const char* csvPathC="") {
  std::string csvPath = csvPathC ? std::string(csvPathC) : std::string();
  if (csvPath.empty()) {
    std::cerr << "Usage: PlotFitParameters(\"results/.../tables/fit_parameters_<mode>_tierOn.csv\")\n";
    return;
  }

  PlotConfig cfg;
  cfg.require_valid_fit_default = true; // physics default

  auto rows = ReadRows(csvPath);
  if (rows.empty()) return;

  std::string outDir = InferOutDir(csvPath);

  // legend exemplars
  auto exPt = ExemplarByPtBin(rows, cfg);
  auto exZ  = ExemplarByZBin(rows, cfg);

  // --- Core physics plots ---
  {
    // M0 vs z in pt bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsZ(rows, cfg, "M0", "M0_err", &presentZ, &presentPt);
    std::string note = MissingZNote(presentZ);
    Plot1D_MultiCurves(outDir+"/M0_vs_z.png", cfg, gs, "M0_vs_z", "z", "M0", false, note, true, exPt);
    DeleteGraphs(gs);
  }
  {
    // A vs z in pt bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsZ(rows, cfg, "A", "A_err", &presentZ, &presentPt);
    std::string note = MissingZNote(presentZ);
    Plot1D_MultiCurves(outDir+"/A_vs_z.png", cfg, gs, "A_vs_z", "z", "A", false, note, true, exPt);
    DeleteGraphs(gs);
  }
  {
    // B vs z in pt bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsZ(rows, cfg, "B", "B_err", &presentZ, &presentPt);
    std::string note = MissingZNote(presentZ);
    Plot1D_MultiCurves(outDir+"/B_vs_z.png", cfg, gs, "B_vs_z", "z", "B", false, note, true, exPt);
    DeleteGraphs(gs);
  }
  {
    // M0 vs pt in z bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsPt(rows, cfg, "M0", "M0_err", &presentZ, &presentPt);
    std::string note = MissingPtNote(presentPt);
    Plot1D_MultiCurves(outDir+"/M0_vs_pt.png", cfg, gs, "M0_vs_pt", "p_{T} (GeV)", "M0", false, note, false, exZ);
    DeleteGraphs(gs);
  }
  {
    // A vs pt in z bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsPt(rows, cfg, "A", "A_err", &presentZ, &presentPt);
    std::string note = MissingPtNote(presentPt);
    Plot1D_MultiCurves(outDir+"/A_vs_pt.png", cfg, gs, "A_vs_pt", "p_{T} (GeV)", "A", false, note, false, exZ);
    DeleteGraphs(gs);
  }
  {
    // B vs pt in z bins
    std::set<double> presentZ; std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsPt(rows, cfg, "B", "B_err", &presentZ, &presentPt);
    std::string note = MissingPtNote(presentPt);
    Plot1D_MultiCurves(outDir+"/B_vs_pt.png", cfg, gs, "B_vs_pt", "p_{T} (GeV)", "B", false, note, false, exZ);
    DeleteGraphs(gs);
  }
  {
    // chi2/ndf vs z in pt bins (log)
    std::set<double> presentZ;
    auto gs = BuildGraphsVsZ(rows, cfg, "chi2_ndf", "", &presentZ, nullptr);
    std::string note = MissingZNote(presentZ);
    Plot1D_MultiCurves(outDir+"/chi2ndf_vs_z.png", cfg, gs, "chi2ndf_vs_z", "z", "#chi^{2}/ndf", true, note, true, exPt);
    DeleteGraphs(gs);
  }
  {
    // chi2/ndf vs pt in z bins (log)
    std::set<std::pair<double,double>> presentPt;
    auto gs = BuildGraphsVsPt(rows, cfg, "chi2_ndf", "", nullptr, &presentPt);
    std::string note = MissingPtNote(presentPt);
    Plot1D_MultiCurves(outDir+"/chi2ndf_vs_pt.png", cfg, gs, "chi2ndf_vs_pt", "p_{T} (GeV)", "#chi^{2}/ndf", true, note, false, exZ);
    DeleteGraphs(gs);
  }

  // --- Diagnostics (kept, but with same top-margin layout) ---
  {
    // dominant weight fraction is in the CSV as dominant_wfrac; not parsed in this macro.
    // If you want it, add it to Row and parser.
  }

  std::cout << "Wrote FitParamPlots to: " << outDir << "\n";
}
