// PlotCoinXsecFromTable.C (rsidis_xs_v5)
//
// Reads table CSV(s) produced by TableCoinXsec.C and generates:
//   - xsec vs phipq plots with pads = pT bins
//   - curves = z (single mode) OR z-settings from group list (group mode)
//   - tiered physics fits: sigma = M0*(1 + A cos(phi) + B cos(2phi)) with tiering
//   - fit_parameters_single.csv / fit_parameters_group.csv
//   - PNG outputs:
//       results/<setting_id>/PNGs/xsec_phipq_z_pt_overlayed_single.png
//       results/<group_id>/PNGs/xsec_phipq_z_pt_overlayed_group.png
//
// Semantics:
//   * Single mode: z is binned (z_lo/z_hi from table), curves are those z-bins.
//   * Group  mode: z is the fixed setting value encoded in curve_label (e.g. z0p36),
//                  curves are those z-settings, labeled cleanly as "z = 0.xx".
//     (No z-range in group mode.)
//
// Fitting philosophy updates (per your request):
//   * Plot points use DEFAULT filter (valid_default etc).
//   * Fits use a stricter, FIT-ONLY filter + robustness guards:
//       - keep only points with rel_err = (xsec_err/xsec) < 0.60 for fit
//       - apply an error floor: sigma_i <- max(sigma_i, 0.25 * median(sigma))
//       - if a single point dominates weights after the floor (dominant_wfrac > 0.60),
//         drop that single most-influential point and refit once
//   * Do NOT draw fit curve unless valid_fit_default==1 (i.e., trustworthy fit).
//
// Run examples:
//   root -l -b -q 'macros/PlotCoinXsecFromTable.C("settings/.../manifest.txt","results","settings")'
//   root -l -b -q 'macros/PlotCoinXsecFromTable.C("groups/.../grp_*.list","results","settings")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <cmath>
#include <cctype>

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

namespace {

static const double kMissing = -999.0;

static std::string NormalizeSlashes(std::string s) {
  std::string out;
  bool prev = false;
  for (char c : s) {
    if (c == '/') {
      if (!prev) out.push_back(c);
      prev = true;
    } else {
      out.push_back(c);
      prev = false;
    }
  }
  if (out.size() > 1 && out.back() == '/') out.pop_back();
  return out;
}

static std::string Dirname(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return p.substr(0, pos);
}

static std::string Basename(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return p;
  return p.substr(pos+1);
}

static bool EndsWith(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

static std::string StripExtension(const std::string& path) {
  std::string base = Basename(path);
  auto pos = base.find_last_of('.');
  if (pos == std::string::npos) return base;
  return base.substr(0, pos);
}

// settingsRoot can be relative ("settings") or absolute ("/.../settings")
static std::string MakeSettingIdFromManifestPath(const std::string& manifestPath,
                                                const std::string& settingsRoot) {
  std::string mdir  = NormalizeSlashes(Dirname(manifestPath));
  std::string sroot = NormalizeSlashes(settingsRoot);

  if (mdir.rfind(sroot + "/", 0) == 0) return mdir.substr(sroot.size() + 1);

  std::string sbase = Basename(sroot); // usually "settings"
  auto pos = mdir.find("/" + sbase + "/");
  if (pos != std::string::npos) return mdir.substr(pos + sbase.size() + 2);
  if (mdir.rfind(sbase + "/", 0) == 0) return mdir.substr(sbase.size() + 1);

  return Basename(mdir);
}

static void MkdirP(const std::string& path) {
  if (path.empty()) return;
  gSystem->mkdir(path.c_str(), true);
}

static std::vector<std::string> SplitCSVLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool inq = false;
  for (size_t i=0;i<line.size();i++) {
    char c = line[i];
    if (c == '"') { inq = !inq; continue; }
    if (c == ',' && !inq) { out.push_back(cur); cur.clear(); continue; }
    cur.push_back(c);
  }
  out.push_back(cur);
  // trim
  for (auto& s : out) {
    auto issp = [](unsigned char c){ return std::isspace(c)!=0; };
    while (!s.empty() && issp((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && issp((unsigned char)s.back())) s.pop_back();
  }
  return out;
}

static bool ToInt(const std::string& s, int& out) {
  if (s.empty()) return false;
  char* endp = nullptr;
  long v = std::strtol(s.c_str(), &endp, 10);
  if (endp == s.c_str()) return false;
  out = (int)v;
  return true;
}

static bool ToLL(const std::string& s, long long& out) {
  if (s.empty()) return false;
  char* endp = nullptr;
  long long v = std::strtoll(s.c_str(), &endp, 10);
  if (endp == s.c_str()) return false;
  out = v;
  return true;
}

static bool ToD(const std::string& s, double& out) {
  if (s.empty()) return false;
  char* endp = nullptr;
  double v = std::strtod(s.c_str(), &endp);
  if (endp == s.c_str()) return false;
  out = v;
  return true;
}

static double WrapPhi(double phi) {
  const double twopi = TMath::TwoPi();
  if (!std::isfinite(phi)) return phi;
  while (phi < 0) phi += twopi;
  while (phi >= twopi) phi -= twopi;
  return phi;
}

// Parse z-setting from curve_label like "z0p36" -> 0.36.
static double ParseZSetting(std::string s) {
  if (s.empty()) return kMissing;
  if (s[0] == 'z' || s[0] == 'Z') s = s.substr(1);
  for (auto& c : s) {
    if (c == 'p') c = '.';
    if (c == 'm') c = '-';
  }
  double z = kMissing;
  if (!ToD(s, z)) return kMissing;
  return z;
}

// median for positive finite values
static double MedianPositive(const std::vector<double>& v) {
  std::vector<double> a;
  a.reserve(v.size());
  for (double x : v) {
    if (std::isfinite(x) && x > 0) a.push_back(x);
  }
  if (a.empty()) return kMissing;
  std::sort(a.begin(), a.end());
  size_t n = a.size();
  if (n % 2 == 1) return a[n/2];
  return 0.5*(a[n/2 - 1] + a[n/2]);
}

// ---------- Fit flag bits ----------
enum FitFlagBits : unsigned int {
  FIT_STATUS_BAD      = 1u << 0,
  NONFINITE_PARAM     = 1u << 1,
  NONFINITE_ERR       = 1u << 2,
  NDF_NONPOSITIVE     = 1u << 3,
  COV_FAILED          = 1u << 4,
  COV_NOT_POSDEF      = 1u << 5,
  CHI2NDF_HUGE        = 1u << 6,
  PROB_TINY           = 1u << 7,
  PHI_UNIFORM_FAIL    = 1u << 8,
  PARAM_AT_LIMIT      = 1u << 9,
  M0_BAD              = 1u << 10,
  M0_RELERR_HUGE      = 1u << 11,
  A_BAD               = 1u << 12,
  A_ERR_BAD           = 1u << 13,
  A_SIG_LOW           = 1u << 14,
  B_BAD               = 1u << 15,
  B_ERR_BAD           = 1u << 16,
  B_SIG_LOW           = 1u << 17,
  DOMINANT_WEIGHT_PT  = 1u << 18,
  NO_FIT_POINTS       = 1u << 19,
  ROBUST_DROPPED_PT   = 1u << 20
};

static void AppendNote(std::string& note, const char* tok) {
  if (!note.empty()) note.push_back(';');
  note += tok;
}

// ---------- Config ----------
struct PlotConfig {
  // plot-point filtering (DEFAULT)
  bool use_valid_default = true;
  bool allow_negative_xsec = false;

  // Fit-only filter + robustness (defaults you approved)
  double fit_relerr_max = 0.60;
  double fit_sigma_floor_frac = 0.25;
  double robust_drop_dom_wfrac_gt = 0.60;

  // fit tier thresholds
  double tier_span0_max = 1.0;          // <1 rad => tier 0
  double tier_span1_max = TMath::Pi();  // <pi => tier 1; else tier 2

  // coverage uniformity
  int phi_nsectors = 8;
  int phi_occ_min_tier1 = 3;
  int phi_occ_min_tier2 = 4;

  // robustness: also hard reject if still dominated after dropping
  double max_weight_frac_reject = 0.80;

  // fit-quality thresholds
  double chi2ndf_max = 10.0;
  double prob_min = 1e-6;

  // "at limit" tolerance for A,B
  double par_limit_tol = 0.98;

  // parameter sanity
  double M0_min = 0.0;
  double M0_relerr_max = 1.0;

  double A_abs_max_plot = 1.2;
  double A_err_max = 1.0;
  double A_sig_min = 1.0;

  double B_abs_max_plot = 1.2;
  double B_err_max = 1.0;
  double B_sig_min = 1.0;

  // labels
  TString xTitle = "#phi_{pq} (rad)";
  TString yTitle = "d#sigma";

  // output
  bool save_png = true;
  bool save_fit_csv = true;
};

// ---------- Row ----------
struct Row {
  std::string mode;
  std::string group_id;
  std::string curve_label;
  std::string setting_id;

  int pt_bin = -1;
  int z_bin = -1;
  int phi_bin = -1;

  double pt_lo=kMissing, pt_hi=kMissing, pt_center=kMissing;
  double z_lo=kMissing, z_hi=kMissing, z_center=kMissing;
  double phi_center=kMissing;

  long long n_sim=0;
  double xsec=kMissing, xsec_err=kMissing;
  double rel_xsec_err=kMissing;
  unsigned int flag_bits=0;
  int valid_default=0;
};

// ---------- Fit record ----------
struct FitOut {
  std::string mode;
  std::string group_id;
  std::string curve_label;
  std::string setting_id;

  int pt_bin=-1;
  int z_bin=-1;
  double pt_lo=kMissing, pt_hi=kMissing, pt_center=kMissing;
  double z_lo=kMissing, z_hi=kMissing, z_center=kMissing;

  int npts=0;          // plotted points
  int npts_fit=0;      // points used in fit (after fit-only filter)
  double phi_min=kMissing, phi_max=kMissing, phi_span=kMissing;
  int phi_occ=0;
  int phi_nsectors=0;
  double dominant_wfrac=kMissing;

  int fit_tier=-1;
  int fit_status=999;
  int cov_status=-1;
  double chi2=kMissing;
  int ndf=-999;
  double prob=kMissing;
  double chi2_ndf=kMissing;

  double M0=kMissing, M0_err=kMissing;
  double A=kMissing, A_err=kMissing;
  double B=kMissing, B_err=kMissing;
  double A_sig=kMissing, B_sig=kMissing;

  unsigned int fit_flag_bits=0;
  int valid_fit_default=0;
  int valid_M0_default=0;
  int valid_A_default=0;
  int valid_B_default=0;
  std::string note;
};

// ---------- CSV reader ----------
static bool ReadTableCSV(const std::string& csvPath, std::vector<Row>& rows) {
  std::ifstream f(csvPath);
  if (!f.is_open()) return false;

  std::string header;
  if (!std::getline(f, header)) return false;
  auto cols = SplitCSVLine(header);

  std::unordered_map<std::string,int> idx;
  for (int i=0;i<(int)cols.size();i++) idx[cols[i]] = i;

  auto getS = [&](const std::vector<std::string>& v, const char* name)->std::string{
    auto it = idx.find(name);
    if (it==idx.end()) return "";
    int i = it->second;
    return (i>=0 && i<(int)v.size()) ? v[i] : "";
  };
  auto getI = [&](const std::vector<std::string>& v, const char* name, int def=-1)->int{
    int out=def;
    auto s = getS(v,name);
    if (!ToInt(s,out)) return def;
    return out;
  };
  auto getLL = [&](const std::vector<std::string>& v, const char* name, long long def=0)->long long{
    long long out=def;
    auto s = getS(v,name);
    if (!ToLL(s,out)) return def;
    return out;
  };
  auto getD = [&](const std::vector<std::string>& v, const char* name, double def=kMissing)->double{
    double out=def;
    auto s = getS(v,name);
    if (!ToD(s,out)) return def;
    return out;
  };

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    auto v = SplitCSVLine(line);
    if (v.size() < 10) continue;

    Row r;
    r.mode       = getS(v,"mode");
    r.group_id   = getS(v,"group_id");
    r.curve_label= getS(v,"curve_label");
    r.setting_id = getS(v,"setting_id");

    r.pt_bin = getI(v,"pt_bin",-1);
    r.z_bin  = getI(v,"z_bin",-1);
    r.phi_bin= getI(v,"phi_bin",-1);

    r.pt_lo = getD(v,"pt_lo",kMissing);
    r.pt_hi = getD(v,"pt_hi",kMissing);
    r.pt_center = getD(v,"pt_center",kMissing);

    r.z_lo = getD(v,"z_lo",kMissing);
    r.z_hi = getD(v,"z_hi",kMissing);
    r.z_center = getD(v,"z_center",kMissing);

    r.phi_center = WrapPhi(getD(v,"phi_center",kMissing));

    r.n_sim = getLL(v,"n_sim",0);
    r.xsec = getD(v,"xsec",kMissing);
    r.xsec_err = getD(v,"xsec_err",kMissing);
    r.rel_xsec_err = getD(v,"rel_xsec_err",kMissing);

    {
      auto s = getS(v,"flag_bits");
      long long tmp=0;
      if (ToLL(s,tmp)) r.flag_bits = (unsigned int)tmp;
    }
    r.valid_default = getI(v,"valid_default",0);

    rows.push_back(r);
  }
  return true;
}

// ---------- point filter (plot) ----------
static bool PassPointDefault(const Row& r, const PlotConfig& cfg) {
  if (!std::isfinite(r.phi_center)) return false;
  if (!std::isfinite(r.xsec) || !std::isfinite(r.xsec_err)) return false;
  if (r.xsec == kMissing || r.xsec_err == kMissing) return false;
  if (r.xsec_err <= 0.0) return false;

  if (cfg.use_valid_default && r.valid_default != 1) return false;
  if (!cfg.allow_negative_xsec && r.xsec <= 0.0) return false;
  return true;
}

// ---------- curve grouping ----------
struct CurveKey {
  int pt_bin=-1;
  int z_bin=-1;           // single mode uses z_bin; group mode uses -1
  double sort_z=kMissing; // single: z_center; group: z_setting
  std::string legend;     // legend label
  std::string setting_id; // source setting id (group)
};

static bool operator<(const CurveKey& a, const CurveKey& b) {
  if (a.pt_bin != b.pt_bin) return a.pt_bin < b.pt_bin;
  if (a.z_bin != b.z_bin) return a.z_bin < b.z_bin;
  if (std::isfinite(a.sort_z) && std::isfinite(b.sort_z) && a.sort_z != b.sort_z) return a.sort_z < b.sort_z;
  if (a.legend != b.legend) return a.legend < b.legend;
  return a.setting_id < b.setting_id;
}

struct CurvePoints {
  CurveKey key;
  std::vector<double> x, y, ey; // plotted points (already default-filtered)
  double pt_lo=kMissing, pt_hi=kMissing, pt_center=kMissing;
  double z_lo=kMissing, z_hi=kMissing, z_center=kMissing; // single: range; group: only z_center
};

// ---------- uniformity ----------
static int PhiOccupancy(const std::vector<double>& phi, int nsectors) {
  if (nsectors <= 0) return 0;
  std::vector<char> occ(nsectors, 0);
  const double twopi = TMath::TwoPi();
  for (double p : phi) {
    if (!std::isfinite(p)) continue;
    double w = WrapPhi(p);
    int idx = (int)std::floor(w / (twopi / nsectors));
    if (idx < 0) idx = 0;
    if (idx >= nsectors) idx = nsectors - 1;
    occ[idx] = 1;
  }
  int cnt=0;
  for (char b : occ) if (b) cnt++;
  return cnt;
}

// dominant weight fraction for a given error vector (after any floors)
static double DominantWeightFrac(const std::vector<double>& ey) {
  double sum = 0.0;
  double mx = 0.0;
  for (double e : ey) {
    if (!std::isfinite(e) || e <= 0) continue;
    double w = 1.0/(e*e);
    sum += w;
    mx = std::max(mx, w);
  }
  if (sum <= 0) return kMissing;
  return mx / sum;
}

// ---------- fit validity ----------
static void ComputeFitValidity(FitOut& fo, const PlotConfig& cfg) {
  auto noteIf = [&](unsigned int bit, const char* tok){
    if (fo.fit_flag_bits & bit) AppendNote(fo.note, tok);
  };

  // covariance bits
  if (fo.cov_status < 0) {
    fo.fit_flag_bits |= COV_FAILED;
  } else if (fo.cov_status == 0) {
    fo.fit_flag_bits |= COV_FAILED;
  } else if (fo.cov_status == 1 || fo.cov_status == 2) {
    fo.fit_flag_bits |= COV_NOT_POSDEF;
  }

  if (fo.fit_status != 0) fo.fit_flag_bits |= FIT_STATUS_BAD;
  if (fo.ndf <= 0) fo.fit_flag_bits |= NDF_NONPOSITIVE;

  auto fin = [](double v){ return std::isfinite(v) && v != kMissing; };

  if (!fin(fo.M0) || (fo.fit_tier>=1 && !fin(fo.A)) || (fo.fit_tier>=2 && !fin(fo.B))) fo.fit_flag_bits |= NONFINITE_PARAM;
  if (!fin(fo.M0_err) || (fo.fit_tier>=1 && !fin(fo.A_err)) || (fo.fit_tier>=2 && !fin(fo.B_err))) fo.fit_flag_bits |= NONFINITE_ERR;

  if (fin(fo.chi2_ndf) && fo.chi2_ndf > cfg.chi2ndf_max) fo.fit_flag_bits |= CHI2NDF_HUGE;
  if (fin(fo.prob) && fo.prob < cfg.prob_min) fo.fit_flag_bits |= PROB_TINY;

  // coverage uniformity
  if (fo.fit_tier == 1) {
    if (fo.phi_occ < cfg.phi_occ_min_tier1) fo.fit_flag_bits |= PHI_UNIFORM_FAIL;
  } else if (fo.fit_tier == 2) {
    if (fo.phi_occ < cfg.phi_occ_min_tier2) fo.fit_flag_bits |= PHI_UNIFORM_FAIL;
  }

  // hard reject still dominated
  if (fin(fo.dominant_wfrac) && fo.dominant_wfrac > cfg.max_weight_frac_reject) fo.fit_flag_bits |= DOMINANT_WEIGHT_PT;

  // at limit
  if (fo.fit_tier >= 1 && fin(fo.A) && std::abs(fo.A) >= cfg.par_limit_tol) fo.fit_flag_bits |= PARAM_AT_LIMIT;
  if (fo.fit_tier >= 2 && fin(fo.B) && std::abs(fo.B) >= cfg.par_limit_tol) fo.fit_flag_bits |= PARAM_AT_LIMIT;

  // M0 sanity
  if (!fin(fo.M0) || fo.M0 <= cfg.M0_min) fo.fit_flag_bits |= M0_BAD;
  if (fin(fo.M0) && fin(fo.M0_err) && fo.M0 > 0) {
    double r = std::abs(fo.M0_err/fo.M0);
    if (r > cfg.M0_relerr_max) fo.fit_flag_bits |= M0_RELERR_HUGE;
  }

  // A sanity
  if (fo.fit_tier >= 1) {
    if (!fin(fo.A) || std::abs(fo.A) > cfg.A_abs_max_plot) fo.fit_flag_bits |= A_BAD;
    if (!fin(fo.A_err) || fo.A_err <= 0 || fo.A_err > cfg.A_err_max) fo.fit_flag_bits |= A_ERR_BAD;
    if (fin(fo.A_sig) && fo.A_sig < cfg.A_sig_min) fo.fit_flag_bits |= A_SIG_LOW;
  }

  // B sanity
  if (fo.fit_tier >= 2) {
    if (!fin(fo.B) || std::abs(fo.B) > cfg.B_abs_max_plot) fo.fit_flag_bits |= B_BAD;
    if (!fin(fo.B_err) || fo.B_err <= 0 || fo.B_err > cfg.B_err_max) fo.fit_flag_bits |= B_ERR_BAD;
    if (fin(fo.B_sig) && fo.B_sig < cfg.B_sig_min) fo.fit_flag_bits |= B_SIG_LOW;
  }

  // notes
  noteIf(NO_FIT_POINTS,     "NO_FIT_POINTS");
  noteIf(ROBUST_DROPPED_PT, "ROBUST_DROPPED_PT");
  noteIf(FIT_STATUS_BAD,    "FIT_STATUS_BAD");
  noteIf(NDF_NONPOSITIVE,   "NDF_NONPOSITIVE");
  noteIf(COV_FAILED,        "COV_FAILED");
  noteIf(COV_NOT_POSDEF,    "COV_NOT_POSDEF");
  noteIf(NONFINITE_PARAM,   "NONFINITE_PARAM");
  noteIf(NONFINITE_ERR,     "NONFINITE_ERR");
  noteIf(CHI2NDF_HUGE,      "CHI2NDF_HUGE");
  noteIf(PROB_TINY,         "PROB_TINY");
  noteIf(PHI_UNIFORM_FAIL,  "PHI_UNIFORM_FAIL");
  noteIf(DOMINANT_WEIGHT_PT,"DOMINANT_WEIGHT_PT");
  noteIf(PARAM_AT_LIMIT,    "PARAM_AT_LIMIT");
  noteIf(M0_BAD,            "M0_BAD");
  noteIf(M0_RELERR_HUGE,    "M0_RELERR_HUGE");
  noteIf(A_BAD,             "A_BAD");
  noteIf(A_ERR_BAD,         "A_ERR_BAD");
  noteIf(A_SIG_LOW,         "A_SIG_LOW");
  noteIf(B_BAD,             "B_BAD");
  noteIf(B_ERR_BAD,         "B_ERR_BAD");
  noteIf(B_SIG_LOW,         "B_SIG_LOW");

  const unsigned int catastrophic = FIT_STATUS_BAD | NONFINITE_PARAM | NONFINITE_ERR | NDF_NONPOSITIVE | COV_FAILED | NO_FIT_POINTS;
  const unsigned int defaultReject = catastrophic | COV_NOT_POSDEF | CHI2NDF_HUGE | PROB_TINY | PHI_UNIFORM_FAIL | DOMINANT_WEIGHT_PT;

  fo.valid_fit_default = ((fo.fit_flag_bits & defaultReject) == 0u) ? 1 : 0;
  fo.valid_M0_default = (fo.valid_fit_default == 1 && (fo.fit_flag_bits & (M0_BAD | M0_RELERR_HUGE)) == 0u) ? 1 : 0;
  fo.valid_A_default  = (fo.valid_fit_default == 1 && fo.fit_tier>=1 && (fo.fit_flag_bits & (A_BAD | A_ERR_BAD | A_SIG_LOW)) == 0u) ? 1 : 0;
  fo.valid_B_default  = (fo.valid_fit_default == 1 && fo.fit_tier>=2 && (fo.fit_flag_bits & (B_BAD | B_ERR_BAD | B_SIG_LOW)) == 0u) ? 1 : 0;
}

// choose tier based on fit-point span and fit-point count, then downgrade by occupancy
static int ChooseTier(int npts_fit, double phi_span_fit, int phi_occ, const PlotConfig& cfg) {
  int tier = 2;
  if (phi_span_fit < cfg.tier_span0_max) tier = 0;
  else if (phi_span_fit < cfg.tier_span1_max) tier = 1;
  else tier = 2;

  if (tier == 2 && npts_fit < 4) tier = (npts_fit >= 3) ? 1 : 0;
  if (tier == 1 && npts_fit < 3) tier = 0;

  // Pre-downgrade by coverage
  if (tier == 2 && phi_occ < cfg.phi_occ_min_tier2) tier = 1;
  if (tier == 1 && phi_occ < cfg.phi_occ_min_tier1) tier = 0;

  return tier;
}

// build fit vectors from plotted points with fit-only filtering and error floor
static void BuildFitVectors(const CurvePoints& cp, const PlotConfig& cfg,
                            std::vector<double>& xf, std::vector<double>& yf, std::vector<double>& ef,
                            bool& didDrop, int& droppedIndex)
{
  xf.clear(); yf.clear(); ef.clear();
  didDrop = false;
  droppedIndex = -1;

  // fit-only relerr cut
  for (size_t i=0;i<cp.x.size();i++) {
    double x = cp.x[i];
    double y = cp.y[i];
    double e = cp.ey[i];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(e)) continue;
    if (y <= 0 || e <= 0) continue;
    double rel = e / y;
    if (std::isfinite(rel) && rel > cfg.fit_relerr_max) continue;
    xf.push_back(x);
    yf.push_back(y);
    ef.push_back(e);
  }

  if (ef.empty()) return;

  // error floor
  double med = MedianPositive(ef);
  if (std::isfinite(med) && med > 0) {
    double floor = cfg.fit_sigma_floor_frac * med;
    for (double& e : ef) if (e < floor) e = floor;
  }

  // robust: if dominated, drop most influential point and refit once
  double dom = DominantWeightFrac(ef);
  if (std::isfinite(dom) && dom > cfg.robust_drop_dom_wfrac_gt && ef.size() >= 4) {
    // find index of max weight = min ef
    size_t imax = 0;
    double best = 1e99;
    for (size_t i=0;i<ef.size();i++) {
      if (ef[i] < best) { best = ef[i]; imax = i; }
    }
    droppedIndex = (int)imax;
    didDrop = true;
    xf.erase(xf.begin()+imax);
    yf.erase(yf.begin()+imax);
    ef.erase(ef.begin()+imax);
  }
}

// ---------- fit ----------
static FitOut FitCurve(const CurvePoints& cp, const PlotConfig& cfg, TF1*& outFunc) {
  FitOut fo;
  fo.curve_label = cp.key.legend;
  fo.setting_id  = cp.key.setting_id;
  fo.pt_bin = cp.key.pt_bin;
  fo.z_bin  = cp.key.z_bin;
  fo.pt_lo = cp.pt_lo; fo.pt_hi = cp.pt_hi; fo.pt_center = cp.pt_center;
  fo.z_lo  = cp.z_lo;  fo.z_hi  = cp.z_hi;  fo.z_center = cp.z_center;

  fo.npts = (int)cp.x.size();
  fo.phi_nsectors = cfg.phi_nsectors;
  outFunc = nullptr;

  // Build fit-only vectors
  std::vector<double> xf, yf, ef;
  bool didDrop=false;
  int droppedIndex=-1;
  BuildFitVectors(cp, cfg, xf, yf, ef, didDrop, droppedIndex);
  fo.npts_fit = (int)xf.size();

  // span/occ for fit points
  if (fo.npts_fit >= 1) {
    double xmin=1e99, xmax=-1e99;
    for (double x : xf) { xmin = std::min(xmin, x); xmax = std::max(xmax, x); }
    fo.phi_min = xmin;
    fo.phi_max = xmax;
    fo.phi_span = (xmax > xmin) ? (xmax - xmin) : 0.0;
    fo.phi_occ = PhiOccupancy(xf, cfg.phi_nsectors);
    fo.dominant_wfrac = DominantWeightFrac(ef);
  }

  if (didDrop) {
    fo.fit_flag_bits |= ROBUST_DROPPED_PT;
  }

  if (fo.npts_fit < 2) {
    fo.fit_flag_bits |= NO_FIT_POINTS;
    fo.fit_status = 999;
    fo.cov_status = -1;
    fo.ndf = -999;
    ComputeFitValidity(fo, cfg);
    return fo;
  }

  // pick tier using fit points
  fo.fit_tier = ChooseTier(fo.npts_fit, fo.phi_span, fo.phi_occ, cfg);

  // compute mean for initial M0
  double ysum=0.0, yMin=1e99, yMax=-1e99;
  for (double y : yf) { ysum += y; yMin = std::min(yMin, y); yMax = std::max(yMax, y); }
  double ymean = ysum / (double)fo.npts_fit;

  // Build graph for fit
  std::vector<double> ex(fo.npts_fit, 0.0);
  TGraphErrors gfit(fo.npts_fit, xf.data(), yf.data(), ex.data(), ef.data());

  TF1* f = nullptr;
  if (fo.fit_tier == 0) {
    f = new TF1(Form("f0_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()), "[0]", fo.phi_min, fo.phi_max);
    f->SetParameters(ymean);
    f->SetParLimits(0, 0.0, 1e9);
  } else if (fo.fit_tier == 1) {
    f = new TF1(Form("f1_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()),
                "[0]*(1 + [1]*cos(x))", fo.phi_min, fo.phi_max);
    double A0 = 0.0;
    if (yMax + yMin > 0) {
      A0 = (yMax - yMin) / (yMax + yMin);
      A0 = std::max(-0.9, std::min(0.9, A0));
    }
    f->SetParameters(ymean, A0);
    f->SetParLimits(0, 0.0, 1e9);
    f->SetParLimits(1, -1.0, 1.0);
  } else {
    f = new TF1(Form("f2_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()),
                "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", fo.phi_min, fo.phi_max);
    double A0 = 0.0;
    if (yMax + yMin > 0) {
      A0 = (yMax - yMin) / (yMax + yMin);
      A0 = std::max(-0.9, std::min(0.9, A0));
    }
    f->SetParameters(ymean, A0, 0.0);
    f->SetParLimits(0, 0.0, 1e9);
    f->SetParLimits(1, -1.0, 1.0);
    f->SetParLimits(2, -1.0, 1.0);
  }

  // Fit
  TFitResultPtr r = gfit.Fit(f, "RQS0");
  fo.fit_status = (int)r;
  if (r.Get()) {
    fo.cov_status = r->CovMatrixStatus();
    fo.chi2 = r->Chi2();
    fo.ndf  = r->Ndf();
    fo.prob = r->Prob();
  } else {
    fo.cov_status = -1;
    fo.chi2 = f->GetChisquare();
    fo.ndf  = f->GetNDF();
    fo.prob = f->GetProb();
  }
  if (fo.ndf > 0) fo.chi2_ndf = fo.chi2 / (double)fo.ndf;

  fo.M0 = f->GetParameter(0);
  fo.M0_err = f->GetParError(0);
  if (fo.fit_tier >= 1) { fo.A = f->GetParameter(1); fo.A_err = f->GetParError(1); }
  if (fo.fit_tier >= 2) { fo.B = f->GetParameter(2); fo.B_err = f->GetParError(2); }

  if (std::isfinite(fo.A) && std::isfinite(fo.A_err) && fo.A_err > 0) fo.A_sig = std::abs(fo.A)/fo.A_err;
  if (std::isfinite(fo.B) && std::isfinite(fo.B_err) && fo.B_err > 0) fo.B_sig = std::abs(fo.B)/fo.B_err;

  auto sane = [&](double v){ return std::isfinite(v) ? v : kMissing; };
  fo.M0 = sane(fo.M0); fo.M0_err = sane(fo.M0_err);
  fo.A  = sane(fo.A);  fo.A_err  = sane(fo.A_err);
  fo.B  = sane(fo.B);  fo.B_err  = sane(fo.B_err);
  fo.chi2 = sane(fo.chi2);
  fo.prob = sane(fo.prob);
  fo.chi2_ndf = sane(fo.chi2_ndf);
  fo.A_sig = sane(fo.A_sig);
  fo.B_sig = sane(fo.B_sig);

  ComputeFitValidity(fo, cfg);

  outFunc = f;
  return fo;
}

// ---------- fit CSV ----------
static void WriteFitHeader(std::ofstream& out) {
  out
    << "mode,group_id,curve_label,setting_id,pt_bin,z_bin,"
    << "pt_lo,pt_hi,pt_center,z_lo,z_hi,z_center,"
    << "npts,npts_fit,phi_min,phi_max,phi_span,phi_occ,phi_nsectors,dominant_wfrac,"
    << "M0,M0_err,A,A_err,B,B_err,A_sig,B_sig,"
    << "chi2,ndf,chi2_ndf,prob,fit_tier,fit_status,cov_status,"
    << "fit_flag_bits,valid_fit_default,valid_M0_default,valid_A_default,valid_B_default,"
    << "note\n";
}

static void WriteFitRow(std::ofstream& out, const FitOut& fo) {
  out
    << fo.mode << ","
    << fo.group_id << ","
    << fo.curve_label << ","
    << fo.setting_id << ","
    << fo.pt_bin << ","
    << fo.z_bin << ","
    << fo.pt_lo << "," << fo.pt_hi << "," << fo.pt_center << ","
    << fo.z_lo << "," << fo.z_hi << "," << fo.z_center << ","
    << fo.npts << "," << fo.npts_fit << ","
    << fo.phi_min << "," << fo.phi_max << "," << fo.phi_span << ","
    << fo.phi_occ << "," << fo.phi_nsectors << "," << fo.dominant_wfrac << ","
    << fo.M0 << "," << fo.M0_err << ","
    << fo.A << "," << fo.A_err << ","
    << fo.B << "," << fo.B_err << ","
    << fo.A_sig << "," << fo.B_sig << ","
    << fo.chi2 << "," << fo.ndf << "," << fo.chi2_ndf << "," << fo.prob << ","
    << fo.fit_tier << "," << fo.fit_status << "," << fo.cov_status << ","
    << fo.fit_flag_bits << ","
    << fo.valid_fit_default << ","
    << fo.valid_M0_default << ","
    << fo.valid_A_default << ","
    << fo.valid_B_default << ","
    << fo.note
    << "\n";
}

} // namespace

// ------------------------------------------------------------
// Main macro entry
// ------------------------------------------------------------
void PlotCoinXsecFromTable(const char* manifestOrGroupPath,
                           const char* resultsRoot = "results",
                           const char* settingsRoot = "settings",
                           bool savePNGs = true,
                           bool saveFitCsv = true)
{
  gROOT->SetBatch(kTRUE);

  if (!manifestOrGroupPath || std::string(manifestOrGroupPath).empty()) {
    std::cerr << "ERROR: manifestOrGroupPath is required\n";
    return;
  }

  const std::string p = NormalizeSlashes(std::string(manifestOrGroupPath));
  const bool isGroup = EndsWith(p, ".list");

  std::string id;          // setting_id or group_id
  std::string resultsDir;
  std::string csvPath;
  std::string fitPath;
  std::string pngPath;

  if (isGroup) {
    id = StripExtension(Basename(p)); // group_id
    resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + id);
    csvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_group.csv";
    fitPath = resultsDir + "/tables/fit_parameters_group.csv";
    pngPath = resultsDir + "/PNGs/xsec_phipq_z_pt_overlayed_group.png";
  } else {
    id = MakeSettingIdFromManifestPath(p, NormalizeSlashes(settingsRoot));
    resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + id);
    csvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_single.csv";
    fitPath = resultsDir + "/tables/fit_parameters_single.csv";
    pngPath = resultsDir + "/PNGs/xsec_phipq_z_pt_overlayed_single.png";
  }

  MkdirP(resultsDir);
  MkdirP(resultsDir + "/tables");
  MkdirP(resultsDir + "/PNGs");

  PlotConfig cfg;
  cfg.save_png = savePNGs;
  cfg.save_fit_csv = saveFitCsv;

  std::vector<Row> rows;
  if (!ReadTableCSV(csvPath, rows)) {
    std::cerr << "ERROR: cannot read table CSV: " << csvPath << "\n";
    return;
  }
  std::cout << "Read rows: " << rows.size() << " from " << csvPath << "\n";

  // Group points into curves (plot filter)
  std::map<CurveKey, CurvePoints> curves;
  std::map<int, std::pair<double,double>> ptEdgesByBin;

  for (const auto& r : rows) {
    if (!PassPointDefault(r, cfg)) continue;

    ptEdgesByBin[r.pt_bin] = {r.pt_lo, r.pt_hi};

    CurveKey k;
    k.pt_bin = r.pt_bin;

    if (isGroup) {
      k.z_bin = -1;
      k.setting_id = r.setting_id;

      double zset = ParseZSetting(r.curve_label);
      k.sort_z = zset;

      if (std::isfinite(zset) && zset != kMissing) {
        k.legend = Form("z = %.2f", zset);
      } else {
        k.legend = r.curve_label.empty() ? r.setting_id : r.curve_label;
      }
    } else {
      k.z_bin = r.z_bin;
      k.setting_id = r.setting_id;
      k.sort_z = r.z_center;
      k.legend = Form("zbin%d", r.z_bin);
    }

    auto& cp = curves[k];
    cp.key = k;
    cp.x.push_back(r.phi_center);
    cp.y.push_back(r.xsec);
    cp.ey.push_back(r.xsec_err);

    cp.pt_lo = r.pt_lo; cp.pt_hi = r.pt_hi; cp.pt_center = r.pt_center;

    if (isGroup) {
      double zset = ParseZSetting(r.curve_label);
      cp.z_center = zset;
      cp.z_lo = kMissing; cp.z_hi = kMissing;
    } else {
      cp.z_lo = r.z_lo; cp.z_hi = r.z_hi; cp.z_center = r.z_center;
    }
  }

  if (curves.empty()) {
    std::cerr << "ERROR: after filtering, 0 curves/points remain.\n";
    return;
  }

  // Sort points in each curve by phi
  for (auto& kv : curves) {
    auto& cp = kv.second;
    std::vector<size_t> idx(cp.x.size());
    for (size_t i=0;i<idx.size();i++) idx[i]=i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return cp.x[a] < cp.x[b]; });
    auto reorder = [&](std::vector<double>& v){
      std::vector<double> tmp; tmp.reserve(v.size());
      for (auto i : idx) tmp.push_back(v[i]);
      v.swap(tmp);
    };
    reorder(cp.x); reorder(cp.y); reorder(cp.ey);
  }

  // pT bins present
  std::vector<int> ptBins;
  for (auto& kv : ptEdgesByBin) ptBins.push_back(kv.first);
  std::sort(ptBins.begin(), ptBins.end());

  const int nPt = (int)ptBins.size();
  int nCols = (nPt <= 2) ? nPt : 2;
  int nRows = (int)std::ceil((double)nPt / (double)nCols);

  // Larger canvas for safe margins
  int width  = 900 * nCols;
  int height = 650 * nRows;

  TCanvas* c = new TCanvas("c_xsec_from_table", "xsec vs phipq (from table)", width, height);
  c->Divide(nCols, nRows, 0.002, 0.002);

  std::ofstream fitcsv;
  if (cfg.save_fit_csv) {
    fitcsv.open(fitPath, std::ios::out);
    WriteFitHeader(fitcsv);
    std::cout << "Fit CSV: " << fitPath << "\n";
  }

  static const Color_t palette[] = { kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7, kCyan+2, kViolet+1 };
  static const int nPal = sizeof(palette)/sizeof(palette[0]);

  std::vector<TF1*> fitFuncs; // keep alive only for drawn fits

  auto ptEdge = [&](int pt_bin)->std::pair<double,double>{
    auto it = ptEdgesByBin.find(pt_bin);
    if (it!=ptEdgesByBin.end()) return it->second;
    return {kMissing,kMissing};
  };

  for (int ip=0; ip<nPt; ip++) {
    c->cd(ip+1);

    gPad->SetGrid(1,1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.18);
    gPad->SetTickx();
    gPad->SetTicky();

    int pt_bin = ptBins[ip];
    auto [ptLo, ptHi] = ptEdge(pt_bin);
    double ptCtr = (std::isfinite(ptLo) && std::isfinite(ptHi)) ? 0.5*(ptLo+ptHi) : kMissing;

    // curves in this pad
    std::vector<CurvePoints*> curvesHere;
    for (auto& kv : curves) if (kv.first.pt_bin == pt_bin) curvesHere.push_back(&kv.second);

    // order by z (sort_z)
    std::sort(curvesHere.begin(), curvesHere.end(),
      [&](const CurvePoints* a, const CurvePoints* b){
        double za = a->key.sort_z, zb = b->key.sort_z;
        if (std::isfinite(za) && std::isfinite(zb) && za != zb) return za < zb;
        return a->key.legend < b->key.legend;
      });

    // y-range based on points
    double yMax = 0.0;
    for (auto* cp : curvesHere)
      for (size_t i=0;i<cp->y.size();i++)
        yMax = std::max(yMax, cp->y[i] + cp->ey[i]);
    if (yMax <= 0.0) yMax = 1.0;

    TH1F* frame = (TH1F*)gPad->DrawFrame(0.0, 0.0, TMath::TwoPi(), 1.20*yMax);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(cfg.xTitle);
    frame->GetYaxis()->SetTitle(cfg.yTitle);

    frame->GetXaxis()->SetTitleSize(0.060);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.050);
    frame->GetYaxis()->SetLabelSize(0.050);

    frame->GetXaxis()->SetTitleOffset(1.25);
    frame->GetYaxis()->SetTitleOffset(1.45);
    frame->GetXaxis()->SetLabelOffset(0.010);
    frame->GetYaxis()->SetLabelOffset(0.010);

    TLegend* leg = new TLegend(0.20, 0.68, 0.58, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.042);
    leg->SetHeader(Form("p_{T} [%.2f, %.2f]", ptLo, ptHi), "C");

    int colorIndex = 0;
    int drawn = 0;

    for (auto* cp : curvesHere) {
      int N = (int)cp->x.size();
      if (N <= 0) continue;

      std::vector<double> ex(N, 0.0);
      TGraphErrors* g = new TGraphErrors(N, cp->x.data(), cp->y.data(), ex.data(), cp->ey.data());

      Color_t col = palette[colorIndex % nPal];
      colorIndex++;

      g->SetMarkerStyle(20 + (drawn % 6));
      g->SetMarkerSize(0.95);
      g->SetLineWidth(2);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->Draw("E1P SAME");

      TF1* f = nullptr;
      FitOut fo = FitCurve(*cp, cfg, f);
      fo.mode = isGroup ? "group" : "single";
      fo.group_id = id;
      fo.pt_center = ptCtr;

      if (isGroup) {
        fo.z_bin = -1;
        fo.z_lo = kMissing;
        fo.z_hi = kMissing;
        fo.z_center = cp->z_center;
      }

      // Always write fit row (even invalid), but only draw fit if trustworthy
      if (cfg.save_fit_csv && fitcsv.is_open()) WriteFitRow(fitcsv, fo);

      if (f) {
        if (fo.valid_fit_default == 1) {
          f->SetLineColor(col);
          f->SetLineWidth(2);
          f->Draw("L SAME");
          fitFuncs.push_back(f);
        } else {
          delete f; // don't draw invalid fit
        }
      }

      if (isGroup) {
        leg->AddEntry(g, cp->key.legend.c_str(), "pe");
      } else {
        leg->AddEntry(g, Form("z #in [%.2f, %.2f]", cp->z_lo, cp->z_hi), "pe");
      }

      drawn++;
    }

    leg->Draw();
  }

  if (cfg.save_fit_csv && fitcsv.is_open()) fitcsv.close();

  if (cfg.save_png) {
    c->SaveAs(pngPath.c_str());
    std::cout << "Saved PNG: " << pngPath << "\n";
  }

  std::cout << "Done.\n";
}
