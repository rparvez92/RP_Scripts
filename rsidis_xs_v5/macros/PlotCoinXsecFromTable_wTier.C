// PlotCoinXsecFromTable_wTier.C (rsidis_xs_v5)
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
//   root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("settings/.../manifest.txt","results","settings")'
//   root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("groups/.../grp_*.list","results","settings")'
//   root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/groups/pass4/pi+sidis/LH2/x0p25/q23p3/thpq2/grp_pass4_piplus_LH2_zOv_x0p25_q23p3_thpq2.list", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/results", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/settings", false, true, true)'
//   root -l -b -q 'macros/PlotCoinXsecFromTable_wTier.C("/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/groups/pass4/pi+sidis/LH2/x0p25/q23p3/thpq2/grp_pass4_piplus_LH2_zOv_x0p25_q23p3_thpq2.list", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/results", "/home/cdaq/users/rparvez/RP_Scripts/rsidis_xs_v5/settings", true, true, true)'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <tuple>
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


  // Fit conditioning (Phase-A stability)
  // Dynamic upper bound for M0: M0_max = min(M0_max_cap, M0_max_scale * max(y))
  double M0_max_scale = 20.0;
  double M0_max_cap   = 1e9;

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

  // diagnostic: draw invalid-fit curves as dashed lines
  // - tierOn mode: only drawn if this is true
  // - tierOff mode: fit curve is always drawn (even if invalid), but dashed if invalid
  bool draw_invalid_fits = true;
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
  std::string attempts;
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

// Like PassPointDefault, but returns a drop-reason token for diagnostics/auditing.
static bool PassPointDefaultReason(const Row& r, const PlotConfig& cfg, std::string& reason) {
  reason.clear();

  if (!std::isfinite(r.phi_center)) { reason="NONFINITE_PHI"; return false; }
  if (r.phi_center == kMissing) { reason="MISSING_PHI"; return false; }

  if (!std::isfinite(r.xsec) || r.xsec == kMissing) { reason="NONFINITE_XSEC"; return false; }
  if (!std::isfinite(r.xsec_err) || r.xsec_err == kMissing) { reason="NONFINITE_XSEC_ERR"; return false; }
  if (r.xsec_err <= 0.0) { reason="NONPOS_XSEC_ERR"; return false; }

  if (cfg.use_valid_default && r.valid_default != 1) { reason="VALID_DEFAULT0"; return false; }
  if (!cfg.allow_negative_xsec && r.xsec <= 0.0) { reason="NONPOS_XSEC"; return false; }

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
  std::vector<long long> row_i; // index into full rows vector for each plotted point
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


// build fit vectors from plotted points with fit-only filtering and error floor.
// Also returns per-plot-point fit exclusion reasons and the list of plot indices
// actually used in the final fit (after robust drop).
static void BuildFitVectors(const CurvePoints& cp, const PlotConfig& cfg,
                            std::vector<double>& xf, std::vector<double>& yf, std::vector<double>& ef,
                            std::vector<int>& usedPlotIdx,
                            std::vector<std::string>& fitExclReasonPlot,
                            bool& didDrop, int& droppedPlotIndex)
{
  xf.clear(); yf.clear(); ef.clear();
  usedPlotIdx.clear();
  didDrop = false;
  droppedPlotIndex = -1;

  fitExclReasonPlot.assign(cp.x.size(), "");

  // fit-only relerr cut
  std::vector<int> candPlotIdx;
  candPlotIdx.reserve(cp.x.size());

  for (size_t i=0;i<cp.x.size();i++) {
    double x = cp.x[i];
    double y = cp.y[i];
    double e = cp.ey[i];

    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(e) || y<=0.0 || e<=0.0) {
      fitExclReasonPlot[i] = "FIT_BADVAL";
      continue;
    }

    const double rel = e / y;
    if (std::isfinite(rel) && rel > cfg.fit_relerr_max) {
      fitExclReasonPlot[i] = "FIT_RELERR_GT_MAX";
      continue;
    }

    candPlotIdx.push_back((int)i);
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
    didDrop = true;
    droppedPlotIndex = candPlotIdx[imax];
    fitExclReasonPlot[droppedPlotIndex] = "FIT_ROBUST_DROPPED";

    xf.erase(xf.begin()+imax);
    yf.erase(yf.begin()+imax);
    ef.erase(ef.begin()+imax);
    candPlotIdx.erase(candPlotIdx.begin()+imax);
  }

  // remaining candidates are used in fit
  usedPlotIdx = candPlotIdx;
  for (int ip : usedPlotIdx) {
    if (fitExclReasonPlot[ip].empty()) fitExclReasonPlot[ip] = "FIT_INCLUDED";
  }
}



// ---------- fit ----------
static std::string AttemptReasonToken(const FitOut& fo) {
  // Priority order roughly matches what you care about in diagnostics
  if (fo.fit_flag_bits & NO_FIT_POINTS) return "NO_PTS";
  if (fo.fit_flag_bits & (COV_NOT_POSDEF | COV_FAILED)) return "COV";
  if (fo.fit_flag_bits & FIT_STATUS_BAD) return "STATUS";
  if (fo.fit_flag_bits & NDF_NONPOSITIVE) return "NDF";
  if (fo.fit_flag_bits & PHI_UNIFORM_FAIL) return "PHI";
  if (fo.fit_flag_bits & DOMINANT_WEIGHT_PT) return "DOM";
  if (fo.fit_flag_bits & (NONFINITE_PARAM | NONFINITE_ERR)) return "NONFIN";
  if (fo.fit_flag_bits & CHI2NDF_HUGE) return "CHI2";
  if (fo.fit_flag_bits & PROB_TINY) return "PROB";
  return "BAD";
}

static FitOut FitCurve(const CurvePoints& cp, const PlotConfig& cfg, bool tierOn,
                       std::vector<int>& usedPlotIdx,
                       std::vector<std::string>& fitExclReasonPlot,
                       TF1*& outFunc)
{
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

  usedPlotIdx.clear();
  fitExclReasonPlot.clear();

  // Build fit-only vectors (and exclusions)
  std::vector<double> xf, yf, ef;
  bool didDrop=false;
  int droppedPlotIndex=-1;
  BuildFitVectors(cp, cfg, xf, yf, ef, usedPlotIdx, fitExclReasonPlot, didDrop, droppedPlotIndex);
  fo.npts_fit = (int)xf.size();

  // span/occ for fit points
  if (fo.npts_fit >= 1) {
    double xmin=1e99, xmax=-1e99;
    for (double x : xf) { xmin = min(xmin, x); xmax = max(xmax, x); }
    fo.phi_min = xmin;
    fo.phi_max = xmax;
    fo.phi_span = (xmax > xmin) ? (xmax - xmin) : 0.0;
    fo.phi_occ = PhiOccupancy(xf, cfg.phi_nsectors);
    fo.dominant_wfrac = DominantWeightFrac(ef);
  }

  if (didDrop) {
    fo.fit_flag_bits |= ROBUST_DROPPED_PT;
  }

  // For tierOff reference: always try tier-2 model (3 params) if possible
  const int forcedTier = 2;

  const int minPointsNeeded = tierOn ? 2 : 3; // tierOff needs >=3 points to constrain 3 params
  if (fo.npts_fit < minPointsNeeded) {
    fo.fit_flag_bits |= NO_FIT_POINTS;
    fo.fit_status = 999;
    fo.cov_status = -1;
    fo.ndf = -999;
    fo.fit_tier = tierOn ? ChooseTier(fo.npts_fit, fo.phi_span, fo.phi_occ, cfg) : forcedTier;
    ComputeFitValidity(fo, cfg);

    // attempt history
    if (!tierOn) fo.attempts = "t2:bad(NO_PTS)";
    return fo;
  }

  // Decide starting tier from fit points (tierOn) or force tier2 (tierOff).
  const int tierStart = tierOn ? ChooseTier(fo.npts_fit, fo.phi_span, fo.phi_occ, cfg) : forcedTier;
  const int minTier   = tierOn ? 0 : forcedTier;

  // compute mean/range for initial guesses
  double ysum=0.0, yMin=1e99, yMax=-1e99;
  double sumw=0.0, sumwy=0.0;
  for (size_t i=0;i<yf.size();i++) {
    const double y = yf[i];
    ysum += y;
    yMin = min(yMin, y);
    yMax = max(yMax, y);
    const double e = ef[i];
    if (std::isfinite(e) && e > 0) {
      const double w = 1.0/(e*e);
      sumw  += w;
      sumwy += w*y;
    }
  }
  const double ymean  = ysum / (double)fo.npts_fit;
  const double ywmean = (sumw > 0.0) ? (sumwy / sumw) : ymean;
  const double ywerr  = (sumw > 0.0) ? (1.0 / std::sqrt(sumw)) : kMissing;

  // Dynamic M0 upper bound to improve Minuit conditioning
  double baseY = yMax;
  if (!std::isfinite(baseY) || baseY <= 0.0)
    baseY = (std::isfinite(ywmean) && ywmean > 0.0) ? ywmean : 1e-6;
  double M0max = cfg.M0_max_scale * baseY;
  if (std::isfinite(cfg.M0_max_cap) && cfg.M0_max_cap > 0.0)
    M0max = std::min(M0max, cfg.M0_max_cap);
  if (!std::isfinite(M0max) || M0max <= cfg.M0_min)
    M0max = cfg.M0_min + std::max(1e-6, 10.0*baseY);

  // Build graph for fit once (data are the same for all tiers)
  std::vector<double> ex(fo.npts_fit, 0.0);
  TGraphErrors gfit(fo.npts_fit, xf.data(), yf.data(), ex.data(), ef.data());

  auto MakeFunc = [&](int tier)->TF1*{
    TF1* f = nullptr;
    if (tier <= 0) {
      f = new TF1(Form("f0_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()), "[0]", fo.phi_min, fo.phi_max);
      f->SetParameters(ywmean);
      f->SetParLimits(0, cfg.M0_min, M0max);
      return f;
    }

    // initial A from peak-to-valley heuristic (bounded)
    double A0 = 0.0;
    if (yMax + yMin > 0) {
      A0 = (yMax - yMin) / (yMax + yMin);
      A0 = max(-0.9, min(0.9, A0));
    }

    if (tier == 1) {
      f = new TF1(Form("f1_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()),
                  "[0]*(1 + [1]*cos(x))", fo.phi_min, fo.phi_max);
      f->SetParameters(ywmean, A0);
      f->SetParLimits(0, cfg.M0_min, M0max);
      f->SetParLimits(1, -1.0, 1.0);
      return f;
    }

    f = new TF1(Form("f2_pt%d_%s", fo.pt_bin, fo.curve_label.c_str()),
                "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", fo.phi_min, fo.phi_max);
    f->SetParameters(ywmean, A0, 0.0);
    f->SetParLimits(0, cfg.M0_min, M0max);
    f->SetParLimits(1, -1.0, 1.0);
    f->SetParLimits(2, -1.0, 1.0);
    return f;
  };

  auto FillFromFunc = [&](FitOut& o, TF1* f, const TFitResultPtr& r){
    o.fit_status = (int)r;
    if (r.Get()) {
      o.cov_status = r->CovMatrixStatus();
      o.chi2 = r->Chi2();
      o.ndf  = r->Ndf();
      o.prob = r->Prob();
    } else {
      o.cov_status = -1;
      o.chi2 = f->GetChisquare();
      o.ndf  = f->GetNDF();
      o.prob = f->GetProb();
    }
    if (o.ndf > 0) o.chi2_ndf = o.chi2 / (double)o.ndf;

    o.M0 = f->GetParameter(0);
    o.M0_err = f->GetParError(0);

    // Phase-A stability: if tier-0 covariance is unreliable, fall back to analytic
    // weighted-mean uncertainty (computed from fit vectors after sigma-floor).
    if (o.fit_tier <= 0 && std::isfinite(ywerr) && ywerr > 0) {
      const bool covBad = (o.cov_status < 2);
      const bool errBad = (!std::isfinite(o.M0_err) || o.M0_err <= 0.0 ||
                           (std::isfinite(o.M0) && o.M0 > 0.0 && (o.M0_err/o.M0) > cfg.M0_relerr_max));
      if (covBad || errBad) {
        o.M0_err = ywerr;
        AppendNote(o.note, "M0ERR_WMEAN");
      }
    }
    if (o.fit_tier >= 1) { o.A = f->GetParameter(1); o.A_err = f->GetParError(1); }
    if (o.fit_tier >= 2) { o.B = f->GetParameter(2); o.B_err = f->GetParError(2); }

    if (std::isfinite(o.A) && std::isfinite(o.A_err) && o.A_err > 0) o.A_sig = abs(o.A)/o.A_err;
    if (std::isfinite(o.B) && std::isfinite(o.B_err) && o.B_err > 0) o.B_sig = abs(o.B)/o.B_err;

    auto sane = [&](double v){ return std::isfinite(v) ? v : kMissing; };
    o.M0 = sane(o.M0); o.M0_err = sane(o.M0_err);
    o.A  = sane(o.A);  o.A_err  = sane(o.A_err);
    o.B  = sane(o.B);  o.B_err  = sane(o.B_err);
    o.chi2 = sane(o.chi2);
    o.prob = sane(o.prob);
    o.chi2_ndf = sane(o.chi2_ndf);
    o.A_sig = sane(o.A_sig);
    o.B_sig = sane(o.B_sig);
  };

  FitOut best = fo;
  best.fit_tier = tierStart;
  best.note.clear();
  best.attempts.clear();
  if (didDrop) best.fit_flag_bits |= ROBUST_DROPPED_PT;

  TF1* bestFunc = nullptr;
  bool foundValid = false;

  std::map<int, FitOut> attempted; // tier -> result

  for (int tier = tierStart; tier >= minTier; --tier) {
    FitOut cand = fo;
    cand.fit_tier = tier;
    cand.note.clear();
    cand.attempts.clear();
    if (didDrop) cand.fit_flag_bits |= ROBUST_DROPPED_PT;

    TF1* f = MakeFunc(tier);
    TFitResultPtr r = gfit.Fit(f, "RQS0");
    FillFromFunc(cand, f, r);
    ComputeFitValidity(cand, cfg);
    attempted[tier] = cand;

    if (!foundValid && cand.valid_fit_default == 1) {
      if (tierOn && tier != tierStart) {
        AppendNote(cand.note, Form("FALLBACK_FROM_TIER%d", tierStart));
      }
      best = cand;
      bestFunc = f;
      foundValid = true;
      break;
    }

    if (!foundValid) {
      // keep best invalid backup: prefer fewer reject bits, then simpler tier, then higher prob
      const unsigned int catastrophic = FIT_STATUS_BAD | NONFINITE_PARAM | NONFINITE_ERR | NDF_NONPOSITIVE | COV_FAILED | NO_FIT_POINTS;
      const unsigned int defaultReject = catastrophic | COV_NOT_POSDEF | CHI2NDF_HUGE | PROB_TINY | PHI_UNIFORM_FAIL | DOMINANT_WEIGHT_PT;

      auto score = [&](const FitOut& o){
        const unsigned int rej = (o.fit_flag_bits & defaultReject);
        const int bitcnt = __builtin_popcount((unsigned int)rej);
        const double p = (std::isfinite(o.prob) ? o.prob : -1.0);
        return std::make_tuple(bitcnt, o.fit_tier, -p); // fewer bits, simpler tier, higher prob
      };

      if (!bestFunc) {
        best = cand;
        bestFunc = f;
      } else if (score(cand) < score(best)) {
        best = cand;
        delete bestFunc;
        bestFunc = f;
      } else {
        delete f;
      }
    } else {
      delete f;
    }
  }

  outFunc = bestFunc;

  // Build attempt history string, e.g. "t2:bad(COV); t1:bad(COV); t0:ok_or_best"
  {
    std::string s;
    auto add = [&](int tier, const std::string& tok){
      if (!s.empty()) s += "; ";
      s += Form("t%d:%s", tier, tok.c_str());
    };

    // prefer order t2 -> t1 -> t0 when present
    for (int t=2; t>=0; --t) {
      auto it = attempted.find(t);
      if (it == attempted.end()) continue;
      if (t == best.fit_tier) {
        add(t, "ok_or_best");
      } else {
        const FitOut& a = it->second;
        if (a.valid_fit_default == 1) add(t, "ok");
        else add(t, std::string("bad(") + AttemptReasonToken(a) + ")");
      }
    }
    best.attempts = s;
  }

  return best;
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
    << "attempts,note\n";
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
    << fo.attempts << ","
    << fo.note
    << "\n";
}

// ---------- parsed points CSV (all rows) ----------
static void WriteParsedHeader(std::ofstream& out) {
  out
    << "mode,group_id,tier_tag,"
    << "curve_label_raw,plot_label,setting_id,"
    << "pt_bin,z_bin,phi_bin,"
    << "pt_lo,pt_hi,pt_center,"
    << "z_lo,z_hi,z_center,"
    << "phi_center,"
    << "xsec,xsec_err,rel_xsec_err,n_sim,flag_bits,valid_default,"
    << "is_plotted,plot_drop_reason,"
    << "is_used_in_fit,fit_excl_reason,"
    << "fit_tier,valid_fit_default,"
    << "y_fit,residual,pull,"
    << "attempts,note\n";
}

static double EvalSigmaModel(double phi, int tier, double M0, double A, double B) {
  if (!std::isfinite(phi) || phi==kMissing) return kMissing;
  if (!std::isfinite(M0) || M0==kMissing) return kMissing;
  const double c1 = std::cos(phi);
  const double c2 = std::cos(2.0*phi);
  if (tier <= 0) return M0;
  if (tier == 1) {
    if (!std::isfinite(A) || A==kMissing) return kMissing;
    return M0*(1.0 + A*c1);
  }
  if (!std::isfinite(A) || A==kMissing) return kMissing;
  if (!std::isfinite(B) || B==kMissing) return kMissing;
  return M0*(1.0 + A*c1 + B*c2);
}

static void WriteParsedRow(std::ofstream& out,
                           const Row& r,
                           long long rowIndex,
                           const std::string& modeTag,
                           const std::string& groupId,
                           const std::string& tierTag,
                           const std::string& plotLabel,
                           bool isPlotted,
                           const std::string& plotDropReason,
                           bool isUsedInFit,
                           const std::string& fitExclReason,
                           const FitOut* fo)
{
  // Fit-based diagnostics (optional)
  int fitTier = -1;
  int validFit = 0;
  double yFit = kMissing;
  double resid = kMissing;
  double pull = kMissing;
  std::string attempts;
  std::string note;

  if (fo) {
    fitTier = fo->fit_tier;
    validFit = fo->valid_fit_default;
    attempts = fo->attempts;
    note = fo->note;

    yFit = EvalSigmaModel(r.phi_center, fo->fit_tier, fo->M0, fo->A, fo->B);
    if (std::isfinite(yFit) && yFit!=kMissing && std::isfinite(r.xsec) && r.xsec!=kMissing) {
      resid = r.xsec - yFit;
      if (std::isfinite(r.xsec_err) && r.xsec_err>0.0) pull = resid / r.xsec_err;
    }
  }

  out
    << modeTag << ","
    << groupId << ","
    << tierTag << ","
    << r.curve_label << ","
    << plotLabel << ","
    << r.setting_id << ","
    << r.pt_bin << ","
    << r.z_bin << ","
    << r.phi_bin << ","
    << r.pt_lo << "," << r.pt_hi << "," << r.pt_center << ","
    << r.z_lo << "," << r.z_hi << "," << r.z_center << ","
    << r.phi_center << ","
    << r.xsec << "," << r.xsec_err << "," << r.rel_xsec_err << "," << r.n_sim << ","
    << r.flag_bits << "," << r.valid_default << ","
    << (isPlotted ? 1 : 0) << ","
    << plotDropReason << ","
    << (isUsedInFit ? 1 : 0) << ","
    << fitExclReason << ","
    << fitTier << ","
    << validFit << ","
    << yFit << ","
    << resid << ","
    << pull << ","
    << attempts << ","
    << note
    << "\n";
}


} // namespace

// ------------------------------------------------------------
// Main macro entry
// ------------------------------------------------------------
void PlotCoinXsecFromTable_wTier(const char* manifestOrGroupPath,
                                 const char* resultsRoot = "results",
                                 const char* settingsRoot = "settings",
                                 bool tierOn = true,
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

  const std::string tierTag = tierOn ? "tierOn" : "tierOff";

  if (isGroup) {
    id = StripExtension(Basename(p)); // group_id
    resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + id);
    csvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_group.csv";
    fitPath = resultsDir + "/tables/fit_parameters_group_" + tierTag + ".csv";
    pngPath = resultsDir + "/PNGs/xsec_phipq_z_pt_overlayed_group_" + tierTag + ".png";
  } else {
    id = MakeSettingIdFromManifestPath(p, NormalizeSlashes(settingsRoot));
    resultsDir = NormalizeSlashes(std::string(resultsRoot) + "/" + id);
    csvPath = resultsDir + "/tables/xsec_phipq_z_pt_overlayed_single.csv";
    fitPath = resultsDir + "/tables/fit_parameters_single_" + tierTag + ".csv";
    pngPath = resultsDir + "/PNGs/xsec_phipq_z_pt_overlayed_single_" + tierTag + ".png";
  }

  const std::string parsedPath = resultsDir + "/tables/parsed_dataPoints_" + std::string(isGroup ? "group" : "single") + "_" + tierTag + ".csv";

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

  for (size_t ir=0; ir<rows.size(); ++ir) {
    const auto& r = rows[ir];

    // Keep pT edges even for rows that won't be plotted (useful for auditing)
    if (r.pt_bin >= 0 && std::isfinite(r.pt_lo) && std::isfinite(r.pt_hi))
      ptEdgesByBin[r.pt_bin] = {r.pt_lo, r.pt_hi};

    // Plot filter (default)
    if (!PassPointDefault(r, cfg)) continue;

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
    cp.row_i.push_back((long long)ir);

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
    {
      std::vector<long long> tmp; tmp.reserve(cp.row_i.size());
      for (auto i : idx) tmp.push_back(cp.row_i[i]);
      cp.row_i.swap(tmp);
    }
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
    fitcsv.open(fitPath, std::ios::out | std::ios::trunc);
    WriteFitHeader(fitcsv);
    fitcsv.flush();
    std::cout << "Fit CSV: " << fitPath << "\n";
  }

  std::ofstream parsedcsv;
  if (cfg.save_fit_csv) {
    parsedcsv.open(parsedPath, std::ios::out | std::ios::trunc);
    WriteParsedHeader(parsedcsv);
    parsedcsv.flush();
    std::cout << "Parsed points CSV: " << parsedPath << "\n";
  }

  // Cache fit outputs and per-point fit-usage reasons for the parsed-data export
  std::map<CurveKey, FitOut> fitByCurve;
  std::map<CurveKey, std::unordered_map<long long, std::string>> fitReasonByCurve;

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
      std::vector<int> usedPlotIdx;
      std::vector<std::string> fitExclReasonPlot;
      FitOut fo = FitCurve(*cp, cfg, tierOn, usedPlotIdx, fitExclReasonPlot, f);
      fo.mode = isGroup ? "group" : "single";
      fo.group_id = id;
      fo.pt_center = ptCtr;

      if (isGroup) {
        fo.z_bin = -1;
        fo.z_lo = kMissing;
        fo.z_hi = kMissing;
        fo.z_center = cp->z_center;
      }

            // Cache fit outputs for parsed points export
      fitByCurve[cp->key] = fo;
      {
        auto& mreason = fitReasonByCurve[cp->key];
        // Map per-plot-point fit inclusion/exclusion token back to original table row index
        for (size_t ipt=0; ipt<cp->row_i.size() && ipt<fitExclReasonPlot.size(); ++ipt) {
          mreason[cp->row_i[ipt]] = fitExclReasonPlot[ipt];
        }
      }

// Always write fit row (even invalid), but only draw fit if trustworthy
      if (cfg.save_fit_csv && fitcsv.is_open()) WriteFitRow(fitcsv, fo);

      if (f) {
        // Drawing policy:
        // - tierOn: draw only valid fits by default; optionally draw invalid as dashed (cfg.draw_invalid_fits)
        // - tierOff: always draw the fit curve (even if invalid) so it's a true reference case
        const bool isValid = (fo.valid_fit_default == 1);
        const bool allowDraw = (!tierOn) || isValid || cfg.draw_invalid_fits;

        if (allowDraw) {
          f->SetLineColor(col);
          f->SetLineWidth(2);
          f->SetLineStyle(isValid ? 1 : 2); // solid=valid, dashed=invalid
          f->Draw("L SAME");
          fitFuncs.push_back(f);
        } else {
          delete f;
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


  // Write parsed (all) data points CSV: includes plotted and non-plotted rows
  if (cfg.save_fit_csv && parsedcsv.is_open()) {
    const std::string modeTag = isGroup ? "group" : "single";

    auto MakeKeyAndLabel = [&](const Row& r)->CurveKey {
      CurveKey k;
      k.pt_bin = r.pt_bin;
      if (isGroup) {
        k.z_bin = -1;
        k.setting_id = r.setting_id;
        double zset = ParseZSetting(r.curve_label);
        k.sort_z = zset;
        if (std::isfinite(zset) && zset != kMissing) k.legend = Form("z = %.2f", zset);
        else k.legend = r.curve_label.empty() ? r.setting_id : r.curve_label;
      } else {
        k.z_bin = r.z_bin;
        k.setting_id = r.setting_id;
        k.sort_z = r.z_center;
        k.legend = Form("zbin%d", r.z_bin);
      }
      return k;
    };

    for (long long ir=0; ir<(long long)rows.size(); ++ir) {
      const Row& r = rows[(size_t)ir];

      std::string plotDropReason;
      const bool isPlotted = PassPointDefaultReason(r, cfg, plotDropReason);

      CurveKey k = MakeKeyAndLabel(r);
      const std::string plotLabel = k.legend;

      // Fit lookup (may not exist if curve had zero plotted points)
      const FitOut* foPtr = nullptr;
      auto itf = fitByCurve.find(k);
      if (itf != fitByCurve.end()) foPtr = &itf->second;

      // Fit inclusion/exclusion token
      bool isUsedInFit = false;
      std::string fitExclReason;
      if (isPlotted) {
        auto itm = fitReasonByCurve.find(k);
        if (itm != fitReasonByCurve.end()) {
          auto itp = itm->second.find(ir);
          if (itp != itm->second.end()) {
            fitExclReason = itp->second;
            isUsedInFit = (fitExclReason == "FIT_INCLUDED");
            if (!isUsedInFit && fitExclReason.empty()) fitExclReason = "FIT_EXCLUDED_UNKNOWN";
          } else {
            fitExclReason = "FIT_EXCLUDED_UNKNOWN";
          }
        } else {
          fitExclReason = "FIT_EXCLUDED_NO_CURVE";
        }
      }

      WriteParsedRow(parsedcsv, r, ir, modeTag, id, tierTag, plotLabel,
                     isPlotted, plotDropReason,
                     isUsedInFit, fitExclReason,
                     foPtr);
    }
    parsedcsv.close();
  }

  if (cfg.save_fit_csv && fitcsv.is_open()) fitcsv.close();

  if (cfg.save_png) {
    c->SaveAs(pngPath.c_str());
    std::cout << "Saved PNG: " << pngPath << "\n";
  }

  std::cout << "Done.\n";
}
