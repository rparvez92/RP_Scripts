#ifndef RSIDIS_XS_V4_BINNING_H
#define RSIDIS_XS_V4_BINNING_H

#include <vector>

namespace rsidis_xs_v4 {

struct Binning {
  std::vector<double> ptEdges;
  std::vector<double> zEdges;
  std::vector<double> phiEdges; // 0..2pi
};

inline double Center(double low, double high) { return 0.5*(low+high); }

} // namespace rsidis_xs_v4

#endif
