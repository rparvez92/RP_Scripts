#ifndef RSIDIS_XS_V4_CSVWRITER_H
#define RSIDIS_XS_V4_CSVWRITER_H

#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>

namespace rsidis_xs_v4 {

class CSVWriter {
public:
  explicit CSVWriter(const std::string& path) : m_out(path) {
    if (!m_out) throw std::runtime_error("Cannot open CSV for write: " + path);
  }

  void WriteHeader(const std::vector<std::string>& cols) {
    for (size_t i=0;i<cols.size();i++) {
      if (i) m_out << ",";
      m_out << cols[i];
    }
    m_out << "\n";
  }

  template <typename... Args>
  void WriteRow(const Args&... args) {
    WriteRowImpl(0, args...);
    m_out << "\n";
  }

private:
  template <typename T, typename... Rest>
  void WriteRowImpl(int idx, const T& v, const Rest&... rest) {
    if (idx) m_out << ",";
    m_out << v;
    WriteRowImpl(idx+1, rest...);
  }
  void WriteRowImpl(int) {}

  std::ofstream m_out;
};

} // namespace rsidis_xs_v4

#endif
