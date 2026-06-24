#ifndef RSIDIS_XS_V4_MANIFESTLITE_H
#define RSIDIS_XS_V4_MANIFESTLITE_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace rsidis_xs_v4 {

struct ManifestInputs {
  std::string treeNameData = "T";
  std::string treeNameSim  = "h10";
  std::vector<std::string> dataFiles;
  std::vector<std::string> simFiles;
  std::string dataGlob;
  std::string simGlob;
};

inline std::vector<std::string> ExtractStringArray(const std::string& text, const std::string& key) {
  std::vector<std::string> out;

  // Use [\s\S] to match across newlines without relying on non-portable dotall flags.
  std::regex re("\"" + key + "\"\\s*:\\s*\\[([\\s\\S]*?)\\]", std::regex_constants::icase);
  std::smatch m;
  if (!std::regex_search(text, m, re)) return out;
  std::string inside = m[1].str();

  std::regex s("\"(.*?)\"");
  for (auto it = std::sregex_iterator(inside.begin(), inside.end(), s);
       it != std::sregex_iterator(); ++it) {
    out.push_back((*it)[1].str());
  }
  return out;
}

inline std::string ExtractString(const std::string& text, const std::string& key) {
  std::regex re("\"" + key + "\"\\s*:\\s*\"(.*?)\"", std::regex_constants::icase);
  std::smatch m;
  if (!std::regex_search(text, m, re)) return "";
  return m[1].str();
}

inline ManifestInputs ReadManifestLite(const std::string& manifestPath) {
  std::ifstream f(manifestPath);
  if (!f) throw std::runtime_error("Cannot open manifest: " + manifestPath);
  std::ostringstream ss;
  ss << f.rdbuf();
  std::string txt = ss.str();

  ManifestInputs in;

  // Tree names (optional)
  std::string tData = ExtractString(txt, "tree_name_data");
  if (tData.empty()) tData = ExtractString(txt, "treeData");
  if (tData.empty()) tData = ExtractString(txt, "tree");
  if (!tData.empty()) in.treeNameData = tData;

  std::string tSim = ExtractString(txt, "tree_name_sim");
  if (tSim.empty()) tSim = ExtractString(txt, "treeSim");
  if (!tSim.empty()) in.treeNameSim = tSim;

  // Files arrays (optional)
  for (auto key : {"data_rootfiles","data_files","dataFiles"}) {
    auto v = ExtractStringArray(txt, key);
    if (!v.empty()) { in.dataFiles = v; break; }
  }
  for (auto key : {"sim_rootfiles","sim_files","simFiles"}) {
    auto v = ExtractStringArray(txt, key);
    if (!v.empty()) { in.simFiles = v; break; }
  }

  // Globs (optional)
  for (auto key : {"data_root_glob","data_glob","dataGlob"}) {
    auto v = ExtractString(txt, key);
    if (!v.empty()) { in.dataGlob = v; break; }
  }
  for (auto key : {"sim_root_glob","sim_glob","simGlob"}) {
    auto v = ExtractString(txt, key);
    if (!v.empty()) { in.simGlob = v; break; }
  }

  return in;
}

} // namespace rsidis_xs_v4

#endif
