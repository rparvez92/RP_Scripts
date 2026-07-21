#ifndef RSIDIS_XS_V4_OUTPUTMANAGER_H
#define RSIDIS_XS_V4_OUTPUTMANAGER_H

#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace rsidis_xs_v4 {

inline bool PathExists(const std::string& path) {
  struct stat st;
  return ::stat(path.c_str(), &st) == 0;
}

inline void MkdirP(const std::string& path) {
  if (path.empty()) return;
  std::string cur;
  for (size_t i = 0; i < path.size(); i++) {
    char c = path[i];
    cur.push_back(c);
    if (c == '/' || i == path.size()-1) {
      if (cur.size() == 1 && cur[0] == '/') continue;
      if (!PathExists(cur)) ::mkdir(cur.c_str(), 0775);
    }
  }
}

inline std::string Dirname(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return ".";
  if (pos == 0) return "/";
  return p.substr(0, pos);
}

inline std::string Basename(const std::string& p) {
  auto pos = p.find_last_of('/');
  if (pos == std::string::npos) return p;
  return p.substr(pos+1);
}

inline std::string NormalizeSlashes(std::string s) {
  std::string out;
  bool prevSlash = false;
  for (char c : s) {
    if (c == '/') {
      if (!prevSlash) out.push_back(c);
      prevSlash = true;
    } else {
      out.push_back(c);
      prevSlash = false;
    }
  }
  if (out.size() > 1 && out.back() == '/') out.pop_back();
  return out;
}

// settingsRoot can be "settings" OR an absolute path ".../settings"
inline std::string MakeSettingIdFromManifestPath(const std::string& manifestPath,
                                                const std::string& settingsRoot = "settings") {
  std::string mdir  = NormalizeSlashes(Dirname(manifestPath));
  std::string sroot = NormalizeSlashes(settingsRoot);

  // Case 1: absolute prefix stripping: <settingsRoot>/<setting_id>
  if (mdir.rfind(sroot + "/", 0) == 0) {
    return mdir.substr(sroot.size() + 1);
  }

  // Case 2: if user passed absolute path but we can match by basename
  std::string sbase = Basename(sroot); // likely "settings"
  auto pos = mdir.find("/" + sbase + "/");
  if (pos != std::string::npos) {
    return mdir.substr(pos + sbase.size() + 2);
  }
  if (mdir.rfind(sbase + "/", 0) == 0) {
    return mdir.substr(sbase.size() + 1);
  }

  // Fallback: leaf folder name
  return Basename(mdir);
}

inline std::string MakeResultsDir(const std::string& resultsRoot,
                                  const std::string& settingId) {
  std::string rr  = NormalizeSlashes(resultsRoot);
  std::string sid = NormalizeSlashes(settingId);
  return NormalizeSlashes(rr + "/" + sid);
}

inline void EnsureSettingOutputDirs(const std::string& resultsDir) {
  MkdirP(resultsDir);
  MkdirP(resultsDir + "/PNGs");
  MkdirP(resultsDir + "/tables");
}

} // namespace rsidis_xs_v4

#endif
