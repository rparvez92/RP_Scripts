// SimToDataMap.h
// Header-only implementation of SimToDataMap() function for mapping simulation variable to data variable
#ifndef SIM_TO_DATA_MAP_H
#define SIM_TO_DATA_MAP_H

#include <string>
#include <map>
#include <stdexcept>

// Link simulation variable name to corresponding data/dummy variable name
inline std::string SimToDataMap(const std::string& simVar) {
    static const std::map<std::string, std::string> stdMap = {
        // HMS Variables
        {"hsdelta", "H.gtr.dp"},
        {"hsytar",  "H.gtr.y"},
        {"hsxptar",  "H.gtr.th"},
        {"hsyptar",  "H.gtr.ph"},
        // SHMS Variables
	{"ssdelta", "P.gtr.dp"},
        {"ssytar",  "P.gtr.y"},
        {"ssxptar",  "P.gtr.th"},
        {"ssyptar",  "P.gtr.ph"},
        // Kinematic Variables
        //{"z", "H.kin.primary.W"},
        {"xbj", "H.kin.primary.x_bj"},
        {"Q2", "H.kin.primary.Q2"},
        {"W", "H.kin.primary.W"},
        {"epsilon", "H.kin.primary.epsilon"},
        {"nu", "H.kin.primary.nu"},
        {"phipq", "P.kin.secondary.ph_xq"},
        {"thetapq", "P.kin.secondary.th_xq"}
    }; //stdMap = Sim To Data Map

    auto it = stdMap.find(simVar); //it=iterator
    if (it != stdMap.end()) {
        return it->second; //gives you the second element of the map
    } else {
        throw std::invalid_argument("No linked variable for simulation variable: " + simVar);
    }
}

#endif // end of include guard
