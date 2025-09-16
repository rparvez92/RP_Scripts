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
        {"hsdelta", "H.gtr.dp"},
        {"hsytar",  "H.gtr.y"},
        {"hsxptar",  "H.gtr.th"},
        {"hsyptar",  "H.gtr.ph"},
        {"xb", "H.kin.x_bj"}, //        {"xb", "H.kin.primary.x_bj"},
        {"q2", "H.kin.Q2"}, //        {"q2", "H.kin.primary.Q2"},
        {"w", "H.kin.W"} //        {"w", "H.kin.primary.W"},
    }; //stdMap = Sim To Data Map

    auto it = stdMap.find(simVar); //it=iterator
    if (it != stdMap.end()) {
        return it->second; //gives you the second element of the map
    } else {
        throw std::invalid_argument("No linked variable for simulation variable: " + simVar);
    }
}

#endif // end of include guard
