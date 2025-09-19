// SimToDataMap.h
// Header-only implementation of SimToDataMap() function for mapping simulation variable to data variable
#ifndef SIM_TO_DATA_MAP_H
#define SIM_TO_DATA_MAP_H

#include <string>
#include <map>
#include <stdexcept>

// link Simulation variable name To corresponding Data-Dummy variable name
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
        {"z", "P.gtr.p/H.kin.primary.nu"},//z needs to be created from other branches
        {"xbj", "H.kin.primary.x_bj"},
        {"Q2", "H.kin.primary.Q2"},
        {"W", "H.kin.primary.W"},
        {"epsilon", "H.kin.primary.epsilon"},
        {"nu", "H.kin.primary.nu"},
        {"thetapq", "P.kin.secondary.th_xq"},
        {"phipq", "P.kin.secondary.ph_xq"}
    }; //stdMap = Sim To Data Map

    auto it1 = stdMap.find(simVar); //it1=iterator1
    if (it1 != stdMap.end()) {
        return it1->second; //gives you the second element of the map
    } else {
        throw std::invalid_argument("No linked variable for simulation variable: " + simVar);
    }
}

// link simulation Branch name to corresponding Physics variable name
inline std::string BranchToPhysicsMap(const std::string& simVar) {
    static const std::map<std::string, std::string> btpMap = {
        // HMS Variables
        {"hsdelta", "HMS Delta"},
        {"hsytar",  "HMS y at Target"},
        {"hsxptar",  "HMS dx/dz"},
        {"hsyptar",  "HMS dy/dz"},
        {"ssdelta", "SHMS Delta"},
        {"ssytar",  "SHMS y at Target"},
        {"ssxptar",  "SHMS dx/dz"},
        {"ssyptar",  "SHMS dy/dz"},
        {"z",  "z"},
        {"xbj",  "Bjorken x"},
        {"Q2",  "Q2"},
        {"W",  "W"},
        {"nu",  "nu"},
        {"epsilon",  "epsilon"},
        {"thetapq",  "thetapq"},
        {"phipq",  "phipq"}
    }; //btpMap = Branch-name To Physics-variable-name Map
    auto it2 = btpMap.find(simVar); //it2=iterator2
    if (it2 != btpMap.end()) {
        return it2->second; //gives you the second element of the map
    } else {
        throw std::invalid_argument("No linked physics variable for branch name: " + simVar);
    }
}

#endif // end of include guard
