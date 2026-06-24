// Mapping.h
// SimToDataMap() and BranchToPhysicsMap() function for mapping

#ifndef MAPPING_H
#define MAPPING_H

#include <string>
#include <map>
#include <stdexcept>

// link Simulation variable name To corresponding Data-Dummy variable name
inline std::string SimToDataMap(const std::string& simVar) {
    static const std::map<std::string, std::string> stdMap = {
        // HMS Variables
        {"hsdelta", "H_gtr_dp"},
        {"hsytar",  "H_gtr_y"},
        {"hsxptar",  "H_gtr_th"},
        {"hsyptar",  "H_gtr_ph"},
        // SHMS Variables
	{"ssdelta", "P_gtr_dp"},
        {"ssytar",  "P_gtr_y"},
        {"ssxptar",  "P_gtr_th"},
        {"ssyptar",  "P_gtr_ph"},
        // Kinematic Variables
        {"z", "z"},//z needs to be created from other branches
        {"xbj", "H_kin_primary_x_bj"},
        {"Q2", "H_kin_primary_Q2"},
        {"W", "H_kin_primary_W"},
        {"epsilon", "H_kin_primary_epsilon"},
        {"nu", "H_kin_primary_nu"},
        {"thetapq", "P_kin_secondary_th_xq"},
        //{"phipq", "P_kin_secondary_ph_xq"}
	{"phipq", "(P_kin_secondary_ph_xq < 0 ? P_kin_secondary_ph_xq + 2*TMath::Pi() : P_kin_secondary_ph_xq)"}//This ensures phipq of both simc and hcana will be mapped to range [0, 2pi].This line maps hcana variable to [0, 2pi]


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
        {"z", "z"},
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
