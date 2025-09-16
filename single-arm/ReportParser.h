#ifndef REPORT_PARSER_H
#define REPORT_PARSER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

// Create a Structure to hold the values
struct ReportValues {
    double charge_mC = 0.0;
    int ps_factor = 0;
    double hms_eff = 0.0;
};

// Define the function
inline ReportValues ParseReportFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open report file: " + filepath);
    }

    // Declare values and line to hold values and line from file
    ReportValues values;
    std::string line;

    // Parse
    while (std::getline(file, line)) {
        std::istringstream iss(line);

        if (line.find("BCM4C Beam Cut Charge") != std::string::npos) {
            std::string label, unit;
            double charge;
            iss >> label >> label >> label >> label >> charge >> unit;
            values.charge_mC = charge / 1000.0;
        }
        else if (line.find("Ps4_factor") != std::string::npos) {
            std::string key, eq;
            int factor;
            iss >> key >> eq >> factor;
            values.ps_factor = factor;
        }
	else if (line.find("E SING FID TRACK EFFIC") != std::string::npos) {
	    size_t colon_pos = line.find(':');
	    if (colon_pos != std::string::npos) {
	        std::istringstream val_stream(line.substr(colon_pos + 1));
	        double eff;
	        val_stream >> eff;
	        values.hms_eff = eff;
	    }
	}
    }

    file.close();
    return values;
}

#endif // REPORT_PARSER_H
