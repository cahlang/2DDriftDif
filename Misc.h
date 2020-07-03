#pragma once

#include <vector>
#include <string>

namespace misc{

	std::vector<double> to_array(const std::string& s); // Separate numbers by commas!
	double single_level_fermi_dist(double fermi_level, double energy);

}