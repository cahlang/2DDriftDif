#include "Misc.h"
#include "Constants.h"

#include <sstream>
#include <boost/lexical_cast.hpp>

namespace misc{


	std::vector<double> to_array(const std::string& string)
	{
		std::vector<double> result;
		std::stringstream ss(string);
		std::string item;
		while (std::getline(ss, item, ',')) result.push_back(boost::lexical_cast<double>(item));

		return result;
	}

	double single_level_fermi_dist(double fermi_level, double energy){
		return 1.0 / (exp(fermi_level - energy) + 1.0);
	}
}