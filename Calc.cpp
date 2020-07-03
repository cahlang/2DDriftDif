#include "DataTypes.h"
#include "Constants.h"
#include "Morphology.h"

#include <cmath>
#include <vector>

namespace calc{

	double bernou(double x){

		if (fabs(x) < 1E-8){
			return 1.0;
		}
		else{
			return x / (exp(x) - 1.0);
		}
	}
	
	MeasuredCurrent outer_circuit_current(Morphology &device_parameters, const PositionDependentParameter &electron_current, const PositionDependentParameter &hole_current, 
		const PositionDependentParameter &displacement_current, const PositionDependentParameter &net_rate){

		MeasuredCurrent current;
		std::vector<double> row_current (net_rate.points_x, 0.0);
		double root_mean_squared = 0.0;
		current.outer_circuit = 0.0;
		current.rms_error = 0.0;

		int number_of_rows = 0;
		for (int i = 0; i < net_rate.points_x; i++){
			for (int j = 0; j < net_rate.points_y; j++){

				int site = i*net_rate.points_y + j;
				int site_x_plus = site + net_rate.points_y;
				if ((device_parameters.is_electrode(site) || device_parameters.is_electrode(site_x_plus))){
					number_of_rows--;
					break;
				}
				// DO SOMETHING ABOUT THIS!
				row_current[i] += electron_current.data[site] + hole_current.data[site] + displacement_current.data[site];
				if (abs(net_rate.data[site_x_plus] * net_rate.spacing_x) >= abs(1E-7*(electron_current.data[site] + hole_current.data[site])) && abs(net_rate.data[site] * net_rate.spacing_x) >= abs(1E-7*(electron_current.data[site] + hole_current.data[site])))
				row_current[i] += (net_rate.data[site_x_plus] - net_rate.data[site]) * net_rate.spacing_x * 0.5;
				


			}
			number_of_rows++;
			current.outer_circuit += row_current[i]/(double)net_rate.points_y;

		}

		current.outer_circuit /= ((double)number_of_rows);
		
		for (int i = 0; i < net_rate.points_x; i++){
			if (row_current[i] != 0.0){
				root_mean_squared += pow(current.outer_circuit*(double)net_rate.points_y - row_current[i], 2.0);
			}
		}

		current.rms_error = pow(root_mean_squared / (double)(number_of_rows*net_rate.points_y), 0.5);

		current.normalization_coef = electron_current.normalization_coef;

		return current;

	}

}
