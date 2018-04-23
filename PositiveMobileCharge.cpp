#include "PositiveMobileCharge.h"
#include "Potential.h"
#include "Calc.h"
#include "Constants.h"
#include "Morphology.h"
#include "Misc.h"

#include <iostream>
#include <cmath>
#include <string>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

void PositiveMobileCharge::initialize(pt::ptree &settings){

	double norm_potential, norm_concentration, norm_current, norm_rate, norm_time;

	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::lenght_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::lenght_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::lenght_norm);
	norm_time = pow(constant::lenght_norm, 2.0) / (norm_potential * constant::mobility_norm);

	concentration.initialize(norm_concentration);
	current_x.initialize(norm_current);
	current_y.initialize(norm_current);
	diff_current_x.initialize(norm_current);
	diff_current_y.initialize(norm_current);
	drift_current_x.initialize(norm_current);
	drift_current_y.initialize(norm_current);
	rate.initialize(norm_rate);
	rate_coef.initialize(1.0 / norm_time);

	c_convergence = misc::to_array(settings.get<std::string>("numerics.hole_convergence_criteria"));
	c_relaxation = misc::to_array(settings.get<std::string>("numerics.hole_relaxation_coefficient"));

	initial_guess = settings.get<int>("numerics.initial_guess", -1);

	std::cout << "Hole concentration initialized." << std::endl;

	return;

}

void PositiveMobileCharge::solve(Morphology &material, const Potential &potential){

	// This function is designed so that the soultion is found independetly of the number of dimensions.
	// The equation to be solved can be written as Ax=b. This functions then applies a sucessive over-relaxation scheme in order to solve the system of equations.

	flag_has_converged = true;

	if (sweep_direction == 1){
		for (int i = 0; i < concentration.points_x; i++){
			for (int j = 0; j < concentration.points_y; j++){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_concentration(i, j, material, potential);
				}
			}
		}
	}
	else if (sweep_direction == -1){
		for (int i = concentration.points_x - 1; i >= 0; i--){
			for (int j = concentration.points_y - 1; j >= 0; j--){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_concentration(i, j, material, potential);
				}
			}
		}
	}
	else{
		/* TODO */
		// Error output
	}
	change_sweep_direction();

	if (material.get_convergence_boundary_condition()){
		for (int site = 0; site < material.cbc_interface.size(); site++){
			concentration.data[material.cbc_interface[site][0]] = concentration.data[material.cbc_interface[site][1]];
		}
	}

	return;

}

void PositiveMobileCharge::calculate_concentration(int i, int j, Morphology &material, const Potential &potential){

	int site, site_x_minus, site_x_plus;
	std::vector<double> a(5, 0.0);
	double b;

	site = i * concentration.points_y + j;

	site_x_minus = site - concentration.points_y;
	site_x_plus = site + concentration.points_y;

	double previous = concentration.data[site];

	if (material.is_electrode_interface(site)){
		get_boundary_condition(site, material.get_electrode_material_interface_number(site), material, potential.electrical);
	}
	else{

		a[0] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_HOMO(site_x_minus) - potential.electrical.data[site] - material.get_HOMO(site))
			+ material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_HOMO(site_x_plus) - potential.electrical.data[site] - material.get_HOMO(site))) + rate_coef.data[site];

		a[1] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_x_minus] - material.get_HOMO(site_x_minus)) * concentration.data[site_x_minus];

		a[2] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_x_plus] - material.get_HOMO(site_x_plus)) * concentration.data[site_x_plus];

		if (concentration.points_y > 1){

			int site_y_minus, site_y_plus;

			if (j != 0)
				site_y_minus = site - 1;
			else
				site_y_minus = site - 1 + concentration.points_y;

			if (j != concentration.points_y - 1)
				site_y_plus = site + 1;
			else
				site_y_plus = site + 1 - concentration.points_y;


			a[0] -= 1.0 / pow(concentration.spacing_y, 2.0) * (material.get_hole_mobility(site, site_y_minus) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_HOMO(site_y_minus) - potential.electrical.data[site] - material.get_HOMO(site))
				+ material.get_hole_mobility(site, site_y_plus) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_HOMO(site_y_plus) - potential.electrical.data[site] - material.get_HOMO(site)));

			a[3] = material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_y_minus] - material.get_HOMO(site_y_minus)) * concentration.data[site_y_minus];

			a[4] = material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_y_plus] - material.get_HOMO(site_y_plus)) * concentration.data[site_y_plus];

		}

		b = -rate.data[site];



		// Now, the next iterative step for the concentrationentration can be calculated.
		concentration.data[site] = (1.0 - c_relaxation[iteration_stage]) * concentration.data[site] + c_relaxation[iteration_stage] * (b - (a[1] + a[2] + a[3] + a[4])) / a[0];
	}

	// Finally, convergence is checked.
	if (fabs(concentration.data[site] - previous) >= c_convergence[iteration_stage] * concentration.data[site]){
		flag_has_converged = false;
	}


}

void PositiveMobileCharge::calculate_current(Morphology &material, const Potential potential){

	int site, site_x_plus;

	for (int i = 0; i < concentration.points_x; i++){
		for (int j = 0; j < concentration.points_y; j++){

			site = i*concentration.points_y + j;
			site_x_plus = site + concentration.points_y;

			if (i < concentration.points_x - 1 && concentration.calculate_this[site] && concentration.calculate_this[site_x_plus]){

				current_x.data[site] = -material.get_hole_mobility(site, site_x_plus) / current_x.spacing_x * (calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_x_plus] - material.get_HOMO(site_x_plus)) * concentration.data[site_x_plus]
					- calc::bernou(potential.electrical.data[site_x_plus] + material.get_HOMO(site_x_plus) - potential.electrical.data[site] - material.get_HOMO(site)) * concentration.data[site]);
			}
			if (current_y.points_y > 1){
				int site_y_plus;

				if (j != current_y.points_y-1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - concentration.points_y;
				if (concentration.calculate_this[site] && concentration.calculate_this[site_y_plus]){
					current_y.data[site] = -material.get_hole_mobility(site, site_y_plus) / current_y.spacing_y * (calc::bernou(potential.electrical.data[site] + material.get_HOMO(site) - potential.electrical.data[site_y_plus] - material.get_HOMO(site_y_plus)) * concentration.data[site_y_plus]
						- calc::bernou(potential.electrical.data[site_y_plus] + material.get_HOMO(site_y_plus) - potential.electrical.data[site] - material.get_HOMO(site)) * concentration.data[site]);
				}
			}
		}

	}

	return;
}

void PositiveMobileCharge::set_boundary_conditions(Morphology &device_properties, const PositionDependentParameter &electric_potential){
	
	for (int site = 0; site < concentration.points_x*concentration.points_y; site++){
		if (device_properties.is_electrode(site)){
			concentration.data[site] = 0.0;
			concentration.calculate_this[site] = false;
		}
	}

	if (device_properties.get_convergence_boundary_condition()){
		for (int site = 0; site < concentration.points_x*concentration.points_y; site++){
			if (device_properties.is_cbc(site)){
				concentration.data[site] = 0.0;
				concentration.calculate_this[site] = false;
			}
		}
	}

	for (int n = 0; n < device_properties.electrode_material_interface.size(); n++){
		ElectrodeMaterialInterface EMinterface = device_properties.get_electrode_material_interface(n);

		for (int m = 0; m < EMinterface.sites.size(); m++){
			int site = EMinterface.sites[m];
			get_boundary_condition(site, n, device_properties, electric_potential);
		}	

	}


	return;


}

void PositiveMobileCharge::get_boundary_condition(int site, int electrode_material_interface_number, Morphology &device_properties, const PositionDependentParameter &electric_potential){

	int electrode_number = device_properties.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
	double work_function = device_properties.get_work_function(electrode_number);

	if (device_properties.get_boundary_condition(electrode_material_interface_number) == 0){
		concentration.data[site] = device_properties.get_DOS(site) * exp(work_function - device_properties.get_HOMO(site));
		concentration.calculate_this[site] = false;
		// This is the thermionic expression, which is simply a constant.
	}
	else if (device_properties.get_boundary_condition(electrode_material_interface_number) == 1){
		concentration.data[site] = device_properties.get_DOS(site) * std::min(exp(work_function + device_properties.get_electrode_potential(electrode_number)
			- device_properties.get_HOMO(site) - electric_potential.data[site]),1.0);
		// This is calculated in each iterative step, insert this part into the solve scheme.
	}
}

void PositiveMobileCharge::set_initial_guess(Morphology &device_properties, const PositionDependentParameter &electric_potential, std::string file_name_template, int file_number){

	if (initial_guess == 0){
		PositionDependentParameter guess_from_file;
		guess_from_file.initialize(concentration.normalization_coef);
		guess_from_file.read_from_file(file_name_template, file_number);
		for (int site = 0; site < concentration.points_x*concentration.points_y; site++){
			if (concentration.calculate_this[site] == true){
				concentration.data[site] = guess_from_file.data[site];
			}
		}
	}
	else if (initial_guess == 1){
		for (int i = 0; i < concentration.points_x; i++){
			for (int j = 0; j < concentration.points_y; j++){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					double N_an = device_properties.get_DOS(site) * exp(device_properties.get_work_function(0) - device_properties.get_HOMO(site));
					double N_cat = device_properties.get_DOS(site) * exp(device_properties.get_work_function(1) - device_properties.get_HOMO(site));

					concentration.data[site] = N_an * pow((N_cat / N_an), ((double)i) / ((double)concentration.points_x - 1.0));

				}
			}
		}

	}
	else if (initial_guess == 2){
		for (int i = 0; i < concentration.points_x; i++){
			for (int j = 0; j < concentration.points_y; j++){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true && !device_properties.is_electrode_interface(site)){
					if (device_properties.get_material_number(device_properties.get_lattice_number(site)) == 0){
						concentration.data[site] = device_properties.get_DOS(site) * std::min(exp(device_properties.get_work_function(1) - device_properties.get_HOMO(site) - electric_potential.data[site]), 1.0);
					}
					else if (device_properties.get_material_number(device_properties.get_lattice_number(site)) == 1){
						concentration.data[site] = device_properties.get_DOS(site) * std::min(exp(device_properties.get_work_function(1) - device_properties.get_HOMO(site) - electric_potential.data[site]), 1.0);
					}
				}
			}
		}
	}
	else if (initial_guess == -1){
		std::cerr << "Initial guess was not read from file." << std::endl;
	}
	else{
		std::cerr << "The chosen initial guess number could not be found" << std::endl;
	}

	if (device_properties.get_convergence_boundary_condition()){
		for (int site = 0; site < device_properties.cbc_interface.size(); site++){
			concentration.data[device_properties.cbc_interface[site][0]] = concentration.data[device_properties.cbc_interface[site][1]];
			concentration.calculate_this[device_properties.cbc_interface[site][0]] = false;
		}
	}

	return;
}

//Ion-related functions

void PositiveMobileCharge::ion_solve(Morphology &material, const Potential &potential){

	// This function is designed so that the soultion is found independetly of the number of dimensions.
	// The equation to be solved can be written as Ax=b. This functions then applies a sucessive over-relaxation scheme in order to solve the system of equations.

	flag_has_converged = true;

	if (sweep_direction == 1){
		for (int i = 0; i < concentration.points_x; i++){
			for (int j = 0; j < concentration.points_y; j++){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_ion_concentration(i, j, material, potential);
				}
			}
		}
	}
	else if (sweep_direction == -1){
		for (int i = concentration.points_x - 1; i >= 0; i--){
			for (int j = concentration.points_y - 1; j >= 0; j--){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_ion_concentration(i, j, material, potential);
				}
			}
		}
	}
	else{
		/* TODO */
		// Error output
	}
	change_sweep_direction();

	if (material.get_convergence_boundary_condition()){
		for (int site = 0; site < material.cbc_interface.size(); site++){
			concentration.data[material.cbc_interface[site][0]] = concentration.data[material.cbc_interface[site][1]];
		}
	}
	for (int i = 1; i < concentration.points_x - 1; i++){
		int site = i*concentration.points_y + concentration.points_y - 1;
		concentration.data[site] = material.get_positive_ion_eq_conc(site);
	}

	return;

}

void PositiveMobileCharge::calculate_ion_concentration(int i, int j, Morphology &material, const Potential &potential){

	int site, site_x_minus, site_x_plus;
	std::vector<double> a(5, 0.0);
	double b;

	site = i * concentration.points_y + j;

	site_x_minus = site - concentration.points_y;
	site_x_plus = site + concentration.points_y;

	double previous = concentration.data[site];



		a[0] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_positive_ion_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] - potential.electrical.data[site])
			+ material.get_positive_ion_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] - potential.electrical.data[site]));

		a[1] = material.get_positive_ion_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_x_minus]) * concentration.data[site_x_minus];

		a[2] = material.get_positive_ion_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_x_plus]) * concentration.data[site_x_plus];

		if (concentration.points_y > 1){

			int site_y_minus, site_y_plus;

			if (j != 0)
				site_y_minus = site - 1;
			else
				site_y_minus = site - 1 + concentration.points_y;

			if (j != concentration.points_y - 1)
				site_y_plus = site + 1;
			else
				site_y_plus = site + 1 - concentration.points_y;


			a[0] -= 1.0 / pow(concentration.spacing_y, 2.0) * (material.get_positive_ion_mobility(site, site_y_minus) * calc::bernou(potential.electrical.data[site_y_minus] - potential.electrical.data[site])
				+ material.get_positive_ion_mobility(site, site_y_plus) * calc::bernou(potential.electrical.data[site_y_plus] - potential.electrical.data[site]));

			a[3] = material.get_positive_ion_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_y_minus]) * concentration.data[site_y_minus];

			a[4] = material.get_positive_ion_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_y_plus]) * concentration.data[site_y_plus];

		}

		b = 0.0;



		// Now, the next iterative step for the concentrationentration can be calculated.
		concentration.data[site] = (1.0 - c_relaxation[iteration_stage]) * concentration.data[site] + c_relaxation[iteration_stage] * (b - (a[1] + a[2] + a[3] + a[4])) / a[0];
	

	// Finally, convergence is checked.
	if (fabs(concentration.data[site] - previous) >= c_convergence[iteration_stage] * concentration.data[site]){
		flag_has_converged = false;
	}


}

void PositiveMobileCharge::calculate_ion_current(Morphology &material, const Potential &potential){

	int site, site_x_plus;

	for (int i = 0; i < concentration.points_x; i++){
		for (int j = 0; j < concentration.points_y; j++){

			site = i*concentration.points_y + j;
			site_x_plus = site + concentration.points_y;

			if (i < concentration.points_x - 1 && concentration.calculate_this[site] && concentration.calculate_this[site_x_plus] && material.positive_ion_transport(site) && material.positive_ion_transport(site_x_plus)){

				current_x.data[site] = -material.get_positive_ion_mobility(site, site_x_plus) / current_x.spacing_x * (calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_x_plus]) * concentration.data[site_x_plus]
					- calc::bernou(potential.electrical.data[site_x_plus] - potential.electrical.data[site]) * concentration.data[site]);
			}
			if (current_y.points_y > 1){
				int site_y_plus;

				if (j != current_y.points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - concentration.points_y;
				if (concentration.calculate_this[site] && concentration.calculate_this[site_y_plus] && material.positive_ion_transport(site) && material.positive_ion_transport(site_y_plus)){
					current_y.data[site] = -material.get_positive_ion_mobility(site, site_y_plus) / current_y.spacing_y * (calc::bernou(potential.electrical.data[site] - potential.electrical.data[site_y_plus]) * concentration.data[site_y_plus]
						- calc::bernou(potential.electrical.data[site_y_plus] - potential.electrical.data[site]) * concentration.data[site]);
				}
			}
		}

	}

	return;
}

void PositiveMobileCharge::set_ion_boundary_conditions(Morphology &device_properties){

	for (int site = 0; site < concentration.points_x*concentration.points_y; site++){
		if (!device_properties.positive_ion_transport(site)){
			concentration.data[site] = 0.0;
			concentration.calculate_this[site] = false;
		}
	}

	for (int i = 1; i < concentration.points_x - 1; i++){
		int site = i*concentration.points_y + concentration.points_y - 1;
		concentration.data[site] = device_properties.get_positive_ion_eq_conc(site);
	}

	return;

}

void PositiveMobileCharge::change_sweep_direction(){
	if (sweep_direction == 1){
		sweep_direction = -1;
	}
	else{
		sweep_direction = 1;
	}
	return;
}

void PositiveMobileCharge::next_iterative_stage(){
	iteration_stage++;
	return;
}

bool PositiveMobileCharge::final_iterative_stage_reached(){

	if (iteration_stage == c_convergence.size() - 1)
		return true;
	else
		return false;

}

void PositiveMobileCharge::reset_iterative_stage(){
	iteration_stage = 0;
	return;
}