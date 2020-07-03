#include "Potential.h"
#include "Constants.h"
#include "Morphology.h"
#include "Misc.h"

#include <iostream>
#include <cmath>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <Eigen/Sparse>


namespace pt = boost::property_tree;
using namespace Eigen;

void Potential::initialize(pt::ptree &settings){

	double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	double norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	double norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;

	chemical.initialize(norm_potential);
	electrical.initialize(norm_potential);
	electrochemical_electron.initialize(norm_potential);
	electrochemical_hole.initialize(norm_potential);

	displacement_current.initialize(norm_current);

	c_convergence = misc::to_array(settings.get<std::string>("numerics.potential_convergence_criteria"));
	c_relaxation = misc::to_array(settings.get<std::string>("numerics.potential_relaxation_coefficient"));

	initial_guess = settings.get<int>("numerics.initial_guess", -1);

	std::cout << "Potential initialized." << std::endl;

	return;

}

void Potential::solve_inverse(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration){

	PositionDependentParameter previous = electrical;

	flag_has_converged = true;

	int size = (electrical.points_x-2)*electrical.points_y;
	SparseMatrix<double> A(size, size);
	A.reserve(VectorXi::Constant(size, 6));
	VectorXd b(size), x(size);

	for (int i = 1; i < electrical.points_x-1; i++){
		for (int j = 0; j < electrical.points_y; j++){

			int site, site_x_minus, site_x_plus, site_y_minus, site_y_plus;

			site = i * electrical.points_y + j; 

			if (electrical.calculate_this[site] == true){

				site_x_minus = site - electrical.points_y;
				site_x_plus = site + electrical.points_y;
				if (j != 0)
					site_y_minus = site - 1;
				else
					site_y_minus = site - 1 + electrical.points_y;
				if (j != electrical.points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - electrical.points_y;

				A.insert(site - electrical.points_y, site - electrical.points_y) = -(material.get_relative_permittivity(site, site_x_minus) + material.get_relative_permittivity(site, site_x_plus)) / pow(electrical.spacing_x, 2.0)
					- (material.get_relative_permittivity(site, site_y_minus) + material.get_relative_permittivity(site, site_y_plus)) / pow(electrical.spacing_y, 2.0) + (electron_concentration.data[site] + hole_concentration.data[site]);
				if (i > 1)
					A.insert(site - electrical.points_y, site_x_minus - electrical.points_y) = material.get_relative_permittivity(site, site_x_minus) / pow(electrical.spacing_x, 2.0);

				if (i < electrical.points_x - 2)
					A.insert(site - electrical.points_y, site_x_plus - electrical.points_y) = material.get_relative_permittivity(site, site_x_plus) / pow(electrical.spacing_x, 2.0);

				A.insert(site - electrical.points_y, site_y_minus - electrical.points_y) = material.get_relative_permittivity(site, site_y_minus) / pow(electrical.spacing_y, 2.0);
				A.insert(site - electrical.points_y, site_y_plus - electrical.points_y) = material.get_relative_permittivity(site, site_y_plus) / pow(electrical.spacing_y, 2.0);

				b(site - electrical.points_y) = electron_concentration.data[site] - hole_concentration.data[site] - material.dopants.data[site] + material.MG_electron_concentration.data[site]
					- material.MG_hole_concentration.data[site] + (electron_concentration.data[site] + hole_concentration.data[site]) * previous.data[site];
				if (i == 1)
					b(site - electrical.points_y) += -material.get_relative_permittivity(site, site_x_minus) * material.get_electrode_potential(0) / pow(electrical.spacing_x, 2.0);
				if (i == electrical.points_x - 2)
					b(site - electrical.points_y) += -material.get_relative_permittivity(site, site_x_plus) * material.get_electrode_potential(1) / pow(electrical.spacing_x, 2.0);

				if (material.negative_ion_transport(site)){
					b(site - electrical.points_y) += negative_ion_concentration.data[site];
				}
				if (material.positive_ion_transport(site)){
					b(site - electrical.points_y) -= positive_ion_concentration.data[site];
				}

			}
		}
	}

	SparseLU<SparseMatrix<double>> solver;
	A.makeCompressed();
	solver.compute(A);

	if (solver.info() != Success) {
		std::cerr << "Factorization failed" << std::endl;
		return;
	}
	x = solver.solve(b);

	for (int site = 0; site < size; site++){
		if (electrical.calculate_this[site + electrical.points_y] == true)
			electrical.data[site + electrical.points_y] = x(site);
		if (fabs(electrical.data[site + electrical.points_y] - previous.data[site + electrical.points_y]) >= c_convergence[iteration_stage] * fabs(electrical.data[site + electrical.points_y]) + 1E-15){
			flag_has_converged = false;

		}

	}
	return;

}

void Potential::solve(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration){

	// This function is designed so that the soultion is found independetly of the number of dimensions.
	// The equation to be solved can be written as Ax=b. This functions then applies a sucessive over-relaxation scheme in order to solve the system of equations.

	flag_has_converged = true;
	
	if (sweep_direction == 1){
		for (int i = 0; i < electrical.points_x; i++){
			for (int j = 0; j < electrical.points_y; j++){
				int site = i*electrical.points_y + j;
				if (electrical.calculate_this[site] == true){
					calculate_electrical_potential(i, j, material, electron_concentration, hole_concentration, negative_ion_concentration, positive_ion_concentration);
				}
			}
		}
	}
	else if (sweep_direction == -1){
		for (int i = electrical.points_x - 1; i >= 0; i--){
			for (int j = electrical.points_y - 1; j >= 0; j--){
				int site = i*electrical.points_y + j;
				if (electrical.calculate_this[site] == true){
					calculate_electrical_potential(i, j, material, electron_concentration, hole_concentration, negative_ion_concentration, positive_ion_concentration);
				}
			}
		}
	}
	change_sweep_direction();

	if (material.get_neumann_zero_boundary()){
		for (int site = 0; site < material.cbc_interface.size(); site++){
			electrical.data[material.cbc_interface[site][0]] = electrical.data[material.cbc_interface[site][1]];
		}
	}


	return;

}

void Potential::solve_one_d(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration){

	PositionDependentParameter previous_solution = electrical;

	flag_has_converged = true;

	std::vector<double> a(electrical.points_x - 2, 0.0), b(electrical.points_x - 2, 0.0), c(electrical.points_x - 2, 0.0), d(electrical.points_x - 2, 0.0);

	for (int i = 0; i < b.size(); i++){

		int site = i + 1;
		int site_x_minus = i;
		int site_x_plus = i + 2;

		if (i != 0)
		a[i] = material.get_relative_permittivity(site, site_x_minus) / pow(electrical.spacing_x, 2.0);

		b[i] = -(material.get_relative_permittivity(site, site_x_minus) + material.get_relative_permittivity(site, site_x_plus)) / pow(electrical.spacing_x, 2.0) - (electron_concentration.data[site] + hole_concentration.data[site]);
		if (i != b.size()-1)
		c[i] = material.get_relative_permittivity(site, site_x_plus) / pow(electrical.spacing_x, 2.0);

		d[i] = electron_concentration.data[site] - hole_concentration.data[site] - material.dopants.data[site] + material.MG_electron_concentration.data[site]
			- material.MG_hole_concentration.data[site] - (electron_concentration.data[site] + hole_concentration.data[site])*previous_solution.data[site];
	}


	std::vector<double> c_prime(electrical.points_x - 2, 0.0), d_prime(electrical.points_x - 2, 0.0);

	int lower_bound = 0;
	int upper_bound = b.size() - 1;

	if (electrical.calculate_this[1] == true)
		d[0] -= electrical.data[0] * material.get_relative_permittivity(1, 0) / pow(electrical.spacing_x, 2.0);
	else{
		d[1] -= electrical.data[1] * material.get_relative_permittivity(2, 1) / pow(electrical.spacing_x, 2.0);
		lower_bound = 1;
	}

	if (electrical.calculate_this[electrical.points_x-2])
		d[d.size() - 1] -= electrical.data[d.size() + 1] * material.get_relative_permittivity(d.size(), d.size() + 1) / pow(electrical.spacing_x, 2.0);
	else{
		d[d.size() - 2] -= electrical.data[d.size()] * material.get_relative_permittivity(d.size() - 1, d.size()) / pow(electrical.spacing_x, 2.0);
		upper_bound = b.size() - 2;
	}

	for (int n = lower_bound; n <= upper_bound; n++){
		if (n == lower_bound)
			c_prime[n] = c[n] / b[n];
		else if (n < upper_bound)
			c_prime[n] = c[n] / (b[n] - a[n] * c_prime[n - 1]);

		if (n == lower_bound)
			d_prime[n] = d[n] / b[n];
		else
			d_prime[n] = (d[n] - a[n] * d_prime[n - 1]) / (b[n] - a[n] * c_prime[n - 1]);
	}

	for (int n = upper_bound; n >= lower_bound; n--){
		int i = n + 1;
		if (n == upper_bound)
			electrical.data[i] = d_prime[n];
		else
			electrical.data[i] = d_prime[n] - c_prime[n] * electrical.data[i + 1];
	}

	for (int i = electrical.points_x - 1; i >= 0; i--){
		for (int j = electrical.points_y - 1; j >= 0; j--){

			int site = i*electrical.points_y + j;

			if (fabs(electrical.data[site] - previous_solution.data[site]) >= c_convergence[iteration_stage] * fabs(electrical.data[site]) + 1E-15){
				flag_has_converged = false;
			}

		}
	}

	return;

}

void Potential::calculate_electrical_potential(int i, int j, Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration){

	int site, site_x_minus, site_x_plus;
	std::vector<double> a(5, 0.0);
	double b = 0.0;

	site = i * electrical.points_y + j;
	site_x_minus = site - electrical.points_y;
	site_x_plus = site + electrical.points_y;

	double previous = electrical.data[site];

	a[0] = -(material.get_relative_permittivity(site, site_x_minus) + material.get_relative_permittivity(site, site_x_plus)) / pow(electrical.spacing_x, 2.0) - (electron_concentration.data[site] + hole_concentration.data[site]); // The Gummel scheme can be added for increased stability (2.0 -> 2.0 + electron con + hole con).

	a[1] = material.get_relative_permittivity(site, site_x_minus) * electrical.data[site_x_minus] / pow(electrical.spacing_x, 2.0);

	a[2] = material.get_relative_permittivity(site, site_x_plus) * electrical.data[site_x_plus] / pow(electrical.spacing_x, 2.0);

	b = electron_concentration.data[site] - hole_concentration.data[site] - material.dopants.data[site] + material.MG_electron_concentration.data[site]
		- material.MG_hole_concentration.data[site] - (electron_concentration.data[site] + hole_concentration.data[site])*electrical.data[site];

	if (material.negative_ion_transport(site)){
		b += negative_ion_concentration.data[site];
	}
	if (material.positive_ion_transport(site)){
		b -= positive_ion_concentration.data[site];
	}

	if (i == 0){
		a[1] = 0.0;
		b += -material.get_relative_permittivity(site, site_x_minus) * material.get_electrode_potential(0) / pow(electrical.spacing_x, 2.0);
	}
	if (i == electrical.points_x - 1){
		a[2] = 0.0;
		b += -material.get_relative_permittivity(site, site_x_plus) * material.get_electrode_potential(1) / pow(electrical.spacing_x, 2.0);
	}

	// If the number of points in the x-direction is larger than 1, the 2D case is assumed. Then the following terms needs to be added:
	if (electrical.points_y > 1){
		int site_y_minus, site_y_plus;
		if (j != 0)
			site_y_minus = site - 1;
		else
			site_y_minus = site - 1 + electrical.points_y;
		if (j != electrical.points_y - 1)
			site_y_plus = site + 1;
		else
			site_y_plus = site + 1 - electrical.points_y;

		a[0] += -(material.get_relative_permittivity(site, site_y_minus) + material.get_relative_permittivity(site, site_y_plus)) / pow(electrical.spacing_y, 2.0);
		a[3] = material.get_relative_permittivity(site, site_y_minus) * electrical.data[site_y_minus] / pow(electrical.spacing_y, 2.0);
		a[4] = material.get_relative_permittivity(site, site_y_plus) * electrical.data[site_y_plus] / pow(electrical.spacing_y, 2.0);

	}

	// Now, the next iterative step for the electrical potential at "site" can be calculated.
	electrical.data[site] = (1.0 - c_relaxation[iteration_stage]) * electrical.data[site] + c_relaxation[iteration_stage] * (b - (a[1] + a[2] + a[3] + a[4])) / a[0];

	// Finally, convergence is checked.
	if (fabs(electrical.data[site] - previous) >= c_convergence[iteration_stage] * fabs(electrical.data[site])){
		flag_has_converged = false;
	}

}

void Potential::calculate_electrochemical_potential(Morphology material, PositionDependentParameter electron_concentration, PositionDependentParameter hole_concentration){

	for (int site = 0; site < electrochemical_electron.points_x*electrochemical_electron.points_y; site++){
		
		if (electrical.calculate_this[site]){
			electrochemical_electron.data[site] = -(material.get_electron_trans_energy(site) + electrical.data[site] - log(electron_concentration.data[site] / material.get_electron_trans_DOS(site)));
			electrochemical_hole.data[site] = -(material.get_hole_trans_energy(site) + electrical.data[site] + log(hole_concentration.data[site] / material.get_hole_trans_DOS(site)));
		}
	}

}

void Potential::calculate_displacement_current(Morphology &material, const PositionDependentParameter previous_potential, const double time_step){

	int site, site_x_plus;

	for (int i = 0; i < electrical.points_x; i++){
		for (int j = 0; j < electrical.points_y; j++){

			site = i*electrical.points_y + j;
			site_x_plus = site + electrical.points_y;

			if (i < electrical.points_x - 1 && !material.is_electrode(site) && !material.is_electrode(site_x_plus)){

				displacement_current.data[site] = -material.get_relative_permittivity(site, site_x_plus) * ((electrical.data[site_x_plus] - electrical.data[site]) / electrical.spacing_x - (previous_potential.data[site_x_plus] - previous_potential.data[site]) / electrical.spacing_x) / time_step;
			}
		}
	}

	return;
}

void Potential::set_boundary_conditions(Morphology device_properties){

	for (int site = 0; site < electrical.points_x*electrical.points_y; site++){
		if (device_properties.is_electrode(site)){
			int electrode_number = device_properties.get_electrode_number(device_properties.get_lattice_number(site));
			electrical.data[site] = device_properties.get_electrode_potential(electrode_number);
			electrical.calculate_this[site] = false;
		}
	}

	if (device_properties.get_neumann_zero_boundary()){
		for (int site = 0; site < electrical.points_x*electrical.points_y; site++){
			if (device_properties.is_cbc(site)){
				electrical.data[site] = 0.0;
				electrical.calculate_this[site] = false;
			}
		}
	}

	for (int n = 0; n < device_properties.electrode_material_interface.size(); n++){
		ElectrodeMaterialInterface EMinterface = device_properties.get_electrode_material_interface(n);
		int electrode_number = device_properties.electrode_material_interface[n].get_electrode_number();

		if (device_properties.get_boundary_condition(n) == 0){
			for (int m = 0; m < EMinterface.sites.size(); m++){

				int site = EMinterface.sites[m];
				electrical.data[site] = device_properties.get_electrode_potential(electrode_number);
				electrical.calculate_this[site] = false;
			}
		}
		else if (device_properties.get_boundary_condition(n) == 1){

			// Nothing here, only the electrode potential is set, the rest is calculated.

		}
		else if (device_properties.get_boundary_condition(n) == -1){
			std::cerr << "Boundary condition was not read from file." << std::endl;
		}
		else{
			std::cerr << "The chosen boundary condition number could not be found" << std::endl;
		}

	}

	return;
}

void Potential::set_initial_guess(Morphology device_properties, std::string file_name_template, int file_number){

	if (initial_guess == 0){
		PositionDependentParameter guess_from_file;
		guess_from_file.initialize(electrical.normalization_coef);
		guess_from_file.read_from_file(file_name_template, file_number);
		for (int site = 0; site < electrical.points_x*electrical.points_y; site++){
			if (electrical.calculate_this[site] == true){
				electrical.data[site] = guess_from_file.data[site];
			}
		}

	}
	else if (initial_guess == 1){
		for (int i = 0; i < electrical.points_x; i++){
			for (int j = 0; j < electrical.points_y; j++){
				int site = i*electrical.points_y + j;
				if (electrical.calculate_this[site] == true){
					electrical.data[site] = electrical.spacing_x * ((double)electrical.points_x - (double)i) / device_properties.get_active_layer_length_x() * (device_properties.get_electrode_potential(0) - device_properties.get_electrode_potential(1));
				}
			}
		}
	}
	else if (initial_guess == 2){
		for (int i = 0; i < electrical.points_x; i++){
			for (int j = 0; j < electrical.points_y; j++){
				int site = i*electrical.points_y + j;
				if (electrical.calculate_this[site] == true){
					if (device_properties.get_material_number(device_properties.get_lattice_number(site)) == 0){
						electrical.data[site] = device_properties.get_work_function(1) - device_properties.get_MG_plus(site);
					}
					else if (device_properties.get_material_number(device_properties.get_lattice_number(site)) == 1){
						electrical.data[site] = device_properties.get_work_function(1) - device_properties.get_MG_minus(site);
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

	if (device_properties.get_neumann_zero_boundary()){
		for (int site = 0; site < device_properties.cbc_interface.size(); site++){
			electrical.data[device_properties.cbc_interface[site][0]] = electrical.data[device_properties.cbc_interface[site][1]];
			electrical.calculate_this[device_properties.cbc_interface[site][0]] = false;
		}
	}

	return;

}

bool Potential::has_converged() {
	return flag_has_converged;
}

void Potential::change_sweep_direction(){
	if (sweep_direction == 1){
		sweep_direction = -1;
	}
	else{
		sweep_direction = 1;
	}
	return;
}

void Potential::next_iterative_stage(){
	iteration_stage++;
	return;
}

bool Potential::final_iterative_stage_reached(){

	if (iteration_stage == c_convergence.size() - 1)
		return true;
	else
		return false;

}

void Potential::reset_iterative_stage(){
	iteration_stage = 0;
	return;
}