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
#include <Eigen/Sparse>

namespace pt = boost::property_tree;
using namespace Eigen;

typedef Triplet<double> Trip;

void PositiveMobileCharge::initialize(pt::ptree &settings){

	double norm_potential, norm_concentration, norm_current, norm_rate, norm_time;

	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::length_norm);
	norm_time = pow(constant::length_norm, 2.0) / (norm_potential * constant::mobility_norm);

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

void PositiveMobileCharge::solve_inverse(Morphology &material, const Potential &potential){

	flag_has_converged = true;
	std::vector<double> previous = concentration.data;

	int size = (concentration.points_x - 2)*concentration.points_y;
	SparseMatrix<double> A(size, size);
	A.reserve(VectorXi::Constant(size, 6));
	VectorXd b(size), x(size);

	std::vector<Trip> tripletList;
	tripletList.reserve(20 * size);

	for (int i = 1; i < concentration.points_x - 1; i++){
		for (int j = 0; j < concentration.points_y; j++){

			int site, site_x_minus, site_x_plus, site_y_minus, site_y_plus;

			site = i * concentration.points_y + j;

			if (concentration.calculate_this[site] == true){

				site_x_minus = site - concentration.points_y;
				site_x_plus = site + concentration.points_y;
				if (j != 0)
					site_y_minus = site - 1;
				else
					site_y_minus = site - 1 + concentration.points_y;
				if (j != concentration.points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - concentration.points_y;

				if (material.is_electrode_interface(site) && !material.surface_recombination_activated(material.get_electrode_material_interface_number(site))){
					// calculates the value of the concentration at site.
					int interface_number = material.get_electrode_material_interface_number(site);

					get_boundary_condition(site, interface_number, material, potential.electrical);

				}
				else{

					if (material.is_electrode_interface(site) && material.surface_recombination_activated(material.get_electrode_material_interface_number(site))){
						int electrode_material_interface_number = material.get_electrode_material_interface_number(site);
						if (material.is_electrode(site_x_minus)){

							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
								- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x));

							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus))));

							if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))){
								int pair_number = material.get_interface_pair(site, site_y_minus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else{
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
							}
							if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))){
								int pair_number = material.get_interface_pair(site, site_y_plus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else{
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
							}

							int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
							double work_function = material.get_work_function(electrode_number);
							double electrode_potential = material.get_electrode_potential(electrode_number);

							b(site - concentration.points_y) = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function - material.get_hole_trans_energy(site)),1.0);
						}
						else if (material.is_electrode(site_x_plus)){
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
								- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus))));

							if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))){
								int pair_number = material.get_interface_pair(site, site_y_minus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else{
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
							}
							if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))){
								int pair_number = material.get_interface_pair(site, site_y_plus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else{
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
							}

							int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
							double work_function = material.get_work_function(electrode_number);
							double electrode_potential = material.get_electrode_potential(electrode_number);

							b(site - concentration.points_y) = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function - material.get_hole_trans_energy(site)), 1.0);

						}
					}
					else{

						if ((material.is_interface(site) && material.is_interface(site_x_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_minus)))){
							int pair_number = material.get_interface_pair(site, site_x_minus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site_x_minus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_x_minus) + potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else{
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus))));
						}

						if ((material.is_interface(site) && material.is_interface(site_x_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_plus)))){
							int pair_number = material.get_interface_pair(site, site_x_plus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else{
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus))));
						}
						if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))){
							int pair_number = material.get_interface_pair(site, site_y_minus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else{
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
						}
						if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))){
							int pair_number = material.get_interface_pair(site, site_y_plus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else{
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
						}
						tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, rate_coef.data[site]));
						b(site - concentration.points_y) = -rate.data[site];
					}



				}

			}
		}
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());

	SparseLU<SparseMatrix<double>>   solver;
	A.makeCompressed();

	solver.compute(A);

	if (solver.info() != Success) {
		std::cerr << "Holecon decomposition failed" << std::endl;
		return;
	}

	x = solver.solve(b);

	for (int site = 0; site < size; site++){
		if (concentration.calculate_this[site + concentration.points_y] == true)
		concentration.data[site + concentration.points_y] = x(site);
		if (fabs(concentration.data[site + concentration.points_y] - previous[site + concentration.points_y]) >= c_convergence[iteration_stage] * concentration.data[site + concentration.points_y]){
			flag_has_converged = false;
		}
	}


	return;


}

void PositiveMobileCharge::solve_inverse(Morphology& material, const Potential& potential, const PositionDependentParameter& previous_time_concentration, double time_step) {

	flag_has_converged = true;
	std::vector<double> previous = concentration.data;

	int size = (concentration.points_x - 2) * concentration.points_y;
	SparseMatrix<double> A(size, size);
	A.reserve(VectorXi::Constant(size, 6));
	VectorXd b(size), x(size);

	std::vector<Trip> tripletList;
	tripletList.reserve(20 * size);

	for (int i = 1; i < concentration.points_x - 1; i++) {
		for (int j = 0; j < concentration.points_y; j++) {

			int site, site_x_minus, site_x_plus, site_y_minus, site_y_plus;

			site = i * concentration.points_y + j;

			if (concentration.calculate_this[site] == true) {

				site_x_minus = site - concentration.points_y;
				site_x_plus = site + concentration.points_y;
				if (j != 0)
					site_y_minus = site - 1;
				else
					site_y_minus = site - 1 + concentration.points_y;
				if (j != concentration.points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - concentration.points_y;

				if (material.is_electrode_interface(site) && !material.surface_recombination_activated(material.get_electrode_material_interface_number(site))) {
					// calculates the value of the concentration at site.
					int interface_number = material.get_electrode_material_interface_number(site);

					get_boundary_condition(site, interface_number, material, potential.electrical);

				}
				else {

					if (material.is_electrode_interface(site) && material.surface_recombination_activated(material.get_electrode_material_interface_number(site))) {
						int electrode_material_interface_number = material.get_electrode_material_interface_number(site);
						if (material.is_electrode(site_x_minus)) {

							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
								- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x));

							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus))));

							if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))) {
								int pair_number = material.get_interface_pair(site, site_y_minus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else {
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
							}
							if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))) {
								int pair_number = material.get_interface_pair(site, site_y_plus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else {
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
							}

							int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
							double work_function = material.get_work_function(electrode_number);
							double electrode_potential = material.get_electrode_potential(electrode_number);

							b(site - concentration.points_y) = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function + electrode_potential - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);
						}
						else if (material.is_electrode(site_x_plus)) {
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
								- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus))));

							if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))) {
								int pair_number = material.get_interface_pair(site, site_y_minus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else {
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
							}
							if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))) {
								int pair_number = material.get_interface_pair(site, site_y_plus);
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
							}
							else {
								tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
								tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
							}

							int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
							double work_function = material.get_work_function(electrode_number);
							double electrode_potential = material.get_electrode_potential(electrode_number);

							b(site - concentration.points_y) = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function + electrode_potential - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);

						}
					}
					else {

						if ((material.is_interface(site) && material.is_interface(site_x_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_minus)))) {
							int pair_number = material.get_interface_pair(site, site_x_minus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site_x_minus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_x_minus) + potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else {
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_minus - concentration.points_y, material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus))));
						}

						if ((material.is_interface(site) && material.is_interface(site_x_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_plus)))) {
							int pair_number = material.get_interface_pair(site, site_x_plus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else {
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_x_plus - concentration.points_y, material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus))));
						}
						if ((material.is_interface(site) && material.is_interface(site_y_minus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus)))) {
							int pair_number = material.get_interface_pair(site, site_y_minus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else {
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_minus - concentration.points_y, material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus))));
						}
						if ((material.is_interface(site) && material.is_interface(site_y_plus) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus)))) {
							int pair_number = material.get_interface_pair(site, site_y_plus);
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0)));
						}
						else {
							tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))));
							tripletList.push_back(Trip(site - concentration.points_y, site_y_plus - concentration.points_y, material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus))));
						}
						tripletList.push_back(Trip(site - concentration.points_y, site - concentration.points_y, rate_coef.data[site] - 1.0 / time_step));
						b(site - concentration.points_y) = -rate.data[site] - previous_time_concentration.data[site] / time_step;
					}



				}

			}
		}
	}

	A.setFromTriplets(tripletList.begin(), tripletList.end());

	SparseLU<SparseMatrix<double>>   solver;
	A.makeCompressed();

	solver.compute(A);

	if (solver.info() != Success) {
		std::cerr << "Holecon decomposition failed" << std::endl;
		return;
	}

	x = solver.solve(b);

	for (int site = 0; site < size; site++) {
		if (concentration.calculate_this[site + concentration.points_y] == true)
			concentration.data[site + concentration.points_y] = x(site);
		if (fabs(concentration.data[site + concentration.points_y] - previous[site + concentration.points_y]) >= c_convergence[iteration_stage] * concentration.data[site + concentration.points_y]) {
			flag_has_converged = false;
		}
	}


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


	if (material.get_neumann_zero_boundary()){
		for (int site = 0; site < material.cbc_interface.size(); site++){
			concentration.data[material.cbc_interface[site][0]] = concentration.data[material.cbc_interface[site][1]];
		}
	}

	return;

}


void PositiveMobileCharge::solve(Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentration, double time_step){

	// This function is designed so that the soultion is found independetly of the number of dimensions.
	// The equation to be solved can be written as Ax=b. This functions then applies a sucessive over-relaxation scheme in order to solve the system of equations.

	flag_has_converged = true;

	if (sweep_direction == 1){
		for (int i = 0; i < concentration.points_x; i++){
			for (int j = 0; j < concentration.points_y; j++){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_concentration(i, j, material, potential, previous_time_concentration.data[site], time_step);
				}
			}
		}
	}
	else if (sweep_direction == -1){
		for (int i = concentration.points_x - 1; i >= 0; i--){
			for (int j = concentration.points_y - 1; j >= 0; j--){
				int site = i*concentration.points_y + j;
				if (concentration.calculate_this[site] == true){
					calculate_concentration(i, j, material, potential, previous_time_concentration.data[site], time_step);
				}
			}
		}
	}
	change_sweep_direction();

	if (material.get_neumann_zero_boundary()){
		for (int site = 0; site < material.cbc_interface.size(); site++){
			concentration.data[material.cbc_interface[site][0]] = concentration.data[material.cbc_interface[site][1]];
		}
	}

	return;
}

void PositiveMobileCharge::solve_one_d(Morphology &material, const Potential &potential){

	flag_has_converged = true;

	std::vector<double> a(concentration.points_x - 2, 0.0), b(concentration.points_x - 2, 0.0), c(concentration.points_x - 2, 0.0), d(concentration.points_x - 2, 0.0);

	for (int i = 0; i < b.size(); i++){

		int site = i + 1;
		int site_x_minus = i;
		int site_x_plus = i + 2;

		if (i == 0){

			double electrode_material_interface_number = material.get_electrode_material_interface_number(site);
			if (concentration.calculate_this[site] == true && material.surface_recombination_activated(electrode_material_interface_number)){

				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
					- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));
				d[i] = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * exp(work_function - material.get_hole_trans_energy(site));
			}

		}
		else if (i == b.size() - 1){

			double electrode_material_interface_number = material.get_electrode_material_interface_number(site);
			if (concentration.calculate_this[site] == true && material.surface_recombination_activated(electrode_material_interface_number)){

				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
				b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site)))
					- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				d[i] = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * exp(work_function - material.get_hole_trans_energy(site));

			}
		}
		else{
			d[i] = -rate.data[site];
			b[i] = rate_coef.data[site];
			if (material.is_interface(site) && material.is_interface(site_x_minus)){
				int pair_number = material.get_interface_pair(site, site_x_minus);
				b[i] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site_x_minus]), 1.0);
				a[i] = material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_minus) + potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);
			}
			else{
				b[i] += -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
			}
			if (material.is_interface(site) && material.is_interface(site_x_plus)){
				int pair_number = material.get_interface_pair(site, site_x_plus);
				b[i] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0);
				c[i] = material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);
			}
			else{
				b[i] += -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));
			}


			/*
			a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
			b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
				+ material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))) + rate_coef.data[site];
			c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));

			*/	
		}

	}

	std::vector<double> c_prime(concentration.points_x - 2, 0.0), d_prime(concentration.points_x - 2, 0.0);

	int lower_bound = 0;
	int upper_bound = b.size() - 1;
	double electrode_material_interface_number = material.get_electrode_material_interface_number(1);
	if (!material.surface_recombination_activated(electrode_material_interface_number)){
		if(concentration.calculate_this[1] == true)
		get_boundary_condition(1, electrode_material_interface_number, material, potential.electrical);

		lower_bound = 1;
		d[1] -= concentration.data[1] * material.get_hole_mobility(2, 1) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[2] + material.get_hole_trans_energy(2) - potential.electrical.data[1] - material.get_hole_trans_energy(1));
	}
	electrode_material_interface_number = material.get_electrode_material_interface_number(concentration.points_x - 2);
	if (!material.surface_recombination_activated(electrode_material_interface_number)){
		if (concentration.calculate_this[concentration.points_x - 2] == true)
		get_boundary_condition(concentration.points_x - 2, electrode_material_interface_number, material, potential.electrical);

		upper_bound = b.size() - 2;
		d[b.size() - 2] -= concentration.data[b.size()] * material.get_hole_mobility(b.size() - 1, b.size()) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[b.size() - 1] + material.get_hole_trans_energy(b.size() - 1) - potential.electrical.data[b.size()] - material.get_hole_trans_energy(b.size()));
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
			concentration.data[i] = d_prime[n];
		else
			concentration.data[i] = d_prime[n] - c_prime[n] * concentration.data[i + 1];

		if (concentration.data[i] > material.get_hole_trans_DOS(i))
			concentration.data[i] = material.get_hole_trans_DOS(i);

		if (concentration.data[i] < 0.0)
			concentration.data[i] = 0.0;
	}

	return;

}

void PositiveMobileCharge::solve_one_d(Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentration, double time_step){

	flag_has_converged = true;

	std::vector<double> a(concentration.points_x - 2, 0.0), b(concentration.points_x - 2, 0.0), c(concentration.points_x - 2, 0.0), d(concentration.points_x - 2, 0.0);

	for (int i = 0; i < b.size(); i++){

		int site = i + 1;
		int site_x_minus = i;
		int site_x_plus = i + 2;

		if (i == 0){

			double electrode_material_interface_number = material.get_electrode_material_interface_number(site);
			if (concentration.calculate_this[site] == true && material.surface_recombination_activated(electrode_material_interface_number)){
			
				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
					-material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));
				d[i] = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * exp(work_function - material.get_hole_trans_energy(site));
			}

		}
		else if (i == b.size() - 1){

			double electrode_material_interface_number = material.get_electrode_material_interface_number(site);
			if (concentration.calculate_this[site] == true && material.surface_recombination_activated(electrode_material_interface_number)){

				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
				b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site)))
					- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				d[i] = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * exp(work_function - material.get_hole_trans_energy(site));

			}
		}
		else{

			d[i] = -rate.data[site] - previous_time_concentration.data[site] / time_step;
			b[i] = rate_coef.data[site] - 1.0 / time_step;
			if (material.is_interface(site) && material.is_interface(site_x_minus)){
				int pair_number = material.get_interface_pair(site, site_x_minus);
				b[i] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site_x_minus]), 1.0);
				a[i] = material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_minus) + potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);
			}
			else{
				b[i] += -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
			}
			if (material.is_interface(site) && material.is_interface(site_x_plus)){
				int pair_number = material.get_interface_pair(site, site_x_plus);
				b[i] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0);
				c[i] = material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0);
			}
			else{
				b[i] += -1.0 / pow(concentration.spacing_x, 2.0) * material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));
			}


			/*
			a[i] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus));
			b[i] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
				+ material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))) + rate_coef.data[site] - 1.0 / time_step;
			c[i] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus));
			d[i] = -rate.data[site] - previous_time_concentration.data[site] / time_step;
			*/
		}


	}

	std::vector<double> c_prime(concentration.points_x - 2, 0.0), d_prime(concentration.points_x - 2, 0.0);

	int lower_bound = 0;
	int upper_bound = b.size() - 1;
	double electrode_material_interface_number = material.get_electrode_material_interface_number(1);
	if (!material.surface_recombination_activated(electrode_material_interface_number)){
		if (concentration.calculate_this[1] == true)
			get_boundary_condition(1, electrode_material_interface_number, material, potential.electrical);

		lower_bound = 1;
		d[1] -= concentration.data[1] * material.get_hole_mobility(2, 1) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[2] + material.get_hole_trans_energy(2) - potential.electrical.data[1] - material.get_hole_trans_energy(1));
	}
	electrode_material_interface_number = material.get_electrode_material_interface_number(concentration.points_x - 2);
	if (!material.surface_recombination_activated(electrode_material_interface_number)){
		if (concentration.calculate_this[concentration.points_x - 2] == true)
			get_boundary_condition(concentration.points_x - 2, electrode_material_interface_number, material, potential.electrical);

		upper_bound = b.size() - 2;
		d[b.size() - 2] -= concentration.data[b.size()] * material.get_hole_mobility(b.size() - 1, b.size()) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[b.size() - 1] + material.get_hole_trans_energy(b.size() - 1) - potential.electrical.data[b.size()] - material.get_hole_trans_energy(b.size()));
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
			concentration.data[i] = d_prime[n];
		else
			concentration.data[i] = d_prime[n] - c_prime[n] * concentration.data[i + 1];

		if (concentration.data[i] > material.get_hole_trans_DOS(i))
			concentration.data[i] = material.get_hole_trans_DOS(i);

		if (concentration.data[i] < 0.0)
			concentration.data[i] = 0.0;

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

	if (material.is_electrode_interface(site) && !material.surface_recombination_activated(material.get_electrode_material_interface_number(site))){
		int interface_number = material.get_electrode_material_interface_number(site);
		get_boundary_condition(site, interface_number, material, potential.electrical);
	}
	else{

		if (material.is_electrode_interface(site) && material.surface_recombination_activated(material.get_electrode_material_interface_number(site))){
			int electrode_material_interface_number = material.get_electrode_material_interface_number(site);
			if (material.is_electrode(site_x_minus)){
				a[0] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))) + rate_coef.data[site]
					- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				a[2] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus)) * concentration.data[site_x_plus];

				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				b = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function - material.get_hole_trans_energy(site)), 1.0);
			}
			else if (material.is_electrode(site_x_plus)){
				a[0] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))) + rate_coef.data[site]
					- material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x;
				a[1] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus)) * concentration.data[site_x_minus];

				int electrode_number = material.electrode_material_interface[electrode_material_interface_number].get_electrode_number();
				double work_function = material.get_work_function(electrode_number);

				b = -material.get_surface_recombination_velocity_hole(electrode_material_interface_number) / concentration.spacing_x * material.get_hole_trans_DOS(site) * std::min(exp(work_function - material.get_hole_trans_energy(site)), 1.0);
			}
		}
		else{

			if ((material.is_interface(site) && material.is_interface(site_x_minus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_minus))){
				int pair_number = material.get_interface_pair(site, site_x_minus);
				a[0] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site_x_minus]), 1.0);
				a[1] += material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_minus) + potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0) * concentration.data[site_x_minus];
			}
			else{
				a[0] += -material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[1] += material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus)) * concentration.data[site_x_minus];
			}
			if ((material.is_interface(site) && material.is_interface(site_x_plus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_plus))){
				int pair_number = material.get_interface_pair(site, site_x_plus);
				a[0] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0);
				a[2] += material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_x * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0) * concentration.data[site_x_plus];
			}
			else{
				a[0] += -material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[2] += material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus)) * concentration.data[site_x_plus];
			}

			a[0] += rate_coef.data[site];
		}
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

			if ((material.is_interface(site) && material.is_interface(site_y_minus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_minus))){
				int pair_number = material.get_interface_pair(site, site_y_minus);
				a[0] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site_y_minus]), 1.0);
				a[3] += material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_minus) + potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0) * concentration.data[site_y_minus];
			}
			else{
				a[0] += -material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[3] += material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus)) * concentration.data[site_y_minus];
			}

			if ((material.is_interface(site) && material.is_interface(site_y_plus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus))){
				int pair_number = material.get_interface_pair(site, site_y_plus);
				a[0] += -material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0);
				a[4] += material.get_interface_hole_transfer_velocity(pair_number) / concentration.spacing_y * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0) * concentration.data[site_y_plus];
			}
			else{
				a[0] += -material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site));
				a[4] += material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus)) * concentration.data[site_y_plus];
			}

		}

		b += -rate.data[site];



		// Now, the next iterative step for the concentrationentration can be calculated.
		concentration.data[site] = (1.0 - c_relaxation[iteration_stage]) * concentration.data[site] + c_relaxation[iteration_stage] * (b - (a[1] + a[2] + a[3] + a[4])) / a[0];
		if (concentration.data[site] > material.get_hole_trans_DOS(site))
			concentration.data[site] = material.get_hole_trans_DOS(site);

		if (concentration.data[site] < 0.0)
			concentration.data[site] = 0.0;
	}

	// Finally, convergence is checked.
	if (fabs(concentration.data[site] - previous) >= c_convergence[iteration_stage] * concentration.data[site]){
		flag_has_converged = false;
	}


}

void PositiveMobileCharge::calculate_concentration(int i, int j, Morphology &material, const Potential &potential, const double previous_concentration, double time_step){

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

		a[0] = -1.0 / pow(concentration.spacing_x, 2.0) * (material.get_hole_mobility(site, site_x_minus) * calc::bernou(potential.electrical.data[site_x_minus] + material.get_hole_trans_energy(site_x_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
			+ material.get_hole_mobility(site, site_x_plus) * calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))) + rate_coef.data[site];

		a[1] = material.get_hole_mobility(site, site_x_minus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_minus] - material.get_hole_trans_energy(site_x_minus)) * concentration.data[site_x_minus];

		a[2] = material.get_hole_mobility(site, site_x_plus) / pow(concentration.spacing_x, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus)) * concentration.data[site_x_plus];

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


			a[0] -= 1.0 / pow(concentration.spacing_y, 2.0) * (material.get_hole_mobility(site, site_y_minus) * calc::bernou(potential.electrical.data[site_y_minus] + material.get_hole_trans_energy(site_y_minus) - potential.electrical.data[site] - material.get_hole_trans_energy(site))
				+ material.get_hole_mobility(site, site_y_plus) * calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site)));

			a[3] = material.get_hole_mobility(site, site_y_minus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_minus] - material.get_hole_trans_energy(site_y_minus)) * concentration.data[site_y_minus];

			a[4] = material.get_hole_mobility(site, site_y_plus) / pow(concentration.spacing_y, 2.0) * calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus)) * concentration.data[site_y_plus];

		}

		a[0] += -1.0 / time_step;

		b = -rate.data[site] - previous_concentration / time_step;


		// Now, the next iterative step for the concentrationentration can be calculated.
		concentration.data[site] = (1.0 - c_relaxation[iteration_stage]) * concentration.data[site] + c_relaxation[iteration_stage] * (b - (a[1] + a[2] + a[3] + a[4])) / a[0];
		if (concentration.data[site] > material.get_hole_trans_DOS(site))
			concentration.data[site] = material.get_hole_trans_DOS(site);

		if (concentration.data[site] < 0.0)
			concentration.data[site] = 0.0;
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

			if (i < concentration.points_x - 1 && !material.is_electrode(site) && !material.is_electrode(site_x_plus)){

				if ((material.is_interface(site) && material.is_interface(site_x_plus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_x_plus))){
					int pair_number = material.get_interface_pair(site, site_x_plus);
					current_x.data[site] = material.get_interface_hole_transfer_velocity(pair_number) * (concentration.data[site] * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site_x_plus]), 1.0)
						- concentration.data[site_x_plus] * std::min(exp(material.get_hole_trans_energy(site_x_plus) + potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0));
				}
				else{
					current_x.data[site] = -material.get_hole_mobility(site, site_x_plus) / current_x.spacing_x * (calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_x_plus] - material.get_hole_trans_energy(site_x_plus)) * concentration.data[site_x_plus]
						- calc::bernou(potential.electrical.data[site_x_plus] + material.get_hole_trans_energy(site_x_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site)) * concentration.data[site]);
				}
			}
			if (current_y.points_y > 1){
				int site_y_plus;

				if (j != current_y.points_y-1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - concentration.points_y;

				if (!material.is_electrode(site) && !material.is_electrode(site_y_plus)){

					if ((material.is_interface(site) && material.is_interface(site_y_plus)) && (material.get_lattice_number(site) != material.get_lattice_number(site_y_plus))){
						int pair_number = material.get_interface_pair(site, site_y_plus);
						current_y.data[site] = material.get_interface_hole_transfer_velocity(pair_number) * (concentration.data[site] * std::min(exp(material.get_hole_trans_energy(site) + potential.electrical.data[site] - material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site_y_plus]), 1.0)
							- concentration.data[site_y_plus] * std::min(exp(material.get_hole_trans_energy(site_y_plus) + potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site) - potential.electrical.data[site]), 1.0));
					}
					else{
						current_y.data[site] = -material.get_hole_mobility(site, site_y_plus) / current_y.spacing_y * (calc::bernou(potential.electrical.data[site] + material.get_hole_trans_energy(site) - potential.electrical.data[site_y_plus] - material.get_hole_trans_energy(site_y_plus)) * concentration.data[site_y_plus]
							- calc::bernou(potential.electrical.data[site_y_plus] + material.get_hole_trans_energy(site_y_plus) - potential.electrical.data[site] - material.get_hole_trans_energy(site)) * concentration.data[site]);
					}
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

	if (device_properties.get_neumann_zero_boundary()){
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
		concentration.data[site] = std::min(device_properties.get_hole_trans_DOS(site) * exp(work_function - device_properties.get_hole_trans_energy(site)), device_properties.get_hole_trans_DOS(site));
		if (!device_properties.surface_recombination_activated(electrode_material_interface_number))
		concentration.calculate_this[site] = false;
		// This is the thermionic expression, which is simply a constant.
	}
	else if (device_properties.get_boundary_condition(electrode_material_interface_number) == 1){
		concentration.data[site] = device_properties.get_hole_trans_DOS(site) * std::min(exp(work_function + device_properties.get_electrode_potential(electrode_number)
			- device_properties.get_hole_trans_energy(site) - electric_potential.data[site]),1.0);
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
					if (device_properties.dopants.data[site] <= -1E23 / concentration.normalization_coef){
						concentration.data[site] = device_properties.dopants.data[site];
					}
					else{
						double N_an = device_properties.get_hole_trans_DOS(site) * exp(device_properties.get_work_function(0) - device_properties.get_hole_trans_energy(site));
						double N_cat = device_properties.get_hole_trans_DOS(site) * exp(device_properties.get_work_function(1) - device_properties.get_hole_trans_energy(site));

						concentration.data[site] = N_an * pow((N_cat / N_an), ((double)i) / ((double)concentration.points_x - 1.0));
					}
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
						concentration.data[site] = device_properties.get_hole_trans_DOS(site) * std::min(exp(device_properties.get_work_function(1) - device_properties.get_hole_trans_energy(site) - electric_potential.data[site]), 1.0);
					}
					else if (device_properties.get_material_number(device_properties.get_lattice_number(site)) == 1){
						concentration.data[site] = device_properties.get_hole_trans_DOS(site) * std::min(exp(device_properties.get_work_function(1) - device_properties.get_hole_trans_energy(site) - electric_potential.data[site]), 1.0);
					}
				}
			}
		}
	}
	else if (initial_guess == 3) {
		for (int i = 0; i < concentration.points_x; i++) {
			for (int j = 0; j < concentration.points_y; j++) {
				int site = i * concentration.points_y + j;
				if (concentration.calculate_this[site] == true && !device_properties.is_electrode_interface(site)) {
					concentration.data[site] = 0.0;
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

	if (material.get_neumann_zero_boundary()){
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
	
		if (concentration.data[site] > material.get_max_ion_concentration(site))
			concentration.data[site] = material.get_max_ion_concentration(site);
		if (concentration.data[site] < 0.0)
			concentration.data[site] = 0.0;

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