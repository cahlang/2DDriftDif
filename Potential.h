#pragma once

#include "DataTypes.h"
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

class Morphology;

class Potential{

	// Stores parameters related to the iterative method
	std::vector <double> c_relaxation;
	std::vector<double> c_convergence;
	int iteration_stage = 0;

	bool flag_has_converged = false;
	double sweep_direction = 1;

	int initial_guess;

public:

	// Stores the electrical potential.
	PositionDependentParameter electrical;
	// Stores the chemical potential, not in use at the moment. Instead, a constant value is assumed for each material.
	PositionDependentParameter chemical;
	// Stores the electrochemical potential for mobile holes and electrons.
	PositionDependentParameter electrochemical_electron;
	PositionDependentParameter electrochemical_hole;

	PositionDependentParameter displacement_current;

	bool has_converged();

	void initialize(pt::ptree &settings);
	void solve(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration); // Calculates the next iterative solution for the electrical potential.
	void solve_inverse(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration);
	void solve_one_d(Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration);
	void set_boundary_conditions(Morphology device_properties);
	void set_initial_guess(Morphology device_properties, std::string file_name_template, int file_number);

	void calculate_electrochemical_potential(Morphology material, PositionDependentParameter electron_concentration, PositionDependentParameter hole_concentration);

	void change_sweep_direction();

	void calculate_electrical_potential(int i, int j, Morphology &material, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration, const PositionDependentParameter &negative_ion_concentration, const PositionDependentParameter &positive_ion_concentration);

	void calculate_displacement_current(Morphology &material, const PositionDependentParameter previous_potential, const double time_step);

	void next_iterative_stage();
	bool final_iterative_stage_reached();
	void reset_iterative_stage();

	double get_electric(int site){
		return electrical.data[site];
	}

};