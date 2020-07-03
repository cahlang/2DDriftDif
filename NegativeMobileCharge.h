#pragma once

#include "DataTypes.h"
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

class Potential;
class Morphology;

class NegativeMobileCharge{

	// Stores parameters related to the iterative method
	std::vector<double> c_relaxation;
	std::vector<double> c_convergence;
	int iteration_stage = 0;

	bool flag_has_converged = false;
	double sweep_direction = 1;

	int initial_guess;

public:

	// Stores the concentration of the negative mobile charge (Typically electrons).
	PositionDependentParameter concentration;
	// Stores the current density in the x- and y-direction, respectively.
	PositionDependentParameter current_x;
	PositionDependentParameter current_y;

	PositionDependentParameter drift_current_x;
	PositionDependentParameter drift_current_y;

	PositionDependentParameter diff_current_x;
	PositionDependentParameter diff_current_y;

	// Stores the net rate generation/recombination NOT dependent on this concentration of charge.
	PositionDependentParameter rate;
	// Stores the net rate generation/recombination dependent on this concentration of charge, divided by this concentration.
	// This distinction is nescessary due to numerical reasons.
	PositionDependentParameter rate_coef;

	// Checks if the calculations of this concentraiton has reached the convergence criteria.
	bool has_converged(){
		return flag_has_converged;
	}

	// Initializes all vectors.
	void initialize(pt::ptree &settings); 

	// Calculates the concentrations of this type of charge for the given potential and current data in rate and rate_coef.
	void solve(Morphology &material, const Potential &potential); 
	void solve_inverse(Morphology &material, const Potential &potential);
	void solve_one_d(Morphology &material, const Potential &potential);
	void solve(Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentraiton, double time_step);
	void solve_one_d(Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentration, double time_step);
	// Substeps of solve.
	void calculate_concentration(int i, int j, Morphology &material, const Potential &potential);
	void calculate_concentration(int i, int j, Morphology &material, const Potential &potential, const double previous_time_concentration, double time_step);
	void calculate_current(Morphology &material, const Potential potential);

	void calculate_concentration_one_d(int i, int j, Morphology &material, const Potential &potential);


	void ion_solve(Morphology &material, const Potential &potential);
	void ion_solve(Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentraiton, double time_step);

	void calculate_ion_concentration(int i, int j, Morphology &material, const Potential &potential);
	void calculate_ion_concentration(int i, int j, Morphology &material, const Potential &potential, const PositionDependentParameter &previous_time_concentraiton, double time_step);
	void calculate_ion_current(Morphology &material, const Potential &potential);

	// Updates the boundary conditions for the concentration for the current calculation.
	void set_boundary_conditions(Morphology &material, const PositionDependentParameter &electric_potential);
	void get_boundary_condition(int site, int electrode_material_interface_number, Morphology &material, const PositionDependentParameter &electric_potential);
	void set_initial_guess(Morphology &device_properties, const PositionDependentParameter &electric_potential, std::string file_name_prototype, int file_number);

	void set_ion_boundary_conditions(Morphology &material);
	void set_ion_initial_guess(Morphology &device_properties, const PositionDependentParameter &electric_potential, std::string file_name_prototype, int file_number);

	// Reads input parameters and settings.
	void read_settings();

	// Functions used by the iterative method
	void change_sweep_direction();
	void next_iterative_stage();
	bool final_iterative_stage_reached();
	void reset_iterative_stage();

};