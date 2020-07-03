// The settings class is used to store input needed by the program.
// Initially, all input is stored in the variables of this class. Once the appropriate class (e.g. Potential) is created, relevant parameters are copied to it.

#pragma once

#include <vector>
#include <array>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;


// Used for storing variables taking on decimal values at each lattice point.
struct PositionDependentParameter{

	void initialize(double norm);
	
	// The number for individual sites should (or will, in the case of output) be separated by blankspaces in the y-direction and new rows in the x-direction.
	// The electrodes are located one point above x=0 and one belove x = sites_x-1 (Or in other words, above the first row in the file and belove the last).

	// Reads the data from a file with a given name (including file extension).
	void read_from_file(std::string file_name);
	// Alternatively, data can be read using a file name template and a number, making it easy to read from a new file for each calculation. 
	// For example, if the data is stored in "potential_0.dat", "potential_1.dat" and so on, the template is "potential" and the number 0, 1...
	void read_from_file(std::string file_name_template, int file_number);
	// Outputs the data to a new file with the given name (overwriting it if it exists).
	void output_data(std::string file_name);
	// Outputs the data to a new file with the given name template and number, see example above, (overwriting it if it exists).
	void output_data(std::string file_name_template, int file_number);

	// The number of points and the distance between points in the specified direction.
	static int points_x, points_y;
	static double spacing_x, spacing_y;

	// Offset from the regular grid, starting at x=0.0, y=0.0 and spacing_x/y between points. For example the current is calculated between points, 
	// and is therefore offset by 0.5*spacing_x and 0.5*spacing_y in the respective directions.
	// If offset[0] == true, the members of data are assumed to be offset by 0.5*spacing_x in the x-direction and if offset[1] == true, they are offset by 0.5*spacing_y in the y-direction.
	std::array<bool, 2> offset{}; // Change to <double, 3> for the 3d case.

	// Stores the values of the parameter.
	std::vector<double> data;

	std::vector<bool> calculate_this;

	// The normalization coefficient (if any) for the parameter.
	double normalization_coef;

};


// Stores a number of material parameters.
struct Material{

	int lattice_number;

	double generation_rate;
	double mobility_electron;
	double mobility_hole; 
	double hole_trans_energy;
	double electron_trans_energy;
	double relative_permittivity;
	double DOS;
	double electron_trans_DOS, hole_trans_DOS;

	double MG_state_DOS;
	double MG_plus_level;
	double MG_minus_level;

	double electron_trap_energy, electron_trap_DOS;
	double hole_trap_energy, hole_trap_DOS;

	double bulk_electron_trap_capture_coef_electron, bulk_electron_trap_capture_coef_hole;
	double bulk_hole_trap_capture_coef_electron, bulk_hole_trap_capture_coef_hole;

	double bulk_reduced_recombination_coef;
	double bulk_bimolecular_recombination_coef;
	double bulk_hole_capture_coef;
	double bulk_electron_capture_coef;

	bool negative_ion_transport;
	bool positive_ion_transport;

	double mobility_negative_ion;
	double mobility_positive_ion;

	double negative_ion_eq_conc;
	double positive_ion_eq_conc;

	double max_ion_concentration;

	double generation_rate_initial;
	double mobility_electron_zero;
	double mobility_hole_zero;

	double field_dependent_mobility_coef;

};


// Stores some of the same material parameters and some specific to interfaces. 
// These can be accessed instead of the normal bulk parameters, for example a different mobility can be specified for movement across the interface.
struct MaterialInterface : Material {

	std::vector<std::array<int, 2>> interface_pairs;

	std::array<int, 2> material_numbers;
	double interface_reduced_recombination_coef;
	double interface_hole_capture_coef;
	double interface_electron_capture_coef;

	double interface_hole_trap_energy;
	double interface_electron_trap_energy;
	double interface_hole_trap_DOS;
	double interface_electron_trap_DOS;

	int hole_trap_material_number;
	int electron_trap_material_number;

	double interface_trap_reduced_recombination_coef;
	double interface_hole_trap_hole_capture_coef;
	double interface_hole_trap_electron_capture_coef;
	double interface_electron_trap_hole_capture_coef;
	double interface_electron_trap_electron_capture_coef;

	double interface_bimolecular_recombination_coef_1, interface_bimolecular_recombination_coef_2;

	double interface_electron_transfer_velocity, interface_hole_transfer_velocity;

};

// Stores some electrode parameters
struct Electrode{

	int lattice_number;

	double work_function;

	double potential;	// This is the applied potential in the current step of the selected measurement.

	double get_work_function(){
		return work_function;
	}

};

// Similarly to material interfaces, separate parameters can be selected for electrode interfaces.
// This structure also keeps track of which sites are located at each pair of one material and one interface. 
// For example, the sites 3, 4, 5 and 6 may be located at the anode/acceptor interface.
struct ElectrodeMaterialInterface : MaterialInterface {

	int electrode_number;
	int material_number;

	std::vector<int> sites;
	void add_electrode_interface_site(int site){
		sites.push_back(site);
		return;
	}
	std::vector<double> spacing;
	void add_electrode_interface_site_spacing(int d){
		spacing.push_back(d);
		return;
	}

	int boundary_condition;

	bool surface_recombination = false;
	double surface_recombination_velocity_electron, surface_recombination_velocity_hole;

	// Currently not in use.
	double eq_hole_concentration;
	double eq_electron_concentration;

	int get_electrode_number(){
		return electrode_number;
	}

};

// Stores parameters related to the selected measurement.
struct Experiment{

	bool experiment_iv, electrostatics, experiment_CELIV, experiment_PhotoCELIV, experiment_transient_photocurrent;
	bool output_time, output_potential, output_current;

	int data_points;
	int current_data_point = 0;

	double potential_step;
	std::vector<double> electrode_potentials;
	std::vector<double> anode_potential;
	std::vector<double> cathode_potential;
	double time_step;
	std::vector<double> time;


	bool iterate_forward = false;

	void initialize(pt::ptree settings);
	void output_data(); 

	void set_current_data_point(int number);
	int get_current_data_point();
	void next_data_point();

	double voltage_rise_speed, voltage_offset, pulse_length;
	double light_delay;
	bool output_full_data;

};


// Stores the current that would be observed in the measurement and the RMS error.
struct MeasuredCurrent{

	double outer_circuit = 0.0;	// The value of the current measured in an outer circuit.
	double rms_error = 0.0;		// Root mean square error value.

	double normalization_coef;	// Copied to here for convenience. 

};

// Stores information related to pair sites. Each pair of adjacent sites of different materials (excluding electrodes) are given a pair number. 
// This is useful for performing some calculation for all interfaces involving variables taking on different values at each site, regardless of the materials at that specific interface.
struct InterfacePair{

	int interface_material;		// The number of the corrsponding material interface properties.
	std::array<int, 2> site{};	// The site numbers of the two sites making up the pair.
	double spacing;				// The spacing between them, this is equal to the grid spacing unless "empty" spaces have been included.

	std::array<int, 2> get_pair_sites(){	// Gets the pair sites for the given pair number.
		return site;
	}

};

