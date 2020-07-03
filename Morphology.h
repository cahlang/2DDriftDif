#include "DataTypes.h"
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

class Morphology{

	std::vector<int> lattice;		// Stores the morphology of the active layer, in the form of material numbers for each site.
	std::vector<bool> is_electrode_lattice;
	std::vector<bool> is_electrode_interface_lattice;
	int points_x, points_y;			// The number of lattice points in each direction
	double spacing_x, spacing_y;	// The distance between lattice points
	std::string morphology_file;	// The name of the morphology file, consisting of material numbers with blankspaces separating numbers in the y-direction and new rows in the x-direction.
	void read_morphology();			// Reads the morphology file
	std::string dopants_file;		// The name of the dopants file, any distribution of dopants can be given to the program
	void read_dopants();			// Reads the dopants file and stores it in "dopants", declared below (under public)

	double active_layer_length_x, active_layer_length_y;			// Stores the length of the active layer in x-direction (perpendicular to electrode surfaces) and y-direction (parallell to electrode surfaces)
																	// Add modes and descriptions for different types of generations and recombinaition here.
	double norm_potential, norm_concentration, norm_current, norm_rate;

	// 

	int initial_guess;


	bool neumann_zero_boundary;
	int neumann_zero_boundary_lattice_number;
	void map_neumann_zero_boundary();
	std::vector<bool> is_neumann_zero_boundary;

	bool light_on_flag = true;
	bool ion_transport = false;
	bool effective_temperature_activated = false;
	bool field_dependent_mobility_activated = false;
	PositionDependentParameter field_dependent_mobility_electron, field_dependent_mobility_hole;

	std::vector<bool> recombination_mode = std::vector<bool>(4);
	// 0: Recombination occurs between charges at the same point. 
	// 1: Recombination between charges at opposite sides of an interface occurs.
	// 2: Recombination between a free charge and a charge in an mid-gap state.

	std::vector<bool> generation_mode = std::vector<bool>(3);
	// 0: Generation occurs everywhere in the material (The rate of generation is set separately for each material). 
	// 1: Generation occurs across interfaces, holes are generated in the material with higher hole_trans_energy and electrons in that with lower electron_trans_energy (The rate of generation is set separately for each interface).

	int generation_profile_type = 0;
	std::vector<double> position_dependent_generation_rate;
	void create_generation_profile(pt::ptree & settings);

public:

	// The mid-gap states at interfaces are called CT states in the simulation, even though no real charge transfer mechanics are included in this simulation.
	PositionDependentParameter MG_electron_concentration;		// Stores the population of negatively charged mid-gap states
	PositionDependentParameter MG_hole_concentration;			// Stores the population of positively charged mid-gap states

	std::vector<struct Electrode> electrode;					// Contains properties for each electrode (typically two)
																// Typically 0 corresponds to the anode and 1 to the cathode

	std::vector<struct Material> material;						// Contains the material properties for each material number
	std::vector<struct MaterialInterface> material_interface;	// Contains the properties of interfaces for each pair of materials
	std::vector<struct ElectrodeMaterialInterface> electrode_material_interface;

																
	std::vector<std::vector <int>> interfaces;					// interfaces[0] corresponds to the negative x-direction, interface[1] to the positive x-direction and so on.
	std::vector<struct InterfacePair> interface_pair_data;		// Contains information about all pairs of adjacent sites with different morphologies.

	std::vector< std::array<int, 2> > cbc_interface;

	PositionDependentParameter dopants;							// Stores the concentration of stationary charge

	// Initializes vectors and reads input parameters
	void initialize(pt::ptree &settings);
	void read_material_data(pt::ptree & settings);
	void read_material_interface_data(pt::ptree &settings);
	void read_electrode_data(pt::ptree &settings);


	void identify_interfaces();							// Automatically detects and stores information on where interfaces are located.
	void determine_interface_pairs();					// Finds all pairs of adjacent sites with different morphology.

	void map_electrodes();

	void set_boundary_conditions(const Experiment &measurement);

	// Calculates the recombination and generation rates for the given values of potential and carrier concentrations. 
	// No field-dependence is included so far, the potential is only needed for calculating SRH recombination via CT states.
	void calculate_generation_rate(const PositionDependentParameter &electron_concentration, PositionDependentParameter &electron_rate, const PositionDependentParameter &hole_concentration, PositionDependentParameter &hole_rate, const PositionDependentParameter &potential, PositionDependentParameter &net_rate);
	void calculate_recombination_rates(const PositionDependentParameter &electron_concentration, PositionDependentParameter &electron_rate, PositionDependentParameter &electron_rate_coef, 
		const PositionDependentParameter &hole_concentration, PositionDependentParameter &hole_rate, PositionDependentParameter &hole_rate_coef, const PositionDependentParameter &potential, PositionDependentParameter &net_rate);

	// Calculates the mid-gap state populations for the given values of potential and carrier concentrations.
	void calculate_MG_state_population(const PositionDependentParameter &electric_potential, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration);
	void set_MG_initial_guess();


	// Introduces temperature-dependence for the mobility and generation rate and adjusts these by an effective temperature
	void effective_temperature();
	bool is_effective_temperature();
	bool is_field_dependent_mobility();
	void calculate_field_dependent_mobility(PositionDependentParameter electric_potential);
	void set_initial_field_dependent_mobility();
	bool initial_values_saved = false;

	// Currently does not do what it is supposed to. Under construction!
	int get_boundary_condition(int interface_number);

	// Get values of the listed parameters, either at the given site or at the interface between the two given sites or for a given pair number.
	// The required input has been selected based on convenience when dealing with the equations. However, the pair number is always known for a pair of sites.
	double get_hole_mobility(int site);
	double get_hole_mobility(int site_one, int site_two);
	double get_hole_mobility_initial(int site_one, int site_two);
	double get_electron_mobility(int site);
	double get_electron_mobility(int site_one, int site_two);
	double get_electron_mobility_initial(int site_one, int site_two);
	double get_interface_hole_transfer_velocity(int pair_number);
	double get_interface_electron_transfer_velocity(int pair_number);

	double get_relative_permittivity(int site);
	double get_relative_permittivity(int site_one, int site_two);
	double get_interface_relative_permittivity(int pair_number);
	double get_generation_rate(int site);
	double get_interface_generation_rate(int pair_number);
	void get_rate_generation(int site, double generation_rate[2], const PositionDependentParameter electrical_potential);
	double get_DOS(int site);
	double get_interface_DOS(int pair_number);
	double get_electron_trans_energy(int site);
	double get_electron_trans_DOS(int site);
	double get_interface_electron_trans_energy(int pair_number);
	double get_hole_trans_energy(int site);
	double get_hole_trans_DOS(int site);
	double get_interface_hole_trans_energy(int pair_number);
	double get_MG_minus(int site);
	double get_MG_plus(int site);
	double get_MG_DOS(int site);
	double get_active_layer_length_x();

	double get_field_dependent_hole_mobility(int site_one, int site_two, double electric_field);
	double get_field_dependent_electron_mobility(int site_one, int site_two, double electric_field);
	double get_field_dependent_mobility_coef(int site);
	double get_interface_field_dependent_mobility_coef(int pair_number);
	double get_field_dependent_mobility_coef(int site_one, int site_two);

	double get_bulk_reduced_recombination_coef(int site);
	double get_bulk_hole_capture_coef(int site);
	double get_bulk_electron_capture_coef(int site);
	double get_bulk_hole_trap_capture_coef_electron(int site);
	double get_bulk_hole_trap_capture_coef_hole(int site);
	double get_bulk_electron_trap_capture_coef_electron(int site);
	double get_bulk_electron_trap_capture_coef_hole(int site);

	double get_bulk_bimolecular_recombination_coef(int site);

	double get_interface_reduced_recombination_coef(int pair_number);
	double get_interface_hole_capture_coef(int pair_number);
	double get_interface_electron_capture_coef(int pair_number);

	double get_interface_bimolecular_recombination_coef_1(int pair_number);
	double get_interface_bimolecular_recombination_coef_2(int pair_number);

	double get_interface_trap_reduced_recombination_coef(int pair_number);
	double get_interface_hole_trap_hole_capture_coef(int pair_number);
	double get_interface_hole_trap_electron_capture_coef(int pair_number);
	double get_interface_electron_trap_hole_capture_coef(int pair_number);
	double get_interface_electron_trap_electron_capture_coef(int pair_number);

	int get_hole_trap_material_number(int pair_number);
	int get_electron_trap_material_number(int pair_number);

	double get_electron_trap_energy(int site);
	double get_electron_trap_DOS(int site);
	double get_hole_trap_energy(int site);
	double get_hole_trap_DOS(int site);

	double get_interface_hole_trap_energy(int pair_number);
	double get_interface_electron_trap_energy(int pair_number);
	double get_interface_hole_trap_DOS(int pair_number);
	double get_interface_electron_trap_DOS(int pair_number);

	double get_electrode_interface_hole_trap_energy(int interface_number);
	double get_electrode_interface_electron_trap_energy(int interface_number);
	double get_electrode_interface_hole_trap_DOS(int interface_number);
	double get_electrode_interface_electron_trap_DOS(int interface_number);

	double get_positive_ion_mobility(int site);
	double get_negative_ion_mobility(int site);
	double get_positive_ion_mobility(int site_one, int site_two);
	double get_negative_ion_mobility(int site_one, int site_two);

	double get_negative_ion_eq_conc(int site);
	double get_positive_ion_eq_conc(int site);

	bool negative_ion_transport(int site);
	bool positive_ion_transport(int site);

	double get_max_ion_concentration(int site);

	int get_material_number(int material_lattice_number);

	int get_electrode_interface_number(int site, int electrode_number);
	int get_electrode_number(int electrode_lattice_number);
	void set_electrode_potential(int electrode_number, double potential);
	double get_electrode_potential(int electrode_number);
	double get_electrode_electron_concentration(int interface_number);
	double get_electrode_hole_concentration(int interface_number);
	double get_work_function(int electrode_number);
	double get_surface_recombination_velocity_electron(int interface_number);
	double get_surface_recombination_velocity_hole(int interface_number);

	double get_electrode_trap_reduced_recombination_coef(int electrode_number);

	int get_pair(int site, int n);
	int get_interface_pair(int site_one, int site_two);

	int get_lattice_number(int site);
	std::array<int,2> get_material_interface_pair(int material_interface_number);
	std::array<int, 2> get_electrode_material_interface_pair(int electrode_material_interface_number);

	ElectrodeMaterialInterface get_electrode_material_interface(int electrode_material_interface_number);
	int get_electrode_material_interface_number(int electrode_number, int site);
	int get_electrode_material_interface_number(int stie);
	// Checks if the given site is at an interface.
	bool is_interface(int site);
	bool is_electrode(int site);
	bool is_electrode_interface(int site);
	bool is_cbc(int site);
	bool surface_recombination_activated(int interface_number);
	bool get_neumann_zero_boundary();
	bool ion_transport_activated();

	bool is_light_on();
	void light_off();
	void light_on();

};