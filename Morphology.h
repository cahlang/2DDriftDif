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

	double active_layer_lenght_x, active_layer_lenght_y;			// Stores the lenght of the active layer in x-direction (perpendicular to electrode surfaces) and y-direction (parallell to electrode surfaces)
																	// Add modes and descriptions for different types of generations and recombinaition here.
	int initial_guess;
	bool convergence_boundary_condition;
	int convergence_boundary_condition_lattice_number;
	void map_cbc();
	std::vector<bool> is_cbc_lattice;

	bool ion_transport = false;

	std::vector<bool> recombination_mode = std::vector<bool>(3);
	// 0: Recombination occurs between charges at the same point. 
	// 1: Recombination between charges at opposite sides of an interface occurs.
	// 2: Recombination between a free charge and a charge in an mid-gap state.

	std::vector<bool> generation_mode = std::vector<bool>(3);
	// 0: Generation occurs everywhere in the material (The rate of generation is set separately for each material). 
	// 1: Generation occurs across interfaces, holes are generated in the material with higher HOMO and electrons in that with lower LUMO (The rate of generation is set separately for each interface).

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

	// Currently does not do what it is supposed to. Under construction!
	int get_boundary_condition(int interface_number);

	// Get values of the listed parameters, either at the given site or at the interface between the two given sites or for a given pair number.
	// The required input has been selected based on convenience when dealing with the equations. However, the pair number is always known for a pair of sites.
	double get_hole_mobility(int site);
	double get_hole_mobility(int site_one, int site_two);
	double get_electron_mobility(int site);
	double get_electron_mobility(int site_one, int site_two);
	double get_relative_permittivity(int site);
	double get_interface_relative_permittivity(int pair_number);
	double get_generation_rate(int site);
	double get_interface_generation_rate(int pair_number);
	void get_rate_generation(int site, double generation_rate[2], const PositionDependentParameter electrical_potential);
	double get_DOS(int site);
	double get_interface_DOS(int pair_number);
	double get_LUMO(int site);
	double get_interface_LUMO(int pair_number);
	double get_HOMO(int site);
	double get_interface_HOMO(int pair_number);
	double get_MG_minus(int site);
	double get_MG_plus(int site);
	double get_MG_DOS(int site);
	double get_active_layer_lenght_x();

	double get_bulk_reduced_recombination_coef(int site);
	double get_bulk_hole_capture_coef(int site);
	double get_bulk_electron_capture_coef(int site);

	double get_interface_reduced_recombination_coef(int pair_number);
	double get_interface_hole_capture_coef(int pair_number);
	double get_interface_electron_capture_coef(int pair_number);

	double get_interface_trap_reduced_recombination_coef(int pair_number);
	double get_interface_trap_hole_capture_coef(int pair_number);
	double get_interface_trap_electron_capture_coef(int pair_number);

	double get_positive_ion_mobility(int site);
	double get_negative_ion_mobility(int site);
	double get_positive_ion_mobility(int site_one, int site_two);
	double get_negative_ion_mobility(int site_one, int site_two);

	double get_negative_ion_eq_conc(int site);
	double get_positive_ion_eq_conc(int site);

	bool negative_ion_transport(int site);
	bool positive_ion_transport(int site);

	int get_material_number(int material_lattice_number);

	int get_electrode_interface_number(int site, int electrode_number);
	int get_electrode_number(int electrode_lattice_number);
	void set_electrode_potential(int electrode_number, double potential);
	double get_electrode_potential(int electrode_number);
	double get_electrode_electron_concentration(int interface_number);
	double get_electrode_hole_concentration(int interface_number);
	double get_work_function(int electrode_number);

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
	bool get_convergence_boundary_condition();
	bool ion_transport_activated();

};