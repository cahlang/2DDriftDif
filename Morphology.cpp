#define BULK_RECOMBINATION 0
#define INTERFACE_RECOMBINATION 1
#define TRAP_ASSISTED_RECOMBINATION_MG 2

#define BULK_GENERATION 0
#define INTERFACE_GENERATION 1

#define SPATIAL_MATERIAL_GAP 0
#define NOT_ALREADY_INTERFACE -1
#define INTERFACES_ELECTRODE_NUMBER -3
#define INTERFACES_CBC_NUMBER -97

#include "Morphology.h"
#include "DataTypes.h"
#include "Constants.h"
#include "Misc.h"

#include <ctime> // Gets a seed for the pseudorandom number generator
#include <cstdlib>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <string>
#include <sstream>
#include <iostream>
namespace pt = boost::property_tree;

#define MAX_NUMBER_MATERIALS 4	// Can be set to anything, but should be kept as small as possible.
#define MAX_NUMBER_MATERIAL_INTERFACES 10

void Morphology::initialize(pt::ptree &settings){

	points_x = settings.get<int>("general.points_x");
	points_y = settings.get<int>("general.points_y", 1);
	active_layer_lenght_x = settings.get<double>("general.lenght_x") / constant::lenght_norm;
	active_layer_lenght_y = settings.get<double>("general.lenght_y") / constant::lenght_norm;
	spacing_x = active_layer_lenght_x / ((double)points_x - 1.0);
	spacing_y = active_layer_lenght_y / ((double)points_y - 1.0);
	morphology_file = settings.get<std::string>("general.morphology_file");
	dopants_file = settings.get<std::string>("general.dopants_file");

	lattice.resize(points_x*points_y, 0);
	interfaces.resize(points_x*points_y, std::vector<int>(4, NOT_ALREADY_INTERFACE));

	read_morphology();

	double norm_potential, norm_concentration, norm_current, norm_rate;
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::lenght_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::lenght_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::lenght_norm);

	dopants.initialize(norm_concentration);
	read_dopants();

	convergence_boundary_condition = settings.get<bool>("general.convergence_boundary_condition", false);
	if (convergence_boundary_condition){
		convergence_boundary_condition_lattice_number = settings.get<int>("general.convergence_boundary_condition_lattice_number");
	}

	MG_electron_concentration.initialize(norm_concentration);
	MG_hole_concentration.initialize(norm_concentration);

	read_material_data(settings);
	read_material_interface_data(settings);
	// identify_interfaces();

	map_electrodes();
	map_cbc();
	determine_interface_pairs();

	int site;

	recombination_mode[0] = settings.get<bool>("general.recombination_mode_0");
	recombination_mode[1] = settings.get<bool>("general.recombination_mode_1");
	recombination_mode[2] = settings.get<bool>("general.recombination_mode_2");

	generation_mode[0] = settings.get<bool>("general.generation_mode_0");
	generation_mode[1] = settings.get<bool>("general.generation_mode_1");

	initial_guess = settings.get<int>("numerics.initial_guess", -1);

	std::cout << "All material parameters read." << std::endl;

	return;

}

void Morphology::read_material_data(pt::ptree & settings){

	double norm_potential, norm_concentration, norm_current, norm_rate;
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::lenght_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::lenght_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::lenght_norm);

	for (int a = 0; a < MAX_NUMBER_MATERIALS; a++){	// Material parameters are stored in the "material" vector. 

		std::ostringstream oss_node;
		oss_node << "material_" << a;
		std::string node_name = oss_node.str();
		boost::optional<pt::ptree&> input_material = settings.get_child_optional(node_name.c_str());

		Material mat;

		if (input_material){

			mat.lattice_number = input_material.get().get<int>("lattice_number");

			mat.mobility_hole = input_material.get().get<double>("mobility_hole") / constant::mobility_norm;
			mat.mobility_electron = input_material.get().get<double>("mobility_electron") / constant::mobility_norm;
			mat.generation_rate = input_material.get().get<double>("generation_rate") / norm_rate;
			mat.HOMO_level = input_material.get().get<double>("HOMO_level") / norm_potential;
			mat.LUMO_level = input_material.get().get<double>("LUMO_level") / norm_potential;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity");
			mat.DOS = input_material.get().get<double>("DOS") / norm_concentration;

			mat.MG_minus_level = input_material.get().get<double>("MG_minus_level", 0.0)/norm_potential;
			mat.MG_plus_level = input_material.get().get<double>("MG_plus_level", 0.0) / norm_potential;
			mat.MG_state_DOS = input_material.get().get<double>("MG_state_DOS", 0.0) / norm_concentration;

			mat.bulk_reduced_recombination_coef = input_material.get().get<double>("bulk_reduced_recombination_coef", 1.0);
			mat.bulk_hole_capture_coef = input_material.get().get<double>("bulk_hole_capture_coef", mat.mobility_hole/mat.relative_permittivity);
			mat.bulk_electron_capture_coef = input_material.get().get<double>("bulk_electron_capture_coef", mat.mobility_electron / mat.relative_permittivity);

			mat.negative_ion_transport = input_material.get().get<bool>("negative_ion_transport", false);
			if (mat.negative_ion_transport){
				mat.mobility_negative_ion = input_material.get().get<double>("mobility_negative_ion")/constant::mobility_norm;
				ion_transport = true;
				mat.negative_ion_eq_conc = input_material.get().get<double>("negative_ion_eq_conc") / norm_concentration;
			}
			mat.positive_ion_transport = input_material.get().get<bool>("positive_ion_transport", false);
			if (mat.positive_ion_transport){
				mat.mobility_positive_ion = input_material.get().get<double>("mobility_positive_ion") / constant::mobility_norm;
				ion_transport = true;
				mat.positive_ion_eq_conc = input_material.get().get<double>("positive_ion_eq_conc") / norm_concentration;
			}

			material.push_back(mat);

		}

	}

	for (int a = 0; a < MAX_NUMBER_MATERIALS; a++){	// Electrode parameters are stored in the "electrode" vector. 

		std::ostringstream oss_node;
		oss_node << "electrode_" << a;
		std::string node_name = oss_node.str();
		boost::optional<pt::ptree&> input_material = settings.get_child_optional(node_name.c_str());

		Electrode mat;

		if (input_material){

			mat.lattice_number = input_material.get().get<int>("lattice_number");

			mat.work_function = input_material.get().get<double>("work_function") / norm_potential;
			// mat.surface_recombination_rate_electron = input_material.get().get<double>("surface_recombination_rate_electron");
			// mat.surface_recombination_rate_hole = input_material.get().get<double>("surface_recombination_rate_hole");
			electrode.push_back(mat);
		
		}

	}

	std::cout << "Material parameters read." << std::endl;

	return;
}

void Morphology::read_material_interface_data(pt::ptree &settings){

	double norm_potential, norm_concentration, norm_current, norm_rate;
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::lenght_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::lenght_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::lenght_norm);

	for (int a = 0; a < MAX_NUMBER_MATERIAL_INTERFACES; a++){	// Material parameters are stored in the "material" vector. 

		std::ostringstream oss_node;
		oss_node << "material_interface_" << a;
		std::string node_name = oss_node.str();
		boost::optional<pt::ptree&> input_material = settings.get_child_optional(node_name.c_str());

		MaterialInterface mat;

		if (input_material){

			mat.mobility_hole = input_material.get().get<double>("mobility_hole") / constant::mobility_norm;
			mat.mobility_electron = input_material.get().get<double>("mobility_electron") / constant::mobility_norm;
			mat.generation_rate = input_material.get().get<double>("generation_rate") / norm_rate;
			mat.HOMO_level = input_material.get().get<double>("HOMO_level") / norm_potential;
			mat.LUMO_level = input_material.get().get<double>("LUMO_level") / norm_potential;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity");
			mat.DOS = input_material.get().get<double>("DOS") / norm_concentration;

			mat.material_numbers[0] = input_material.get().get<int>("material_1");
			mat.material_numbers[1] = input_material.get().get<int>("material_2");

			mat.interface_reduced_recombination_coef = input_material.get().get<double>("interface_reduced_recombination_coef", 1.0);
			mat.interface_hole_capture_coef = input_material.get().get<double>("interface_hole_capture_coef", mat.mobility_hole / mat.relative_permittivity);
			mat.interface_electron_capture_coef = input_material.get().get<double>("interface_electron_capture_coef", mat.mobility_electron / mat.relative_permittivity);

			mat.interface_trap_reduced_recombination_coef = input_material.get().get<double>("interface_trap_reduced_recombination_coef", 1.0);
			mat.interface_trap_hole_capture_coef = input_material.get().get<double>("interface_trap_hole_capture_coef", mat.mobility_hole / mat.relative_permittivity);
			mat.interface_trap_electron_capture_coef = input_material.get().get<double>("interface_trap_electron_capture_coef", mat.mobility_electron / mat.relative_permittivity);

			material_interface.push_back(mat);

		}



	}

	for (int n = 0; n < MAX_NUMBER_MATERIAL_INTERFACES; n++){

		std::ostringstream oss_node;
		oss_node << "electrode_interface_" << n;
		std::string node_name = oss_node.str();
		boost::optional<pt::ptree&> input_material = settings.get_child_optional(node_name.c_str());

		ElectrodeMaterialInterface mat;

		if (input_material){

			mat.material_number = input_material.get().get<int>("material_number");
			mat.electrode_number = input_material.get().get<int>("electrode_number");

			mat.mobility_hole = input_material.get().get<double>("mobility_hole", material[mat.material_number].mobility_hole * constant::mobility_norm) / constant::mobility_norm;
			mat.mobility_electron = input_material.get().get<double>("mobility_electron", material[mat.material_number].mobility_electron * constant::mobility_norm) / constant::mobility_norm;
			mat.generation_rate = input_material.get().get<double>("generation_rate", material[mat.material_number].generation_rate * norm_rate) / norm_rate;
			mat.HOMO_level = input_material.get().get<double>("HOMO_level", material[mat.material_number].HOMO_level * norm_potential) / norm_potential;
			mat.LUMO_level = input_material.get().get<double>("LUMO_level", material[mat.material_number].LUMO_level * norm_potential) / norm_potential;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity", material[mat.material_number].relative_permittivity);
			mat.DOS = input_material.get().get<double>("DOS", material[mat.material_number].DOS * norm_concentration) / norm_concentration;

			mat.interface_reduced_recombination_coef = input_material.get().get<double>("interface_reduced_recombination_coef", 1.0);
			mat.interface_hole_capture_coef = input_material.get().get<double>("interface_hole_capture_coef", mat.mobility_hole / mat.relative_permittivity);
			mat.interface_electron_capture_coef = input_material.get().get<double>("interface_electron_capture_coef", mat.mobility_electron / mat.relative_permittivity);

			mat.interface_trap_reduced_recombination_coef = input_material.get().get<double>("interface_trap_reduced_recombination_coef", 1.0);
			mat.interface_trap_hole_capture_coef = input_material.get().get<double>("interface_trap_hole_capture_coef", mat.mobility_hole / mat.relative_permittivity);
			mat.interface_trap_electron_capture_coef = input_material.get().get<double>("interface_trap_electron_capture_coef", mat.mobility_electron / mat.relative_permittivity);

			mat.boundary_condition = input_material.get().get<int>("boundary_condition");

			electrode_material_interface.push_back(mat);

		}



	}

	std::cout << "Interface parameters read." << std::endl;

	return;

}

void Morphology::map_cbc(){

	is_cbc_lattice.resize(points_x*points_y, false);

	for (int i = 0; i < points_x*points_y; i++){

		for (int j = 0; j < is_cbc_lattice.size(); j++){
			if (lattice[i] == convergence_boundary_condition_lattice_number){
				is_cbc_lattice[i] = true;
			}
		}

	}

	std::cout << "Cbc mapped." << std::endl;

	return;
}

void Morphology::map_electrodes(){

	is_electrode_lattice.resize(points_x*points_y, false);

	for (int i = 0; i < points_x*points_y; i++){

		for (int j = 0; j < electrode.size(); j++){
			if (lattice[i] == electrode[j].lattice_number){
				is_electrode_lattice[i] = true;
			}
		}

	}

	std::cout << "Electrodes mapped." << std::endl;

	return;
}

void Morphology::calculate_recombination_rates(const PositionDependentParameter &electron_concentration, PositionDependentParameter &electron_rate, PositionDependentParameter &electron_rate_coef, 
	const PositionDependentParameter &hole_concentration, PositionDependentParameter &hole_rate, PositionDependentParameter &hole_rate_coef, const PositionDependentParameter &electric_potential, PositionDependentParameter &net_rate){

	if (recombination_mode[BULK_RECOMBINATION]){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){

				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){
					double R = -get_bulk_reduced_recombination_coef(site) * (get_bulk_hole_capture_coef(site) + get_bulk_electron_capture_coef(site));
					electron_rate_coef.data[site] += R * hole_concentration.data[site];
					hole_rate_coef.data[site] += R * electron_concentration.data[site];
				}
			}
		}
	}

	if (recombination_mode[INTERFACE_RECOMBINATION]){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

			double R = -get_interface_reduced_recombination_coef(pair_number) * (get_interface_hole_capture_coef(pair_number) + get_interface_electron_capture_coef(pair_number));

					electron_rate_coef.data[site[0]] += R * hole_concentration.data[site[1]];
					hole_rate_coef.data[site[1]] += R * electron_concentration.data[site[0]];

					electron_rate_coef.data[site[1]] += R * hole_concentration.data[site[0]];
					hole_rate_coef.data[site[0]] += R * electron_concentration.data[site[1]];
			
		}
	}

	if (recombination_mode[TRAP_ASSISTED_RECOMBINATION_MG]){

		for (int n = 0; n < electrode_material_interface.size(); n++){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				int site = electrode_material_interface[n].sites[m];
				int electrode_number = electrode_material_interface[n].get_electrode_number();

				if (electrode_number == 0 && get_MG_plus(site) != 0.0){

					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_MG_plus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_plus(site)
						- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
					double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

					electron_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * hole_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);
					hole_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * electron_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);

				}
				else if (electrode_number == 1 && get_MG_minus(site) != 0.0){
					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_MG_minus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_minus(site)
						- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
					double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

					electron_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * hole_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);
					hole_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * electron_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);

				}

			}
		}

		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

			if (site[0] < points_y*(points_x - 2) && site[0] >= points_y * 2 && site[1] < points_y*(points_x - 2) && site[1] >= points_y * 2){

				// Calculate the fermi level at before contact, we only consider CT- states below this and CT+ above this.
				double fermi_level = (std::min(get_HOMO(site[0]), get_HOMO(site[1])) + std::max(get_LUMO(site[0]), get_LUMO(site[1]))) / 2;


				// We assume that the CT+ state is at site 0 and CT- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
				if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
					std::swap(site[0], site[1]);
				}

				if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){

					double electron_con = electron_concentration.data[site[1]];
					double electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_LUMO(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

					double hole_con = hole_concentration.data[site[0]];
					double hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

					double C_n = get_interface_trap_electron_capture_coef(pair_number);
					double C_p = get_interface_trap_hole_capture_coef(pair_number);

					electron_rate_coef.data[site[1]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * hole_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate_coef.data[site[0]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * electron_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

					electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_LUMO(site[1])) + 1.0), 1.0);

					hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

					electron_rate_coef.data[site[1]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * hole_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate_coef.data[site[0]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[0]) * C_n * C_p * electron_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

				}
				else if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){

				}
				else{
					std::cerr << "Something strange occured while calculating CT recombination rates." << std::endl;
				}
			}

		}
	}

	for (int site = 0; site<points_x*points_y; site++){

		net_rate.data[site] += hole_rate_coef.data[site] * hole_concentration.data[site] - electron_rate_coef.data[site] * electron_concentration.data[site];
	}

	return;
}

void Morphology::calculate_generation_rate(const PositionDependentParameter &electron_concentration, PositionDependentParameter &electron_rate, const PositionDependentParameter &hole_concentration, 
	PositionDependentParameter &hole_rate, const PositionDependentParameter &electric_potential, PositionDependentParameter &net_rate){


	if (recombination_mode[BULK_RECOMBINATION]){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){
				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){
					double intrinsic_rate = get_bulk_reduced_recombination_coef(site) * (get_bulk_hole_capture_coef(site) + get_bulk_electron_capture_coef(site)) * pow(get_DOS(site), 2.0) * exp(get_LUMO(site) - get_HOMO(site));
					electron_rate.data[site] += intrinsic_rate;
					hole_rate.data[site] += intrinsic_rate;
				}
			}
		}
	}

	if (recombination_mode[INTERFACE_RECOMBINATION]){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

			if (get_LUMO(site[0]) < get_LUMO(site[1]) || get_HOMO(site[0]) < get_HOMO(site[1])){

				double intrinsic_rate = get_interface_reduced_recombination_coef(pair_number) * (get_interface_hole_capture_coef(pair_number) + get_interface_electron_capture_coef(pair_number)) *
					pow(get_interface_DOS(pair_number), 2.0) * exp(get_LUMO(site[1]) - get_HOMO(site[0]) - (electric_potential.data[site[0]] - electric_potential.data[site[1]]));

				hole_rate.data[site[0]] += intrinsic_rate;
				electron_rate.data[site[1]] += intrinsic_rate;
			}
			else if (get_LUMO(site[0]) > get_LUMO(site[1]) || get_HOMO(site[0]) > get_HOMO(site[1])){

				double intrinsic_rate = get_interface_reduced_recombination_coef(pair_number) * (get_interface_hole_capture_coef(pair_number) + get_interface_electron_capture_coef(pair_number)) *
					pow(get_interface_DOS(pair_number), 2.0) * exp(get_LUMO(site[0]) - get_HOMO(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]]));

				hole_rate.data[site[1]] += intrinsic_rate;
				electron_rate.data[site[0]] += intrinsic_rate;
			}
		}
	}

	if (recombination_mode[TRAP_ASSISTED_RECOMBINATION_MG]){

		for (int n = 0; n < electrode_material_interface.size(); n++){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				int site = electrode_material_interface[n].sites[m];
				int electrode_number = electrode_material_interface[n].get_electrode_number();

				double intrinsic_concentration_squared = pow(get_DOS(site), 2.0) * exp(get_LUMO(site) - get_HOMO(site));


				if (electrode_number == 0 && get_MG_plus(site) != 0.0){

					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_MG_plus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_plus(site)
						- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
					double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

					electron_rate.data[site] += get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * intrinsic_concentration_squared / (electron_con + electron_con_zero + hole_con + hole_con_zero);
					hole_rate.data[site] += get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * intrinsic_concentration_squared / (electron_con + electron_con_zero + hole_con + hole_con_zero);

				}
				else if (electrode_number == 1 && get_MG_minus(site) != 0.0){
					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_MG_minus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_minus(site)
						- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
					double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

					electron_rate.data[site] += get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * intrinsic_concentration_squared / (electron_con + electron_con_zero + hole_con + hole_con_zero);
					hole_rate.data[site] += get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * intrinsic_concentration_squared / (electron_con + electron_con_zero + hole_con + hole_con_zero);

				}

			}
		}

		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();


			if (site[0] < points_y*(points_x - 2) && site[0] >= points_y * 2 && site[1] < points_y*(points_x - 2) && site[1] >= points_y * 2){

				// Calculate the fermi level at before contact, we only consider CT- states below this and CT+ above this.
				double fermi_level = (std::min(get_HOMO(site[0]), get_HOMO(site[1])) + std::max(get_LUMO(site[0]), get_LUMO(site[1]))) / 2;


				// We assume that the CT+ state is at site 0 and CT- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
				if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
					std::swap(site[0], site[1]);
				}

				if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){

					double electron_con = electron_concentration.data[site[1]];
					double electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_LUMO(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

					double hole_con = hole_concentration.data[site[0]];
					double hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

					double intrinsic_concentration_squared = pow(get_interface_DOS(pair_number), 2.0) * exp(get_LUMO(site[1]) - get_HOMO(site[0]) - (electric_potential.data[site[0]] - electric_potential.data[site[1]]));

					double C_n = get_interface_trap_electron_capture_coef(pair_number);
					double C_p = get_interface_trap_hole_capture_coef(pair_number);

					electron_rate.data[site[1]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate.data[site[0]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[0]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

					electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_LUMO(site[1])) + 1.0), 1.0);

					hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

					electron_rate.data[site[1]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate.data[site[0]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[0]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
				}
				else if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){

				}
				else{
					std::cerr << "Something strange occured while calculating CT recombination rates." << std::endl;
				}
			}

		}
	}

	if (generation_mode[BULK_GENERATION]){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){
				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){
					double bulk_generation_rate = get_generation_rate(site);
					electron_rate.data[site] += bulk_generation_rate;
					hole_rate.data[site] += bulk_generation_rate;
				}
			}
		}
	}

	if (generation_mode[INTERFACE_GENERATION]){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
			if (get_LUMO(site[0]) < get_LUMO(site[1]) || get_HOMO(site[0]) < get_HOMO(site[1])){
				hole_rate.data[site[0]] += get_interface_generation_rate(pair_number);
				electron_rate.data[site[1]] += get_interface_generation_rate(pair_number);
			}
			else if (get_LUMO(site[0]) > get_LUMO(site[1]) || get_HOMO(site[0]) > get_HOMO(site[1])){
				hole_rate.data[site[1]] += get_interface_generation_rate(pair_number);
				electron_rate.data[site[0]] += get_interface_generation_rate(pair_number);
			}
			else{
				std::cerr << "Unable to determine donor/acceptor for pair number" << pair_number << std::endl;
			}
		}
	}

	for (int site = 0; site < points_x*points_y; site++){
		net_rate.data[site] += hole_rate.data[site] - electron_rate.data[site];
	}

	return;

}

void Morphology::set_MG_initial_guess(){

	if (initial_guess == 2){

		for (int n = 0; n < electrode_material_interface.size(); n++){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				int site = electrode_material_interface[n].sites[m];
				int electrode_number = electrode_material_interface[n].get_electrode_number();

				if (electrode_number == 0 && get_MG_plus(site) != 0.0){
					MG_hole_concentration.data[site] = 0.5*get_MG_DOS(site);
				}
				if (electrode_number == 1 && get_MG_minus(site) != 0.0){
					MG_electron_concentration.data[site] = 0.5*get_MG_DOS(site);
				}

			}
		}
	}

	for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
		std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

		if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
			std::swap(site[0], site[1]);
		}

		if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){
			MG_hole_concentration.data[site[0]] = 0.5*get_MG_DOS(site[0]);
			MG_electron_concentration.data[site[1]] = 0.5*get_MG_DOS(site[1]);
		}
	}

}

void Morphology::calculate_MG_state_population(const PositionDependentParameter &electric_potential, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration){


	for (int n = 0; n < electrode_material_interface.size(); n++){
		for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
			int site = electrode_material_interface[n].sites[m];
			int electrode_number = electrode_material_interface[n].get_electrode_number();
			/*
			if (electrode_number == 0){

				MG_hole_concentration.data[site] = get_MG_DOS(site) * std::min(1.0/(exp((get_MG_plus(site) - get_work_function(electrode_number)
					+ (get_ele
					rode_potential(electrode_number) - electric_potential.data[site])))+1.0),1.0);
			}
			if (electrode_number == 1){
				MG_electron_concentration.data[site] = get_MG_DOS(site) * std::min(1.0/(exp((get_work_function(electrode_number) - get_MG_minus(site)
					+ (get_electrode_potential(electrode_number) - electric_potential.data[site])))+1.0),1.0);
			}*/

			if (electrode_number == 0 && get_MG_plus(site) != 0.0){

				double fermi_level = get_work_function(electrode_number) + get_electrode_potential(electrode_number);
				double energy_MG_plus = get_MG_plus(site) + electric_potential.data[site];

				double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * (1.0-misc::single_level_fermi_dist(fermi_level, energy_MG_plus));
				double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * misc::single_level_fermi_dist(fermi_level, energy_MG_plus);
				double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
				double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

				MG_hole_concentration.data[site] = get_MG_DOS(site) * (hole_con + electron_con_zero) / (electron_con + electron_con_zero + hole_con + hole_con_zero);
				// MG_hole_concentration.data[site] = get_MG_DOS(site) * (1-misc::single_level_fermi_dist(fermi_level, energy_MG_plus));
			}
			else if (electrode_number == 1 && get_MG_minus(site) != 0.0){

				double fermi_level = get_work_function(electrode_number) + get_electrode_potential(electrode_number);
				double energy_MG_minus = get_MG_minus(site) + electric_potential.data[site];

				double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_MG_minus(site) - get_work_function(electrode_number)
					+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
				double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_minus(site)
					- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
				double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
				double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

				MG_electron_concentration.data[site] = get_MG_DOS(site) * (electron_con + hole_con_zero) / (electron_con + electron_con_zero + hole_con + hole_con_zero);
			}

		}
	}

				// Reset all mid-gap state concentrations
	for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
		std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
		if (site[0] < points_y*(points_x - 1) && site[0] >= points_y && site[1] < points_y*(points_x - 1) && site[1] >= points_y){
			MG_electron_concentration.data[site[0]] = 0.0;
			MG_electron_concentration.data[site[1]] = 0.0;
			MG_hole_concentration.data[site[0]] = 0.0;
			MG_hole_concentration.data[site[1]] = 0.0;
		}
	}

	for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
		std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

		if (site[0] < points_y*(points_x - 1) && site[0] >= points_y && site[1] < points_y*(points_x - 1) && site[1] >= points_y){

			// Calculate the fermi level at before contact, we only consider E- states below this and E+ above this.
			double fermi_level = (std::min(get_HOMO(site[0]), get_HOMO(site[1])) + std::max(get_LUMO(site[0]), get_LUMO(site[1]))) / 2;


			// We assume that the E+ state is at site 0 and E- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
			if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
				std::swap(site[0], site[1]);
			}

			if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){
				double electron_con = electron_concentration.data[site[1]];
				double electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_LUMO(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

				double hole_con = hole_concentration.data[site[0]];
				double hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

				double C_n = get_electron_mobility(site[0], site[1]) / get_interface_relative_permittivity(pair_number);
				double C_p = get_hole_mobility(site[0], site[1]) / get_interface_relative_permittivity(pair_number);

				MG_hole_concentration.data[site[0]] = get_MG_DOS(site[0]) * (C_p * hole_con + C_n * electron_con_zero) / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

				electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_LUMO(site[1])) + 1.0), 1.0);

				hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_HOMO(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

				MG_electron_concentration.data[site[1]] = get_MG_DOS(site[1]) * (C_n * electron_con + C_p * hole_con_zero) / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

			}
			else if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
				// No mid-gap states, no error message needed.
			}
			else{
				std::cerr << "Something strange occured while calculating mid-gap state population." << std::endl;
			}

		}

	}
		
	

	return;

}

void Morphology::read_morphology(){

	std::ifstream file(morphology_file.c_str());
	int site = 0;

	while (file >> lattice[site]){
		site++;
			// Add error handling
	}

	file.close();

	std::cout << "Morphology read." << std::endl;

	return;
}

void Morphology::read_dopants(){

	std::ifstream file(dopants_file.c_str());
	int site = 0;

	double dopant_concentration;

	while (file >> dopant_concentration){
		dopants.data[site] = dopant_concentration/dopants.normalization_coef;
		site++;
		// Add error handling
	}

	file.close();

	std::cout << "Dopants read." << std::endl;

	return;
}

void Morphology::determine_interface_pairs(){


	int interface_pair_number = 1;

	int number_of_material_interfaces = material_interface.size();
	int number_of_electrode_material_interfaces = electrode_material_interface.size();



	if (convergence_boundary_condition){

		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){

				int site, site_x_minus, site_x_plus;

				site = i * points_y + j;

				site_x_minus = site - points_y;
				site_x_plus = site + points_y;

				if (i > 0){

					int x_minus_neigh = site_x_minus;
					while (x_minus_neigh > 0 && get_lattice_number(x_minus_neigh) == SPATIAL_MATERIAL_GAP)
						x_minus_neigh -= points_y;

					if ((is_cbc(site) || is_cbc(x_minus_neigh)) && !(is_cbc(site) && is_cbc(x_minus_neigh)) && (!is_electrode(site) && !is_electrode(x_minus_neigh))){

						std::array<int, 2> cbc_sites;

						if (is_cbc(site)){
							cbc_sites[0] = site;
							cbc_sites[1] = x_minus_neigh;
						}
						else if (is_cbc(x_minus_neigh)){
							cbc_sites[0] = x_minus_neigh;
							cbc_sites[1] = site;
						}

						cbc_interface.push_back(cbc_sites);

					}
				}

				if (i < points_x - 1){
					int x_plus_neigh = site_x_plus;
					while (x_plus_neigh < points_x*points_y && get_lattice_number(x_plus_neigh) == SPATIAL_MATERIAL_GAP)
						x_plus_neigh += points_y;

					if ((is_cbc(site) || is_cbc(x_plus_neigh)) && !(is_cbc(site) && is_cbc(x_plus_neigh)) && (!is_electrode(site) && !is_electrode(x_plus_neigh))){

						std::array<int, 2> cbc_sites;

						if (is_cbc(site)){
							cbc_sites[0] = site;
							cbc_sites[1] = x_plus_neigh;
						}
						else if (is_cbc(x_plus_neigh)){
							cbc_sites[0] = x_plus_neigh;
							cbc_sites[1] = site;
						}

						cbc_interface.push_back(cbc_sites);

					}
				}
				if (points_y > 1){

					int site_y_minus, site_y_plus;

					if (j != 0)
						site_y_minus = site - 1;
					else
						site_y_minus = site - 1 + points_y;

					if (j != points_y - 1)
						site_y_plus = site + 1;
					else
						site_y_plus = site + 1 - points_y;

					int y_minus_neigh = site_y_minus;

					while (get_lattice_number(y_minus_neigh) == SPATIAL_MATERIAL_GAP){

						if (y_minus_neigh%points_x != 0){
							y_minus_neigh -= 1;
						}
						else{
							y_minus_neigh += points_y - 1;
						}
					}

					if ((is_cbc(site) || is_cbc(y_minus_neigh)) && !(is_cbc(site) && is_cbc(y_minus_neigh)) && (!is_electrode(site) && !is_electrode(y_minus_neigh))){

						std::array<int, 2> cbc_sites;

						if (is_cbc(site)){
							cbc_sites[0] = site;
							cbc_sites[1] = y_minus_neigh;
						}
						else if (is_cbc(y_minus_neigh)){
							cbc_sites[0] = y_minus_neigh;
							cbc_sites[1] = site;
						}

						cbc_interface.push_back(cbc_sites);
					}
				
						int y_plus_neigh = site_y_plus;

						while (get_lattice_number(y_plus_neigh) == SPATIAL_MATERIAL_GAP){

							if (y_plus_neigh%points_x != 0){
								y_plus_neigh += 1;
							}
							else{
								y_plus_neigh -= points_y - 1;
							}
						}

						if ((is_cbc(site) || is_cbc(y_plus_neigh)) && !(is_cbc(site) && is_cbc(y_plus_neigh)) && (!is_electrode(site) && !is_electrode(y_plus_neigh))){

							std::array<int, 2> cbc_sites;

							if (is_cbc(site)){
								cbc_sites[0] = site;
								cbc_sites[1] = y_plus_neigh;
							}
							else if (is_cbc(y_plus_neigh)){
								cbc_sites[0] = y_plus_neigh;
								cbc_sites[1] = site;
							}

							cbc_interface.push_back(cbc_sites);

						}

					
				}

			}
		}
	}
	
	for (int i = 0; i < points_x; i++){
		for (int j = 0; j < points_y; j++){

			int site, site_x_minus, site_x_plus;

			site = i * points_y + j;

			site_x_minus = site - points_y;
			site_x_plus = site + points_y;

			if (i > 0){

				int x_minus_neigh = site_x_minus;
				while (x_minus_neigh > 0 && get_lattice_number(x_minus_neigh) == SPATIAL_MATERIAL_GAP)
					x_minus_neigh -= points_y;

				if ((is_electrode(site) || is_electrode(x_minus_neigh)) && (!is_cbc(site) && !is_cbc(x_minus_neigh))){

					if (is_electrode(site) && is_electrode(x_minus_neigh)){
						if (get_lattice_number(site) != get_lattice_number(x_minus_neigh)){
							std::cerr << "Two neighbouring sites have different electrode numbers. This is currently not implemented and may cause problems." << std::endl;
							interfaces[site][0] = INTERFACES_ELECTRODE_NUMBER;
							interfaces[x_minus_neigh][1] = INTERFACES_ELECTRODE_NUMBER;
						}
					}
					else{
						int electrode_number, material_number, electrode_site, material_site;
						if (is_electrode(site)){
							electrode_site = site;
							material_site = x_minus_neigh;
							electrode_number = get_electrode_number(get_lattice_number(site));
							material_number = get_material_number(get_lattice_number(x_minus_neigh));
						}
						else if (is_electrode(x_minus_neigh)){
							electrode_site = x_minus_neigh;
							material_site = site;
							electrode_number = get_electrode_number(get_lattice_number(electrode_site));
							material_number = get_material_number(get_lattice_number(site));
						}

						interfaces[site][0] = INTERFACES_ELECTRODE_NUMBER;
						interfaces[x_minus_neigh][1] = INTERFACES_ELECTRODE_NUMBER;

						for (int n = 0; n < number_of_electrode_material_interfaces; n++){

							bool correct_interface_material = (electrode_number == electrode_material_interface[n].electrode_number && material_number == electrode_material_interface[n].material_number);

							if (correct_interface_material){
								electrode_material_interface[n].add_electrode_interface_site(material_site);
							}
						}

					}
				}
			}
			
			if (i < points_x - 1){
				int x_plus_neigh = site_x_plus;
				while (x_plus_neigh < points_x*points_y && get_lattice_number(x_plus_neigh) == SPATIAL_MATERIAL_GAP)
					x_plus_neigh += points_y;

				if ((is_electrode(site) || is_electrode(x_plus_neigh)) && (!is_cbc(site) && !is_cbc(x_plus_neigh))){

					if (is_electrode(site) && is_electrode(x_plus_neigh)){
						if (get_lattice_number(site) != get_lattice_number(x_plus_neigh)){
							std::cerr << "Two neighbouring sites have different electrode numbers. This is currently not implemented and may cause problems." << std::endl;
							interfaces[site][1] = INTERFACES_ELECTRODE_NUMBER;
							interfaces[x_plus_neigh][0] = INTERFACES_ELECTRODE_NUMBER;
						}
					}
					else{
						int electrode_number, material_number, electrode_site, material_site;
						if (is_electrode(site)){
							electrode_site = site;
							material_site = x_plus_neigh;
							electrode_number = get_electrode_number(get_lattice_number(site));
							material_number = get_material_number(get_lattice_number(x_plus_neigh));
						}
						else if (is_electrode(x_plus_neigh)){
							electrode_site = x_plus_neigh;
							material_site = site;
							electrode_number = get_electrode_number(get_lattice_number(x_plus_neigh));
							material_number = get_material_number(get_lattice_number(site));
						}

						interfaces[site][1] = INTERFACES_ELECTRODE_NUMBER;
						interfaces[x_plus_neigh][0] = INTERFACES_ELECTRODE_NUMBER;

						for (int n = 0; n < number_of_electrode_material_interfaces; n++){

							bool correct_interface_material = (electrode_number == electrode_material_interface[n].electrode_number && material_number == electrode_material_interface[n].material_number);

							if (correct_interface_material){
								electrode_material_interface[n].add_electrode_interface_site(material_site);
							}
						}

					}
				}
			}
			if (points_y > 1){

				int site_y_minus, site_y_plus;

				if (j != 0)
					site_y_minus = site - 1;
				else
					site_y_minus = site - 1 + points_y;

				if (j != points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - points_y;

				int y_minus_neigh = site_y_minus;

				while (get_lattice_number(y_minus_neigh) == SPATIAL_MATERIAL_GAP){

					if (y_minus_neigh%points_x != 0){
						y_minus_neigh -= 1;
					}
					else{
						y_minus_neigh += points_y - 1;
					}
				}

				if ((is_electrode(site) || is_electrode(y_minus_neigh)) && (!is_cbc(site) && !is_cbc(y_minus_neigh))){

					if (is_electrode(site) && is_electrode(y_minus_neigh)){
						if (get_lattice_number(site) != get_lattice_number(y_minus_neigh)){
							std::cerr << "Two neighbouring sites have different electrode numbers. This is currently not implemented and may cause problems." << std::endl;
							interfaces[site][2] = INTERFACES_ELECTRODE_NUMBER;
							interfaces[y_minus_neigh][3] = INTERFACES_ELECTRODE_NUMBER;
						}
					}
					else{
						int electrode_number, material_number, electrode_site, material_site;
						if (is_electrode(site)){
							electrode_site = site;
							material_site = y_minus_neigh;
							electrode_number = get_electrode_number(get_lattice_number(site));
							material_number = get_material_number(get_lattice_number(y_minus_neigh));
						}
						else if (is_electrode(y_minus_neigh)){
							electrode_site = y_minus_neigh;
							material_site = site;
							electrode_number = get_electrode_number(get_lattice_number(y_minus_neigh));
							material_number = get_material_number(get_lattice_number(site));
						}

						interfaces[site][2] = INTERFACES_ELECTRODE_NUMBER;
						interfaces[y_minus_neigh][3] = INTERFACES_ELECTRODE_NUMBER;

						for (int n = 0; n < number_of_electrode_material_interfaces; n++){

							bool correct_interface_material = (electrode_number == electrode_material_interface[n].electrode_number && material_number == electrode_material_interface[n].material_number);

							if (correct_interface_material){
								electrode_material_interface[n].add_electrode_interface_site(material_site);
	
							}
						}

					}
				}

				int y_plus_neigh = site_y_plus;

				while (get_lattice_number(y_plus_neigh) == SPATIAL_MATERIAL_GAP){

					if (y_plus_neigh%points_x != 0){
						y_plus_neigh += 1;
					}
					else{
						y_plus_neigh -= points_y - 1;
					}
				}

				if ((is_electrode(site) || is_electrode(y_plus_neigh)) && (!is_cbc(site) && !is_cbc(y_plus_neigh))){

					if ((is_electrode(site) && is_electrode(y_plus_neigh))){
						if (get_lattice_number(site) != get_lattice_number(y_plus_neigh)){
							std::cerr << "Two neighbouring sites have different electrode numbers. This is currently not implemented and may cause problems." << std::endl;
							interfaces[site][3] = INTERFACES_ELECTRODE_NUMBER;
							interfaces[y_plus_neigh][2] = INTERFACES_ELECTRODE_NUMBER;
						}
					}
					else{
						int electrode_number, material_number, electrode_site, material_site;
						if (is_electrode(site)){
							electrode_site = site;
							material_site = y_plus_neigh;
							electrode_number = get_electrode_number(get_lattice_number(site));
							material_number = get_material_number(get_lattice_number(y_plus_neigh));
						}
						else if (is_electrode(y_plus_neigh)){
							electrode_site = y_plus_neigh;
							material_site = site;
							electrode_number = get_electrode_number(get_lattice_number(y_plus_neigh));
							material_number = get_material_number(get_lattice_number(site));
						}

						interfaces[site][3] = INTERFACES_ELECTRODE_NUMBER;
						interfaces[y_plus_neigh][2] = INTERFACES_ELECTRODE_NUMBER;

						for (int n = 0; n < number_of_electrode_material_interfaces; n++){

							bool correct_interface_material = (electrode_number == electrode_material_interface[n].electrode_number && material_number == electrode_material_interface[n].material_number);

							if (correct_interface_material){
								electrode_material_interface[n].add_electrode_interface_site(material_site);
							}
						}

					}
				}


			}


		}
	}

	std::cout << "Electrode interfaces located" << std::endl;

	for (int i = 0; i < points_x; i++){
		for (int j = 0; j < points_y; j++){

			int site, site_x_minus, site_x_plus;

			site = i * points_y + j;

			if (i > 0){
				site_x_minus = site - points_y;

				// If there is a spatial gap between materials, we look for the next point in the x-direction which is located on the other side of the gap.
				int x_minus_neigh = site_x_minus;
				while (get_lattice_number(x_minus_neigh) == SPATIAL_MATERIAL_GAP)
					x_minus_neigh -= points_y;

				if (get_lattice_number(site) != get_lattice_number(x_minus_neigh) && interfaces[site][0] == NOT_ALREADY_INTERFACE && (!is_cbc(site) && !is_cbc(x_minus_neigh))){

					int site_material = get_material_number(get_lattice_number(site));
					int neigh_material = get_material_number(get_lattice_number(x_minus_neigh));

					// Check that the two "neighbouring" sites are in different materials and store their pair number in a way that makes it easy to relate a site to it's pair number.
					if (site_material != neigh_material){

						interfaces[site][0] = interface_pair_number;
						interfaces[x_minus_neigh][1] = interface_pair_number;

						// Gather and store the different parameters of the pair.
						InterfacePair pair;

						pair.site[0] = x_minus_neigh;
						pair.site[1] = site;

						for (int n = 0; n < number_of_material_interfaces; n++){

							// Get the two material numbers of this material interface.
							std::array<int, 2> material_numbers = get_material_interface_pair(n);

							// Check if the material interface numbers match our pair.
							bool correct_interface_material = (site_material == material_numbers[0] && neigh_material == material_numbers[1]) ||
								(neigh_material == material_numbers[0] && site_material == material_numbers[1]);

							// Store the material interface number in the pair data.
							if (correct_interface_material){
								pair.interface_material = n;
								break;
							}
							else if (n == number_of_material_interfaces - 1){
								std::cerr << "No material properties could be found for the interface between material number " << lattice[site] << " and " << lattice[x_minus_neigh] << std::endl;
							}
						}
						// Calculate the spacing between the pair sites, unless there is a spatial gap, this is simply the regular spacing between sites.
						pair.spacing = (double)(site - x_minus_neigh) * spacing_x;

						// Add the pair data to a vector storing all the pairs.
						interface_pair_data.push_back(pair);

						// Pick the next pair number since the current one has been used.
						interface_pair_number++;

					}


				}
			}
			if (i < points_x - 1){
				site_x_plus = site + points_y;

				int x_plus_neigh = site_x_plus;

				while (get_lattice_number(x_plus_neigh) == SPATIAL_MATERIAL_GAP)
					x_plus_neigh += points_y;

				if (get_lattice_number(site) != get_lattice_number(x_plus_neigh) && interfaces[site][1] == NOT_ALREADY_INTERFACE && (!is_cbc(site) && !is_cbc(x_plus_neigh))){

					int site_material = get_material_number(get_lattice_number(site));
					int neigh_material = get_material_number(get_lattice_number(x_plus_neigh));

					if (site_material != neigh_material){
						interfaces[site][1] = interface_pair_number;
						interfaces[x_plus_neigh][0] = interface_pair_number;
					}

					InterfacePair pair;

					pair.site[0] = site;
					pair.site[1] = x_plus_neigh;

					for (int n = 0; n < number_of_material_interfaces; n++){

						std::array<int, 2> material_numbers = get_material_interface_pair(n);

						bool correct_interface_material = (site_material == material_numbers[0] && neigh_material == material_numbers[1]) ||
							(neigh_material == material_numbers[0] && site_material == material_numbers[1]);

						if (correct_interface_material){
							pair.interface_material = n;
							break;
						}
						else if (n == number_of_material_interfaces - 1){
							std::cerr << "No material properties could be found for the interface between material number " << lattice[site] << " and " << lattice[x_plus_neigh] << std::endl;
						}
					}

					pair.spacing = (double)(x_plus_neigh - site) * spacing_x;

					interface_pair_data.push_back(pair);

					interface_pair_number++;

				}
			}
			if (points_y > 1){
				
				int site_y_minus, site_y_plus;

				if (j != 0)
					site_y_minus = site - 1;
				else
					site_y_minus = site - 1 + points_y;

				if (j != points_y - 1)
					site_y_plus = site + 1;
				else
					site_y_plus = site + 1 - points_y;

					int y_minus_neigh = site_y_minus;

				while (get_lattice_number(y_minus_neigh) == SPATIAL_MATERIAL_GAP){

					if (y_minus_neigh%points_x != 0){
						y_minus_neigh -= 1;
					}
					else{
						y_minus_neigh += points_y - 1;
					}
				}

				if (get_lattice_number(site) != get_lattice_number(y_minus_neigh) && interfaces[site][2] == NOT_ALREADY_INTERFACE && (!is_cbc(site) && !is_cbc(y_minus_neigh))){

					int site_material = get_material_number(get_lattice_number(site));
					int neigh_material = get_material_number(get_lattice_number(y_minus_neigh));

					if (site_material != neigh_material){
						interfaces[site][2] = interface_pair_number;
						interfaces[y_minus_neigh][3] = interface_pair_number;

						InterfacePair pair;

						pair.site[0] = y_minus_neigh;
						pair.site[1] = site;

						for (int n = 0; n < number_of_material_interfaces; n++){

							std::array<int, 2> material_numbers = get_material_interface_pair(n);

							bool correct_interface_material = (site_material == material_numbers[0] && neigh_material == material_numbers[1]) ||
								(neigh_material == material_numbers[0] && site_material == material_numbers[1]);

							if (correct_interface_material){
								pair.interface_material = n;
								break;
							}
							else if (n == number_of_material_interfaces - 1){
								std::cerr << "No material properties could be found for the interface between material number " << lattice[site] << " and " << lattice[y_minus_neigh] << std::endl;
							}
						}

						pair.spacing = (double)(site - y_minus_neigh) * spacing_y; 

						interface_pair_data.push_back(pair);

						interface_pair_number++;

					}


				}

				int y_plus_neigh = site_y_plus;

				while (get_lattice_number(y_plus_neigh) == SPATIAL_MATERIAL_GAP){

					if (y_plus_neigh%points_x != points_y-1){
						y_plus_neigh += 1;
					}
					else{
						y_plus_neigh -= points_y - 1;
					}
				}

				if (get_lattice_number(site) != get_lattice_number(y_plus_neigh) && interfaces[site][3] == NOT_ALREADY_INTERFACE && (!is_cbc(site) && !is_cbc(y_plus_neigh))){

					int site_material = get_material_number(get_lattice_number(site));
					int neigh_material = get_material_number(get_lattice_number(y_plus_neigh));

					if (site_material != neigh_material){
						interfaces[site][3] = interface_pair_number;
						interfaces[y_plus_neigh][2] = interface_pair_number;

						InterfacePair pair;

						pair.site[0] = site;
						pair.site[1] = y_plus_neigh;

						for (int n = 0; n < number_of_material_interfaces; n++){

							std::array<int, 2> material_numbers = get_material_interface_pair(n);

							bool correct_interface_material = (site_material == material_numbers[0] && neigh_material == material_numbers[1]) ||
								(neigh_material == material_numbers[0] && site_material == material_numbers[1]);

							if (correct_interface_material){
								pair.interface_material = n;
								break;
							}
							else if (n == number_of_material_interfaces - 1){
								std::cerr << "No material properties could be found for the interface between material number " << lattice[site] << " and " << lattice[y_plus_neigh] << std::endl;
							}
						}

						pair.spacing = (double)(site - y_plus_neigh) * spacing_y; 

						interface_pair_data.push_back(pair);

						interface_pair_number++;

					}
				}
			}
		}
	}

	std::cout << "Organic interfaces located" << std::endl;

	is_electrode_interface_lattice.resize(points_x*points_y, false);
	for (int n = 0; n < electrode_material_interface.size(); n++){
		for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
			int site = electrode_material_interface[n].sites[m];
			is_electrode_interface_lattice[site] = true;
		}
	}


	return;

}

int Morphology::get_boundary_condition(int interface_number){
	return electrode_material_interface[interface_number].boundary_condition;
}

double Morphology::get_relative_permittivity(int site){

	int material_number;
	if (is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].relative_permittivity;

}

double Morphology::get_interface_relative_permittivity(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].relative_permittivity;
}

double Morphology::get_DOS(int site){

	int material_number;
	if (is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].DOS;

}

double Morphology::get_interface_DOS(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].DOS;
}

double Morphology::get_LUMO(int site){

	int material_number;
	if (is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].LUMO_level;

}

double Morphology::get_interface_LUMO(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].LUMO_level;
}

double Morphology::get_HOMO(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].HOMO_level;

}

double Morphology::get_interface_HOMO(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].HOMO_level;
}

double Morphology::get_generation_rate(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].generation_rate;
}

double Morphology::get_interface_generation_rate(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].generation_rate;
}

double Morphology::get_electron_mobility(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].mobility_electron;

}

double Morphology::get_electron_mobility(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_electron;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_electron;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_electron;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_electron;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_electron;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_electron;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_electron;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_electron;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_hole_mobility(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].mobility_hole;

}

double Morphology::get_hole_mobility(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_hole;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_hole;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_hole;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_hole;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_hole;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_hole;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_hole;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_hole;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_MG_minus(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].MG_minus_level;
}

double Morphology::get_MG_plus(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].MG_plus_level;
}

double Morphology::get_MG_DOS(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].MG_state_DOS;
}

double Morphology::get_bulk_reduced_recombination_coef(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].bulk_reduced_recombination_coef;
}

double Morphology::get_interface_reduced_recombination_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_reduced_recombination_coef;
}

double::Morphology::get_interface_trap_reduced_recombination_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_trap_reduced_recombination_coef;
}

double Morphology::get_bulk_hole_capture_coef(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].bulk_hole_capture_coef;
}

double Morphology::get_bulk_electron_capture_coef(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].bulk_electron_capture_coef;
}

double Morphology::get_interface_electron_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_capture_coef;
}

double Morphology::get_interface_hole_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_capture_coef;
}

double Morphology::get_interface_trap_electron_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_trap_electron_capture_coef;
}

double Morphology::get_interface_trap_hole_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_trap_hole_capture_coef;
}

double Morphology::get_negative_ion_mobility(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].mobility_negative_ion;

}

double Morphology::get_negative_ion_mobility(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_negative_ion;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_negative_ion;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_negative_ion;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_negative_ion;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_negative_ion;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_negative_ion;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_negative_ion;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_negative_ion;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_positive_ion_mobility(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].mobility_positive_ion;

}

double Morphology::get_positive_ion_mobility(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_positive_ion;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_positive_ion;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_positive_ion;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_positive_ion;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_positive_ion;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_positive_ion;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_positive_ion;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_positive_ion;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_negative_ion_eq_conc(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else if (is_electrode(site)){
		return 0.0;
	}
	else{
		material_number = get_material_number(lattice[site]);
	}
	return material[material_number].negative_ion_eq_conc;


}

double Morphology::get_positive_ion_eq_conc(int site){

	int material_number;
	if (convergence_boundary_condition && is_cbc(site)){
		int neigh_site;
		for (int n = 0; n < cbc_interface.size(); n++){
			if (cbc_interface[n][0] == site){
				neigh_site = cbc_interface[n][1];
				break;
			}
			if (n == cbc_interface.size() - 1){
				std::cerr << "Could not find cbc neighbour for site " << site << "." << std::endl;
			}
		}
		material_number = get_material_number(lattice[neigh_site]);
	}
	else if (is_electrode(site)){
		return 0.0;
	}
	else{
		material_number = get_material_number(lattice[site]);
	}

	return material[material_number].positive_ion_eq_conc;


}

int Morphology::get_electrode_number(int electrode_lattice_number){

	int electrode_number;
	for (int n = 0; n < electrode.size(); n++){
		if (electrode[n].lattice_number == electrode_lattice_number)
			return n;
	}
	std::cerr << "No electrode found for lattice number " << electrode_lattice_number << "." << std::endl;
	return NAN;
}

int Morphology::get_material_number(int material_lattice_number){

	int material_number;
	if (convergence_boundary_condition && convergence_boundary_condition_lattice_number == material_lattice_number){
		return convergence_boundary_condition_lattice_number;
	}
	for (int n = 0; n < material.size(); n++){
		if (material[n].lattice_number == material_lattice_number)
			return n;
	}
	std::cerr << "No semiconductor found for lattice number " << material_lattice_number << "." << std::endl;
	return NAN;
}

void Morphology::set_electrode_potential(int electrode_number, double potential){
	electrode[electrode_number].potential = potential;
	return;
}

double Morphology::get_electrode_potential(int electrode_number){
	return electrode[electrode_number].potential;
}

double Morphology::get_electrode_electron_concentration(int interface_number){
	return electrode_material_interface[interface_number].eq_electron_concentration;
}

double Morphology::get_electrode_hole_concentration(int interface_number){
	return electrode_material_interface[interface_number].eq_hole_concentration;
}


double Morphology::get_electrode_trap_reduced_recombination_coef(int interface_number){
	return electrode_material_interface[interface_number].interface_reduced_recombination_coef;
}

double Morphology::get_work_function(int electrode_number){
	return electrode[electrode_number].work_function;
}

double Morphology::get_active_layer_lenght_x(){
	return active_layer_lenght_x;
}

bool Morphology::is_interface(int site){

	bool interface = false;

	for (int n = 0; n < 4; n++){
		if (interfaces[site][n] != -1){
			interface = true;
		}
	}
	return interface;
}

int Morphology::get_interface_pair(int site_one, int site_two){

	for (int n = 0; n < interface_pair_data.size(); n++){
		if ((site_one == interface_pair_data[n].site[0] && site_two == interface_pair_data[n].site[1]) || (site_one == interface_pair_data[n].site[1] && site_two == interface_pair_data[n].site[0]))
			return n;
	}

	std::cerr << "No pair found for the selected sites." << std::endl;
	return 0;
}

int Morphology::get_pair(int site, int pair_number){
	if (site == interface_pair_data[pair_number].site[0]){
		return interface_pair_data[pair_number].site[1];
	}
	else if (site == interface_pair_data[pair_number].site[1]){
		return interface_pair_data[pair_number].site[0];
	}
	else{
		std::cerr << "Mismatch between pair number and site encountered." << std::endl;
		return 0;
	}
}

int Morphology::get_lattice_number(int site){
	return lattice[site];
}

std::array<int,2> Morphology::get_material_interface_pair(int interface_material_number){
	return material_interface[interface_material_number].material_numbers;
}

std::array<int, 2> Morphology::get_electrode_material_interface_pair(int electrode_material_interface_number){
	std::array<int, 2> pair = { electrode_material_interface[electrode_material_interface_number].electrode_number, electrode_material_interface[electrode_material_interface_number].material_number };
	return pair;
}

ElectrodeMaterialInterface Morphology::get_electrode_material_interface(int number){
	return electrode_material_interface[number];
}

int Morphology::get_electrode_material_interface_number(int electrode_number, int site){

	for (int n = 0; n < electrode_material_interface.size(); n++){
		bool electrode_found, site_found;
		electrode_found = electrode_material_interface[n].electrode_number == electrode_number;
		if (electrode_found){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				site_found = electrode_material_interface[n].sites[m] == site;
				if (site_found){
					return n;
				}
			}
		}
	}

	std::cerr << "No material found between the given electrode and site" << std::endl;
	return -1;

}

int Morphology::get_electrode_material_interface_number(int site){

	for (int n = 0; n < electrode_material_interface.size(); n++){

		bool site_found;
		for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
			site_found = electrode_material_interface[n].sites[m] == site;
			if (site_found){
				return n;
			}
		}

	}

	std::cerr << "No interface was found between site number " << site << "and any of the electrodes." << std::endl;
	return -1;

}

bool Morphology::is_electrode(int site){
	return is_electrode_lattice[site];
}

bool Morphology::negative_ion_transport(int site){

	if (is_electrode(site) || is_cbc(site)){
		return false;
	}
	else{
		int material_number = get_material_number(lattice[site]);
		return material[material_number].negative_ion_transport;
	}

}

bool Morphology::positive_ion_transport(int site){

	if (is_electrode(site) || is_cbc(site)){
		return false;
	}
	else{
		int material_number = get_material_number(lattice[site]);
		return material[material_number].positive_ion_transport;
	}

}

bool Morphology::is_electrode_interface(int site){
	return is_electrode_interface_lattice[site];
}

bool Morphology::is_cbc(int site){
	if (!convergence_boundary_condition){
		return false;
	}
	else{
		return is_cbc_lattice[site];
	}
}

bool Morphology::get_convergence_boundary_condition(){
	return convergence_boundary_condition;
}

bool Morphology::ion_transport_activated(){
	return ion_transport;
}