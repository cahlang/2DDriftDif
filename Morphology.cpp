#define BULK_RECOMBINATION 0
#define INTERFACE_RECOMBINATION 1
#define TRAP_ASSISTED_RECOMBINATION_MG 3
#define BULK_TRAP_ASSISTED_RECOMBINATION 2
#define INTERFACE_TRAP_ASSISTED_RECOMBINATION 4

#define BULK_GENERATION 0
#define INTERFACE_GENERATION 1

#define UNIFORM_GENERATION_PROFILE 0
#define EXPONENTIAL_FROM_ANODE_GENERATION_PROFILE 1
#define EXPONENTIAL_FROM_CATHODE_GENERATION_PROFILE 2

#define SPATIAL_MATERIAL_GAP 0
#define NOT_ALREADY_INTERFACE -1
#define INTERFACES_ELECTRODE_NUMBER -3
#define INTERFACES_CBC_NUMBER -97

#include "Morphology.h"
#include "DataTypes.h"
#include "Constants.h"
#include "Misc.h"

#include <cstdlib>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <string>
#include <sstream>
#include <iostream>
namespace pt = boost::property_tree;

#define MAX_NUMBER_MATERIALS 10	// Can be set to anything, but should be kept as small as possible.
#define MAX_NUMBER_MATERIAL_INTERFACES 40

void Morphology::initialize(pt::ptree &settings){

	points_x = settings.get<int>("general.points_x");
	points_y = settings.get<int>("general.points_y", 1);
	active_layer_length_x = settings.get<double>("general.length_x") / constant::length_norm;
	active_layer_length_y = settings.get<double>("general.length_y") / constant::length_norm;
	spacing_x = active_layer_length_x / ((double)points_x);
	spacing_y = active_layer_length_y / ((double)points_y);
	morphology_file = settings.get<std::string>("general.morphology_file");

	effective_temperature_activated = settings.get<bool>("general.effective_temperature",false);
	field_dependent_mobility_activated = settings.get<bool>("general.field_dependent_mobility", false);

	recombination_mode[0] = settings.get<bool>("general.recombination_mode_0");
	recombination_mode[1] = settings.get<bool>("general.recombination_mode_1");
	recombination_mode[2] = settings.get<bool>("general.recombination_mode_2");
	recombination_mode[3] = settings.get<bool>("general.recombination_mode_3");
	recombination_mode[4] = settings.get<bool>("general.recombination_mode_4");

	generation_mode[0] = settings.get<bool>("general.generation_mode_0");
	generation_mode[1] = settings.get<bool>("general.generation_mode_1");



	initial_guess = settings.get<int>("numerics.initial_guess", -1);

	lattice.resize(points_x*points_y, 0);
	interfaces.resize(points_x*points_y, std::vector<int>(4, NOT_ALREADY_INTERFACE));

	read_morphology();

	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::length_norm);

	dopants.initialize(norm_concentration);
	std::string no_dopants = "ignore";
	dopants_file = settings.get<std::string>("general.dopants_file", no_dopants);
	if (!(dopants_file == no_dopants))
	read_dopants();

	neumann_zero_boundary = settings.get<bool>("general.neumann_zero_boundary", false);
	if (neumann_zero_boundary){
		neumann_zero_boundary_lattice_number = settings.get<int>("general.neumann_zero_boundary_lattice_number");
	}

	MG_electron_concentration.initialize(norm_concentration);
	MG_hole_concentration.initialize(norm_concentration);

	if (field_dependent_mobility_activated){
		field_dependent_mobility_hole.initialize(constant::mobility_norm);
		field_dependent_mobility_electron.initialize(constant::mobility_norm);
	}
	read_material_data(settings);
	read_material_interface_data(settings);
	// identify_interfaces();

	map_electrodes();
	map_neumann_zero_boundary();
	determine_interface_pairs();

	generation_profile_type = settings.get<int>("generation.generation_type", 0);
	if (generation_profile_type != UNIFORM_GENERATION_PROFILE){
		create_generation_profile(settings);
	}


	std::cout << "All material parameters read." << std::endl;

	return;

}

void Morphology::read_material_data(pt::ptree & settings){

	double norm_potential, norm_concentration, norm_current, norm_rate;
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::length_norm);

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
			if (generation_mode[0] && generation_profile_type == UNIFORM_GENERATION_PROFILE)
				mat.generation_rate = input_material.get().get<double>("generation_rate") / norm_rate;

			mat.hole_trans_energy = input_material.get().get<double>("hole_transport_energy") / norm_potential;
			mat.hole_trans_DOS = input_material.get().get<double>("hole_transport_DOS") / norm_concentration;
			mat.electron_trans_energy = input_material.get().get<double>("electron_transport_energy") / norm_potential;
			mat.electron_trans_DOS = input_material.get().get<double>("electron_transport_DOS") / norm_concentration;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity");

			mat.MG_minus_level = input_material.get().get<double>("MG_minus_level", 0.0)/norm_potential;
			mat.MG_plus_level = input_material.get().get<double>("MG_plus_level", 0.0) / norm_potential;
			mat.MG_state_DOS = input_material.get().get<double>("MG_state_DOS", 0.0) / norm_concentration;

			mat.bulk_reduced_recombination_coef = input_material.get().get<double>("bulk_reduced_recombination_coef", 1.0);
		
			if (recombination_mode[BULK_TRAP_ASSISTED_RECOMBINATION]){
				mat.hole_trap_energy = input_material.get().get<double>("hole_trap_energy", 0.0) / norm_potential;
				mat.hole_trap_DOS = input_material.get().get<double>("hole_trap_DOS",0.0) / norm_concentration;
				mat.electron_trap_energy = input_material.get().get<double>("electron_trap_energy",0.0) / norm_potential;
				mat.electron_trap_DOS = input_material.get().get<double>("electron_trap_DOS", 0.0) / norm_concentration;

				mat.bulk_hole_trap_capture_coef_hole = input_material.get().get<double>("bulk_hole_trap_capture_coef_hole", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.bulk_hole_trap_capture_coef_electron = input_material.get().get<double>("bulk_hole_trap_capture_coef_electron", constant::elementary_charge * mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.bulk_electron_trap_capture_coef_hole = input_material.get().get<double>("bulk_electron_trap_capture_coef_hole", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.bulk_electron_trap_capture_coef_electron = input_material.get().get<double>("bulk_electron_trap_capture_coef_electron", constant::elementary_charge * mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);

			}

			mat.bulk_bimolecular_recombination_coef = input_material.get().get<double>("bulk_bimolecular_recombination_coef", constant::elementary_charge * (mat.mobility_hole + mat.mobility_electron) * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);

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

			if (is_field_dependent_mobility())
				mat.field_dependent_mobility_coef = input_material.get().get<double>("field_dependent_mobility_coef") / (sqrt(constant::length_norm / norm_potential));
			

			if (mat.negative_ion_transport || mat.positive_ion_transport)
				mat.max_ion_concentration = input_material.get().get<double>("max_ion_concentration") / norm_concentration;

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
			electrode.push_back(mat);
		
		}

	}

	std::cout << "Material parameters read." << std::endl;

	return;
}

void Morphology::read_material_interface_data(pt::ptree &settings){

	double norm_potential, norm_concentration, norm_current, norm_rate;
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::length_norm);
	double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

	for (int a = 0; a < MAX_NUMBER_MATERIAL_INTERFACES; a++){	// Material parameters are stored in the "material" vector. 

		std::ostringstream oss_node;
		oss_node << "material_interface_" << a;
		std::string node_name = oss_node.str();
		boost::optional<pt::ptree&> input_material = settings.get_child_optional(node_name.c_str());

		MaterialInterface mat;

		if (input_material){

			mat.mobility_hole = input_material.get().get<double>("mobility_hole") / constant::mobility_norm;
			mat.mobility_electron = input_material.get().get<double>("mobility_electron") / constant::mobility_norm;
			if(generation_mode[1])
				mat.generation_rate = input_material.get().get<double>("generation_rate") / norm_rate;

			mat.hole_trans_energy = input_material.get().get<double>("hole_transport_energy") / norm_potential;
			mat.hole_trans_DOS = input_material.get().get<double>("hole_transport_DOS") / norm_concentration;
			mat.electron_trans_energy = input_material.get().get<double>("electron_transport_energy") / norm_potential;
			mat.electron_trans_DOS = input_material.get().get<double>("electron_transport_DOS") / norm_concentration;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity");

			mat.interface_electron_transfer_velocity = input_material.get().get<double>("electron_transfer_velocity") / (constant::length_norm / norm_time);
			mat.interface_hole_transfer_velocity = input_material.get().get<double>("hole_transfer_velocity") / (constant::length_norm / norm_time);

			mat.material_numbers[0] = input_material.get().get<int>("material_1");
			mat.material_numbers[1] = input_material.get().get<int>("material_2");

			mat.interface_reduced_recombination_coef = input_material.get().get<double>("interface_reduced_recombination_coef", 1.0);
			mat.interface_hole_capture_coef = input_material.get().get<double>("interface_hole_capture_coef", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_electron_capture_coef = input_material.get().get<double>("interface_electron_capture_coef", constant::elementary_charge *  mat.mobility_electron  * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);

			mat.interface_bimolecular_recombination_coef_1 = input_material.get().get<double>("interface_bimolecular_recombination_coef_1", 
				constant::elementary_charge * (material[mat.material_numbers[0]].mobility_hole + material[mat.material_numbers[1]].mobility_electron) * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_bimolecular_recombination_coef_2 = input_material.get().get<double>("interface_bimolecular_recombination_coef_2",
				constant::elementary_charge * (material[mat.material_numbers[1]].mobility_hole + material[mat.material_numbers[0]].mobility_electron) * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			if (recombination_mode[INTERFACE_TRAP_ASSISTED_RECOMBINATION]){
				mat.electron_trap_material_number = input_material.get().get<int>("electron_trap_material_number");
				mat.hole_trap_material_number = input_material.get().get<int>("hole_trap_material_number");
				mat.interface_hole_trap_DOS = input_material.get().get<double>("interface_hole_trap_DOS", 0.0) / norm_concentration;
				mat.interface_electron_trap_DOS = input_material.get().get<double>("interface_electron_trap_DOS", 0.0) / norm_concentration;
				mat.interface_hole_trap_energy = input_material.get().get<double>("interface_hole_trap_energy", 0.0) / norm_potential;
				mat.interface_electron_trap_energy = input_material.get().get<double>("interface_electron_trap_energy", 0.0) / norm_potential;
				mat.interface_trap_reduced_recombination_coef = input_material.get().get<double>("interface_trap_reduced_recombination_coef", 1.0);
				mat.interface_hole_trap_hole_capture_coef = input_material.get().get<double>("interface_hole_trap_hole_capture_coef", constant::elementary_charge *  mat.mobility_hole  * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.interface_hole_trap_electron_capture_coef = input_material.get().get<double>("interface_hole_trap_electron_capture_coef", constant::elementary_charge *  mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.interface_electron_trap_hole_capture_coef = input_material.get().get<double>("interface_electron_trap_hole_capture_coef", constant::elementary_charge *  mat.mobility_hole  * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
				mat.interface_electron_trap_electron_capture_coef = input_material.get().get<double>("interface_electron_trap_electron_capture_coef", constant::elementary_charge *  mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			}
			if (is_field_dependent_mobility())
				mat.field_dependent_mobility_coef = input_material.get().get<double>("field_dependent_mobility_coef") / (sqrt(constant::length_norm / norm_potential));

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
			mat.hole_trans_energy = input_material.get().get<double>("hole_transport_energy", material[mat.material_number].hole_trans_energy * norm_potential) / norm_potential;
			mat.electron_trans_energy = input_material.get().get<double>("electron_transport_energy", material[mat.material_number].electron_trans_energy * norm_potential) / norm_potential;
			mat.relative_permittivity = input_material.get().get<double>("relative_permittivity", material[mat.material_number].relative_permittivity);
			mat.electron_trans_DOS = input_material.get().get<double>("electron_transport_DOS", material[mat.material_number].electron_trans_DOS * norm_concentration) / norm_concentration;
			mat.hole_trans_DOS = input_material.get().get<double>("electron_transport_DOS", material[mat.material_number].hole_trans_DOS * norm_concentration) / norm_concentration;

			mat.interface_reduced_recombination_coef = input_material.get().get<double>("interface_reduced_recombination_coef", 1.0);
			mat.interface_hole_capture_coef = input_material.get().get<double>("interface_hole_capture_coef", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_electron_capture_coef = input_material.get().get<double>("interface_electron_capture_coef", constant::elementary_charge * mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);

			mat.interface_trap_reduced_recombination_coef = input_material.get().get<double>("interface_trap_reduced_recombination_coef", 1.0);
			mat.interface_electron_trap_hole_capture_coef = input_material.get().get<double>("interface_electron_trap_hole_capture_coef", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_electron_trap_electron_capture_coef = input_material.get().get<double>("interface_electron_trap_electron_capture_coef", constant::elementary_charge * mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_hole_trap_hole_capture_coef = input_material.get().get<double>("interface_hole_trap_hole_capture_coef", constant::elementary_charge * mat.mobility_hole * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);
			mat.interface_hole_trap_electron_capture_coef = input_material.get().get<double>("interface_hole_trap_electron_capture_coef", constant::elementary_charge * mat.mobility_electron * constant::mobility_norm / (mat.relative_permittivity * constant::permittivity)) / (constant::elementary_charge * constant::mobility_norm / constant::permittivity);

			mat.surface_recombination = input_material.get().get<bool>("surface_recombination", false);
			if (mat.surface_recombination){
				mat.surface_recombination_velocity_electron = input_material.get().get<double>("surface_recombination_velocity_electron") / (constant::length_norm / norm_time);
				mat.surface_recombination_velocity_hole = input_material.get().get<double>("surface_recombination_velocity_hole") / (constant::length_norm / norm_time);
			}
			mat.boundary_condition = input_material.get().get<int>("boundary_condition");

			electrode_material_interface.push_back(mat);

		}



	}

	std::cout << "Interface parameters read." << std::endl;

	return;

}

void Morphology::map_neumann_zero_boundary(){

	is_neumann_zero_boundary.resize(points_x*points_y, false);

	for (int i = 0; i < points_x*points_y; i++){

		for (int j = 0; j < is_neumann_zero_boundary.size(); j++){
			if (lattice[i] == neumann_zero_boundary_lattice_number){
				is_neumann_zero_boundary[i] = true;
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
					double R = -get_bulk_reduced_recombination_coef(site) * (get_bulk_bimolecular_recombination_coef(site));
					electron_rate_coef.data[site] += R * hole_concentration.data[site];
					hole_rate_coef.data[site] += R * electron_concentration.data[site];
				}
			}
		}
	}

	if (recombination_mode[INTERFACE_RECOMBINATION]){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
			double R;
			
			if (get_material_number(lattice[site[0]]) == material_interface[interface_pair_data[pair_number].interface_material].material_numbers[0]){
				R = -get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_1(pair_number);
				hole_rate_coef.data[site[0]] += R * electron_concentration.data[site[1]];
				electron_rate_coef.data[site[1]] += R * hole_concentration.data[site[0]];
				R = -get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_2(pair_number);
				electron_rate_coef.data[site[0]] += R * hole_concentration.data[site[1]];
				hole_rate_coef.data[site[1]] += R * electron_concentration.data[site[0]];
			}
			else{
				R = -get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_2(pair_number);
				hole_rate_coef.data[site[0]] += R * electron_concentration.data[site[1]];
				electron_rate_coef.data[site[1]] += R * hole_concentration.data[site[0]];
				R = -get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_1(pair_number);
				electron_rate_coef.data[site[0]] += R * hole_concentration.data[site[1]];
				hole_rate_coef.data[site[1]] += R * electron_concentration.data[site[0]];
			}

		}
	}
	if (recombination_mode[BULK_TRAP_ASSISTED_RECOMBINATION]){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){

				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){

					double n_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_electron_trap_energy(site)),1.0);
					double p_1 = get_hole_trans_DOS(site)*std::min(exp(get_electron_trap_energy(site) - get_hole_trans_energy(site)),1.0);

					double R = - get_electron_trap_DOS(site) * get_bulk_electron_trap_capture_coef_electron(site) * get_bulk_electron_trap_capture_coef_hole(site)
						 / (get_bulk_electron_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_electron_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					electron_rate_coef.data[site] += R * hole_concentration.data[site];
					hole_rate_coef.data[site] += R * electron_concentration.data[site];

					n_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_hole_trap_energy(site)), 1.0);
					p_1 = get_hole_trans_DOS(site)*std::min(exp(get_hole_trap_energy(site) - get_hole_trans_energy(site)), 1.0);

					R = - get_hole_trap_DOS(site) * get_bulk_hole_trap_capture_coef_electron(site) * get_bulk_hole_trap_capture_coef_hole(site)
						/ (get_bulk_hole_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_hole_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					electron_rate_coef.data[site] += R * hole_concentration.data[site];
					hole_rate_coef.data[site] += R * electron_concentration.data[site];

				}

			}
		}

	}
	if (recombination_mode[INTERFACE_TRAP_ASSISTED_RECOMBINATION]){

		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
			int trap_site, opposite_site;
			if (get_material_number(get_lattice_number(site[0])) == get_hole_trap_material_number(pair_number)){
				trap_site = site[0];
				opposite_site = site[1];
			}
			else{
				trap_site = site[1];
				opposite_site = site[0];
			}

			double p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_hole_trap_energy(pair_number))
				+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site]);

			double n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
				get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			double n = electron_concentration.data[trap_site];
			double p = hole_concentration.data[opposite_site];
		
			double C_n = get_interface_hole_trap_electron_capture_coef(pair_number);
			double C_p = get_interface_hole_trap_hole_capture_coef(pair_number);

			double R = -get_interface_hole_trap_DOS(pair_number) * C_n * C_p / (C_n * (n + n_1) + C_p * (p + p_1));

			hole_rate_coef.data[opposite_site] += n * R;

			electron_rate_coef.data[trap_site] += p * R;

			if (get_material_number(get_lattice_number(site[0])) == get_electron_trap_material_number(pair_number)){
				trap_site = site[0];
				opposite_site = site[1];
			}
			else{
				trap_site = site[1];
				opposite_site = site[0];
			}

			p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_electron_trap_energy(pair_number))
				+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site]);

			n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
				get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			n = electron_concentration.data[trap_site] + electron_concentration.data[opposite_site];
			p = hole_concentration.data[trap_site] + hole_concentration.data[opposite_site];

			C_n = get_interface_electron_trap_electron_capture_coef(pair_number);
			C_p = get_interface_electron_trap_hole_capture_coef(pair_number);

			R = -get_interface_electron_trap_DOS(pair_number) * C_n * C_p / (C_n * (n + n_1) + C_p * (p + p_1));

			if(hole_concentration.data[trap_site] >= hole_concentration.data[opposite_site])
				hole_rate_coef.data[trap_site] += n * R;
			else
				hole_rate_coef.data[opposite_site] += n * R;

			if(electron_concentration.data[trap_site] >= electron_concentration.data[opposite_site])
				electron_rate_coef.data[trap_site] += p * R;
			else
				electron_rate_coef.data[opposite_site] += p * R;
			/*
			hole_rate_coef.data[trap_site] += (hole_concentration.data[trap_site] / p) * n * R;
			hole_rate_coef.data[opposite_site] += (hole_concentration.data[opposite_site] / p) * n * R;

			electron_rate_coef.data[trap_site] += (electron_concentration.data[trap_site] / n) * p * R;
			electron_rate_coef.data[opposite_site] += (electron_concentration.data[opposite_site] / n) * p * R;
			*/
		}

	}
	/*if (recombination_mode[TRAP_ASSISTED_RECOMBINATION_MG]){

		for (int n = 0; n < electrode_material_interface.size(); n++){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				int site = electrode_material_interface[n].sites[m];
				int electrode_number = electrode_material_interface[n].get_electrode_number();

				if (electrode_number == 0 && get_MG_plus(site) != 0.0){

					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_hole_trans_DOS(site) * std::min(1.0 / (exp(get_MG_plus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_electron_trans_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_plus(site)
						- (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
					double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

					electron_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * hole_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);
					hole_rate_coef.data[site] += -get_electrode_trap_reduced_recombination_coef(n) * get_MG_DOS(site) * electron_con / (electron_con + electron_con_zero + hole_con + hole_con_zero);

				}
				else if (electrode_number == 1 && get_MG_minus(site) != 0.0){
					double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_hole_trans_DOS(site) * std::min(1.0 / (exp(get_MG_minus(site) - get_work_function(electrode_number)
						+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
					double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_electron_trans_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_minus(site)
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
				double fermi_level = (std::min(get_hole_trans_energy(site[0]), get_hole_trans_energy(site[1])) + std::max(get_electron_trans_energy(site[0]), get_electron_trans_energy(site[1]))) / 2;


				// We assume that the CT+ state is at site 0 and CT- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
				if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
					std::swap(site[0], site[1]);
				}

				if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){

					double electron_con = electron_concentration.data[site[1]];
					double electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_electron_trans_energy(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

					double hole_con = hole_concentration.data[site[0]];
					double hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

					double C_n = get_interface_trap_electron_capture_coef(pair_number);
					double C_p = get_interface_trap_hole_capture_coef(pair_number);

					electron_rate_coef.data[site[1]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * hole_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate_coef.data[site[0]] += - get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * electron_con / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

					electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_electron_trans_energy(site[1])) + 1.0), 1.0);

					hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

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
	}*/

	for (int site = 0; site<points_x*points_y; site++){

		//hole_rate.data[site] += hole_rate_coef.data[site] * hole_concentration.data[site];
		//electron_rate.data[site] += electron_rate_coef.data[site] * electron_concentration.data[site];
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
					double intrinsic_rate = get_bulk_reduced_recombination_coef(site) * (get_bulk_bimolecular_recombination_coef(site)) * get_electron_trans_DOS(site) * get_hole_trans_DOS(site) * exp(get_electron_trans_energy(site) - get_hole_trans_energy(site));
					electron_rate.data[site] += intrinsic_rate;
					hole_rate.data[site] += intrinsic_rate;
				}
			}
		}
	}

	if (recombination_mode[INTERFACE_RECOMBINATION]){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();

			if (get_material_number(lattice[site[0]]) == material_interface[interface_pair_data[pair_number].interface_material].material_numbers[0]){
				double R = get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_1(pair_number);
				double intrinsic_concentration = get_electron_trans_DOS(site[1]) * get_hole_trans_DOS(site[0]) * exp(get_electron_trans_energy(site[1]) - get_hole_trans_energy(site[0]));
				hole_rate.data[site[0]] += R * intrinsic_concentration;
				electron_rate.data[site[1]] += R * intrinsic_concentration;
				R = get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_2(pair_number);
				intrinsic_concentration = get_electron_trans_DOS(site[0]) * get_hole_trans_DOS(site[1]) * exp(get_electron_trans_energy(site[0]) - get_hole_trans_energy(site[1]));
				electron_rate.data[site[0]] += R * intrinsic_concentration;
				hole_rate.data[site[1]] += R * intrinsic_concentration;
			}
			else{
				double R = get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_2(pair_number);
				double intrinsic_concentration = get_electron_trans_DOS(site[1]) * get_hole_trans_DOS(site[0]) * exp(get_electron_trans_energy(site[1]) - get_hole_trans_energy(site[0]));
				hole_rate.data[site[0]] += R * intrinsic_concentration;
				electron_rate.data[site[1]] += R * intrinsic_concentration;
				R = get_interface_reduced_recombination_coef(pair_number) * get_interface_bimolecular_recombination_coef_1(pair_number);
				intrinsic_concentration = get_electron_trans_DOS(site[0]) * get_hole_trans_DOS(site[1]) * exp(get_electron_trans_energy(site[0]) - get_hole_trans_energy(site[1]));
				electron_rate.data[site[0]] += R * intrinsic_concentration;
				hole_rate.data[site[1]] += R * intrinsic_concentration;
			}
		}
	}
	if (recombination_mode[BULK_TRAP_ASSISTED_RECOMBINATION]){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){

				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){

					double n_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_electron_trap_energy(site)), 1.0);
					double p_1 = get_hole_trans_DOS(site)*std::min(exp(get_electron_trap_energy(site) - get_hole_trans_energy(site)), 1.0);

					double R = get_electron_trap_DOS(site) * get_bulk_electron_trap_capture_coef_electron(site) * get_bulk_electron_trap_capture_coef_hole(site)
						/ (get_bulk_electron_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_electron_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					double intrinsic_concentration =  get_electron_trans_DOS(site) * get_hole_trans_DOS(site) * exp(get_electron_trans_energy(site) - get_hole_trans_energy(site));

					electron_rate.data[site] += R * intrinsic_concentration;
					hole_rate.data[site] += R * intrinsic_concentration;

					n_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_hole_trap_energy(site)), 1.0);
					p_1 = get_hole_trans_DOS(site)*std::min(exp(get_hole_trap_energy(site) - get_hole_trans_energy(site)), 1.0);

					R = get_hole_trap_DOS(site) * get_bulk_hole_trap_capture_coef_electron(site) * get_bulk_hole_trap_capture_coef_hole(site)
						/ (get_bulk_hole_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_hole_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					electron_rate.data[site] += R * intrinsic_concentration;
					hole_rate.data[site] += R * intrinsic_concentration;

				}

			}
		}

	}
	if (recombination_mode[INTERFACE_TRAP_ASSISTED_RECOMBINATION]){

		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
			int trap_site, opposite_site;
			if (get_material_number(get_lattice_number(site[0])) == get_hole_trap_material_number(pair_number)){
				trap_site = site[0];
				opposite_site = site[1];
			}
			else{
				trap_site = site[1];
				opposite_site = site[0];
			}

			double p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_hole_trap_energy(pair_number))
				+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site]);

			double n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
				get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			double n = electron_concentration.data[trap_site] + electron_concentration.data[opposite_site];
			double p = hole_concentration.data[trap_site] + hole_concentration.data[opposite_site];

			double C_n = get_interface_hole_trap_electron_capture_coef(pair_number);
			double C_p = get_interface_hole_trap_hole_capture_coef(pair_number);

			double R = get_interface_hole_trap_DOS(pair_number) * C_n * C_p / (C_n * (n + n_1) + C_p * (p + p_1));

			
			double  intrinsic_concentration = get_electron_trans_DOS(trap_site) * get_hole_trans_DOS(opposite_site) * exp(get_electron_trans_energy(trap_site) + electric_potential.data[trap_site] - get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			hole_rate.data[opposite_site] += intrinsic_concentration * R;
			electron_rate.data[trap_site] += intrinsic_concentration * R;

			intrinsic_concentration = get_electron_trans_DOS(opposite_site) * get_hole_trans_DOS(trap_site) * exp(get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site] - get_hole_trans_energy(trap_site) + electric_potential.data[trap_site]);

			hole_rate.data[trap_site] += intrinsic_concentration * R;
			electron_rate.data[opposite_site] += intrinsic_concentration * R;

			if (get_material_number(get_lattice_number(site[0])) == get_electron_trap_material_number(pair_number)){
				trap_site = site[0];
				opposite_site = site[1];
			}
			else{
				trap_site = site[1];
				opposite_site = site[0];
			}

			p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_electron_trap_energy(pair_number))
				+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site]);

			n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
				get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			n = electron_concentration.data[trap_site] + electron_concentration.data[opposite_site];
			p = hole_concentration.data[trap_site] + hole_concentration.data[opposite_site];

			C_n = get_interface_electron_trap_electron_capture_coef(pair_number);
			C_p = get_interface_electron_trap_hole_capture_coef(pair_number);

			R = get_interface_electron_trap_DOS(pair_number) * C_n * C_p / (C_n * (n + n_1) + C_p * (p + p_1));

			intrinsic_concentration = get_electron_trans_DOS(trap_site) * get_hole_trans_DOS(opposite_site) * exp(get_electron_trans_energy(trap_site) + electric_potential.data[trap_site] - get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

			hole_rate.data[opposite_site] += intrinsic_concentration * R;
			electron_rate.data[trap_site] += intrinsic_concentration * R;

			intrinsic_concentration = get_electron_trans_DOS(opposite_site) * get_hole_trans_DOS(trap_site) * exp(get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site] - get_hole_trans_energy(trap_site) + electric_potential.data[trap_site]);

			hole_rate.data[trap_site] += intrinsic_concentration * R;
			electron_rate.data[opposite_site] += intrinsic_concentration * R;
		}

	}
	/*if (recombination_mode[BULK_TRAP_ASSISTED_RECOMBINATION]){

		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){

				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){

					double intrinsic_concentration = get_hole_trans_DOS(site) * get_electron_trans_DOS(site) * exp(get_electron_trans_energy(site) - get_hole_trans_energy(site));

					double p_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_electron_trap_energy(site)), 1.0);
					double n_1 = get_hole_trans_DOS(site)*std::min(exp(get_electron_trap_energy(site) - get_hole_trans_energy(site)), 1.0);

					double R = get_bulk_electron_trap_capture_coef_electron(site) * get_bulk_electron_trap_capture_coef_hole(site)
						/ (get_bulk_electron_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_electron_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					electron_rate.data[site] += R * intrinsic_concentration;
					hole_rate.data[site] += R * intrinsic_concentration;

					p_1 = get_electron_trans_DOS(site)*std::min(exp(get_electron_trans_energy(site) - get_hole_trap_energy(site)), 1.0);
					n_1 = get_hole_trans_DOS(site)*std::min(exp(get_hole_trap_energy(site) - get_hole_trans_energy(site)), 1.0);

					R = get_bulk_hole_trap_capture_coef_electron(site) * get_bulk_hole_trap_capture_coef_hole(site)
						/ (get_bulk_hole_trap_capture_coef_electron(site) * (electron_concentration.data[site] + n_1) + get_bulk_hole_trap_capture_coef_hole(site) * (hole_concentration.data[site] + p_1));

					electron_rate.data[site] += R * intrinsic_concentration;
					hole_rate.data[site] += R * intrinsic_concentration;

				}

			}
		}
	}
	if (recombination_mode[TRAP_ASSISTED_RECOMBINATION_MG]){

		for (int n = 0; n < electrode_material_interface.size(); n++){
			for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
				int site = electrode_material_interface[n].sites[m];
				int electrode_number = electrode_material_interface[n].get_electrode_number();

				double intrinsic_concentration_squared = get_electron_trans_DOS(site) * get_hole_trans_DOS(site) * exp(get_electron_trans_energy(site) - get_hole_trans_energy(site));


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
				double fermi_level = (std::min(get_hole_trans_energy(site[0]), get_hole_trans_energy(site[1])) + std::max(get_electron_trans_energy(site[0]), get_electron_trans_energy(site[1]))) / 2;


				// We assume that the CT+ state is at site 0 and CT- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
				if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
					std::swap(site[0], site[1]);
				}

				if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){

					double electron_con = electron_concentration.data[site[1]];
					double electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_electron_trans_energy(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

					double hole_con = hole_concentration.data[site[0]];
					double hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

					double intrinsic_concentration_squared = pow(get_interface_DOS(pair_number), 2.0) * exp(get_electron_trans_energy(site[1]) - get_hole_trans_energy(site[0]) - (electric_potential.data[site[0]] - electric_potential.data[site[1]]));

					double C_n = get_interface_trap_electron_capture_coef(pair_number);
					double C_p = get_interface_trap_hole_capture_coef(pair_number);

					electron_rate.data[site[1]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[1]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));
					hole_rate.data[site[0]] += get_interface_trap_reduced_recombination_coef(pair_number) * get_MG_DOS(site[0]) * C_n * C_p * intrinsic_concentration_squared / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

					electron_con_zero = get_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_electron_trans_energy(site[1])) + 1.0), 1.0);

					hole_con_zero = get_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

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
	}*/

	if (generation_mode[BULK_GENERATION] && is_light_on()){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){
				int site = j + i * points_y;
				if (!is_electrode(site) && !is_cbc(site)){
					if (generation_profile_type == UNIFORM_GENERATION_PROFILE){
						double bulk_generation_rate = get_generation_rate(site);
						electron_rate.data[site] += bulk_generation_rate;
						hole_rate.data[site] += bulk_generation_rate;
					}
					else if (generation_profile_type == EXPONENTIAL_FROM_ANODE_GENERATION_PROFILE || generation_profile_type == EXPONENTIAL_FROM_CATHODE_GENERATION_PROFILE){
						electron_rate.data[site] += position_dependent_generation_rate[site];
						hole_rate.data[site] += position_dependent_generation_rate[site];
					}
				}
			}
		}
	}

	if (generation_mode[INTERFACE_GENERATION] && is_light_on()){
		for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
			std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
			if (get_electron_trans_energy(site[0]) < get_electron_trans_energy(site[1]) || get_hole_trans_energy(site[0]) < get_hole_trans_energy(site[1])){
				hole_rate.data[site[0]] += get_interface_generation_rate(pair_number);
				electron_rate.data[site[1]] += get_interface_generation_rate(pair_number);
			}
			else if (get_electron_trans_energy(site[0]) > get_electron_trans_energy(site[1]) || get_hole_trans_energy(site[0]) > get_hole_trans_energy(site[1])){
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

	for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
		std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
		int trap_site, opposite_site;
		if (get_material_number(get_lattice_number(site[0])) == get_hole_trap_material_number(pair_number)){
			trap_site = site[0];
			opposite_site = site[1];
		}
		else{
			trap_site = site[1];
			opposite_site = site[0];
		}

		MG_hole_concentration.data[trap_site] = 0.5 * get_interface_hole_trap_DOS(pair_number);

		if (get_material_number(get_lattice_number(site[0])) == get_electron_trap_material_number(pair_number)){
			trap_site = site[0];
			opposite_site = site[1];
		}
		else{
			trap_site = site[1];
			opposite_site = site[0];
		}

		MG_electron_concentration.data[trap_site] = 0.5 * get_interface_electron_trap_DOS(pair_number);

	}
}


void Morphology::calculate_MG_state_population(const PositionDependentParameter &electric_potential, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration){


	for (int pair_number = 0; pair_number < interface_pair_data.size(); pair_number++){
		std::array<int, 2> site = interface_pair_data[pair_number].get_pair_sites();
		int trap_site, opposite_site;
		if (get_material_number(get_lattice_number(site[0])) == get_hole_trap_material_number(pair_number)){
			trap_site = site[0];
			opposite_site = site[1];
		}
		else{
			trap_site = site[1];
			opposite_site = site[0];
		}

		double p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_hole_trap_energy(pair_number))
			+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site]);

		double n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
			get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_hole_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

		double n = electron_concentration.data[trap_site] + electron_concentration.data[opposite_site];
		double p = hole_concentration.data[trap_site] + hole_concentration.data[opposite_site];

		double C_n = get_interface_hole_trap_electron_capture_coef(pair_number);
		double C_p = get_interface_hole_trap_hole_capture_coef(pair_number);

		MG_hole_concentration.data[trap_site] = get_interface_hole_trap_DOS(pair_number) * (C_p * p + C_n * n_1) / (C_n * (n + n_1) + C_p * (p + p_1));

		if (get_material_number(get_lattice_number(site[0])) == get_electron_trap_material_number(pair_number)){
			trap_site = site[0];
			opposite_site = site[1];
		}
		else{
			trap_site = site[1];
			opposite_site = site[0];
		}

		p_1 = get_hole_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_hole_trans_energy(trap_site), get_interface_electron_trap_energy(pair_number))
			+ get_hole_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_hole_trans_energy(opposite_site) + electric_potential.data[opposite_site], get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site]);

		n_1 = get_electron_trans_DOS(trap_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number), get_electron_trans_energy(trap_site)) +
			get_electron_trans_DOS(opposite_site) * misc::single_level_fermi_dist(get_interface_electron_trap_energy(pair_number) + electric_potential.data[trap_site], get_electron_trans_energy(opposite_site) + electric_potential.data[opposite_site]);

		n = electron_concentration.data[trap_site] + electron_concentration.data[opposite_site];
		p = hole_concentration.data[trap_site] + hole_concentration.data[opposite_site];

		C_n = get_interface_electron_trap_electron_capture_coef(pair_number);
		C_p = get_interface_electron_trap_hole_capture_coef(pair_number);

		MG_electron_concentration.data[trap_site] = get_interface_electron_trap_DOS(pair_number) * (C_n * n + C_p * p_1) / (C_n * (n + n_1) + C_p * (p + p_1));

	}

}

void Morphology::create_generation_profile(pt::ptree & settings){

	position_dependent_generation_rate.resize(points_x*points_y, 0.0);
	
	norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	norm_concentration = constant::permittivity * norm_potential / (constant::elementary_charge * pow(constant::length_norm, 2.0));
	norm_current = constant::mobility_norm * constant::boltzmann * constant::temperature * norm_concentration / constant::length_norm;
	norm_rate = norm_current / (constant::elementary_charge * constant::length_norm);


	double photon_flux = settings.get<double>("generation.photon_flux")/(norm_rate*constant::length_norm);
	double alpha = settings.get<double>("generation.alpha") / (1.0 / constant::length_norm);
	int absorber_material_number = settings.get<int>("generation.absorber_material_number");


	std::vector<int> absorber_position;
	absorber_position.resize(points_y * 2);

	for (int j = 0; j < points_y; j++){
		bool position_zero_found = false;
		for (int i = 1; i < points_x-1; i++){

			int site = i*points_y + j;
			if (get_material_number(get_lattice_number(site)) == absorber_material_number && !position_zero_found){
				absorber_position[2*j] = i;
				position_zero_found = true;
			}
			if (position_zero_found && (get_material_number(get_lattice_number(site)) != absorber_material_number)){
				absorber_position[2*j + 1] = i-1;
				break;
			}
			if (i == points_x - 2)
				absorber_position[2 * j + 1] = i;
		}
	}

	for (int j = 0; j < points_y; j++){
		for (int i = 1; i < points_x - 1; i++){
			int site = i*points_y + j;

			if (i >= absorber_position[2 * j] && i <= absorber_position[2 * j + 1]){
				if (generation_profile_type == EXPONENTIAL_FROM_ANODE_GENERATION_PROFILE)
					position_dependent_generation_rate[site] = photon_flux*alpha*std::exp(-alpha * ((double)i - (double)absorber_position[2 * j]) * spacing_x);

				if (generation_profile_type == EXPONENTIAL_FROM_CATHODE_GENERATION_PROFILE){
					position_dependent_generation_rate[site] = photon_flux*alpha*std::exp(-alpha * ((double)absorber_position[2 * j + 1] - (double)i) * spacing_x);

				}
			}
		}
	}

	return;

}

/*
void Morphology::calculate_MG_state_population(const PositionDependentParameter &electric_potential, const PositionDependentParameter &electron_concentration, const PositionDependentParameter &hole_concentration){


	for (int n = 0; n < electrode_material_interface.size(); n++){
		for (int m = 0; m < electrode_material_interface[n].sites.size(); m++){
			int site = electrode_material_interface[n].sites[m];
			int electrode_number = electrode_material_interface[n].get_electrode_number();
			
		//	if (electrode_number == 0){
		//
			//	MG_hole_concentration.data[site] = get_MG_DOS(site) * std::min(1.0/(exp((get_MG_plus(site) - get_work_function(electrode_number)
			//		+ (get_ele
			//		rode_potential(electrode_number) - electric_potential.data[site])))+1.0),1.0);
		//	}
		//	if (electrode_number == 1){
		//		MG_electron_concentration.data[site] = get_MG_DOS(site) * std::min(1.0/(exp((get_work_function(electrode_number) - get_MG_minus(site)
		//			+ (get_electrode_potential(electrode_number) - electric_potential.data[site])))+1.0),1.0);
		//	}

			if (electrode_number == 0 && get_MG_plus(site) != 0.0){

				double fermi_level = get_work_function(electrode_number) + get_electrode_potential(electrode_number);
				double energy_MG_plus = get_MG_plus(site) + electric_potential.data[site];

				double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_electron_trans_DOS(site) * (1.0-misc::single_level_fermi_dist(fermi_level, energy_MG_plus));
				double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_hole_trans_DOS(site) * misc::single_level_fermi_dist(fermi_level, energy_MG_plus);
				double electron_con = get_electron_mobility(site) / get_relative_permittivity(site) * electron_concentration.data[site];
				double hole_con = get_hole_mobility(site) / get_relative_permittivity(site) * hole_concentration.data[site];

				MG_hole_concentration.data[site] = get_MG_DOS(site) * (hole_con + electron_con_zero) / (electron_con + electron_con_zero + hole_con + hole_con_zero);
				// MG_hole_concentration.data[site] = get_MG_DOS(site) * (1-misc::single_level_fermi_dist(fermi_level, energy_MG_plus));
			}
			else if (electrode_number == 1 && get_MG_minus(site) != 0.0){

				double fermi_level = get_work_function(electrode_number) + get_electrode_potential(electrode_number);
				double energy_MG_minus = get_MG_minus(site) + electric_potential.data[site];

				double electron_con_zero = get_hole_mobility(site) / get_relative_permittivity(site) * get_electron_trans_DOS(site) * std::min(1.0 / (exp(get_MG_minus(site) - get_work_function(electrode_number)
					+ (electric_potential.data[site] - get_electrode_potential(electrode_number))) + 1.0), 1.0);
				double hole_con_zero = get_electron_mobility(site) / get_relative_permittivity(site) * get_hole_trans_DOS(site) * std::min(1.0 / (exp(get_work_function(electrode_number) - get_MG_minus(site)
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
			double fermi_level = (std::min(get_hole_trans_energy(site[0]), get_hole_trans_energy(site[1])) + std::max(get_electron_trans_energy(site[0]), get_electron_trans_energy(site[1]))) / 2;


			// We assume that the E+ state is at site 0 and E- at site 1, if they happen to be the other way around, we swap site 0 and site 1.
			if (get_MG_plus(site[0]) == 0.0 && get_MG_minus(site[1]) == 0.0){
				std::swap(site[0], site[1]);
			}

			if (get_MG_plus(site[0]) != 0.0 && get_MG_minus(site[1]) != 0.0){
				double electron_con = electron_concentration.data[site[1]];
				double electron_con_zero = get_electron_trans_DOS(site[1]) * std::min(1.0 / (exp(get_MG_plus(site[0]) - get_electron_trans_energy(site[1]) + (electric_potential.data[site[0]] - electric_potential.data[site[1]])) + 1.0), 1.0);

				double hole_con = hole_concentration.data[site[0]];
				double hole_con_zero = get_hole_trans_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_plus(site[0])) + 1.0), 1.0);

				double C_n = get_electron_mobility(site[0], site[1]) / get_interface_relative_permittivity(pair_number);
				double C_p = get_hole_mobility(site[0], site[1]) / get_interface_relative_permittivity(pair_number);

				MG_hole_concentration.data[site[0]] = get_MG_DOS(site[0]) * (C_p * hole_con + C_n * electron_con_zero) / (C_n * (electron_con + electron_con_zero) + C_p * (hole_con + hole_con_zero));

				electron_con_zero = get_electron_trans_DOS(site[1]) * std::min(1.0 / (exp(get_MG_minus(site[1]) - get_electron_trans_energy(site[1])) + 1.0), 1.0);

				hole_con_zero = get_hole_trans_DOS(site[0]) * std::min(1.0 / (exp(get_hole_trans_energy(site[0]) - get_MG_minus(site[1]) - (electric_potential.data[site[1]] - electric_potential.data[site[0]])) + 1.0), 1.0);

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
*/
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



	if (neumann_zero_boundary){

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

void Morphology::effective_temperature(){

	//double effective_temperature = pow((pow(constant::temperature, 2.0) + pow(0.67*(get_electrode_potential(0) - get_electrode_potential(1))*constant::temperature * 1.5E-9 / (get_active_layer_length_x()*constant::length_norm), 2.0)), 0.5);
	double electric_field = abs((get_electrode_potential(0) - get_electrode_potential(1))*norm_potential / (get_active_layer_length_x()*constant::length_norm));

	for (int i = 0; i < material.size(); i++){
		if (!initial_values_saved){
			material[i].mobility_electron_zero = material[i].mobility_electron;
			material[i].mobility_hole_zero = material[i].mobility_hole; 
			material[i].generation_rate_initial = material[i].generation_rate;
		}
		material[i].mobility_electron = material[i].mobility_electron_zero*pow(10.0, (constant::temperature - 300.0) / 50.0);
		material[i].mobility_hole = material[i].mobility_hole_zero*pow(10.0, (constant::temperature - 300.0) / 50.0);
		material[i].generation_rate = material[i].generation_rate_initial / (1.0 + 1.205E-11*pow(constant::temperature, 2.0)*exp(2809.0 / constant::temperature)*exp(-0.4776*pow(electric_field,0.5)/constant::temperature));
	}
	for (int i = 0; i < material_interface.size(); i++){
		if (!initial_values_saved){
			material_interface[i].mobility_electron_zero = material_interface[i].mobility_electron;
			material_interface[i].mobility_hole_zero = material_interface[i].mobility_hole;
			material_interface[i].generation_rate_initial = material_interface[i].generation_rate;
		}
		material_interface[i].mobility_electron = material_interface[i].mobility_electron_zero*pow(10.0, (constant::temperature - 300.0) / 50.0);
		material_interface[i].mobility_hole = material_interface[i].mobility_hole_zero*pow(10.0, (constant::temperature - 300.0) / 50.0);
		material_interface[i].generation_rate = material_interface[i].generation_rate_initial / (1.0 + 1.205E-11*pow(constant::temperature, 2.0)*exp(2809.0 / constant::temperature)*exp(-0.4776*pow(electric_field, 0.5) / constant::temperature));
	}
	std::cout << 1.0 / (1.0 + 1.205E-11*pow(constant::temperature, 2.0)*exp(2809.0 / constant::temperature)*exp(-0.4776*pow(electric_field, 0.5) / constant::temperature)) << std::endl;
	initial_values_saved = true;
}

bool Morphology::is_effective_temperature(){
	return effective_temperature_activated;
}

bool Morphology::is_field_dependent_mobility(){
	return field_dependent_mobility_activated;
}

void Morphology::calculate_field_dependent_mobility(PositionDependentParameter electrical_potential){

	for (int i = 0; i < points_x-1; i++){
		for (int j = 0; j < points_y; j++){
			int site = i * points_y + j;
			int site_x_plus = site + points_y;
			double electric_field_x = (electrical_potential.data[site] - electrical_potential.data[site_x_plus]) / spacing_x;
			if (electric_field_x <= 0.0){
				field_dependent_mobility_electron.data[site] = get_electron_mobility_initial(site, site_x_plus) * exp(get_field_dependent_mobility_coef(site, site_x_plus) * sqrt(abs(electric_field_x)));
				field_dependent_mobility_hole.data[site] = get_hole_mobility_initial(site, site_x_plus) * exp(get_field_dependent_mobility_coef(site, site_x_plus) * sqrt(abs(electric_field_x)));
			}
			else{
				field_dependent_mobility_electron.data[site] = get_electron_mobility_initial(site, site_x_plus);
				field_dependent_mobility_hole.data[site] = get_hole_mobility_initial(site, site_x_plus);
			}
		}
	}


}

void Morphology::set_initial_field_dependent_mobility(){

	for (int i = 0; i < material.size(); i++){
		if (!initial_values_saved){
			material[i].mobility_electron_zero = material[i].mobility_electron;
			material[i].mobility_hole_zero = material[i].mobility_hole;
		}
	}
	for (int i = 0; i < material_interface.size(); i++){
		if (!initial_values_saved){
			material_interface[i].mobility_electron_zero = material_interface[i].mobility_electron;
			material_interface[i].mobility_hole_zero = material_interface[i].mobility_hole;
		}
	}
	initial_values_saved = true;

	for (int i = 0; i < points_x-1; i++){
		for (int j = 0; j < points_y; j++){
			int site = i * points_y + j;
			int site_x_plus = site + points_y;
			field_dependent_mobility_electron.data[site] = get_electron_mobility_initial(site, site_x_plus);
			field_dependent_mobility_hole.data[site] = get_hole_mobility_initial(site, site_x_plus);
		}
	}

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

double Morphology::get_relative_permittivity(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].relative_permittivity;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].relative_permittivity;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].relative_permittivity;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].relative_permittivity;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].relative_permittivity;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].relative_permittivity;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].relative_permittivity;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].relative_permittivity;
	}

	else{
		std::cerr << "Relative permittivity can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}

}

double Morphology::get_field_dependent_mobility_coef(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].field_dependent_mobility_coef;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].field_dependent_mobility_coef;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].field_dependent_mobility_coef;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].field_dependent_mobility_coef;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].field_dependent_mobility_coef;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].field_dependent_mobility_coef;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].field_dependent_mobility_coef;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].field_dependent_mobility_coef;
	}

	else{
		std::cerr << "Field-dependent mobility coeficient can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}

}

double Morphology::get_hole_trans_DOS(int site){

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
	return material[material_number].hole_trans_DOS;

}

double Morphology::get_electron_trans_DOS(int site){

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
	return material[material_number].electron_trans_DOS;

}


double Morphology::get_interface_DOS(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].DOS;
}

double Morphology::get_electron_trans_energy(int site){

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
	return material[material_number].electron_trans_energy;

}

double Morphology::get_interface_field_dependent_mobility_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].field_dependent_mobility_coef;
}

double Morphology::get_interface_electron_trans_energy(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].electron_trans_energy;
}

double Morphology::get_hole_trans_energy(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].hole_trans_energy;

}

double Morphology::get_interface_hole_trans_energy(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].hole_trans_energy;
}

double Morphology::get_generation_rate(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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

	if (is_field_dependent_mobility() && ((site_one - site_two) % points_y == 0)){
		if (site_one < site_two)
			return field_dependent_mobility_electron.data[site_one];
		else
			return field_dependent_mobility_electron.data[site_two];
	}
	else if (is_cbc(site_one)){
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

double Morphology::get_interface_electron_transfer_velocity(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_transfer_velocity;
}

double Morphology::get_electron_mobility_initial(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_electron_zero;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_electron_zero;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_electron_zero;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_electron_zero;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_electron_zero;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_electron_zero;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_electron_zero;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_electron_zero;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_hole_mobility(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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

	if (is_field_dependent_mobility() && ((site_one - site_two) % points_y == 0)){
		if (site_one < site_two)
			return field_dependent_mobility_hole.data[site_one];
		else
			return field_dependent_mobility_hole.data[site_two];
	}
	else if (is_cbc(site_one)){
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

double Morphology::get_interface_hole_transfer_velocity(int pair_number){

	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_transfer_velocity;
}

double Morphology::get_hole_mobility_initial(int site_one, int site_two){

	if (is_cbc(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_hole_zero;
	}
	else if (is_cbc(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_hole_zero;
	}
	else if (is_electrode(site_one)){
		return material[get_material_number(lattice[site_two])].mobility_hole_zero;
	}
	else if (is_electrode(site_two)){
		return material[get_material_number(lattice[site_one])].mobility_hole_zero;
	}
	else if (lattice[site_one] == lattice[site_two]){
		return material[get_material_number(lattice[site_one])].mobility_hole_zero;
	}
	else if (site_one < 0 || site_one > points_x*points_y){
		return material[get_material_number(lattice[site_two])].mobility_hole_zero;
	}
	else if (site_two < 0 || site_two > points_x*points_y){
		return material[get_material_number(lattice[site_one])].mobility_hole_zero;
	}
	else if (is_interface(site_one)){
		return material_interface[interface_pair_data[get_interface_pair(site_one, site_two)].interface_material].mobility_hole_zero;
	}

	else{
		std::cerr << "Mobility can't be determined between two sites of different materials unless they make up an interface pair." << std::endl;
		std::cerr << "The current sites are " << site_one << site_two << std::endl;
		return NAN;
	}
}

double Morphology::get_MG_minus(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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

double Morphology::get_interface_trap_reduced_recombination_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_trap_reduced_recombination_coef;
}

double Morphology::get_bulk_hole_capture_coef(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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

double Morphology::get_bulk_electron_trap_capture_coef_electron(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].bulk_electron_trap_capture_coef_electron;
}

double Morphology::get_bulk_electron_trap_capture_coef_hole(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].bulk_electron_trap_capture_coef_hole;
}

double Morphology::get_bulk_hole_trap_capture_coef_electron(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].bulk_hole_trap_capture_coef_electron;
}

double Morphology::get_bulk_hole_trap_capture_coef_hole(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].bulk_hole_trap_capture_coef_hole;
}


double Morphology::get_bulk_bimolecular_recombination_coef(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].bulk_bimolecular_recombination_coef;
}

double Morphology::get_electron_trap_energy(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].electron_trap_energy;
}

double Morphology::get_hole_trap_energy(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].hole_trap_energy;
}

double Morphology::get_electron_trap_DOS(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].electron_trap_DOS;
}

double Morphology::get_hole_trap_DOS(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].hole_trap_DOS;
}


double Morphology::get_interface_electron_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_capture_coef;
}

double Morphology::get_interface_hole_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_capture_coef;
}

double Morphology::get_interface_bimolecular_recombination_coef_1(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_bimolecular_recombination_coef_1;
}

double Morphology::get_interface_bimolecular_recombination_coef_2(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_bimolecular_recombination_coef_2;
}

double Morphology::get_interface_electron_trap_electron_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_trap_electron_capture_coef;
}

double Morphology::get_interface_electron_trap_hole_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_trap_hole_capture_coef;
}

double Morphology::get_interface_hole_trap_electron_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_trap_electron_capture_coef;
}

double Morphology::get_interface_hole_trap_hole_capture_coef(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_trap_hole_capture_coef;
}

double Morphology::get_interface_hole_trap_energy(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_trap_energy;
}

double Morphology::get_interface_electron_trap_energy(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_trap_energy;
}

double Morphology::get_interface_hole_trap_DOS(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_hole_trap_DOS;
}

double Morphology::get_interface_electron_trap_DOS(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].interface_electron_trap_DOS;
}

int Morphology::get_hole_trap_material_number(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].hole_trap_material_number;
}

int Morphology::get_electron_trap_material_number(int pair_number){
	return material_interface[interface_pair_data[pair_number].interface_material].electron_trap_material_number;
}


double Morphology::get_negative_ion_mobility(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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
	if (neumann_zero_boundary && is_cbc(site)){
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


double Morphology::get_max_ion_concentration(int site){

	int material_number;
	if (neumann_zero_boundary && is_cbc(site)){
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
	return material[material_number].max_ion_concentration;


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
	if (neumann_zero_boundary && neumann_zero_boundary_lattice_number == material_lattice_number){
		return neumann_zero_boundary_lattice_number;
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

double Morphology::get_surface_recombination_velocity_electron(int interface_number){
	return electrode_material_interface[interface_number].surface_recombination_velocity_electron;
}

double Morphology::get_surface_recombination_velocity_hole(int interface_number){
	return electrode_material_interface[interface_number].surface_recombination_velocity_hole;
}


double Morphology::get_active_layer_length_x(){
	return active_layer_length_x;
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
	if (!neumann_zero_boundary){
		return false;
	}
	else{
		return is_neumann_zero_boundary[site];
	}
}

bool Morphology::get_neumann_zero_boundary(){
	return neumann_zero_boundary;
}

bool Morphology::ion_transport_activated(){
	return ion_transport;
}

bool Morphology::surface_recombination_activated(int interface_number){
	return electrode_material_interface[interface_number].surface_recombination;
}

bool Morphology::is_light_on(){
	return light_on_flag;
}


void Morphology::light_off(){
	light_on_flag = false;
}

void Morphology::light_on(){
	light_on_flag = true;
}