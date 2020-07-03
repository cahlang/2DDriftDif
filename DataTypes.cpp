#include "DataTypes.h"
#include "Constants.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#define MAX_NUMBER_ELECTRODES 10

namespace pt = boost::property_tree;

int PositionDependentParameter::points_x;
int PositionDependentParameter::points_y;
double PositionDependentParameter::spacing_x;
double PositionDependentParameter::spacing_y;

void PositionDependentParameter::read_from_file(std::string file_name){

	std::ifstream from_file;
	from_file.open(file_name.c_str(), std::ios::in);

	for (int i = 0; i < points_x; i++){
		for (int j = 0; i < points_y; j++){
			if (from_file.eof()){
				std::cerr << "Reached end of " << file_name << " before all data was read." << std::endl;
				break;
			}
			int site = i * points_y + j;
			from_file >> data[site];
		}
	}

	if (!from_file.eof()){
		std::cerr << "All data was read before reaching the end of " << file_name << "." << std::endl;
	}

	return;
}

void PositionDependentParameter::read_from_file(std::string file_name_template, int file_number){

	std::string file_name;

	std::ostringstream oss_file;
	oss_file << file_name_template << "_" << file_number << ".dat";

	file_name = oss_file.str();

	std::ifstream from_file;
	from_file.open(file_name.c_str());

	double input;

	if (from_file.is_open()){
		for (int i = 0; i < points_x; i++){
			for (int j = 0; j < points_y; j++){
				if (from_file.eof()){
					std::cerr << "Reached end of " << file_name << " before all data was read." << std::endl;
					break;
				}
				int site = i * points_y + j;
				from_file >> input;
				data[site] = input / normalization_coef;
			}
		}

		if (from_file >> input){
			std::cerr << "All data was read before reaching the end of " << file_name << "." << std::endl;
		}
	}
	else{
		std::cerr << "Could not open " << file_name << "." << std::endl;
	}
	return;
}


void PositionDependentParameter::output_data(std::string file_name){

	std::ofstream to_file;
	to_file.open(file_name.c_str(), std::ios::ate);

	int site;

	for (int i = 0; i < points_x; i++){
		for (int j = 0; j < points_y; j++){
			site = j + i*points_y;
			to_file << data[site]*normalization_coef << " ";
		}
		to_file << "\n";
	}

	to_file.close();

	return;
}

void PositionDependentParameter::output_data(std::string file_name_template, int file_number){

	std::string file_name;

	std::ostringstream oss_file;
	oss_file << file_name_template << "_" << file_number << ".dat";

	file_name = oss_file.str();

	std::ofstream to_file;
	to_file.open(file_name.c_str(), std::ios::ate);

	int site;

	for (int i = 0; i < points_x; i++){
		for (int j = 0; j < points_y; j++){
			site = j + i*points_y;
			to_file << data[site] * normalization_coef << " ";
		}
		to_file << "\n";
	}

	to_file.close();

	return;
}

void PositionDependentParameter::initialize(double norm){

	calculate_this.resize(points_x*points_y, true);
	data.resize(points_x*points_y, 0.0);
	normalization_coef = norm;

	return;

}

void Experiment::initialize(pt::ptree settings){

	double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
	experiment_iv = settings.get<bool>("experiment.iv",false);
	electrostatics = settings.get<bool>("experiment.electrostatics",false);
	experiment_CELIV = settings.get<bool>("experiment.CELIV", false);
	experiment_PhotoCELIV = settings.get<bool>("experiment.PhotoCELIV", false);
	experiment_transient_photocurrent = settings.get<bool>("experiment.transient_photocurrent", false);

	output_full_data = settings.get<bool>("general.output_full_data", true);

	if (experiment_iv){

		data_points = settings.get<int>("experiment.data_points");
		anode_potential.resize(data_points, 0.0);
		cathode_potential.resize(data_points, 0.0);
		time.resize(data_points, 0.0);

		anode_potential[0] = settings.get<double>("experiment.applied_voltage") / norm_potential;
		cathode_potential[0] = 0.0;

		potential_step = settings.get<double>("experiment.potential_step") / norm_potential;

		iterate_forward = settings.get<bool>("numerics.iterate_forward", false);


		for (int i = 1; i < data_points; i++){
			anode_potential[i] = anode_potential[0] + (double)i * potential_step;
		}
	}
	else if (electrostatics){
		for (int a = 0; a < MAX_NUMBER_ELECTRODES; a++){
			double electrode_potential;
			std::ostringstream oss_node;
			oss_node << "experiment.electrode_" << a;
			std::string node_name = oss_node.str();
			electrode_potential = settings.get<double>(node_name.c_str(), 0.0) / norm_potential;
			electrode_potentials.push_back(electrode_potential);

		}
	}
	else if (experiment_CELIV){

		double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
		double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

		voltage_offset = settings.get<double>("experiment.offset_voltage", 0.0)/norm_potential;
		voltage_rise_speed = settings.get<double>("experiment.voltage_rise_speed")/(norm_potential/norm_time);
		time_step = settings.get<double>("experiment.time_step")/norm_time;
		pulse_length = settings.get<double>("experiment.pulse_length")/norm_time;

		data_points = settings.get<double>("experiment.data_points");

		anode_potential.resize(data_points, 0.0);
		cathode_potential.resize(data_points, 0.0);
		time.resize(data_points, 0.0);

		iterate_forward = settings.get<bool>("numerics.iterate_forward", true);

		for (int i = 0; i < data_points; i++){

			time[i] = (double)i * time_step;
			if (time[i] <= pulse_length){
				anode_potential[i] = -voltage_offset - voltage_rise_speed * time[i];
			}
			else{
				anode_potential[i] = -voltage_offset;
			}
		}

	}
	else if (experiment_PhotoCELIV){

		double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
		double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

		voltage_offset = settings.get<double>("experiment.offset_voltage", 0.0) / norm_potential;
		voltage_rise_speed = settings.get<double>("experiment.voltage_rise_speed") / (norm_potential / norm_time);
		time_step = settings.get<double>("experiment.time_step") / norm_time;
		pulse_length = settings.get<double>("experiment.pulse_length") / norm_time;
		light_delay = settings.get<double>("experiment.light_delay") / norm_time;

		data_points = settings.get<double>("experiment.data_points");

		anode_potential.resize(data_points, 0.0);
		cathode_potential.resize(data_points, 0.0);
		time.resize(data_points, 0.0);

		iterate_forward = settings.get<bool>("numerics.iterate_forward", true);

		for (int i = 0; i < data_points; i++){

			time[i] = (double)i * time_step;
			if (time[i] <= pulse_length + light_delay && time[i] > light_delay){
				anode_potential[i] = -voltage_offset - voltage_rise_speed * time[i];
			}
			else{
				anode_potential[i] = -voltage_offset;
			}
		}

	}
	else if (experiment_transient_photocurrent){

		double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
		double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

		pulse_length = settings.get<double>("experiment.pulse_length") / norm_time;
		time_step = settings.get<double>("experiment.time_step") / norm_time;
		voltage_offset = settings.get<double>("experiment.offset_voltage", 0.0) / norm_potential;

		data_points = settings.get<double>("experiment.data_points");

		anode_potential.resize(data_points, 0.0);
		cathode_potential.resize(data_points, 0.0);
		time.resize(data_points, 0.0);

		for (int i = 0; i < data_points; i++){

			time[i] = (double)i * time_step;

			anode_potential[i] = voltage_offset;

		}

	}

	std::cout << "Measurement initialized." << std::endl;

	return;

}

void Experiment::next_data_point(){
	current_data_point++;
	return;
}
void Experiment::set_current_data_point(int number){
	current_data_point = number;
	return;
}
int Experiment::get_current_data_point(){
	return current_data_point;
}