#include "Control.h"
#include "DataTypes.h"
#include "Potential.h"
#include "NegativeMobileCharge.h"
#include "PositiveMobileCharge.h"
#include "Morphology.h"
#include "Constants.h"
#include "Calc.h"

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <string>
#include <sstream>
#include <algorithm>

namespace pt = boost::property_tree;

double constant::length_norm;
double constant::temperature;
double constant::mobility_norm;
int constant::max_iteration_step;
bool mid_gap_states;

namespace control{

	bool read_settings(const std::string &file_name, Morphology &morphology, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, Experiment &measurement){

		std::cout << "Reading settings..." << std::endl;

		pt::ptree settings;

		pt::read_ini(file_name, settings);

		PositionDependentParameter::points_x = settings.get<int>("general.points_x");
		PositionDependentParameter::points_y = settings.get<int>("general.points_y", 1);


		constant::length_norm = settings.get<double>("constants.length_normalization");
		constant::mobility_norm = settings.get<double>("constants.mobility_normalization");
		constant::temperature = settings.get<double>("constants.temperature");

		PositionDependentParameter::spacing_x = settings.get<double>("general.length_x") / ((double)PositionDependentParameter::points_x) / constant::length_norm;
		PositionDependentParameter::spacing_y = settings.get<double>("general.length_y") / ((double)PositionDependentParameter::points_y) / constant::length_norm;

		mid_gap_states = settings.get<bool>("general.mid_gap_states", false);

		constant::max_iteration_step = settings.get<int>("numerics.max_iteration_step", 1E6);

		morphology.initialize(settings);

		if (morphology.ion_transport_activated()){
			negative_ion.initialize(settings);
			positive_ion.initialize(settings);
		}	

		potential.initialize(settings);
		electron.initialize(settings);
		hole.initialize(settings);

		measurement.initialize(settings);

		std::cout << "All settings sucessfully read" << std::endl;

		return true;

	}

	bool solver(Morphology &material, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, 
		NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, MeasuredCurrent &outer_circuit_current, PositionDependentParameter &net_rate){

		bool has_converged = false;
		bool has_converged_inner = false;
		int counter = 0;

		while (!has_converged && counter <= constant::max_iteration_step){

			determine_rates(material, net_rate, electron, hole, potential.electrical);

			if (electron.concentration.points_y == 1)
				potential.solve_one_d(material, electron.concentration, hole.concentration, negative_ion.concentration, positive_ion.concentration);
			else
				potential.solve_inverse(material, electron.concentration, hole.concentration, negative_ion.concentration, positive_ion.concentration);


			if (electron.concentration.points_y == 1){
				electron.solve_one_d(material, potential);
				hole.solve_one_d(material, potential);
			}
			else{
				electron.solve_inverse(material, potential);
				hole.solve_inverse(material, potential);
			}

			if (mid_gap_states)
				material.calculate_MG_state_population(potential.electrical, electron.concentration, hole.concentration);

			if (material.ion_transport_activated()){
				negative_ion.ion_solve(material, potential);
				positive_ion.ion_solve(material, potential);
			}

			

			counter++;

			if (potential.has_converged() && electron.has_converged() && hole.has_converged()){
				has_converged = true;
				if (!potential.final_iterative_stage_reached() && !electron.final_iterative_stage_reached() && !hole.final_iterative_stage_reached()){
					has_converged = false;
					electron.next_iterative_stage();
					hole.next_iterative_stage();
					potential.next_iterative_stage();
				}
			}
			else if (counter % 1000 == 0){

				if (potential.has_converged())
					printf("The solution for the potential has been found...\n");

				if (electron.has_converged() && electron.concentration.points_y != 1)
					printf("The solution for the electron concentration has been found...\n");

				if (hole.has_converged() && hole.concentration.points_y != 1)
					printf("The solution for the hole concentration has been found...\n");
				if (material.ion_transport_activated() && negative_ion.has_converged()){
					printf("The solution for the negative ion concentration has been found...\n");
				}
				if (material.ion_transport_activated() && positive_ion.has_converged()){
					printf("The solution for the positive ion concentration has been found...\n");
				}
			}

			
		}

		determine_rates(material, net_rate, electron, hole, potential.electrical);

		electron.calculate_current(material, potential);
		hole.calculate_current(material, potential);
		if (material.ion_transport_activated()){
			negative_ion.calculate_ion_current(material, potential);
			positive_ion.calculate_ion_current(material, potential);
		}


		potential.calculate_electrochemical_potential(material, electron.concentration, hole.concentration);

		PositionDependentParameter negative_charge_carrier_current, positive_charge_carrier_current;
		negative_charge_carrier_current = electron.current_x;
		positive_charge_carrier_current = hole.current_x;

		if (material.ion_transport_activated()){
			for (int site = 0; site < electron.current_x.data.size(); site++){

				if (material.negative_ion_transport(site))
					negative_charge_carrier_current.data[site] += negative_ion.current_x.data[site];

				if (material.positive_ion_transport(site))
					positive_charge_carrier_current.data[site] += positive_ion.current_x.data[site];

			}
		}

		outer_circuit_current = calc::outer_circuit_current(material, negative_charge_carrier_current, positive_charge_carrier_current, potential.displacement_current, net_rate);

		printf("Number of iterations required: %d\n", counter);

		potential.reset_iterative_stage();
		electron.reset_iterative_stage();
		hole.reset_iterative_stage();

		return true;

	}

	bool solver(Morphology &material, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, MeasuredCurrent &outer_circuit_current, PositionDependentParameter &net_rate,
		PositionDependentParameter &previous_holecon_concentration, PositionDependentParameter &previous_electroncon_concentration, PositionDependentParameter &previous_neg_ion_con, PositionDependentParameter &previous_pos_ion_con, PositionDependentParameter &previous_potential, double time_step){

		bool has_converged = false;
		int counter = 0;

		while (!has_converged && counter <= constant::max_iteration_step){

			determine_rates(material, net_rate, electron, hole, potential.electrical);
			if (electron.concentration.points_y == 1){
				potential.solve_one_d(material, electron.concentration, hole.concentration, negative_ion.concentration, positive_ion.concentration);
			}
			else{
				potential.solve_inverse(material, electron.concentration, hole.concentration, negative_ion.concentration, positive_ion.concentration);
			}
		

			if (electron.concentration.points_y == 1){
				electron.solve_one_d(material, potential, previous_electroncon_concentration, time_step);
				hole.solve_one_d(material, potential, previous_holecon_concentration, time_step);
			}
			else{
				electron.solve_inverse(material, potential, previous_electroncon_concentration, time_step);
				hole.solve_inverse(material, potential, previous_holecon_concentration, time_step);
			}
			
		//	if (mid_gap_states)
		//		material.calculate_MG_state_population(potential.electrical, electron.concentration, hole.concentration);

			if (material.ion_transport_activated()){
				negative_ion.ion_solve(material, potential);
				positive_ion.ion_solve(material, potential);
			}

			counter++;

			if (potential.has_converged() && electron.has_converged() && hole.has_converged()){
				has_converged = true;
				if (!potential.final_iterative_stage_reached() && !electron.final_iterative_stage_reached() && !hole.final_iterative_stage_reached()){
					has_converged = false;
					electron.next_iterative_stage();
					hole.next_iterative_stage();
					potential.next_iterative_stage();
				}
			}
			else if (counter % 1000 == 0){

				if (potential.has_converged())
					printf("The solution for the potential has been found...\n");

				if (electron.has_converged())
					printf("The solution for the electron concentration has been found...\n");

				if (hole.has_converged())
					printf("The solution for the hole concentration has been found...\n");
				if (material.ion_transport_activated() && negative_ion.has_converged()){
					printf("The solution for the negative ion concentration has been found...\n");
				}
				if (material.ion_transport_activated() && positive_ion.has_converged()){
					printf("The solution for the positive ion concentration has been found...\n");
				}
			}


		}

		determine_rates(material, net_rate, electron, hole, potential.electrical);

		electron.calculate_current(material, potential);
		hole.calculate_current(material, potential);
		potential.calculate_displacement_current(material, previous_potential, time_step);
		if (material.ion_transport_activated()){
			negative_ion.calculate_ion_current(material, potential);
			positive_ion.calculate_ion_current(material, potential);
		}


		potential.calculate_electrochemical_potential(material, electron.concentration, hole.concentration);

		PositionDependentParameter negative_charge_carrier_current, positive_charge_carrier_current;
		negative_charge_carrier_current = electron.current_x;
		positive_charge_carrier_current = hole.current_x;

		if (material.ion_transport_activated()){
			for (int site = 0; site < electron.current_x.data.size(); site++){

				if (material.negative_ion_transport(site))
					negative_charge_carrier_current.data[site] += negative_ion.current_x.data[site];

				if (material.positive_ion_transport(site))
					positive_charge_carrier_current.data[site] += positive_ion.current_x.data[site];

			}
		}

		outer_circuit_current = calc::outer_circuit_current(material, negative_charge_carrier_current, positive_charge_carrier_current, potential.displacement_current, net_rate);

		printf("Number of iterations required: %d\n", counter);

		potential.reset_iterative_stage();
		electron.reset_iterative_stage();
		hole.reset_iterative_stage();

		return true;
	}

	void determine_rates(Morphology &material, PositionDependentParameter &net_rate, NegativeMobileCharge &electron, PositiveMobileCharge &hole, PositionDependentParameter &potential){

		std::fill(electron.rate.data.begin(), electron.rate.data.end(), 0.0);
		std::fill(hole.rate.data.begin(), hole.rate.data.end(), 0.0);
		std::fill(electron.rate_coef.data.begin(), electron.rate_coef.data.end(), 0.0);
		std::fill(hole.rate_coef.data.begin(), hole.rate_coef.data.end(), 0.0);
		std::fill(net_rate.data.begin(), net_rate.data.end(), 0.0);

		material.calculate_generation_rate(electron.concentration, electron.rate, hole.concentration, hole.rate, potential, net_rate);
		material.calculate_recombination_rates(electron.concentration, electron.rate, electron.rate_coef, hole.concentration, hole.rate, hole.rate_coef, potential, net_rate);

		return;

	}

	void run_measurement(Morphology &morphology, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole,
		NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, Experiment &measurement){

		if (measurement.experiment_iv){

			measurement.set_current_data_point(0);

			std::cout << "Initiating calculations" << std::endl;

			while (measurement.current_data_point < measurement.data_points){

				int i = measurement.get_current_data_point();

				MeasuredCurrent current;

				PositionDependentParameter net_rate;
				net_rate.initialize(electron.rate.normalization_coef);

				morphology.set_electrode_potential(0, measurement.anode_potential[i] - (morphology.get_work_function(0)-morphology.get_work_function(1)));
				morphology.set_electrode_potential(1, measurement.cathode_potential[i]);

				if (morphology.is_effective_temperature()){
					morphology.effective_temperature();
				}

				potential.set_boundary_conditions(morphology);
				if (!measurement.iterate_forward || i == 0)
					potential.set_initial_guess(morphology, "initial_guess\\potential", i);
				electron.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					electron.set_initial_guess(morphology, potential.electrical, "initial_guess\\electroncon", i);
				hole.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					hole.set_initial_guess(morphology, potential.electrical, "initial_guess\\holecon", i);
				if (!measurement.iterate_forward || i == 0)
					morphology.set_MG_initial_guess();
				if (morphology.ion_transport_activated()){
					negative_ion.set_ion_boundary_conditions(morphology);
					positive_ion.set_ion_boundary_conditions(morphology);
				}

				std::cout << "Boundary conditions and initial guess set." << std::endl;
				std::cout << "Solving for data point " << i << "..." << std::endl;

				solver(morphology, potential, electron, hole, negative_ion, positive_ion, current, net_rate);
				if (measurement.output_full_data){
					potential.electrical.output_data("potential", i);
					electron.concentration.output_data("electroncon", i);
					hole.concentration.output_data("holecon", i);
					electron.current_x.output_data("ecurrentx", i);
					hole.current_x.output_data("hcurrentx", i);
					if (potential.electrical.points_y > 1){
						hole.current_y.output_data("hcurrenty", i);
						electron.current_y.output_data("ecurrenty", i);
					}
					potential.electrochemical_electron.output_data("electron_fermi_level", i);
					potential.electrochemical_hole.output_data("hole_fermi_level", i);
					electron.rate.output_data("erate", i);
					hole.rate.output_data("hrate", i);
					if (mid_gap_states){
						morphology.MG_electron_concentration.output_data("MG_minus", i);
						morphology.MG_hole_concentration.output_data("MG_plus", i);
					}
					if (morphology.ion_transport_activated()){
						negative_ion.concentration.output_data("negative_ion_con", i);
						positive_ion.concentration.output_data("positive_ion_con", i);
						negative_ion.current_x.output_data("negative_ion_currentx", i);
						positive_ion.current_x.output_data("positive_ion_currentx", i);
					}
					net_rate.output_data("rate", i);
				}
				std::ofstream fout("iv.dat", std::ios::app);
				fout << (measurement.anode_potential[i] - measurement.cathode_potential[i]) * constant::boltzmann * constant::temperature / constant::elementary_charge << "	" << std::scientific << current.outer_circuit * current.normalization_coef << "	" << current.rms_error * current.normalization_coef << std::endl;
				fout.close();

				measurement.next_data_point();

			}

		}
		else if (measurement.electrostatics){

			measurement.set_current_data_point(0);

			std::cout << "Initiating calculations" << std::endl;

			MeasuredCurrent discard;

			PositionDependentParameter net_rate;
			net_rate.initialize(electron.rate.normalization_coef);

			for (int i = 0; i < morphology.electrode.size();i++)
				morphology.set_electrode_potential(i, measurement.electrode_potentials[i] - (morphology.get_work_function(i) - morphology.get_work_function(0)));

			potential.set_boundary_conditions(morphology);
			//potential.set_initial_guess(morphology, "ignore", 0);
			electron.set_boundary_conditions(morphology, potential.electrical);
			//electron.set_initial_guess(morphology, potential.electrical, "ignore",0);
			hole.set_boundary_conditions(morphology, potential.electrical);
			//hole.set_initial_guess(morphology, potential.electrical, "ignore",0);
			//morphology.set_MG_initial_guess();
			if (morphology.ion_transport_activated()){
				negative_ion.set_ion_boundary_conditions(morphology);
				positive_ion.set_ion_boundary_conditions(morphology);
			}

			std::cout << "Boundary conditions and initial guess set." << std::endl;
			std::cout << "Solving..." << std::endl;

			solver(morphology, potential, electron, hole, negative_ion, positive_ion, discard, net_rate);

				potential.electrical.output_data("potential.dat");
				electron.concentration.output_data("electroncon.dat");
				hole.concentration.output_data("holecon.dat");
				electron.current_x.output_data("ecurrentx.dat");
				hole.current_x.output_data("hcurrentx.dat");
				potential.electrochemical_electron.output_data("electron_fermi_level.dat");
				potential.electrochemical_hole.output_data("hole_fermi_level.dat");
				electron.rate.output_data("erate.dat");
				hole.rate.output_data("hrate.dat");
				if (mid_gap_states){
					morphology.MG_electron_concentration.output_data("MG_minus.dat");
					morphology.MG_hole_concentration.output_data("MG_plus.dat");
				}
				if (morphology.ion_transport_activated()){
					negative_ion.concentration.output_data("negative_ion_con.dat");
					positive_ion.concentration.output_data("positive_ion_con.dat");
					negative_ion.current_x.output_data("negative_ion_currentx.dat");
					positive_ion.current_x.output_data("positive_ion_currentx.dat");
				}
				net_rate.output_data("rate.dat");
			
		}
		else if (measurement.experiment_CELIV || measurement.experiment_PhotoCELIV){

			measurement.set_current_data_point(0);

			std::cout << "Initiating calculations" << std::endl;

			PositionDependentParameter previous_electron_concentration, previous_hole_concentration, previous_neg_ion_concentration, previous_pos_ion_concentration, previous_potential;

			double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
			double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

			while (measurement.current_data_point < measurement.data_points){

				int i = measurement.get_current_data_point();

				MeasuredCurrent current;

				PositionDependentParameter net_rate;
				net_rate.initialize(electron.rate.normalization_coef);

				morphology.set_electrode_potential(0, measurement.anode_potential[i] - (morphology.get_work_function(0) - morphology.get_work_function(1)));
				morphology.set_electrode_potential(1, measurement.cathode_potential[i]);

				if (morphology.is_effective_temperature()){
					morphology.effective_temperature();
				}

				if (morphology.is_field_dependent_mobility() && i == 0){
					morphology.set_initial_field_dependent_mobility();
				}
				else if (morphology.is_field_dependent_mobility()){
					morphology.calculate_field_dependent_mobility(previous_potential);
				}

				if (measurement.experiment_PhotoCELIV && i == 1)
					morphology.light_on();
				else
					morphology.light_off();

				potential.set_boundary_conditions(morphology);
				if (!measurement.iterate_forward || i == 0)
					potential.set_initial_guess(morphology, "initial_guess\\potential", i);
				electron.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					electron.set_initial_guess(morphology, potential.electrical, "initial_guess\\electroncon", i);
				hole.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					hole.set_initial_guess(morphology, potential.electrical, "initial_guess\\holecon", i);
				if (!measurement.iterate_forward || i == 0)
					morphology.set_MG_initial_guess();
				if (morphology.ion_transport_activated()){
					negative_ion.set_ion_boundary_conditions(morphology);
					positive_ion.set_ion_boundary_conditions(morphology);
				}

				std::cout << "Boundary conditions and initial guess set." << std::endl;
				std::cout << "Solving for data point " << i << "..." << std::endl;
				
				if (i == 0){
					solver(morphology, potential, electron, hole, negative_ion, positive_ion, current, net_rate);
				}
				else{
					solver(morphology, potential, electron, hole, negative_ion, positive_ion, current, net_rate, 
						previous_hole_concentration, previous_electron_concentration, previous_neg_ion_concentration, previous_pos_ion_concentration, previous_potential, measurement.time_step);
				}
				if (measurement.output_full_data){
					potential.electrical.output_data("potential", i);
					electron.concentration.output_data("electroncon", i);
					hole.concentration.output_data("holecon", i);
					electron.current_x.output_data("ecurrentx", i);
					hole.current_x.output_data("hcurrentx", i);
					potential.displacement_current.output_data("dispcurrentx", i);
					potential.electrochemical_electron.output_data("electron_fermi_level", i);
					potential.electrochemical_hole.output_data("hole_fermi_level", i);
					electron.rate.output_data("erate", i);
					hole.rate.output_data("hrate", i);
					if (mid_gap_states){
						morphology.MG_electron_concentration.output_data("MG_minus", i);
						morphology.MG_hole_concentration.output_data("MG_plus", i);
					}
					if (morphology.ion_transport_activated()){
						negative_ion.concentration.output_data("negative_ion_con", i);
						positive_ion.concentration.output_data("positive_ion_con", i);
						negative_ion.current_x.output_data("negative_ion_currentx", i);
						positive_ion.current_x.output_data("positive_ion_currentx", i);
					}
					net_rate.output_data("rate", i);
				}
				std::ofstream fout("CELIV_transient.dat", std::ios::app);
				fout << measurement.time[i]*norm_time << "	" << std::scientific << - current.outer_circuit * current.normalization_coef << "	" << current.rms_error * current.normalization_coef << std::endl;
				fout.close();

				measurement.next_data_point();

				previous_electron_concentration = electron.concentration;
				previous_hole_concentration = hole.concentration;
				previous_potential = potential.electrical;
				if (morphology.ion_transport_activated()){
					previous_neg_ion_concentration = negative_ion.concentration;
					previous_pos_ion_concentration = positive_ion.concentration;
				}

			}

		}
		else if (measurement.experiment_transient_photocurrent){

			measurement.set_current_data_point(0);

			std::cout << "Initiating calculations" << std::endl;

			PositionDependentParameter previous_electron_concentration, previous_hole_concentration, previous_neg_ion_concentration, previous_pos_ion_concentration, previous_potential;

			double norm_potential = constant::boltzmann * constant::temperature / constant::elementary_charge;
			double norm_time = pow(constant::length_norm, 2) / (norm_potential*constant::mobility_norm);

			morphology.light_off();

			while (measurement.current_data_point < measurement.data_points){

				int i = measurement.get_current_data_point();

				MeasuredCurrent current;

				PositionDependentParameter net_rate;
				net_rate.initialize(electron.rate.normalization_coef);

				morphology.set_electrode_potential(0, measurement.anode_potential[i] - (morphology.get_work_function(0) - morphology.get_work_function(1)));
				morphology.set_electrode_potential(1, measurement.cathode_potential[i]);

				if (morphology.is_effective_temperature()){
					morphology.effective_temperature();
				}

				potential.set_boundary_conditions(morphology);
				if (!measurement.iterate_forward || i == 0)
					potential.set_initial_guess(morphology, "initial_guess\\potential", i);
				electron.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					electron.set_initial_guess(morphology, potential.electrical, "initial_guess\\electroncon", i);
				hole.set_boundary_conditions(morphology, potential.electrical);
				if (!measurement.iterate_forward || i == 0)
					hole.set_initial_guess(morphology, potential.electrical, "initial_guess\\holecon", i);
				if (!measurement.iterate_forward || i == 0)
					morphology.set_MG_initial_guess();
				if (morphology.ion_transport_activated()){
					negative_ion.set_ion_boundary_conditions(morphology);
					positive_ion.set_ion_boundary_conditions(morphology);
				}

				std::cout << "Boundary conditions and initial guess set." << std::endl;
				std::cout << "Solving for data point " << i << "..." << std::endl;

				if (i == 0){
					solver(morphology, potential, electron, hole, negative_ion, positive_ion, current, net_rate);
					morphology.light_on();
				}
				else{
					solver(morphology, potential, electron, hole, negative_ion, positive_ion, current, net_rate,
						previous_hole_concentration, previous_electron_concentration, previous_neg_ion_concentration, previous_pos_ion_concentration, previous_potential, measurement.time_step);
				}
				if (measurement.current_data_point == (int)(measurement.pulse_length / measurement.time_step))
					morphology.light_off();

				if (measurement.output_full_data){
					potential.electrical.output_data("potential", i);
					electron.concentration.output_data("electroncon", i);
					hole.concentration.output_data("holecon", i);
					electron.current_x.output_data("ecurrentx", i);
					hole.current_x.output_data("hcurrentx", i);
					potential.displacement_current.output_data("dispcurrentx", i);
					potential.electrochemical_electron.output_data("electron_fermi_level", i);
					potential.electrochemical_hole.output_data("hole_fermi_level", i);
					electron.rate.output_data("erate", i);
					hole.rate.output_data("hrate", i);
					if (mid_gap_states){
						morphology.MG_electron_concentration.output_data("MG_minus", i);
						morphology.MG_hole_concentration.output_data("MG_plus", i);
					}
					if (morphology.ion_transport_activated()){
						negative_ion.concentration.output_data("negative_ion_con", i);
						positive_ion.concentration.output_data("positive_ion_con", i);
						negative_ion.current_x.output_data("negative_ion_currentx", i);
						positive_ion.current_x.output_data("positive_ion_currentx", i);
					}
					net_rate.output_data("rate", i);
				}
				std::ofstream fout("transient_photocurrent.dat", std::ios::app);
				fout << measurement.time[i] * norm_time << "	" << std::scientific << -current.outer_circuit * current.normalization_coef << "	" << current.rms_error * current.normalization_coef << std::endl;
				fout.close();

				measurement.next_data_point();

				previous_electron_concentration = electron.concentration;
				previous_hole_concentration = hole.concentration;
				previous_potential = potential.electrical;
				if (morphology.ion_transport_activated()){
					previous_neg_ion_concentration = negative_ion.concentration;
					previous_pos_ion_concentration = positive_ion.concentration;
				}

			}

		}

	}



	void initialize_calculation(Morphology &morphology, Potential &potential, NegativeMobileCharge &electroncon, PositiveMobileCharge &holecon, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, Experiment &measurement){

		// potential.set_boundary_conditions(morphology, measurement.anode_potential[i], measurement.cathode_potential[i]);
		// electroncon.set_boundary_conditions(morphology, potential.electrical);
		// holecon.set_boundary_conditions(morphology, potential.electrical);

		// potential.set_initial_guess();
		// electroncon.set_initial_guess();
		// holecon.set_initial_guess();


		return;
	}


}

