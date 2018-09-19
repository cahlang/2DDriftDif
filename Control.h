
#pragma once

#include "DataTypes.h"

#include <string>

class Morphology;
class Potential;
class NegativeMobileCharge;
class PositiveMobileCharge;
class StationaryCharge;

// This namespace contains functions for performing sets of tasks, such as reading all input and distributing the parameters to their appropriate variables.

namespace control{

	// Initializes all vectors and reads all settings.
	bool read_settings(const std::string &file_name, Morphology &morphology, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, Experiment &measurement);

	// Calculates the electron and hole concentrations, population of CT states and the electrical potential for the next iterative step.
	bool solver(Morphology &material, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, MeasuredCurrent &outer_circuit_current, PositionDependentParameter &net_rate);
	bool solver(Morphology &material, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, MeasuredCurrent &outer_circuit_current, PositionDependentParameter &net_rate,
		 PositionDependentParameter &previous_holecon, PositionDependentParameter &previous_electroncon, PositionDependentParameter &previous_neg_ion_con, PositionDependentParameter &previous_pos_ion_con, PositionDependentParameter &previous_potential, double time_step);

	// Set boundary conditions and initial guess.
	void initiate_calculation(Morphology &material, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, MeasuredCurrent &outer_circuit_current, PositionDependentParameter &net_rate, Experiment &measurement);

	// Calculates the generation and recombination rates for the given concentrations and electrical potential.
	void determine_rates(Morphology &morph, PositionDependentParameter &net_rate, NegativeMobileCharge &electron, PositiveMobileCharge &hole, PositionDependentParameter &potential);

	// Simulates the selected measurement. Currently only IV-curves are implemented, but it is quite simple to build more complex ones once the time dependence is implemented.
	void run_measurement(Morphology &morphology, Potential &potential, NegativeMobileCharge &electron, PositiveMobileCharge &hole, NegativeMobileCharge &negative_ion, PositiveMobileCharge &positive_ion, Experiment &measurement);


}