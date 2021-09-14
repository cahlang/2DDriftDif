#pragma once

#include "DataTypes.h"

class Morphology;

namespace calc{

	// Calculates the value of the so-called Bernoulli equation for the given value.
	double bernou(double x);

	// Calculates the current measured in an outer circuit (for example that observed in an IV-measurement) and the RMS error.
	MeasuredCurrent outer_circuit_current(Morphology &device_parameters, const PositionDependentParameter &electron_current, const PositionDependentParameter &hole_current, const PositionDependentParameter &displacement_current, const PositionDependentParameter &net_rate, const PositionDependentParameter &recombination_current_x, const PositionDependentParameter& recombination_current_y);

}