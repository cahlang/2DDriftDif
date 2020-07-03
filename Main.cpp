#include <iostream>
#include <mpi.h>
#include <boost/property_tree/ptree.hpp>

#include "Control.h" // Contains functions for running various experiments. Use them in order to keep main() organized and clear. The functions are wrapped in a namespace, syntax: control::function().
#include "DataTypes.h"
#include "Potential.h"
#include "NegativeMobileCharge.h"
#include "PositiveMobileCharge.h"
#include "Morphology.h"

int main(int argc, char *argv[]){

	// Configure parallelization, division of the lattice and ghost exchange must be implemeted!
	// MPI_Init(&argc, &argv);

	std::string file_name = "Settings.ini";

	Morphology morphology;

	Potential potential;
	NegativeMobileCharge electron;
	PositiveMobileCharge hole;

	NegativeMobileCharge negative_ion;
	PositiveMobileCharge positive_ion;

	Experiment experiment;

	if (control::read_settings(file_name, morphology, potential, electron, hole, negative_ion, positive_ion, experiment)){
		printf("Settings read\n");
	}

	control::run_measurement(morphology, potential, electron, hole, negative_ion, positive_ion, experiment);



	// MPI_Finalize();

	return 0;
}