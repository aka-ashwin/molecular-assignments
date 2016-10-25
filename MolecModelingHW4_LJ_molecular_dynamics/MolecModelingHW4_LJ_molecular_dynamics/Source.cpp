//COMMENTED OUT - test MD by checking conservation of total energy
////test: get energies over time for rho* = .5, T* = 2.


////initialize arrays
//posn positions_Etest[num_particles];
//posn velocities_Etest[num_particles];
//posn forces_Etest[num_particles];

////posn init_absolute_displacement[num_particles];
////posn absolute_displacement[num_particles];

//double KEs_Etest[num_prod_cycles_test];
//double PEs_Etest[num_prod_cycles_test];

////normalization decision - choose epsilon, sigma 'already' 1 in our units -> input r_cutoff, rho, T must be scaled accordingly
//double epsilon = 1.;
//double sigma = 1.;


////input parameters
//double r_cutoff = 3.;
//double rho_star = .5;
//double T_star = 2.;
//double dt = .005;
//double mass = 1.;

////initialize mass array from input mass - kind of clunky to use array in this case, but not terrible I think
//double masses[num_particles];
//for (int imass = 0; imass < num_particles; imass++) {
//	masses[imass] = mass;
//}


////derive lattice parameters, initialize lattice
//double V = num_particles/ rho_star;
//double sidelen = pow(V, 1 / 3.);
//double boxdims[3] = { sidelen, sidelen, sidelen };

//double approx_parts_per_side = pow(num_particles, 1 / 3.);
//int num_parts_per_side = ceil(approx_parts_per_side);
//
////use ceil so that there are enough lattice positions if num_particles is not a cube

////if num particles is perfect cube, make sure rounding error doesn't add extra particle per side:
//int alt_num_parts_per_side = floor(approx_parts_per_side);
//if (alt_num_parts_per_side*alt_num_parts_per_side*alt_num_parts_per_side == num_particles) {
//	num_parts_per_side = alt_num_parts_per_side;
//}

//fill_simple_symm_cubic_lattice(positions_Etest, num_parts_per_side, num_particles, sidelen);


////initialize velocities
//fill_Boltzmann_velocities(velocities_Etest, num_particles, random, masses, T_star * epsilon);

////loop until desired number of steps (incl initialization step) have been completed
//for (int i = 0; i < num_prod_cycles_test; i++) {
//	//update energy arrays
//	PEs_Etest[i] = LJ_energy_w_cutoff(positions_Etest, num_particles, true, epsilon, sigma, boxdims, r_cutoff, false);
//	KEs_Etest[i] = get_KE(velocities_Etest, masses, num_particles);

//	//write lattice to file
//	print_lattice_to_file("test_positions", positions_Etest, num_particles, "x, y, z", true, "xyz");

//	//time-evolve system if not at last half-step
//	if (i == num_prod_cycles_test - 1) { continue; }
//	else {
//		update_Verlet(positions_Etest, velocities_Etest, num_particles, boxdims, forces_Etest, masses, false, default_absolute_disp_array, true, 3.0, epsilon, sigma, dt);
//	}
//}

//
////gather KE/N, PE/N, E/N into single big array
////kind of kludgy but I'm using an array of posns (structs w 3 doubles each) for this - easier to write to file w current methods
//posn energies[num_prod_cycles_test];
//for (int i = 0; i < num_prod_cycles_test; i++) {
//	energies[i].x = KEs_Etest[i] / (double) num_particles;
//	energies[i].y = PEs_Etest[i] / (double) num_particles;
//	energies[i].z = energies[i].x + energies[i].y;
//}

//
////write KE, PE, E to file

//string energy_filename = "energies_from_LJ_MD_test_fixed_" + to_string(num_particles) +"_particles";
//string E_header = "KE, PE, E";
//
//
//print_lattice_to_file(energy_filename, energies, num_prod_cycles_test, E_header, true, "tsv");



//calculate diffusion constant via <(delta r^2)(t)>



//test: get energies over time for rho* = .5, T* = 2.
