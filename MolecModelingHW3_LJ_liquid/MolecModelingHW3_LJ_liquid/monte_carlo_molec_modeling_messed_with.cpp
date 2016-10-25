#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <ctime>
#include "C:/Users/Ashwin/Desktop/Molec Modeling/random_mars.h"
#include <math.h>
//#include <chplot.h>
using namespace std;



// need:
// Part 1 : to write initial lattice to .dat file
// and make an image of it - use VMD?




// Part 1: initialize and print 
//initialize 3D lattice; will store points as array of double arrays (= a double**), with their x,y,z positions

//dumb global initialization - remove
struct posn {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	//initialization function
	posn(double xin = 0.0, double yin = 0.0, double zin = 0.0) {
		x = xin;
		y = yin;
		z = zin;
	}
	//setting function
	void set_all(double xin, double yin, double zin) {
		x = xin;
		y = yin;
		z = zin;
	}
	void scale(double scalefactor) {
		x = x*scalefactor;
		y = y*scalefactor;
		z = z*scalefactor;
	}

	//get length = (x^2 + y^2 + z^2)^(1/2) of posn
	double length() {
		return sqrt(x*x + y*y + z*z);
	}

	//normalize vector (divide by length) - return unit vector in same direction
	void normalize() {
		scale(length());
	}

	
};

posn add_posns(posn first, posn second, double first_scaling = 1.0, double second_scaling = 1.0) {
	double xsum = first_scaling * first.x + second_scaling * second.x;
	double ysum = first_scaling * first.y + second_scaling * second.y;
	double zsum = first_scaling * first.z + second_scaling * second.z;
	posn output_posn = posn(xsum, ysum, zsum);
	return output_posn;
}

double dot_posns(posn first, posn second, double scale = 1.0) {
	return scale*(first.x*second.x + first.y*second.y + first.z*second.z);
}



//fill in simple cubic lattice (3D)
//starts at 0.0 for each dimension
//assumes input lattice of num_cols*num_rows*num_planes particles
//note: MODIFIES INPUT LATTICE
void fill_simple_cubic_lattice(posn* input_lattice, int num_cols, int num_rows, int num_planes, int num_particles = 0,  double boxlen_x = 1.0, double boxlen_y = 1.0, double boxlen_z = 1.0) {
	// lattice parameters
	double xstep = boxlen_x / num_cols;
	double ystep = boxlen_y / num_rows;
	double zstep = boxlen_z / num_planes;

	if (num_particles == 0) { num_particles = num_cols*num_rows*num_planes; }

	// loop through and fill lattice
	int lattice_index = 0;
	double curr_x = -1*xstep;
	double curr_y = -1*ystep;
	double curr_z = -1*zstep;

	for (int row_index = 0; row_index < num_rows; row_index++) {
		curr_x += xstep;

		for (int col_index = 0; col_index < num_cols; col_index++) {
			curr_y += ystep;

			for (int plane_index = 0; plane_index < num_planes; plane_index++) {
				curr_z += zstep;
				input_lattice[lattice_index].set_all(curr_x, curr_y, curr_z);
				lattice_index++;
				//if have already filled input lattice (to desired amount), end loop!
				if (lattice_index >= num_particles) { break; }
			}
			curr_z = -zstep;
		}
		curr_y = -ystep;
		
	}
}

void fill_simple_symm_cubic_lattice(posn* input_lattice, int num_parts_per_side, int num_particles = 0, double boxlen= 1.0) {
	fill_simple_cubic_lattice(input_lattice, num_parts_per_side, num_parts_per_side, num_parts_per_side, num_particles, boxlen, boxlen, boxlen);
}



// print lattice: currently just prints to console

void print_lattice_to_console(posn* positions, int num_particles) {
	for (int i = num_particles - 1; i >= 0; i--) {
		cout << i;
		cout << ": ";
		cout << positions[i].x;
		cout << ",";
		cout << positions[i].y;
		cout << ",";
		cout << positions[i].z;
		cout << "\n";
	}
	//code to stall console at this point
	//std::cout << "Press ENTER to continue...";
	//std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

}

void print_double_array_to_file(string filename, double* input_array, int num_entries, string filehead = "", bool has_header = false, string filetype = "csv") {
	string separator = "\t";
	if (filetype == "csv") { separator = ", "; }
	
	string filepath = "outputs/";
	string full_filename = filepath + (filename + ("." + filetype));
	ofstream myfile;
	myfile.open(full_filename, ios::app);
	//if has header, include header line (blank by default)
	if (has_header) {
		myfile << filehead + "\n";
	}
	//loop through particles and write their coordinates, separated by desired separator
	//
	int currInd;
	for (int i = num_entries - 1; i >= 0; i--) {
		currInd = num_entries - (i + 1);
		myfile << input_array
		myfile << "\n";
	}
	myfile.close();


}

void print_lattice_to_file(string filename, posn* positions, int num_particles, string filehead = "index, x, y, z", bool has_header = false, string filetype = "csv", string atomtype_if_xyz = "Ar") {
	string separator = "\t";
	if (filetype == "csv") { separator = ", "; }

	string filepath = "outputs/";
	string full_filename = filepath + (filename + ("." + filetype));
	ofstream myfile;
	myfile.open(full_filename, ios::app);
	//if xyz, print number of particles first
	if (filetype == "xyz") { myfile << num_particles; myfile << "\n"; }
	//if has header, of if is in xyz format which requires a header line, include header line (blank by default)
	if (has_header | filetype == "xyz") {
		myfile << filehead + "\n";
	}
	//loop through particles and write their coordinates, separated by desired separator
	//
	int currInd;
	for (int i = num_particles - 1; i >= 0; i--) {
		currInd = num_particles - (i + 1);
		//if xyz, include atom type
		myfile << atomtype_if_xyz;
		myfile << separator;
		//then write atom coordinates, spaced by separator 
		myfile << positions[currInd].x;
		myfile << separator;
		myfile << positions[currInd].y;
		myfile << separator;
		myfile << positions[currInd].z;
		myfile << "\n";
	}
	myfile.close();


}



// Part 2: calculate energy

//helper function: given particle i and j's positions, gets distance  particle i to nearest image of particle j
//assumes box length in x and y dims is 1

//calculate 3d vector length (or 2d, if z is not specified)
double vector_length(double x, double y, double z = 0.0) {
	return sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
}

//calculate length of a posn structure
double posn_length(posn input_posn) {
	double x = input_posn.x;
	double y = input_posn.y;
	double z = input_posn.z;
	return vector_length(x, y, z);
}

//apply periodic boundary conditions for a single dimension of position; keeps position within [0,boxlen)
double apply_pbcs_1d(double x, double boxlen = 1.0) {
	double pbc_x = x;

	//use of 'while' rather than 'if' allows for proposed x to be more than 1 boxlen out of the box. 
	while (pbc_x >= boxlen) { pbc_x -= 1.0*boxlen; }
	while(pbc_x < 0.0) { pbc_x += 1.0*boxlen; }

	return pbc_x;
}


double default_boxdims[3] = { 1.0, 1.0, 1.0 };

//apply periodic boundary conditions in 3D
posn apply_pbcs_3d(posn input_vector, double* boxdims = default_boxdims) {
	double pbc_x = apply_pbcs_1d(input_vector.x, boxdims[0]);
	double pbc_y = apply_pbcs_1d(input_vector.y, boxdims[1]);
	double pbc_z = apply_pbcs_1d(input_vector.z, boxdims[2]);

	posn pbc_vec = posn(pbc_x, pbc_y, pbc_z);
	return pbc_vec;
}


double minimum_image_distance_1d(double xfrom, double xto, double boxlen = 1.0) {
	double dist_1d = xto - xfrom;
	double min_dist = dist_1d;

	//distance can be in [-.5 boxlen, .5 boxlen)

	if (dist_1d > .5*boxlen) { min_dist -= 1.0*boxlen; }
	else if (dist_1d < -.5*boxlen) { min_dist += 1.0*boxlen; }

	return min_dist;
}



double minimum_image_distance(posn particlefrom, posn particleto, double* boxdims = default_boxdims) {
	
	double xfrom = particlefrom.x;
	double yfrom = particlefrom.y;
	double zfrom = particlefrom.z;
	
	double xto = particleto.x;
	double yto = particleto.y;
	double zto = particleto.z;

	double xboxlen = boxdims[0];
	double yboxlen = boxdims[1];
	double zboxlen = boxdims[2];


	double min_xdist = minimum_image_distance_1d(xfrom, xto, xboxlen);
	double min_ydist = minimum_image_distance_1d(yfrom, yto, yboxlen);
	double min_zdist = minimum_image_distance_1d(zfrom, zto, zboxlen);

	//note {min_xdist, min_ydist, min_zdist} = minimum-image version of the vector r_to - r_from

	double minim_distance = vector_length(min_xdist, min_ydist, min_zdist);
	return minim_distance;
}

posn minimum_image_form_of_displacement(posn input_displacement, double*box_dims = default_boxdims) {
	posn minim_displacement = input_displacement;

	minim_displacement.x = minimum_image_distance_1d(0., minim_displacement.x, box_dims[0]);
	minim_displacement.y = minimum_image_distance_1d(0., minim_displacement.y, box_dims[1]);
	minim_displacement.z = minimum_image_distance_1d(0., minim_displacement.z, box_dims[2]);

	return minim_displacement;
}

posn minimum_image_displacement_posns(posn particle_from, posn particle_to, double* box_dims = default_boxdims) {
	return minimum_image_form_of_displacement(add_posns(particle_to, particle_from, 1.0, -1.0), box_dims);

}

//hard-sphere distance calculations

bool no_hardsph_overlap_w_i(posn* positions, int num_particles, double sph_radius, int i) {
	bool no_overlap = true;

	posn i_posn = positions[i];
	posn j_posn;
	double ij_dist;

	for (int j = 0; j < num_particles; j++) {
		if (j != i) {
			j_posn = positions[j];
			ij_dist = minimum_image_distance(i_posn, j_posn);
			if (ij_dist < sph_radius*2.0) {
				no_overlap = false;
				break;
			}
		}
	}
	return no_overlap;
}



//Lennard-Jones energy- and force-calculation functions
// if has no cutoff, r_cutoff is not used, so it is given a default value for this use case
// energy does not include U_tail = (8./3) pi rho*((1/3)r_cutoff^(-9) - r_cutoff^(-3))
//	because this requires additional input of rho, is inefficient for our fixed-rho sims, and is irrelevant for our Monte Carlo steps.

//note: by default uses LJ coordinates with epsilon, sigma = 1

//LJ force calculation
//force on i from j = 4*epsilon * (12 sigma^12 r_ij^(-13) - 6 sigma^6 r_ij^(-7)) r_ji hat
//where r_ji hat = unit vector in (r_i - r_j direction)
//can be given (r/sigma)^(-6), (r/sigma)
posn LJ_force_j_on_i(posn part_i, posn part_j, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2, bool scaled_distance_inv_and_pow_given = false, double scaled_distance_inv = 1.0, double scaled_distance_invsixthpow = 1.0) {
	posn force_vec = posn(0., 0., 0.);

	//get vector between positions 
	posn j_minus_i_minim = minimum_image_displacement_posns(part_j, part_i, boxdims);

	if (!scaled_distance_inv_and_pow_given) {
		scaled_distance_inv = sigma/j_minus_i_minim.length();
		scaled_distance_invsixthpow = pow(scaled_distance_inv, 6.0);
	}


	//if distance less than cutoff radius, calculate force
	if ((1./scaled_distance_inv) < scaled_r_cutoff) {

		//magnitude of force = 4 * epsilon * [term]
		//[term] = 12 sigma^12/r^13 - 6 sigma^6 / r^7
		// = 6 sigma^6/r^7 (2 sigma^6/r^6 - 1)
		// = 6 (r*^(-6)) * (1/r*)/sigma * (2r*^(-6) - 1)
		double F_magnitude = 4 * epsilon;
		F_magnitude *= 6 * scaled_distance_invsixthpow * (scaled_distance_invsixthpow*scaled_distance_inv/sigma) * (2 * scaled_distance_invsixthpow - 1.0);

		//set force vector equal to unit vector in r_j - r_i direction
		force_vec = j_minus_i_minim;
		force_vec.normalize();
		
		//scale force vector by magnitude of force
		//note that sign of magnitude is correct because it is large (pushing j away from i) when distance is small, and vice-versa.
		force_vec.scale(F_magnitude);


	}
	return force_vec;
}

posn default_position_array[1] = { posn(0,0,0) };



double LJ_energy_w_cutoff_ithpart(posn* positions, int num_particles, bool has_cutoff, int i, bool use_below_only = true, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, double r_cutoff = .2, bool calc_force_too = false, posn* force_array = default_position_array) {
	//note energy starts at 0 rather than U_tail
	double energy = 0.0;
	posn net_force = posn(0.0, 0.0, 0.0);

	posn i_posn = positions[i];
	posn j_posn;

	//initialize variables for use in loop
	double scaled_dist_inv;
	double scaled_dist_inv6th;
	double U_ij;

	//set upper bound for loop over other particles' indices: = i by default, = number of particles otherwise
	int j_ub_notreached = i;
	if (!use_below_only) { j_ub_notreached = num_particles; }


	// loop through all desired j-values
	for (int j = j_ub_notreached - 1; j >= 0; j--) {
		//if j == i, skip this iteration
		if (j == i) { continue; }

		//else, get j's position...
		j_posn = positions[j];
		//note that minimum_image_displacement always gives the vector difference i_posn - j_posn
		double min_distance = minimum_image_displacement_posns(i_posn, j_posn, boxdims).length();
		//if no cutoff, or if have cutoff and distance is less than it..
		if (!has_cutoff | (has_cutoff & (min_distance < r_cutoff))) {
			//calculate pairwise Lennard-Jones energy
			scaled_dist_inv = sigma / min_distance;
			scaled_dist_inv6th = pow(scaled_dist_inv, 6.);
			U_ij = 4 * epsilon*scaled_dist_inv6th*(scaled_dist_inv6th - 1.0);
			energy += U_ij;
			if (calc_force_too) {
				net_force = add_posns(net_force, LJ_force_j_on_i(i_posn, j_posn, epsilon, sigma, boxdims, has_cutoff, r_cutoff, true, scaled_dist_inv, scaled_dist_inv6th));			
			}
		}
	}

	//if calculating force, updated desired element of force array
	if(calc_force_too){ force_array[i] = add_posns(force_array[i], net_force); }

	return energy;
}



double LJ_energy_w_cutoff(posn* positions, int num_particles, bool has_cutoff, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, double r_cutoff = .2, bool calc_force_too = false, posn* force_array = default_position_array) {
	//note energy starts at 0 rather than U_tail
	double energy = 0.0;

	//initialize variables for use in loop

	// loop through all i values > 1; note we don't need i = 0, since the 1-0 interaction is covered by i = 1, j = 0.
	// get energy of i's interactions with all particles j for j < i; add to total energy
	for (int i = num_particles - 1; i > 0; i--) {
		energy += LJ_energy_w_cutoff_ithpart(positions, num_particles, has_cutoff, i, true, epsilon, sigma, boxdims, r_cutoff, calc_force_too, force_array);
	}
	return energy;
}

//Part 3 - Monte Carlo implementation

//propose a step: select an atom uniformly at random, generate e_x and e_y uniformly at random from [-1,1]
//and move it by dx = (d - d_0)*ex and dy = (d - d_0)*ey, where d = 1/14 and d_0 is an input value
//respecting periodic boundary conditions

//hard sphere MC step
//note- assumes box is 1 x 1 x 1
void MC_update_hardsph(posn* positions, int num_particles, int part_index, posn rand_move, double sph_radius) {
	posn old_posn = positions[part_index];

	posn proposed_posn = add_posns(positions[part_index], rand_move);
	proposed_posn = apply_pbcs_3d(proposed_posn);

	positions[part_index] = proposed_posn;

	bool valid_move = no_hardsph_overlap_w_i(positions, num_particles, sph_radius, part_index);

	if (!valid_move) {
		positions[part_index] = old_posn;
	}

}

void MC_cycle_hardsph(posn* positions, int num_particles, RanMars* rng, double sph_radius_d0, double scalefac= 1.0) {
	for (int i = num_particles; i > 0; i--) {

		double x_move = .5 * ((2 * rng->uniform()) - 1);
		double y_move = .5 * ((2 * rng->uniform()) - 1);
		double z_move = .5 * ((2 * rng->uniform()) - 1);

		posn proposed_move = posn(x_move, y_move, z_move);

		proposed_move.scale(scalefac);

		//choose particle to move
		//generate uniform double between 0 and num_particles; 
		//for each particle i from 0 to num_particles-1 there is probability 1/num_particles of being between i and i + 1
		//so casting the result to an int gives us a uniform integer between 0 and num_particles
		int part_index = (int)(num_particles * rng->uniform());

		MC_update_hardsph(positions, num_particles, part_index, proposed_move, sph_radius_d0);
	}
}

void MC_loop_hardsph(posn* positions, int num_particles, double density, int num_prod_cycles, double sph_radius_d_0, bool write_to_file, string filename = "hardsph") {
	int	num_equil_cycles = 5;
	


	//initialize RNG
	RanMars *random;
	//seed with time
	unsigned int seed = unsigned int(time(NULL));
	int real_seed = abs((int)seed) % 900000000;
	random = new RanMars(real_seed);


	//initialize lattice

	//get n = # of particles per lattice side; want cubic lattice s.t. n = N^(1/3) (or smallest integer above this)
	int lattice_side_numparts = ceil(pow((double)num_particles, 1. / 3.));


	//set box-dimension array = num_particles_per_side*[lattice_param, lattice_param, lattice_param]
	double volume = num_particles / density;
	double lattice_param = pow(volume, 1. / 3.);
	double box_sidelen = lattice_param * lattice_side_numparts;
	double box_dimensions[3] = { box_sidelen, box_sidelen, box_sidelen };

	//fill provided positions array with particles initialized on a simple cubic lattice fitting the input params
	//note that lattice must be initialized before this function is called
	fill_simple_symm_cubic_lattice(positions, lattice_side_numparts, num_particles, lattice_param);

	//for now, will scale proposal moves by .1 * volume^(1/3) factor
	double scalefac = .0000001*lattice_param;

	//equilibration cycles 
	for (int i = 0; i < num_equil_cycles; i++) {
		MC_cycle_hardsph(positions, num_particles, random, sph_radius_d_0, scalefac); 
		if (write_to_file) {
			print_lattice_to_file(filename, positions, num_particles, "equilibration cycle", true, "xyz", "Ar");
		}
	}

	//production cycles
	for (int j = 0; j < num_prod_cycles; j++) {
		MC_cycle_hardsph(positions, num_particles, random, sph_radius_d_0, scalefac);
		if (write_to_file) {
			print_lattice_to_file(filename, positions, num_particles, "production cycle", true, "xyz", "Ar");
		}
	}


}


//constants used in MC step: parameters of energy equation, plus input random move
//will leave random number generation to MC cycle function (and seeding to MC loop function)
void MC_update_LJ(posn* positions, int part_index, posn rand_move, double acceptance_threshold, double* box_dims, double temperature_scaled, int num_particles, bool LJ_cutoff, double epsilon = 1.0, double sigma = 1.0, double r_cutoff = .2) {
	
	posn old_part = positions[part_index];

	posn proposed_posn = add_posns(old_part, rand_move);

	//apply pbcs; note that dimensions can equal 0 but not the box length
	//potential issue: assumes proposal displacements are small enough that they never move particles more than a box length
	//a reasonable assumption here and for most realistic displacements, but not good to have as a general case.
	

	proposed_posn = apply_pbcs_3d(proposed_posn, box_dims);

	//get difference in energy, = difference in particle i's total interaction energy
	double init_i_energy = LJ_energy_w_cutoff_ithpart(positions, num_particles, LJ_cutoff, part_index, false, epsilon, sigma, box_dims, r_cutoff);
	//update i's position and get its new interaction energies
	positions[part_index] = proposed_posn;
	double final_i_energy = LJ_energy_w_cutoff_ithpart(positions, num_particles, LJ_cutoff, part_index, false, epsilon, sigma, box_dims, r_cutoff);

	//get delta_E = E_f - E_i
	double delta_E = final_i_energy - init_i_energy;

	//if delta_E <= 0, pi factor is greater than 1 -> always accept 
	// in this case nothing further needs to be done, since positions[i] is already the proposal position

	//if delta_E > 0, have to calculate acceptance step:

	if (delta_E > 0.0) {
		//boltzmann exponent = delta_E/kT; scaled_temperature is T* = kT/epsilon, so kT = epsilon times T*, so boltz exponent = delta_E / (epsilon T*)
		double boltz_exponent = -1. * delta_E / (epsilon*temperature_scaled);
		double boltz_factor = exp(boltz_exponent);

		//compare boltzmann factor to threshold X (~ Unif(0,1)); if >= X, accept and do nothing; otherwise reject and revert
		if (boltz_factor < acceptance_threshold) {
			positions[part_index] = old_part;
		}
	}
}


//box_dims, temperature_scaled, num_particles, LJ_cutoff, epsilon, sigma, r_cutoff);

// Monte Carlo cycle: run update process n times, where n is the number of particles
void MC_cycle_LJ(posn* positions, int num_particles, RanMars *rng, double* box_dims, double temperature_scaled, bool use_LJ_cutoff, double epsilon, double sigma, double r_cutoff = .2) {
	for (int i = num_particles - 1; i > 0; i--) {
		//generate x and y move RVs; RanMars.uniform gives unif(0,1) and we want unif(-1,1), so we multiply by 2 and subtract 1.
		//for our hard-sphere model, we relied on scaling the move in the individual MC steps
		//here, we will do scaling here, by a fixed ratio .5*relevant box dim, since the ranges from [.5,1] and [-1,-.5] in box-side units are made redundant by PBCs

		double x_move = .5*box_dims[0] * ((2 * rng->uniform()) - 1);
		double y_move = .5*box_dims[1] * ((2 * rng->uniform()) - 1);
		double z_move = .5*box_dims[2] * ((2 * rng->uniform()) - 1);

		posn proposed_move = posn(x_move, y_move, z_move);
		

		//choose particle to move
		//generate uniform double between 0 and num_particles; 
		//for each particle i from 0 to num_particles-1 there is probability 1/num_particles of being between i and i + 1
		//so casting the result to an int gives us a uniform integer between 0 and num_particles
		int part_index = (int)(num_particles * rng->uniform());

		//finally, generate acceptance threshold uniformly between 0,1
		double accept_threshold = rng->uniform();


		MC_update_LJ(positions,part_index, proposed_move, accept_threshold, box_dims, temperature_scaled, num_particles, use_LJ_cutoff, epsilon, sigma, r_cutoff);
	}
}







//Parts 3 on: running loops

double default_input_double_array[1] = { 0.0 };

// Monte Carlo loop: initialize lattice, run a few equilibration cycles, then run production cycles
void MC_loop_LJ(int num_prod_cycles, posn* positions, int num_particles, double density_scaled, double temperature_scaled, double sigma, double epsilon, bool scale_to_LJ_params = true, bool use_LJ_cutoff = true, double r_cutoff_LJ_scaled = .2, bool print_all_positions_to_file = true, string posn_filename = "irrelevant", bool calc_energies = true, double* energies = default_input_double_array, bool calc_force_too = false, posn* force_array = default_position_array) {
	int	num_equil_cycles = 5;

	//initialize RNG
	RanMars *random;
	//seed with time
	unsigned int seed = unsigned int(time(NULL));
	seed = abs((int)seed)% 900000000;
	random = new RanMars(seed);


	//initialize lattice

	//get n = # of particles per lattice side; want cubic lattice s.t. n = N^(1/3) (or smallest integer above this)
	int lattice_side_numparts = ceil(pow((double) num_particles, 1. / 3.));
	
	//get lattice parameter alpha: want density = N/(alpha^3), so alpha = (N/density)^(1/3)
	//density scaled = density*sigma^3, so we'll get back regular density first
	// (there's probably a less clunky way to do this, but I didn't want to mess around with factors of sigma
	// and wanted to leave open the ability to use non-scaled coordinates)
	double density = density_scaled / (pow(sigma, 3.));
	double lattice_param = pow(num_particles / density, 1. / 3.);

	//if using scaled coordinates x* = x/sigma, etc., then scale lattice param by sigma
	//and pass sigma* = 1.0, epsilon* = 1.0 to cycling functions
	double sigma_to_pass = sigma;
	double epsilon_to_pass = epsilon;
	if (scale_to_LJ_params) { 
		lattice_param *= 1. / (sigma); 
		sigma_to_pass = 1.0;
		epsilon_to_pass = 1.0;
	}
	//also pass appropriately scaled cutoff radius
	double r_cutoff_to_pass = r_cutoff_LJ_scaled * sigma;
	if (scale_to_LJ_params) { r_cutoff_to_pass = r_cutoff_LJ_scaled; }

	//set box-dimension array = num_particles_per_side*[lattice_param, lattice_param, lattice_param]
		//gives side lengths of box in units we are using--either scaled or unscaled ones
	double box_sidelen = lattice_param * lattice_side_numparts;
	double box_dimensions[3] = { box_sidelen, box_sidelen, box_sidelen };

	//fill provided positions array with particles initialized on a simple cubic lattice fitting the input params
	//note that lattice must be initialized before this function is called
	fill_simple_symm_cubic_lattice(positions, lattice_side_numparts, num_particles, lattice_param);

	//equilibration cycles 
	for (int i = 0; i < num_equil_cycles; i++) {
		MC_cycle_LJ(positions, num_particles, random, box_dimensions, temperature_scaled, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, r_cutoff_to_pass);
		if (print_all_positions_to_file) {
			print_lattice_to_file(posn_filename, positions, num_particles, "equilibration cycle", true, "xyz", "Ar");
		}
	}

	//production cycles
	for (int j = 0; j < num_prod_cycles; j++) {
		MC_cycle_LJ(positions, num_particles, random, box_dimensions, temperature_scaled, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, r_cutoff_to_pass);
		
		if (print_all_positions_to_file) {
			print_lattice_to_file(posn_filename, positions, num_particles, "production cycle", true, "xyz", "Ar");
		}


		//if calculating energies, get current energy and put it in the energies array argument
		if (calc_energies) {
			energies[j] = LJ_energy_w_cutoff(positions, num_particles, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, box_dimensions, r_cutoff_to_pass, calc_force_too, force_array);
		}
	}
}



// 
void main() {
	constexpr int num_prod_cycles_hardsph = 100;
	constexpr int num_prod_cycles_LJ = 100;

	constexpr int num_particles = 1000;
	
	//still working on getting lattice initialization done in a way that allows for user input to console
	//(I know there are some fancy non-fixed array-like things ('vectors'?) but I don't want to use them until I'm sure I can make them work)
	//(And they seem inefficient for our current purposes, though obviously necessary when we start changing N)
	//failed code below
	//cout << "Please enter the number of particles desired";
	//cin >> num_particles;



	//hard sphere cycles
	
	//cout << "Running hard spheres.\n";
	//posn positions_sph[num_particles];

	////test: run hard spheres in 3D and get xyz file
	//MC_loop_hardsph(positions_sph, num_particles, 1000.0, num_prod_cycles_hardsph, .001, true, "testloop_posns_hard_sphere");

	//cout << "done!";

	
	//Lennard-Jones cycles

	//first- LJ test and energy grabbing

	//run w scaled density = 1., scaled temperature = 2.
	//note we will print all lattices to file for visualization
	//and calculate energy, but not force just yet
	
	//initialize position and energy list
	double energies[num_prod_cycles_LJ];
	posn positions_LJ[num_particles];

	double scaled_r_cutoff = 3.0;
	MC_loop_LJ(num_prod_cycles_LJ, positions_LJ, num_particles, 1., 2., 1.0, 1.0, true, true, scaled_r_cutoff, true, "testloop_posns_LJ", true, energies);
	
	int i;
	for (i = 0; i < num_prod_cycles_LJ; i++) {
		cout << energies[i];
		cout << "\n";
	}
	string dummy_input;
	cin >> dummy_input;
}