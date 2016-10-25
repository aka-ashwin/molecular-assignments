#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <ctime>
#include "C:/Users/Ashwin/Desktop/Molec Modeling/random_mars.h"
#include <math.h>
//#include <chplot.h>
using namespace std;




// Part 1: initialize and print 
//initialize 3D lattice; will store points as array of posns (struct with 3 doubles: x,y,z)

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
		double len = sqrt(x*x + y*y + z*z);
		x = x / len;
		y = y / len;
		z = z / len;
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

void clear_posn_array(posn* array, int num_posns) {
	for (int i = 0; i < num_posns; i++) {
		array[i].set_all(0., 0., 0.);
	}
}






//--

struct intposn {
	int x = 0;
	int y = 0;
	int z = 0;

	//initialize
	intposn(int xin = 0, int yin = 0, int zin = 0) {
		x = xin;
		y = yin;
		z = zin;
	}

	//set
	void set_all(int xin, int yin, int zin) {
		x = xin;
		y = yin;
		z = zin;
	}


	int* to_array() {
		int arr[3] = { x,y,z };
		return arr;
	}

	void from_array(int* arr_in) {
		x = arr_in[0];
		y = arr_in[1];
		z = arr_in[2];
	}

};


//fill in simple cubic lattice (3D)
//starts at 0.0 for each dimension
//assumes input lattice of num_cols*num_rows*num_planes particles
//note: MODIFIES INPUT LATTICE
void fill_simple_cubic_lattice(posn* input_lattice, int num_cols, int num_rows, int num_planes, int num_particles = 0, double boxlen_x = 1.0, double boxlen_y = 1.0, double boxlen_z = 1.0) {
	// lattice parameters
	double xstep = boxlen_x / num_cols;
	double ystep = boxlen_y / num_rows;
	double zstep = boxlen_z / num_planes;

	if (num_particles == 0) { num_particles = num_cols*num_rows*num_planes; }

	// loop through and fill lattice
	int lattice_index = 0;
	double curr_x = -1 * xstep;
	double curr_y = -1 * ystep;
	double curr_z = -1 * zstep;

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

void fill_simple_symm_cubic_lattice(posn* input_lattice, int num_parts_per_side, int num_particles = 0, double boxlen = 1.0) {
	fill_simple_cubic_lattice(input_lattice, num_parts_per_side, num_parts_per_side, num_parts_per_side, num_particles, boxlen, boxlen, boxlen);
}


void fill_Boltzmann_velocities(posn* velocities, int num_particles, RanMars* rng, double* masses, double kT) {
	//at temperature T, velocity j-component Boltzmann distrib = Gaussian(0, kT/m)
	// can draw from Gaussian(0,1) and multiply by kT/m to get the proper value

	double total_momentum_first_assigned[3] = { 0.,0.,0., };

	double curr_velocity[3];
	double curr_mass;

	double curr_velocity_component;

	for (int i = 0; i < num_particles; i++) {
		curr_mass = masses[i];
		for (int j = 0; j < 3; j++) {
			//generate velocity from Maxwell-Boltzmann distrib
			curr_velocity_component = (rng->gaussian())*sqrt(kT / curr_mass);
			//assign array element
			curr_velocity[j] = curr_velocity_component;

			//add component of momentum to corresponding component of total velocity
			total_momentum_first_assigned[j] += curr_mass*curr_velocity_component;
		}

		velocities[i] = posn(curr_velocity[0], curr_velocity[1], curr_velocity[2]);
	}

	//make sure total velocity of system = 0
	//go through particles and subtract original avg momentum of system from each one

	//get (orig) avg momentum
	double avg_momentum_first_assigned[3] = { total_momentum_first_assigned[0] / num_particles, total_momentum_first_assigned[1] / num_particles, total_momentum_first_assigned[2] / num_particles };

	//go thru and subtract it
	for (int i = 0; i < num_particles; i++) {
		curr_mass = masses[i];
		velocities[i].x -= avg_momentum_first_assigned[0] / curr_mass;
		velocities[i].y -= avg_momentum_first_assigned[1] / curr_mass;
		velocities[i].z -= avg_momentum_first_assigned[2] / curr_mass;
	}




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


//writing to files


string default_filepath = "outputs_absposn_fixed/";

double default_double_arr[1] = { 0. };

void print_double_array_to_file(string filename, double* input_array, int num_entries, string filehead = "", bool has_header = false, string filetype = "csv", bool prefix_2nd_array = false, double*second_array = default_double_arr) {
	string separator = "\t";
	if (filetype == "csv") { separator = ", "; }

	string full_filename = default_filepath + (filename + ("." + filetype));
	ofstream myfile;
	myfile.open(full_filename, ios::app);
	//if has header, include header line (blank by default)
	if (has_header) {
		myfile << filehead + "\n";
	}
	//loop through particles and write their coordinates, separated by desired separator
	//
	for (int i = 0; i < num_entries; i++) {
		if (prefix_2nd_array) {
			myfile << second_array[i];
			myfile << separator;
		}
		myfile << input_array[i];
		myfile << "\n";
	}
	myfile.close();


}


void print_lattice_to_file(string filename, posn* positions, int num_particles, string filehead = "index, x, y, z", bool has_header = false, string filetype = "csv", string atomtype_if_xyz = "Ar") {
	string separator = "\t";
	if (filetype == "csv") { separator = ", "; }

	string full_filename = default_filepath + (filename + ("." + filetype));
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
		if (filetype == "xyz") {
			myfile << atomtype_if_xyz;
			myfile << separator;
		}
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
	while (pbc_x < 0.0) { pbc_x += 1.0*boxlen; }

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



//Lennard-Jones force, pressure, and energy calculation functions
// if has no cutoff, r_cutoff is not used, so it is given a default value for this use case
// energy does not include U_tail = (8./3) pi rho*((1/3)r_cutoff^(-9) - r_cutoff^(-3))
//	because this requires additional input of rho, is inefficient for our fixed-rho sims, and is irrelevant for our Monte Carlo steps.

//note: by default uses LJ coordinates with epsilon, sigma = 1

//LJ force calculation - seems to yield incorrect results when used for pressure calculation
//force on i from j = 4*epsilon * (12 sigma^12 r_ij^(-13) - 6 sigma^6 r_ij^(-7)) r_ji hat
//where r_ji hat = unit vector in (r_i - r_j direction)
//can be given (r/sigma)^(-6), (r/sigma)
posn LJ_force_j_on_i(posn part_i, posn part_j, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2, bool scaled_distance_inv_and_pow_given = false, double scaled_distance_inv = 1.0, double scaled_distance_invsixthpow = 1.0) {
	posn force_vec = posn(0., 0., 0.);

	//get vector between positions 
	posn i_minus_j_minim = minimum_image_displacement_posns(part_j, part_i, boxdims);

	if (!scaled_distance_inv_and_pow_given) {
		scaled_distance_inv = sigma / i_minus_j_minim.length();
		scaled_distance_invsixthpow = pow(scaled_distance_inv, 6.0);
	}


	//if distance less than cutoff radius, calculate force
	if ((1. / scaled_distance_inv) < scaled_r_cutoff) {

		//magnitude of force = 4 * epsilon * [term]
		//[term] = 12 sigma^12/r^13 - 6 sigma^6 / r^7
		// = 6 sigma^6/r^7 (2 sigma^6/r^6 - 1)
		// = 6 (r*^(-6)) * (1/r*)/sigma * (2r*^(-6) - 1)
		double F_magnitude = 4. * epsilon;
		F_magnitude *= 6 * scaled_distance_invsixthpow * (scaled_distance_inv / sigma) * (2 * scaled_distance_invsixthpow - 1.0);

		//set force vector equal to unit vector in r_i - r_j direction
		force_vec = i_minus_j_minim;
		force_vec.normalize();

		//scale force vector by magnitude of force
		//note that sign of magnitude is correct because it is positive (pushing i away from j) when distance is small, and vice-versa.
		force_vec.scale(F_magnitude);


	}
	return force_vec;
}

posn LJ_force_b_on_a_redone(posn part_a, posn part_b, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2) {
	posn f_ba = posn(0., 0., 0.);
	
	double disp[3];
	disp[0] = part_a.x - part_b.x;
	disp[1] = part_a.y - part_b.y;
	disp[2] = part_a.z - part_b.z;


	//apply pbcs: disp is a double array that contains x, y, z displacement
	for (int dim = 0; dim < 3; dim++) {
		if (disp[dim] > .5*boxdims[dim]) { disp[dim] -= boxdims[dim]; }
		else if (disp[dim] < -.5*boxdims[dim]) { disp[dim] += boxdims[dim]; }
	}

	//get magnitude of displacement, r_ij = (xdisp^2 + ydisp^2 + zdisp^2)^(1/2)
	double r_ij = sqrt(disp[0] * disp[0] + disp[1] * disp[1] + disp[2] * disp[2]);

	//if r_ij greater than cutoff, return default (0) force
	if (use_cutoff & (r_ij > scaled_r_cutoff)) {
		return f_ba;
	}
	
	//otherwise, calculate force
	else {
		//have displacement between particles i and j (a and b)
		//now get its inverse and its inverse 6th power for use in force magnitude calculation

		double r_ij_inv = 1. / r_ij;
		double r_ij_inv6th = pow(r_ij_inv, 6.);


		//magnitude of F = 4*epsilon * 6 x^6 * x/sigma (2*x^6 - 1), where x = (sigma/r_ij)   
		// since F is in ij direction, F = (/r_i/ - /r_j/) direction unit-vector * 24*epsilon * (x^7 /sigma) * (2*x^6 - 1)
		//								= (/r_i/ - /r_j/) * [24 * epsilon * (x^8/sigma) * (2*x^6 - 1)]
		//													[magnitude of F divided by r_ij]
		double F_magnitude_div_rij = 24. * epsilon;
		F_magnitude_div_rij *= (r_ij_inv6th*r_ij_inv*r_ij_inv / sigma) * (2.*r_ij_inv6th - 1.);


		//scale displacement array ( = /r_i/ - /r_j/) by |F|/r_ij to get |F| (/r_ij/ direction) = /F/
		f_ba = posn(disp[0] * F_magnitude_div_rij, disp[1] * F_magnitude_div_rij, disp[2] * F_magnitude_div_rij);

		return f_ba;
	}
}

//fill an array of num_particles particles with the forces on the particles
void fill_LJ_forces(posn* force_array, posn* positions, int num_particles, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2) {
	posn i_posn;
	posn j_posn;

	posn force_on_i_from_lesser_js;
	posn f_j_on_i;

	//clear previous forces
	for (int i = 0; i < num_particles; i++) {
		force_array[i].set_all(0., 0., 0.);
	}

	for (int i = 1; i < num_particles; i++) {
		i_posn = positions[i];
		force_on_i_from_lesser_js.set_all(0., 0., 0.);
		for (int j = 0; j < i; j++) {
			j_posn = positions[j];
			//add F_{j on i} to F_i
			f_j_on_i = LJ_force_j_on_i(i_posn, j_posn, epsilon, sigma, boxdims, true, scaled_r_cutoff);

			//subtract F_j_on_i (= add F_i_on_j) to F_j
			force_array[j] = add_posns(force_array[j], f_j_on_i, 1.0, -1.0);

			//and add it to F_i
			force_on_i_from_lesser_js = add_posns(force_on_i_from_lesser_js, f_j_on_i);
		}

		//write force to array
		force_array[i] = add_posns(force_array[i], force_on_i_from_lesser_js);
	}

}

void fill_LJ_forces_new(posn* force_array, posn* positions, int num_particles, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2) {
	posn i_posn;
	posn j_posn;

	posn force_on_i_from_lesser_js;
	posn f_j_on_i;

	//clear previous forces
	for (int i = 0; i < num_particles; i++) {
		force_array[i].set_all(0., 0., 0.);
	}

	for (int i = 1; i < num_particles; i++) {
		i_posn = positions[i];
		force_on_i_from_lesser_js.set_all(0., 0., 0.);
		for (int j = 0; j < i; j++) {
			j_posn = positions[j];
			//add F_{j on i} to F_i
			f_j_on_i = LJ_force_b_on_a_redone(i_posn, j_posn, epsilon, sigma, boxdims, true, scaled_r_cutoff);

			//subtract F_j_on_i (= add F_i_on_j) to F_j
			force_array[j] = add_posns(force_array[j], f_j_on_i, 1.0, -1.0);

			//and add it to F_i
			force_on_i_from_lesser_js = add_posns(force_on_i_from_lesser_js, f_j_on_i);
		}

		//write force to array
		force_array[i] = add_posns(force_array[i], force_on_i_from_lesser_js);
	}

}

//copy and paste of above: just changed inputs so that almost-always-default inputs epsilon and sigma are at the end of the list
void fill_LJ_forces_redone_changedinputs(posn* force_array, posn* positions, int num_particles, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2, double epsilon = 1.0, double sigma = 1.0) {
	posn i_posn;
	posn j_posn;

	posn force_on_i_from_lesser_js;
	posn f_j_on_i;

	//clear previous forces before starting force-calculation loop
	for (int i = 0; i < num_particles; i++) {
		force_array[i].set_all(0., 0., 0.);
	}

	for (int i = 1; i < num_particles; i++) {
		i_posn = positions[i];
		force_on_i_from_lesser_js.set_all(0., 0., 0.);
		for (int j = 0; j < i; j++) {
			j_posn = positions[j];
			//add F_{j on i} to F_i
			f_j_on_i = LJ_force_b_on_a_redone(i_posn, j_posn, epsilon, sigma, boxdims, true, scaled_r_cutoff);

			//subtract F_j_on_i (= add F_i_on_j) to F_j
			force_array[j] = add_posns(force_array[j], f_j_on_i, 1.0, -1.0);

			//and add it to F_i
			force_on_i_from_lesser_js = add_posns(force_on_i_from_lesser_js, f_j_on_i);
		}

		//write force to array
		force_array[i] = add_posns(force_array[i], force_on_i_from_lesser_js);
	}

}


//VIRIAL PRESSURE calculations

//calculate tail term of virial pressure 
// = (16/3) pi rho^2 (v^3 - v)
//where rho is density, u is 1/r_cutoff, v = u^3
//note we will calculate the scaled pressure, using scaled density and implicitly scaled coordinates
//can later add option to get non-scaled pressure in this function

//my math library doesn't seem to have pi as a constant, so I'm putting it in manually
double virial_pressure_tail(double scaled_density, double scaled_r_cutoff) {
	double rcut_inv_3 = pow(scaled_r_cutoff, -3.);
	double scaled_tail_term = (16. / 3.) * (3.1415926536) * scaled_density*scaled_density*rcut_inv_3*((2. / 3.)*rcut_inv_3*rcut_inv_3 - 1.);

	return scaled_tail_term;
}

//calculate virial pressure force term = 1/3V * sum(F_i dot r_i) from force array
double virial_pressure_from_forces(posn* positions, posn* forces, int num_particles, double volume = 1.0) {
	double virial_pressure_value = 0.0;
	for (int i = 0; i < num_particles; i++) {
		virial_pressure_value += dot_posns(positions[i], forces[i]);
	}

	virial_pressure_value *= 1 / (3. * volume);
	return virial_pressure_value;
}

//virial pressure calculation, no forces needed
//however, requires posn[num_particles] array input, which will be filled with forces
double virial_pressure(posn* positions, int num_particles, posn* force_holder, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2, double scaled_density = 1.0) {

	fill_LJ_forces(force_holder, positions, num_particles, epsilon, sigma, boxdims, use_cutoff, scaled_r_cutoff);
	double volume = boxdims[0] * boxdims[1] * boxdims[2];

	double virial_pressure = virial_pressure_from_forces(positions, force_holder, num_particles, volume);

	//if use cutoff, add tail term = (16/3)pi rho^2((2/3)u^9 - u^3) = (16/3) pi rho^2 v ((2/3)v^2 - 1)
	//note we will calculate the scaled pressure, using scaled density and implicitly scaled coordinates
	//can later add option to get non-scaled pressure in this function
	//note that epsilon and sigma inputs = epsilon and sigma as used in force calculation, so they will be scaled if the sim is scaled
	//and hence cannot be used to get unscaled pressure
	if (use_cutoff) {
		double scaled_tail_term = virial_pressure_tail(scaled_density, scaled_r_cutoff);
		virial_pressure += scaled_tail_term;
	}

	return virial_pressure;
}

//new virial pressure method - old ones were giving me strange overall values, so I'm going to write a full method for it that hopefully will work
double virial_pressure_redone(posn* positions, int num_particles, posn* force_holder, double epsilon = 1.0, double sigma = 1.0, double* boxdims = default_boxdims, bool use_cutoff = true, double scaled_r_cutoff = .2, double scaled_density = 1.0) {
	double P_virial_no_tail = 0.0;

	//P_virial_no_tail = tail term + (1/(3V)) * sum(F_i dot r_i)
	//this sum = sum_{all i != j, including repeats)(F_j on i dot r_i)
	// = sum(all pairs i != j, no repeats)(F_j on i dot (r_i - r_j)
	posn i_posn;
	posn j_posn;

	//displacement vector - not using as a posn for now
	double disp[3];


	double r_ij;
	double r_ij_inv6th;
	double P_contribution;

	for (int i = 1; i < num_particles; i++) {
		i_posn = positions[i];

		for (int j = i - 1; j < i; j++) {
			j_posn = positions[j];
			disp[0] = i_posn.x - j_posn.x;
			disp[1] = i_posn.y - j_posn.y;
			disp[2] = i_posn.z - j_posn.z;


			//apply pbcs: disp is a double array that contains x, y, z displacement
			for (int dim = 0; dim < 3; dim++) {
				if (disp[dim] > .5*boxdims[dim]) { disp[dim] -= boxdims[dim]; }
				else if (disp[dim] < -.5*boxdims[dim]) { disp[dim] += boxdims[dim]; }
			}

			//get magnitude of displacement, r_ij = (xdisp^2 + ydisp^2 + zdisp^2)^(1/2)
			r_ij = sqrt(disp[0] * disp[0] + disp[1] * disp[1] + disp[2] * disp[2]);

			//if r_ij greater than cutoff, skip:
			if (use_cutoff & (r_ij > scaled_r_cutoff)) { continue; }


			//and get its inverse 6th power, for use in pressure contribution calculation
			r_ij_inv6th = pow(r_ij, -6.);


			//magnitude of F = 4*epsilon * 6 x^6 * x/sigma (2*x^6 - 1), where x = (sigma/r_ij)   
			//(below uses /v/ to indicate vectors, v for the same vector's magnitude)
			//this magnitude * /F's direction/, dotted into /r_ij/, is the ij contribution to the pressure
			//note /F/ is in /r_i/ - /r_j/ = /r_ij/ direction
			//so the ij scalar contribution to the pressure is simply the scalar r_ij times the magnitude of F
			// = 4*epsilon * (x^6 /sigma) * (2*x^6 - 1)
			P_contribution = 24. * epsilon;
			P_contribution *= (r_ij_inv6th / sigma) * (2.*r_ij_inv6th - 1.);

			P_virial_no_tail += P_contribution;

		}
	}


	//force-based virial pressure term must be scaled by 1/(3V)
	double volume = boxdims[0] * boxdims[1] * boxdims[2];
	P_virial_no_tail *= 1. / (3. * volume);

	return P_virial_no_tail;
}


posn default_position_array[1] = { posn(0,0,0) };

//LJ ENERGY calculations

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
	posn f_j_on_i;


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
				//add F_ji to force on i, subtract it from force on j
				f_j_on_i = LJ_force_j_on_i(i_posn, j_posn, epsilon, sigma, boxdims, has_cutoff, r_cutoff, true, scaled_dist_inv, scaled_dist_inv6th);
				net_force = add_posns(net_force, f_j_on_i);
				force_array[j] = add_posns(force_array[j], f_j_on_i, 1.0, -1.0);
			}
		}
	}

	//if calculating force, updated desired element of force array
	if (calc_force_too) { force_array[i] = add_posns(force_array[i], net_force); }

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

//Monte Carlo implementation

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

void MC_cycle_hardsph(posn* positions, int num_particles, RanMars* rng, double sph_radius_d0, double scalefac = 1.0) {
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
	double box_dimensions[3] = { lattice_param, lattice_param, lattice_param };

	//fill provided positions array with particles initialized on a simple cubic lattice fitting the input params
	//note that lattice must be initialized before this function is called
	fill_simple_symm_cubic_lattice(positions, lattice_side_numparts, num_particles, lattice_param);

	//for now, will scale proposal moves by .1 * volume^(1/3) factor
	double scalefac = .1*lattice_param;

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
void MC_cycle_LJ(posn* positions, int num_particles, RanMars *rng, double* box_dims, double temperature_scaled, bool use_LJ_cutoff, double epsilon, double sigma, double r_cutoff = .2, double rand_move_scaler = 1.0) {
	for (int i = num_particles - 1; i > 0; i--) {
		//generate x and y move RVs; RanMars.uniform gives unif(0,1) and we want unif(-1,1), so we multiply by 2 and subtract 1.
		//for our hard-sphere model, we relied on scaling the move in the individual MC steps
		//here, we will do scaling here, by a fixed ratio .5*relevant box dim, since the ranges from [.5,1] and [-1,-.5] in box-side units are made redundant by PBCs

		double x_move = rand_move_scaler * .5*box_dims[0] * ((2 * rng->uniform()) - 1);
		double y_move = rand_move_scaler * .5*box_dims[1] * ((2 * rng->uniform()) - 1);
		double z_move = rand_move_scaler * .5*box_dims[2] * ((2 * rng->uniform()) - 1);

		posn proposed_move = posn(x_move, y_move, z_move);


		//choose particle to move
		//generate uniform double between 0 and num_particles; 
		//for each particle i from 0 to num_particles-1 there is probability 1/num_particles of being between i and i + 1
		//so casting the result to an int gives us a uniform integer between 0 and num_particles

		int part_index = (int)(num_particles * rng->uniform()) % num_particles;

		//I think it might be possible for corner-cases involving floating point arithmetic to yield num_particles here
		//so this while loop makes sure it is in the proper region

		while ((part_index < 0) | (part_index >= num_particles)) {
			part_index = (int)(num_particles * rng->uniform()) % num_particles;
		}

		//finally, generate acceptance threshold uniformly between 0,1
		double accept_threshold = rng->uniform();


		MC_update_LJ(positions, part_index, proposed_move, accept_threshold, box_dims, temperature_scaled, num_particles, use_LJ_cutoff, epsilon, sigma, r_cutoff);
	}
}


//Molecular Dynamics implementation
intposn default_absolute_disp = intposn(0, 0, 0);
intposn default_absolute_disp_array[1] = { default_absolute_disp };


//update functions

//apply basic Verlet algorithm for a single particle + force
posn default_posn = posn(0.,0.,0.);
int default_boxes_away[3] = { 0,0,0 };

posn posn_post_Verlet_given_velocity(posn curr_posn, posn* velocity_pointer, posn force, double* box_dims, double mass = 1.0, double dt = .005, bool update_velocity = true, bool update_absolute_posn = false, intposn* absolute_posn = default_absolute_disp_array) {
	//Verlet algorithm: x(t+dt) = x(t) + (x(t) - x(t - dt)) + (1/m) F(t)dt^2
	//note (x(t) - x(t - dt))/dt is serving as velocity here; we can as well write x(t+dt) = x(t) + v(t)dt + (1/m) F(t)dt^2
	//default dt = .005 will be the usual delta-time we use for lennard-jones particles in scaled coordinates
	//that is, this is meant to be used as dt = .005 in scaled-time units (which will be the case if posns + force + mass are scaled)

	//first let's take 2x(t) - x(t-dt) = x(t) + v(t)dt
	//posn new_posn = add_posns(curr_posn,velocity,1,dt);
	posn velocity = velocity_pointer[0];

	double curr_x = curr_posn.x;
	double curr_y = curr_posn.y;
	double curr_z = curr_posn.z;

	double new_x = curr_posn.x + velocity.x * dt;
	double new_y = curr_posn.y + velocity.y * dt;
	double new_z = curr_posn.z + velocity.z * dt;

	//then add F * (dt^2)/m
	//new_posn = add_posns(new_posn, force, 1.0, dt*dt / mass);

	new_x += force.x * (dt*dt) / mass;
	new_y += force.y * (dt*dt) / mass;
	new_z += force.z * (dt*dt) / mass;

	posn new_posn = posn(new_x, new_y, new_z);

	//apply periodic boundary conditions: keep in box of dimensions = box_dims

	//create array to hold position elements- easier to loop through than posn
	double new_posn_array[3] = { new_posn.x, new_posn.y, new_posn.z };
	//update velocity (if input says so): note that velocity = direction of move ignoring application of PBCs
	posn v_out = posn((new_x - curr_x) / dt, (new_y - curr_y) / dt, (new_z - curr_z) / dt) ;

	velocity_pointer[0] = v_out;


	double curr_boxlen;
	int disp_from_center_box[3] = { 0,0,0 };

	//loop through array and apply periodic boundary conditions to each component: keep within [0, L_i)
	for (int i = 0; i < 3; i++){
		curr_boxlen = box_dims[i];

		if (new_posn_array[i] >= curr_boxlen) { 
			new_posn_array[i] -= curr_boxlen;
			//if keeping track of absolute position (for autocorrelation fn), update that as well
			if (update_absolute_posn) {
				disp_from_center_box[i] += 1;
			}
		}
		else if (new_posn_array[i] < 0.0) {
			new_posn_array[i] += curr_boxlen;
			//if keeping track of absolute position (for autocorrelation fn), update that as well
			if (update_absolute_posn) {
				disp_from_center_box[i] -= 1;
			}
		}
	}


	//pull back updates from array into the posn
	new_posn = posn(new_posn_array[0], new_posn_array[1], new_posn_array[2]);

	//and into absolute posn
	if (update_absolute_posn) {
		absolute_posn->x += disp_from_center_box[0];
		absolute_posn->y += disp_from_center_box[1];
		absolute_posn->z += disp_from_center_box[2];
	}


	return new_posn;
}


posn posn_post_Verlet_given_old_posn(posn curr_posn, posn old_posn, posn force, double* box_dims, double mass = 1.0, double dt = .005, bool update_absolute_posn = false, intposn* absolute_posn= default_absolute_disp_array) {

	//same as above, but have to get v from old_posn, curr_posn:
	//v = (x(t) - x(t-dt))/dt
	posn velocity = add_posns(curr_posn, old_posn, 1. / dt, -1. / dt);
	posn vel_pointer[1] = { velocity };

	return posn_post_Verlet_given_velocity(curr_posn, vel_pointer, force, box_dims, mass, dt, false, update_absolute_posn, absolute_posn);
} 



//verlet update given forces: updates positions, velocities, and optionally absolute displacements
void update_Verlet_given_forces(posn* positions, posn* velocities, int num_particles, double*boxdims, posn* forces, double* masses, bool update_absolute_posn = false, intposn* absolute_disps = default_absolute_disp_array, double epsilon = 1.0, double sigma = 1.0, double dt = .005) {
	posn old_posn;
	posn velocity_pointer[1];

	intposn curr_abs_posn[1];
	for (int i = 0; i < num_particles; i++) {
		old_posn = positions[i];

		velocity_pointer[0] = velocities[i];

		//get new position
		//if updating absolute posn, create pointer, pass it, get back value; otherwise don't:
		if(update_absolute_posn){
			curr_abs_posn[0] = absolute_disps[i];
			positions[i] = posn_post_Verlet_given_velocity(positions[i], velocity_pointer, forces[i], boxdims, masses[i], dt, true, true, curr_abs_posn);
			absolute_disps[i].set_all(curr_abs_posn[0].x, curr_abs_posn[0].y, curr_abs_posn[0].z);
		}
		else{ 
			positions[i] = posn_post_Verlet_given_velocity(positions[i], velocity_pointer, forces[i], boxdims, masses[i], dt, true, false);
		}
		//pull pointer's new pointed-to value back into velocity array 
		velocities[i] = velocity_pointer[0];
	}
}

//full verlet update function: calcs forces and updates
void update_Verlet(posn* positions, posn* velocities, int num_particles, double*boxdims, posn* forces, double* masses, bool update_absolute_posn = false, intposn* absolute_disps = default_absolute_disp_array, bool use_r_cutoff = true, double scaled_r_cutoff = 3., double epsilon = 1.0, double sigma = 1.0, double dt = .005) {
	//update forces
	fill_LJ_forces_redone_changedinputs(forces, positions, num_particles, boxdims, use_r_cutoff, scaled_r_cutoff);
	//update positions, velocities, and (optionally) absolute positions
	update_Verlet_given_forces(positions, velocities, num_particles, boxdims, forces, masses, update_absolute_posn, absolute_disps, epsilon, sigma, dt);
}


//calculate kinetic energy (placed here because only relevant for molecular dynamics)

double get_KE(posn* velocities, double* masses, int num_particles) {
	double KE = 0.;
	for (int i = 0; i < num_particles; i++) {
		//KE += .5 * mv^2
		KE += .5* masses[i] * (velocities[i].x*velocities[i].x + velocities[i].y*velocities[i].y + velocities[i].z*velocities[i].z);
	}
	return KE;
}


//get total distance^2 from current posn, init posn, and absolute posn (of box particle is currently in) + box dimensions

double get_delta_r_sq(posn curr_posn, posn init_posn, intposn absolute_posn, double* box_dims) {
	posn displacement = curr_posn;
	
	//get (box-crossing vector) multiplied by corresponding box-dimensions elements
	displacement.x += absolute_posn.x * box_dims[0] - init_posn.x;
	displacement.y += absolute_posn.y * box_dims[1] - init_posn.y;
	displacement.z += absolute_posn.z * box_dims[2] - init_posn.z;


	double delta_r_sq = displacement.x*displacement.x + displacement.y*displacement.y + displacement.z*displacement.z;

	return delta_r_sq;
}

double get_avg_delta_r_sq(posn* positions, int num_particles, posn* init_positions, intposn* absolute_positions, double* box_dims) {
	double total_delta_r_sq = 0.;

	for (int i = 0; i < num_particles; i++) {
		total_delta_r_sq += get_delta_r_sq(positions[i], init_positions[i], absolute_positions[i], box_dims);
	}

	double avg_delta_r_sq = total_delta_r_sq / num_particles;
	return avg_delta_r_sq;
}



//Parts 3 on: running loops

double default_input_double_array[1] = { 0.0 };

// Monte Carlo loop: initialize lattice, run a few equilibration cycles, then run production cycles
void MC_loop_LJ(int num_prod_cycles, posn* positions, int num_particles, double density_scaled, double temperature_scaled, double sigma, double epsilon, bool scale_to_LJ_params = true, bool use_LJ_cutoff = true, double r_cutoff_LJ_scaled = .2, bool print_all_positions_to_file = true, string posn_filename = "irrelevant", bool calc_energies = true, double* energies = default_input_double_array, bool also_calc_force_and_pressure = false, posn* force_array = default_position_array, double* virial_pressures = default_input_double_array, double rand_move_scalefac = 1.0) {
	int	num_equil_cycles = 5;

	//initialize RNG
	RanMars *random;
	//seed with time
	unsigned int seed = unsigned int(time(NULL));
	seed = abs((int)seed) % 900000000;
	random = new RanMars(seed);


	//initialize lattice

	//get n = # of particles per lattice side; want cubic lattice s.t. n = N^(1/3) (or smallest integer above this)
	int lattice_side_numparts = ceil(pow((double)num_particles, 1. / 3.));

	//get lattice sidelength L: density = N/(L^3), so L = (N/density)^(1/3)
	//density scaled = density*sigma^3, so we'll get back regular density first
	// (there's probably a less clunky way to do this, but I didn't want to mess around with factors of sigma
	// and wanted to leave open the ability to use non-scaled coordinates)
	double density = density_scaled / (pow(sigma, 3.));
	double box_sidelen = pow(num_particles / density, 1. / 3.);

	//if using scaled coordinates x* = x/sigma, etc., then scale lattice param by sigma
	//and pass sigma* = 1.0, epsilon* = 1.0 to cycling functions
	double sigma_to_pass = sigma;
	double epsilon_to_pass = epsilon;
	if (scale_to_LJ_params) {
		box_sidelen *= 1. / (sigma);
		sigma_to_pass = 1.0;
		epsilon_to_pass = 1.0;
	}
	//also pass appropriately scaled cutoff radius
	double r_cutoff_to_pass = r_cutoff_LJ_scaled * sigma;
	if (scale_to_LJ_params) { r_cutoff_to_pass = r_cutoff_LJ_scaled; }

	//set box-dimension array = [sidelen, sidelen, sidelen]
	//gives side lengths of box in units we are using--either scaled or unscaled ones
	double box_dimensions[3] = { box_sidelen, box_sidelen, box_sidelen };

	//fill provided positions array with particles initialized on a simple cubic lattice fitting the input params
	//note that lattice must be initialized before this function is called
	fill_simple_symm_cubic_lattice(positions, lattice_side_numparts, num_particles, box_sidelen);


	//equilibration cycles 
	for (int i = 0; i < num_equil_cycles; i++) {
		MC_cycle_LJ(positions, num_particles, random, box_dimensions, temperature_scaled, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, r_cutoff_to_pass, rand_move_scalefac);
		if (print_all_positions_to_file) {
			print_lattice_to_file(posn_filename, positions, num_particles, "equilibration cycle", true, "xyz", "Ar");
		}
	}

	//production cycles

	//initialize variables for force+pressure calc if necessary
	//cheap, so would rather calculate them always than accidentally not in some cases, which could happen if I add significantly to this code later	
	double scaled_volume = box_dimensions[0] * box_dimensions[1] * box_dimensions[2];
	double virial_pressure = 0.;


	//note - changed to use new force-pressure calc; doesn't involve calculating forces during energy calc step
	for (int j = 0; j < num_prod_cycles; j++) {
		MC_cycle_LJ(positions, num_particles, random, box_dimensions, temperature_scaled, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, r_cutoff_to_pass, rand_move_scalefac);

		if (print_all_positions_to_file) {
			print_lattice_to_file(posn_filename, positions, num_particles, "production cycle", true, "xyz", "Ar");
		}


		//if calculating energies, get current energy and put it in the energies array argument
		if (calc_energies) {
			energies[j] = LJ_energy_w_cutoff(positions, num_particles, use_LJ_cutoff, epsilon_to_pass, sigma_to_pass, box_dimensions, r_cutoff_to_pass, also_calc_force_and_pressure, force_array);
			if (also_calc_force_and_pressure) {
				virial_pressure = virial_pressure_redone(positions, num_particles, force_array, epsilon, sigma, box_dimensions, use_LJ_cutoff, r_cutoff_to_pass, scaled_volume);
				virial_pressure += virial_pressure_tail(density_scaled, r_cutoff_to_pass);
				virial_pressures[j] = virial_pressure;

			}
		}
	}
}



// 
void main() {
	constexpr int num_prod_cycles_test = 300;
	constexpr int num_prod_cycles = 5000;

	constexpr int num_particles = 343;

	//initialize RNG
	RanMars *random;
	//seed with time
	unsigned int seed = unsigned int(time(NULL));
	seed = abs((int)seed) % 900000000;
	random = new RanMars(seed);



	//initialize arrays
	posn positions[num_particles];
	posn velocities[num_particles];
	posn forces[num_particles];

	posn init_positions[num_particles];
	intposn absolute_displacement[num_particles];

	
	//normalization decision - choose epsilon, sigma 'already' 1 in our units -> input r_cutoff, rho, T must be scaled accordingly
	double epsilon = 1.;
	double sigma = 1.;
	//COMMENTED OUT - test MD by checking conservation of total energy
	//test: get energies over time for rho* = .5, T* = 2.


	//initialize arrays
	posn positions_Etest[num_particles];
	posn velocities_Etest[num_particles];
	posn forces_Etest[num_particles];

	//posn init_absolute_displacement[num_particles];
	//posn absolute_displacement[num_particles];

	double KEs_Etest[num_prod_cycles_test];
	double PEs_Etest[num_prod_cycles_test];


	//input parameters
	double r_cutoff = 3.;
	double rho_star = .5;
	double T_star = 2.;
	double dt = .005;
	double mass = 1.;

	//initialize mass array from input mass - kind of clunky to use array in this case, but not terrible I think
	double masses[num_particles];
	for (int imass = 0; imass < num_particles; imass++) {
		masses[imass] = mass;
	}


	//derive lattice parameters, initialize lattice
	double V = num_particles/ rho_star;
	double sidelen = pow(V, 1 / 3.);
	double boxdims[3] = { sidelen, sidelen, sidelen };

	double approx_parts_per_side = pow(num_particles, 1 / 3.);
	int num_parts_per_side = ceil(approx_parts_per_side);
	
	//use ceil so that there are enough lattice positions if num_particles is not a cube

	//if num particles is perfect cube, make sure rounding error doesn't add extra particle per side:
	int alt_num_parts_per_side = floor(approx_parts_per_side);
	if (alt_num_parts_per_side*alt_num_parts_per_side*alt_num_parts_per_side == num_particles) {
		num_parts_per_side = alt_num_parts_per_side;
	}

	fill_simple_symm_cubic_lattice(positions_Etest, num_parts_per_side, num_particles, sidelen);


	//initialize velocities
	fill_Boltzmann_velocities(velocities_Etest, num_particles, random, masses, T_star * epsilon);

	//loop until desired number of steps (incl initialization step) have been completed
	for (int i = 0; i < num_prod_cycles_test; i++) {
		//update energy arrays
		PEs_Etest[i] = LJ_energy_w_cutoff(positions_Etest, num_particles, true, epsilon, sigma, boxdims, r_cutoff, false);
		KEs_Etest[i] = get_KE(velocities_Etest, masses, num_particles);

		//write lattice to file
		print_lattice_to_file("test_positions", positions_Etest, num_particles, "x, y, z", true, "xyz");

		//time-evolve system if not at last half-step
		if (i == num_prod_cycles_test - 1) { continue; }
		else {
			update_Verlet(positions_Etest, velocities_Etest, num_particles, boxdims, forces_Etest, masses, false, default_absolute_disp_array, true, 3.0, epsilon, sigma, dt);
		}
	}

	
	//gather KE/N, PE/N, E/N into single big array
	//kind of kludgy but I'm using an array of posns (structs w 3 doubles each) for this - easier to write to file w current methods
	posn energies[num_prod_cycles_test];
	for (int i = 0; i < num_prod_cycles_test; i++) {
		energies[i].x = KEs_Etest[i] / (double) num_particles;
		energies[i].y = PEs_Etest[i] / (double) num_particles;
		energies[i].z = energies[i].x + energies[i].y;
	}

	
	//write KE, PE, E to file

	string energy_filename = "energies_from_LJ_MD_test_fixed_" + to_string(num_particles) +"_particles";
	string E_header = "KE, PE, E";
	
	
	print_lattice_to_file(energy_filename, energies, num_prod_cycles_test, E_header, true, "tsv");



	//calculate diffusion constant via <(delta r^2)(t)>

	//input parameters
	//COMMENTED OUT - msd calculation with rho variation 
	constexpr int num_rhos = 8;
	double rho_stars[num_rhos] = {.1, .2, .3, .4, .5, .7, .9, 1.1};
	double r_cutoff = 3.;
	//note: temporarily upped T (and lowered t in compensation) to test vs literature
	double T_star = 2.5;
	double dt = .005;
	double mass = 1.;

	//initialize mass array from input mass - kind of clunky to use array in this case, but not terrible I think
	double masses[num_particles];
	for (int imass = 0; imass < num_particles; imass++) {
		masses[imass] = mass;
	}

	double avg_delta_rsq_t[num_prod_cycles];

	//loop thru densities
	for (int i_rho = 0; i_rho < num_rhos; i_rho++) {
		//clear absolute positions
		for (int i_disp = 0; i_disp < num_particles; i_disp++) {
			absolute_displacement[i_disp] = intposn(0, 0, 0);
		}

		double rho_star = rho_stars[i_rho];
		
		//derive lattice parameters, initialize lattice
		double V = num_particles / rho_star;
		double sidelen = pow(V, 1 / 3.);
		double boxdims[3] = { sidelen, sidelen, sidelen };

		double approx_parts_per_side = pow(num_particles, 1 / 3.);
		int num_parts_per_side = ceil(approx_parts_per_side);

		//use ceil so that there are enough lattice positions if num_particles is not a cube

		//if num particles is perfect cube, make sure rounding error doesn't add extra particle per side:
		int alt_num_parts_per_side = floor(approx_parts_per_side);
		if (alt_num_parts_per_side*alt_num_parts_per_side*alt_num_parts_per_side == num_particles) {
			num_parts_per_side = alt_num_parts_per_side;
		}

		//initialize lattice, and do same for lattice of init posns
		//(since lattice initialization is (currently) deterministic, this shouldn't present any problems.)
		fill_simple_symm_cubic_lattice(positions, num_parts_per_side, num_particles, sidelen);
		fill_simple_symm_cubic_lattice(init_positions, num_parts_per_side, num_particles, sidelen);


		//initialize velocities
		fill_Boltzmann_velocities(velocities, num_particles, random, masses, T_star * epsilon);

		//loop until desired number of steps (incl initialization step) have been completed
		for (int i = 1; i < num_prod_cycles; i++) {
			//write lattice to file
			print_lattice_to_file("positions_"+ to_string(int(rho_star*100.)), positions, num_particles, "x, y, z", true, "xyz");

			//update <(delta_r)^2>
			avg_delta_rsq_t[i] = get_avg_delta_r_sq(positions, num_particles, init_positions, absolute_displacement, boxdims);

			//write delta_r_sq to file
			double baby_array[1] = { avg_delta_rsq_t[i] };
			double baby_time_array[1] = { (i-1)*dt };
			print_double_array_to_file("r_sq_t_for_rho_pt_" + to_string(int(rho_star*100.)), baby_array, 1, "", false, "dat", true, baby_time_array);


			//time-evolve system if not at last half-step
			if (i == num_prod_cycles - 1) { continue; }
			else {
				update_Verlet(positions, velocities, num_particles, boxdims, forces, masses, true, absolute_displacement, true, 3.0, epsilon, sigma, dt);
			}

			////at desired steps, resample velocities
			//if (i % 100 == 0 & i > 0) {
			//	fill_Boltzmann_velocities(velocities, num_particles, random, masses, T_star*epsilon);

			//}			
		}
		
	}









	//-----------------------------------------------------------------------
	//commented out Monte Carlo stuff - might be useful source for loop-over-data methods
	//still working on getting lattice initialization done in a way that allows for user input to console
	//(I know there are some fancy non-fixed array-like things ('vectors'?) but I don't want to use them until I'm sure I can make them work)
	//(And they seem inefficient for our current purposes, though obviously necessary when we start changing N)
	//failed code below
	//cout << "Please enter the number of particles desired";
	//cin >> num_particles;


	//hard sphere cycles
	////COMMENT OUT HERE - hard spheres

	//cout << "Running hard spheres.\n";
	//posn positions_sph[num_particles];

	////test: run hard spheres in 3D and get xyz file
	//MC_loop_hardsph(positions_sph, num_particles, 343.0, num_prod_cycles_hardsph, .001, true, "testloop_posns_hard_sphere");

	//cout << "done!";


	//Lennard-Jones cycles

	//first- LJ test and energy grabbing

	//run w scaled density = 1., scaled temperature = 2.
	//note we will print all lattices to file for visualization
	//and calculate energy, but not force just yet

	//initialize position and energy list

	////COMMENT OUT FROM HERE - energy calc
	//	double energies[num_prod_cycles_LJ];
	//	posn positions_LJ[num_particles];
	//
	//	//scale down random moves to reduce rejection/increase convergence
	//	double rand_move_scale_energycalc = .1;
	//	
	//	double scaled_r_cutoff = 3.0;
	//	MC_loop_LJ(num_prod_cycles_LJ, positions_LJ, num_particles, 1., 2., 1.0, 1.0, true, true, scaled_r_cutoff, true, "testloop_posns_LJ", true, energies,false, default_position_array, default_input_double_array, rand_move_scale_energycalc);
	//	
	//	int i;
	//	for (i = 0; i < num_prod_cycles_LJ; i++) {
	//		cout << energies[i];
	//		cout << "\n";
	//	}
	//	//save all cycles' Es to file, with summary stats (avg, stdev)
	//	//calc summary energy values (avg, stdev) to save to file
	//	
	//	double E_avg = 0.0;
	//	double Esq_avg = 0.0;
	//	for (int index = 0; index < num_prod_cycles_LJ; index++) {
	//		E_avg += energies[index] / num_prod_cycles_LJ;
	//		Esq_avg += energies[index]* energies[index] / num_prod_cycles_LJ;
	//	}
	//	double E_variance = Esq_avg - E_avg*E_avg;
	//	double stdev_E_per_particle = (1/num_prod_cycles_LJ) * sqrt(E_variance / num_prod_cycles_LJ);
	//	
	//	string energy_filename = "energies_from_LJ_energy_test_" + to_string(num_particles) +"_particles";
	//	string E_header = "T* = 2.0, avg energy = " + to_string(E_avg) + "\n or " + to_string(E_avg/num_particles) + " per particle,\n stdev energy per particle= " + to_string(stdev_E_per_particle) +"\n" ;
	//	print_double_array_to_file(energy_filename, energies, num_prod_cycles_LJ, E_header, true);
	//


	////COMMENTED OUT FOR NOW - density variation
	//posn positions_LJ_density[num_particles];
	//posn force_holder[num_particles];

	//const int num_densities_to_test = 3;
	//double scaled_densities[num_densities_to_test] = { .1, .5, .9 };

	//double virial_pressure_array[num_prod_cycles_LJ];
	//double energy_holder_array[num_prod_cycles_LJ];


	//double scaled_r_cutoff_rhos = 3.0;

	//double curr_scaled_density;
	//double curr_scaled_kinetic_pressure;
	//double scaled_T = 2.;

	////want to scale down random moves at higher rho to keep acceptance rate reasonable
	//double rand_move_scales[3] = { 1.0, .7, .2 };
	//double curr_rand_scalefac;


	//for (int i_rho = 0; i_rho < num_densities_to_test; i_rho++) {
	//	curr_scaled_density = scaled_densities[i_rho];
	//	curr_scaled_kinetic_pressure = curr_scaled_density * scaled_T;

	//	curr_rand_scalefac = rand_move_scales[i_rho];

	//	string position_filename = "testloop_posns_LJ_" + to_string(num_prod_cycles_LJ) + "_" + to_string(i_rho);



	//	MC_loop_LJ(num_prod_cycles_LJ, positions_LJ_density, num_particles, curr_scaled_density, scaled_T, 1.0, 1.0, true, true, scaled_r_cutoff_rhos, true, position_filename, true, energy_holder_array, true, force_holder, virial_pressure_array, curr_rand_scalefac);

	//	//write pressures to file
	//	string header = "Scaled density = ";
	//	header += to_string(curr_scaled_density);
	//	header += ": scaled kinetic pressure = ";
	//	header += to_string(curr_scaled_kinetic_pressure);

	//	string filename = "LJ_pressures_rho_variation";

	//	print_double_array_to_file(filename, virial_pressure_array, num_prod_cycles_LJ, header, true);

	//	//get average of F dot rs = P_virial_no_tail
	//	double curr_P_virial = 0.0;
	//	for (int j_LJ = 0; j_LJ < num_prod_cycles_LJ; j_LJ++) {
	//		curr_P_virial += virial_pressure_array[j_LJ];
	//	}
	//	curr_P_virial = curr_P_virial / num_prod_cycles_LJ;

	//	//get variance of P
	//	double curr_variance_P = 0.0;
	//	for (int j_LJ = 0; j_LJ < num_prod_cycles_LJ; j_LJ++) {
	//		curr_variance_P += pow(virial_pressure_array[j_LJ] - curr_P_virial, 2.0);
	//	}
	//	//variance = E[(x - xbar)^2] = SUM[(x-xbar)^2]/n
	//	curr_variance_P = curr_variance_P / num_prod_cycles_LJ;

	//	//stdev P = stdev mean = sqrt(var/n)
	//	double stdev_P = sqrt(curr_variance_P / num_prod_cycles_LJ);


	//	double curr_P_scaled = curr_P_virial + curr_scaled_kinetic_pressure;

	//	double curr_Ps[4] = { curr_P_scaled, stdev_P, curr_P_virial, curr_scaled_kinetic_pressure };

	//	string pressure_filename = "LJ_pressures_rhovariation_longrun_justpressures";

	//	string P_header = "rho = " + to_string(scaled_densities[i_rho]) + ", num cycles = " + to_string(num_prod_cycles_LJ);
	//	P_header += ": P, stdev, virial, kinetic";

	//	print_double_array_to_file(pressure_filename, curr_Ps, 4, P_header, true);

	//	//write energies to file - just so I have more energies to look at

	//	double E_avg = 0.0;
	//	double Esq_avg = 0.0;
	//	for (int index = 0; index < num_prod_cycles_LJ; index++) {
	//		E_avg += energy_holder_array[index] / num_prod_cycles_LJ;
	//		Esq_avg += energy_holder_array[index] * energy_holder_array[index] / num_prod_cycles_LJ;
	//	}
	//	double E_variance = Esq_avg - E_avg*E_avg;
	//	double stdev_E_avg = sqrt(E_variance / num_prod_cycles_LJ);

	//	string energy_filename = "energies_from_rho_variation_w_" + to_string(num_particles);
	//	string E_header = "T* = " + to_string(scaled_T) + ", avg energy = " + to_string(E_avg) + "or " + to_string(E_avg / num_particles) + " per particle, stdev energy avg= " + to_string(stdev_E_avg); ;
	//	print_double_array_to_file(energy_filename, energy_holder_array, num_prod_cycles_LJ, E_header, true);


	//	//clear force array before using in next loop
	//	clear_posn_array(force_holder, num_particles);
	//}

	//Argon comparison
	//based on literature values in Gilgen, R., R. Kleinrahm, and W. Wagner. "Measurement and correlation of the (pressure, density, temperature) relation of argon I. The homogeneous gas and liquid regions in the temperature range from 90 K to 340 K at pressures up to 12 MPa." The Journal of Chemical Thermodynamics 26, no. 4 (1994): 383-398.
	//(newer and seem somewhat more precise (error in density ~ .01 to .1%, error in T,P < .001%

	//point to approximate: T = 325.00 K, rho = 167.129 kg/m^3, P = 11.00868 MPa 
	//						T* = 2.72651006711, rho* = .09946, P*= .2641 
	//						note P* = P*_kinetic - .00707869127; P*_virial = - .00707869127
	// this point was chosen for having the second-highest temperature studied
	// (since at low temperatures the LJ simulation tends not to accept very frequently at all)
	// (and the highest temperature studied may have been one with more error than the others, since it was the last possible T to study
	// It also had the second-highest rho and P out of the points studied at this T
	// which allows us to use the adjacent points to get an expected bounding range
	// these adjacent points are about 10% of rho (and 10% of P) away, so we will test variations
	// within this range


	////COMMENT OUT FROM HERE - Ar cycles
	//	double scaled_T_Ar = 2.72651006711;
	//	double scaled_rho_Ar = .09946;
	//
	//	//coarser initial search
	//	//constexpr int num_rho_perturbations = 5;	
	//	//double scaled_rho_Ar_perturbations[num_rho_perturbations] = { -.01, -.005, 0.0, .005, .01};
	//	
	//	//found that leaving rho*_lit unchanged gave by far the closest P* to P*_lit
	//	//however, this value was very slightly (~.5%) above the literature value
	//	//since P* was fairly consistently linear in rho* at roughly the same % rate of change
	//	//(eg a 5% increase in rho* as performed above led to about a 6% increase in P*)
	//	//I further tested values about .5% below the literature value
	//	constexpr int num_rho_perturbations = 2;	
	//	double scaled_rho_Ar_perturbations[num_rho_perturbations] = { -.0005, -.0004};
	//
	//
	//
	//	posn positions_Ar[num_particles];
	//	posn force_holder_Ar[num_particles];
	//
	//	double virial_pressure_array_Ar[num_prod_cycles_Ar];
	//	double energy_holder_array_Ar[num_prod_cycles_Ar];
	//	double scaled_r_cutoff_Ar = 3.;
	//
	//
	//	double curr_scaled_rho;
	//	double curr_scaled_P_kinetic;
	//
	//	for (int i_Ar = 0; i_Ar < num_rho_perturbations; i_Ar++) {
	//		curr_scaled_rho = scaled_rho_Ar + scaled_rho_Ar_perturbations[i_Ar];
	//		curr_scaled_P_kinetic = curr_scaled_rho*scaled_T_Ar;
	//		
	//		MC_loop_LJ(num_prod_cycles_Ar, positions_Ar, num_particles, curr_scaled_rho, scaled_T_Ar, 1.0, 1.0, true, true, scaled_r_cutoff_Ar, true, "Ar_posns_LJ" + to_string(i_Ar), true, energy_holder_array_Ar, true, force_holder_Ar, virial_pressure_array_Ar);
	//		//write pressures to file
	//		string header = "Scaled density = ";
	//		header += to_string(curr_scaled_rho);
	//		header += ": scaled kinetic pressure = ";
	//		header += to_string(curr_scaled_P_kinetic);
	//
	//		string filename = "LJ_pressures_Ar_approximation_" + to_string(num_prod_cycles_Ar) + "run";
	//
	//		print_double_array_to_file(filename, virial_pressure_array_Ar, num_prod_cycles_Ar, header, true);
	//
	//		//get average of F dot rs = P_virial_no_tail
	//		double curr_P_virial = 0.0;
	//		for (int j_Ar = 0; j_Ar < num_prod_cycles_Ar; j_Ar++) {
	//			curr_P_virial += virial_pressure_array_Ar[j_Ar];
	//		}
	//		curr_P_virial = curr_P_virial / num_prod_cycles_Ar;
	//
	//		//get variance of P
	//		double curr_variance_P = 0.0;
	//		for (int j_Ar = 0; j_Ar < num_prod_cycles_Ar;  j_Ar++) {
	//			curr_variance_P += pow(virial_pressure_array_Ar[j_Ar] - curr_P_virial, 2.0);
	//		}
	//		//variance = E[(x - xbar)^2] = SUM[(x-xbar)^2]/n
	//		curr_variance_P = curr_variance_P / num_prod_cycles_Ar;
	//
	//		//stdev P = stdev mean = sqrt(var/n)
	//		double stdev_P = sqrt(curr_variance_P / num_prod_cycles_Ar);
	//
	//
	//		double curr_P_scaled = curr_P_virial + curr_scaled_P_kinetic;
	//
	//		double curr_Ps[4] = { curr_P_scaled, stdev_P, curr_P_virial, curr_scaled_P_kinetic };
	//
	//		string pressure_filename = "LJ_pressures_Ar_approximation_longrun_justpressures";
	//
	//		string P_header = "rho perturbation = " + to_string(scaled_rho_Ar_perturbations[i_Ar]) + ", num cycles = " + to_string(num_prod_cycles_Ar);
	//		P_header += ": P, stdev, virial, kinetic";
	//
	//		print_double_array_to_file(pressure_filename, curr_Ps, 4, P_header, true);
	//
	//		//clear force array before using in next loop
	//		clear_posn_array(force_holder_Ar, num_particles);
	//	}


	//	string dummy_input;
	//	cin >> dummy_input;

}