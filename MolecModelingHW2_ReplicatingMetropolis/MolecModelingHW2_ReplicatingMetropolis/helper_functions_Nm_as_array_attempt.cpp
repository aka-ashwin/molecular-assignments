#include <iostream>
#include <limits>
#include <string>
#include <ctime>
#include "C:/Users/Ashwin/Desktop/Molec Modeling/random_mars.h"
using namespace std;



// need:
// Part 1 : to write initial lattice to .dat file
		// and make an image of it - use VMD?



// function for initializing arrays as per Metropolis et al


// Part 1: initialize and print (and get image)
//initialize lattice; will use xs and ys as globally-accessed posn array
// (I tried to use a struct that wrapped around the arrays and pass it, but I couldn't get it to work)
// (because only pointers can be passed, but passing pointers rather than arrays seems to cause errors with array lookup
struct posn {
	double x, y;
};


int num_particles = 224;
int num_rows = 16;
int num_cols = 14;
posn positions[224];

//fill in lattice: requires global num_particles, num_rows, num_cols integers
// as well as global positions[num_particles] array.
// assumes 2D
void fill_Metropolis_lattice() {
	// lattice parameters
	double xstep = 1.0/(float)num_cols;
	double offset = 1.0/(2.0 * (float)num_cols);
	double ystep = 1.0 / (float) num_rows;

	// loop through and fill lattice
	int lattice_index = 0;
	double curr_x = 0.0;
	double curr_y = 0.0;
	for(int row_index = 0; row_index < num_rows; row_index++) {
		curr_x = 0.0;
		if(row_index % 2 == 1) { curr_x += offset; }

		for (int col_index = 0; col_index < num_cols; col_index++) {
			positions[lattice_index].x = curr_x;
			positions[lattice_index].y = curr_y;
			lattice_index++;
			curr_x += xstep;
		}
		curr_y += ystep;
	}
}

// print lattice: currently just prints to console

void print_lattice() {
	for (int i = num_particles-1; i >= 0; i--) {
		cout << i;
		cout << ": ";
		cout << positions[i].x;
		cout << ",";
		cout << positions[i].y;
		cout << "\n";
	}

	//std::cout << "Press ENTER to continue...";
	//std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

}




// Part 2: calculate energy

//helper function: given particle i and j's positions, gets distance  particle i to nearest image of particle j
//assumes box length in x and y dims is 1

double vector_length(double x, double y) {
	return sqrt(pow(x, 2.0) + pow(y, 2.0));
}

double minimum_image_distance(posn particlefrom, posn particleto) {
	double xfrom = particlefrom.x;
	double yfrom = particlefrom.y;
	double xto = particleto.x;
	double yto = particleto.y;
	double xdist = xto - xfrom;
	double ydist = yto - yfrom;

	double min_xdist = xdist;
	double min_ydist = ydist;

	if (xdist > .5) { min_xdist -= 1.0; }
	else if(xdist < -.5) { min_xdist += 1.0; }

	if (ydist > .5) { min_ydist -= 1.0; }
	else if (ydist < -.5) { min_ydist += 1.0; }

	double minim_distance = vector_length(min_xdist, min_ydist);
	return minim_distance;
}

//energy-calculation function: since energy is either 0 (no hard spheres overlap) or infinite (some overlap), we just have a boolean
	//which is true if the configuration is acceptable (energy 0), and false if not (energy infinite)
	//we can also stop looping through as soon as we have found an infinite-energy pair
// doesn't take input; just acts on xs, ys
bool get_configuration_validity(double sph_radius) {
	bool config_is_valid = true;
	posn i_posn;
	posn j_posn;
	// loop through all i values > 1; note we don't need i = 0, since the 1-0 interaction is covered by i = 1, j = 0.
	for (int i = num_particles-1; i >0; i--) {
		i_posn = positions[i];
		// loop through all j-values below i; this prevents double-counting
		for (int j = i-1; j >= 0; j--) {
			j_posn = positions[j];
			double min_distance = minimum_image_distance(i_posn, j_posn);
			if (min_distance <= sph_radius) { 
				config_is_valid = false; 
				break;
			}

		}

	
	}
	return config_is_valid;
}

//Part 3 - Monte Carlo implementation

//propose a step: select an atom uniformly at random, generate e_x and e_y uniformly at random from [-1,1]
//and move it by dx = (d - d_0)*ex and dy = (d - d_0)*ey, where d = 1/14 and d_0 is an input value
//respecting periodic boundary conditions

//constants used in MC step: box length, parameter d
// note box length should be used in periodic BC step of 'energy calculation'; currently is not
double box_len = 1.0;
double d = 1. / 14.;
//will leave random number generation to MC cycle function (and seeding to MC loop function)
void MC_update(int part_index, double x_move_rand, double y_move_rand, double d_0) {
	double alpha = d - d_0;
	
	posn old_part = positions[part_index];

	double curr_part_x = old_part.x;
	double curr_part_y = old_part.y;

	double proposed_x = curr_part_x + alpha*x_move_rand;
	double proposed_y = curr_part_y + alpha*y_move_rand;

	//apply pbcs; note that x = 0.0 or y = 0.0 are currently valid; I thought since they were starting positions they probably ought to be
	//potential issue: assumes proposal displacements are small enough that they never move particles more than a box length
	//a reasonable assumption here and for most realistic displacements, but not good to have as a general case.
	if (proposed_x >= box_len) { proposed_x -= box_len; }
	else if (proposed_x < 0) { proposed_x += box_len; }
	if (proposed_y >= box_len) { proposed_y -= box_len; }
	else if (proposed_y < 0) { proposed_y += box_len; }

	posn proposed_part;
	proposed_part.x = proposed_x;
	proposed_part.y = proposed_y;

	//test proposition; note we are using d_0 = 1/14 as the sphere radius
	positions[part_index] = proposed_part;
	bool config_validity = get_configuration_validity(d_0);
	//if configuration is valid, keep it (this is a very simple acceptance step, since our proposal generation and the energy landscape are so simple)
	//if invalid, discard it and take the particle back to its old position
	if (!config_validity) {
		positions[part_index] = old_part;
	}
}


// Monte Carlo cycle: run update process n times, where n is the number of particles
void MC_cycle(RanMars *rng, double d_0) {
	for (int i = num_particles-1; i > 0; i--) {
		//generate x and y move RVs; RanMars.uniform gives unif(0,1) and we want unif(-1,1), so we multiply by 2 and subtract 1.
		double x_move = (2 * rng->uniform()) - 1;
		double y_move = (2 * rng->uniform()) - 1;

		//choose particle to move
		//generate uniform double between 0 and num_particles; 
			//for each particle i from 0 to num_particles-1 there is probability 1/num_particles of being between i and i + 1
			//so casting the result to an int gives us a uniform integer between 0 and num_particles
		int part_index = (int)(num_particles * rng->uniform());

		MC_update(part_index, x_move, y_move, d_0);
	}
}




//Parts 6-7 helper: updating N_m

struct Nm_arr {
	int counts[64];
};


Nm_arr default_Nm_arr;

//making a global variable for now; can't get function-passing to work properly (turns array to pointer, which causes problems)
int Nms[64];

void update_N_ms(Nm_arr local_Nms, double d_0, double bin_size) {


	// note bin_size = "deltaA^2" of notes * 1/pi; this makes the expressions a little tidier 
	posn particle_i;
	posn particle_j;

	double sq_dist_scaling = 1.0 / bin_size;
	double sq_dist_const = pow(d_0, 2.);


	for (int i = num_particles-1; i > 0; i--) {
		particle_i = positions[i];
		for (int j = i - 1; j >= 0; j--) {
			particle_j = positions[j];
			//get sq dist
			double ij_dist_sq = pow(minimum_image_distance(particle_i, particle_j), 2.0);
			
			//bin sq dist
			// - .5 ensures that int-casting creates bins such that e.g. bin 0 counts values from d_0^2 to d_0^2 + 1 increment
			// rather than d_0^2 to d_0^2 + .5, as would otherwise be the case.
			int ij_dist_index = (int)(((ij_dist_sq -sq_dist_const)*sq_dist_scaling) - .5);
			if (ij_dist_index <= 63) {
				Nms[ij_dist_index] += 1;
			}
			

		}


	}

}


double get_histogram_count_moment(Nm_arr local_Nms, int moment) {
	double moment_value = 0.0;
	int total_entries_count = 0;
	for (int i = 63; i >= 0; i--) {
		moment_value += Nms[i] * pow(Nms[i], moment);
		total_entries_count += Nms[i];
	}
	moment_value = moment_value / total_entries_count;
	return moment_value;
}

double get_histogram_bin_moment(Nm_arr local_Nms, int moment) {
	double moment_value = 0.0;
	int total_entries_count = 0;
	for (int i = 63; i >= 0; i--) {
		moment_value += Nms[i] * pow(i, moment);
		total_entries_count += Nms[i];
	}
	moment_value = moment_value / total_entries_count;
	return moment_value;
}


double get_histogram_bin_count_correl(Nm_arr local_Nms) {
	double prod_mean_value = 0.0;
	double bin_mean_value = 0.0;
	double count_mean_value = 0.0;

	int total_entries_count = 0;
	for (int i = 63; i >= 0; i--) {
		prod_mean_value += pow(Nms[i],2)*i;
		bin_mean_value += Nms[i] * i;
		count_mean_value += pow(Nms[i], 2);
		total_entries_count += Nms[i];
	}
	double avging_coeff = 1. / ((float)total_entries_count);

	prod_mean_value *= avging_coeff;
	bin_mean_value *= avging_coeff;
	count_mean_value *= avging_coeff;


	double cov = prod_mean_value - bin_mean_value * count_mean_value;

	double bin_variance = get_histogram_bin_moment(Nms, 2);
	double correl_value = cov / bin_variance;
	return correl_value;
}

//simple linear regression: optimal slope b = correlation, optimal intercept a = mean of y - b*(mean of x) 
//the way I'm indexing the bins, my bin indices are equal to 1 plus "m" from the notes
//so I need N_{index = -1/2}, while the intercept I'm calculating should be N_{index = 0)
//so N_needed = intercept -.5*slope = count_mean - b*(bin_mean + .5)


double get_histogram_min_edge_estimate(Nm_arr local_Nms) {


	double slope = get_histogram_bin_count_correl(Nms);
	double bin_mean = get_histogram_bin_moment(Nms, 1);
	double count_mean = get_histogram_count_moment(Nms, 1);
	double intercept = count_mean - slope * bin_mean;

	double min_edge_estimate = intercept - slope * .5;

	return min_edge_estimate;
}


//Parts 3 on: running loops

// Monte Carlo loop: initialize lattice, run a few equilibration cycles, then run production cycles
void MC_loop(double d_0, int num_prod_cycles, bool calc_Nms = false, Nm_arr local_Nms = default_Nm_arr, double Nm_bin_size = 1.0, bool print_during_production = false) {
	int	num_equil_cycles = 5;

	//initialize RNG
	RanMars *random;
	//seed with time
	random = new RanMars(1005);


	//initialize lattice
	fill_Metropolis_lattice();

	//equilibration cycles 
	for (int i = 0; i < num_equil_cycles; i++) {
		MC_cycle(random, d_0);
	}

	//production cycles
	for (int j = 0; j < num_prod_cycles; j++) {
		MC_cycle(random, d_0);
		if (print_during_production) {
			print_lattice();
		}
		if (calc_Nms) {
			update_N_ms(local_Nms, d_0, Nm_bin_size);
		}
	}
}




// 
void main() {
	int	num_prod_cycles = 5;

	////run w test nus; print out final positions for each sim
	//int test_nus[] = {2,5,7};
	//double curr_d_0;
	//for(int i = 0; i < 3; i++) {
	//	curr_d_0 = d*(1.0 - pow(2.0, test_nus[i]-8.));
	//	cout << "test production cycle: for d_0 = ";
	//	cout << curr_d_0;
	//	cout << "\n";
	//	MC_loop(curr_d_0, num_prod_cycles);
	//	print_lattice();

	//run w/ nu = 5, calculate N_m
	//initialize Nms

	Nm_arr local_Nms;

	//global Nms hack
	for (int i = 63; i >= 0; i--) {
		local_Nms.counts[i] = 0;
		Nms[i] = 0;
		
	}
	//fix K
	float K = 1.5;

	//loop w/ nu = 5
	int	curr_nu = 5;
	//calc d_0 
	double curr_d_0 = d*(1.0 - pow(2.0, curr_nu - 8 ));
	//calc area increment; will use "delta_A^2" from assignment divided by pi, = (K^2 -1)d_0^2/64
	double scaled_dA_squared = (double) (pow(K, 2) - 1)*pow(curr_d_0, 2) / 64.;
	
	//run cycles
	MC_loop(curr_d_0, num_prod_cycles, true, local_Nms,scaled_dA_squared);


	//calculations from Nms

	double N_critical = get_histogram_min_edge_estimate(Nms);

	//get nbar analogue PA/NkT = N_critical * 64/(N^2 * (K^2 - 1)), where N = number of particles
	double nbar_analogue = N_critical * 64. / (pow(num_particles, 2) * (pow(K, 2) - 1));

	double mean_count = get_histogram_count_moment(Nms,1);
	double mean_bin = get_histogram_bin_moment(Nms, 1);
	double var_bin = get_histogram_bin_moment(Nms, 2);
	
	cout << "mean count = ";
	cout << mean_count;

	cout << "nbar analogue = ";
	cout << nbar_analogue;

	cout << "Nms = ";
	cout << Nms;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	


}