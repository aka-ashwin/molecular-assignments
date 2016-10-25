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



// function for initializing arrays as per Metropolis et al


// Part 1: initialize and print (and get image)
//initialize lattice; will use xs and ys as globally-accessed posn array
// (I tried to use a struct that wrapped around the arrays and pass it, but I couldn't get it to work)
// (because only pointers can be passed, but passing pointers rather than arrays seems to cause errors with array lookup)
struct posn {
	double x, y, z;
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

void print_lattice_to_console() {
	for (int i = num_particles-1; i >= 0; i--) {
		cout << i;
		cout << ": ";
		cout << positions[i].x;
		cout << ",";
		cout << positions[i].y;
		cout << ",";
		cout << positions[i].z;
		cout << "\n";
	}

	//std::cout << "Press ENTER to continue...";
	//std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

}

void print_lattice_to_file(string filename, string filehead = "index, x, y", bool has_header = false, string filetype = "xyz") {
	string separator = ", ";
	string filepath = "outputs/" ;
	string full_filename = filepath + (filename + ("." + filetype));
	ofstream myfile;
	myfile.open(full_filename, ios::out);
	if (has_header | filetype == "xyz") {
		myfile << filehead + "\n";
	}
	int currInd;
	posn currPosn;
	for (int i = num_particles-1; i >= 0; i--) {
		currInd = num_particles - (i + 1);
		currPosn = positions[currInd];
		if (filetype == "xyz") {
			myfile << "Ar";
			myfile << separator;
		}
		myfile << currPosn.x;
		myfile << separator;
		myfile << currPosn.y;
		myfile << separator;
		myfile << currPosn.z;
		myfile << "\n";
	}
	myfile.close();


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

//making a global variable for now; can't get function-passing to work properly (turns array to pointer, which causes problems)
Nm_arr Nms;

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
				Nms.counts[ij_dist_index] += 1;
			}
			

		}


	}

}


double get_histogram_count_moment(Nm_arr local_Nms, int moment) {
	double moment_value = 0.0;
	int total_entries_count = 0;
	for (int i = 63; i >= 0; i--) {
		moment_value += Nms.counts[i] * pow(Nms.counts[i], moment);
		total_entries_count += Nms.counts[i];
	}
	moment_value = moment_value / total_entries_count;
	return moment_value;
}

double get_histogram_bin_moment(Nm_arr local_Nms, int moment) {
	double moment_value = 0.0;
	int total_entries_count = 0;
	for (int i = 63; i >= 0; i--) {
		moment_value += Nms.counts[i] * pow(i, moment);
		total_entries_count += Nms.counts[i];
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
		prod_mean_value += pow(Nms.counts[i],2)*i;
		bin_mean_value += Nms.counts[i] * i;
		count_mean_value += pow(Nms.counts[i], 2);
		total_entries_count += Nms.counts[i];
	}
	double avging_coeff = 1. / ((float)total_entries_count);

	prod_mean_value *= avging_coeff;
	bin_mean_value *= avging_coeff;
	count_mean_value *= avging_coeff;


	double cov = prod_mean_value - bin_mean_value * count_mean_value;

	double bin_variance = get_histogram_bin_moment(Nms, 2) - pow(get_histogram_bin_moment(Nms,1),2);
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

Nm_arr default_Nm_arr;
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
			print_lattice_to_console();
		}
		if (calc_Nms) {
			update_N_ms(local_Nms, d_0, Nm_bin_size);
		}
	}
}




// 
void main() {
	int	num_prod_cycles = 64;


	//initialize and print lattice
	cout << "building initial lattice...\n\n";
	fill_Metropolis_lattice();
	print_lattice_to_file("init", "initial lattice", false, "dat");

	//run w test nus; print out final positions for each sim, and whether they contain overlapping particles
	int test_nus[] = {2,5,7};
	double curr_d_0;
	string base_fname = "testcycles_nu_";
	string base_header = "test cycles: " + to_string(num_prod_cycles) + "production cycles with nu = ";
	string curr_nu_str;
	for (int i = 0; i < 3; i++) {
		curr_nu_str = to_string(test_nus[i]);
		curr_d_0 = d*(1.0 - pow(2.0, test_nus[i] - 8.));

		cout << "test production cycle: for d_0 = ";
		cout << curr_d_0;
		cout << "\n";
		MC_loop(curr_d_0, num_prod_cycles);

		bool lattice_validity = get_configuration_validity(curr_d_0);

		if (lattice_validity) { cout << "Success! Final position has no particles overlapping.\n"; }
		else { cout << "Failure! Something went wrong - final position has at least one pair of particles less than d_0 away from each other.\n"; }

		print_lattice_to_file(base_fname + curr_nu_str, base_header + curr_nu_str);
	}
	//run w/ nu = 5, calculate N_m
	//initialize Nms

	cout << "\nCalculating a single (PA/NkT) - 1 by regressing N_m on m \n";

	Nm_arr local_Nms;

	
	for (int i = 63; i >= 0; i--) {
		local_Nms.counts[i] = 0;
		Nms.counts[i] = 0;
		
	}
	//fix K
	float K = 1.5;

	//loop w/ nu = 5
	float curr_nu = 5;
	//calc d_0 
	curr_d_0 = d*(1.0 - pow(2.0, curr_nu - 8 ));
	//calc area increment; will use "delta_A^2" from assignment divided by pi, = (K^2 -1)d_0^2/64
	double scaled_dA_squared = (double) (pow(K, 2) - 1)*pow(curr_d_0, 2) / 64.;
	
	//run cycles
	MC_loop(curr_d_0, num_prod_cycles, true, local_Nms,scaled_dA_squared);


	//calculations from Nms



	double N_critical = get_histogram_min_edge_estimate(Nms);

	//get nbar analogue PA/NkT = N_critical * 64/(N * (K^2 - 1)), where N = number of particles
	//note: using N instead of N^2 because I am calculated N_m per cycle, not per update
	//I am also dividing by num_prod_cycles; 
		//since all derivations are linear in the N_m, for convenience I average over production cycles here
		//rather than earlier
	double nbar_analogue = N_critical * 64. / (num_prod_cycles * pow(num_particles, 1) * (pow(K, 2) - 1));

	//print N_m results; need to convert to using a function that takes an array and prints it to file.
	string filepath = "outputs/";
	string filetype = "dat";
	string filename = "Nms_for_ms_to_get_nbar_analogue";
	int num_Nm_bins = 16;

	string full_filename = filepath +  (filename + ("." + filetype));
	fstream myfile;
	myfile.open(full_filename, ios::out);
	//print header
	myfile << "data for nu=5: nbar analogue ((PA/Nkt)-1) calculated from this was ";
	myfile << nbar_analogue;
	myfile << "\n m, N_m \n";

	//loop through N_m array and print m, N_m for each m.
	int currInd;
	int currCount;
	for (int i = num_Nm_bins - 1; i >= 0; i--) {
		currInd = num_Nm_bins - (i + 1);
		currCount = Nms.counts[currInd];
		myfile << currInd;
		myfile << ", ";
		myfile << currCount;
		myfile << "\n";
	}
	myfile.close();
	
	//old writing process: prints to console
	//cout << "slope = ";
	//cout << get_histogram_bin_count_correl(Nms);
	//cout << "\n";

	//cout << "nbar analogue = ";
	//cout << nbar_analogue;
	//cout << "\n";

	//cout << "Nms = ";
	//cout << "\n";
	//for (int m = 63; m >= 0; m--) {
	//	cout << m;
	//	cout << ": ";
	//	cout << Nms.counts[m];
	//	cout << "\n";
	//}
	
	//Part 7 - loop over nu-values, using appropriate K-values, and get PA/NkT - 1; get graph of this versus (A/A_0) -1
	cout << "Calculating multiple (PA/NkT) - 1 values to compare them with (A/A_0) -1 \n\n";


	const int num_nus = 7;
	float nu_list[num_nus] = { 2.,4.,5.5,6.0,6.25,6.5,7.0 };
	float K_list[num_nus] = { 1.2, 1.3, 1.5, 1.6, 1.7, 1.8, 1.95 };

	double nbar_analogues[num_nus];

	double A_ratios[num_nus];
	//loop through nus; get corresponding A/A0s and simulated nbar analogues to plot
	for (int i = 0; i < 7; i++) {
		curr_nu = nu_list[i];
		K = K_list[i];

		//set up loop constants and run loops; same as for previous fixed-nu example.
		//calc d_0 
		double curr_d_0 = d*(1.0 - pow(2.0, curr_nu - 8));
		//calc area increment; will use "delta_A^2" from assignment divided by pi, = (K^2 -1)d_0^2/64
		double scaled_dA_squared = (double)(pow(K, 2) - 1)*pow(curr_d_0, 2) / 64.;

		//run cycles
		MC_loop(curr_d_0, num_prod_cycles, true, local_Nms, scaled_dA_squared);


		//calculations from Nms
		double N_critical = get_histogram_min_edge_estimate(Nms);

		//get nbar analogue PA/NkT = N_critical * 64/(N * (K^2 - 1)), where N = number of particles
		double nbar_analogue = N_critical * 64. / (num_prod_cycles * pow(num_particles, 1) * (pow(K, 2) - 1));
		nbar_analogues[i] = nbar_analogue;

		//get A/A_0 = 1/(sqrt(3) d_0^2 N /2)
		double A_ratio_inv = sqrt(3) * pow(curr_d_0, 2) * num_particles / 2.;
		A_ratios[i] = (1.0 / A_ratio_inv);


		//clear Nms
		for (int i = 63; i >= 0; i--) {
			Nms.counts[i] = 0;
		}
	}

		
	//print table of values

	string finalpart_filename = "nbar_analogues_vs_A_ratios";
	string full_finalpart_filename = filepath + (finalpart_filename +("." + filetype));
	myfile.open(full_finalpart_filename, ios::out);
	//print header
	myfile << "nu, (A/A_0) -1, (PA/Nkt) -1 \n";

	//loop through nu values and print nu, (A/A_0) -1, (PA/Nkt) -1  for each nu.

	float current_nu_in_loop;
	double current_A_ratio;
	double current_nbar_analogue;

	for (int i = num_nus-1; i >= 0; i--) {
		currInd = num_nus - (i + 1);
		current_nu_in_loop = nu_list[currInd];
		current_A_ratio = A_ratios[currInd];
		current_nbar_analogue = nbar_analogues[currInd];
		myfile << current_nu_in_loop;
		myfile << ", ";
		//need to subtract 1 to get A/A_0 -1
		myfile << current_A_ratio-1;
		myfile << ", ";
		myfile << current_nbar_analogue;
		myfile << "\n";
	}
	myfile.close();

	////old output - prints to file
	//for (int i = 0; i < 7; i++) {
	//	cout << "A/A_0: ";
	//	cout << A_ratios[i];
	//	cout << "\n";

	//	cout << "PA/NkT -1: ";
	//	cout << nbar_analogues[i];
	//	cout << "\n";
	//}




	// line to keep console open once program is done; keeps console output visible 
	std::cout << "Press ENTER to close...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}