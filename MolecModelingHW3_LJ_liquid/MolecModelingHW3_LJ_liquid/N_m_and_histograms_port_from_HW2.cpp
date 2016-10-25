//Parts 6-7 helper: updating N_m
#include <math.h>


struct Nm_arr {
	int counts[64];
};

//making a global variable for now; can't get function-passing to work properly (turns array to pointer, which causes problems)
Nm_arr Nms;

void update_N_ms(posn* positions, int num_particles, Nm_arr local_Nms, double d_0, double bin_size) {


	// note bin_size = "deltaA^2" of notes * 1/pi; this makes the expressions a little tidier 
	posn particle_i = posn(0., 0., 0.);
	posn particle_j = posn(0., 0., 0.);

	double sq_dist_scaling = 1.0 / bin_size;
	double sq_dist_const = pow(d_0, 2.);


	for (int i = num_particles - 1; i > 0; i--) {
		particle_i = positions[i];
		for (int j = i - 1; j >= 0; j--) {
			particle_j = positions[j];
			//get sq dist
			double ij_dist_sq = pow(minimum_image_distance(particle_i, particle_j), 2.0);

			//bin sq dist
			// - .5 ensures that int-casting creates bins such that e.g. bin 0 counts values from d_0^2 to d_0^2 + 1 increment
			// rather than d_0^2 to d_0^2 + .5, as would otherwise be the case.
			int ij_dist_index = (int)(((ij_dist_sq - sq_dist_const)*sq_dist_scaling) - .5);
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
		prod_mean_value += pow(Nms.counts[i], 2)*i;
		bin_mean_value += Nms.counts[i] * i;
		count_mean_value += pow(Nms.counts[i], 2);
		total_entries_count += Nms.counts[i];
	}
	double avging_coeff = 1. / ((float)total_entries_count);

	prod_mean_value *= avging_coeff;
	bin_mean_value *= avging_coeff;
	count_mean_value *= avging_coeff;


	double cov = prod_mean_value - bin_mean_value * count_mean_value;

	double bin_variance = get_histogram_bin_moment(Nms, 2) - pow(get_histogram_bin_moment(Nms, 1), 2);
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
