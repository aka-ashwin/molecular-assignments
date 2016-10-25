#include <iostream>
#include <limits>
using namespace std;

// function for initializing arrays as per Metropolis et al

//lattice parameters

struct lattice {
	double* xs;
	double* ys;
};


lattice fill_Metropolis_lattice() {
	int num_particles = 224;
	int num_rows = 16;
	int num_cols = 14;
	double xstep = 1.0/(float)num_cols;
	double offset = 1.0/(2.0 * (float)num_cols);
	double ystep = 1.0 / (float) num_rows;
	cout << xstep;
	cout << '\n';
	cout << offset;
	cout << '\n';
	cout << ystep;
	cout << '\n';

	//initialize and fill lattice
	double xs[224];
	double ys[224];

	int lattice_index = 0;
	double curr_x = 0.0;
	double curr_y = 0.0;
	for (int row_index = 0; row_index < num_rows; row_index++) {
		curr_x = 0.0;
		if(row_index % 2 == 1) { curr_x += offset; }

		for (int col_index = 0; col_index < num_cols; col_index++) {
			xs[lattice_index] = curr_x;
			ys[lattice_index] = curr_y;
			lattice_index++;
			cout << lattice_index;
			cout << ": ";
			cout << curr_x;
			cout << ", ";
			cout << curr_y;
			cout << "\n";
			curr_x += xstep;
		}
		curr_y += ystep;
	}
	//make lattice to wrap around values and return
	lattice metro_lattice;
	metro_lattice.xs = xs;
	metro_lattice.ys = ys;
	return metro_lattice;
}


void main() {


	lattice aleph_test = fill_Metropolis_lattice();

	cout << "main data \n";
	for (int i = 100; i >= 0; i--) {
		cout << i;
		cout << ": ";
		cout << aleph_test.ys[i];
		cout << "\n";
	}

	std::cout << "Press ENTER to continue...";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

}

