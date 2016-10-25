#include <stdio.h>
#include <iostream>
using namespace std;

struct posn3d {
	double x, y, z;
};

void initialize_posns(posn3d* input_lattice, int numpts) {
	for (int i = 0; i < numpts; i++) {
		posn3d currpt;
		currpt.x = .1 * (i^2);
		currpt.y = 0;
		currpt.z = 0;
		input_lattice[i] = currpt;
	}
}


int sum_intarr(int* intlist, int numelems) {
	int curr_int;
	int intsum = 0;
	for (int i = 0; i < numelems; i++) {
		curr_int = intlist[i];
		intsum += curr_int;
	}
	return intsum;

}

int main()
{
	const int num_posns = 5;
	posn3d posn_lattice[num_posns];

	initialize_posns(posn_lattice, num_posns);


	for (int i = 0; i < num_posns; i++) {
		cout << (posn_lattice[i].x);
		cout << "\n";
	}


	//int a[3] = { 1,2,3 };
	//int i;
	//int asum = sum_intarr(a, 3);
	//cout << asum;

	char jab;
	cin >> jab;
	
	
	return 0;
}
