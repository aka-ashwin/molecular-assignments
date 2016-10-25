#include <iostream>
#include <string>
#include "time.h"
#include "stdlib.h"
#include <fstream>

using namespace std;

//function delcarations
void print_box(int, int, int, int);
void clear_screen();
void pause(double);

//==================================================
// Method to Edit
//==================================================
int main() {
	int boxlength; //length of box
	int x, vx; //x position and velocity 
	int y, vy; //y position and velocity 
	int nsteps; //total length of simulation
	int step; //current simulation step
	ofstream myfile("position.dat", ios::out);

	//WRITE CODE HERE
	//set up simulation parameters
	boxlength = 20;
	nsteps = 100;
	//initialize particle
	x = 10;
	vx = 1;
	y = 5;
	vy = 1;
	//loop simulation
	int newx;
	int newy;
	for (step = 0; step < nsteps; step++) {
		//numerical integration
		newx = x + vx;
		newy = y + vy;

		//boundary conditions: 
		//if update would take particle outside of bounds,
		//	instead take it to opposide side of box

		x = newx % boxlength;
		y = newy % boxlength;


		//write stepnumber, position to file
		myfile << step << ", " << x << ", " << y << endl;

		//print box
		print_box(boxlength, step, x, y);

		//
	}




	myfile.close();
	return 0;
}

//==================================================
// DO NOT EDIT CODE BELOW THIS LINE
//==================================================

//function definitions
void print_box(int boxl, int step, int x, int y) {
	// function that prints the 1-dimensional 'box'
	// input arguements:
	//   boxl - length of the box
	//   step - timestep of simulation
	//   x - x position 
	//   y - y position 

	if ((x < 0) || (x >= boxl) || (y<0) || (y >= boxl)) {
		cout << "Error! Particle position is outside of box! Boundary Conditions are likely incorrect." << endl;
		exit(1);
	}

	int i, j;
	clear_screen();
	cout << "Step " << step << endl;
	for (j = 0; j< boxl; j++) {
		if (j == 0) {
			for (i = -1; i <= boxl; i++)
				cout << "-";
			cout << endl;
		}
		cout << "|";
		for (i = 0; i < boxl; i++) {
			if ((i == x) && (j == y)) {
				cout << "x";
			}
			else {
				cout << " ";
			}
		}
		cout << "|" << endl;
		if (j == boxl - 1) {
			for (i = -1; i <= boxl; i++)
				cout << "-";
			cout << endl;
		}
	}
	pause(0.1);
}
void clear_screen() {
	//function that 'clears' the screen by printing 100 new lines
	cout << string(100, '\n');
}
void pause(double s) {
	//function that pauses the code for s seconds
	//otherwise, the output is written to the screen too quickly
	clock_t ti, t;
	ti = clock();
	double diff = 0;
	while (diff < s) {
		t = clock() - ti;
		diff = (float)t / CLOCKS_PER_SEC;
	}
}
