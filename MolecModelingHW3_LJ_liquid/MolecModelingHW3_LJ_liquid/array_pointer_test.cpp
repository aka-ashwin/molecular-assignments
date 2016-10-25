#include <time.h>
#include <iostream>
#include "C:/Users/Ashwin/Desktop/Molec Modeling/random_mars.h"
using namespace std;

void main() {
	int i;
	unsigned int seed = unsigned int(time(NULL));
	cout << seed;
	cout << "\n";
	seed = abs((int)seed) % 900000000;
	cout << seed;
	cout << "\n";
	RanMars* random;
	random = new RanMars(seed);

	cin >> i;
}