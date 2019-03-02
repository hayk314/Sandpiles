#include <iostream>
#include "output.h"
#include "Abel.h"


using namespace std;
using namespace outputFunctions;
using namespace Abel;


int main()
{
	unsigned int n;
	cout << "Please write the number of grains on the origin" << endl;
	cin >> n;

	//double t = Abel::moveStandard(n);
	double t = Abel::moveStandard_1Step(n);

	cout << "done in " << t << "seconds" << endl;

	char _end_;
	cout << endl << "Enter any character to quit";
	cin >> _end_;

	return 1;
}