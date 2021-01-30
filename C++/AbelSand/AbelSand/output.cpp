#include <iostream>
#include <fstream>
#include <string>

#include "output.h"

using namespace output_functions;
using namespace std;

//template<typename T_array>
BoxCoord output_functions::TrimmedArray(bool ** a, int sz1, int sz2)
{
    // trim any zero row/columns from the borders of a and return the resulting rectangle

	output_functions::BoxCoord x;
	x.i1 = -1; 
	x.i2 = -1; 
	x.j1 = -1; 
	x.j2 = -1;

	bool q = false;

	for (int i = 0; i < sz1; ++i){
		q = false;
		for (int j = 0; j < sz2; ++j)
			if (a[i][j] != 0){
				q = true;
				x.i1 = i;
				break;
			}
		if (q == true)
			break;
	}

	// i2
	for (int i = sz1 - 1; i >= 0; --i){
		q = false;
		for (int j = 0; j < sz2; ++j)
			if (a[i][j] != 0){
				q = true;
				x.i2 = i;
				break;
			}
		if (q == true)
			break;
	}
	
	// j1
	for (int j = 0; j < sz2; ++j){
		q = false;
		for (int i = 0; i < sz1; ++i)
			if (a[i][j] != 0){
				q = true;
				x.j1 = j;
				break;
			}
		if (q == true)
			break;
	}
	
	// j2
	for (int j = sz2 - 1; j >= 0; --j){
		q = false;
		for (int i = 0; i < sz1; ++i)
			if (a[i][j] != 0){
				q = true;
				x.j2 = j;
				break;
			}
		if (q == true)
			break;
	}
	

	return x;

}

//template<typename T_array>
void output_functions::ArrayToCSV(unsigned int ** z_lat, bool ** visited, int i1, int i2, int j1, int j2, const char* filename)
{
	// given the box [i1, i2]x[j1,j2] in the array T_array, saves the box into a csv file

	ofstream csv_file;
	string row_data = "";

	csv_file.open(filename);
	
	for (int i = i1; i <= i2; ++i){
		row_data = "";
		for (int j = j1; j <= j2 - 1; ++j)
			if (visited[i][j] == true)
				row_data += std::to_string(static_cast <unsigned long long>(z_lat[i][j])) + ",";
			else
				row_data += "-1,";

		if (visited[i][j2] == true)
			row_data += std::to_string(static_cast <unsigned long long>(z_lat[i][j2])) + "\n";
		else
			row_data +=  "-1\n";

		csv_file << row_data;
	}

	csv_file.close();
}

template<typename T_array>
void output_functions::ArrayToCSV(T_array ** a, int sz1, int sz2, const char* filename)
{
	// outputs the array a into csv file with the given filename

	ofstream csv_file;
	string row_data = "";

	csv_file.open(filename);

	for (int i = 0; i < sz1; ++i){
		row_data = "";
		for (int j = 0; j < sz2 - 1; ++j)
			row_data = row_data + std::to_string(static_cast <unsigned long long>(a[i][j])) + ",";

		row_data = row_data + std::to_string(static_cast <unsigned long long>(a[i][sz2 - 1])) + "\n";
		
		csv_file << row_data;
	}

	csv_file.close();
}

void output_functions::ArrayToPPM(unsigned int ** z_lat, bool ** visited, int i1, int i2, int j1, int j2, const char * filename)
{
	// output the subarray [i1, i2]x[j1,j2] into a ppm (raw image) file 

	ofstream ppm_file;
	string row_data = "";

	ppm_file.open(filename);

	ppm_file << (string("P3\n") + std::to_string(i2 - i1 + 1) + " " + std::to_string(j2 - j1 + 1) + "\n255\n" );

	for (int i = i1; i <= i2; ++i){
		row_data = "";
		for (int j = j1; j <= j2; ++j){
			if (z_lat[i][j] == 1)
				row_data += "255 128 255 ";
			else if (z_lat[i][j] == 2)
				row_data += "255 0 0 ";
			else if (z_lat[i][j] == 3)	
				row_data += "0 128 255 ";
			else
				row_data += "0 0 0 ";				
		}
		ppm_file << (row_data + "\n");
	}
	ppm_file.close();
}