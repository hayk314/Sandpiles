#include <iostream>
#include <fstream>
#include <string>

#include "output.h"

using namespace outputFunctions;
using namespace std;

//template<typename T_array>
Box_coord outputFunctions::TrimmedArray(bool ** A, int sz1, int sz2)
{
    // trim any zero row/columns from the borders of A and return the resulting rectangle

	outputFunctions::Box_coord X;
	X.i1 = -1; X.i2 = -1; X.j1 = -1; X.j2 = -1;

	bool q = false;

	for (int i = 0; i < sz1; i++)
	{
		q = false;
		for (int j = 0; j < sz2; j++)
			if (A[i][j] != 0)
			{
				q = true;
				X.i1 = i;
				break;
			}
		if (q == true)
			break;
	}

	// i2
	for (int i = sz1 - 1; i >= 0; i--)
	{
		q = false;
		for (int j = 0; j < sz2; j++)
			if (A[i][j] != 0)
			{
				q = true;
				X.i2 = i;
				break;
			}
		if (q == true)
			break;
	}


	// j1
	for (int j = 0; j < sz2; j++)
	{
		q = false;
		for (int i = 0; i < sz1; i++)
			if (A[i][j] != 0)
			{
				q = true;
				X.j1 = j;
				break;
			}
		if (q == true)
			break;
	}


	// j2
	for (int j = sz2 - 1; j >= 0; j--)
	{
		q = false;
		for (int i = 0; i < sz1; i++)
			if (A[i][j] != 0)
			{
				q = true;
				X.j2 = j;
				break;
			}
		if (q == true)
			break;
	}



	return X;

}

//template<typename T_array>
void outputFunctions::Array_toCSV(unsigned int ** Z_lat, bool ** Visited, int i1, int i2, int j1, int j2, const char* fileName)
{
	// given the box [i1, i2]x[j1,j2] in the array T_array, saves the box into a csv file

	ofstream csv_File;
	string rowData = "";

	csv_File.open(fileName);


	for (int i = i1; i <= i2; i++)
	{
		rowData = "";
		for (int j = j1; j <= j2 - 1; j++)
			if (Visited[i][j] == true)
				rowData += std::to_string(static_cast <unsigned long long>(Z_lat[i][j])) + ",";
			else
				rowData += "-1,";

		if (Visited[i][j2] == true)
			rowData += std::to_string(static_cast <unsigned long long>(Z_lat[i][j2])) + "\n";
		else
			rowData +=  "-1\n";

		csv_File << rowData;
	}

	csv_File.close();
}

template<typename T_array>
void outputFunctions::Array_toCSV(T_array ** A, int sz1, int sz2, const char* fileName)
{
	// outputs the array A into csv file with the given fileName

	ofstream csv_File;
	string rowData = "";

	csv_File.open(_Fname);

	for (int i = 0; i < sz1; i++)
	{
		rowData = "";
		for (int j = 0; j < sz2 - 1; j++)
			rowData = rowData + std::to_string(static_cast <unsigned long long>(A[i][j])) + ",";

		rowData = rowData + std::to_string(static_cast <unsigned long long>(A[i][sz2 - 1])) + "\n";


		csv_File << rowData;
	}

	csv_File.close();
}

void outputFunctions::Array_toPPM(unsigned int ** Z_lat, bool ** Visited, int i1, int i2, int j1, int j2, const char * fileName)
{
	// output the subarray [i1, i2]x[j1,j2] into a ppm (raw image) file 

	ofstream ppm_File;
	string rowData = "";

	ppm_File.open(fileName);

	ppm_File << (string("P3\n") + std::to_string(i2 - i1 + 1) + " " + std::to_string(j2 - j1 + 1) + "\n255\n" );

	for (int i = i1; i <= i2; i++)
	{
		rowData = "";
		for (int j = j1; j <= j2; j++)
		{
			if (Z_lat[i][j] == 1)
				rowData += "255 128 255 ";
			else if (Z_lat[i][j] == 2)
				rowData += "255 0 0 ";
			else if (Z_lat[i][j] == 3)	
				rowData += "0 128 255 ";
			else
				rowData += "0 0 0 ";				
		}
		ppm_File << (rowData + "\n");
	}
	ppm_File.close();

}
