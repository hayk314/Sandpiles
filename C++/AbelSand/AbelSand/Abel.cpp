#include <iostream>
#include <vector>
#include <time.h>
#include <string>

#include "Abel.h"
#include "output.h"

using namespace std;
using namespace Abel;
using namespace outputFunctions;

double Abel::moveStandard( unsigned int N)
{
		// sandpile stabilization with N particles at the origin in a standard way, 
	    // with preallocated stack of moving vertices

		unsigned int** Z_lat = new unsigned int*[N_size];      // models the standard lattice Z^2
		bool** V_sites = new bool*[N_size];                    // vertices of Z^2 which were visited during the process
		bool** Move_ = new bool*[N_size];                      // vertices of Z^2 which are already in the walking stack
		unsigned int** Odometer = new unsigned int*[N_size];   // total number of topplings of a given vertex of Z^2


		for (int k = 0; k < N_size; k++)
		{
			Z_lat[k] = new unsigned int[N_size];
			V_sites[k] = new bool[N_size];
			Odometer[k] = new unsigned int[N_size];
			Move_[k] = new bool[N_size];

			for (int i = 0; i<N_size; i++)
			{
				Z_lat[k][i] = 0;  Move_[k][i] = false; V_sites[k][i] = false;
			}
		}

		L_coord* walking_ = new L_coord[N_size*N_size]; //pre-allocated stack of vertices which need to be toppled
		
		V_sites[N_size2][N_size2] = true;
		Move_[N_size2][N_size2] = true;

		Z_lat[N_size2][N_size2] = N;
		walking_[0].x = N_size2;
		walking_[0].y = N_size2;


		unsigned int x = 0, y = 0;
		int top_ = 0;

		unsigned int lx = 0, ly = 0;
		unsigned int d = 0;

		clock_t t1, t2;
		t1 = clock();

		unsigned int max_of_top = 0;
		unsigned long N_of_Moves = 0;


		while (top_ >= 0)
		{
			N_of_Moves += 1;

			if (max_of_top<top_) { max_of_top = top_; }

			x = walking_[top_].x;      y = walking_[top_].y;
			
			top_--; Move_[x][y] = false;
			
			d = (Z_lat[x][y] >> 2);
			Z_lat[x][y] = Z_lat[x][y] - (d << 2);

			for (int k = 0; k < 4; k++)
			{
				lx = x + dx[k];
				ly = y + dy[k];

				V_sites[lx][ly] = true;
				Z_lat[lx][ly] += d;

				if (Move_[lx][ly] == false && Z_lat[lx][ly] >= 4)
				{
					walking_[++top_].x = lx; 
					walking_[top_].y = ly; 
					Move_[lx][ly] = true;
				}
			}	

		}

		t2 = clock();

		outputFunctions::Box_coord B;
		B = TrimmedArray(V_sites, N_size, N_size);
		Array_toCSV(Z_lat, V_sites, B.i1, B.i2, B.j1, B.j2, ("Abel" + std::to_string(N) + ".csv").c_str());


		// clean-ups

		for (int k = 0; k<N_size; k++)
		{
			delete[] Z_lat[k];
			delete[] V_sites[k];
			delete[] Odometer[k];
			delete[] Move_[k];
		}

		delete[] Z_lat; delete[] V_sites; delete[] Odometer; delete[] Move_; delete[] walking_;


		cout << "Maximal number in the stack was " << max_of_top << endl;
		cout << "Number of moves in the main loop was " << N_of_Moves << endl;

		return ((float)(t2)-float(t1))*0.001;
	
}

double Abel::moveStandard_1Step(unsigned int N)
{

	unsigned int** Z_lat = new unsigned int*[N_size];
	bool** V_sites = new bool*[N_size];
	bool** Move = new bool*[N_size]; 
	unsigned int** Odometer = new unsigned int*[N_size];


	for (int k = 0; k<N_size; k++)
	{
		Z_lat[k] = new unsigned int[N_size];
		V_sites[k] = new bool[N_size];
		Odometer[k] = new unsigned int[N_size];
		Move[k] = new bool[N_size];

		for (int i = 0; i < N_size; i++)
		{
			Z_lat[k][i] = 0;  V_sites[k][i] = false; Odometer[k][i] = 0; Move[k][i] = false;
		}
	}

	L_coord* walking = new L_coord[N_size*N_size]; //pre-allocated stack

	V_sites[N_size2][N_size2] = true;

	int top = -1;
	unsigned int x = 0, y = 0, lx = 0, ly = 0;


	clock_t t1, t2;
	t1 = clock();

	unsigned int max_of_top = 0;
	unsigned long N_of_Moves = 0;

	for (unsigned int i = N ; i--; )
	{
		if (++Z_lat[N_size2][N_size2] >= 4)
		{
			walking[++top].x = N_size2; 
			walking[top].y = N_size2;
		}

		while (top >= 0)
		{
			N_of_Moves += 1;

			if (max_of_top < top) { max_of_top = top; }

			x = walking[top].x;      
			y = walking[top].y;

			Z_lat[x][y] = Z_lat[x][y] - 4;
			if (Z_lat[x][y] < 4)
			{
				top --; 
				Move[x][y] = false;
			}


			for (int k = 0; k < 4; k++)
			{
				lx = x + dx[k];
				ly = y + dy[k];
				
				V_sites[lx][ly] = true;
				Z_lat[lx][ly] ++;

				if (Move[lx][ly] == false && Z_lat[lx][ly] >= 4)
				{
					walking[++top].x = lx;
					walking[top].y = ly;
					Move[lx][ly] = true;
				}
			}
		}

	}




	// now running the process until stabilisation


	t2 = clock();

	cout << endl;

	cout << "Maximal number in the stack was " << max_of_top << endl;
	cout << "Number of moves in the main loop was " << N_of_Moves << endl;



	Box_coord B = TrimmedArray(V_sites, N_size, N_size);
	Array_toCSV (Z_lat, V_sites, B.i1, B.i2, B.j1, B.j2, ("Abel" +   std::to_string(N) + ".csv").c_str()  );
	Array_toPPM(Z_lat, V_sites, B.i1, B.i2, B.j1, B.j2, ("Abel" + std::to_string(N) + ".ppm").c_str());

	// clean-ups

	for (int k = 0; k<N_size; k++)
	{
		delete[] Z_lat[k];
		delete[] V_sites[k];
		delete[] Odometer[k];
		delete[] Move[k];
	}

	delete[] Z_lat; delete[] V_sites; delete[] Odometer; delete[] Move; delete[] walking;




	return ((double)(t2)-double(t1))*0.001;
}