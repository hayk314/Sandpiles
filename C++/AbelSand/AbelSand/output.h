//using namespace std;

namespace outputFunctions
{
	struct Box_coord
	{
		int i1, i2, j1, j2;
	};

	Box_coord TrimmedArray(bool** A, int sz1, int sz2);

	void Array_toCSV(unsigned int** Z_lat, bool** Visited, int i1, int i2, int j1, int j2, const char* fileName);
	
	template <typename T_array>
	void Array_toCSV(T_array** A, int sz1, int sz2, const char* fileName);


	void Array_toPPM(unsigned int ** Z_lat, bool ** Visited, int i1, int i2, int j1, int j2, const char* fileName);

}
