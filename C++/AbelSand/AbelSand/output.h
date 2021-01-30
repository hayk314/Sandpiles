//using namespace std;

namespace output_functions
{
	struct BoxCoord
	{
		int i1, i2, j1, j2;
	};

	BoxCoord TrimmedArray(bool** a, int sz1, int sz2);

	void ArrayToCSV(unsigned int** z_lat, bool** visited, int i1, int i2, int j1, int j2, const char* filename);
	
	template <typename T_array>
	void ArrayToCSV(T_array** a, int sz1, int sz2, const char* filename);
	
	void ArrayToPPM(unsigned int ** z_lat, bool ** visited, int i1, int i2, int j1, int j2, const char* filename);
}
