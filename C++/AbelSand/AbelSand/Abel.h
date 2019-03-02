namespace Abel
{
	const int N_size = 1025;   //513;
	const int N_size2 = 512;   //256;  // (N-1)/2;

	const int dx[4] = { 1,0,-1,0 };
	const int dy[4] = { 0,1,0,-1 };

	struct L_coord
	{
		unsigned int x;
		unsigned int y;
	};

	double moveStandard(unsigned int N);        // performs Abelian sandpile toppling 1 by 1 
	double moveStandard_1Step(unsigned int N);  // performs toppling of Abelian sandpile by larger mass

}