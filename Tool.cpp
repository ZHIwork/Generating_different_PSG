#include "Tool.h"



Tool::Tool()
{
	//std::cout << "Tool constructor" << std::endl;
}


Tool::~Tool()
{
	//std::cout << "Tool deconstructor" << std::endl;
}


//convert coordinate from Direct to Cartesian
void Tool::directtoCart(const T(&direct)[3], T(&cart)[3])
{
	int m = 1, n = 3, k = 3;
	double alpha = 1, beta = 0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, alpha, direct, k, latticeVector, n, beta, cart, n);
}




//convert coordinate from Cartesian to Direct
void Tool::carttoDirect(const T(&cart)[3], T(&direct)[3])
{
	int m = 1, n = 3, k = 3;
	double alpha = 1, beta = 0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, alpha, cart, k, inverseLatticeVector, n, beta, direct, n);
}


//get distance between two sets of Cartesian coordinates  
T Tool::distanceinCart(const T(&a)[3], const T(&b)[3])
{
	return sqrt((a[0] - b[0])*(a[0] - b[0]) + (a[1] - b[1])*(a[1] - b[1]) + (a[2] - b[2])*(a[2] - b[2]));
}

//get angle between one center atom and two bonded atoms using Cartesian coordinates
T Tool::angleinCart(const T(&center)[3], const T(&bonded1)[3], const T(&bonded2)[3], const T &c1, const T &c2)
{
	//cosine
	T costheta =((center[0] - bonded1[0])*(center[0] - bonded2[0]) + (center[1] - bonded1[1])*(center[1] - bonded2[1]) + (center[2] - bonded1[2])*(center[2] - bonded2[2])) / (c1*c2);
	return costheta;
}


//get real direct shortest coordinates between center atom and bonded atom 
void Tool::directShortestCoor(const T (&center)[3], const T (&comp)[3], T (&shortest)[3])
{
	//for x
	T xcomp = fabs(center[0] - comp[0]);//x in comp not change
	shortest[0] = comp[0];
	T abs_xOsubs1 = fabs(center[0] - (comp[0] - 1));
	T abs_xOadds1 = fabs(center[0] - (comp[0] + 1));
	if (abs_xOsubs1 < xcomp)//x in comp -1
	{
		shortest[0] = comp[0] - 1;
		xcomp = abs_xOsubs1;
		if (abs_xOadds1 < xcomp)//x in comp +1
			shortest[0] = comp[0] + 1;

	}
	else
	{
		if (abs_xOadds1 < xcomp)//x in comp +1
			shortest[0] = comp[0] + 1;
	}

	//for y
	T ycomp = fabs(center[1] - comp[1]);//y in comp not change
	shortest[1] = comp[1];
	T abs_yOsubs1 = fabs(center[1] - (comp[1] - 1));
	T abs_yOadds1 = fabs(center[1] - (comp[1] + 1));
	if (abs_yOsubs1 < ycomp)//y in comp -1
	{
		shortest[1] = comp[1] - 1;
		ycomp = abs_yOsubs1;
		if (abs_yOadds1 < ycomp)//y in comp +1
			shortest[1] = comp[1] + 1;

	}
	else
	{
		if (abs_yOadds1 < ycomp)//y in comp +1
			shortest[1] = comp[1] + 1;
	}

	//for z
	T zcomp = fabs(center[2] - comp[2]);//z in comp not change
	shortest[2] = comp[2];
	T abs_zOsubs1 = fabs(center[2] - (comp[2] - 1));
	T abs_zOadds1 = fabs(center[2] - (comp[2] + 1));
	if (abs_zOsubs1 < zcomp)//z in comp -1
	{
		shortest[2] = comp[2] - 1;
		zcomp = abs_zOsubs1;
		if (abs_zOadds1 < zcomp)//z in comp +1
			shortest[2] = comp[2] + 1;

	}
	else
	{
		if (abs_zOadds1 < zcomp)//z in comp +1
			shortest[2] = comp[2] + 1;
	}
};


bool Tool::MC_probability(const T &old_energy, const T &new_energy)
{
	//difference of energy
	T deltE = new_energy - old_energy;
	//probability
	T comp = exp(-deltE / (temperature*BoltzmannConstant));
	if (comp < 1.0)
	{
		//get another random number
		long double randomNum = generator_lb(0.0, 1.0);
		if (randomNum < comp)
			return 1;
		else
			return 0;
	}
	else
		return 1;
}


T Tool::generator_lb(T begin, T end)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(begin, end);
	return dis(gen);
}


void Tool::copyGraph(const int(&a)[][129], int(&b)[][129], const int &row)
{
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < 128; ++j)
			b[i][j] = a[i][j];
	}
}


void Tool::copyCoordinates(const T a[], T b[], const int &N)
{
	T *temp = new T[N];
	for (int i = 0; i < N; ++i)
	{
		temp[i] = a[i];
		b[i] = temp[i];
	}
	delete[] temp;
}
