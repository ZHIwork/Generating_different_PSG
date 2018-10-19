#include "BFGS.h"

EnergyForce BFGS_EF;
Tool BFGS_Tool;

BFGS::BFGS()
{
	//std::cout << "call BFGS constructor" << std::endl;
	//initialization
	dimension =  193 * 3;
	//initialize inverse of B0 
	inverseB0 = new T[dimension*dimension];
	setZero(inverseB0, dimension*dimension);
	for (int i = 0; i < dimension; i++)
		inverseB0[i*dimension + i] = 1.0e-15 * 1.0e-15;

	//initialize x  
	x_before = new T[dimension];
	setZero(x_before, dimension);
	x_after = new T[dimension];
	setZero(x_after, dimension);

	//initialize inverseB
	inverseB_before = new T[dimension*dimension];
	setZero(inverseB_before, dimension*dimension);
	inverseB_after = new T[dimension*dimension];
	setZero(inverseB_after, dimension*dimension);

	//initialize g
	g_before = new T[dimension];
	setZero(g_before, dimension);
	g_after = new T[dimension];
	setZero(g_after, dimension);

	//initialize S
	S_before = new T[dimension];
	setZero(S_before, dimension);

	//initialize theta and y
	theta = new T[dimension];
	setZero(theta, dimension);
	y = new T[dimension];
	setZero(y, dimension);

	//initialize other parameters 
	Alpha = 1;
	Beta = 0.1;//between 0 to 1; Should be quit small;
	Gamma = 0.9;//between beta to 1; Should be quit large;
	t = 0.5;//between 0 to 1;

	//initialize temporary array for calculation only
	temp2Dimension = new T[dimension*dimension];
	setZero(temp2Dimension, dimension*dimension);
	tempDimension = new T[dimension];
	setZero(tempDimension, dimension);
	tempOne = new T[1];
	setZero(tempOne, 1);
	I = new T[dimension*dimension];
	setZero(I, dimension*dimension);
	for (int i = 0; i < dimension; i++)
		I[i*dimension + i] = 1.0;

}


BFGS::~BFGS()
{
	//free memory of inverseB0
	free1Darray(inverseB0);

	//free memory of x
	free1Darray(x_before);
	free1Darray(x_after);

	//free memory of inverseB
	free1Darray(inverseB_before);
	free1Darray(inverseB_after);

	//free memory of g
	free1Darray(g_before);
	free1Darray(g_after);

	//free memory of S
	free1Darray(S_before);

	//free memory of theta and y
	free1Darray(theta);
	free1Darray(y);

	//free temperory array
	free1Darray(temp2Dimension);
	free1Darray(tempDimension);
	free1Darray(I);
	delete tempOne;
}



//Line Search
void BFGS::LineSearch(T &totalEnergy, const int doubleO[2][2])
{

	//set x_before and inverseB_before
	convertDirecttoCartesianCoordinates(atomCoordinate, x_before, dimension);//assign x_before(Cartesian) with all atoms coordinates(Direct)
	assignValuefor1Darray(inverseB0, inverseB_before, dimension*dimension);

	/*initialize inverse Hessian approximation*/
	//get acceleration g_before
	BFGS_EF.obtainEnergyForce(atomCoordinate, g_before, totalEnergy, doubleO);

	//get S_before = -B0^-1 * g0
	int m = 193 * 3 , n = 1, k = 193 * 3 ;
	T a = -1, b = 0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, a, inverseB_before, k, g_before, n, b, S_before, n);

	//calculate x_after = x_before + alpha * S_before.
	numMultiply(Alpha, S_before, tempDimension, dimension);
	add(x_before, tempDimension, x_after, dimension);
	//get g_after according to x_after
	convertandLimitCarttoDirectCoordinates(x_after, atomCoordinate, dimension);//convert Cartesian to Direct and limit value between 0 to 1
	BFGS_EF.obtainEnergyForce(atomCoordinate, g_after, totalEnergy, doubleO);

	//std::cout << "First cut energy: "<<totalEnergy <<std::endl;

	//calculate theta = x_after - x_before
	subtract(x_after, x_before, theta, dimension);
	
	//calculate y = g_after - g_before
	subtract(g_after, g_before, y, dimension);

	//get Inverse Hessian approximation 
	numMultiply(specialMultiply(y, theta, dimension) / specialMultiply(y, y, dimension), I, inverseB_before, dimension*dimension);

	T BFGSOldEnergy = totalEnergy;

	//start iterations
	int num = 0;
	while (num < 100) 
	{
		//get S_before = -inverseB_before * g_before
		m = 193 * 3;
		n = 1;
		k = 193 * 3;
		a = -1;
		b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, inverseB_before, k, g_before, n, b, S_before, n);

		//get value of alpha
		Alpha = GetAlpha(Alpha, S_before, g_before, totalEnergy, doubleO);

		//calculate x_after = x_before + alpha * S_before.
		numMultiply(Alpha, S_before, tempDimension, dimension);
		add(x_before, tempDimension, x_after, dimension);
		//get g_after according to x_after
		convertandLimitCarttoDirectCoordinates(x_after, atomCoordinate, dimension);//convert Cartesian to Direct and limit value between 0 to 1
		BFGS_EF.obtainEnergyForce(atomCoordinate, g_after, totalEnergy, doubleO);
		

		std::cout << num << "  energy:" << totalEnergy << std::endl;

		//set limit and stop
		/*if ((BFGSOldEnergy - totalEnergy) < 1e-5)
			break;
		BFGSOldEnergy = totalEnergy;*/


		//calculate theta = x_after - x_before
		subtract(x_after, x_before, theta, dimension);

		//calculate y = g_after - g_before
		subtract(g_after, g_before, y, dimension);


		//Best BFGS update
		//Common Factor 1 / [(theta)^T * y]
		m = 1;
		n = 1;
		k = 193 * 3;
		a = 1;
		b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, theta, k, y, n, b, tempOne, n);
		T CommonFactor = 1.0 / tempOne[0];

		//firstItem1 = 1.0 + {[(y)^T * inverseB_before] * y} * CommonFactor 
		m = 1;
		n = 193 * 3;
		k = 193 * 3;
		a = 1;
		b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, y, k, inverseB_before, n, b, tempDimension, n);//[(y)^T * inverseB_before]

		m = 1;
		n = 1;
		k = 193 * 3;
		a = 1;
		b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, tempDimension, k, y, n, b, tempOne, n);//{[(y)^T * inverseB_before] * y}
		T firstItem1 = 1.0 + (tempOne[0]) * CommonFactor;//1.0 + {[(y)^T * inverseB_before] * y} * CommonFactor 

		//firstItem = [theta * (theta)^T] * CommonFactor * firstItem1 
		m = 193 * 3, n = 193 * 3, k = 1;
		a = 1, b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, theta, k, theta, n, b, temp2Dimension, n);//[theta * (theta)^T] * CommonFactor
		T* firstItem;
		firstItem = new	T[193 * 3 * 193 * 3];//define firstItem in 2D dimension
		numMultiply(CommonFactor*firstItem1, temp2Dimension, firstItem, dimension*dimension);//[theta * (theta)^T] * CommonFactor * firstItem1 
		//secondItem1 =[theta * (y)^T] * inverseB_before
		T *secondItem1;
		secondItem1 = new T[193 * 3 * 193 * 3];//define secondItem1 in 2D dimension
		m = 193 * 3, n = 193 * 3, k = 1;
		a = 1, b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, theta, k, y, n, b, temp2Dimension, n);
		m = 193 * 3, n =193 * 3, k = 193 * 3;
		a = 1, b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, temp2Dimension, k, inverseB_before, n, b, secondItem1, n);

		

		//secondItem2 = (inverseB_before * y) * (theta)^T
		T* secondItem2;
		secondItem2 = new T[193 * 3 * 193 * 3];//define secondItem2 in 2D dimension
		m = 193 * 3, n = 1, k = 193 * 3;
		a = 1, b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, inverseB_before, k, y, n, b, tempDimension, n);//(inverseB_before * y)
		m = 193 * 3, n = 193 * 3, k = 1;
		a = 1, b = 0;
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			m, n, k, a, tempDimension, k, theta, n, b, secondItem2, n);//(inverseB_before * y) * (theta)^T
		//secondItem = (secondItem1 + secondItem2) * CommonFactor
		T* secondItem;
		secondItem = new T[193 * 3 * 193 * 3];//define secondItem in 2D dimension
		add(secondItem1, secondItem2, temp2Dimension, dimension*dimension);//(secondItem1 + secondItem2)
		numMultiply(CommonFactor, temp2Dimension, secondItem, dimension*dimension);//(secondItem1 + secondItem2) * CommonFactor

		//get inverseB_after = inverseB_before + firstItem - secondItem
		add(inverseB_before, firstItem, temp2Dimension, dimension*dimension);
		subtract(temp2Dimension, secondItem, inverseB_after, dimension*dimension);


		//prepare new cycle parameters
		assignValuefor1Darray(x_after, x_before, dimension);//assign x_after to x_before
		assignValuefor1Darray(inverseB_after, inverseB_before, dimension*dimension);//assign inverseB_after to inverseB_before
		assignValuefor1Darray(g_after, g_before, dimension);//assign g_after to g_before
		num++;
		

		//free memory
		delete[] firstItem;
		delete[] secondItem1;
		delete[] secondItem2;
		delete[] secondItem;
	}

	//assign x_before to atomCoordinate
	convertandLimitCarttoDirectCoordinates(x_before, atomCoordinate, dimension);
};


//update Alpha in each iteration using Backtracking Line Search
T BFGS::GetAlpha(const T &alpha, const T* Sk, const T* gk, const T &totalEnergy, const int doubleO[2][2])
{
	
	T newAlpha = alpha;//set old alpha as the new alpha if it does not satisfy the condition
	T EnergyX = 0;//define a new total energy

	//calculate x_before + alpha * S_before (left hand)
	numMultiply(alpha, Sk, tempDimension, dimension);
	add(tempDimension, x_before, XCart, dimension);//assign the result to new coordiantes arrray X

	//calculate  alpha * beta * (gk)^T * S_before  into temp (right hand)
	T F_righthand_secondterm = 0;
	F_righthand_secondterm = specialMultiply(gk, Sk, dimension) * alpha * Beta;
	
	//Curvature condition items
	convertandLimitCarttoDirectCoordinates(XCart, XDirect, dimension);
	BFGS_EF.obtainEnergyForce(XDirect, gX, EnergyX, doubleO);
	//left hand
	T Cur_lefthand[1];
	int m = 1, n = 1, k = 193 * 3;
	double a = 1, b = 0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, a, gX, k, Sk, n, b, Cur_lefthand, n);

	//right hand
	T Cur_righthand[1];
	m = 1, n = 1, k = 193 * 3;
	a = Gamma, b = 0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		m, n, k, a, gk, k, Sk, n, b, Cur_righthand, n);

	//determine equation
	if (EnergyX <= (totalEnergy + F_righthand_secondterm) &&
		Cur_lefthand > Cur_righthand)//Armijo condition and Curvature condition
	{
		//update new alpha
		newAlpha = t * alpha;
	}
	return newAlpha;
};


//initialize a 1D array to 0
void BFGS::setZero(T* a, const int &N)
{
	for (int i = 0; i < N; i++)
		a[i] = 0.0;
}



//assign value for a 1D array by another 1D array
void BFGS::assignValuefor1Darray(const T* a, T* b, const int &N)
{
	T *temp = new T[N];
	for (int i = 0; i < N; i++)
	{
		temp[i] = a[i];
		b[i] = temp[i];
	}
	delete[] temp;
}


//add two array
void BFGS::add(const T* a, const T* b, T c[], const int &N)
{
	for (int i = 0; i < N; i++)
		c[i] = a[i] + b[i];
}


//substract two array
void BFGS::subtract(const T* a, const T* b, T* c, const int &N)
{
	for (int i = 0; i < N; i++)
		c[i] = a[i] - b[i];
}


//number multiplication
void BFGS::numMultiply(const T &num, const T* a, T* b, const int &N)
{
	for (int i = 0; i < N; ++i)
		b[i] = num * a[i];
}


//special multiplication (a)T * b = a_number instead of matrix
T BFGS::specialMultiply(const T* a, const T* b, const int &N)
{
	T sum = 0;
	for (int i = 0; i < N; i++)
		sum += a[i] * b[i];
	return sum;
}


//convert atom coordinates from Direct to Cartesian
void BFGS::convertDirecttoCartesianCoordinates(const T* coorDirect, T* coorCartesian, const int &N)
{
	T tempDirect[3], tempCartesian[3];
	for (int i = 0; i < N; i=i+3)
	{
		//get direct coordinates
		tempDirect[0] = coorDirect[i];//x
		tempDirect[1] = coorDirect[i + 1];//y
		tempDirect[2] = coorDirect[i + 2];//z
		
		//convert direct to cartesian
		BFGS_Tool.directtoCart(tempDirect, tempCartesian);

		//assign tempCartesian value to real coorCartesian
		coorCartesian[i] = tempCartesian[0];//x
		coorCartesian[i+1] = tempCartesian[1];//y
		coorCartesian[i+2] = tempCartesian[2];//z
		
	}
}


//convert atom coordiantes from Cartesian to Direct and limit it between 0 to 1
void BFGS::convertandLimitCarttoDirectCoordinates(const T* coorCartesian, T* coorDirect, const int &N)
{
	T tempDirect[3], tempCartesian[3];
	for (int i = 0; i < N; i = i + 3)
	{
		//get Cartesian coordinates
		tempCartesian[0] = coorCartesian[i];//x
		tempCartesian[1] = coorCartesian[i + 1];//y
		tempCartesian[2] = coorCartesian[i + 2];//z

		//convert direct to cartesian
		BFGS_Tool.carttoDirect(tempCartesian, tempDirect);

		//limit x, y, z between 0 to 1
		for (int j = 0; j < 3; ++j)
		{
			if ((int)tempDirect[j] >= 1.0)
				tempDirect[j] -= (int)tempDirect[j];
			else if ((int)tempDirect[j] < 0.0)
				tempDirect[j] -= ((int)tempDirect[j] - 1.0);
		}

		//assign tempCartesian value to real coorCartesian
		coorDirect[i] = tempDirect[0];//x
		coorDirect[i + 1] = tempDirect[1];//y
		coorDirect[i + 2] = tempDirect[2];//x
	}
}


void BFGS::limitCoor(T a[], const int &N)
{
	for (int i = 0; i < N; ++i)
	{
	
		while (a[i] >= 1.0) {
			std::cout << i << "     " << a[i] << std::endl;
			a[i] -= 1.0;
			
		}
		while (a[i] < 0.0) {
			std::cout << i << "     " << a[i] << std::endl; 
			a[i] += 1.0;
			
		}
	}
}


//free mempry of 1D array
void BFGS::free1Darray(T *array)
{
	delete[] array;
}