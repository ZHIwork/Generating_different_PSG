#include "CONSTANT.h"

//IO interface
const std::string importFileName = "PSG.vasp";//input file name.
const std::string exportFileName = "final.vasp";//output file name.


//simulation parameters
const int converganceTimes = 100;//when minimize energy, the number of times to converge.
const int cutBondsTimes = 3;//when simulate amorphous system, the number of times to cut bonds.


//Environment parameters
const T temperature = 6000.0;//temperature that used to form amorphous system.


//Constant that uses in equations
const T pi = 3.141592653589793;//pi constant

const T BoltzmannConstant = 8.617330350e-5;//Boltzmann Constant

const T massSi = 46.82e-27;//mass of Si
const T massP = 51.43e-27;//mass of P
const T massO = 26.67e-27;//mass of O

//bond length
const T b0Si = 1.55148;//Si-O
const T b0P = 1.578;//P-O
const T b0PDoubleBond = 1.458;//P=O

//two body force constant
const T kbSi = 26.96;
const T kbP = 28.00;//(guess value)
const T kbPDoubleBond = 35.00;//(guess value)

//three body force constant
const T kthetaSi = 1.685;
const T kthetaO = 0.58;
const T kthetaP = 1.80;//(guess value)
const T kthetaPDoubleBond = 2.00;//(guess value)

//angle
const T thetaSi = 109.47;//O-Si-O
const T thetaO = 120.0;//P-O-Si or Si-O-Si
const T thetaP = 103.0;//O-P-O
const T thetaPDoubleBond = 114.7;//O=P-O
const T cos0_Si = cos(thetaSi * pi / 180.0);
const T cos0_O = cos(thetaO * pi / 180.0);
const T cos0_P = cos(thetaP * pi / 180.0);
const T cos0_PDoubleBond = cos(thetaPDoubleBond * pi / 180.0);

const T repulsiveCoeforSiSiorSiPorPP = 8.0;
const T repulsiveCoeforOO = 6.0;



//POSCAR CONSTANT
T latticeConstant = 0.0; //universal scaling factor (lattice constant),
						 //which is used to scale all lattice vectors and all atomic coordinates. 
T latticeVector[3 * 3]{ 0 }; //three lattice vectors defining the unit cell of the system are given

T inverseLatticeVector[3 * 3]{ 0.072025928,	0,	0,
0.002807564,	0.070502637,	0,
0.003485514,	0.000127587,	0.071120349
};//inverse lattice vector

const int atomType = 3; //the quantity of atom's types (need to be correct when use different POSCAR)
std::string atomName[atomType]; //all atom's names from POSCARS
const int atomQuantity[atomType]={129, 62, 2}; //the quantity of each type of atoms (O, Si and P orders)
const int totalAtomQuantity = 193; // the total number atoms (need to be correct when use different POSCAR)

T *atomCoordinate = new T[totalAtomQuantity * 3];//store all coordinates of x, y, z into this dynamic array. 
void freeAtomCoordinate() //free memory
{
	delete[] atomCoordinate;
}

T *atomOldCoordinate = new T[totalAtomQuantity * 3];//store all coordinates of x, y, z into this dynamic array. 
void freeAtomOldCoordinate() //free memory
{
	delete[] atomOldCoordinate;
}

T *atomAcceleration = new T[totalAtomQuantity * 3];//store all acceleration coordinates of x, y, z into this dynamic array. 
void freeAtomAcceleration() //free memory
{
	delete[] atomAcceleration;
}

int startPoint[atomType] = {0, 129*3, 191*3};//start point for different type atoms, need to be changed when use different POSCAR


//For Group storing
int graph[64][129]{ 0 };//first one is the number of Si and P, second is the number of O.
int oldGraph[64][129]{ 0 };//first one is the number of Si and P, second is the number of O.
int groupSi[62][4];//array to store all bonded O with each Si 
int groupO[129][2];//array to store all bonded Si with each O
int groupP[2][5];//array to store all bonded O with each P


//energy parameters
T totalEnergy = 0;//total energy of whole system
const T frequency = 1e-15;
T oldEnergy = 0;
T infiniteNumber = 1e8;
