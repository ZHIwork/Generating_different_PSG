#ifndef CONSTANT_H_
#define CONSTANT_H_
#include <string>
#include <cmath>

//define type
typedef double T;

//IO interface
extern const std::string importFileName;
extern const std::string exportFileName;


//simulation parameters
extern const int converganceTimes;
extern const int cutBondsTimes;


//Environment parameters
extern const T temperature;


//Constant that uses in equations
extern const T pi;

extern const T BoltzmannConstant;

extern const T massSi;
extern const T massP;
extern const T massO;

extern const T b0Si;
extern const T b0P;
extern const T b0PDoubleBond;
extern const T kbSi;
extern const T kbP;
extern const T kbPDoubleBond;
extern const T kthetaSi;
extern const T kthetaO;
extern const T kthetaP;
extern const T kthetaPDoubleBond;
extern const T thetaSi;
extern const T thetaO;
extern const T thetaP;
extern const T thetaPDoubleBond;
extern const T cos0_Si;
extern const T cos0_O;
extern const T cos0_P;
extern const T cos0_PDoubleBond;
extern const T repulsiveCoeforSiSiorSiPorPP;
extern const T repulsiveCoeforOO;



//POSCAR CONSTANT
extern T latticeVector[];
extern T inverseLatticeVector[];
extern T latticeConstant;
extern const int atomType;
extern std::string atomName[];
extern const int atomQuantity[];
extern const int totalAtomQuantity;
extern T *atomCoordinate;
extern void freeAtomCoordinate();
extern T *atomOldCoordinate;
extern void freeAtomOldCoordinate();
extern T *atomAcceleration;
extern void freeAtomAcceleration();
extern int startPoint[];


//For Group storing
extern int graph[][129];//number 129 need to be changed when changed POSCAR(also in Group.cpp need to be revised)
extern int oldGraph[][129];
extern int groupSi[][4];
extern int groupO[][2];
extern int groupP[][5];


//
extern T totalEnergy;
extern const T frequency;
extern T oldEnergy;
extern T infiniteNumber;

#endif /* CONSTANT_H_ */





