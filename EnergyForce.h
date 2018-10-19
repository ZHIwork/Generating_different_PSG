#pragma once
#include "CONSTANT.h"
#include "Tool.h"

class EnergyForce
{
private:
	T angular_energy_Si, angular_energy_O, angular_energy_P;
	T radial_energySiandO, radial_energyPandO;
	T repulsive_energy_SiSiorSiPorPP, repulsive_energy_O_O;
	T directCenterCoorSi[3], directBondedCoorOwithSi[4][3], realDirectBondedCoorOwithSi[4][3];
	T directCenterCoorP[3], directBondedCoorOwithP[5][3], realDirectBondedCoorOwithP[5][3];
	T directCenterCoorO[3], directBondedCoorSiwithO[2][3], realDirectBondedCoorSiwithO[2][3];
	T cartesianCenterCoorSi[3], realCartesianDirectBondedCoorOwithSi[4][3];
	T cartesianCenterCoorP[3], realCartesianDirectBondedCoorOwithP[5][3];
	T cartesianCenterCoorO[3], realCartesianDirectBondedCoorSiwithO[2][3];
	int orderObondedwithSi[4], orderSibondedwithO[2];
	int orderObondedwithP[5], orderPbondedwithO[2];
	T Rf_O1_on_Si[3], Rf_O2_on_Si[3], Rf_O3_on_Si[3], Rf_O4_on_Si[3];
	T Rf_O1_on_P[3], Rf_O2_on_P[3], Rf_O3_on_P[3], Rf_O4_on_P[3], Rf_O5_on_P[3];
	T Rf_on_Si[3], Rf_on_P[3];
	T fSiO1O2_on_O1[3], fSiO1O2_on_O2[3], fSiO1O2_on_Si[3];
	T fSiO1O3_on_O1[3], fSiO1O3_on_O3[3], fSiO1O3_on_Si[3];
	T fSiO1O4_on_O1[3], fSiO1O4_on_O4[3], fSiO1O4_on_Si[3];
	T fSiO2O3_on_O2[3], fSiO2O3_on_O3[3], fSiO2O3_on_Si[3];
	T fSiO2O4_on_O2[3], fSiO2O4_on_O4[3], fSiO2O4_on_Si[3];
	T fSiO3O4_on_O3[3], fSiO3O4_on_O4[3], fSiO3O4_on_Si[3];

	T fPO1O2_on_O1[3], fPO1O2_on_O2[3], fPO1O2_on_P[3];
	T fPO1O3_on_O1[3], fPO1O3_on_O3[3], fPO1O3_on_P[3];
	T fPO1O4_on_O1[3], fPO1O4_on_O4[3], fPO1O4_on_P[3];
	T fPO1O5_on_O1[3], fPO1O5_on_O5[3], fPO1O5_on_P[3];
	T fPO2O3_on_O2[3], fPO2O3_on_O3[3], fPO2O3_on_P[3];
	T fPO2O4_on_O2[3], fPO2O4_on_O4[3], fPO2O4_on_P[3];
	T fPO2O5_on_O2[3], fPO2O5_on_O5[3], fPO2O5_on_P[3];
	T fPO3O4_on_O3[3], fPO3O4_on_O4[3], fPO3O4_on_P[3];
	T fPO3O5_on_O3[3], fPO3O5_on_O5[3], fPO3O5_on_P[3];
	T fPO4O5_on_O4[3], fPO4O5_on_O5[3], fPO4O5_on_P[3];

	T fOSi1Si2_on_Si1[3], fOSi1Si2_on_Si2[3], fOSi1Si2_on_O[3];

	T directRepSiSi[3], directRepSiSishortest[3], cartesianRepSiSishortest[3], repSiSilength, frepSiSi[3];
	T directRepSiP[3], directRepSiPshortest[3], cartesianRepSiPshortest[3], repSiPlength, frepSiP[3];
	T directRepPP[3], directRepPPshortest[3], cartesianRepPPshortest[3], repPPlength, frepPP[3];
	T directRepOO[3], directRepOOshortest[3], cartesianRepOOshortest[3], repOOlength, frepOO[3];
public:
	EnergyForce();
	~EnergyForce();
	void obtainEnergyForce(const T* coordinates, T* acceleration, T& total_energy, const int doubleO[2][2]);
};

