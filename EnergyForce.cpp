#include "EnergyForce.h"


//call Tool.h
Tool EFTool;


EnergyForce::EnergyForce()
{
	

	
}


EnergyForce::~EnergyForce()
{
	//std::cout << "EnergyForce deconstructor" << std::endl;
}


void EnergyForce::obtainEnergyForce(const T* coordinates, T* acceleration, T& total_energy, const int doubleO[2][2])
{
	//initialize energy parameters
	angular_energy_Si = 0;
	angular_energy_O = 0;
	angular_energy_P = 0;
	radial_energySiandO = 0;
	radial_energyPandO = 0;
	repulsive_energy_SiSiorSiPorPP = 0;
	repulsive_energy_O_O = 0;

	//initialize acceleration dynamic array to 0 
	for (int i = 0; i < totalAtomQuantity * 3; ++i)
		acceleration[i] = 0.0;

	/* energy for Si */
	//calculate length between each two bounded atoms (Si and O)
	int a = 1;
	for (int i = 0; i < atomQuantity[1]; ++i) 
	{
		//store Si position using atomCoordiante(need to know start point of Si atom)
		directCenterCoorSi[0] = coordinates[i * 3 + startPoint[1]];
		directCenterCoorSi[1] = coordinates[i * 3 + 1 + startPoint[1]];
		directCenterCoorSi[2] = coordinates[i * 3 + 2 + startPoint[1]];

		//store bonded O position 
		for (int j = 0; j < 4; ++j)
		{
			orderObondedwithSi[j] = groupSi[i][j];//get bonded O orders

			//get coordinates of each bonded O
			directBondedCoorOwithSi[j][0] = coordinates[orderObondedwithSi[j] * 3 + startPoint[0]];//x
			directBondedCoorOwithSi[j][1] = coordinates[orderObondedwithSi[j] * 3 + 1 + startPoint[0]];//y
			directBondedCoorOwithSi[j][2] = coordinates[orderObondedwithSi[j] * 3 + 2 + startPoint[0]];//z
		}

		//get real direct shortest coordinates of bonded O with center Si
		for (int j = 0; j < 4; ++j)
			EFTool.directShortestCoor(directCenterCoorSi, directBondedCoorOwithSi[j], realDirectBondedCoorOwithSi[j]);


		//--------------------------------------------------------------------------------------------------------------

		//convert Direct to Cartesian coordinates
		EFTool.directtoCart(directCenterCoorSi, cartesianCenterCoorSi);//convert center Si coordiantes from Direct to Cartesian
		for (int j = 0; j < 4; ++j)
			EFTool.directtoCart(realDirectBondedCoorOwithSi[j], realCartesianDirectBondedCoorOwithSi[j]);


		/*std::cout << realDirectBondedCoorOwithSi[0][0] << "    "
			<< realDirectBondedCoorOwithSi[0][1] << "    "
			<< realDirectBondedCoorOwithSi[0][2] << std::endl;*/


		//length between center Si and bonded O
		T length_Si_O1 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0]);
		T length_Si_O2 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1]);
		T length_Si_O3 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[2]);
		T length_Si_O4 = EFTool.distanceinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[3]);


		//std::cout << length_Si_O1 << "    " << length_Si_O2 << "    " << length_Si_O3 << "    " << length_Si_O4 << std::endl;

		//angle
		T cosSi_O1_O2 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[1], length_Si_O1, length_Si_O2);
		T cosSi_O1_O3 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[2], length_Si_O1, length_Si_O3);
		T cosSi_O1_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[0], realCartesianDirectBondedCoorOwithSi[3], length_Si_O1, length_Si_O4);
		T cosSi_O2_O3 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1], realCartesianDirectBondedCoorOwithSi[2], length_Si_O2, length_Si_O3);
		T cosSi_O2_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[1], realCartesianDirectBondedCoorOwithSi[3], length_Si_O2, length_Si_O4);
		T cosSi_O3_O4 = EFTool.angleinCart(cartesianCenterCoorSi, realCartesianDirectBondedCoorOwithSi[2], realCartesianDirectBondedCoorOwithSi[3], length_Si_O3, length_Si_O4);




		//--------------------------------------------------------------------------------------------------------------

		/*calculate radial energy and radial force
		 * and store forces to  */

		 //coefficient
		T diff_Si_O1 = length_Si_O1 - b0Si;
		T diff_Si_O2 = length_Si_O2 - b0Si;
		T diff_Si_O3 = length_Si_O3 - b0Si;
		T diff_Si_O4 = length_Si_O4 - b0Si;

		//radial energy
		radial_energySiandO += 0.5 * kbSi *
			(diff_Si_O1*diff_Si_O1 + diff_Si_O2 * diff_Si_O2 +
				diff_Si_O3 * diff_Si_O3 + diff_Si_O4 * diff_Si_O4);

		//radial force for Si
		T fcoe = kbSi * diff_Si_O1 / length_Si_O1;
		Rf_O1_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][0] - cartesianCenterCoorSi[0]);//force on x by O1
		Rf_O1_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][1] - cartesianCenterCoorSi[1]);//force on y by O1
		Rf_O1_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[0][2] - cartesianCenterCoorSi[2]);//force on z by O1

		fcoe = kbSi * diff_Si_O2 / length_Si_O2;
		Rf_O2_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][0] - cartesianCenterCoorSi[0]);//force on x by O2
		Rf_O2_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][1] - cartesianCenterCoorSi[1]);//force on y by O2
		Rf_O2_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[1][2] - cartesianCenterCoorSi[2]);//force on z by O2

		fcoe = kbSi * diff_Si_O3 / length_Si_O3;
		Rf_O3_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][0] - cartesianCenterCoorSi[0]);//force on x by O3
		Rf_O3_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][1] - cartesianCenterCoorSi[1]);//force on y by O3
		Rf_O3_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[2][2] - cartesianCenterCoorSi[2]);//force on z by O3

		fcoe = kbSi * diff_Si_O4 / length_Si_O4;
		Rf_O4_on_Si[0] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][0] - cartesianCenterCoorSi[0]);//force on x by O4
		Rf_O4_on_Si[1] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][1] - cartesianCenterCoorSi[1]);//force on y by O4
		Rf_O4_on_Si[2] = fcoe * (realCartesianDirectBondedCoorOwithSi[3][2] - cartesianCenterCoorSi[2]);//force on z by O4

		//sum all the radial forces applied on Si
		Rf_on_Si[0] = Rf_O1_on_Si[0] + Rf_O2_on_Si[0] + Rf_O3_on_Si[0] + Rf_O4_on_Si[0];
		Rf_on_Si[1] = Rf_O1_on_Si[1] + Rf_O2_on_Si[1] + Rf_O3_on_Si[1] + Rf_O4_on_Si[1];
		Rf_on_Si[2] = Rf_O1_on_Si[2] + Rf_O2_on_Si[2] + Rf_O3_on_Si[2] + Rf_O4_on_Si[2];

		//store radial force into acceleration(need to know startpoint of Si)
		acceleration[i * 3 + startPoint[1]] += Rf_on_Si[0] / massSi;
		acceleration[i * 3 + 1 + startPoint[1]] += Rf_on_Si[1] / massSi;
		acceleration[i * 3 + 2 + startPoint[1]] += Rf_on_Si[2] / massSi;

		//store radial force into acceleration which just add minus sign to Si forces(need to know startpoint of O)
		acceleration[orderObondedwithSi[0] * 3 + startPoint[0]] -= Rf_O1_on_Si[0] / massO;//O1, x
		acceleration[orderObondedwithSi[0] * 3 + 1 + startPoint[0]] -= Rf_O1_on_Si[1] / massO;//O1, y
		acceleration[orderObondedwithSi[0] * 3 + 2 + startPoint[0]] -= Rf_O1_on_Si[2] / massO;//O1, z

		acceleration[orderObondedwithSi[1] * 3 + startPoint[0]] -= Rf_O2_on_Si[0] / massO;//O2, x
		acceleration[orderObondedwithSi[1] * 3 + 1 + startPoint[0]] -= Rf_O2_on_Si[1] / massO;//O2, y
		acceleration[orderObondedwithSi[1] * 3 + 2 + startPoint[0]] -= Rf_O2_on_Si[2] / massO;//O2, z

		acceleration[orderObondedwithSi[2] * 3 + startPoint[0]] -= Rf_O3_on_Si[0] / massO;//O3, x
		acceleration[orderObondedwithSi[2] * 3 + 1 + startPoint[0]] -= Rf_O3_on_Si[1] / massO;//O3, y
		acceleration[orderObondedwithSi[2] * 3 + 2 + startPoint[0]] -= Rf_O3_on_Si[2] / massO;//O3, z

		acceleration[orderObondedwithSi[3] * 3 + startPoint[0]] -= Rf_O4_on_Si[0] / massO;//O4, x
		acceleration[orderObondedwithSi[3] * 3 + 1 + startPoint[0]] -= Rf_O4_on_Si[1] / massO;//O4, y
		acceleration[orderObondedwithSi[3] * 3 + 2 + startPoint[0]] -= Rf_O4_on_Si[2] / massO;//O4, z



		//--------------------------------------------------------------------------------------------------------------

		/*calculate angular energy and angular force
		 * and store forces to acceleration */

		//coefficient
		T diff_Si_O1_O2 = cosSi_O1_O2 - cos0_Si;
		T diff_Si_O1_O3 = cosSi_O1_O3 - cos0_Si;
		T diff_Si_O1_O4 = cosSi_O1_O4 - cos0_Si;
		T diff_Si_O2_O3 = cosSi_O2_O3 - cos0_Si;
		T diff_Si_O2_O4 = cosSi_O2_O4 - cos0_Si;
		T diff_Si_O3_O4 = cosSi_O3_O4 - cos0_Si;

		//angular energy for Si
		angular_energy_Si += 0.5 * kthetaSi * b0Si * b0Si *
			(diff_Si_O1_O2*diff_Si_O1_O2 + diff_Si_O1_O3 * diff_Si_O1_O3 +
				diff_Si_O1_O4 * diff_Si_O1_O4 + diff_Si_O2_O3 * diff_Si_O2_O3 +
				diff_Si_O2_O4 * diff_Si_O2_O4 + diff_Si_O3_O4 * diff_Si_O3_O4);

		//angular force for Si
		fcoe = -1.0 * kthetaSi * b0Si * b0Si;

		/* force calculation (Si, O1 and O2) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O2
			fSiO1O2_on_O1[m] = fcoe * ((diff_Si_O1_O2 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O1_O2 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O2
			fSiO1O2_on_O2[m] = fcoe * ((diff_Si_O1_O2 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O2 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
			fSiO1O2_on_Si[m] = -(fSiO1O2_on_O1[m] + fSiO1O2_on_O2[m]);


		/* force calculation (Si, O1 and O3) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O3
			fSiO1O3_on_O1[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O1_O3 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O3
			fSiO1O3_on_O3[m] = fcoe * ((diff_Si_O1_O3 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O3 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O3 on Si
			fSiO1O3_on_Si[m] = -(fSiO1O3_on_O1[m] + fSiO1O3_on_O3[m]);


		/* force calculation (Si, O1 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O1 and O4
			fSiO1O4_on_O1[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O1) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O1_O4 * (realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1));
		for (int m = 0; m < 3; m++)//force between Si O1 and O4
			fSiO1O4_on_O4[m] = fcoe * ((diff_Si_O1_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[0][m] - cartesianCenterCoorSi[m]) / length_Si_O1 
				- cosSi_O1_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
			fSiO1O4_on_Si[m] = -(fSiO1O4_on_O1[m] + fSiO1O4_on_O4[m]);


		/* force calculation (Si, O2 and O3) */
		for (int m = 0; m < 3; m++)//force between Si O2 and O3
			fSiO2O3_on_O2[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O2_O3 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//force between Si O2 and O3
			fSiO2O3_on_O3[m] = fcoe * ((diff_Si_O2_O3 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O2_O3 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O3 on Si
			fSiO2O3_on_Si[m] = -(fSiO2O3_on_O2[m] + fSiO2O3_on_O3[m]);


		/* force calculation (Si, O2 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O2 and O4
			fSiO2O4_on_O2[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O2) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O2_O4 * (realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2));
		for (int m = 0; m < 3; m++)//force between Si O2 and O4
			fSiO2O4_on_O4[m] = fcoe * ((diff_Si_O2_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[1][m] - cartesianCenterCoorSi[m]) / length_Si_O2 
				- cosSi_O2_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
			fSiO2O4_on_Si[m] = -(fSiO2O4_on_O2[m] + fSiO2O4_on_O4[m]);


		/* force calculation (Si, O3 and O4) */
		for (int m = 0; m < 3; m++)//force between Si O3 and O4
			fSiO3O4_on_O3[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O3) *
			((realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4 
				- cosSi_O3_O4 * (realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3));
		for (int m = 0; m < 3; m++)//force between Si O3 and O4
			fSiO3O4_on_O4[m] = fcoe * ((diff_Si_O3_O4 / length_Si_O4) *
			((realCartesianDirectBondedCoorOwithSi[2][m] - cartesianCenterCoorSi[m]) / length_Si_O3 
				- cosSi_O3_O4 * (realCartesianDirectBondedCoorOwithSi[3][m] - cartesianCenterCoorSi[m]) / length_Si_O4));
		for (int m = 0; m < 3; m++)//final calculation for force between Si O3 and O4 on Si
			fSiO3O4_on_Si[m] = -(fSiO3O4_on_O3[m] + fSiO3O4_on_O4[m]);

		//store angular force into forceAcceleration(need to know startpoint of Si)
		for (int m = 0; m < 3; m++)
			acceleration[i * 3 + m + startPoint[1]] += (fSiO1O2_on_Si[m] + fSiO1O3_on_Si[m] + fSiO1O4_on_Si[m] +
				fSiO2O3_on_Si[m] + fSiO2O4_on_Si[m] + fSiO3O4_on_Si[m]) / massSi;


		//store angular force into forceAcceleration(need to know startpoint of O)
		for (int m = 0; m < 3; m++) {
			acceleration[orderObondedwithSi[0] * 3 + m + startPoint[0]] += (fSiO1O2_on_O1[m] + fSiO1O3_on_O1[m] + fSiO1O4_on_O1[m]) / massO;
			acceleration[orderObondedwithSi[1] * 3 + m + startPoint[0]] += (fSiO1O2_on_O2[m] + fSiO2O3_on_O2[m] + fSiO2O4_on_O2[m]) / massO;
			acceleration[orderObondedwithSi[2] * 3 + m + startPoint[0]] += (fSiO1O3_on_O3[m] + fSiO2O3_on_O3[m] + fSiO3O4_on_O3[m]) / massO;
			acceleration[orderObondedwithSi[3] * 3 + m + startPoint[0]] += (fSiO1O4_on_O4[m] + fSiO2O4_on_O4[m] + fSiO3O4_on_O4[m]) / massO;
		}

		//---------------------------------------------------------------------------------------
		
		
		// calculate repulsive energy and force between Si and Si
		for (int k = a; k < atomQuantity[1]; ++k)
		{
			//store all other Si coordiante in Direct(need to know startpoint of Si)
			directRepSiSi[0] = coordinates[k * 3 + startPoint[1]];
			directRepSiSi[1] = coordinates[k * 3 + 1 + startPoint[1]];
			directRepSiSi[2] = coordinates[k * 3 + 2 + startPoint[1]];

			//get real direct coor of Si near the center Si
			EFTool.directShortestCoor(directCenterCoorSi, directRepSiSi, directRepSiSishortest);

			//convert direct to cartesian 
			EFTool.directtoCart(directRepSiSishortest, cartesianRepSiSishortest);

			//get length between two Si atoms using Cartesian	
			repSiSilength = EFTool.distanceinCart(cartesianCenterCoorSi, cartesianRepSiSishortest);

			if (repSiSilength <= 3.10)
			{
				//calculate energy
				repulsive_energy_SiSiorSiPorPP += 0.5 * repulsiveCoeforSiSiorSiPorPP * (3.10 - repSiSilength)*(3.10 - repSiSilength)*(3.10 - repSiSilength);

				//calculate force
				for (int m = 0; m < 3; ++m)
				{
					frepSiSi[m] = -3 * repulsiveCoeforSiSiorSiPorPP *
						(3.10 - repSiSilength)*(3.10 - repSiSilength) / repSiSilength *
						(cartesianRepSiSishortest[m] - cartesianCenterCoorSi[m]);
					//store force in acceleration (need to know startpoint of Si)
					acceleration[i * 3 + m + startPoint[1]] += frepSiSi[m] / massSi;
					acceleration[k * 3 + m + startPoint[1]] -= frepSiSi[m] / massSi;
				}
			}
		}
		++a;


		// calculate repulsive energy and force between Si and P
		for (int k = 0; k < atomQuantity[2]; ++k)
		{
			//store all P coordiante in Direct(need to know startpoint of P)
			directRepSiP[0] = coordinates[k * 3 + startPoint[2]];
			directRepSiP[1] = coordinates[k * 3 + 1 + startPoint[2]];
			directRepSiP[2] = coordinates[k * 3 + 2 + startPoint[2]];

			//get real direct coor of P near the center Si
			EFTool.directShortestCoor(directCenterCoorSi, directRepSiP, directRepSiPshortest);

			//convert direct to cartesian 
			EFTool.directtoCart(directRepSiPshortest, cartesianRepSiPshortest);

			//get length between P and Si atoms using Cartesian	
			repSiPlength = EFTool.distanceinCart(cartesianCenterCoorSi, cartesianRepSiPshortest);

			if (repSiPlength <= 3.10)
			{
				//calculate energy
				repulsive_energy_SiSiorSiPorPP += 0.5 * repulsiveCoeforSiSiorSiPorPP * 
					(3.10 - repSiPlength)*(3.10 - repSiPlength)*(3.10 - repSiPlength);

				//calculate force
				for (int m = 0; m < 3; ++m)
				{
					frepSiP[m] = -3 * repulsiveCoeforSiSiorSiPorPP *
						(3.10 - repSiPlength)*(3.10 - repSiPlength) / repSiPlength *
						(cartesianRepSiPshortest[m] - cartesianCenterCoorSi[m]);
					//store force in acceleration (need to know startpoint of P)
					acceleration[i * 3 + m + startPoint[1]] += frepSiP[m] / massSi;
					acceleration[k * 3 + m + startPoint[2]] -= frepSiP[m] / massP;
				}
			}
		}
		


	}

	/*----------------------------------------------------------------------------------------------------------------------------------------



				energy for O



	-----------------------------------------------------------------------------------------------------------------------------------------*/

	//calculate angular energy of O
	int b = 1;
	for (int i = 0; i < atomQuantity[0]; ++i)
	{
		if (i != doubleO[0][0] && i != doubleO[0][1])
		{
			//store O position(need to know the startpoint of O)
			directCenterCoorO[0] = coordinates[i * 3 + startPoint[0]];
			directCenterCoorO[1] = coordinates[i * 3 + 1 + startPoint[0]];
			directCenterCoorO[2] = coordinates[i * 3 + 2 + startPoint[0]];

			//get order of Si or P bonded with center O
			for (int j = 0; j < 2; ++j)
			{
				orderSibondedwithO[j] = groupO[i][j];
				for (int k = 0; k < 3; ++k)//get direct coordinates of Si(need to know the startpoint of Si)
					directBondedCoorSiwithO[j][k] = coordinates[orderSibondedwithO[j] * 3 + k + startPoint[1]];
			}

			//get real direct coordinate of Si that bonded with O
			for (int j = 0; j < 2; ++j)
				EFTool.directShortestCoor(directCenterCoorO, directBondedCoorSiwithO[j], realDirectBondedCoorSiwithO[j]);

			//convert coordinates from Direct to Cartesian
			EFTool.directtoCart(directCenterCoorO, cartesianCenterCoorO);//for center O
			for (int j = 0; j < 2; ++j)//for bonded Si
				EFTool.directtoCart(realDirectBondedCoorSiwithO[j], realCartesianDirectBondedCoorSiwithO[j]);

			//length
			T length_O_Si1 = EFTool.distanceinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[0]);
			T length_O_Si2 = EFTool.distanceinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[1]);

			//angle
			T cosO_Si1_Si2 = EFTool.angleinCart(cartesianCenterCoorO, realCartesianDirectBondedCoorSiwithO[0], realCartesianDirectBondedCoorSiwithO[1], length_O_Si1, length_O_Si2);

			//calculate angular energy
			T diff_O_Si1_Si2 = cosO_Si1_Si2 - cos0_O;
			T fcoe = 0.0;
			//all atoms are P
			if (orderSibondedwithO[0] >= atomQuantity[1] && orderSibondedwithO[1] >= atomQuantity[1])
			{
				angular_energy_O += 0.5 * kthetaO * b0P * b0P * diff_O_Si1_Si2 * diff_O_Si1_Si2;
				fcoe = -1.0 * kthetaO * b0P * b0P;
			}
			//all atoms are Si
			else if (orderSibondedwithO[0] < atomQuantity[1] && orderSibondedwithO[1] < atomQuantity[1])
			{
				angular_energy_O += 0.5 * kthetaO * b0Si * b0Si * diff_O_Si1_Si2 * diff_O_Si1_Si2;
				fcoe = -1.0 * kthetaO * b0Si * b0Si;
			}
			//one is Si and another one is P
			else
			{
				angular_energy_O += 0.5 * kthetaO * b0P * b0Si * diff_O_Si1_Si2 * diff_O_Si1_Si2;
				fcoe = -1.0 * kthetaO * b0Si * b0P;
			}


			//force calculation (Si, O2 and O4)
			for (int m = 0; m < 3; m++)//force between O Si1 and Si2
				fOSi1Si2_on_Si1[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si1) *
				((realCartesianDirectBondedCoorSiwithO[1][m] - cartesianCenterCoorO[m]) / length_O_Si2
					- cosO_Si1_Si2 * (realCartesianDirectBondedCoorSiwithO[0][m] - cartesianCenterCoorO[m]) / length_O_Si1));
			for (int m = 0; m < 3; m++)//force between Si O2 and O4
				fOSi1Si2_on_Si2[m] = fcoe * ((diff_O_Si1_Si2 / length_O_Si2) *
				((realCartesianDirectBondedCoorSiwithO[0][m] - cartesianCenterCoorO[m]) / length_O_Si1
					- cosO_Si1_Si2 * (realCartesianDirectBondedCoorSiwithO[1][m] - cartesianCenterCoorO[m]) / length_O_Si2));
			for (int m = 0; m < 3; m++)//final calculation for force between Si O2 and O4 on Si
				fOSi1Si2_on_O[m] = -(fOSi1Si2_on_Si1[m] + fOSi1Si2_on_Si2[m]);

			//store angular force into atomAccelelation(need to know the startpoint of O)
			for (int m = 0; m < 3; m++)
				acceleration[i * 3 + m + startPoint[0]] += fOSi1Si2_on_O[m] / massO;




			//store angular force into atomAccelelation(need to know the startpoint of Si)
			//all atoms are P
			if (orderSibondedwithO[0] >= atomQuantity[1] && orderSibondedwithO[1] >= atomQuantity[1])
			{
				for (int m = 0; m < 3; m++)
				{
					acceleration[orderSibondedwithO[0] * 3 + m + startPoint[2]] += fOSi1Si2_on_Si1[m] / massP;
					acceleration[orderSibondedwithO[1] * 3 + m + startPoint[2]] += fOSi1Si2_on_Si2[m] / massP;
				}
			}
			//all atoms are Si
			else if (orderSibondedwithO[0] < atomQuantity[1] && orderSibondedwithO[1] < atomQuantity[1])
			{
				for (int m = 0; m < 3; m++)
				{
					acceleration[orderSibondedwithO[0] * 3 + m + startPoint[1]] += fOSi1Si2_on_Si1[m] / massSi;
					acceleration[orderSibondedwithO[1] * 3 + m + startPoint[1]] += fOSi1Si2_on_Si2[m] / massSi;
				}
			}
			//first one is Si and second one is P
			else if (orderSibondedwithO[0] < atomQuantity[1] && orderSibondedwithO[1] >= atomQuantity[1])
			{
				for (int m = 0; m < 3; m++)
				{
					acceleration[orderSibondedwithO[0] * 3 + m + startPoint[1]] += fOSi1Si2_on_Si1[m] / massSi;
					acceleration[orderSibondedwithO[1] * 3 + m + startPoint[2]] += fOSi1Si2_on_Si2[m] / massP;
				}
			}
			//first one is P and second one is Si
			else
			{
				for (int m = 0; m < 3; m++)
				{
					acceleration[orderSibondedwithO[0] * 3 + m + startPoint[1]] += fOSi1Si2_on_Si1[m] / massP;
					acceleration[orderSibondedwithO[1] * 3 + m + startPoint[2]] += fOSi1Si2_on_Si2[m] / massSi;
				}
			}


		}


		//calculate repulsive energy between O and O
		for (int j = b; j < atomQuantity[0]; ++j)//get direct coordinates of center O
		{
			directRepOO[0] = coordinates[j * 3 + startPoint[0]];
			directRepOO[1] = coordinates[j * 3 + 1 + startPoint[0]];
			directRepOO[2] = coordinates[j * 3 + 2 + startPoint[0]];

			//find real direct coordinates of O near center O
			EFTool.directShortestCoor(directCenterCoorO, directRepOO, directRepOOshortest);

			//convert coordiantes from Direct to Cartesian for RepOO
			EFTool.directtoCart(directRepOOshortest, cartesianRepOOshortest);

			//get length
			repOOlength = EFTool.distanceinCart(cartesianCenterCoorO, cartesianRepOOshortest);

			if (repOOlength <= 2.533)
			{
				//calculate energy
				repulsive_energy_O_O += 0.5*repulsiveCoeforOO*(2.533 - repOOlength)*(2.533 - repOOlength)*(2.533 - repOOlength);

				//calculate force for O and O
				for (int m = 0; m < 3; ++m)
				{
					frepOO[m] = -3 * repulsiveCoeforOO * (2.533 - repOOlength)*(2.533 - repOOlength) / repOOlength *
						(cartesianRepOOshortest[m] - cartesianCenterCoorO[m]);
					//store acceleration to acceleration(need to know the startpoint of O)
					acceleration[i * 3 + m + startPoint[0]] += frepOO[m] / massO;
					acceleration[j * 3 + m + startPoint[0]] -= frepOO[m] / massO;
				}
			}
		}
		++b;

	}



	/*----------------------------------------------------------------------------------------------------------------------------------------
	


				energy for P



	-----------------------------------------------------------------------------------------------------------------------------------------*/
	int c = 1;
	for (int i = 0; i < atomQuantity[2]; ++i)
	{
		
		//P has 3 single bonds with 1 double bond 
		if (doubleO[1][i] == 4)
		{
			//store P position using atomCoordiante(need to know start point of P atom)
			directCenterCoorP[0] = coordinates[i * 3 + startPoint[2]];
			directCenterCoorP[1] = coordinates[i * 3 + 1 + startPoint[2]];
			directCenterCoorP[2] = coordinates[i * 3 + 2 + startPoint[2]];

			//store bonded O position with P 
			for (int j = 0; j < 4; ++j)
			{
				orderObondedwithP[j] = groupP[i][j];//get bonded O orders

				//get coordinates of each bonded O
				directBondedCoorOwithP[j][0] = coordinates[orderObondedwithP[j] * 3 + startPoint[0]];//x
				directBondedCoorOwithP[j][1] = coordinates[orderObondedwithP[j] * 3 + 1 + startPoint[0]];//y
				directBondedCoorOwithP[j][2] = coordinates[orderObondedwithP[j] * 3 + 2 + startPoint[0]];//z
			}

			//get real direct shortest coordinates of bonded O with center Si
			for (int j = 0; j < 4; ++j)
				EFTool.directShortestCoor(directCenterCoorP, directBondedCoorOwithP[j], realDirectBondedCoorOwithP[j]);


			//--------------------------------------------------------------------------------------------------------------

			//convert Direct to Cartesian coordinates
			EFTool.directtoCart(directCenterCoorP, cartesianCenterCoorP);//convert center P coordiantes from Direct to Cartesian
			for (int j = 0; j < 4; ++j)//convert bonded O from Direct to Cartsesian 
				EFTool.directtoCart(realDirectBondedCoorOwithP[j], realCartesianDirectBondedCoorOwithP[j]);

			//length between center P and bonded O
			T length_P_O1 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0]);
			T length_P_O2 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1]);
			T length_P_O3 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[2]);
			T length_P_O4 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[3]);//double bond

			//angle
			T cosP_O1_O2 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[1], length_P_O1, length_P_O2);
			T cosP_O1_O3 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[2], length_P_O1, length_P_O3);
			T cosP_O1_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[3], length_P_O1, length_P_O4);
			T cosP_O2_O3 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1], realCartesianDirectBondedCoorOwithP[2], length_P_O2, length_P_O3);
			T cosP_O2_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1], realCartesianDirectBondedCoorOwithP[3], length_P_O2, length_P_O4);
			T cosP_O3_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[2], realCartesianDirectBondedCoorOwithP[3], length_P_O3, length_P_O4);

			//--------------------------------------------------------------------------------------------------------------

			/*calculate radial energy and radial force
			 * and store forces to  */

			 //coefficient
			T diff_P_O1 = length_P_O1 - b0P;
			T diff_P_O2 = length_P_O2 - b0P;
			T diff_P_O3 = length_P_O3 - b0P;
			T diff_P_O4 = length_P_O4 - b0PDoubleBond;//double bond

			//radial energy
			radial_energyPandO += 0.5 * kbP *(diff_P_O1*diff_P_O1 + diff_P_O2 * diff_P_O2 + diff_P_O3 * diff_P_O3)
				+ 0.5 * kbPDoubleBond * (diff_P_O4 * diff_P_O4);

			//radial force for P
			T fcoe = kbP * diff_P_O1 / length_P_O1;
			Rf_O1_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[0][0] - cartesianCenterCoorP[0]);//force on x by O1
			Rf_O1_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[0][1] - cartesianCenterCoorP[1]);//force on y by O1
			Rf_O1_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[0][2] - cartesianCenterCoorP[2]);//force on z by O1

			fcoe = kbP * diff_P_O2 / length_P_O2;
			Rf_O2_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[1][0] - cartesianCenterCoorP[0]);//force on x by O2
			Rf_O2_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[1][1] - cartesianCenterCoorP[1]);//force on y by O2
			Rf_O2_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[1][2] - cartesianCenterCoorP[2]);//force on z by O2

			fcoe = kbP * diff_P_O3 / length_P_O3;
			Rf_O3_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[2][0] - cartesianCenterCoorP[0]);//force on x by O3
			Rf_O3_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[2][1] - cartesianCenterCoorP[1]);//force on y by O3
			Rf_O3_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[2][2] - cartesianCenterCoorP[2]);//force on z by O3

			fcoe = kbPDoubleBond * diff_P_O4 / length_P_O4;
			Rf_O4_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[3][0] - cartesianCenterCoorP[0]);//force on x by O4
			Rf_O4_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[3][1] - cartesianCenterCoorP[1]);//force on y by O4
			Rf_O4_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[3][2] - cartesianCenterCoorP[2]);//force on z by O4


			//sum all the radial forces applied on P
			Rf_on_P[0] = Rf_O1_on_P[0] + Rf_O2_on_P[0] + Rf_O3_on_P[0] + Rf_O4_on_P[0];
			Rf_on_P[1] = Rf_O1_on_P[1] + Rf_O2_on_P[1] + Rf_O3_on_P[1] + Rf_O4_on_P[1];
			Rf_on_P[2] = Rf_O1_on_P[2] + Rf_O2_on_P[2] + Rf_O3_on_P[2] + Rf_O4_on_P[2];

			//store radial force into acceleration(need to know startpoint of P)
			acceleration[i * 3 + startPoint[2]] += Rf_on_P[0] / massP;
			acceleration[i * 3 + 1 + startPoint[2]] += Rf_on_P[1] / massP;
			acceleration[i * 3 + 2 + startPoint[2]] += Rf_on_P[2] / massP;

			//store radial force into acceleration for O which just add minus sign to P forces(need to know startpoint of O)
			acceleration[orderObondedwithP[0] * 3 + startPoint[0]] -= Rf_O1_on_P[0] / massO;//O1, x
			acceleration[orderObondedwithP[0] * 3 + 1 + startPoint[0]] -= Rf_O1_on_P[1] / massO;//O1, y
			acceleration[orderObondedwithP[0] * 3 + 2 + startPoint[0]] -= Rf_O1_on_P[2] / massO;//O1, z

			acceleration[orderObondedwithP[1] * 3 + startPoint[0]] -= Rf_O2_on_P[0] / massO;//O2, x
			acceleration[orderObondedwithP[1] * 3 + 1 + startPoint[0]] -= Rf_O2_on_P[1] / massO;//O2, y
			acceleration[orderObondedwithP[1] * 3 + 2 + startPoint[0]] -= Rf_O2_on_P[2] / massO;//O2, z

			acceleration[orderObondedwithP[2] * 3 + startPoint[0]] -= Rf_O3_on_P[0] / massO;//O3, x
			acceleration[orderObondedwithP[2] * 3 + 1 + startPoint[0]] -= Rf_O3_on_P[1] / massO;//O3, y
			acceleration[orderObondedwithP[2] * 3 + 2 + startPoint[0]] -= Rf_O3_on_P[2] / massO;//O3, z

			acceleration[orderObondedwithP[3] * 3 + startPoint[0]] -= Rf_O4_on_P[0] / massO;//O4, x
			acceleration[orderObondedwithP[3] * 3 + 1 + startPoint[0]] -= Rf_O4_on_P[1] / massO;//O4, y
			acceleration[orderObondedwithP[3] * 3 + 2 + startPoint[0]] -= Rf_O4_on_P[2] / massO;//O4, z

			//--------------------------------------------------------------------------------------------------------------

			/*calculate angular energy and angular force
			 * and store forces to acceleration */

			 //coefficient
			T diff_P_O1_O2 = cosP_O1_O2 - cos0_P;
			T diff_P_O1_O3 = cosP_O1_O3 - cos0_P;
			T diff_P_O1_O4 = cosP_O1_O4 - cos0_PDoubleBond;//double bond
			T diff_P_O2_O3 = cosP_O2_O3 - cos0_P;
			T diff_P_O2_O4 = cosP_O2_O4 - cos0_PDoubleBond;//double bond
			T diff_P_O3_O4 = cosP_O3_O4 - cos0_PDoubleBond;//double bond

			//angular energy for Si
			angular_energy_P += 0.5 * kthetaP * b0P * b0P *
				(diff_P_O1_O2*diff_P_O1_O2 + diff_P_O1_O3 * diff_P_O1_O3 + diff_P_O2_O3 * diff_P_O2_O3) +
				0.5 * kthetaPDoubleBond * b0PDoubleBond * b0P *
				(diff_P_O1_O4 * diff_P_O1_O4 + diff_P_O2_O4 * diff_P_O2_O4 + diff_P_O3_O4 * diff_P_O3_O4);

			//angular force for P
			fcoe = -1.0 * kthetaP * b0P * b0P;

			/* force calculation (P, O1 and O2) */
			for (int m = 0; m < 3; m++)//force between P O1 and O2
				fPO1O2_on_O1[m] = fcoe * ((diff_P_O1_O2 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O1_O2 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O2
				fPO1O2_on_O2[m] = fcoe * ((diff_P_O1_O2 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O2 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
				fPO1O2_on_P[m] = -(fPO1O2_on_O1[m] + fPO1O2_on_O2[m]);


			/* force calculation (P, O1 and O3) */
			for (int m = 0; m < 3; m++)//force between P O1 and O3
				fPO1O3_on_O1[m] = fcoe * ((diff_P_O1_O3 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O1_O3 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O3
				fPO1O3_on_O3[m] = fcoe * ((diff_P_O1_O3 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O3 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//final calculation for force between P O1 and O3 on P
				fPO1O3_on_P[m] = -(fPO1O3_on_O1[m] + fPO1O3_on_O3[m]);

			/* force calculation (P, O2 and O3) */
			for (int m = 0; m < 3; m++)//force between P O2 and O3
				fPO2O3_on_O2[m] = fcoe * ((diff_P_O2_O3 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O2_O3 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//force between P O2 and O3
				fPO2O3_on_O3[m] = fcoe * ((diff_P_O2_O3 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O2_O3 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//final calculation for force between P O2 and O3 on P
				fPO2O3_on_P[m] = -(fPO2O3_on_O2[m] + fPO2O3_on_O3[m]);


			fcoe = -1.0 * kthetaPDoubleBond * b0PDoubleBond * b0P;

			/* force calculation (P, O1 and O4) */
			for (int m = 0; m < 3; m++)//force between P O1 and O4
				fPO1O4_on_O1[m] = fcoe * ((diff_P_O1_O4 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O1_O4 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O4
				fPO1O4_on_O4[m] = fcoe * ((diff_P_O1_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O1 and O2 on P
				fPO1O4_on_P[m] = -(fPO1O4_on_O1[m] + fPO1O4_on_O4[m]);

			/* force calculation (P, O2 and O4) */
			for (int m = 0; m < 3; m++)//force between P O2 and O4
				fPO2O4_on_O2[m] = fcoe * ((diff_P_O2_O4 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O2_O4 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//force between P O2 and O4
				fPO2O4_on_O4[m] = fcoe * ((diff_P_O2_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O2_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O2 and O4 on P
				fPO2O4_on_P[m] = -(fPO2O4_on_O2[m] + fPO2O4_on_O4[m]);

			/* force calculation (P, O3 and O4) */
			for (int m = 0; m < 3; m++)//force between P O3 and O4
				fPO3O4_on_O3[m] = fcoe * ((diff_P_O3_O4 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O3_O4 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//force between P O3 and O4
				fPO3O4_on_O4[m] = fcoe * ((diff_P_O3_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O3_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O3 and O4 on P
				fPO3O4_on_P[m] = -(fPO3O4_on_O3[m] + fPO3O4_on_O4[m]);


			//store angular force into forceAcceleration(need to know startpoint of P)
			for (int m = 0; m < 3; m++)
				acceleration[i * 3 + m + startPoint[2]] += (fPO1O2_on_P[m] + fPO1O3_on_P[m] + fPO1O4_on_P[m]
					+ fPO2O3_on_P[m] + fPO2O4_on_P[m] + fPO3O4_on_P[m]) / massP;


			//store angular force into forceAcceleration(need to know startpoint of O)
			for (int m = 0; m < 3; m++) {
				acceleration[orderObondedwithP[0] * 3 + m + startPoint[0]] += (fPO1O2_on_O1[m] + fPO1O3_on_O1[m] + fPO1O4_on_O1[m] ) / massO;
				acceleration[orderObondedwithP[1] * 3 + m + startPoint[0]] += (fPO1O2_on_O2[m] + fPO2O3_on_O2[m] + fPO2O4_on_O2[m] ) / massO;
				acceleration[orderObondedwithP[2] * 3 + m + startPoint[0]] += (fPO1O3_on_O3[m] + fPO2O3_on_O3[m] + fPO3O4_on_O3[m] ) / massO;
				acceleration[orderObondedwithP[3] * 3 + m + startPoint[0]] += (fPO1O4_on_O4[m] + fPO2O4_on_O4[m] + fPO3O4_on_O4[m] ) / massO;
			}
		}
		
		//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


		//P has 5 bonds
		else if (doubleO[1][i] == 5)
		{
			//store P position using atomCoordiante(need to know start point of P atom)
			directCenterCoorP[0] = coordinates[i * 3 + startPoint[2]];
			directCenterCoorP[1] = coordinates[i * 3 + 1 + startPoint[2]];
			directCenterCoorP[2] = coordinates[i * 3 + 2 + startPoint[2]];

			//store bonded O position with P 
			for (int j = 0; j < 5; ++j)
			{
				orderObondedwithP[j] = groupP[i][j];//get bonded O orders

				//get coordinates of each bonded O
				directBondedCoorOwithP[j][0] = coordinates[orderObondedwithP[j] * 3 + startPoint[0]];//x
				directBondedCoorOwithP[j][1] = coordinates[orderObondedwithP[j] * 3 + 1 + startPoint[0]];//y
				directBondedCoorOwithP[j][2] = coordinates[orderObondedwithP[j] * 3 + 2 + startPoint[0]];//z
			}

			//get real direct shortest coordinates of bonded O with center Si
			for (int j = 0; j < 5; ++j)
				EFTool.directShortestCoor(directCenterCoorP, directBondedCoorOwithP[j], realDirectBondedCoorOwithP[j]);


			//--------------------------------------------------------------------------------------------------------------

			//convert Direct to Cartesian coordinates
			EFTool.directtoCart(directCenterCoorP, cartesianCenterCoorP);//convert center P coordiantes from Direct to Cartesian
			for (int j = 0; j < 5; ++j)
				EFTool.directtoCart(realDirectBondedCoorOwithP[j], realCartesianDirectBondedCoorOwithP[j]);

			//length between center P and bonded O
			T length_P_O1 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0]);
			T length_P_O2 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1]);
			T length_P_O3 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[2]);
			T length_P_O4 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[3]);
			T length_P_O5 = EFTool.distanceinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[4]);


			//std::cout << length_Si_O1 << "    " << length_Si_O2 << "    " << length_Si_O3 << "    " << length_Si_O4 << std::endl;

			//angle
			T cosP_O1_O2 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[1], length_P_O1, length_P_O2);
			T cosP_O1_O3 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[2], length_P_O1, length_P_O3);
			T cosP_O1_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[3], length_P_O1, length_P_O4);
			T cosP_O1_O5 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[0], realCartesianDirectBondedCoorOwithP[4], length_P_O1, length_P_O5);
			T cosP_O2_O3 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1], realCartesianDirectBondedCoorOwithP[2], length_P_O2, length_P_O3);
			T cosP_O2_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1], realCartesianDirectBondedCoorOwithP[3], length_P_O2, length_P_O4);
			T cosP_O2_O5 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[1], realCartesianDirectBondedCoorOwithP[4], length_P_O2, length_P_O5);
			T cosP_O3_O4 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[2], realCartesianDirectBondedCoorOwithP[3], length_P_O3, length_P_O4);
			T cosP_O3_O5 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[2], realCartesianDirectBondedCoorOwithP[4], length_P_O3, length_P_O5);
			T cosP_O4_O5 = EFTool.angleinCart(cartesianCenterCoorP, realCartesianDirectBondedCoorOwithP[3], realCartesianDirectBondedCoorOwithP[4], length_P_O4, length_P_O5);



			//--------------------------------------------------------------------------------------------------------------

			/*calculate radial energy and radial force
			 * and store forces to  */

			 //coefficient
			T diff_P_O1 = length_P_O1 - b0P;
			T diff_P_O2 = length_P_O2 - b0P;
			T diff_P_O3 = length_P_O3 - b0P;
			T diff_P_O4 = length_P_O4 - b0P;
			T diff_P_O5 = length_P_O5 - b0P;

			//radial energy
			radial_energyPandO += 0.5 * kbP *
				(diff_P_O1*diff_P_O1 + diff_P_O2 * diff_P_O2 +
					diff_P_O3 * diff_P_O3 + diff_P_O4 * diff_P_O4 +
					diff_P_O5 * diff_P_O5);

			//radial force for P
			T fcoe = kbP * diff_P_O1 / length_P_O1;
			Rf_O1_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[0][0] - cartesianCenterCoorP[0]);//force on x by O1
			Rf_O1_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[0][1] - cartesianCenterCoorP[1]);//force on y by O1
			Rf_O1_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[0][2] - cartesianCenterCoorP[2]);//force on z by O1

			fcoe = kbP * diff_P_O2 / length_P_O2;
			Rf_O2_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[1][0] - cartesianCenterCoorP[0]);//force on x by O2
			Rf_O2_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[1][1] - cartesianCenterCoorP[1]);//force on y by O2
			Rf_O2_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[1][2] - cartesianCenterCoorP[2]);//force on z by O2

			fcoe = kbP * diff_P_O3 / length_P_O3;
			Rf_O3_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[2][0] - cartesianCenterCoorP[0]);//force on x by O3
			Rf_O3_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[2][1] - cartesianCenterCoorP[1]);//force on y by O3
			Rf_O3_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[2][2] - cartesianCenterCoorP[2]);//force on z by O3

			fcoe = kbP * diff_P_O4 / length_P_O4;
			Rf_O4_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[3][0] - cartesianCenterCoorP[0]);//force on x by O4
			Rf_O4_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[3][1] - cartesianCenterCoorP[1]);//force on y by O4
			Rf_O4_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[3][2] - cartesianCenterCoorP[2]);//force on z by O4

			fcoe = kbP * diff_P_O5 / length_P_O5;
			Rf_O5_on_P[0] = fcoe * (realCartesianDirectBondedCoorOwithP[4][0] - cartesianCenterCoorP[0]);//force on x by O5
			Rf_O5_on_P[1] = fcoe * (realCartesianDirectBondedCoorOwithP[4][1] - cartesianCenterCoorP[1]);//force on y by O5
			Rf_O5_on_P[2] = fcoe * (realCartesianDirectBondedCoorOwithP[4][2] - cartesianCenterCoorP[2]);//force on z by O5

			//sum all the radial forces applied on Si
			Rf_on_P[0] = Rf_O1_on_P[0] + Rf_O2_on_P[0] + Rf_O3_on_P[0] + Rf_O4_on_P[0] + Rf_O5_on_P[0];
			Rf_on_P[1] = Rf_O1_on_P[1] + Rf_O2_on_P[1] + Rf_O3_on_P[1] + Rf_O4_on_P[1] + Rf_O5_on_P[1];
			Rf_on_P[2] = Rf_O1_on_P[2] + Rf_O2_on_P[2] + Rf_O3_on_P[2] + Rf_O4_on_P[2] + Rf_O5_on_P[2];

			//store radial force into acceleration(need to know startpoint of P)
			acceleration[i * 3 + startPoint[2]] += Rf_on_P[0] / massP;
			acceleration[i * 3 + 1 + startPoint[2]] += Rf_on_P[1] / massP;
			acceleration[i * 3 + 2 + startPoint[2]] += Rf_on_P[2] / massP;

			//store radial force into acceleration for O which just add minus sign to P forces(need to know startpoint of O)
			acceleration[orderObondedwithP[0] * 3 + startPoint[0]] -= Rf_O1_on_P[0] / massO;//O1, x
			acceleration[orderObondedwithP[0] * 3 + 1 + startPoint[0]] -= Rf_O1_on_P[1] / massO;//O1, y
			acceleration[orderObondedwithP[0] * 3 + 2 + startPoint[0]] -= Rf_O1_on_P[2] / massO;//O1, z

			acceleration[orderObondedwithP[1] * 3 + startPoint[0]] -= Rf_O2_on_P[0] / massO;//O2, x
			acceleration[orderObondedwithP[1] * 3 + 1 + startPoint[0]] -= Rf_O2_on_P[1] / massO;//O2, y
			acceleration[orderObondedwithP[1] * 3 + 2 + startPoint[0]] -= Rf_O2_on_P[2] / massO;//O2, z

			acceleration[orderObondedwithP[2] * 3 + startPoint[0]] -= Rf_O3_on_P[0] / massO;//O3, x
			acceleration[orderObondedwithP[2] * 3 + 1 + startPoint[0]] -= Rf_O3_on_P[1] / massO;//O3, y
			acceleration[orderObondedwithP[2] * 3 + 2 + startPoint[0]] -= Rf_O3_on_P[2] / massO;//O3, z

			acceleration[orderObondedwithP[3] * 3 + startPoint[0]] -= Rf_O4_on_P[0] / massO;//O4, x
			acceleration[orderObondedwithP[3] * 3 + 1 + startPoint[0]] -= Rf_O4_on_P[1] / massO;//O4, y
			acceleration[orderObondedwithP[3] * 3 + 2 + startPoint[0]] -= Rf_O4_on_P[2] / massO;//O4, z

			acceleration[orderObondedwithP[4] * 3 + startPoint[0]] -= Rf_O5_on_P[0] / massO;//O5, x
			acceleration[orderObondedwithP[4] * 3 + 1 + startPoint[0]] -= Rf_O5_on_P[1] / massO;//O5, y
			acceleration[orderObondedwithP[4] * 3 + 2 + startPoint[0]] -= Rf_O5_on_P[2] / massO;//O5, z


			//--------------------------------------------------------------------------------------------------------------

			/*calculate angular energy and angular force
			 * and store forces to acceleration */

			 //coefficient
			T diff_P_O1_O2 = cosP_O1_O2 - cos0_P;
			T diff_P_O1_O3 = cosP_O1_O3 - cos0_P;
			T diff_P_O1_O4 = cosP_O1_O4 - cos0_P;
			T diff_P_O1_O5 = cosP_O1_O5 - cos0_P;
			T diff_P_O2_O3 = cosP_O2_O3 - cos0_P;
			T diff_P_O2_O4 = cosP_O2_O4 - cos0_P;
			T diff_P_O2_O5 = cosP_O2_O5 - cos0_P;
			T diff_P_O3_O4 = cosP_O3_O4 - cos0_P;
			T diff_P_O3_O5 = cosP_O3_O5 - cos0_P;
			T diff_P_O4_O5 = cosP_O4_O5 - cos0_P;

			//angular energy for Si
			angular_energy_P += 0.5 * kthetaP * b0P * b0P *
				(diff_P_O1_O2*diff_P_O1_O2 + diff_P_O1_O3 * diff_P_O1_O3 +
					diff_P_O1_O4 * diff_P_O1_O4 + diff_P_O1_O5 * diff_P_O1_O5 +
					diff_P_O2_O3 * diff_P_O2_O3 + diff_P_O2_O4 * diff_P_O2_O4 +
					diff_P_O2_O5 * diff_P_O2_O5 + diff_P_O3_O4 * diff_P_O3_O4 +
					diff_P_O3_O5 * diff_P_O3_O5 + diff_P_O4_O5 * diff_P_O4_O5);

			//angular force for P
			fcoe = -1.0 * kthetaP * b0P * b0P;

			/* force calculation (P, O1 and O2) */
			for (int m = 0; m < 3; m++)//force between P O1 and O2
				fPO1O2_on_O1[m] = fcoe * ((diff_P_O1_O2 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O1_O2 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O2
				fPO1O2_on_O2[m] = fcoe * ((diff_P_O1_O2 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O2 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//final calculation for force between Si O1 and O2 on Si
				fPO1O2_on_P[m] = -(fPO1O2_on_O1[m] + fPO1O2_on_O2[m]);


			/* force calculation (P, O1 and O3) */
			for (int m = 0; m < 3; m++)//force between P O1 and O3
				fPO1O3_on_O1[m] = fcoe * ((diff_P_O1_O3 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O1_O3 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O3
				fPO1O3_on_O3[m] = fcoe * ((diff_P_O1_O3 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O3 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//final calculation for force between P O1 and O3 on P
				fPO1O3_on_P[m] = -(fPO1O3_on_O1[m] + fPO1O3_on_O3[m]);


			/* force calculation (P, O1 and O4) */
			for (int m = 0; m < 3; m++)//force between P O1 and O4
				fPO1O4_on_O1[m] = fcoe * ((diff_P_O1_O4 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O1_O4 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O4
				fPO1O4_on_O4[m] = fcoe * ((diff_P_O1_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O1 and O2 on P
				fPO1O4_on_P[m] = -(fPO1O4_on_O1[m] + fPO1O4_on_O4[m]);

			/* force calculation (P, O1 and O5) */
			for (int m = 0; m < 3; m++)//force between P O1 and O5
				fPO1O5_on_O1[m] = fcoe * ((diff_P_O1_O5 / length_P_O1) *
				((realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5
					- cosP_O1_O5 * (realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1));
			for (int m = 0; m < 3; m++)//force between P O1 and O5
				fPO1O5_on_O5[m] = fcoe * ((diff_P_O1_O5 / length_P_O5) *
				((realCartesianDirectBondedCoorOwithP[0][m] - cartesianCenterCoorP[m]) / length_P_O1
					- cosP_O1_O5 * (realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5));
			for (int m = 0; m < 3; m++)//final calculation for force between P O1 and O2 on P
				fPO1O5_on_P[m] = -(fPO1O5_on_O1[m] + fPO1O5_on_O5[m]);


			/* force calculation (P, O2 and O3) */
			for (int m = 0; m < 3; m++)//force between P O2 and O3
				fPO2O3_on_O2[m] = fcoe * ((diff_P_O2_O3 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O2_O3 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//force between P O2 and O3
				fPO2O3_on_O3[m] = fcoe * ((diff_P_O2_O3 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O2_O3 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//final calculation for force between P O2 and O3 on P
				fPO2O3_on_P[m] = -(fPO2O3_on_O2[m] + fPO2O3_on_O3[m]);


			/* force calculation (P, O2 and O4) */
			for (int m = 0; m < 3; m++)//force between P O2 and O4
				fPO2O4_on_O2[m] = fcoe * ((diff_P_O2_O4 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O2_O4 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//force between P O2 and O4
				fPO2O4_on_O4[m] = fcoe * ((diff_P_O2_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O2_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O2 and O4 on P
				fPO2O4_on_P[m] = -(fPO2O4_on_O2[m] + fPO2O4_on_O4[m]);

			/* force calculation (P, O2 and O5) */
			for (int m = 0; m < 3; m++)//force between P O2 and O5
				fPO2O5_on_O2[m] = fcoe * ((diff_P_O2_O5 / length_P_O2) *
				((realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5
					- cosP_O2_O5 * (realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2));
			for (int m = 0; m < 3; m++)//force between P O2 and O5
				fPO2O5_on_O5[m] = fcoe * ((diff_P_O2_O5 / length_P_O5) *
				((realCartesianDirectBondedCoorOwithP[1][m] - cartesianCenterCoorP[m]) / length_P_O2
					- cosP_O2_O5 * (realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5));
			for (int m = 0; m < 3; m++)//final calculation for force between P O2 and O5 on P
				fPO2O5_on_P[m] = -(fPO2O5_on_O2[m] + fPO2O5_on_O5[m]);


			/* force calculation (P, O3 and O4) */
			for (int m = 0; m < 3; m++)//force between P O3 and O4
				fPO3O4_on_O3[m] = fcoe * ((diff_P_O3_O4 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O3_O4 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//force between P O3 and O4
				fPO3O4_on_O4[m] = fcoe * ((diff_P_O3_O4 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O3_O4 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//final calculation for force between P O3 and O4 on P
				fPO3O4_on_P[m] = -(fPO3O4_on_O3[m] + fPO3O4_on_O4[m]);

			/* force calculation (P, O3 and O5) */
			for (int m = 0; m < 3; m++)//force between P O3 and O5
				fPO3O5_on_O3[m] = fcoe * ((diff_P_O3_O5 / length_P_O3) *
				((realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5
					- cosP_O3_O5 * (realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3));
			for (int m = 0; m < 3; m++)//force between P O3 and O5
				fPO3O5_on_O5[m] = fcoe * ((diff_P_O3_O5 / length_P_O5) *
				((realCartesianDirectBondedCoorOwithP[2][m] - cartesianCenterCoorP[m]) / length_P_O3
					- cosP_O3_O5 * (realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5));
			for (int m = 0; m < 3; m++)//final calculation for force between P O3 and O5 on P
				fPO3O5_on_P[m] = -(fPO3O5_on_O3[m] + fPO3O5_on_O5[m]);

			/* force calculation (P, O4 and O5) */
			for (int m = 0; m < 3; m++)//force between P O4 and O5
				fPO4O5_on_O4[m] = fcoe * ((diff_P_O4_O5 / length_P_O4) *
				((realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5
					- cosP_O4_O5 * (realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4));
			for (int m = 0; m < 3; m++)//force between P O4 and O5
				fPO4O5_on_O5[m] = fcoe * ((diff_P_O4_O5 / length_P_O5) *
				((realCartesianDirectBondedCoorOwithP[3][m] - cartesianCenterCoorP[m]) / length_P_O4
					- cosP_O4_O5 * (realCartesianDirectBondedCoorOwithP[4][m] - cartesianCenterCoorP[m]) / length_P_O5));
			for (int m = 0; m < 3; m++)//final calculation for force between P O4 and O5 on P
				fPO4O5_on_P[m] = -(fPO4O5_on_O4[m] + fPO4O5_on_O5[m]);

			//store angular force into forceAcceleration(need to know startpoint of P)
			for (int m = 0; m < 3; m++)
				acceleration[i * 3 + m + startPoint[2]] += (fPO1O2_on_P[m] + fPO1O3_on_P[m] + fPO1O4_on_P[m] + fPO1O5_on_P[m]
					+ fPO2O3_on_P[m] + fPO2O4_on_P[m] + fPO2O5_on_P[m] + fPO3O4_on_P[m] + fPO3O5_on_P[m] + fPO4O5_on_P[m]) / massP;


			//store angular force into forceAcceleration(need to know startpoint of O)
			for (int m = 0; m < 3; m++) {
				acceleration[orderObondedwithP[0] * 3 + m + startPoint[0]] += (fPO1O2_on_O1[m] + fPO1O3_on_O1[m] + fPO1O4_on_O1[m] + fPO1O5_on_O1[m]) / massO;
				acceleration[orderObondedwithP[1] * 3 + m + startPoint[0]] += (fPO1O2_on_O2[m] + fPO2O3_on_O2[m] + fPO2O4_on_O2[m] + fPO2O5_on_O2[m]) / massO;
				acceleration[orderObondedwithP[2] * 3 + m + startPoint[0]] += (fPO1O3_on_O3[m] + fPO2O3_on_O3[m] + fPO3O4_on_O3[m] + fPO3O5_on_O3[m]) / massO;
				acceleration[orderObondedwithP[3] * 3 + m + startPoint[0]] += (fPO1O4_on_O4[m] + fPO2O4_on_O4[m] + fPO3O4_on_O4[m] + fPO4O5_on_O4[m]) / massO;
				acceleration[orderObondedwithP[4] * 3 + m + startPoint[0]] += (fPO1O5_on_O5[m] + fPO2O5_on_O5[m] + fPO3O5_on_O5[m] + fPO4O5_on_O5[m]) / massO;
			}
		}
		//--------------------------------------------------------------------------------------------------------------------------------------------------------------

		// calculate repulsive energy and force between P and P
		for (int k = c; k < atomQuantity[2]; ++k)
		{
			//store another P coordiante in Direct(need to know startpoint of P)
			directRepPP[0] = coordinates[k * 3 + startPoint[2]];
			directRepPP[1] = coordinates[k * 3 + 1 + startPoint[2]];
			directRepPP[2] = coordinates[k * 3 + 2 + startPoint[2]];

			//get real direct coor of P near the center P
			EFTool.directShortestCoor(directCenterCoorP, directRepPP, directRepPPshortest);

			//convert direct to cartesian 
			EFTool.directtoCart(directRepPPshortest, cartesianRepPPshortest);

			//get length between two P atoms using Cartesian	
			repPPlength = EFTool.distanceinCart(cartesianCenterCoorP, cartesianRepPPshortest);

			if (repPPlength <= 3.10)
			{
				//calculate energy
				repulsive_energy_SiSiorSiPorPP += 0.5 * repulsiveCoeforSiSiorSiPorPP * 
					(3.10 - repPPlength)*(3.10 - repPPlength)*(3.10 - repPPlength);

				//calculate force
				for (int m = 0; m < 3; ++m)
				{
					frepPP[m] = -3 * repulsiveCoeforSiSiorSiPorPP *
						(3.10 - repPPlength)*(3.10 - repPPlength) / repPPlength *
						(cartesianRepPPshortest[m] - cartesianCenterCoorP[m]);
					//store force in acceleration (need to know startpoint of P)
					acceleration[i * 3 + m + startPoint[2]] += frepPP[m] / massP;
					acceleration[k * 3 + m + startPoint[2]] -= frepPP[m] / massP;
				}
			}
		}
		++c;


	}


	//--------------------------------------------------------------------------------------------------------------

	/*std::cout << "radial energy: " << radial_energySiandO << std::endl;
	std::cout << "Si angular energy: " << angular_energy_Si << std::endl;
	std::cout << "O angular energy: " << angular_energy_O << std::endl;
	std::cout << "Si repulsive energy: " << repulsive_energy_SiSiorSiPorPP << std::endl;
	std::cout << "O repulsive energy: " << repulsive_energy_O_O << std::endl;*/

	total_energy = radial_energySiandO + radial_energyPandO + angular_energy_Si + angular_energy_O +
		repulsive_energy_O_O + repulsive_energy_SiSiorSiPorPP;



};
