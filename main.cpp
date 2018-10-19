//============================================================================
// Name        : PSG
// Author      : Zhi Song
// Version     : 1.0.0
// Copyright   : Your copyright notice
// Description : Create code to generate PSG system base on some
//				 crystal structure in C++
//============================================================================

#include <iostream>
#include "CONSTANT.h"
#include "IO.h"
#include "Group.h"	
#include "EnergyForce.h"
#include "BFGS.h"
#include "Tool.h"
#include <stdlib.h> 


int main() {

	

	int doubleO[2][2] = { {99,104}, {4,4} };

	//declare a struct contains energy and MPI rank
	struct{
		T energy;
		int rank;
	}totalEnergy, lowEnergy;

	

	//initialize totalEnergy
	totalEnergy.energy = 0;
	totalEnergy.rank = 0;
	lowEnergy.energy = 0;
	lowEnergy.rank = 0;

	//read file
	IO mainIO;
	mainIO.readfromFile();


	//first group atom
	Group mainGroup;
	mainGroup.Graph();

	graph[atomQuantity[1]][doubleO[0][0]] = 2;
	graph[atomQuantity[1]+1][doubleO[0][1]] = 2;

	mainGroup.groupfromGraph(doubleO);

	std::cout << "P1: " << groupP[0][0] << "     " << groupP[0][1] << "     " << groupP[0][2] << "     "
		<< groupP[0][3] << "     " << groupP[0][4] << "     " << std::endl;
	std::cout << "P2: " << groupP[1][0] << "     " << groupP[1][1] << "     " << groupP[1][2] << "     "
		<< groupP[1][3] << "     " << groupP[1][4] << "     " << std::endl;

	//call constructers 
	EnergyForce mainEF;
	Tool mainTool;

	//set up initial old energy
	mainEF.obtainEnergyForce(atomCoordinate, atomAcceleration, oldEnergy, doubleO);
	std::cout << "Old Energy: "<<oldEnergy <<std::endl;	

	//go into the cutting bonds process
	for (int i = 0; i < cutBondsTimes; ++i) {

		//store parameters
		mainTool.copyGraph(graph, oldGraph, atomQuantity[1]);//store graph to old graph

		//cut bonds using graph
		mainGroup.cutBond(graph, doubleO);

		//group from graph
		mainGroup.groupfromGraph(doubleO);

		std::cout << "P1: " << groupP[0][0] << "     " << groupP[0][1] << "     " << groupP[0][2] << "     "
			<< groupP[0][3] << "     " << groupP[0][4] << "     "<<std::endl;
		std::cout << "P2: " << groupP[1][0] << "     " << groupP[1][1] << "     " << groupP[1][2] << "     "
			<< groupP[1][3] << "     " << groupP[1][4] << "     " << std::endl;

		mainIO.output(1);
		//do linesearch (BFGS)
		BFGS mainBFGS;
		mainBFGS.LineSearch(totalEnergy.energy, doubleO);
		std::cout << totalEnergy.energy << std::endl;
		mainIO.output(0);

		//MC determination 
		//auto p = mainTool.MC_probability(oldEnergy, totalEnergy.energy);

		//std::cout<< "Old energy: " << oldEnergy << "    New Energy: " << totalEnergy.energy << std::endl;

		//if p==0, energy equals infinteNumber which will not be accepted
		//if (p == 0)
			//totalEnergy.energy = infiniteNumber;
	}
	system("Pause");
	return 0;
}



