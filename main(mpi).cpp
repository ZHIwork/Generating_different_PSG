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
#include "mpi.h"


int main(int argc, char** argv) {


	int acceptedTimes = 0;
	int doubleO[2][2] = { {-1,-1}, {5,4} };
	

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
	mainGroup.groupfromGraph(doubleO);


	//call constructers 
	EnergyForce mainEF;
	Tool mainTool;

	//set up initial old energy
	mainEF.obtainEnergyForce(atomCoordinate, atomAcceleration, oldEnergy, doubleO);
	std::cout << "Old Energy: "<<oldEnergy <<std::endl;	

	//declare MPI parameters
	int size, rank;
    MPI_Status status;

    //initialize MPI
    MPI_Init(NULL, NULL);
    // Get the number of processes
    size = MPI::COMM_WORLD.Get_size();
    // Get the rank of the processes
    rank = MPI::COMM_WORLD.Get_rank();


	//go into the cutting bonds process
	for (int i = 0; i < cutBondsTimes; ++i) {

		//auto start = chrono::steady_clock::now();
		
		//=====================================================================================================================
		if (rank == 0)//master rank 0
		{
			//cout<<"Times: "<<i<<"  Total energy: "<<total_energy.energy<<" lowenergy: "<<lowenergy.energy<<"\n";
			totalEnergy.energy = infiniteNumber;//set an infinite large number here
			totalEnergy.rank = rank;//store master number into totalEnergy.rank
		}
		else //slave rank
		{

			totalEnergy.rank = rank;//store each slave number into totalEnergy.rank

			//store parameters
			mainTool.copyGraph(graph, oldGraph, atomQuantity[1]);//store graph to old graph
			

			//cut bonds using graph
			mainGroup.cutBond(graph, doubleO);

			//group from graph
			mainGroup.groupfromGraph(doubleO);
			
			mainIO.output(1);
			//do linesearch (BFGS)
			BFGS mainBFGS;
			mainBFGS.LineSearch(totalEnergy.energy, doubleO);
			std::cout << totalEnergy.energy <<std::endl;
			//mainIO.output(rank);
			
			//MC determination 
			//auto p = mainTool.MC_probability(oldEnergy, totalEnergy.energy);

			//std::cout<< "Old energy: " << oldEnergy << "    New Energy: " << totalEnergy.energy << std::endl;

			//if p==0, energy equals infinteNumber which will not be accepted
			//if (p == 0)
				//totalEnergy.energy = infiniteNumber;
		}

		//=====================================================================================================================
		//obtain all total energy from each node and find the lowest one and its rank
		MPI_Allreduce(&totalEnergy, &lowEnergy, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

		//=====================================================================================================================
		//if rank number = rank number of lowest energy and this rank is not master rank 0,
		//send all coordinates and graph to master rank
		if (rank == lowEnergy.rank && rank != 0)
		{
		
			//std::cout<<"Send data" << "   "<<rank<<std::endl;

			//send graph to master rank, tag=0
			MPI::COMM_WORLD.Send(&graph, 64 * 129, MPI::INT, 0, 0);

			//send coordinates to master rank, tag=1
			MPI::COMM_WORLD.Send(&atomCoordinate, totalAtomQuantity * 3, MPI::DOUBLE, 0, 1);
		}

		//auto end = chrono::steady_clock::now();
		//auto diff = end - start;
		//cout << chrono::duration <double, milli>(diff).count() / 1000.0 << " s" << endl;

		//=====================================================================================================================
		//when in the master rank 0 and rank number of lowest energy is not from master rank,
		//master rank receive all coordiantes and graph from slave ranks
		if (rank == 0 && lowEnergy.rank != 0) 
		{
			//receive graph, tag = 0
			MPI::COMM_WORLD.Recv(&graph, 64 * 129, MPI::INT, MPI::ANY_SOURCE, 0);

			//receive coordiates, tag = 1	
			MPI::COMM_WORLD.Recv(&atomCoordinate, totalAtomQuantity * 3, MPI::DOUBLE, MPI::ANY_SOURCE, 1);
	
			//store lowest energy to oldEnergy
			oldEnergy = lowEnergy.energy;

			std::cout << "Accpeted Energy: " << oldEnergy << std::endl;

			//accept new structure and increment accept times 
			++acceptedTimes;

			mainIO.output(acceptedTimes);
		}
		
	
		MPI::COMM_WORLD.Bcast(&graph, 64 * 129, MPI::DOUBLE, 0);
        	MPI::COMM_WORLD.Bcast(&atomCoordinate, totalAtomQuantity * 3, MPI::DOUBLE, 0);
		MPI::COMM_WORLD.Bcast(&oldEnergy, 1, MPI::DOUBLE, 0);
    	}

    	MPI::COMM_WORLD.Barrier();
    	if(rank == 0)
    	{
			//output file
			mainIO.output(0);
           	//std::cout << "Accepted times are: "<<acceptedTimes<<std::endl;
    	}

    	//Finalize the MPI environment. No more MPI calls can be made after this
    	MPI::Finalize();


	return 0;
}



