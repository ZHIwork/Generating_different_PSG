#include "IO.h"



IO::IO()
{
	//std::cout << "IO constructor" << std::endl;
}


IO::~IO()
{
	//std::cout << "IO deconstructor" << std::endl;
}


void IO::readfromFile()
{	
	std::ifstream inFile;
	inFile.open(importFileName);
	if (!inFile.is_open()) // failed to open file
	{
		std::cout << "Could not open the file " << importFileName << std::endl;
	}
	else {
		std::string temp;//store useless line
		T lc;//get lattice constant from file
		T lv;//get each lattice vector coordinate from file
		std::string type;//get types of atoms 
		//int Quantity;//get quantity of each type atoms
		T Coor;//get coordinates from file
		for(int i=0; i<totalAtomQuantity+8; i++)
		{
			if (i == 1)//Second line of POSCAR
			{
				inFile >> lc;
				latticeConstant = lc;
			}
			else if (i > 1 && i <= 4)//Third to Fifth lines of POSCAR
			{
				inFile >> lv;
				latticeVector[(i - 2) * 3] = lv;
				inFile >> lv;
				latticeVector[(i - 2) * 3 + 1] = lv;
				inFile >> lv;
				latticeVector[(i - 2) * 3 + 2] = lv;

			}
			else if (i == 5) //Sixth line of POSCAR
			{
				for (int a = 0; a < atomType; a++)
				{
					inFile >> type;
					atomName[a] = type;
				}
			}
			/*else if (i == 6) //Seventh line of POSCAR
			{
				for( int a=0; a<atomType; a++)
				{
					inFile >> _atomQuantity;
					atomQuantity[a] = _atomQuantity;
				}
			}*/
			else if (i >= 8) //Nineth line to the end 
			{
				inFile >> Coor;
				atomCoordinate[(i - 8) * 3] = Coor;
				inFile >> Coor;
				atomCoordinate[(i - 8) * 3 + 1] = Coor;
				inFile >> Coor;
				atomCoordinate[(i - 8) * 3 + 2] = Coor;
			}
			getline(inFile, temp);
		};

		//get inverse lattice vector(need to be revised Not use mannul steps)
		/*for (int i = 0; i < 3*3; ++i)
		{
			inverseLatticeVector[i] = 0.0;
		}
		for (int i = 0; i < 3; ++i)
		{
			inverseLatticeVector[i * 3 + i] = 0.06977393246;
		}*/
	};

	inFile.close();
};



void IO::output(const int &n)
{
	std::ifstream inFile;
	inFile.open(importFileName);
	if (!inFile.is_open()) // failed to open file
	{
		std::cout << "Could not open the file " << importFileName << std::endl;
	}
	else
	{
		std::string temp;//store line
		std::string num = std::to_string(n);
		std::ofstream Savefile(exportFileName + num);//create a file to output data
		Savefile.setf(std::ios::fixed, std::ios::floatfield);
		Savefile.precision(12);
		for (int i = 0; i < 8; i++)
		{
			getline(inFile, temp);
			Savefile << temp;
			if (i != 7)
				Savefile << std::endl;
		}
		Savefile << std::endl;
		Savefile.clear();

		//output position
		for (int i = 0; i < totalAtomQuantity * 3; i=i+3)
		{

			Savefile << atomCoordinate[i];
			Savefile << "    ";
			Savefile << atomCoordinate[i + 1];
			Savefile << "    ";
			Savefile << atomCoordinate[i + 2];
			Savefile << "    ";
			Savefile << std::endl;
			
		}
	}
}
