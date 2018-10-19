#pragma once
#include <fstream>
#include <iostream>
#include "CONSTANT.h"



class IO
{
public:
	IO();
	~IO();

public:
	void readfromFile();
	void output(const int &n);
};

