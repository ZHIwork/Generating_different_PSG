#pragma once

#include <iostream>
#include <random>
#include <cmath>
#include "CONSTANT.h"
#include "mkl.h"

class Tool
{
public:
	Tool();
	~Tool();

public:
	void directtoCart(const T(&direct)[3], T(&cart)[3]);
	void carttoDirect(const T(&cart)[3], T(&direct)[3]);
	T distanceinCart(const T(&a)[3], const T(&b)[3]);
	T angleinCart(const T(&center)[3], const T(&bonded1)[3], const T(&bonded2)[3], const T &c1, const T &c2);
	void directShortestCoor(const T(&center)[3], const T(&comp)[3], T(&shortest)[3]);
	bool MC_probability(const T &old_energy, const T &new_energy);
	T generator_lb(T begin, T end);
	void copyGraph(const int(&a)[][129], int(&b)[][129], const int &row);
	void copyCoordinates(const T a[], T b[], const int &N);
};
