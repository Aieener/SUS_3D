/*
* hardrods.cpp
* 3-D Rods
* Author: Yuding Ai 
* Date: 2015.06.08
*/

#include "hardrods.h"
#include <iostream>
#include <vector>
#include <string>

HR::HR()
{
	X = 0;
	Y= 0;
	Z = 0;
	Length = 1;
	Orientation = 0;
}

HR::HR(int x, int y, int z, int l, int o)
{
	X = x;
	Y= y;
	Z = z;
	Length = l;
	Orientation = o;
}

int HR::getX() const
{
	return X;
}

int HR::getY() const
{
	return Y;
}

int HR::getZ() const
{
	return Z;
}

int HR::getLength() const
{
	return Length;
}

int HR::getOrientation() const
{
	return Orientation;
}

