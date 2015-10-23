// Boxgen.cpp
// initial box config generator
// Author: Yuding Ai
// Data : July, 23, 2015

#include "Boxgen.h"

Boxgen::Boxgen(int Lx,int Ly, int Lz,int Orien,int Len)
{
	lx = Lx;
	ly = Ly;
	lz = Lz;
	o = Orien;
	length = Len;
	int x,y,z;

	if (o == VER) // Ver
	{
		for ( int i = 0; i < Len; i++)
		{
			for (int j = 0; j < Len; j++)
			{
				x=lx + j;
				y=ly + 0;
				z=lz + i;
				HR rv(x,y,z,length,o);
				BVRodlist.push_back(rv);
			}
		}
	    BHRodlist.clear(); // the list storage the Horizantal Rods;
	    BURodlist.clear(); // the list storage the Up Rods;
	}
	else if (o == HOR) // Hor
	{
		for ( int i = 0; i < Len; i++)
		{
			for (int j = 0; j < Len; j++)
			{
				x=lx + 0;
				y=ly + j;
				z=lz + i;
				HR rh(x,y,z,length,o);
				BHRodlist.push_back(rh);
			}
		}
		BVRodlist.clear(); // the list storage the Vertical Rods;
	    BURodlist.clear(); // the list storage the Up Rods;
	}
	else if (o == UP) // Up
	{
		for ( int i = 0; i < Len; i++)
		{
			for (int j = 0; j < Len; j++)
			{
				x=lx + i;
				y=ly + j;
				z=lz + 0;
				HR ru(x,y,z,length,o);
				BURodlist.push_back(ru);
			}
		}
		BVRodlist.clear(); // the list storage the Vertical Rods;
	    BHRodlist.clear(); // the list storage the Horizantal Rods;
	}

	else
	{
		BVRodlist.clear(); // the list storage the Vertical Rods;
	    BHRodlist.clear(); // the list storage the Horizantal Rods;
	    BURodlist.clear(); // the list storage the Up Rods;
	}
}

int Boxgen::getLength() const
{
	return length;
}

int Boxgen::getOri() const
{
	return o;
}


array<int,3> Boxgen::getPos() const
{
	array<int,3> bar;
	bar[0] =lx;
	bar[1] =ly;
	bar[2] =lz;
	return bar;
}

const vector<HR>& Boxgen::getBVRodlist() const
{
	return BVRodlist;
}

const vector<HR>& Boxgen::getBHRodlist() const
{
	return BHRodlist;
}

const vector<HR>& Boxgen::getBURodlist() const
{
	return BURodlist;
}