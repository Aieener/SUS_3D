// Planegen.cpp
// initial Plane config generator
// Author: Yuding Ai
// Data : July, 27, 2015

#include "Planegen.h"

Planegen::Planegen(int Plane, int Lay, int L, int Orien,int Len)
{
	plane = Plane;
    l = L;
    lay = Lay;
	o = Orien;
	length = Len;

	int x,y,z;
	if (o == 0) // VerRod, lying in the "y" direction thus can only construct for "z" plane and "x" plane
	{
		for ( int i = 0; i < l; i++)
		{
			for (int j = 0; j < l/length; j++)
			{
				if(plane == Z_)
				{
					x=i;
					y=j*length;
					z=lay;
					HR rv(x,y,z,length,o);
					PVRodlist.push_back(rv);
				}

				else if (plane == X_)
				{
					x=lay;
					y=j*length;
					z=i;
					HR rv(x,y,z,length,o);
					PVRodlist.push_back(rv);
				}
			}
		}
	    PHRodlist.clear(); // the list storage the Horizantal Rods;
	    PURodlist.clear(); // the list storage the Up Rods;
	}
	else if (o == 1) // HorRod, lying in the "x" direction thus can only construct for "y" plane and "z" plane
	{
		for ( int i = 0; i < l; i++)
		{
			for (int j = 0; j < l/length; j++)
			{
				if(plane == Y_)
				{
					x=j*length;
					y=lay;
					z=i;
					HR rh(x,y,z,length,o);
					PVRodlist.push_back(rh);
				}

				else if (plane == Z_)
				{
					x=j*length;
					y=i;
					z=lay;
					HR rh(x,y,z,length,o);
					PHRodlist.push_back(rh);
				}
			}
		}
		PVRodlist.clear(); // the list storage the Vertical Rods;
	    PURodlist.clear(); // the list storage the Up Rods;
	}
	else if (o == 2) // UpRod,  lying in the "z" direction thus can only construct for "x" plane and "y" plane
	{
		for ( int i = 0; i < l; i++)
		{
			for (int j = 0; j < l/length; j++)
			{
				if(plane == X_)
				{
					x=lay;
					y=i;
					z=j*length;
					HR ru(x,y,z,length,o);
					PURodlist.push_back(ru);
				}

				else if (plane == Y_)
				{
					x=i;
					y=lay;
					z=j*length;
					HR ru(x,y,z,length,o);
					PURodlist.push_back(ru);
				}
			}
		}
		PVRodlist.clear(); // the list storage the Vertical Rods;
	    PHRodlist.clear(); // the list storage the Horizantal Rods;
	}

	else
	{
		PVRodlist.clear(); // the list storage the Vertical Rods;
	    PHRodlist.clear(); // the list storage the Horizantal Rods;
	    PURodlist.clear(); // the list storage the Up Rods;
	}
}

int Planegen::getLength() const
{
	return length;
}

int Planegen::getOri() const
{
	return o;
}

const vector<HR>& Planegen::getPVRodlist() const
{
	return PVRodlist;
}

const vector<HR>& Planegen::getPHRodlist() const
{
	return PHRodlist;
}

const vector<HR>& Planegen::getPURodlist() const
{
	return PURodlist;
}

