// cells.cpp
// 3-D Rods
// Author: Yuding Ai
// Date: 2015.07.15


#include "cells.h"
extern const string EXC_INVALID_DESC = "Not a valid description of cells!";

Cells::Cells()
{
	Boxlist.clear();
	Planelist.clear();
	N0 = 1;
	N1 = 1;
	N2 = 1;
	size = 1;
	//initialize my cell
	arr = new Square[1];
	arr[0] = Square();
}


Cells::Cells(int X, int Y, int Z, int init,int length)
{
	Boxlist.clear(); // the list of subBox
	Planelist.clear(); // the list of subplane
	N0 = X; // length
	N1 = Y; // width
    N2 = Z; // height
    size = (N0)*(N1)*(N2);
    arr = new Square[size];
    if (init == EMPTY)
    {
    	//initialize my cell with empty config
		for(int i = 0; i < size; i++)
		{
			arr[i] = Square();
		}
    }

    else if (init == BOX)
	{
		//initialize my cell with filling up with subboxes config

    	//firstly: initialize my cell with full config
		for(int i = 0; i < size; i++)
		{
			arr[i] = Square(1);
		}
		//then: randomly filling up it with subboxes
		srand (time(NULL));
    	for(int i = 0; i < X/length; i++)
    	{   
			for ( int j = 0; j < Y/length; j++)
			{
				for (int k = 0; k < Z/length; k++)
				{
					int Or; //pick a random config of subbox to add into Box
					Or = rand()%3; 
					int l,m,n;
					l = i*length;
					m = j*length;
					n = k*length;
					Boxgen subox(l,m,n,Or,length);
					Boxlist.push_back(subox);
				}
			}
    	}
    }

    else if (init == PLANE)
    {
    	//initialize my cell with filling up with subplane config
    	//firstly: initialize my cell with full config
		for(int i = 0; i < size; i++)
		{
			arr[i] = Square(1);
		}

		//then: randomly filling up it with subplanes
		srand (time(NULL));
		int Pl; // pick an plane first
		Pl = rand()%3;

		if (Pl == X_)
		{
			int Or;// Pick a random Orientation: in the case of X_plane, one can only choose from VER or UP
			for (int i = 0; i < N0; i ++)
			{				
				Or = 2*(rand()%2);// Or can be either 0 --> Ver or 2 --> UP;				
				Planegen subplane(X_, i, N0, Or, length);
				// Planegen(int Plane, int Lay, int L, int Orien,int Len); 
				Planelist.push_back(subplane);
			}

		}

		else if (Pl == Y_)
		{
			int Or;// Pick a random Orientation: in the case of Y_plane, one can only choose from HOR or UP
			for (int i = 0; i< N1; i++)
			{
				Or = rand()%2 +1;// Or can be either 1 --> HOR or 2 --> UP;
				Planegen subplane(Y_, i, N1, Or, length);
				Planelist.push_back(subplane);
			}			
		}

		else if (Pl == Z_)
		{
			int Or;// Pick a random Orientation: in the case of Z_plane, one can only choose from VER or HOR
			for (int i = 0; i<N2; i++)
			{
				Or = rand()%2; // Or can be either 0 --> VER or 1 --> HOR;
				Planegen subplane(Z_, i, N2, Or, length);
				Planelist.push_back(subplane);
			}
		}
    }
}

/*
*  Destructor
*/
Cells::~Cells()
{   
	delete [] arr;	
}

// *** Getters *** //
int Cells::getN0() const
{
	return N0;
}
int Cells::getN1() const
{
	return N1;
}

int Cells::getN2() const
{
	return N2;
}
int Cells::getSize() const
{
	return (N0)*(N1)*(N2);
}

const vector<Boxgen>& Cells::getBoxlist() const
{
	return Boxlist;
}

const vector<Planegen>& Cells::getPlanelist() const
{
	return Planelist;
}
// *** Other Functionality *** //

int Cells::getIdx( int x, int y, int z) const
{
	return  x + N0*y + N1*N0*z;
}

Square& Cells::getSquare( int x, int y, int z) const
{
	if (x >= N0 || y >= N1 || z >= N2)
	{
		throw EXC_INVALID_DESC;
	}

	int idx = getIdx(x,y,z);
	return arr[idx];
}

