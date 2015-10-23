// cells.h	
// 3-D Rods
// Author: Yuding Ai
// Date: 2015.06.05

#ifndef CELLS_H
#define CELLS_H
#include <vector>
#include "hardrods.h"
#include "square.h"
#include "Boxgen.h"
#include "Planegen.h"
using namespace std;

// template< int N0, int N1, int N2> inline const int rten_idx(const int i0, const int i1, const int i2)
// {
// 	return i2 + N2*(i1 + N1*i0); 
// }

// template< int N0, int N1, int N2> inline const int rten_idx(const int i0, const int i1, const int i2)
// {
// 	return i2 + N2*(i1 + N1*i0); 
// }

class Cells
{
	private:
		int N0;         // Weith; 
		int N1;         // Length; 
		int N2;         // Height;
		int size;       // the total size of the grid
		Square * arr;   // linearized array of Squares
		std::vector<Boxgen> Boxlist; // the list storage the subboxes;
		std::vector<Planegen> Planelist; // the list storage the subboxes;


	public:
		enum INIT {EMPTY,BOX,PLANE};
		enum PLANE_DIRECTION {X_,Y_,Z_};
		enum BOX_DIRECTION {VER,HOR,UP};
		Cells();

		/*
		*  Constructor builds a cells from an input stream
		*  @param      int r - row index
		*  @param      int c - column index
		*  initialize all states to 0;
		*/
		Cells(int X,int Y, int Z, int init,int length);


		/*
		*  Destructor
		*/
		~Cells();

		// *** Getters *** //
		int getN0() const;
		int getN1() const;
		int getN2() const;
		int getSize() const;
		int getIdx( int x, int y, int z) const;
		const vector<Boxgen>& getBoxlist() const;
		const vector<Planegen>& getPlanelist() const;



		// *** Other Functionality *** //
		/*
		*  @function      getSquare
		*  @param         int i - row index, in range 0 through numRows - 1
		*  @param         int j - column index, in range 0 through numCols - 1
		*  @return        Square& - return a reference to the (i,j)-Square stored
		*                 in this cells
		*/
		Square& getSquare(int x, int y,int z) const;

	};


#endif /* CELLS_H */