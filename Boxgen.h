// Boxgen.h
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "square.h"
#include "cells.h"
#include "hardrods.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <deque>
#include <array>
using namespace std;


// vector1.insert( vector1.end(), vector2.begin(), vector2.end() );


#ifndef BOXGEN_H
#define BOXGEN_H
class Boxgen
{
    private:
    	//data members;
        std::vector<HR> BVRodlist; // the list storage the Vertical Rods;
        std::vector<HR> BHRodlist; // the list storage the Horizantal Rods;
        std::vector<HR> BURodlist; // the list storage the Up Rods;
        int o; // the orientation of the box 
        /*
        0 ---------------> VER filling up with Verticle Rods;
        1 ---------------> HOR filling up with Horizontal Rods;
        2 ---------------> UP  filling up with Up Rods;
        3 and other nonsense # ---------------> empty box;
        */
    	int length; // length of the rod
        int lx,ly,lz; // the location of the subbox to be assign to
       
    public:
        //Constructer: initialize the oritation of our box
        //and then assighn the corresponding rod pos into a list.
        //when we add a subbox into our Box/Cells, we assign the rod list into the real Rodlist;
        enum BOX_DIRECTION {VER,HOR,UP};

    	Boxgen(int Lx,int Ly, int Lz, int Orien,int Len); 

    	// ********* Getters********//
        int getLength() const;
        array<int,3> getPos() const;
        int getOri() const;

        const vector<HR>& getBVRodlist() const;
        const vector<HR>& getBHRodlist() const;
        const vector<HR>& getBURodlist() const;

        // ******** Setters ******//

    	// ******** Other Functianality *******//

};

#endif /* BOXGEN_H */