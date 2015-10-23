// Planegen.h
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


#ifndef PLANEGEN_H
#define PLANEGEN_H
class Planegen
{
    private:
    	//data members;
        std::vector<HR> PVRodlist; // the list storage the Vertical Rods;
        std::vector<HR> PHRodlist; // the list storage the Horizantal Rods;
        std::vector<HR> PURodlist; // the list storage the Up Rods;
        int o; // the orientation of the plane 
        /*
        0 ---------------> filling up with Verticle Rods;
        1 ---------------> filling up with Horizontal Rods;
        2 ---------------> filling up with Up Rods;
        3 and other nonsense # ---------------> empty box;
        */
    	int length; // length of the rod
        int plane; // the plane idex:
        /*
         plane is a string index indicates with plane: 
         x ------------- y,z plane
         y ------------- z,x plane
         z ------------- x,y plane
        */
        int lay,l; // the # of the layer of the plane and the size of the BOX: here we assume a square box is lxlxl 

    public:
        //Constructer: initialize the oritation of our box
        //and then assighn the corresponding rod pos into a list.
        //when we add a subbox into our Box/Cells, we assign the rod list into the real Rodlist;
        enum PLANE_DIRECTION {X_,Y_,Z_};

    	Planegen(int Plane, int Lay, int L, int Orien,int Len); 

    	// ********* Getters********//
        int getLength() const;
        int getOri() const;

        const vector<HR>& getPVRodlist() const;
        const vector<HR>& getPHRodlist() const;
        const vector<HR>& getPURodlist() const;

        // ******** Setters ******//

    	// ******** Other Functianality *******//

};

#endif /* PLANEGEN_H */