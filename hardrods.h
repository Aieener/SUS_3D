/*
* hardrods.h
* 3-D  rod
* Author: Yuding Ai 
* Date: 2015.06.08
*/
#ifndef HARDRODS_H
#define HARDRODS_H
#include <iostream>
#include <vector>
#include <string>

class HR
{
	private:
		int X; // the x index;
		int Y; // the y index;
		int Z; // the z index;
		int Orientation; // the orientation of the rod; 3 options:  ver,hor, up
		int Length; // the length of the rod

	public:
		//class constants;
		enum ORI {VER,HOR,UP};

		/*
		* Default constructor
		* @para     none
		* @set X = Y = Orientation = 0; Length = 1
		*/
		HR();

		/*
		* constructor
		* @para     int x, y, L, orientation
		*/
		HR(int x, int y, int z, int l, int o);

		//**** Getters *****//
		int getX() const;
		int getY() const;
		int getZ() const;
		int getLength() const;
		int getOrientation() const;
};

#endif /* HARDRODS_H */
