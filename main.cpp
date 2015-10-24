// My main file for SUS_3D
// Data: July27 2015
// Authur: Yuding Ai
#include "MC.h"
#include <iostream>
using namespace std;

int main()
{
	stringstream t;
	double start = clock();

	// ======================= Plotting the final config ========================
	vector<HR> VR,HR,UR;
	MC m(1E6L,8,64,64,64,1);
	m.MCSUS();
	// ======================= end of simulation, print out the time =======
	double end = clock();
	cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<" sec."<<endl;
	t <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC) <<" sec."<<endl;
	ofstream myfile ("time.txt");
	string datat = t.str();
	myfile << datat;
	myfile.close();	
	return 0;
}