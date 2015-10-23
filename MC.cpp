/*
* MC.cpp
* Simulation of 2-D Rods By GCMC with SUS, No energy
* Author: Yuding Ai
* Date: 2015.10.19
* *************************** MC implementation ********************************
* This simulation follows Successive Umbrella sampling without energy.
* ******************************************************************************
*
*/


#include "MC.h"

MC::MC(long int ST, int LEN, int N0, int N1, int N2, double Z)
{
	VRodlist.clear(); // the list storage the Vertical Rods;
    HRodlist.clear(); // the list storage the Horizantal Rods;
    URodlist.clear(); // the list storage the Up Rods;
	n0 = N0;          // length of the box
	n1 = N1;          // weight of the box
	n2 = N2;          // hight  of the box
	length = LEN;     // length of rod
	step = ST;        // MC steps
	z = Z;            // activity z = exp(beta*miu)
	nh=nv=nu=dh=dv=du=ah=av=au=0;            // initialize all my number count into 0; ie. nh stands for number of horizontal rods
	                                         // dh stands for the times of deletion for for horizontal rods and ah is for addition.
}

const vector<HR>& MC::getVRodlist() const
{
	return VRodlist;
}

const vector<HR>& MC::getHRodlist() const
{
	return HRodlist;
}

const vector<HR>& MC::getURodlist() const
{
	return URodlist;
}


void MC::Add(Cells &s,double &prob,double &proba)
{
	int x,y,z,o; // pick a random position and orientation for the HR to be added;
	x = rand()%n0; // x range[0,n0-1]
	y = rand()%n1; // y range[0,n1-1]
	z = rand()%n2; // z range[0,n2-1]
	o = rand()%3  ; // 0 range {0,1,2}

	if(s.getSquare(x,y,z).isEmpty()) // if it's open, try to do Addition;
	{
		HR rod(x,y,z,length,o);

		if(o == 0)
		{
		//======================== Vertical ===============================
			// the vertical case
			int counter = 0;

			for (int j = 0; j < length-1; j++)
			{
				// check if the vertical space is wide open
				if(s.getSquare(x,(y+j+1)%n1,z).isOccupied())
				{
					counter++;
					break;
				}						
			}
			// cout <<"counter = "<< counter<< endl;
			if (counter == 0)
			{
				if(prob<=proba)
				{
					// Do addition;
					// push the new rod into the VRodlist;
					VRodlist.push_back(rod);
					av++;
					nv++;// accumulate the # of ver rod;
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{	
						s.getSquare(x,(y+i)%n1,z).setStatus(1);
					}					
				}		
			}									
		}

		else if(o == 1)
		{
        //======================= Horizontal  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare((x+1+j)%n0,y,z).isOccupied())
				{
					counter++;
					break;
				}							
			}
			if (counter == 0)
			{
				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					HRodlist.push_back(rod);
					ah++;
					nh++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare((x+i)%n0,y,z).setStatus(1);
					}
				}
			}
		}
		else 
		{
        //======================= Up  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare(x,y,(z+j+1)%n2).isOccupied())
				{
					counter++;
					break;
				}							
			}
			if (counter == 0)
			{
				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					URodlist.push_back(rod);
					au++;
					nu++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare(x,y,(z+i)%n2).setStatus(1);
					}
				}
			}							
		}
    }
}

void MC::Del(Cells &s,double &prob,double &probd, double &size)
{

	if(nv + nh + nu >0)// make sure there are Vertical rod;
	{
		int indx; // pick a random index from the Rodlist;
		indx = rand()%int(nv+nh+nu);

		//remove Rodlist[indx];
		int x,y,z;// the position of the target on the cells;

		if(prob <= probd)
		{					
			if(indx < nv)
			{
				x = VRodlist[indx].getX();
				y = VRodlist[indx].getY();
				z = VRodlist[indx].getZ();	
				// --------------------- it's a vertical rod -----------------------

				for(int i = 0; i<VRodlist[indx].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,(y+i)%n1,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				VRodlist.erase(VRodlist.begin() + indx);
				nv--;// substract the # of ver rod;
				dv++;
			}
			
			else if (indx < nv + nh)
			{
				x = HRodlist[indx - nv].getX();
				y = HRodlist[indx - nv].getY();
				z = HRodlist[indx - nv].getZ();	
				// --------------------- it's a horizontal rod -----------------------

				for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare((x+i)%n0,y,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				HRodlist.erase(HRodlist.begin() + indx - nv);
				nh--;// substract the # of ver rod;
				dh++;
			}
			else
			{
				x = URodlist[indx - nv - nh].getX();
				y = URodlist[indx - nv - nh].getY();
				z = URodlist[indx - nv - nh].getZ();	
				// --------------------- it's a up rod -----------------------
				for(int i = 0; i<URodlist[indx-nv-nh].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,y,(z+i)%n2).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				URodlist.erase(URodlist.begin()+indx-nv-nh);
				nu--;// substract the # of ver rod;
				du++;

			}
		}
	}	
}


array<double,10000> MC::MCSUS()
{
	Cells s(n0,n1,n2,EMPTY,length);  //  setting the lattice;
	//==========================================================  declare instance variables ============================================================= //
	stringstream sh;
	sh.precision(20);
	double addordel;           // the prob to decide either add or del;
	double probd,proba;      // the acceptance prob of addition and deletion; 
	double prob;               // the prob to decide either accept add/del;
	double aaccp,daccp;      // the acceptance probabilities: 
	double V = double(n0*n1*n2);    // the total lattice size
	double K = double(length); //
    array<double,10000> WF;

	srand(time(NULL));
	long int i = 0; // counter for each window
	double w = 1.0; // a counter that keep track of the index of window
	double fu,fl; // occurrence counter


	//================================Start my MC simulation=================================
	while (w <= 0.8*V/K) // while loop terminate until finish the last window; window[V/K]; !!!rods: 0.8*V/K
	{
		i = 0;  // initialize my step counter; 
		fl = fu = 0; // initialize my occurrence counter;
		while(i < step )
		// Simulation for each window with "step" amount of step
		{
			i++;
			// generate a random probability to decide either add or del;
			addordel = rand()%2;
		    double size = nv+nh+nu;			
			
			prob = ((double) rand() / (RAND_MAX)); 

			aaccp = (z*V)/((size+1.0)*K)*(exp(WF[int(size+1)] - WF[int(size)]));
			daccp = (size*K)/(z*V)*(exp(WF[int(size-1)] - WF[int(size)]));	

			probd = min(1.0,daccp);
			proba = min(1.0,aaccp);

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				if (nv+nh+nu < w  ) // only try to add if the size is below the upper window
				{
					//try to Do Addition;
					Add(s,prob,proba);
				}
			}
			// ============================Deletion=============================
			else 
			{
				if (nv+nh+nu > w - 1) // only try to delete if the size is above the lower window size
				{
					//Do deletion;
					Del(s,prob,probd,size);	
				}			
			}	

			if (nv+nh+nu == w) // only update the fu when addition succeessfully.
			{
				fu++; // if at the upper window, update fu
			}
			else if (nv+nh+nu == w-1 )  // only update the fu when deletion succeessfully.
			{
				fl++;//if at the lower window, update fl
			}	
		}


		// =======================  if fu and fl != 0 Update the upper window side ================================
        if (fu!=0 && fl != 0)
        {
		    WF[w] = WF[w] + log(fu/fl);
	        // "linearly extrapolate" for WF[w+1] by using W[w] and WF[w-1]
	        WF[w+1] = WF[w];

	        cout << fl<<"  "<<fu <<" " <<nv <<"  "<<WF[w+1]<<endl;
			// ======================= Print out the data into terminal =============================================		
			cout <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;
			// initial config determine the intial value of fu and fl
		    w++; // switch into the next window
        }
        // else reset fu and fl and repeat the simulation
	}

	for(int i = 0; i< 0.8*V/K+1; i++) // V/K multiply 0.8 for rods
	{
		sh<<WF[i]<<endl;
	}

	ofstream myfile ("SUSWeight_function.txt");
	myfile.precision(20);
	string data2 = sh.str();
	myfile << data2;
	myfile.close();

	return WF;  
}


void MC::plot(const vector<HR>& VRodlist, const vector<HR>& HRodlist,const vector<HR>& URodlist)
{
	stringstream stv,sth,stu;

	for (int i = 0; i < VRodlist.size(); i++)
	{
		stv<< VRodlist[i].getX() << "   "<< VRodlist[i].getY() << "   "<< VRodlist[i].getZ()<<endl;
	}

	ofstream myfilev ("3dplotv.txt");
	string datav = stv.str();
	myfilev << datav;
	myfilev.close();

	for (int j = 0; j < HRodlist.size(); j++)
	{
		sth<< HRodlist[j].getX() << "   "<< HRodlist[j].getY() << "   "<< HRodlist[j].getZ()<<endl;
	}

	ofstream myfileh ("3dploth.txt");
	string datah = sth.str();
	myfileh << datah;
	myfileh.close();

	for (int k = 0; k < URodlist.size(); k++)
	{
		stu<< URodlist[k].getX() << "   "<< URodlist[k].getY() << "   "<< URodlist[k].getZ()<<endl;
	}

	ofstream myfileu ("3dplotu.txt");
	string datau = stu.str();
	myfileu << datau;
	myfileu.close();

}





