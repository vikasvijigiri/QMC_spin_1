#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <filesystem>


int lx, ly, Ns, Nm, Np, Nb[3], np1, S;
int  nbins, isteps, iter;
long n1, Lc;
int **JHsites, **JQsites, **JBsites, *JBsgnx, *JBsgny, *JHsgnx, *JHsgny, *JQsgnx, *JQsgny;
double J[3];
double Beta;
int pf, tf, pf1, tf1;

void
SSEvariables::declare_variables()
{
	// Read input variables from file.
	//std::cout << "path is   " << std::filesystem::current_path() << std::endl;
	//std::cout <<  std::filesystem::exists("input_param.dat") << std::endl;
	std::string INPUT_FILE = "input_param.dat";
	std::ifstream vars(INPUT_FILE);
	vars >> lx;
	vars >> J[0]; 			// Heisenberg AFM interaction strength
	vars >> J[1]; 
	vars >> J[2]; 			// Biquad AFM interaction strength	
	vars >> Beta;
	vars >> iter;
	vars >> isteps;
	S = 1;
	
	
	// Initialize global variables, spin = 1;
	ly = lx;			
	Ns = lx * ly;				// No. of Fundamental spins.
	Np = Ns;						// No. of projector bonds.
	Nm = 2 * S * Ns;		// No. of Mini-spins
			
	if (lx == 2) {
		// For single plaquette on square lattice.
		Nb[0] = 16;				// H
		Nb[1] = 8;				// QQ
		Nb[2] = 8;				// Biquad	
	}
	else {	
		// For PBC square lattice. Minimum size is 3x3	
		Nb[0] = 8*Ns;			// No. of Heisenberg bonds
		Nb[1] = 8*Ns; 		// No. of Q_Q bonds
		Nb[2] = 4*Ns;			// No. of biquad bonds	
	}
}


void SSEvariables::lattice_sites()
{

	
	pf = Nb[0]/2, tf = 4, pf1 = Nb[2]/2, tf1 = 2;

	if (lx == 2){
		// Chain lattice for 2 x 2 square lattice (plaquette)
		
		// QQ bonds.
		JQsgnx = new int [Nb[1]];
		JQsgny = new int [Nb[1]];
		JQsites = new int *[Nb[1]];		
		for (int i = 0; i < Nb[1]; i++) JQsites[i] = new int[8];
		JQsites[0][0] = 0;
		JQsites[0][1] = 1;
		JQsites[0][2] = 4;
		JQsites[0][3] = 5;		// aa
		JQsites[0][4] = 2;
		JQsites[0][5] = 3;
		JQsites[0][6] = 6;
		JQsites[0][7] = 7;
		JQsgnx[0] = 1;
		JQsgny[0] = 0; 		
		//
		JQsites[1][0] = 0;
		JQsites[1][1] = 5;
		JQsites[1][2] = 1;
		JQsites[1][3] = 4;		// ab
		JQsites[1][4] = 2;
		JQsites[1][5] = 3;
		JQsites[1][6] = 6;
		JQsites[1][7] = 7;
		JQsgnx[1] = 1;
		JQsgny[1] = 0; 			
		//
		JQsites[2][0] = 0;
		JQsites[2][1] = 1;
		JQsites[2][2] = 4;
		JQsites[2][3] = 5;		// ba
		JQsites[2][4] = 2;
		JQsites[2][5] = 7;
		JQsites[2][6] = 3;
		JQsites[2][7] = 6;	
		JQsgnx[2] = 1;
		JQsgny[2] = 0; 			
		//
		JQsites[3][0] = 0;
		JQsites[3][1] = 5;
		JQsites[3][2] = 1;
		JQsites[3][3] = 4;		// bb
		JQsites[3][4] = 2;
		JQsites[3][5] = 7;
		JQsites[3][6] = 3;
		JQsites[3][7] = 6;
		JQsgnx[3] = 1;
		JQsgny[3] = 0; 			
		//***************
		JQsites[4][0] = 0;
		JQsites[4][1] = 2;
		JQsites[4][2] = 4;
		JQsites[4][3] = 6;		// aa
		JQsites[4][4] = 1;
		JQsites[4][5] = 3;
		JQsites[4][6] = 5;
		JQsites[4][7] = 7;
		JQsgnx[4] = 0;
		JQsgny[4] = 1; 			
		//
		JQsites[5][0] = 0;
		JQsites[5][1] = 2;
		JQsites[5][2] = 4;
		JQsites[5][3] = 6;		// ab
		JQsites[5][4] = 1;
		JQsites[5][5] = 7;
		JQsites[5][6] = 3;
		JQsites[5][7] = 5;
		JQsgnx[5] = 0;
		JQsgny[5] = 1; 			
		//
		JQsites[6][0] = 0;
		JQsites[6][1] = 6;
		JQsites[6][2] = 2;
		JQsites[6][3] = 4;		// ba
		JQsites[6][4] = 1;
		JQsites[6][5] = 3;
		JQsites[6][6] = 5;
		JQsites[6][7] = 7;	
		JQsgnx[6] = 0;
		JQsgny[6] = 1; 				
		//
		JQsites[7][0] = 0;
		JQsites[7][1] = 6;
		JQsites[7][2] = 2;
		JQsites[7][3] = 4;		// bb
		JQsites[7][4] = 1;
		JQsites[7][5] = 7;
		JQsites[7][6] = 3;
		JQsites[7][7] = 5;
		JQsgnx[7] = 0;
		JQsgny[7] = 1; 					

		JHsgnx = new int [Nb[0]];
		JHsgny = new int [Nb[0]];
		JHsites = new int *[Nb[0]];
		for (int i = 0; i < Nb[0]; ++i) JHsites[i] = new int[2];

		JHsites[0][0] = 0;		// H-bonds,  horizontal above
		JHsites[0][1] = 1;
		JHsites[1][0] = 0;		// H-bonds,  right-aligned
		JHsites[1][1] = 5;
		JHsites[2][0] = 1;		// H-bonds,  left-aligned
		JHsites[2][1] = 4;
		JHsites[3][0] = 4;		// H-bonds,  horizontal below
		JHsites[3][1] = 5;	
		JHsgnx[0] = 1;
		JHsgnx[1] = 1;
		JHsgnx[2] = 1;
		JHsgnx[3] = 1; 
		JHsgny[0] = 0;
		JHsgny[1] = 0;
		JHsgny[2] = 0;
		JHsgny[3] = 0; 										

		JHsites[4][0] = 2;		// H-bonds,  horizontal above
		JHsites[4][1] = 3;
		JHsites[5][0] = 3;		// H-bonds,  right-aligned
		JHsites[5][1] = 6;
		JHsites[6][0] = 2;		// H-bonds,  left-aligned
		JHsites[6][1] = 7;
		JHsites[7][0] = 6;		// H-bonds,  horizontal below
		JHsites[7][1] = 7;			
		JHsgnx[4] = 1;
		JHsgnx[5] = 1;
		JHsgnx[6] = 1;
		JHsgnx[7] = 1; 
		JHsgny[4] = 0;
		JHsgny[5] = 0;
		JHsgny[6] = 0;
		JHsgny[7] = 0; 				

		JHsites[8][0] = 0;		// H-bonds,  horizontal above
		JHsites[8][1] = 2;
		JHsites[9][0] = 2;		// H-bonds,  right-aligned
		JHsites[9][1] = 4;
		JHsites[10][0] = 0;		// H-bonds,  left-aligned
		JHsites[10][1] = 6;
		JHsites[11][0] = 4;		// H-bonds,  horizontal below
		JHsites[11][1] = 6;						
		JHsgnx[8] = 0;
		JHsgnx[9] = 0;
		JHsgnx[10] = 0;
		JHsgnx[11] = 0; 
		JHsgny[8] = 1;
		JHsgny[9] = 1;
		JHsgny[10] = 1;
		JHsgny[11] = 1; 

		
		JHsites[12][0] = 1;		// H-bonds,  horizontal above
		JHsites[12][1] = 3;
		JHsites[13][0] = 1;		// H-bonds,  right-aligned
		JHsites[13][1] = 7;
		JHsites[14][0] = 3;		// H-bonds,  left-aligned
		JHsites[14][1] = 5;
		JHsites[15][0] = 5;		// H-bonds,  horizontal below
		JHsites[15][1] = 7;
		JHsgnx[12] = 0;
		JHsgnx[13] = 0;
		JHsgnx[14] = 0;
		JHsgnx[15] = 0; 
		JHsgny[12] = 1;
		JHsgny[13] = 1;
		JHsgny[14] = 1;
		JHsgny[15] = 1; 		

	
		JBsgnx = new int [Nb[2]];
		JBsgny = new int [Nb[2]];		
		JBsites = new int*[Nb[2]];
		for (int i = 0; i < Nb[2]; i++) JBsites[i] = new int[4];
		JBsites[0][0] = 0;
		JBsites[0][1] = 1;
		JBsites[0][2] = 4;
		JBsites[0][3] = 5;	
		JBsgnx[0] = 1;
		JBsgny[0] = 0;		
		//
		JBsites[1][0] = 0;
		JBsites[1][1] = 5;
		JBsites[1][2] = 1;
		JBsites[1][3] = 4;	
		JBsgnx[1] = 1;
		JBsgny[1] = 0;				
		
		JBsites[2][0] = 2;
		JBsites[2][1] = 3;
		JBsites[2][2] = 6;
		JBsites[2][3] = 7;
		JBsgnx[2] = 1;
		JBsgny[2] = 0;				
		//
		JBsites[3][0] = 2;
		JBsites[3][1] = 7;
		JBsites[3][2] = 3;
		JBsites[3][3] = 6;
		JBsgnx[3] = 1;
		JBsgny[3] = 0;	
		//***************
		//
		JBsites[4][0] = 0;
		JBsites[4][1] = 2;
		JBsites[4][2] = 4;
		JBsites[4][3] = 6;	
		JBsgnx[4] = 0;
		JBsgny[4] = 1;			
		
		//
		JBsites[5][0] = 0;
		JBsites[5][1] = 6;
		JBsites[5][2] = 2;
		JBsites[5][3] = 4;		// bb
		JBsgnx[5] = 0;
		JBsgny[5] = 1;		
		
		//
		JBsites[6][0] = 1;
		JBsites[6][1] = 3;
		JBsites[6][2] = 5;
		JBsites[6][3] = 7;	
		JBsgnx[6] = 0;
		JBsgny[6] = 1;				
		//
		JBsites[7][0] = 1;
		JBsites[7][1] = 7;
		JBsites[7][2] = 3;
		JBsites[7][3] = 5;
		JBsgnx[7] = 0;
		JBsgny[7] = 1;				
		
	
		
	}		
	else {
		// Square lattice. 
		// Heisenberg bonds.

		// Square lattice. 
		// Heisenberg bonds.
		JHsgnx = new int [Nb[0]];
		JHsgny = new int [Nb[0]];
		JHsites = new int *[Nb[0]];
		for (int i = 0; i < Nb[0]; ++i)JHsites[i] = new int[2];
		// QQ bonds.
		JQsgnx = new int [Nb[1]];
		JQsgny = new int [Nb[1]];
		JQsites = new int *[Nb[1]];
		for (int i = 0; i < Nb[1]; i++)JQsites[i] = new int[8];
		// Biquad bonds.
		JBsgnx = new int [Nb[2]];
		JBsgny = new int [Nb[2]];		
		JBsites = new int*[Nb[2]];
		for (int i = 0; i < Nb[2]; ++i)JBsites[i] = new int[4];
	
		//
		//int N1 = Nb[0]/Ns;
		//int N2 = Nb[1]/Ns;
		
		for (int y1 = 0; y1 < ly; ++y1)
		{
			for (int x1 = 0; x1 < lx; ++x1)
			{
				// ***************************************************************************************
				// 												Heisenberg sites
				// ***************************************************************************************		
				int s1 = x1 + y1 * lx;
				int x2 = (x1 + 1) % lx;
				int s2 = x2 + y1 * lx;
			
				// Heisenberg sites
				// Horizontal bonds
				JHsites[4*s1][0] = s1;//(s1 < s2) ? s1 : s2;		// H bonds,  horizontal above
				JHsites[4*s1][1] = s2;//(s1 < s2) ? s2 : s1;
			
				JHsites[4*s1+1][0] = s1 + Nm / 2;		// H-bonds,  slanted right
				JHsites[4*s1+1][1] = s2;

				JHsites[4*s1+2][0] = s1;		// H-bonds,  slanted left
				JHsites[4*s1+2][1] = s2 + Nm / 2;

				JHsites[4*s1+3][0] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 	// H-bonds,  horizontal below
				JHsites[4*s1+3][1] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 
				
				JHsgnx[4*s1]   = 1;
				JHsgnx[4*s1+1] = 1;
				JHsgnx[4*s1+2] = 1;
				JHsgnx[4*s1+3] = 1; 
				
				JHsgny[4*s1]   = 0;
				JHsgny[4*s1+1] = 0;
				JHsgny[4*s1+2] = 0;
				JHsgny[4*s1+3] = 0; 	

        
								
				// vertical bonds
				x2 = x1;
				int y2 = (y1+1) % ly;
				s2 = x2+y2*lx;
				
				
				JHsites[4*Ns+4*s1][0] = s1;//(s1 < s2) ? s1 : s2;		// H bonds,  horizontal above
				JHsites[4*Ns+4*s1][1] = s2;//(s1 < s2) ? s2 : s1;
				
				JHsites[4*Ns+4*s1+1][0] = s1 + Nm / 2;;		// H-bonds,  slanted right
				JHsites[4*Ns+4*s1+1][1] = s2;


				JHsites[4*Ns+4*s1+2][0] = s1;		// H-bonds,  slanted left
				JHsites[4*Ns+4*s1+2][1] = s2 + Nm / 2;
				
				JHsites[4*Ns+4*s1+3][0] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;	// H-bonds,  horizontal below
				JHsites[4*Ns+4*s1+3][1] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 
				JHsgnx[4*Ns+4*s1]   = 0;
				JHsgnx[4*Ns+4*s1+1] = 0;
				JHsgnx[4*Ns+4*s1+2] = 0;
				JHsgnx[4*Ns+4*s1+3] = 0; 
				
				JHsgny[4*Ns+4*s1]   = 1;
				JHsgny[4*Ns+4*s1+1] = 1;
				JHsgny[4*Ns+4*s1+2] = 1;
				JHsgny[4*Ns+4*s1+3] = 1; 


				// ***************************************************************************************
				//												Q of Q bonds
				// ***************************************************************************************	
				s1 = x1 + y1 * lx;
				x2 = (x1 + 1) % lx; 
				s2 = x2 + y1 * lx;
				
				int s3 = (s1 + lx)%Ns;
				int s4 = (s2 + lx)%Ns;
				
				// QoQ bonds			
				JQsites[4*s1][0] = s1;//(s1 < s2) ? s1 : s2;
				JQsites[4*s1][1] = s2;//(s1 < s2) ? s2 : s1;
				JQsites[4*s1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 
				JQsites[4*s1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[4*s1][4] = s3;//(s3 < s4) ? s3 : s4;
				JQsites[4*s1][5] = s4;//(s3 < s4) ? s4 : s3;
				JQsites[4*s1][6] = s3 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[4*s1][7] = s4 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 
				
				//
				JQsites[4*s1+1][0] = s1;//(s1 < s2) ? s1 : s2;
				JQsites[4*s1+1][1] = s2;//(s1 < s2) ? s2 : s1;
				JQsites[4*s1+1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 
				JQsites[4*s1+1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;

				JQsites[4*s1+1][4] = s3;
				JQsites[4*s1+1][5] = s4 + Nm/2;
				JQsites[4*s1+1][6] = s3 + Nm/2;
				JQsites[4*s1+1][7] = s4;
				
				//
				JQsites[4*s1+2][0] = s1;
				JQsites[4*s1+2][1] = s2 + Nm/2;
				JQsites[4*s1+2][2] = s1 + Nm/2;
				JQsites[4*s1+2][3] = s2;

				JQsites[4*s1+2][4] = s3;//(s3 < s4) ? s3 : s4;
				JQsites[4*s1+2][5] = s4;//(s3 < s4) ? s4 : s3;
				JQsites[4*s1+2][6] = s3 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2;
				JQsites[4*s1+2][7] = s4 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 	
				
				//
				JQsites[4*s1+3][0] = s1;
				JQsites[4*s1+3][1] = s2 + Nm/2;
				JQsites[4*s1+3][2] = s1 + Nm/2;
				JQsites[4*s1+3][3] = s2;

				JQsites[4*s1+3][4] = s3;
				JQsites[4*s1+3][5] = s4 + Nm/2;
				JQsites[4*s1+3][6] = s3 + Nm/2;
				JQsites[4*s1+3][7] = s4;				
				JQsgnx[4*s1]   = 1;
				JQsgnx[4*s1+1] = 1;
				JQsgnx[4*s1+2] = 1;
				JQsgnx[4*s1+3] = 1; 
				
				JQsgny[4*s1]   = 0;
				JQsgny[4*s1+1] = 0;
				JQsgny[4*s1+2] = 0;
				JQsgny[4*s1+3] = 0; 
				s1 = x1 + y1 * lx;
				s2 = (s1 + lx)%Ns;
				
				s3 = x2 + y1 * lx;
				s4 = (s3 + lx)%Ns;
				
				// std::cout <<  "   is  " << s1 << "    " << s2 << "    " << s3 << "    " << s4 << std::endl;
				JQsites[4*Ns+4*s1][0] = s1;//(s1 < s2) ? s1 : s2;
				JQsites[4*Ns+4*s1][1] = s2;//(s1 < s2) ? s2 : s1;
				JQsites[4*Ns+4*s1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JQsites[4*Ns+4*s1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[4*Ns+4*s1][4] = s3;//(s3 < s4) ? s3 : s4;
				JQsites[4*Ns+4*s1][5] = s4;//(s3 < s4) ? s4 : s3;
				JQsites[4*Ns+4*s1][6] = s3 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[4*Ns+4*s1][7] = s4 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 
				
				//
				JQsites[4*Ns+4*s1+1][0] = s1;//(s1 < s2) ? s1 : s2;
				JQsites[4*Ns+4*s1+1][1] = s2;//(s1 < s2) ? s2 : s1;
				JQsites[4*Ns+4*s1+1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JQsites[4*Ns+4*s1+1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[4*Ns+4*s1+1][4] = s3;
				JQsites[4*Ns+4*s1+1][5] = s4 + Nm/2;
				JQsites[4*Ns+4*s1+1][6] = s3 + Nm/2;
				JQsites[4*Ns+4*s1+1][7] = s4;
				
				//
				JQsites[4*Ns+4*s1+2][0] = s1;
				JQsites[4*Ns+4*s1+2][1] = s2 + Nm/2;
				JQsites[4*Ns+4*s1+2][2] = s1 + Nm/2;
				JQsites[4*Ns+4*s1+2][3] = s2;

				JQsites[4*Ns+4*s1+2][4] = s3;//(s3 < s4) ? s3 : s4;
				JQsites[4*Ns+4*s1+2][5] = s4;//(s3 < s4) ? s4 : s3;
				JQsites[4*Ns+4*s1+2][6] = s3 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[4*Ns+4*s1+2][7] = s4 + Nm/2;//(s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2;	
				
				//
				JQsites[4*Ns+4*s1+3][0] = s1;
				JQsites[4*Ns+4*s1+3][1] = s2 + Nm/2;
				JQsites[4*Ns+4*s1+3][2] = s1 + Nm/2;
				JQsites[4*Ns+4*s1+3][3] = s2;

				JQsites[4*Ns+4*s1+3][4] = s3;
				JQsites[4*Ns+4*s1+3][5] = s4 + Nm/2;
				JQsites[4*Ns+4*s1+3][6] = s3 + Nm/2;
				JQsites[4*Ns+4*s1+3][7] = s4;
				
				JQsgnx[4*Ns+4*s1]   = 0;
				JQsgnx[4*Ns+4*s1+1] = 0;
				JQsgnx[4*Ns+4*s1+2] = 0;
				JQsgnx[4*Ns+4*s1+3] = 0; 
				
				JQsgny[4*Ns+4*s1]   = 1;
				JQsgny[4*Ns+4*s1+1] = 1;
				JQsgny[4*Ns+4*s1+2] = 1;
				JQsgny[4*Ns+4*s1+3] = 1; 

				// ***************************************************************************************
				//			Biqaduratic interaction bonds 
				// ***************************************************************************************
				// Biquad sites
				// Horizontal bonds
				s1 = x1 + y1 * lx;
				x2 = (x1 + 1) % lx;
				s2 = x2 + y1 * lx;				

				JBsites[2*s1][0] = s1;//
				JBsites[2*s1][1] = s2;//
				JBsites[2*s1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;//
				JBsites[2*s1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;//
				JBsites[2*s1+1][0] = s1;
				JBsites[2*s1+1][1] = s2 + Nm/2;
				JBsites[2*s1+1][2] = s1 + Nm/2;
				JBsites[2*s1+1][3] = s2;
				JBsgnx[2*s1]   = 1;
				JBsgnx[2*s1+1] = 1;
				JBsgny[2*s1]   = 0;
				JBsgny[2*s1+1] = 0;

				// Vertical bonds				
				x2 = x1;
				y2 = (y1+1) % ly;
				s2 = x2+y2*lx;

				//std::cout << JBsites[4*s1][0] << "   " << JBsites[4*s1][1] << "   " << JBsites[4*s1][2] << "   " << JBsites[4*s1][3] << std::endl;
				//std::cout << 2*s1 << "  " << 2*s1+1 << "  " << 2*Ns+2*s1 << "  " << 2*Ns+2*s1+1 << std::endl;
				
								
				JBsites[2*Ns+2*s1][0] = s1;//(s1 < s2) ? s1 : s2;//
				JBsites[2*Ns+2*s1][1] = s2;//(s1 < s2) ? s2 : s1;//
				JBsites[2*Ns+2*s1][2] = s1 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;//s
				JBsites[2*Ns+2*s1][3] = s2 + Nm/2;//(s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;//
				JBsites[2*Ns+2*s1+1][0] = s1;
				JBsites[2*Ns+2*s1+1][1] = s2 + Nm/2;
				JBsites[2*Ns+2*s1+1][2] = s1 + Nm/2;
				JBsites[2*Ns+2*s1+1][3] = s2;
				JBsgnx[2*Ns+2*s1]   = 0;
				JBsgnx[2*Ns+2*s1+1] = 0;
				JBsgny[2*Ns+2*s1]   = 1;
				JBsgny[2*Ns+2*s1+1] = 1;						
					
			}
		}
	}	
}





/*

				// Heisenberg sites
				// Horizontal bonds
				JHsites[4*s1][0] = (s1 < s2) ? s1 : s2;		// H bonds,  horizontal above
				JHsites[4*s1][1] = (s1 < s2) ? s2 : s1;
			
				JHsites[4*s1+1][0] = s2;		// H-bonds,  slanted right
				JHsites[4*s1+1][1] = s1 + Nm / 2;

				JHsites[4*s1+2][0] = s1;		// H-bonds,  slanted left
				JHsites[4*s1+2][1] = s2 + Nm / 2;

				JHsites[4*s1+3][0] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 	// H-bonds,  horizontal below
				JHsites[4*s1+3][1] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 
				
				// vertical bonds
				x2 = x1;
				int y2 = (y1+1) % ly;
				s2 = x2+y2*lx;
				
				
				JHsites[4*Ns+4*s1][0] = (s1 < s2) ? s1 : s2;		// H bonds,  horizontal above
				JHsites[4*Ns+4*s1][1] = (s1 < s2) ? s2 : s1;
				
				JHsites[4*Ns+4*s1+1][0] = s2;		// H-bonds,  slanted right
				JHsites[4*Ns+4*s1+1][1] = s1 + Nm / 2;


				JHsites[4*Ns+4*s1+2][0] = s1;		// H-bonds,  slanted left
				JHsites[4*Ns+4*s1+2][1] = s2 + Nm / 2;
				
				JHsites[4*Ns+4*s1+3][0] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;	// H-bonds,  horizontal below
				JHsites[4*Ns+4*s1+3][1] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 


				// ***************************************************************************************
				//												Q of Q bonds
				// ***************************************************************************************	
				s1 = x1 + y1 * lx;
				x2 = (x1 + 1) % lx; 
				s2 = x2 + y1 * lx;
				
				int s3 = (s1 + lx)%Ns;
				int s4 = (s2 + lx)%Ns;
				
				// QoQ bonds			
				JQsites[N2*s1][0] = (s1 < s2) ? s1 : s2;
				JQsites[N2*s1][1] = (s1 < s2) ? s2 : s1;
				JQsites[N2*s1][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 
				JQsites[N2*s1][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[N2*s1][4] = (s3 < s4) ? s3 : s4;
				JQsites[N2*s1][5] = (s3 < s4) ? s4 : s3;
				JQsites[N2*s1][6] = (s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[N2*s1][7] = (s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 
				
				//
				JQsites[N2*s1+1][0] = (s1 < s2) ? s1 : s2;
				JQsites[N2*s1+1][1] = (s1 < s2) ? s2 : s1;
				JQsites[N2*s1+1][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2; 
				JQsites[N2*s1+1][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;

				JQsites[N2*s1+1][4] = s3;
				JQsites[N2*s1+1][5] = s4 + Nm/2;
				JQsites[N2*s1+1][6] = s4;
				JQsites[N2*s1+1][7] = s3 + Nm/2;
				
				//
				JQsites[N2*s1+2][0] = s1;
				JQsites[N2*s1+2][1] = s2 + Nm/2;
				JQsites[N2*s1+2][2] = s2;
				JQsites[N2*s1+2][3] = s1 + Nm/2;

				JQsites[N2*s1+2][4] = (s3 < s4) ? s3 : s4;
				JQsites[N2*s1+2][5] = (s3 < s4) ? s4 : s3;
				JQsites[N2*s1+2][6] = (s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2;
				JQsites[N2*s1+2][7] = (s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 	
				
				//
				JQsites[N2*s1+3][0] = s1;
				JQsites[N2*s1+3][1] = s2 + Nm/2;
				JQsites[N2*s1+3][2] = s2;
				JQsites[N2*s1+3][3] = s1 + Nm/2;

				JQsites[N2*s1+3][4] = s3;
				JQsites[N2*s1+3][5] = s4 + Nm/2;
				JQsites[N2*s1+3][6] = s4;
				JQsites[N2*s1+3][7] = s3 + Nm/2;				

				s1 = x1 + y1 * lx;
				s2 = (s1 + lx)%Ns;
				
				s3 = x2 + y1 * lx;
				s4 = (s3 + lx)%Ns;
				
				// std::cout <<  "   is  " << s1 << "    " << s2 << "    " << s3 << "    " << s4 << std::endl;
				JQsites[N2*s1+4][0] = (s1 < s2) ? s1 : s2;
				JQsites[N2*s1+4][1] = (s1 < s2) ? s2 : s1;
				JQsites[N2*s1+4][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JQsites[N2*s1+4][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[N2*s1+4][4] = (s3 < s4) ? s3 : s4;
				JQsites[N2*s1+4][5] = (s3 < s4) ? s4 : s3;
				JQsites[N2*s1+4][6] = (s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[N2*s1+4][7] = (s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2; 
				
				//
				JQsites[N2*s1+5][0] = (s1 < s2) ? s1 : s2;
				JQsites[N2*s1+5][1] = (s1 < s2) ? s2 : s1;
				JQsites[N2*s1+5][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JQsites[N2*s1+5][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2; 

				JQsites[N2*s1+5][4] = s3;
				JQsites[N2*s1+5][5] = s4 + Nm/2;
				JQsites[N2*s1+5][6] = s4;
				JQsites[N2*s1+5][7] = s3 + Nm/2;
				
				//
				JQsites[N2*s1+6][0] = s1;
				JQsites[N2*s1+6][1] = s2 + Nm/2;
				JQsites[N2*s1+6][2] = s2;
				JQsites[N2*s1+6][3] = s1 + Nm/2;

				JQsites[N2*s1+6][4] = (s3 < s4) ? s3 : s4;
				JQsites[N2*s1+6][5] = (s3 < s4) ? s4 : s3;
				JQsites[N2*s1+6][6] = (s3 + Nm/2 < s4 + Nm/2) ? s3 + Nm/2 : s4 + Nm/2; 
				JQsites[N2*s1+6][7] = (s3 + Nm/2 < s4 + Nm/2) ? s4 + Nm/2 : s3 + Nm/2;	
				
				//
				JQsites[N2*s1+7][0] = s1;
				JQsites[N2*s1+7][1] = s2 + Nm/2;
				JQsites[N2*s1+7][2] = s2;
				JQsites[N2*s1+7][3] = s1 + Nm/2;

				JQsites[N2*s1+7][4] = s3;
				JQsites[N2*s1+7][5] = s4 + Nm/2;
				JQsites[N2*s1+7][6] = s4;
				JQsites[N2*s1+7][7] = s3 + Nm/2;


				// ***************************************************************************************
				//			Biqaduratic interaction bonds 
				// ***************************************************************************************
				// Biquad sites
				// Horizontal bonds
				s1 = x1 + y1 * lx;
				x2 = (x1 + 1) % lx;
				s2 = x2 + y1 * lx;				

				JBsites[2*s1][0] = (s1 < s2) ? s1 : s2;
				JBsites[2*s1][1] = (s1 < s2) ? s2 : s1;
				JBsites[2*s1][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JBsites[2*s1][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;
				JBsites[2*s1+1][0] = s1;
				JBsites[2*s1+1][1] = s2 + Nm/2;
				JBsites[2*s1+1][2] = s2;
				JBsites[2*s1+1][3] = s1 + Nm/2;				


				// Vertical bonds				
				x2 = x1;
				y2 = (y1+1) % ly;
				s2 = x2+y2*lx;

				//std::cout << JBsites[4*s1][0] << "   " << JBsites[4*s1][1] << "   " << JBsites[4*s1][2] << "   " << JBsites[4*s1][3] << std::endl;
				//std::cout << 2*Ns+2*s1+1 << std::endl;
				
								
				JBsites[2*Ns+2*s1][0] = (s1 < s2) ? s1 : s2;
				JBsites[2*Ns+2*s1][1] = (s1 < s2) ? s2 : s1;
				JBsites[2*Ns+2*s1][2] = (s1 + Nm/2 < s2 + Nm/2) ? s1 + Nm/2 : s2 + Nm/2;
				JBsites[2*Ns+2*s1][3] = (s1 + Nm/2 < s2 + Nm/2) ? s2 + Nm/2 : s1 + Nm/2;
				JBsites[2*Ns+2*s1+1][0] = s1;
				JBsites[2*Ns+2*s1+1][1] = s2 + Nm/2;
				JBsites[2*Ns+2*s1+1][2] = s2;
				JBsites[2*Ns+2*s1+1][3] = s1 + Nm/2;	
				
				//std::cout << s1 << "    " << JHsites[4*s1][0] << "  " << JHsites[4*s1][1] << "    " << JHsites[4*s1+1][0] << "  " << JHsites[4*s1+1][1] << "    " <<JHsites[4*s1+2][0] << "  " << JHsites[4*s1+2][1] << "    " << JHsites[4*s1+3][0] << "  " << JHsites[4*s1+3][1] << std::endl;
				//std::cout << s1 << "    " << JHsites[4*Ns+4*s1][0] << "  " << JHsites[4*Ns+4*s1][1] << "    " << JHsites[4*Ns+4*s1+1][0] << "  " << JHsites[4*Ns+4*s1+1][1] << "    " <<JHsites[4*Ns+4*s1+2][0] << "  " << JHsites[4*Ns+4*s1+2][1] << "    " << JHsites[4*Ns+4*s1+3][0] << "  " << JHsites[4*Ns+4*s1+3][1] << std::endl;				
				//std::cout << 	JQsites[N2*s1+4][0]	<< "  " << JQsites[N2*s1+4][1] << "  " << JQsites[N2*s1+4][2] << "  " << JQsites[N2*s1+4][3] << "      " << JQsites[N2*s1+4][4] << "  " << JQsites[N2*s1+4][5] << "  " << JQsites[N2*s1+4][6] << "  " << JQsites[N2*s1+4][7] << "  " << std::endl;  	
				//std::cout << JBsites[2*s1][0] << "   " << JBsites[2*s1][1] << "   " << JBsites[2*s1][2] << "   " << JBsites[2*s1][3] << std::endl;

*/
