#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"
//#include "../include/writeresults.hpp"
#include "../include/SSEupdates.hpp"
#include <iostream>
#include <iomanip>
//#include <utility>
//#include <cmath>

using namespace std;


void
SSEupdates::observables_processing(int bins, int iters)//, writeresults& write)
{
  double V, x, y;
  // VBS order params, no.1
  //***********************
  double V_realx1 = 0.;
  double V_realx2 = 0.;
  double V_realy1 = 0.;
  double V_realy2 = 0.;

  // VBS order params, no.2
  //***********************
  double B_realx1 = 0.;
  double B_realx2 = 0.;
  double B_realy1 = 0.;
  double B_realy2 = 0.;
  
  // Neel order params.
  //***********************  
  double N_real3 = 0.;
  double N_real2 = 0.;
  double N_real1 = 0.;  
  	
  
  // std::cout << "vikas " << std::endl;
  // Neel and VBS Measurements  
  if (lx==2){
 		for (int j = 0; j < Ns; ++j) { 
		  V = B[j]/iters;
		  x = psite[j][0];
		  y = psite[j][1];
		  //std::cout << B[j] << std::endl;
		  
		  N_real1 += cos(x * (pi + 2. * pi / lx) + y * pi) * V;
		  N_real2 += cos(x * pi + y * (pi + 2. * pi / lx)) * V;
		  N_real3 += cos(x * pi + y * pi) * V;    
    }  
    //*****************
		int i, j;
	 	i = 0;
	 	j = 0; 
	 	x = 0;
	 	y = 0; 
		V = ( Cx[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		V_realx2 += cos(x * pi) * V; 
		std::cout  << "(0,0)  " << V/(Ns*Ns) << "   " << Cx[i][j]  << "  " << n1 << std::endl;		
		V = ( Dx[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		B_realx2 += cos(x * pi) * V; 
	 	//*************
	 	i = 1;
	 	j = 0; 
	 	x = 0;
	 	y = 1; 
		V = ( Cx[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		V_realx2 += cos(x * pi) * V;
		std::cout  << "(1,0)  " << V/(Ns*Ns) << "   " << Cx[i][j] << "  " << n1 << std::endl;			
		V = ( Dx[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		B_realx2 += cos(x * pi) * V;  
	 	//*************
	 	i = 1;
	 	j = 1; 
	 	x = 0;
	 	y = 0; 
		V = ( Cx[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		V_realx2 += cos(x * pi) * V;
		std::cout  << "(1,1)  " << V/(Ns*Ns) << "   " << Cx[i][j] << "  " << n1 << std::endl;			
		V = ( Dx[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		B_realx2 += cos(x * pi) * V;  
		//###########################################################3
		i = 0;
	 	j = 0; 
	 	x = 0;
	 	y = 0; 
		V = ( Cy[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		V_realy2 += cos(y * pi) * V;
		V = ( Dy[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		B_realy2 += cos(y * pi) * V;  
	 	//*************
	 	i = 1;
	 	j = 0; 
	 	x = 1;
	 	y = 0; 
		V = ( Cy[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		V_realy2 += cos(y * pi) * V;
		V = ( Dy[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		B_realy2 += cos(y * pi) * V;   
	 	//*************
	 	i = 1;
	 	j = 1; 
	 	x = 0;
	 	y = 0; 
		V = ( Cy[i][j]) / (Beta * Beta * J[0] * J[0])/iters;
		V_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		V_realy2 += cos(y * pi) * V;  
		V = ( Dy[i][j]) / (Beta * Beta * J[2] * J[2])/iters;
		B_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		B_realy2 += cos(y * pi) * V;  
  }
  else{
		for (int j = 0; j < Ns; ++j) {
		  //V = B[j]/Lc/iters;
		  V = B[j]/iters;
		  x = psite[j][0];
		  y = psite[j][1];
		  //std::cout << B[j] << std::endl;
		  
		  N_real1 += cos(x * (pi + 2. * pi / lx) + y * pi) * V;
		  N_real2 += cos(x * pi + y * (pi + 2. * pi / lx)) * V;
		  N_real3 += cos(x * pi + y * pi) * V;    
	 	
		  for (int i = j; i < Ns; ++i) {
		    x = float(psite[j][0] - psite[i][0]);
		    y = float(psite[j][1] - psite[i][1]);
		         
		    // VBS no.1
		    V = Cx[i][j] / (Beta * Beta * J[0] * J[0])/iters; Cx[i][j]=0.;
		    V_realx1 += cos(x * pi + y * (2. * pi / lx)) * V; 
		    V_realx2 += cos(x * pi) * V;
		//std::cout << i << "  " << j  << "   " << Cx[i][j]  << std::endl;					    
		    V = Cy[i][j] / (Beta * Beta * J[0] * J[0])/iters; Cy[i][j]=0.;
		    V_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		    V_realy2 += cos(y * pi) * V;
		//std::cout << i << "  " << j  << "   " << Cy[i][j]  << std::endl;					    		    

		    // VBS no.2
		    V = Dx[i][j] / (Beta * Beta * J[2] * J[2])/iters; Dx[i][j]=0.;
		    B_realx1 += cos(x * pi + y * (2. * pi / lx)) * V;
		    B_realx2 += cos(x * pi) * V;
		//std::cout << i << "  " << j  << "   " << Dx[i][j]  << std::endl;					    		    
		    V = Dy[i][j] / (Beta * Beta * J[2] * J[2])/iters; Dy[i][j]=0.;
		    B_realy1 += cos(x * (2. * pi / lx) + y * pi) * V;
		    B_realy2 += cos(y * pi) * V;
		//std::cout << i << "  " << j  << "   " << Dy[i][j]  << std::endl;					    		    		    
		  }
		}
  }
  //std::cout << "Bin : " << bins << " val  : " << N_real3 << std::endl << std::endl;
  //std::cout << std::endl;

  // Neel  order parameter.
  O_N1 = (N_real1);
  O_N2 = (N_real2);
  O_N3 = (N_real3);
  O_N = O_N3/float(Ns);
	R_N = 1. - (O_N1+O_N2)/(2.*O_N3);
	//std:: cout << O_N1 << "   " << O_N2 << "  " << O_N3 << std::endl;		


  
  // VBS  order parameter. 1st param
  O_V1 = (V_realy1);
  O_V2 = (V_realx1);
  O_V3 = (V_realy2);
  O_V4 = (V_realx2);
	O_V = (O_V3 + O_V4)/(2.*float(Ns*Ns));	  
	R_V = 1. - (O_V1/O_V3+ O_V2/O_V4)/2.;
	//std:: cout << O_V1 << "   " << O_V2 << "  " << O_V3 << "   " << O_V4  << std::endl;	

  // VBS  order parameter. 2nd param
  O_B1 = (B_realy1);
  O_B2 = (B_realx1);
  O_B3 = (B_realy2);
  O_B4 = (B_realx2);
	O_B = (O_B3 + O_B4)/(2.*float(Ns*Ns));  
	R_B = 1. - (O_B1/O_B3 + O_B2/O_B4)/2.;
	//std:: cout << O_B1 << "   " << O_B2 << "  " << O_B3 << "   " << O_B4  << std::endl;	


	// Normalizing
	enrg  /= iters;
	mag   /= iters;
	mag_square   /= iters;
	mag_four   /= iters;
	stiff /= iters;
	
	enrg   = -enrg / (Beta * Ns);
	mag    = mag;// / float(Ns);
	mag_square = mag_square;// / float(Ns * Ns); 
	mag_four = mag_four;// / float(Ns * Ns * Ns * Ns); 
	stiff  = stiff / (Beta * Ns);
}

/*
std::pair<double, double>  SSEupdates::Jackniffe(double *x, double x_tot, int nbins){
    // Do the jackknife estimates
    double *x_jack = new double [nbins];
    for (int i=0; i<nbins; i++)
        x_jack[i] = (x_tot - x[i]) / float(nbins - 1);
    //double x_av = x_tot / nbins; // do the overall averages
    //double g_av  = x_av;
    double g_jack_av = 0.; double g_jack_err = 0.0;
    double dg;
    for (int i=0; i<nbins; i++){ // do the final jackknife averages
        dg = x_jack[i];
        g_jack_av += dg;
        g_jack_err += dg*dg;
    }    
    g_jack_av /= nbins;
    g_jack_err /= nbins;
    g_jack_err = sqrt((nbins - 1) * std::abs(g_jack_err - g_jack_av*g_jack_av));
    //print(" Overall average is %8.4f"% g_av)
    //print(" Jackknife average is %8.4f +/- %6.4f"% (g_jack_av, g_jack_err))
    return std::make_pair(g_jack_av, g_jack_err);
}    
*/

void
SSEupdates::Initiate_observables()
{

  enrg  = 0.;
  stiff = 0.;
  mag   = 0.;
  mag_square = 0.;
  mag_four = 0.;
	//***********
	Cx = new double *[Ns];
	Cy = new double *[Ns];
	/*
	Cxx = new double *[Ns];
	Cyy = new double *[Ns];
	Dxx = new double *[Ns];
	Dyy = new double *[Ns];		
	*/		
	B = new double  [Ns];
	Dx = new double *[Ns];
	Dy = new double *[Ns];	

	//***********
	for (int k=0; k<Ns; k++){
		Cx[k] = new double [Ns];
		Cy[k] = new double [Ns];
		Dx[k] = new double [Ns];
		Dy[k] = new double [Ns];
		/*
		Cxx[k] = new double [Ns];
		Cyy[k] = new double [Ns];
		Dxx[k] = new double [Ns];
		Dyy[k] = new double [Ns];		
		*/
	}	
	
	//***********
	for (int i=0; i<Ns; i++){
	 	B[i] = 0.;
	} 	
			

	
	for (int i=0; i<Ns; i++){
		for (int j=0; j<Ns; j++){
			Cx[i][j] = 0.; Cy[i][j] = 0.;
			Dx[i][j] = 0.; Dy[i][j] = 0.;
			//Cxx[i][j] = 0.; Cyy[i][j] = 0.;
			//Dxx[i][j] = 0.; Dyy[i][j] = 0.;		
		} 
	}		
	
	
	// Position vectors  
	int z1=0, y1=0;//, g1=0, g2=0;

	//NL = 2 * Ns;
	psite = new int *[Ns];
	
	for (int k=0; k<Ns; k++){
		psite[k] = new int [2];
	}

	for (int k = 1; k <= Ns; k++)
	{
		//Block 1(Bottom left)
		if (k % lx <= int((lx) / 2) && k <= int((Ns / 2) + lx))
		{
			if (k % lx != 0)
			{
				z1 = k % lx - 1;
				y1 = int((k - 1) / lx);
				//g1 = z1 + y1;
				//g2 = y1 - z1;
			}
			else
			{
				z1 = -k % lx - 1;
				y1 = int((k - 1) / lx);
				//g1 = z1 + y1;
				//g2 = y1 - z1;
			}
		}
		//Block 2(Bottom right)
		else if (k % lx > int((lx) / 2) && k <= int((Ns / 2) + lx))
		{
			z1 = k % lx - 1 - lx;
			y1 = int((k - 1) / lx);
			//g1 = z1 + y1;
			//g2 = y1 - z1;
		}
		//Block 3(top left)
		else if (k % lx <= int((lx) / 2) && k >= int((Ns / 2) + lx))
		{
			z1 = k % lx - 1;
			y1 = int((k - 1) / lx) - lx;
			//g1 = z1 + y1;
			//g2 = y1 - z1;
		}
		//Block 4(top right)
		else if (k % lx > int((lx) / 2) && k >= int((Ns / 2) + lx))
		{
			z1 = k%lx - 1 - lx;
			y1 = int((k - 1) / lx) - lx;
			//g1 = z1 + y1;
			//g2 = y1 - z1;
		}
		psite[k - 1][0] = z1;//g1 - 1;
		psite[k - 1][1] = y1;//g2 - 1;
	}
	/*
	//std::cout 	<< " Initialized " << std::endl;
	for (int k=0; k<Ns; k++){
		std::cout << "psite is  "<< k << "       " << psite[k][0] << "  " << psite[k][1] << std::endl;
	}
	*/	

}
