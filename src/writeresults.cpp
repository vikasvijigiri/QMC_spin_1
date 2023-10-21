#include "../include/SSEvariables.hpp"
#include "../include/writeresults.hpp"
#include <sstream>
//#include <string>
#include <fstream>
#include <iostream>
//#include <mpi.h>


using namespace std;

//***************************************************************************************
void writeresults::output(double enrg_local, double mag_local, double mag_square_local, double mag_four_local,  double stiff_local, double O_N_local, double R_N_local,  
double O_V_local,  double R_V_local,  double O_B_local,  double R_B_local)
{
	/*
	double *temp_global  = NULL;
	double *enrg_global  = NULL;
	double *mag_global   = NULL;
	double *mag_square_global   = NULL;
	double *stiff_global = NULL;
	double *O_N_global   = NULL;
	double *O_N1_global  = NULL;
	double *O_N2_global  = NULL;
	double *R_N_global   = NULL;
			
	double *O_V_global   = NULL;
	double *O_V1_global  = NULL;
	double *O_V2_global  = NULL;
	double *O_V3_global  = NULL;
	double *R_V_global   = NULL;
	
	double *O_B_global   = NULL;
	double *O_B1_global  = NULL;
	double *O_B2_global  = NULL;
	double *O_B3_global  = NULL;
	double *R_B_global   = NULL;
	
	temp_global  = (double*) malloc(numprocs* sizeof(double));
	enrg_global  = (double*) malloc(numprocs* sizeof(double));
	mag_global   = (double*) malloc(numprocs* sizeof(double));
	mag_square_global   = (double*) malloc(numprocs* sizeof(double));
	stiff_global = (double*) malloc(numprocs* sizeof(double));
	O_N_global   = (double*) malloc(numprocs* sizeof(double));
	O_N1_global  = (double*) malloc(numprocs* sizeof(double));
	O_N2_global  = (double*) malloc(numprocs* sizeof(double));
	R_N_global   = (double*) malloc(numprocs* sizeof(double));	
			
	O_V_global   = (double*) malloc(numprocs* sizeof(double));
	O_V1_global  = (double*) malloc(numprocs* sizeof(double));
	O_V2_global  = (double*) malloc(numprocs* sizeof(double));
	O_V3_global  = (double*) malloc(numprocs* sizeof(double));
	R_V_global   = (double*) malloc(numprocs* sizeof(double));	
				
	O_B_global   = (double*) malloc(numprocs* sizeof(double));
	O_B1_global  = (double*) malloc(numprocs* sizeof(double));
	O_B2_global  = (double*) malloc(numprocs* sizeof(double));
	O_B3_global  = (double*) malloc(numprocs* sizeof(double));
	R_B_global   = (double*) malloc(numprocs* sizeof(double));			


	MPI_Gather(&temp_local,  1, MPI_DOUBLE, temp_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	MPI_Gather(&enrg_local,  1, MPI_DOUBLE, enrg_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&mag_local,   1, MPI_DOUBLE, mag_global,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&mag_square_local,   1, MPI_DOUBLE, mag_square_global,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&stiff_local, 1, MPI_DOUBLE, stiff_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&O_N_local,   1, MPI_DOUBLE, O_N_global,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&O_N1_local,  1, MPI_DOUBLE, O_N1_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&O_N2_local,  1, MPI_DOUBLE, O_N2_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Gather(&O_V_local,   1, MPI_DOUBLE, O_V_global,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	MPI_Gather(&O_V1_local,  1, MPI_DOUBLE, O_V1_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	MPI_Gather(&O_V2_local,  1, MPI_DOUBLE, O_V2_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);	
	MPI_Gather(&O_V3_local,  1, MPI_DOUBLE, O_V3_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);				
	
	MPI_Gather(&O_B_local,   1, MPI_DOUBLE, O_B_global,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);		
	MPI_Gather(&O_B1_local,  1, MPI_DOUBLE, O_B1_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);		
	MPI_Gather(&O_B2_local,  1, MPI_DOUBLE, O_B2_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);		
	MPI_Gather(&O_B3_local,  1, MPI_DOUBLE, O_B3_global,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);					
	*/
		
	//if (myrank == 0)
	//{
		//std::cout << "  " << std::endl;
		/*
		std::cout << "O(N) is     ------------->   " << 1- O_N_global[0]/(J[0]) << std::endl;
		std::cout << "O(V) is     ------------->   " << 1- O_V_global[0]/(J[0]) << std::endl;
		std::cout << "(1 - 0.5*(S0.S1 + S2.S3)) is " << C << std::endl;
		std::cout << "(S0.s1)(S2.S3) is ------->   " << O_B_global[0] + C - 1/(J[0]) << std::endl;
		*/
		/*	
		double C1 = 1. - O_N_global[0]/(J[0]);
		double C2 = 1. - O_V_global[0]/(J[0]);

		std::cout << "<(S0.S1)>     ------------->   " << C2 << std::endl;
		std::cout << "<(S2.S3)>     ------------->   " << C1 << std::endl;
		std::cout << "<(S0.S1)(S2.S3)> ---------->   " << O_B_global[0]/(J[0]*J[0]) << std::endl;	
		std::cout << "<(Sz0*Sz1)(Sz2*Sz3)> ------>   " << O_B1_global[0]  << std::endl;	
		*/	
		
		// Build the filename for the given process

	  ofstream file;
    file.open("output_seeded_data_input_details.txt");
			file << "J: " << J[0] << '\n';
			file << "QoQ: " << J[1] << '\n';
			file << "Biquad: " << J[2] << '\n';
			file << "T: " << 1/Beta << '\n';
			file << "L: " << lx << '\n';
			file << "Eq steps: " << isteps << '\n';
			file << "No. of bins: " << nbins << '\n';
			file << "Iter^n inside bins: " << iter << '\n';
			//file << "Energy: " << enrg_global[0]*Ns << '\n';
			file << '\n';	
			//file << "      -----------------------------------------------------------------------------------------------------------------------------------------------------------------";
			file << '\n';
			file << "         J_QQ     ";				
			file << "       Enrg       ";
			file << "     mag       ";
			file << "    mag^2       ";
			file << "    stif       ";
			file << "    O_N        ";
			file << "    R_N      ";
			file << "     O_V      ";
			file << "      R_V      ";
			file << "      O_B      ";
			file << "      R_B      ";
			file << '\n';	
			//file << "      -----------------------------------------------------------------------------------------------------------------------------------------------------------------;
			file << '\n';	
			file << '\n';		
    file.close();
    
		// Append to the same file 
	  ofstream files;
    files.open("output_seeded_data_avg.txt", ios::app);
    
    if (!files){
				std::cout << "No file found. Error in file creation! " << std::endl;
		}			
    else {
    			// Else append to the file
						files << setw(15) << std::scientific << enrg_local/* *Ns */;
						//files << setw(15) << std::scientific << mag_local/* *Ns */;
						files << setw(15) << std::scientific << mag_square_local/* *Ns */;
						files << setw(15) << std::scientific << mag_four_local/* *Ns */;						
						files << setw(15) << std::scientific << stiff_local/* *Ns */;
						files << setw(15) << std::scientific << O_N_local; //O_N_local/4.;
						files << setw(15) << std::scientific << R_N_local;
						files << setw(15) << std::scientific << O_V_local;
						files << setw(15) << std::scientific << R_V_local;
						files << setw(15) << std::scientific << O_B_local;//O_B_local;
						files << setw(15) << std::scientific << R_B_local;//R_B_local;						
						files << '\n';
        files.close();
    }
}


//***************************************************************************************
template < typename A>
	std::string writeresults::to_string_with_precision(A a_value, int n)
	{
		std::ostringstream out;
		out << std::setprecision(n) << a_value;
		return out.str();
	}
//***************************************************************************************
