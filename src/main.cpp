#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"
#include "../include/SSEupdates.hpp"			//Headers
//#include "SSEobservables.hpp"
#include "../include/writeresults.hpp"
//#include <mpi.h>
#include <chrono>
#include <iostream>

static ran rann;

int main() {
  SSEvariables vars; // Creating all necessary instances
  SSEupdates update;
  //SSEobservables obs;
  writeresults write;

  //MPI_Init(NULL, NULL);
  //int myrank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  //int numprocs;
  //MPI_Comm_size(MPI_COMM_WORLD, &numprocs);	
	//std::cout << rann() << std::endl;
  vars.declare_variables();
  //vars.set_temperatures();

  SSElattice * lattice; //Declare spin array of SSE class.

  lattice = new SSElattice[Nm];

  for (int xc = 0; xc < Nm; ++xc) {
    lattice[xc].set_S((rann() <= 0.5) ? -1 : 1); // Start with random configuration
  }

  std::cout << "J_Heisenberg: " << J[0] << std::endl;
  std::cout << "J_QoQ: " << J[1] << std::endl;
  std::cout << "J_Biquad: " << J[2] << std::endl;
  std::cout << "Beta: " << Beta << std::endl;
  std::cout << "L: " << lx << std::endl;
  std::cout << "Ms iters: " << iter << std::endl;
  std::cout << "Eq iters: " << isteps << std::endl;
  //std::cout << "nbins: " << nbins << std::endl;

  vars.lattice_sites();
  update.initvrtx_dirloop();
  update.initialize(lattice);
  update.weights();

  /* 
  //Print the initial random configuration
  for (int xc = 0; xc < Nm; ++xc) {
   std::cout << "xc : " <<  lattice[xc].S() << '\n'; // Start with random configuration
  }
  */
  
  // Start measuring time
  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < isteps; ++i) {
    // Equilibration  to check how much time does it take to complete.  
    update.mcstep(i, lattice);
    update.checkl();
  }

  // Start measuring time
  auto t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast < std::chrono::nanoseconds > (t2 - t1);
  double time = (elapsed.count() * 1e-9) / 60.;
  std::cout << "Total Equilibration time ~ " << time << ", minutes" << std::endl;

  update.Initiate_observables();

  //for (int j = 0; j < nbins; ++j) {
    // Measurements
    for (int i = 0; i < iter; ++i) {
      //update.mcstep(lattice);
      update.mcstep_measurement(lattice);
      //update.neel_measure(lattice);
    }
    update.observables_processing(nbins, iter); //, write);
  	write.output(update.enrg, update.mag, update.mag_square, update.mag_four, update.stiff, update.O_N, update.R_N,
		  					 update.O_V, update.R_V, update.O_B, update.R_B);        
  //}
  t1 = std::chrono::high_resolution_clock::now();
  elapsed = std::chrono::duration_cast < std::chrono::nanoseconds > (t1 - t2);
  time = (elapsed.count() * 1e-9) / 60.;
  std::cout << "Total Measurements time: " << time <<  " minutes." << std::endl;
  //if (myrank == 0)
  //{
  //std::cout << "Total Actual Execution time = " << timeeq + timems << ", minutes" << std::endl;
  //}

  delete[] lattice;
  //MPI_Finalize();
  return 0;
}
