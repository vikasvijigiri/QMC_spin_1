#ifndef _SSEUPDATES_HPP_DEFINED_
#define _SSEUPDATES_HPP_DEFINED_

class SSEupdates {
      	public: 
      	SSEupdates(){}

  // ********************************************************************
  // 		Update related structures 
  // ********************************************************************

	long **str, *pstr;
  int *frst, *last; 
 	static const int nvx=6;
 	const double pi = 2. * acos(0.);
	double prob_in, prob_rm, cum_prob[3]; 	
	
	int Nd[3], No[3], type[3], flip[5][2];
	
 	double amax;
 	double wgt[2][2];
 	double awgt[2][2];
 	double dwgt[2][2];
	
 	int vxoper[nvx+1];
 	int legvx[2][2][2][2];
 	int vxleg[4][nvx+1];
 	int vxnew[4][4][nvx+1];
 	double vxprb[4][4][nvx+1];
 	double vxp[4][4][nvx+1];

	
 	int nl;
	
	void weights();
	void initvrtx_dirloop();

	void mcstep(int, SSElattice* );
	void mcstep_measurement(SSElattice* );
	void checkl();
	
	
	void diag_update(int, SSElattice * );				// Method 1
	//void diag_update_IJ(int, SSElattice * ); // Method 2
	//void diag_update_IQ(int, SSElattice * ); // Method 2
	void looper(SSElattice* );
	void diag_update_measurement(SSElattice * ); // Method 1
	void looper_measurement(SSElattice *  );
  void initialize(SSElattice * );
  
  
  // Uncomment this if you want to debug something.	
	/*
	void print_state_string(SSElattice * );	
	void print_state(SSElattice * );	
	void print_op_string(int );	
	void print_proj_op_string(int );
	bool check_states(SSElattice *, int * );
	int* copy_state(SSElattice * );
  */
  
  
  
  // ********************************************************************
  // 		Measurements related data 
  // ********************************************************************
  
  
  // Observables
  // local Variables
  double enrg; //*enrg_bin, 
  double O_N, O_N1, O_N2, O_N3, R_N;
  double O_V, O_V1, O_V2, O_V3, O_V4, R_V;
  double O_B, O_B1, O_B2, O_B3, O_B4, R_B;
  double *B, **Cx, **Cy, **Dx, **Dy; // Off-diagonal correlations
  //double **Cxx, **Cyy, **Dxx, **Dyy; // Off-diagonal correlations  

  double mag, mag_square, mag_four, stiff;
  
  // Lattice position vectors  	
  int **psite;
      	
  // Functions
  void observables_processing(int , int );//, writeresults& );
  std::pair<double, double>  Jackniffe(double *, double , int );
  void Initiate_observables();
  void neel_measure(SSElattice* );

};
#endif

