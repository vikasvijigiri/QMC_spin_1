#include "../include/SSEvariables.hpp"
#include "../include/SSElattice.hpp"
#include "../include/SSEupdates.hpp"
#include <iostream>

static ran rann;
using namespace std;

//********************************************************************************************


void
SSEupdates::checkl()
{
	long Lc_new = long(n1 + n1 / 3);

	if (Lc_new > Lc)
	{
		int ns = 4;
		long **tmp;
		tmp = new long *[ns];
		for (long i = 0; i < ns; ++i)
			tmp[i] = new long[Lc_new];
		
		// Store in a temporary array..
		for (int j = 0; j < ns; ++j)
		{
			for (long i = 0; i < Lc_new; ++i)
			{
				if (i < Lc)
				{
					tmp[j][i] = str[j][i];
				}
				else
				{
					tmp[j][i] = 0;
				}
			}
		}
		
		// Delete old string.
		for (int i = 0; i < ns; ++i)
		{
			delete[] str[i];
		}
		delete[] str;
		
		// Allot values to newly created arrays.
		Lc = Lc_new;
		str = new long *[ns];
		for (int i = 0; i < ns; ++i) str[i] = new long[Lc];

		for (int j = 0; j < ns; ++j)
		{
			for (long i = 0; i < Lc; ++i)
			{
				str[j][i] = tmp[j][i];
			}
		}

		for (int i = 0; i < ns; ++i)
		{
			delete[] tmp[i];
		}
		delete[] tmp;
    //long l = Np * Lc;
	 	//cout << "Cutoff length is   " << Np << "  " << Lc << "   " << l << endl;
	}
}


//********************************************************************************************

void
SSEupdates::initvrtx_dirloop()
{

	int v0, v1, v2, v3, v4, v5;

	int vx, ic, oc, vxn;
	int st[4];

	// Define legvx, vxoper, vxwgt, and vxleg
	v0 = 1;
	v1 = v0 + 1;
	v2 = v0 + 2;
	v3 = v0 + 3;
	v4 = v0 + 4;
	v5 = v0 + 5;

	for (int t2 = 0; t2 < 2; ++t2)
	{
		for (int t1 = 0; t1 < 2; ++t1)
		{
			for (int s2 = 0; s2 < 2; ++s2)
			{
				for (int s1 = 0; s1 < 2; ++s1)
				{
					legvx[s1][s2][t1][t2] = 0;
				}
			}
		}
	}
	// ss1 ss2 ---> tt1 tt2
	legvx[1][0][0][1] = v0;	//offdiagonal
	legvx[0][1][1][0] = v1;	//offdiagonal
	legvx[1][0][1][0] = v2;	//diagonal
	legvx[0][1][0][1] = v3;	//diagonal
	legvx[0][0][0][0] = v4;	//diagonal
	legvx[1][1][1][1] = v5;	//diagonal

	vxoper[v0] = 2;	//offdiagonal
	vxoper[v1] = 2;	//offdiagonal
	vxoper[v2] = 1;	//diagonal
	vxoper[v3] = 1;	//diagonal
	vxoper[v4] = 1;	//diagonal
	vxoper[v5] = 1;	//diagonal      

	for (int t2 = 0; t2 < 2; ++t2)
	{
		for (int t1 = 0; t1 < 2; ++t1)
		{
			for (int s2 = 0; s2 < 2; ++s2)
			{
				for (int s1 = 0; s1 < 2; ++s1)
				{
					vx = legvx[s1][s2][t1][t2];
					if (vx != 0)
					{
						vxleg[0][vx] = s1;
						vxleg[1][vx] = s2;
						vxleg[2][vx] = t1;
						vxleg[3][vx] = t2;
					}
				}
			}
		}
	}

	// Create vxnew links from input channel (ic) ---> output channel (oc).
	for (vx = 1; vx <= nvx; ++vx)
	{
		for (ic = 0; ic < 4; ++ic)
		{
			for (oc = 0; oc < 4; ++oc)
			{
				st[0] = vxleg[0][vx];
				st[1] = vxleg[1][vx];
				st[2] = vxleg[2][vx];
				st[3] = vxleg[3][vx];
				st[ic] = 1 - st[ic];
				st[oc] = 1 - st[oc];
				vxn = legvx[st[0]][st[1]][st[2]][st[3]];
				if (vxn != 0)
				{
					vxnew[oc][ic][vx] = vxn;
				}
			}
		}
	}

	// Initialize cumulative probabilities saving in vxprob.
	for (vx = 1; vx <= nvx; ++vx)
	{
		for (ic = 0; ic < 4; ++ic)
		{
			for (oc = 0; oc < 4; ++oc)
			{
				vxprb[oc][ic][vx] = 0;
			}
		}
	}

	//*************************** Projection operator updates *******************************
	int s1, s2, t1, t2;

	s1 = 1;
	s2 = 1;
	t1 = 1;
	t2 = 1;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 2;
	vxprb[oc][ic][vx] = 0.5;
	oc = 3;
	vxprb[oc][ic][vx] = 0.5;

	ic = 1;
	oc = 2;
	vxprb[oc][ic][vx] = 0.5;
	oc = 3;
	vxprb[oc][ic][vx] = 0.5;

	ic = 2;
	oc = 1;
	vxprb[oc][ic][vx] = 0.5;
	oc = 0;
	vxprb[oc][ic][vx] = 0.5;

	ic = 3;
	oc = 0;
	vxprb[oc][ic][vx] = 0.5;
	oc = 1;
	vxprb[oc][ic][vx] = 0.5;

	//*************************            
	s1 = 0;
	s2 = 0;
	t1 = 0;
	t2 = 0;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 2;
	vxprb[oc][ic][vx] = 0.5;
	oc = 3;
	vxprb[oc][ic][vx] = 0.5;

	ic = 1;
	oc = 2;
	vxprb[oc][ic][vx] = 0.5;
	oc = 3;
	vxprb[oc][ic][vx] = 0.5;

	ic = 2;
	oc = 1;
	vxprb[oc][ic][vx] = 0.5;
	oc = 0;
	vxprb[oc][ic][vx] = 0.5;

	ic = 3;
	oc = 0;
	vxprb[oc][ic][vx] = 0.5;
	oc = 1;
	vxprb[oc][ic][vx] = 0.5;

	//####################################################

	//*************************            
	s1 = 0;
	s2 = 1;
	t1 = 1;
	t2 = 0;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 3;
	vxprb[oc][ic][vx] = 1.0;

	ic = 1;
	oc = 2;
	vxprb[oc][ic][vx] = 1.0;

	ic = 2;
	oc = 1;
	vxprb[oc][ic][vx] = 1.0;

	ic = 3;
	oc = 0;
	vxprb[oc][ic][vx] = 1.0;

	//*************************            
	s1 = 0;
	s2 = 1;
	t1 = 0;
	t2 = 1;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 2;
	vxprb[oc][ic][vx] = 1.0;

	ic = 1;
	oc = 3;
	vxprb[oc][ic][vx] = 1.0;

	ic = 2;
	oc = 0;
	vxprb[oc][ic][vx] = 1.0;

	ic = 3;
	oc = 1;
	vxprb[oc][ic][vx] = 1.0;

	//####################################################

	//*************************            
	s1 = 1;
	s2 = 0;
	t1 = 0;
	t2 = 1;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 3;
	vxprb[oc][ic][vx] = 1.0;

	ic = 1;
	oc = 2;
	vxprb[oc][ic][vx] = 1.0;

	ic = 2;
	oc = 1;
	vxprb[oc][ic][vx] = 1.0;

	ic = 3;
	oc = 0;
	vxprb[oc][ic][vx] = 1.0;

	//*************************            
	s1 = 1;
	s2 = 0;
	t1 = 1;
	t2 = 0;
	vx = legvx[s1][s2][t1][t2];

	ic = 0;
	oc = 2;
	vxprb[oc][ic][vx] = 1.0;

	ic = 1;
	oc = 3;
	vxprb[oc][ic][vx] = 1.0;

	ic = 2;
	oc = 0;
	vxprb[oc][ic][vx] = 1.0;

	ic = 3;
	oc = 1;
	vxprb[oc][ic][vx] = 1.0;

}

//********************************************************************************************
void
SSEupdates::weights()
{
	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			wgt[s1][s2] = 2.0 * 0.25 * pow(-1, 1 + s1) *pow(-1, 1 + s2);
			if (wgt[s1][s2] > amax)
				amax = wgt[s1][s2];
		}
	}

	// amax[j]+=1.0f;
	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			wgt[s1][s2] = amax - wgt[s1][s2];
		}
	}

	for (int s2 = 0; s2 < 2; ++s2)
	{
		for (int s1 = 0; s1 < 2; ++s1)
		{
			awgt[s1][s2] = wgt[s1][s2];
			if (awgt[s1][s2] > 1e-6)
			{
				dwgt[s1][s2] = 1.0 / awgt[s1][s2];
			}
			else
			{
				dwgt[s1][s2] = 1.e6;
			}
		}
	}
	
	//std::cout << awgt[0][1] << "  " << awgt[1][0] << std::endl; 

 	// Which operator to insert,  probabilities...
 	prob_in = 0.0;
 	prob_rm = 0.0;
 	cum_prob[0] = 0.5 * J[0] * Nb[0] *Beta;
 	cum_prob[1] = 0.5 * 0.5 * J[1] * Nb[1] *Beta;
 	cum_prob[2] = 0.5 * J[2] * Nb[2] *Beta;
 	for (int i=0; i<3; i++){
 		//cum_prob[i] = pow(2.0, (i == 1) ? 2 : 0)* pow(0.5,type[i]) * J[i] * Nb[i] *Beta;
 		//cout << "i is " << i  << "  "<< cum_prob[i] << endl;
 		prob_in += cum_prob[i];
 	}	
 	for (int i=0; i<3; i++){
 		
 		cum_prob[i] /= prob_in;
 		//cout << "Normalized cumulative prob "<< cum_prob[i] << endl;
 	}
 	prob_rm = 1.0/prob_in;	

	// vertex type after flipping
	//  1 <---> 2
	flip[1][0] = 2;
	flip[1][0] = 2;	
	
	flip[2][0] = 1;	
	flip[2][0] = 1;
	
	//  1 <---> 3
	flip[1][1] = 3;	
	flip[1][1] = 3;	

	flip[3][1] = 1;	
	flip[3][1] = 1;	

	//  2 <---> 4
	flip[2][1] = 4;	
	flip[2][1] = 4;	
	
	flip[4][1] = 2;	
	flip[4][1] = 2;			

	//  3 <---> 4
	flip[4][0] = 3;	
	flip[4][0] = 3;	
	
	flip[3][0] = 4;	
	flip[3][0] = 4;		
}


void
SSEupdates::initialize(SSElattice *lattice)
{
	Lc = 30; // Initial cutoff length
	nl = 3;  // Initial no. of loops
	n1 = 0;  // Initial diagonal operators
	int ss1, ss2, ns = 4;



	// H-bond part
	str = new long *[ns];
	pstr = new long [Np];
	for (int i = 0; i < ns; i++)
	{
		str[i] = new long[Lc];
	}
				
	for (int i = 0; i < ns; i++)
	{
		for (long j = 0; j < Lc; j++)
		{
			str[i][j] = 0;
		}
	}
	/*
	// Projector part
	for (int j = 0; j < Np; j++)
	{
		pstr[j] = 0;
		//lattice[j].set_S(pow(-1,j));
		lattice[j].set_S((rann() < 0.5) ? 1 : -1);
		lattice[j + Np].set_S(-lattice[j].S());
	}
	*/
	np1 = 0;
	for (int j = 0; j < Np; j++)
	{
		np1 += 1;
		ss1 = (1 + lattice[j].S()) / 2;
		ss2 = (1 + lattice[j + Np].S()) / 2;

		// Template projector bonds at every bond and at each operator slice placed after Lc.
		pstr[j] = legvx[ss1][ss2][ss1][ss2];
	}
	// Template S-bonds
	frst = new int[Nm];
	last = new int[Nm];
	
	for (int i = 0; i < Nm; i++)
	{
		frst[i] = -1;
		last[i] = -1;
	}
				
     	
	type[0] = 1;
	type[1] = 4;
	type[2] = 2;
	
	Nd[0] = 0;
	Nd[1] = 0;
	Nd[2] = 0;
	
}
