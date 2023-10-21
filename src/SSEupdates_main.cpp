#include "../include/SSEvariables.hpp"
#include "../include/SSElattice.hpp"
#include "../include/SSEupdates.hpp"
#include <iostream>

static ran rann;
using namespace std;

void
SSEupdates::mcstep(int w, SSElattice * lattice) {
  diag_update(w, lattice);
  looper(lattice);
}

void
SSEupdates::mcstep_measurement(SSElattice * lattice) {
  diag_update_measurement(lattice);
  looper_measurement(lattice);
}

void
SSEupdates::diag_update(int w, SSElattice * lattice) {
  int b, k = 0, o, op, ss1, ss2;
  double p, r, cp;
  /*
  int dum_latt[128];

  //int bn[128];
	for (int i = 0; i < Nm; ++i) {
	    //fn = 0;
			dum_latt[i] = lattice[i].S();
			//std::cout << " up site: " << i << "   " << lattice[i].S() <<"   " << dum_latt[i] << "          " << n1 << std::endl;			
	}	
	*/
	
	/*
	int wc=1000000;
	int* lattice_prev;
	// Seperate loop to insert the diagonal QQ.
	if(w<wc){	
		//cout << "\niter: " << w << ", initial : ";
		//print_state(lattice);
		lattice_prev=copy_state(lattice); 
	}
	*/
  //No = 0;
  for (long i = 0; i < Lc; ++i) {
    o = str[0][i] * str[1][i];
    if (o == 0) {
      r = rann();
      cp = 0.0;
      for (k = 0; k < 3; k++) {
        cp += cum_prob[k];
        if (cp > r and cum_prob[k] > 1e-6) { // J or Q or Biquad?
          b = int(rann() * Nb[k]);
          p = 1.0;
          if (k == 0) {
            ss1 = (1 + lattice[JHsites[b][0]].S()) / 2;
            ss2 = (1 + lattice[JHsites[b][1]].S()) / 2;
            p = p * awgt[ss1][ss2];
            //p = p * Beta* Nb[0] * J[0] * 0.5  / float(Lc-n1);
          } else if (k == 1) {
            for (int j = 0; j < 4; ++j) {
              ss1 = (1 + lattice[JQsites[b][2 * j]].S()) / 2;
              ss2 = (1 + lattice[JQsites[b][2 * j + 1]].S()) / 2;
              p = p * awgt[ss1][ss2];
            }
          } else if (k == 2) {
            for (int j = 0; j < 2; ++j) {
              ss1 = (1 + lattice[JBsites[b][2 * j]].S()) / 2;
              ss2 = (1 + lattice[JBsites[b][2 * j + 1]].S()) / 2;
              p = p * awgt[ss1][ss2];
            }
          }
          //p = p * Beta *Nb[1] * J[1] * 0.5 * 0.5 / float(Lc-n1);
          p = p * prob_in / float(Lc - n1);
          if (rann() <= p) // Bond
          {
            str[0][i] = 1; // Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to QoQ and Bi and H);  1 for diag 2 for off-diag
            str[1][i] = 1; // Off-digonal (or On-diag) operator within the 2nd biquad bond (relevant to QoQ));  1 for diag 2 for off-diag
            str[2][i] = k; // Which type of operator, 0 for H, 1 for  QQ, 2 for Biquad;
            str[3][i] = b + 1; // Bond-index;
            Nd[k] += 1;
            n1 += 1;
          }
          break;
        }
      }
    } else if (o == 1) {
      k = str[2][i];
      p = prob_rm * float(Lc - n1 + 1);
      if (rann() <= p) {
        Nd[k] -= 1;
        str[0][i] = 0;
        str[1][i] = 0;
        str[2][i] = 0;
        str[3][i] = 0;
        n1 -= 1;
      }
    } else {
      b = str[3][i] - 1;
      k = str[2][i];
      if (k == 0) {
        lattice[JHsites[b][0]].flip();
        lattice[JHsites[b][1]].flip();
      } else if (k == 1) {
        for (int j = 0; j < 2; j++) {
          op = str[j][i];
          if (op == 2) {
            lattice[JQsites[b][0 + 4 * j]].flip();
            lattice[JQsites[b][1 + 4 * j]].flip();
          } else if (op == 3) {
            lattice[JQsites[b][2 + 4 * j]].flip();
            lattice[JQsites[b][3 + 4 * j]].flip();
          } else if (op == 4) {
            lattice[JQsites[b][0 + 4 * j]].flip();
            lattice[JQsites[b][1 + 4 * j]].flip();
            lattice[JQsites[b][2 + 4 * j]].flip();
            lattice[JQsites[b][3 + 4 * j]].flip();
          }
        }
      } else if (k == 2) {
        op = str[0][i];
        if (op == 2) {
          lattice[JBsites[b][0]].flip();
          lattice[JBsites[b][1]].flip();
        } else if (op == 3) {
          lattice[JBsites[b][2]].flip();
          lattice[JBsites[b][3]].flip();
        } else if (op == 4) {
          lattice[JBsites[b][0]].flip();
          lattice[JBsites[b][1]].flip();
          lattice[JBsites[b][2]].flip();
          lattice[JBsites[b][3]].flip();
        }
      }
    }
  //std::cout << " time: " << i << "  op str : " <<  str[0][i] << " " << str[1][i] << "  " <<str[2][i] << " bond : " << str[3][i] - 1 << " site: " <<  JHsites[b][0] << "," << JHsites[b][1]<< endl;    
  }
  /*
	for (int i = 0; i < Nm; ++i) {
			//std::cout << " dn site: " << i << "   " << lattice[i].S() <<"   " << dum_latt[i] << "          " << No << std::endl;				
			if (dum_latt[i] != lattice[i].S()){
          printf("wrong vertex ");
          exit(0);				
			}
	}
	*/	
	/*
  for (int i = 0; i < Np; ++i)
	{
		op = vxoper[pstr[i]];	
		if (op == 2)
		{
			//b = pstr[1][i] - 1;
			//s1 = i; //Jsites[b][0];
			//s2 = i+Np; //Jsites[b][1];
			lattice[i].flip();
			lattice[i+Np].flip();
		}

		if(w<wc){
				print_proj_op_string(i);	
				print_state(lattice);
		}	

	}
		if(w<wc){
		bool check=check_states(lattice,lattice_prev);
		cout << "\niter: " << w << ", final   : ";
		//print_state(lattice);
		cout << " " << std::boolalpha << check;
		cout<<endl;
		}
	*/			
}

void
SSEupdates::looper(SSElattice * lattice) {
  long o, v, vi, v0, v1, v2;
  int vx, ic, oc;
  int s1, s2;
  int b, k;
  double r, p;

  //int a = ((Nd[0]+No[0]) + 4 * (Nd[1]+No[1]) + 2 * (Nd[2]+No[2]) + np1);
  int a = ((Nd[0]) + 4 * (Nd[1]) + 2 * (Nd[2]) + np1);
  //std::cout << 4*a << std::endl;
  //int a = (4*Lc + np1);

  long * X = new long[4 * a];
  long * visited_vertex = new long[Lc + Np];
  long * lpos = new  long[a];

  //std::fill_n(frst, Nm, -1);
  //std::fill_n(last, Nm, -1);
  std::fill_n(visited_vertex, Lc, 0);
  
  // Linked list construction
  v0 = 0;
  for (long i = 0; i < Lc; ++i) {
    o = str[0][i];
    if (o != 0) {
         
      b = str[3][i] - 1;
      v = str[2][i];
  //if (o == 1)std::cout << " dn vikas time: " << i << "  op str : " <<  str[0][i] << " " << str[1][i] << "  " <<str[2][i] << " bond : " << b << "  " << " site " << JHsites[b][0]<<"," << JHsites[b][1] << endl; 
      for (int j = 0; j < type[v]; j++) {
        if (v == 0) {
          s1 = JHsites[b][0];
          s2 = JHsites[b][1];
        } else if (v == 1) {
          s1 = JQsites[b][2 * j];
          s2 = JQsites[b][2 * j + 1];
        } else if (v == 2) {
          s1 = JBsites[b][2 * j];
          s2 = JBsites[b][2 * j + 1];
        } else {
          printf("wrong vertex ");
          exit(0);
        }
        lpos[long(v0 / 4)] = Np * i + j;
        v1 = last[s1];
        v2 = last[s2];

        if (v1 != -1) {
          X[v1] = v0;
          X[v0] = v1;
        //std::cout << "sankalp bhaiyya " << i << "   " << v0 << ","<< X[v0] << std::endl;              
        } else {
          frst[s1] = v0;
        }
        if (v2 != -1) {
          X[v2] = v0 + 1;
          X[v0 + 1] = v2;
        } else {
          frst[s2] = v0 + 1;
        }
        last[s1] = v0 + 2;
        last[s2] = v0 + 3;
        v0 += 4;         
      }
    }
  }
  
  
  // Projector part
  for (int i = 0; i < Np; ++i) {
    o = pstr[i];
    if (o != 0) {
      s1 = i;
      s2 = i + Np;
      lpos[long(v0 / 4)] = Np * (i + Lc);
      v1 = last[s1];
      v2 = last[s2];

      if (v1 != -1) {
        X[v1] = v0;
        X[v0] = v1;
      } else {
        frst[s1] = v0;
      }
      if (v2 != -1) {
        X[v2] = v0 + 1;
        X[v0 + 1] = v2;
      } else {
        frst[s2] = v0 + 1;
      }
      last[s1] = v0 + 2;
      last[s2] = v0 + 3;
      v0 += 4;
    }
  }

  // PBC loops 
  for (int k = 0; k < Nm; ++k) {
    v1 = frst[k];
    if (v1 != -1) {
      v2 = last[k];
      X[v2] = v1;
      X[v1] = v2;
    }
  }
	//std::cout << "nl: " << nl << " v0: " << v0 << " a: " << a  << "  " << No[0] << "  " << No[1]  << "  " << No[2] << "  " << Np << "  " << Lc << "  " << Np*(Np+Lc) << std::endl;
  /*
  for (int i=0; i<v0; i++){
    std::cout<< "love " << i << "   " << X[i] << " --> " << std::endl;
  }
	*/
  

	nl = 5*lx;
  // Main loop 
  long i;
  long l;
  for (int j = 0; j < nl; j++) {
    vi = long(rann() * v0);
    v1 = vi;
    l = lpos[long(v1 / 4)];
    //std::cout << l << std::endl;
    i = long(l / Np);
    while (1) {
      ic = v1 % 4;
      if (i < Lc) {
        oc = ic ^ 1;
        k = str[2][i];
        //std::cout <<"Naba " << i << "  " << k << "  " << ic << "  " << oc  << std::endl;        
        if (k == 0) {
          str[0][i] = (1 + 2) - str[0][i];
        } else if (k == 1) {
          v = long((l % Np) / 2); // Index of which bond out of those 4
					//std::cout << v << "   " << i << "    " << (l % Np) % 2 << "    " << l << std::endl;
          str[v][i] = flip[str[v][i]][(l % Np) % 2];
        } else if (k == 2) {
          v = (l % Np) % 2;
          str[0][i] = flip[str[0][i]][v];
        }
      } else {
        r = rann();
        vx = pstr[i - Lc];
        p = 0;
        for (oc = 0; oc < 4; ++oc) {
          p += vxprb[oc][ic][vx];
          if (r < p) {
            pstr[i - Lc] = vxnew[oc][ic][vx];
            break;
          }
        }
      }
      v2 = 4 * (long(v1 / 4)) + oc;
      visited_vertex[i] += 1;
      /*
      if (v1 == v2 && v2 == vi)
      {
      	break;
      }
      */
      v = X[v2];
      if (v == vi) {
        break;
      } else {
        v1 = v;
        l = lpos[long(v1 / 4)];
        i = long(l / Np);
        //std::cout << "updated i: " << i << std::endl;
      }
    }
    //std::cout << std::endl;
  }
	/*
  int number_visited = 0;
  for (int j = 0; j < Lc; j++) {
    //i = int(lpos[j] / Np);
    k = str[2][j];
    number_visited += visited_vertex[j] / type[k];
    //visited_vertex[i] = 0;
  }
	
  if (number_visited  < n1) {
    nl += int(nl / 3);
  }
  */
  // cout << "nl is " << nl << endl; 
  // Update spin configuration here
  for (int j = 0; j < Nm; ++j) {
    // Valid only where there is projector operators.
    if (last[j] != -1) {
      long i = long(lpos[long(last[j] / 4)] / Np) - Lc;
      o = pstr[i];
      b = last[j] % 4;
      lattice[j].set_S(2 * vxleg[b][o] - 1);
    } else {
      if (rann() < 0.5) lattice[j].flip();
    }
    last[j] = -1;
    frst[j] = -1;
  }
  delete[] lpos;
  delete[] visited_vertex;
  delete[] X;
}




void
SSEupdates::diag_update_measurement(SSElattice * lattice) {
  int b, k, o, ss1, ss2;//, spn;
  double p, r, cp;
  //double b1[4], a1[4];
  // Neel order params
  //double N_real1 = 0.;
  //double N_real2 = 0.;
  //double N_real3 = 0.;

  //double N_imag1 = 0.;
  //double N_imag2 = 0.;
  //double N_imag3 = 0.;
  double mg = 0.;
  //double mag_abs = 0.;
  //double mag_sq = 0.;
  //double mag_fo = 0.;
  //double mag_fo = 0.;
  //double O_N_abs;
  //double O_N1_abs;
  //double O_N2_abs;
  long jjx = 0, jjy = 0;  


  
  
  mg = 0.;
  // Measurement of some diagonal quantities.
  for(int x = 0; x < Ns; x++){
			B[x] += (lattice[0].S() + lattice[Ns].S())/2.*(lattice[x].S() + lattice[Ns+x].S())/2.;       
			//B[x] += Bx[x];
			mg += 0.5*(lattice[x].S() + lattice[x+Ns].S()) * pow(-1, ((x % Ns) % lx) + (x % Ns) / lx);	
	}
	mg = mg/float(Ns);
	
	
	     
  for (long i = 0; i < Lc; ++i) {
    o = str[0][i] * str[1][i];
    if (o == 0) {
      r = rann();
      cp = 0.0;
      for (k = 0; k < 3; k++) {
        cp += cum_prob[k];
        if (cp > r and cum_prob[k] > 1e-7) { // J or Q?
          b = int(rann() * Nb[k]);
          p = 1.0;
          if (k == 0) {
            ss1 = (1 + lattice[JHsites[b][0]].S()) / 2;
            ss2 = (1 + lattice[JHsites[b][1]].S()) / 2;
            p = p * awgt[ss1][ss2];
          } else if (k == 1) {
            for (int j = 0; j < 4; ++j) {
              ss1 = (1 + lattice[JQsites[b][2 * j]].S()) / 2;
              ss2 = (1 + lattice[JQsites[b][2 * j + 1]].S()) / 2;
              p = p * awgt[ss1][ss2];
            }
          } else if (k == 2) {
            for (int j = 0; j < 2; ++j) {
              ss1 = (1 + lattice[JBsites[b][2 * j]].S()) / 2;
              ss2 = (1 + lattice[JBsites[b][2 * j + 1]].S()) / 2;
              p = p * awgt[ss1][ss2];
            }
          }
          p = p * prob_in / float(Lc - n1);
          if (rann() <= p) // Bond
          {
            str[0][i] = 1; // Off-digonal (or On-diag) operator within the 1st biquad bond (relevant to QoQ and Bi and H); 
            str[1][i] = 1; // Off-digonal (or On-diag) operator within the 2nd biquad bond (relevant to QoQ)); 
            str[2][i] = k; // Which type of operator, I, H, QQ, Biquad;
            str[3][i] = b + 1; // Bond-index;
            Nd[k] += 1;
            n1 += 1;
          }
          break;
        }
      }
    } else if (o == 1) {
      k = str[2][i];
      p = prob_rm * float(Lc - n1 + 1);
      if (rann() <= p) {
        Nd[k] -= 1;
        str[0][i] = 0;
        str[1][i] = 0;
        str[2][i] = 0;
        str[3][i] = 0;
        n1 -= 1;
      }
    } else {
      b = str[3][i] - 1;
      k = str[2][i];
      if (k == 0) { // Encountered off-digonal Heisenberg bond
					ss1 = JHsites[b][0];
					ss2 = JHsites[b][1];		
					//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();																		
					lattice[ss1].flip();
					lattice[ss2].flip();
					//if (spn != 0) {
						jjx+=JHsgnx[b]*lattice[ss2].S();
						jjy+=JHsgny[b]*lattice[ss2].S();
				  //}	
					//mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);			          				
      } else if (k == 1) { // Encountered off-digonal QoQ bond
		      for (int j = 0; j < 2; j++) {
		        o = str[j][i];
		        if (o == 2) {
		        	ss1 = JQsites[b][0 + 4 * j];
		        	ss2 = JQsites[b][1 + 4 * j];
							//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();		        	
		          lattice[ss1].flip();
		          lattice[ss2].flip();
		          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
							//if (spn != 0) {
								jjx+=JQsgnx[b]*lattice[ss2].S();
								jjy+=JQsgny[b]*lattice[ss2].S();
							//}			          
		        } else if (o == 3) {
		        	ss1 = JQsites[b][2 + 4 * j];
		        	ss2 = JQsites[b][3 + 4 * j];	
							//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();		        		        
		          lattice[ss1].flip();
		          lattice[ss2].flip();
		          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
							//if (spn != 0) {
								jjx+=JQsgnx[b]*lattice[ss2].S();
								jjy+=JQsgny[b]*lattice[ss2].S();
							//}					          
		        } else if (o == 4) {
		        	ss1 = JQsites[b][0 + 4 * j];
		        	ss2 = JQsites[b][1 + 4 * j];
							//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();		        				
		          lattice[ss1].flip();
		          lattice[ss2].flip();
		          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
							//if (spn != 0) {
								jjx+=JQsgnx[b]*lattice[ss2].S();
								jjy+=JQsgny[b]*lattice[ss2].S();
							//}			
		        	ss1 = JQsites[b][2 + 4 * j];
		        	ss2 = JQsites[b][3 + 4 * j];		          
							//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();
		          lattice[ss1].flip();
		          lattice[ss2].flip();
		          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
							//if (spn != 0) {
								jjx+=JQsgnx[b]*lattice[ss2].S();
								jjy+=JQsgny[b]*lattice[ss2].S();
							//}						          
		        }
		      }
      } else if (k == 2) { // Encountered off-digonal Bi-Bi bond
        o = str[0][i];
        if (o == 2) {
        	ss1 = JBsites[b][0];
        	ss2 = JBsites[b][1];
					//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();
          lattice[ss1].flip();
          lattice[ss2].flip();
          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
					//if (spn != 0) {
						jjx+=JBsgnx[b]*lattice[ss2].S();
						jjy+=JBsgny[b]*lattice[ss2].S();
					//}				          
        } else if (o == 3) {
        	ss1 = JBsites[b][2];
        	ss2 = JBsites[b][3];
					//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();
          lattice[ss1].flip();
          lattice[ss2].flip();
          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
					//if (spn != 0) {
						jjx+=JBsgnx[b]*lattice[ss2].S();
						jjy+=JBsgny[b]*lattice[ss2].S();
					//}		          
        } else if (o == 4) {
        	ss1 = JBsites[b][0];
        	ss2 = JBsites[b][1];
					//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();
          lattice[ss1].flip();
          lattice[ss2].flip();
          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
					//if (spn != 0) {
						jjx+=JBsgnx[b]*lattice[ss2].S();
						jjy+=JBsgny[b]*lattice[ss2].S();
					//}				
        	ss1 = JBsites[b][2];
        	ss2 = JBsites[b][3];  
					//spn = lattice[ss2].S()+lattice[(ss2 < Ns) ? ss2+Ns : ss2%Ns].S();        	        
          lattice[ss1].flip();
          lattice[ss2].flip();
          //mg += 2. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
					//if (spn != 0) {
						jjx+=JBsgnx[b]*lattice[ss2].S();
						jjy+=JBsgny[b]*lattice[ss2].S();
					//}		          
        }
      }
  	}
  //mag_abs += mg; 
  //mag_sq  += pow(mg, 2);  // Uncomment when you are averaging over imaginary time slice.  
  //mag_fo  += pow(mg, 4);   
  }  	
  //mag_abs /= float(Lc);
  //mag_sq  /= float(Lc);
  //mag_fo  /= float(Lc);
  //mag_sq  = pow(mg, 2);
  //mag_abs = abs(mg);  
  
  /*
  mag  += mag_abs;
  mag_square += mag_sq;
  mag_four += mag_fo; 
  */
  
  mag  += mg;
  mag_square += pow(mg, 2);
  mag_four += pow(mg, 4);    

  // Avg Stiffness along x and y;
  stiff += 0.5 * float(pow(jjx, 2) + pow(jjy, 2));
  // Avg Energy 
  enrg += float(n1);
}

// ****************************************************************************
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//  				Second loop starts here 

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
// ****************************************************************************	

void
SSEupdates::looper_measurement(SSElattice * lattice) {
  int o, v, vi, v0, v1, v2, vx, ic, oc;
  int s1, s2;
  int b, k;
  double r, p;

  //int a = ((Nd[0]+No[0]) + 4 * (Nd[1]+No[1]) + 2 * (Nd[2]+No[2]) + np1);
  int a = ((Nd[0]) + 4 * (Nd[1]) + 2 * (Nd[2]) + np1);
  int * X = new int[4 * a];
  long * lpos = new long[a];
  long * flag = new long[Lc];

  //std::fill_n(frst, Nm, -1);
  //std::fill_n(last, Nm, -1);

  // Linked list construction
  v0 = 0;
  for (long i = 0; i < Lc; ++i) {
    o = str[0][i];
    flag[i] = 0;
    if (o != 0) {
      b = str[3][i] - 1;
      v = str[2][i];
      //std::cout << i  << "   " << v << "   " << b << std::endl;              
      for (int j = 0; j < type[v]; j++) {
        if (v == 0) {
          s1 = JHsites[b][0];
          s2 = JHsites[b][1];
        } else if (v == 1) {
          s1 = JQsites[b][2 * j];
          s2 = JQsites[b][2 * j + 1];
        } else if (v == 2) {
          s1 = JBsites[b][2 * j];
          s2 = JBsites[b][2 * j + 1];
        } else {
          printf("wrong vertex ");
          exit(0);
        }
        lpos[int(v0 / 4)] = Np * i + j;
        v1 = last[s1];
        v2 = last[s2];

        if (v1 != -1) {
          X[v1] = v0;
          X[v0] = v1;
        } else {
          frst[s1] = v0;
        }
        if (v2 != -1) {
          X[v2] = v0 + 1;
          X[v0 + 1] = v2;
        } else {
          frst[s2] = v0 + 1;
        }
        last[s1] = v0 + 2;
        last[s2] = v0 + 3;
        v0 += 4;
      }
    }
  }
  //std::cout << std::endl;
  int  j, vb1, vb2, b2, b3, b4, b5 ;
  long i, ii;
	long l;
  // ****************** Measurements ****************** 

  // Correlators here
  for (long jj = 0; jj < v0/4; ++jj) {
  		ii = jj*4;
      //if (ii % 4 == 0) {
      l = lpos[jj];
      i = long(l / Np);
      if (flag[i] == 0) {
        o = str[0][i] * str[1][i];
        flag[i] = 1;
        if (o != 0) {
          k = str[2][i];
          //<H-Bond.H-Bond>
          if (k == 0) { // k = 0 --> H, 1--> QoQ, 2--> Bi
            b = str[3][i] - 1;
            vb2 = long(b / (4*Ns));
            b3 = long(b / 4);
            b5 = long((b % (4*Ns)) / 4);
            l = lpos[long(((ii + 4) % v0) / 4)];
            j = long(l / Np) % Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 0) { // k = 0 --> H, 1--> QoQ, 2--> Bi
                b = str[3][j] - 1;
                vb1 = long(b / ( 4*Ns ));
                b2 = long(b / 4);
                b4 = long((b % ( 4*Ns )) / 4);
                if (vb1 == 0 and vb2 == 0) {
                	//std::cout << b3 << " H " << b2  << std::endl; 
                  Cx[b2][b3] += n1-1;
                  if (b2 != b3) Cx[b3][b2] += n1-1;
                }
                if (vb1 == 1 and vb2 == 1) {
                	////std::cout << b5 << " V " << b4  << std::endl;                   
                  Cy[b5][b4] += n1-1;
                  if (b5 != b4) Cy[b4][b5] += n1-1;
                }
              }
            }
          }
          //<Biquad.Biquad>
          else if (k == 2) {
            b = str[3][i] - 1;
            vb2 = long(b / (2*Ns ));
            b3 = long(b / 2);
            b5 = long((b % (2*Ns )) / 2);
            //B[int(b / 2)] += 1.;
            l = lpos[long(((ii + 8) % v0) / 4)];
            j = long(l / Np) % Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 2) { // k = 0 --> H, 1--> QoQ, 2--> Bi
                b = str[3][j] - 1;
                vb1 = long(b / (2*Ns ));
                b2 = long(b / 2);
                b4 = long((b % (2*Ns )) / 2);
                if (vb1 == 0 and vb2 == 0) { 	// Along x. So, b2, b3 are horizontal bonds
                  Dx[b2][b3] += n1-1;
                  if (b2 != b3) Dx[b3][b2] += n1-1;
                }
                if (vb1 == 1 and vb2 == 1) {	// Along y. So, b4, b5 are vertical bonds
                  Dy[b4][b5] += n1-1;
                  if (b4 != b5) Dy[b5][b4] += n1-1;
                }
              }
            }
          }
        }
      }
      //}
  }
  
  /*
  for (int j =0; j<Ns; j++){
		for (int i =j; i<Ns; i++){
			Cxx[i][j] += (n1-1)*Cx[i][j];
			Dxx[i][j] += (n1-1)*Dx[i][j];		
			Cyy[i][j] += (n1-1)*Cy[i][j];
			Dyy[i][j] += (n1-1)*Dy[i][j];		
			Cx[i][j] = 0;Cy[i][j] = 0;Dx[i][j] = 0;Dy[i][j] = 0;
		}
  }
  */
  
  // Update part
  // Projector part
  for (int i = 0; i < Np; ++i) {
    o = pstr[i];
    if (o != 0) {
      s1 = i;
      s2 = i + Np;
      lpos[int(v0 / 4)] = Np * (i + Lc);
      v1 = last[s1];
      v2 = last[s2];

      if (v1 != -1) {
        X[v1] = v0;
        X[v0] = v1;
      } else {
        frst[s1] = v0;
      }
      if (v2 != -1) {
        X[v2] = v0 + 1;
        X[v0 + 1] = v2;
      } else {
        frst[s2] = v0 + 1;
      }
      last[s1] = v0 + 2;
      last[s2] = v0 + 3;
      v0 += 4;
    }
  }


 
  

  // PBC, can be exploited to calculate observables.
  for (int j = 0; j < Nm; ++j) {
    v1 = frst[j];
    if (v1 != -1) {
      v2 = last[j];
      X[v2] = v1;
      X[v1] = v2;
    }
  }
  

    
  // Main loop 
  for (int j = 0; j < nl; j++) {
    vi = int(rann() * v0);
    v1 = vi;
    l = lpos[long(v1 / 4)];
    i = long(l / Np);
    while (1) {
      ic = v1 % 4;
      if (i < Lc) {
        oc = ic ^ 1;
        k = str[2][i];
        if (k == 0) {
          str[0][i] = (1 + 2) - str[0][i];
        } else if (k == 1) {
          v = long((l % Np) / 2); // Index of which bond out of those 4.
          str[v][i] = flip[str[v][i]][(l % Np) % 2];
        } else if (k == 2) {
        	v = (l % Np) % 2;
          str[0][i] = flip[str[0][i]][v];
        }
      } else {
        r = rann();
        vx = pstr[i - Lc];
        p = 0;
        for (oc = 0; oc < 4; ++oc) {
          p += vxprb[oc][ic][vx];
          if (r < p) {
            pstr[i - Lc] = vxnew[oc][ic][vx];
            break;
          }
        }
      }
      v2 = 4 * (long(v1 / 4)) + oc;
      v = X[v2];
      if (v == vi) break;
      else {
        v1 = v;
        l = lpos[long(v1 / 4)];
        i = long(l / Np);
      }
    }
  }

  // Update spin configuration here
  for (int j = 0; j < Nm; ++j) {
    // Valid only where there is projector operators.
    if (last[j] != -1) {
      int i = long(lpos[long(last[j] / 4)] / Np) - Lc;
      o = pstr[i];
      b = last[j] % 4;
      lattice[j].set_S(2 * vxleg[b][o] - 1);
    } else {
      if (rann() < 0.5) lattice[j].flip();
    }
    frst[j] = -1;
    last[j] = -1;
  }

  delete[] lpos;
  delete[] X;
  delete[] flag;
  
  
  
}


/*
				if (int(b/(2*Ns)) == 0 and  int(b1/(2*Ns)) == 0){
           Cx[(b - b%4)/2][(b1 - b1%4)/2] += 1.;
           if ((b - b%4)/2 != (b1 - b1%4)/2){
           		Cx[(b1 - b1%4)/2][(b - b%4)/2] += 1.;
           }
        }
        if (int(b1/(2*Ns)) == 1 and  int(b/(2*Ns)) == 1){
           Cy[int((b1% (2*Ns))  / 4)][int((b%  (2*Ns))  / 4)] += 1.;
           if (int((b1% (2*Ns))  / 4) != int((b% (2*Ns))  / 4)){ 
	          	Cy[int((b%  (2*Ns))  / 4)][int((b1% (2*Ns))  / 4)] += 1.;
	         }
				} 
*/


        //if (lattice[ss1].S()+lattice[ss2].S() == 0) {
        /*
        if (int(b / (4 * Ns)) == 0){
				jjx += lattice[ss2].S();
        }
        else if (int(b / (4 * Ns)) == 1) {
          jjy += lattice[ss2].S();
        }
        */
        //}
        /*        
        else {
          jjx += lattice[ss2].S()*(xy(1,s1)+xy(1,s2));
          jjy += lattice[ss2].S()*(xy(2,s1)+xy(2,s2));
        }
        */
        //jj[(ss2%Ns)/Ns] += (lattice[ss2].S());
        //jj[(ss2%Ns)/Ns] += (lattice[ss2].S());
        
        
/*

		          for (int z=0; z<2; z++){
				        ss1 = JBsites[b][z];
				        lattice[ss1].flip();
				          V = (lattice[ss1].S());
					x = psite[ss1 % Ns][0];
				        	y = psite[ss1 % Ns][1];
					        
					W1 += cos(x * (pi + 2. * pi / lx) + y * pi) * V;
					W2 += cos(x * pi + y * (pi + 2. * pi / lx)) * V;
					W3 += cos(x * pi + y * pi) * V;
					d *= ss1%Ns;				        

				       	if (z%2==0) mg += 4. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);			          
			}
			if (d == 0){
				V = ((lattice[0].S() + lattice[0 + Ns].S())/2.);
				N_real3 = 0.;//-N_real3 + V*W3;
				N_real2 = 0.;//-N_real2 + V*W2;
				N_real1 = 0.;//-N_real1 + V*W1;
			}
			else {
				V = ((lattice[0].S() + lattice[0 + Ns].S())/2.); 
				N_real3 = N_real3 + V*W3;
				N_real2 = N_real2 + V*W2;
				N_real1 = N_real1 + V*W1;
			}

*/


/*


				x = psite[ss1 % Ns][0];
				y = psite[ss1 % Ns][1];

				V = 2.*(lattice[ss1].S());//+lattice[ss1%Ns].S());
				N_real1 += cos(x * (pi + 2. * pi / lx) + y * pi) * V;
				N_real2 += cos(x * pi + y * (pi + 2. * pi / lx)) * V;
				N_real3 += cos(x * pi + y * pi) * V;

				N_imag1 += sin(x * (pi + 2. * pi / lx) + y * pi) * V;
				N_imag2 += sin(x * pi + y * (pi + 2. * pi / lx)) * V;
				N_imag3 += sin(x * pi + y * pi) * V;
*/    


		      /*
		      	if (fm[b] == -2){
			      	for(int x = 0; x < Ns; x++){
			      		By[x] =  (lattice[0].S() + lattice[Ns].S())/2.*(lattice[x].S() + lattice[Ns+x].S())/2.;
			      	}
				for (int z=0; z<4; z++){
					ss1 = JBsites[b][z];
					lattice[ss1].flip();
					if (z%2==0) mg += 4. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
				}
			      	for(int x = 0; x < Ns; x++){
			      		//V = Bx[x];
			      		Bx[x] += By[x] - (lattice[0].S() + lattice[Ns].S())/2.*(lattice[x].S() + lattice[Ns+x].S())/2.; 
			      		//Bx[x] += V;	
			      	}
		      	}
		      	else {
			      	W2 = (lattice[0].S() + lattice[Ns].S())/2.;
				for (int z=0; z<4; z++){
					ss1 = JBsites[b][z];
					lattice[ss1].flip();
					//V = Bx[ss1%Ns];
					Bx[ss1%Ns] += W2*lattice[ss1].S();
					//Bx[ss1%Ns] += V;
					
					if (z%2==0) mg += 4. * (lattice[ss1].S()) * pow(-1, ((ss1 % Ns) % lx) + (ss1 % Ns) / lx);	
			}
			}
			*/		



/*

  //double W2; //W1, W2, W3;
  //int d;
  for (int j = 0; j < Ns; ++j) {
    //Bx[j] += (lattice[j].S() + lattice[j + Ns].S());
    V =  (lattice[j].S() + lattice[j + Ns].S())/2.;
    x = psite[j][0];
    y = psite[j][1];
    N_real1 += cos(x * (pi + 2. * pi / lx) + y * pi) * V;
    N_real2 += cos(x * pi + y * (pi + 2. * pi / lx)) * V;
    N_real3 += cos(x * pi + y * pi) * V;

    //std::cout << "site " <<  j  << "  " << cos(x * pi + y * pi)  << "  " << V << "  " << N_real3 << std::endl;

    N_imag1 += sin(x * (pi + 2. * pi / lx) + y * pi) * V;
    N_imag2 += sin(x * pi + y * (pi + 2. * pi / lx)) * V;
    //N_imag3 += sin(x * pi + y * pi) * V;

    mag += (lattice[j].S() + lattice[j + Ns].S()) * pow(-1, (j % lx) + j / lx);
  }
  O_N_abs  = pow(N_real3, 2);
  O_N1_abs = (pow(N_real1, 2) + pow(N_imag1, 2));
  O_N2_abs = (pow(N_real2, 2) + pow(N_imag2, 2)); 
  //for (int j = 0; j < Ns; j++) {
    //Bx[j] = (lattice[0].S() + lattice[Ns].S())/2.*(lattice[j].S() + lattice[j + Ns].S())/2.;
    //B[j] += Bx[j];
    //cout << Bx[j] << "  " << (lattice[0].S() + lattice[Ns].S())/2. << "  " << (lattice[j].S() + lattice[j + Ns].S())/2. << endl;
  //}
  

*/			 


/*

      if (flag[i] == 0) {
        o = str[0][i] * str[1][i];
        flag[i] = 1;
        if (o != 0) {
          k = str[2][i];
          //<H-Bond.H-Bond>
          if (k == 0) { // k = 0 --> H, 1--> QoQ, 2--> Bi
            b = str[3][i] - 1;
            b1 = b;
            vb2 = int(b1 / (4 * Ns));
            b3 = int(b1 / 2);
            b5 = int((b1 % (4 * Ns)) / 2);
            //if (vb2 == 0) Bx[b3] += 1.;
            //if (vb2 == 1) By[b5] += 1.;
            l = lpos[int(((ii + 4) % v0) / 4)];
            j = int(l / Np) % Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 0) { // k = 0 --> H, 1--> QoQ, 2--> Bi
                b = str[3][j] - 1;
                vb1 = int(b / (4 * Ns));
                b2 = int(b / 2);
                b4 = int((b % (4 * Ns)) / 2);
                if (vb1 == 0 and vb2 == 0) {
                  Cx[b2][b3] += 1.;
                  if (b2 != b3) Cx[b3][b2] += 1.;
                }
                if (vb1 == 1 and vb2 == 1) {
                  Cy[b5][b4] += 1.;
                  if (b5 != b4) Cy[b4][b5] += 1.;
                }
              }
            }
          }
          //<Biquad.Biquad>
          else if (k == 2) {
            b = str[3][i] - 1;
            b1 = b;
            vb2 = int(b1 / (2 * Ns));
            b3 = int(b1 / 2);
            b5 = int((b1 % (2 * Ns)) / 2);
            //B[int(b / 2)] += 1.;
            l = lpos[int(((ii + 8) % v0) / 4)];
            j = int(l / Np) % Lc;
            o = str[0][j] * str[1][j];
            if (o != 0) {
              k = str[2][j];
              if (k == 2) { // k = 0 --> H, 1--> QoQ, 2--> Bi
                b = str[3][j] - 1;
                vb1 = int(b / (2 * Ns));
                b2 = int(b / 2);
                b4 = int((b % (2 * Ns)) / 2);
                if (vb1 == 0 and vb2 == 0) { // Along x. So, b2, b3 are horizontal bonds
                  Dx[b2][b3] += 1.;
                  if (b2 != b3) Dx[b3][b2] += 1.;
                }
                if (vb1 == 1 and vb2 == 1) {	// Along y. So, b4, b5 are vertical bonds
                  Dy[b4][b5] += 1.;
                  if (b4 != b5) Dy[b5][b4] += 1.;
                }
              }
            }
          }
        }
      }   
      
      */
