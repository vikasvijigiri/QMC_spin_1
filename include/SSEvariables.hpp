#ifndef _SSEVARIABLES_HPP_DEFINED_
#define _SSEVARIABLES_HPP_DEFINED_
#include <random>
#include <ctime>
#include <chrono>
using namespace std::chrono;

#pragma once

extern	int lx, ly, Ns, Nm, Np, Nb[3], np1, S;
extern	double J[3];
extern long n1, Lc;
extern	int nbins, isteps, iter;
extern double Beta;
extern int **JHsites, **JQsites, **JBsites,  *JBsgnx, *JBsgny, *JHsgnx, *JHsgny, *JQsgnx, *JQsgny;
extern int pf, tf, pf1, tf1;


class SSEvariables {
	public:
      	void lattice_sites( );
      	void declare_variables();
      	//void set_temperatures();

};

class ran {
  private:
    std::mt19937 mt;
    std::uniform_real_distribution < double > dist;

  public:
    ran(double lower = 0.0, double upper = 1.0): mt(std::chrono::high_resolution_clock::now().time_since_epoch().count()), dist(lower, upper) {}
  double operator()() {
    return dist(mt);
  }
};
#endif
