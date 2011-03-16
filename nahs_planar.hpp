#ifndef NAHS_PLANAR_H
#define NAHS_PLANAR_H

#include <iostream>
#include <fftw3.h>


class nahs_planar {

public:
	nahs_planar (int N, double dr, double* R, double Delta);


public:
	void set_mu(double* mu);
	double* get_mu();

	void set_rho(double* rhob);

	double* calc_n();

	double* calc_c1();
	double get_c1(int s, int i) {return c1[s*N + i];}

	double calc_F();
	double get_F();

	void hardwall();
	void load_rho();

	double *get_rho(int s) {return &rho[s*N];}
	const double get_rho(int s, int i) {return rho[s*N + i];}


	double n(int s, int g, int i) {return n_[(s*4 + g)*N +i];}

	int get_left(int s) {return left[s];}
	int get_right(int s) {return right[s];}

	double &set_rho(int s,int i) {return rho[s*N + i];}

	
private:

	double &p(int s, int g, int i) {
//return p_[(s*4 + g)*N +i]
return p_[i*8 + s*4 + g]
;}

	double &omega(int i) {return omega_[i];}

	double &K(int a, int b, int i) {
		return K_[(a + b*(b+1)/2)*NK + i];
	}
	
	double &K(int j, int i) {
		return K_[j*NK + i];
	}

	int left[2], right[2];	

private:
	int N, NK;
	double dr, Delta, R12, F, Omega;

	double R[2], mu[2];


	double* r;
	double* rho;
	double* c1;
	fftw_complex* rhok;

	double* wr, *wk;
	double* n_;
	double* p_;
	double* K_;

	double* omega_;

	double *FFTW_dble;
	fftw_complex *FFTW_cmplx;

	fftw_plan FFTWplanF, FFTWplanB;


};
	
#endif
