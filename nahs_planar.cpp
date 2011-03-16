#include<iostream>
#include<fstream>
#include<cmath>

//#include <fftw3.h>

#include"nahs_planar.hpp"

using namespace std;

nahs_planar::nahs_planar(int N, double dr, double R[2], double Delta) {
	double *H, *D0, *D1; // Heavside, Dirac and Dirac' functions.

	// Get raw parameters
	this->N = N;
	this->dr = dr;
	this->R[0] = R[0];
	this->R[1] = R[1];
	this->Delta = Delta;



	// Calc and set derived parameters
	R12 = Delta * (R[0] + R[1]);
	NK = int(2.0*R12/dr) + 10;

	int range = int(2*(R[0]>R[1] ? R[0] : R[1])/dr) + NK;
	left[0] = range;
	left[1] = range + int((R[1] - R[0] + 0.3)/dr);

	right[0] = N - range;
	right[1] = N - range;




	// Allocate arrays
	r = (double*) malloc(N*sizeof(double));

	rho = (double*) malloc(2*N*sizeof(double));
    c1 = (double*) malloc(2*N*sizeof(double));
	rhok = (fftw_complex*) malloc(2*(N/2+1)*sizeof(fftw_complex));

	wr = (double*) malloc(2*4*N*sizeof(double));
	wk = (double*) malloc(2*4*(N/2+1)*sizeof(double));

    n_ = (double*) malloc(2*4*N*sizeof(double));
    p_ = (double*) malloc(2*4*N*sizeof(double));
	
	K_ = (double*) malloc(10*NK*sizeof(double));

	omega_ = (double*) malloc(N*sizeof(double));


	// Define FFTW arrays
	FFTW_dble = (double*) fftw_malloc(N*sizeof(double));
	FFTW_cmplx = (fftw_complex*) fftw_malloc((N/2+1)*sizeof(fftw_complex));

    FFTWplanF = fftw_plan_dft_r2c_1d(N, FFTW_dble, FFTW_cmplx,FFTW_ESTIMATE);
    FFTWplanB = fftw_plan_dft_c2r_1d(N, FFTW_cmplx, FFTW_dble,FFTW_ESTIMATE);


	// Initialise arrays 
    for (int i = 0 ; i < N; ++i) {
		r[i] = dr*i;
	}

	// Set weight functions

	H = (double*) malloc(N*sizeof(double));
	D0 = (double*) malloc(N*sizeof(double));
	D1 = (double*) malloc(N*sizeof(double));



	for (int s = 0; s < 2; ++s) {
		// First make Heaviside and Dirac Functions
		for (int i = 0; i < N; ++i) {
			H[i] = 0.0; D0[i] = 0.0;
			if (r[i] < R[s]) H[i] = 1.0;
			if ( fabs(r[i] - R[s]) < dr/2.0) {
				H[i] = 0.5;
				D0[i] = 1.0/dr;
			}
		} // i loop
		
		// Then make a Dirac p from that function.
		// The gradient loops around so that hopefully R = 0 works ??
		D1[0] = (D0[1]-D0[N-1])/2.0/dr;
		for (int i = 1; i < N-1; ++i) {
			D1[i] = (D0[i+1]-D0[i-1])/2.0/dr;
		}
		D1[N-1] = (D0[0]-D0[N-2])/2.0/dr;

		// Assign the weighted densities.
		for (int i = 0; i < N/2+1; ++i) {
    		wr[s*4*N + i]       = 0.5*D0[i] + 0.25*R[s]*D1[i];
    		wr[(s*4 + 1)*N + i] = 0.25*(H[i] + R[s]*D0[i]);
    		wr[(s*4 + 2)*N + i] = 6.2831853071795862*R[s]*H[i];
    		wr[(s*4 + 3)*N + i] = 3.1415926535897931*H[i]*
									(R[s]*R[s]-r[i]*r[i]);
		}

		// and then mirror
		for (int i = 1; i < N/2; ++i) {
			wr[s*4*N + N-i]       = wr[s*4*N + i];
    		wr[(s*4 + 1)*N + N-i] = wr[(s*4 + 1)*N + i];
    		wr[(s*4 + 2)*N + N-i] = wr[(s*4 + 2)*N + i];
    		wr[(s*4 + 3)*N + N-i] = wr[(s*4 + 3)*N + i];
		}


		for (int j = 0; j < 4; ++j) {
			for (int i = 0; i < N; ++i) {
				FFTW_dble[i] = wr[(s*4 + j)*N + i];
			}
			
			fftw_execute(FFTWplanF);
		
			for (int i = 0; i < (N/2 + 1); ++i) {
				wk[(s*4 + j)*(N/2 + 1) + i] = FFTW_cmplx[i][0]*dr/N;
			}
		}
			
                    
	} // end s loop

	free(H);
	free(D0);
	free(D1);

	double *D2, *D3, *D4; // Heavside, Dirac'' and Dirac''' functions.
	H = (double*) malloc(NK*sizeof(double));
	D0 = (double*) malloc(NK*sizeof(double));
	D1 = (double*) malloc(NK*sizeof(double));
	D2 = (double*) malloc(NK*sizeof(double));
	D3 = (double*) malloc(NK*sizeof(double));
	D4 = (double*) malloc(NK*sizeof(double));

	for (int i = 0; i < NK; ++i) {
		H[i] = 0.0; D0[i] = 0.0;
		if (r[i] < R12) H[i] = 1.0;
		if ( fabs(r[i] - R12) < dr/2.0) {
			H[i] = 0.5;
			D0[i] = 1.0/dr;
		}
	} // i loop
		
	// Then make a Dirac p from that function.
	// The gradient loops around so that hopefully R = 0 works ??
	D1[0] = (D0[1]-D0[NK-1])/2.0/dr;
	for (int i = 1; i < NK-1; ++i) {
		D1[i] = (D0[i+1]-D0[i-1])/2.0/dr;
	}
	D1[NK-1] = (D0[0]-D0[NK-2])/2.0/dr;

	D2[0] = (D1[1]-D1[NK-1])/2.0/dr;
	for (int i = 1; i < NK-1; ++i) {
		D2[i] = (D1[i+1]-D1[i-1])/2.0/dr;
	}
	D2[NK-1] = (D1[0]-D1[NK-2])/2.0/dr;

	D3[0] = (D2[1]-D2[NK-1])/2.0/dr;
	for (int i = 1; i < NK-1; ++i) {
		D3[i] = (D2[i+1]-D2[i-1])/2.0/dr;
	}
	D3[NK-1] = (D2[0]-D2[NK-2])/2.0/dr;

	D4[0] = (D3[1]-D3[NK-1])/2.0/dr;
	for (int i = 1; i < NK-1; ++i) {
		D4[i] = (D3[i+1]-D3[i-1])/2.0/dr;
	}
	D4[NK-1] = (D3[0]-D3[NK-2])/2.0/dr;

	for (int i = 0; i < NK/2; ++i) {
		K_[i]       = 3.1415926*H[i]*(R12*R12-r[i]*r[i]);
        K_[NK + i]   = 6.2831853071795862*R12*H[i];
        K_[3*NK + i] = 0.25*(H[i] + R12*D0[i]);
        K_[6*NK + i] = 0.5*D0[i] + 0.25*R12*D1[i];
 
        K_[2*NK + i] = 6.2831853071795862*(H[i]+R12*D0[i]);
        K_[4*NK + i] = 0.5*D0[i]-0.25*R12*D1[i];
        K_[7*NK + i] = 0.25*(-D1[i]-R12*D2[i]);

        K_[5*NK + i] = 0.0099471839432434591*(-3.0*D1[i]+R12*D2[i]);
        K_[8*NK + i] = 0.0099471839432434591*R12*D3[i];

        K_[9*NK + i] = 0.0099471839432434591*(3.0*D3[i]+R12*D4[i]);
	}

	for (int j = 0; j < 10; ++j) {
		for (int i = 1; i < NK/2; ++i) {
			K(j,NK-i) = K(j,i);
		}
	}

	for (int i = 0; i < 10*NK; ++i) {
		K_[i] *= dr;
	}

	for (int j = 0; j < 10; ++j) {
	double a=0.0;
		for (int i = 0; i < NK; ++i) {
			a += K_[j*NK + i];
		}
	cout << a << endl;
	}

	free(H);
	free(D0);
	free(D1);
	free(D2);
	free(D3);
	free(D4);

	cout << "Just saving..";

	// Output the weight functions in real and Fourier space.
	ofstream newfile("results/wk.dat");	
	for (int i=0; i<(N/2 + 1); i++) {
		newfile << r[i];
			for (int j = 0; j < 8; ++j) {
				newfile << " " << wk[j*(N/2+1) + i];
			}
		newfile << endl;
	} // end i loop
	newfile.close();

	newfile.open("results/wr.dat");	
	for (int i=0; i<N; i++) {
		newfile << r[i];
			for (int j = 0; j < 8; ++j) {
				newfile << " " << wr[j*N + i];
			}
		newfile << endl;
	} // end i loop
	newfile.close();

	newfile.open("results/K.dat");	
	for (int i=0; i<NK; i++) {
		newfile << r[i];
			for (int j = 0; j < 10; ++j) {
				newfile << " " << K_[j*NK + i];
			}
		newfile << endl;
	} // end i loop
	newfile.close();

	cout << "done" << endl;
}



void nahs_planar::set_rho(double* rhob) {
	for (int s = 0; s < 2; ++s) {
		for (int i = 0; i < N; ++i) {
			rho[s*N + i] = rhob[s];
		} // i loop
	} // s loop
}

void nahs_planar::load_rho() {
	double dum;
	ifstream rhofile("results/rho_keep6.dat");		
		for (int i =0; i < right[0]; ++i) {
			rhofile >> dum >> rho[i] >> rho[N + i] >> dum >> dum;
//			cout << rho[i] << " " << rho[N+i] << endl;
		}
	rhofile.close();
}

double* nahs_planar::calc_n() {
	for (int s = 0; s < 2; ++s) {
		for (int i = 0; i < N; ++i) {
			FFTW_dble[i] = rho[s*N + i];
		}
		
		fftw_execute(FFTWplanF);
	
		for (int i = 0; i < (N/2 + 1); ++i) {
			rhok[s*(N/2 + 1) + i][0] = FFTW_cmplx[i][0];
			rhok[s*(N/2 + 1) + i][1] = FFTW_cmplx[i][1];
		}

		for (int j = 0; j < 4; ++j) {
			for (int i = 0; i < N/2 + 1; ++i) {
				FFTW_cmplx[i][0] = wk[(s*4 + j)*(N/2 + 1) + i] * 
									rhok[s*(N/2 + 1) + i][0];
				FFTW_cmplx[i][1] = wk[(s*4 + j)*(N/2 + 1) + i] * 
									rhok[s*(N/2 + 1) + i][1];
			}
			
			fftw_execute(FFTWplanB);
		
			for (int i = 0; i < N; ++i) {
				n_[(s*4 + j)*N + i] = FFTW_dble[i];
			}
		} 
	}

/*	ofstream newfile("results/n.dat");	
	for (int i=0; i<N; i++) {
		newfile << r[i];
			for (int j = 0; j < 8; ++j) {
				newfile << " " << n_[j*N + i];
			}
		newfile << endl;
	} // end i loop
	newfile.close();
*/

	return n_;
}

double nahs_planar::calc_F() {
	F = 0;

	for (int s = 0; s < 2; ++s) {
		for (int i = 0; i < N; ++i) {
			F += rho[s*N + i] * (log(rho[s*N + i]) - 1.0 + mu[s]);
		}
	}

	F *= dr;

	int Nko2 = NK/2;			// Must match python code
	int off0 = left[0];
	int off1 = right[1];

#include"Phi_code.c"
	F += sum;

	
	return F;
}

double nahs_planar::get_F() {
	return F;
}

double* nahs_planar::calc_c1() {
	calc_n();
	
	for (int i = 0; i < 2*N; ++i) {
		c1[i] = 0.0;	// c1 does not actually contain ideal part 
				   		// log(rho[s*N + i]);
	}

	for (int i = 0; i < 8*N; ++i) {
		p_[i] = 0.0;
	}


	int Nko2 = NK/2;			// Must match python code

#include"phi_code.c"


/*	ofstream newfile("results/p.dat");	
	for (int i=0; i<N; i++) {
		newfile << r[i];
			for (int j = 0; j < 8; ++j) {
				newfile << " " << p_[j*N + i];
			}
		newfile << endl;
	} // end i loop
	newfile.close();
*/


	for (int s = 0; s < 2; ++s) {
		for (int j = 0; j < 4; ++j) {
    		for (int i = 0; i < N; ++i) {
    			FFTW_dble[i] = p(s,j,i);
    		}
		
			fftw_execute(FFTWplanF);
	
			for (int i = 0; i < N/2 + 1; ++i) {
				FFTW_cmplx[i][0] *= wk[(s*4 + j)*(N/2 + 1) + i];
				FFTW_cmplx[i][1] *= wk[(s*4 + j)*(N/2 + 1) + i];
			}
			
			fftw_execute(FFTWplanB);
		
			for (int i = 0; i < N; ++i) {
				c1[s*N + i] += FFTW_dble[i];
			}
		} 
	}

	return c1;
}

void nahs_planar::set_mu(double* mu) {
	this->mu[0] = mu[0];
	this->mu[1] = mu[1];
}

double* nahs_planar::get_mu() {
	return mu;
}

void nahs_planar::hardwall() {
	for (int i = 0; i < left[0]; ++i) {
		rho[i] = 0.0;
	}
	for (int i = 0; i < left[1]; ++i) {
		rho[N + i] = 0.0;
	}
}
