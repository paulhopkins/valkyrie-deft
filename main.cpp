#include<iostream>
#include<fstream>
#include<cmath>

#include"nahs_planar.hpp"

using namespace std;

int main(int argc, char* argv) {
	int N = 2048;
	double dr = 0.02;
	double R[2] = {0.5,1.0};
	double mu[2] = {0.0,0.0};
	double alpha = 0.08;

	nahs_planar *myDFT = new nahs_planar(N, dr, R, 0.2);

	cout << "Initialised"<< endl; 

	double rhob[2] = {0.01,0.01};
	rhob[0] = 4.223539288213692364e-05;
	rhob[1] = 1.275666657617313959e-01;

//	rhob[0] = 8.195416335694483134e-03;
//	rhob[1] = 9.004472879821673326e-02;
//	rhob[0] = 8.312891860807745600e-03;
//	rhob[1] =  8.984551185446439103e-02;
//	rhob[0] = 1.546224012985698837e-02;
//	rhob[1] =  8.033441245529272801e-02;

	myDFT->set_rho(rhob);

	cout << "rho set";

	myDFT->calc_c1();

	cout << "c1 calced" << endl;


	ofstream newfile("results/c1.dat");	
	
	for (int i =0; i < N; ++i) {
		newfile << i << " " << myDFT->get_c1(0,i) << " " << 
			myDFT->get_c1(1,i) << " " << myDFT->get_rho(0,i) << " " << 
				myDFT->get_rho(1,i) << endl;
	}


	newfile.close();





	int left[2], right[2];
	double rho_new[2*N], c1_max;

	left[0] = myDFT->get_left(0);
	right[0] = myDFT->get_right(0);
	left[1] = myDFT->get_left(1);
	right[1] = myDFT->get_right(1);


	mu[0] = log(rhob[0]) + myDFT->get_c1(0,N/2);
	mu[1] = log(rhob[1]) + myDFT->get_c1(1,N/2);

	cout << "Omega" << (myDFT->calc_F()/N*dr) << endl;


//	myDFT->load_rho();

	cout << "alldone" << endl;

	myDFT->hardwall();

	for (unsigned int loop = 0; loop < 1<<30; ++loop) {

		for (unsigned int inner = 0; inner < 21; ++inner) {
			myDFT->calc_c1();
			c1_max = 0.0;

			for (int s = 0; s < 2; ++s) {
				for (int i = left[s]; i < right[s]; ++i) {

				rho_new[s*N + i] = mu[s]-myDFT->get_c1(s,i);			
				c1_max = fabs(rho_new[s*N + i] - log(myDFT->get_rho(s,i)) ) > c1_max ? fabs(rho_new[s*N + i]-log(myDFT->get_rho(s,i)) ) : c1_max;
				rho_new[s*N + i] = exp(rho_new[s*N + i]);

				myDFT->set_rho(s,i) = (1.0-alpha)*myDFT->get_rho(s,i) 
		+ alpha*rho_new[s*N + i];
				}
			}
			
		}

		cout << loop << " Error: " << c1_max << endl;

		newfile.open("results/rho_cur.dat");		
		for (int i =0; i < N; ++i) {
			newfile << (i-left[0])*dr << " " << myDFT->get_rho(0,i) << " " << 
				myDFT->get_rho(1,i) << " " << 
				mu[0]-myDFT->get_c1(0,i) - log(myDFT->get_rho(0,i)) << " "
				<< mu[1]-myDFT->get_c1(1,i)-log(myDFT->get_rho(1,i)) << endl;
		}
		newfile.close();
	}
			
	
	myDFT->set_mu(mu);



	cout << "hardwall";

	

/*
	cout << "running" << myDFT->get_rho()[0] << endl;
   	cout << myDFT->get_F() << endl;
*/
	cout << endl;
	delete myDFT;

}
