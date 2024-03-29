

int j,abs_imj;
double a, b, eta, eta2, eta3, eta4,eta5,eta6;
double oeta,oeta2,oeta3,oeta4,oeta5,oeta6;
double opieta, opieta2, opieta3, opieta4;
double o9p9eta3,o9p9eta4, o9p9eta5,o9p9eta6;
double n0202, n020202, n1212, n121212;
double n00,n01,n02,n03;
double n10,n11,n12,n13;


for (int i=0; i<N; ++i){

	n00 = n(0,0,i); n01 = n(0,1,i);
	n02 = n(0,2,i); n03 = n(0,3,i);
	n0202 = n02*n02;
	n020202 = n0202*n02;

   for (int imj=-Nko2; imj<Nko2; ++imj) {

 
        j = (i + imj + N) % N;
        abs_imj = abs(imj);
		n10 = n(1,0,j); n11 = n(1,1,j);
		n12 = n(1,2,j); n13 = n(1,3,j);

        eta = 1.0 - n03 - n13;

		n1212 = n12*n12;
		n121212 = n1212*n12;

        eta2 = eta*eta;
        eta3 = eta2*eta;
        eta4 = eta3*eta;                   
        eta5 = eta4*eta;
        eta6 = eta5*eta;

        oeta = 1.0/eta;
        oeta2 = 1.0/eta2;
        oeta3 = 1.0/eta3;
        oeta4 = 1.0/eta4;
        oeta5 = 1.0/eta5;
        oeta6 = 1.0/eta6;

        opieta = 0.31830988618379069*oeta;
        opieta2 = 0.31830988618379069*oeta2;
        opieta3 = 0.31830988618379069*oeta3;
        opieta4 = 0.31830988618379069*oeta4;

        o9p9eta3 = 0.10132118364233778*oeta3;
        o9p9eta4 = 0.10132118364233778*oeta4;
        o9p9eta5 = 0.10132118364233778*oeta5;
        o9p9eta6 = 0.10132118364233778*oeta6;


           
        a = K(0,0,abs_imj); 

//        if(a!=0.0)
{



        	p(0,0,i) += a*(n10*oeta + 0.0833333333333333*(n121212)*opieta3 + n11*n12*oeta2);

        	p(0,1,i) += a*(0.25*n02*(n121212)*opieta4 + 2.0*n11*n02*n12*oeta3+ n10*n02*oeta2);

        	p(0,2,i) += a*((0.75*n11*n12*(n0202)*0.31830988618379069 + 0.25*n01*(n121212)*0.31830988618379069)*oeta4+ (2.0*n01*n11*n12 + 0.25*n10*(n0202)*0.31830988618379069)*oeta3+ 0.125*(n0202)*(n121212)*o9p9eta5 + n10*n01*oeta2);

        	p(0,3,i) += a*((n01*n02*(n121212)*0.31830988618379069 + n11*n12*(n020202)*0.31830988618379069)*oeta5+ (6.0*n01*n11*n02*n12 + 0.25*n00*(n121212)*0.31830988618379069 + 0.25*n10*(n020202)*0.31830988618379069)*oeta4+ (2.0*n00*n11*n12 + 2.0*n10*n01*n02)*oeta3+ 0.208333333333333*(n020202)*(n121212)*o9p9eta6 + n00*n10*oeta2);

        	p(1,0,j) += a*(n00*oeta + 0.0833333333333333*(n020202)*opieta3 + n01*n02*oeta2);
			 
        	p(1,1,j) += a*(0.25*n12*(n020202)*opieta4 + 2.0*n01*n02*n12*oeta3+ n00*n12*oeta2);

        	p(1,2,j) += a*((0.75*n01*n02*(n1212)*0.31830988618379069 + 0.25*n11*(n020202)*0.31830988618379069)*oeta4+ (2.0*n01*n11*n02 + 0.25*n00*(n1212)*0.31830988618379069)*oeta3+ 0.125*(n020202)*(n1212)*o9p9eta5 + n00*n11*oeta2);

        	p(1,3,j) += a*((n01*n02*(n121212)*0.31830988618379069 + n11*n12*(n020202)*0.31830988618379069)*oeta5+ (6.0*n01*n11*n02*n12 + 0.25*n00*(n121212)*0.31830988618379069 + 0.25*n10*(n020202)*0.31830988618379069)*oeta4+ (2.0*n00*n11*n12 + 2.0*n10*n01*n02)*oeta3+ 0.208333333333333*(n020202)*(n121212)*o9p9eta6 + n00*n10*oeta2);

        }


            
        a = K(0,1,abs_imj); 

//        if(a!=0.0)
{



        	p(0,0,i) += a*(n11*oeta + 0.125*(n1212)*opieta2);

        	p(0,1,i) += a*((0.25*n02*(n1212)*0.31830988618379069 + 0.0833333333333333*(n121212)*0.31830988618379069)*oeta3+ (n11*n02 + n11*n12)*oeta2+ n10*oeta);

        	p(0,2,i) += a*((0.09375*(n0202)*(n1212)*0.10132118364233778+ 0.0625*(n02)*(n121212)*0.10132118364233778)*oeta4+ (0.5*n11*n12*(n02)*0.31830988618379069 + 0.25*n01*(n1212)*0.31830988618379069 + 0.25*n11*(n0202)*0.31830988618379069)*oeta3+ (n01*n11 + 0.25*n10*(n02)*0.31830988618379069)*oeta2);

        	p(0,3,i) += a*((0.125*(n020202)*(n1212)*0.10132118364233778+ 0.125*(n0202)*(n121212)*0.10132118364233778)*oeta5+ (0.75*n01*n02*(n1212)*0.31830988618379069 + 0.75*n11*n12*(n0202)*0.31830988618379069 + 0.25*n01*(n121212)*0.31830988618379069 + 0.25*n11*(n020202)*0.31830988618379069)*oeta4+ (2.0*n01*n11*n02 + 2.0*n01*n11*n12 + 0.25*n00*(n1212)*0.31830988618379069 + 0.25*n10*(n0202)*0.31830988618379069)*oeta3+ (n00*n11 + n10*n01)*oeta2);

        	p(1,0,j) += a*(n01*oeta + 0.125*(n0202)*opieta2);

        	p(1,1,j) += a*((0.25*n12*(n0202)*0.31830988618379069 + 0.0833333333333333*(n020202)*0.31830988618379069)*oeta3+ (n01*n02 + n01*n12)*oeta2+ n00*oeta);

        	p(1,2,j) += a*((0.0625*(n020202)*(n12)*0.10132118364233778+ 0.09375*(n0202)*(n1212)*0.10132118364233778)*oeta4+ (0.5*n01*n02*(n12)*0.31830988618379069 + 0.25*n01*(n1212)*0.31830988618379069 + 0.25*n11*(n0202)*0.31830988618379069)*oeta3+ (n01*n11 + 0.25*n00*(n12)*0.31830988618379069)*oeta2);

        	p(1,3,j) += a*((0.125*(n020202)*(n1212)*0.10132118364233778+ 0.125*(n0202)*(n121212)*0.10132118364233778)*oeta5+ (0.75*n01*n02*(n1212)*0.31830988618379069 + 0.75*n11*n12*(n0202)*0.31830988618379069 + 0.25*n01*(n121212)*0.31830988618379069 + 0.25*n11*(n020202)*0.31830988618379069)*oeta4+ (2.0*n01*n11*n02 + 2.0*n01*n11*n12 + 0.25*n00*(n1212)*0.31830988618379069 + 0.25*n10*(n0202)*0.31830988618379069)*oeta3+ (n00*n11 + n10*n01)*oeta2);
        }

            
        a = K(0,2,abs_imj); 

        if(a!=0.0){



        	p(0,0,i) += a*(n12*oeta);

        	p(0,1,i) += a*(n02*n12*oeta2);

        	p(0,2,i) += a*((0.25*n12*(n0202)*0.31830988618379069 + 0.0833333333333333*(n121212)*0.31830988618379069)*oeta3+ (n01*n12 + n11*n12)*oeta2+ n10*oeta);

        	p(0,3,i) += a*((0.25*n02*(n121212)*0.31830988618379069 + 0.25*n12*(n020202)*0.31830988618379069)*oeta4+ (2.0*n01*n02*n12 + 2.0*n11*n02*n12)*oeta3+ (n00*n12 + n10*n02)*oeta2);

        	p(1,0,j) += a*(n02*oeta);

        	p(1,1,j) += a*(n02*n12*oeta2);

        	p(1,2,j) += a*((0.25*n02*(n1212)*0.31830988618379069 + 0.0833333333333333*(n020202)*0.31830988618379069)*oeta3+ (n01*n02 + n11*n02)*oeta2+ n00*oeta);

        	p(1,3,j) += a*((0.25*n02*(n121212)*0.31830988618379069 + 0.25*n12*(n020202)*0.31830988618379069)*oeta4+ (2.0*n01*n02*n12 + 2.0*n11*n02*n12)*oeta3+ (n00*n12 + n10*n02)*oeta2);
        }

            
        a = K(0,3,abs_imj); 

        if(a!=0.0){



        	p(0,0,i) += a*(-log(eta));

        	p(0,1,i) += a*(n02*oeta);

        	p(0,2,i) += a*(n01*oeta + 0.125*(n0202)*opieta2);

        	p(0,3,i) += a*((n00 + n10)*oeta + (0.0833333333333333*(n020202)*0.31830988618379069 + 0.0833333333333333*(n121212)*0.31830988618379069)*oeta3+ (n01*n02 + n11*n12)*oeta2);

        	p(1,0,j) += a*(-log(eta));

        	p(1,1,j) += a*(n12*oeta);

        	p(1,2,j) += a*(n11*oeta + 0.125*(n1212)*opieta2);

        	p(1,3,j) += a*((n00 + n10)*oeta + (0.0833333333333333*(n020202)*0.31830988618379069 + 0.0833333333333333*(n121212)*0.31830988618379069)*oeta3+ (n01*n02 + n11*n12)*oeta2);
        }

            
        a = K(1,1,abs_imj); 

//        if(a!=0.0)
{



        	p(0,1,i) += a*(n11*oeta + 0.125*(n1212)*opieta2);

        	p(0,2,i) += a*(0.0625*(n02)*(n1212)*o9p9eta3 + 0.25*n11*(n02)*opieta2);

        	p(0,3,i) += a*((0.25*n01*(n1212)*0.31830988618379069 + 0.25*n11*(n0202)*0.31830988618379069)*oeta3+ 0.09375*(n0202)*(n1212)*o9p9eta4 + n01*n11*oeta2);

        	p(1,1,j) += a*(n01*oeta + 0.125*(n0202)*opieta2);

        	p(1,2,j) += a*(0.0625*(n0202)*(n12)*o9p9eta3 + 0.25*n01*(n12)*opieta2);

        	p(1,3,j) += a*((0.25*n01*(n1212)*0.31830988618379069 + 0.25*n11*(n0202)*0.31830988618379069)*oeta3+ 0.09375*(n0202)*(n1212)*o9p9eta4 + n01*n11*oeta2);
        }

            
        a = K(1,2,abs_imj); 

        if(a!=0.0){



        	p(0,1,i) += a*(n12*oeta);

        	p(0,2,i) += a*((0.25*n12*(n02)*0.31830988618379069 + 0.125*(n1212)*0.31830988618379069)*oeta2+ n11*oeta);

        	p(0,3,i) += a*((0.25*n02*(n1212)*0.31830988618379069 + 0.25*n12*(n0202)*0.31830988618379069)*oeta3+ (n01*n12 + n11*n02)*oeta2);

        	p(1,1,j) += a*(n02*oeta);

        	p(1,2,j) += a*((0.25*n02*(n12)*0.31830988618379069 + 0.125*(n0202)*0.31830988618379069)*oeta2+ n01*oeta);

        	p(1,3,j) += a*((0.25*n02*(n1212)*0.31830988618379069 + 0.25*n12*(n0202)*0.31830988618379069)*oeta3+ (n01*n12 + n11*n02)*oeta2);
        }

            
        a = K(1,3,abs_imj); 

        if(a!=0.0){
			b = a*(-log(eta));


        	p(0,1,i) += b;

        	p(0,2,i) += a*(0.25*(n02)*opieta);

        	p(0,3,i) += a*((n01 + n11)*oeta + (0.125*(n0202)*0.31830988618379069 + 0.125*(n1212)*0.31830988618379069)*oeta2);

        	p(1,1,j) += b;

        	p(1,2,j) += a*(0.25*(n12)*opieta);

        	p(1,3,j) += a*((n01 + n11)*oeta + (0.125*(n0202)*0.31830988618379069 + 0.125*(n1212)*0.31830988618379069)*oeta2);
        }

            
        a = K(2,2,abs_imj); 

        if(a!=0.0){



        	p(0,2,i) += a*(n12*oeta);

        	p(0,3,i) += a*(n02*n12*oeta2);

        	p(1,2,j) += a*(n02*oeta);

        	p(1,3,j) += a*(n02*n12*oeta2);
        }

            
        a = K(2,3,abs_imj); 

        if(a!=0.0){
			b = a*(-log(eta));


        	p(0,2,i) += b;

        	p(0,3,i) += a*((n02 + n12)*oeta);

        	p(1,2,j) += b;

        	p(1,3,j) += a*((n02 + n12)*oeta);
        }

            
        a = K(3,3,abs_imj); 

        if(a!=0.0){
			b = a*(-log(eta));


        	p(0,3,i) += b;

        	p(1,3,j) += b;
	
        }


	}    // imj_loop
                
}     // i_loop
