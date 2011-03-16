#line 55 "planar_nahs.py" 
        int j;
        double a, eta, b;
        double sum = 0.0;
            
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(0,0, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,0,i)*n(1,0,j)/(eta) + 0.0416666666666667*(n(0,2,i)*n(0,2,i)*n(0,2,i))*(n(1,2,j)*n(1,2,j)*n(1,2,j))/((9.869604401089358)*((eta)*(eta)*(eta)*(eta)*(eta))) + 0.25*n(0,1,i)*n(0,2,i)*(n(1,2,j)*n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta)*(eta)*(eta))) + 0.25*n(1,1,j)*n(1,2,j)*(n(0,2,i)*n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta)*(eta)*(eta))) + 0.0833333333333333*n(0,0,i)*(n(1,2,j)*n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta)*(eta))) + 0.0833333333333333*n(1,0,j)*(n(0,2,i)*n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta)*(eta))) + 2*n(0,1,i)*n(1,1,j)*n(0,2,i)*n(1,2,j)/((eta)*(eta)*(eta)) + n(0,0,i)*n(1,1,j)*n(1,2,j)/((eta)*(eta))+ n(1,0,j)*n(0,1,i)*n(0,2,i)/((eta)*(eta)));
                        sum += b;
                        omega(i) += b;
     p(0,0,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(0,1, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,0,i)*n(1,1,j)/(eta) + n(1,0,j)*n(0,1,i)/(eta) + 0.03125*(n(0,2,i)*n(0,2,i))*(n(1,2,j)*n(1,2,j)*n(1,2,j))/((9.869604401089358)*((eta)*(eta)*(eta)*(eta))) + 0.03125*(n(0,2,i)*n(0,2,i)*n(0,2,i))*(n(1,2,j)*n(1,2,j))/((9.869604401089358)*((eta)*(eta)*(eta)*(eta))) + 0.0833333333333333*n(0,1,i)*(n(1,2,j)*n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta)*(eta))) + 0.0833333333333333*n(1,1,j)*(n(0,2,i)*n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta)*(eta))) + 0.25*n(0,1,i)*n(0,2,i)*(n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta)*(eta))) + 0.25*n(1,1,j)*n(1,2,j)*(n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta)*(eta))) + n(0,1,i)*n(1,1,j)*n(0,2,i)/((eta)*(eta))+ n(0,1,i)*n(1,1,j)*n(1,2,j)/((eta)*(eta))+ 0.125*n(0,0,i)*(n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta))) + 0.125*n(1,0,j)*(n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta))));
                        sum += b;
                        omega(i) += b;
     p(0,1,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(0,2, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,0,i)*n(1,2,j)/(eta) + n(1,0,j)*n(0,2,i)/(eta) + 0.0833333333333333*n(0,2,i)*(n(1,2,j)*n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta)*(eta))) + 0.0833333333333333*n(1,2,j)*(n(0,2,i)*n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta)*(eta))) + n(0,1,i)*n(0,2,i)*n(1,2,j)/((eta)*(eta))+ n(1,1,j)*n(0,2,i)*n(1,2,j)/((eta)*(eta)));
                        sum += b;
                        omega(i) += b;
     p(0,2,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(0,3, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(-n(0,0,i)*log(eta) - n(1,0,j)*log(eta) + n(0,1,i)*n(0,2,i)/(eta) + n(1,1,j)*n(1,2,j)/(eta) + 0.0416666666666667*(n(0,2,i)*n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta))) + 0.0416666666666667*(n(1,2,j)*n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta))));
                        sum += b;
                        omega(i) += b;
     p(0,3,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(1,1, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,1,i)*n(1,1,j)/(eta) + 0.03125*(n(0,2,i)*n(0,2,i))*(n(1,2,j)*n(1,2,j))/((9.869604401089358)*((eta)*(eta)*(eta))) + 0.125*n(0,1,i)*(n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta))) + 0.125*n(1,1,j)*(n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta))));
                        sum += b;
                        omega(i) += b;
     p(1,1,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(1,2, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,1,i)*n(1,2,j)/(eta) + n(1,1,j)*n(0,2,i)/(eta) + 0.125*n(0,2,i)*(n(1,2,j)*n(1,2,j))/(3.14159265358979*((eta)*(eta))) + 0.125*n(1,2,j)*(n(0,2,i)*n(0,2,i))/(3.14159265358979*((eta)*(eta))));
                        sum += b;
                        omega(i) += b;
     p(1,2,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(1,3, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(-n(0,1,i)*log(eta) - n(1,1,j)*log(eta) + 0.125*(n(0,2,i)*n(0,2,i))/(3.14159265358979*(eta)) + 0.125*(n(1,2,j)*n(1,2,j))/(3.14159265358979*(eta)));
                        sum += b;
                        omega(i) += b;
     p(1,3,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(2,2, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,2,i)*n(1,2,j)/(eta));
                        sum += b;
                        omega(i) += b;
     p(2,2,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(2,3, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(-n(0,2,i)*log(eta) - n(1,2,j)*log(eta));
                        sum += b;
                        omega(i) += b;
     p(2,3,i) = b;
                        }
                    
                    }
               }
                   
               for (int imj=-Nko2; imj<Nko2; ++imj) {
                   a = K(3,3, abs(imj));
                   if(a!=0.0){
                    for (int i=off0; i<off1; ++i){
                        j = (i - imj + N) % N;
                        eta = 1.0 - n(0,3,i) - n(1,3,j);
                        b = a*(n(0,3,i) + n(1,3,j) + (eta)*log(eta));
                        sum += b;
                        omega(i) += b;
     p(3,3,i) = b;
                        }
                    
                    }
               }
                             
