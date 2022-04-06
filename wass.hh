#ifndef WASS_H
#define WASS_H

#include "misfit.hh"
#include <iostream>
#include <cassert>

using namespace std;
template<class T> using Ctn=vector<T>;


template< class T >
class Wass : public Misfit<T>{
public:
    Ctn<T> cdf(Ctn<T> f, Ctn<T> x){
        //precondition error checking
        assert( f.size() == x.size() );
    
        int N = f.size();
        Ctn<T> F(0,N);
        for(int i = 0; i < N - 1; i++){
            std::cout << i << " ";
            T df = f[i+1] - f[i];
            T dx = x[i+1] - x[i];
            assert( dx > 0 );
            F[i+1] = F[i] + df/dx;
            std::cout << i+1 << endl;
        }
        return F;
    }
    
    Ctn<T> quantile(Ctn<T> F, 
        Ctn<T> x, 
        Ctn<T> p, 
        T eps){
         //precondition error checking
        // assert( F.size() == x.size() )
         //looping vars
         int P = p.size(); 
         int N = x.size();
         int i_x = 0;
    
         //initialize Quantile
         Ctn<T> Q(0, P);
    
         for(int i_p = 0; i_p < P; i_p++){
             while( i_x < N ) {
                 if( F[i_x] <= p[i_p] and p[i_p] <= F[i_x+1] ){
                     T df = F[i_x+1] - F[i_x];
                     if( df < eps ){
                         //left bin -- division might cause problems
                         Q[i_p] = F[i_x];
                     }
                     else{ //else interpolate between left and right
                         Q[i_p] = (p - F[i_x]) / df;
                     }
                     if( i_x != N-1 ) i_x++;
                 }
             }
         }
         return Q;
    }
    
    T eval(Ctn<T> m){
        std::cout << "Wasserstein dummy eval\n";
        return 0.0;
    }
};

#endif
