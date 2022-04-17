#ifndef WASS_SPLIT_H
#define WASS_SPLIT_H

#include "wass.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSplit : public Wass<T>{
public:
    //inherit superclass constructor
    using Wass<T>::Wass;

    //implement virtual renormalization routine from superclass
    valarray<valarray<T> > renormalize(valarray<T> f){
        //splitting maps to two pdfs
        valarray<valarray<T> > g(f,2);

        //tmp var for efficiency
        int N = f.size();
 
        assert( N == this->get_nt() );
       
        T pos_int = 0.0;
        T neg_int = 0.0; 
        //split between positive and negative parts
        for(int i = 0; i < N; i++){
            if( g[0][i] < 0 ) g[0][i] = 0.0;
            if( g[1][i] > 0 ) g[1][i] = 0.0;
            if( i >= 1 ){
                T dt = this->t[i] - this->t[i-1];
                pos_int += 0.5 * dt * (g[0][i] + g[0][i-1]);
                neg_int += 0.5 * dt * (g[1][i] + g[1][i-1]);
            }
        }
  
        T tol=1e-20;
        if( pos_int >= tol ){
            T pos_int_inv = 1.0 / pos_int;
            for(int i = 0; i < N; i++) g[0][i] = pos_int_inv * g[0][i];
        }
        else{
            for(int i = 0; i < N; i++) g[0][i] = 0.0;
        }

        if( neg_int >= tol ){
            T neg_int_inv = 1.0 / neg_int;
            for(int i = 0; i < N; i++) g[1][i] = neg_int_inv * g[1][i];
        }
        else{
            for(int i = 0; i < N; i++) g[1][i] = 0.0;
        }
        return g;
    }
};

#endif
