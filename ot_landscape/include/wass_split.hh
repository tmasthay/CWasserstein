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
 
        cerr << "N = " << N << endl;
        cerr << "nt = " << this->nt << endl;

        assert( N == this->nt or N == this->num_rec * this->nt );
       
        int num_rec_local = N / this->nt;

        for(int i_r = 0; i_r < num_rec_local; i_r++){
            T pos_int = 0.0;
            T neg_int = 0.0;

            //split between positive and negative parts
            for(int i_t = 0; i_t < this->nt; i_t++){
                int i_gbl = i_t + i_r * this->nt;
                if( g[0][i_gbl] < 0 ) g[0][i_gbl] = 0.0;
                if( g[1][i_gbl] > 0 ) g[1][i_gbl] = 0.0;
                if( i_t >= 1 ){
                    T dt = this->t[i_t] - this->t[i_t-1];
                    pos_int += 0.5 * dt * (g[0][i_gbl] + g[0][i_gbl-1]);
                    neg_int += 0.5 * dt * (g[1][i_gbl] + g[1][i_gbl-1]);
                }
            }
  
            T tol=1e-20;
            if( pos_int >= tol ){
                T pos_int_inv = 1.0 / pos_int;
                for(int i_t = 0; i_t < this->nt; i_t++){
                    int i_gbl = i_t + i_r * this->nt;
                    g[0][i_gbl] = pos_int_inv * g[0][i_gbl];
                }
            }
            else{
                for(int i_t = 0; i_t < this->nt; i_t++) {
                    int i_gbl = i_t + i_r * this->nt;
                    g[0][i_gbl] = 0.0;
                }
            }

            if( neg_int >= tol ){
                T neg_int_inv = 1.0 / neg_int;
                for(int i_t = 0; i_t < this->nt; i_t++){
                    int i_gbl = i_t + i_r * this->nt;
                    g[1][i_gbl] = neg_int_inv * g[1][i_gbl];
                }
            }
            else{
                for(int i_t = 0; i_t < this->nt; i_t++){
                    int i_gbl = i_t + i_r * this->nt;
                    g[1][i_gbl] = 0.0;
                }
            }
        }
        return g;
    }
};

#endif
