#ifndef WASS_SLICER_H
#define WASS_SLICER_H

#include "wass.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSlicer : public Wass<T>{
protected:
    int num_dists = 1;
    T tol = 0.0;
public:
    //inherit superclass constructor
    using Wass<T>::Wass;

    virtual T renorm_op(T x) = 0;

    void set_dists(int y) { this->num_dists = y; }
    void set_tol(T y) { this->tol = y; }

    //implement virtual renormalization routine from superclass
    valarray<valarray<T> > renormalize(valarray<T> f){
        //splitting maps to two pdfs
        valarray<valarray<T> > g(f,this->num_dists);

        //tmp var for efficiency
        int N = f.size();
 
        assert( N == this->nt or N == this->num_rec * this->nt );
       
        int num_rec_local = N / this->nt;

        for(int i_dists = 0; i_dists < num_dists; i_dists++){
            for(int i_r = 0; i_r < num_rec_local; i_r++){
                valarray<T> C(0.0, num_dists);
    
                //split between positive and negative parts
                for(int i_t = 0; i_t < this->nt; i_t++){
                    int i_gbl = i_t + i_r * this->nt;
                    g[i_dists][i_gbl] = renorm_op(g[i_dists][i_gbl]);
                    if( i_t >= 1 ){
                        T dt = this->t[i_t] - this->t[i_t-1];
                        C[i_dists] += 0.5 * dt * (g[i_dists][i_gbl] + g[i_dists][i_gbl-1]);
                    }
                }
            
                if( C[i_dists] >= this->tol ){
                    T C_inv = 1.0 / C[i_dists];
                    for(int i_t = 0; i_t < this->nt; i_t++){
                        int i_gbl = i_t + i_r * this->nt;
                        g[i_dists][i_gbl] = C_inv * g[i_dists][i_gbl];
                    }
                }
                else{
                    for(int i_t = 0; i_t < this->nt; i_t++){
                        int i_gbl = i_t + i_r * this->nt;
                        g[i_dists][i_gbl] = 0.0;
                    }
                }
            }
        }
        return g;
    }
};

#endif
