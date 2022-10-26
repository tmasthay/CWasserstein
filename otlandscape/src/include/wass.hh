#ifndef WASS_H
#define WASS_H

#include "misfit.hh"
#include <iostream>
#include <cassert>
#include <vector>

using namespace std;
//template<class T> using Ctn=vector<T>;


template< class T >
class Wass : public Misfit<T>{
protected:
    valarray<T> data;
    valarray< valarray<T> > data_ren;
    valarray< valarray<T> > quantile_var;
    valarray<T> t;
    valarray<T> p;
    int num_rec;
    int curr_trace;
    int nt;
    int np;
    T eps=1e-10;
    bool init_ren=false;
    bool init_quantile=false;
public:
    Wass(const &valarray<T> data, 
        const &valarray<T> t, 
        const &valarray<T> p){ 
        //set container class vars
        this->data = data;
        this->t = t;
        this->p = p;

        //set counter vars
        this->nt = t.size();
        this->np = p.size();
    }

    //pure virtual function
    virtual valarray< valarray<T> > renormalize(const &valarray<T> f) = 0;

    valarray<T> cdf(const &valarray<T> f){
        //precondition error checking
        assert( f.size() == t.size() );

        //set size parameter for small speedup
        int N = f.size() + 1;

        //initialize cdf
        valarray<T> F(0.0,N);

        //integrate f to recover F
        for(int i = 1; i < N; i++){
            T dt = t[i] - t[i-1];
            F[i] = F[i-1] + dt * f[i-1];
        }

        //normalize so it's a valid cdf
        for(int i = 0; i < N; i++){
           F[i] /= F[N-1];
        }
        return F;
    }
    
    //quantile function 
    valarray<T> quantile(const &valarray<T> F){
         //precondition error checking
         assert( F.size() == nt );
    
         //initialize Quantile
         valarray<T> Q(0.0, np);
  
         //initialize loop variable
         int i_t = 0;
         //loop over all probability values p[i_p])
         for(int i_p = 0; i_p < np; i_p++){
             while(i_t < nt - 1){
                 if( F[i_t] <= p[i_p] and p[i_p] <= F[i_t+1] ){
                     T df = F[i_t+1] - F[i_t];
                     if( df < eps ){
                         //left bin -- division might cause problems
                         Q[i_p] = t[i_t];
                     }
                     else{ //else interpolate between left and right
                         T alpha = (p[i_p] - F[i_t]) / df;
                         Q[i_p] = (1-alpha) * t[i_t] 
                             + alpha * t[i_t+1];
                     }
                     break;
                 }
                 i_t++;
             }
         }
         return Q;
    }
 
    valarray< valarray<T> > cdf_multi(
        const &valarray< valarray<T> >){
        valarray< valarray<T> > F = valarray< valarray<T> >(f);
        for(int i = 0; i < F.size(); i++){
            F[i] = cdf(f[i]);
        }
        return F;
    }

    valarray< valarray<T> > quantile_multi(
        const &valarray< valarray<T> >){
        valarray< valarray<T> > F = valarray< valarray<T> >(f);
        for(int i = 0; i < F.size(); i++){
            F[i] = quantile(f[i]);
        }
        return F;
    }
    
    //implement virtual evaluation function -- workhorse function
    T eval(const valarray<T> &m){
       //renormalize synthetic data
       const valarray<valarray<T> > m_ren(renormalize(m));

       //renormalize data if needed
       if( not init_ren ){
           data_ren = renormalize(data);
           init_ren = true;
       }

       
       //sanity check dimensions
       assert(m_ren.size() == data_ren.size());

       //num probability densities 
       int num_splits = m_ren.size();

       if( not init_quantile ){
           init_quantile=true;
           //come back and implement so that it's not redundant
           quantile_var = quantile_multi(data);
       }

       valarray< valarray<T> > q = quantile_multi(m_ren);

       assert(q.size() == quantile_var.size());
       assert(q[0].size() == np and quantile_var[0].size() == np);

       T sum = 0.0;
       for(int i = 0; i < q.size(); i++){
           for(int i_p = 0; i_p < np - 1; i_p++){
               T dp = p[i_p+1] - p[i_p];
               T dqr = quantile_var[i_p+1] - q[i_p+1];
               T dql = quantile_var[i_p] - q[i_p];
               sum += 0.5 * dp * (dqr*dqr + dql*dql);
           }
       }
    return sum;
};

#endif
