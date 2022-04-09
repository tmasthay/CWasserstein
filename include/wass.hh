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
private:
    valarray<T> data;
    valarray<valarray<T> > data_ren;
    valarray<T> quantile_var;
    valarray<T> t;
    valarray<T> p;
    int num_rec;
    int nt;
    int np;
    T eps=1e-10;
    bool init_ren=false;
    bool init_quantile=false;
public:
    Wass(valarray<T> data, valarray<T> t, valarray<T> p, int num_rec){
        //set container class vars
        this->data = data;
        this->t = t;
        this->p = p;

        //set counter vars
        this->nt = t.size();
        this->num_rec = num_rec;
        this->np = p.size();

        //sanity check
        assert(data.size() == num_rec * nt);
    }
    valarray<T> cdf(valarray<T> f){
        //precondition error checking
        assert( f.size() == t.size() );

        //set size parameter for small speedup
        int N = f.size();

        //initialize cdf
        valarray<T> F(N,0);

        //integrate f to recover F
        for(int i = 0; i < N - 1; i++){
            T avgf = 0.5 * (f[i+1] + f[i]);
            T dt = t[i+1] - t[i];
            assert( dt > 0 );
            F[i+1] = F[i] + avgf * dt;
        }
 
        //normalize so it's a valid cdf
        for(int i = 0; i < N; i++){
            F[i] /= F[N-1];
        }
        return F;
    }
    
    //quantile function 
    valarray<T> quantile(valarray<T> F){
         //precondition error checking
         assert( F.size() == nt );
    
         //initialize Quantile
         valarray<T> Q(np, 0);
    
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
   
    
    virtual valarray<valarray<T> > renormalize(valarray<T> f) = 0; // pure virtual
    
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
       assert(m_ren[0].size() == num_rec * nt);

       //num_splits = k <--> renormalization routine outputs k prob. densities
       int num_splits = m_ren.size();

       //total_sum for return value
       T total_sum = 0.0;

       if( not init_quantile ){
           init_quantile=true;
           //come back and implement so that it's not redundant
       }

       for(int i_s = 0; i_s < num_splits; i_s++){
           for(int i_r = 0; i_r < num_rec; i_r++){
               //define slice parameters 
               size_t lengths[] = {nt};
               size_t strides[]= {1};

               //define slice object for i_s-th split and i_r-th trace, zero-indexed
               gslice trace_indices(i_r * nt,
                   valarray<size_t>(lengths, 1),
                   valarray<size_t>(strides, 1));

               //extract data along slice
               const valarray<valarray<T> > m_ren_tmp(m_ren);
               const valarray<valarray<T> > d_ren_tmp(data_ren);

               valarray<T> m_trace = m_ren_tmp[i_s][trace_indices];
               valarray<T> d_trace = d_ren_tmp[i_s][trace_indices];

               //compute quantile functions --> note: we should really only compute quantile for
               //    d_q once! --> come back to this.
               valarray<T> f_q(quantile(cdf(m_trace)));
               valarray<T> d_q(quantile(cdf(d_trace)));
                
               //print some info for debugging -- superfluous
               for(int i_p = 0; i_p < np; i_p++){
//                   cout << "(" << f_q[i_p]) < "," << d_q[i_p]) << ")" << endl;
                   cout << "(" << p[i_p] << "," << f_q[i_p] << ", " << 
                       d_q[i_p] << ")" << endl;
               }
 
               //perform sanity checks
               assert( f_q.size() == np );
               //assert( d_q.size() == np );
  
               //sum up squared difference of quantiles
               for(int i_p = 0; i_p < np - 1; i_p++){
                   T dp = p[i_p+1] - p[i_p];
                   T dqr = f_q[i_p+1] - d_q[i_p+1];
                   T dql = f_q[i_p] - d_q[i_p];
                   total_sum += 0.5 * dp * (dqr*dqr + dql*dql);
               }
           }
       }
       return total_sum;
    }
};

#endif
