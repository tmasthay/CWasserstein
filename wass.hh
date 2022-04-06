#ifndef WASS_H
#define WASS_H

#include "misfit.hh"
#include <iostream>
#include <cassert>
#include <vector>

using namespace std;
template<class T> using Ctn=vector<T>;


template< class T >
class Wass : public Misfit<T>{
private:
    Ctn<T> data;
    Ctn<Ctn<T>> data_ren;
    Ctn<T> t;
    Ctn<T> p;
    int num_rec;
    int nt;
    int np;
    T eps=1e-10;
    bool init_ren=false;
public:
    Wass(Ctn<T> data, Ctn<T> t, Ctn<T> p, int num_rec){
        //set container class vars
        this->data = data;
        this->t.assign(t.begin(), t.end());
        this->p.assign(p.begin(), p.end());

        //set counter vars
        this->nt = t.size();
        this->num_rec = num_rec;
        this->np = p.size();

        //sanity check
        assert(data.size() == num_rec * nt);
    }
    Ctn<T> cdf(Ctn<T> f){
        //precondition error checking
        assert( f.size() == t.size() );

        //set size parameter for small speedup
        int N = f.size();

        //initialize cdf
        Ctn<T> F(N,0);

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
    Ctn<T> quantile(Ctn<T> F){
         //precondition error checking
         assert( F.size() == nt );
    
         //initialize Quantile
         Ctn<T> Q(np, 0);
    
         //loop over all probability values p.at(i_p)
         for(int i_p = 0; i_p < np; i_p++){
             for(int i_t = 0; i_t < nt - 1; i_t++){
                 if( F.at(i_t) <= p.at(i_p) and p.at(i_p) <= F.at(i_t+1) ){
                     T df = F[i_t+1] - F[i_t];
                     if( df < eps ){
                         //left bin -- division might cause problems
                         Q.at(i_p) = F.at(i_t);
                     }
                     else{ //else interpolate between left and right
                         T alpha = (p.at(i_p) - F.at(i_t)) / df;
                         Q.at(i_p) = (1-alpha) * F.at(i_t) + alpha * F.at(i_t+1);
                     }
                 }
             }
         }
         return Q;
    }
   
    
    virtual Ctn<Ctn<T>> renormalize(Ctn<T> f) = 0; // pure virtual
    
    //implement virtual evaluation function -- workhorse function
    T eval(Ctn<T> m){
       //renormalize synthetic data
       Ctn<Ctn<T>> m_ren(renormalize(m));

       //renormalize data if needed
       if( not init_ren ){
           data_ren = renormalize(data);
           init_ren = true;
       }

       //sanity check dimensions
       assert(m_ren.size() == data_ren.size());
       assert(m_ren.at(0).size() == num_rec * nt);

       //num_splits = k <--> renormalization routine outputs k prob. densities
       int num_splits = m_ren.size();
       cerr << "numsplits = " << num_splits << endl;
       cerr << "num_rec = " << num_rec << endl;

       //total_sum for return value
       T total_sum = 0.0;
       for(int i_s = 0; i_s < num_splits; i_s++){
           for(int i_r = 0; i_r < num_rec; i_r++){
               //get proper index slices
               typename Ctn<T>::iterator start = m_ren.at(i_s).begin() + i_r * nt;
               typename Ctn<T>::iterator end = start + nt;

               //recover quantile functions
               Ctn<T> f_q(quantile(cdf(Ctn<T>(start,end))));
               Ctn<T> d_q(quantile(cdf(Ctn<T>(start,end))));
 
               assert( f_q.size() == np );
               assert( d_q.size() == np );
  
               //sum up squared difference of quantiles
               for(int i_p = 0; i_p < np - 1; i_p++){
                   T dp = p.at(i_p+1) - p.at(i_p);
                   T dqr = f_q.at(i_p+1) - d_q.at(i_p+1);
                   T dql = f_q.at(i_p) - d_q.at(i_p);
                   total_sum += 0.5 * dp * (dqr*dqr + dql*dql);
               }
           }
       }
       return total_sum;
    }
};

#endif
