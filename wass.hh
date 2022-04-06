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
public:
    Wass(Ctn<T> data, Ctn<T> t, Ctn<T> p, int num_rec, int nt){
        this->data = data;
        assert(data.size() == num_rec * nt);
        this->t.assign(t.begin(), t.end());
        this->p.assign(p.begin(), p.end());
        this->num_rec = num_rec;
        this->nt = nt;
        this->data_ren.swap(this->renormalize(this->data));
    }
    Ctn<T> cdf(Ctn<T> f){
        //precondition error checking
        assert( f.size() == t.size() );

        int N = f.size();
        Ctn<T> F(N,0);
        for(int i = 0; i < N - 1; i++){
            T avgf = 0.5 * (f[i+1] + f[i]);
            T dt = t[i+1] - t[i];
            assert( dt > 0 );
            F[i+1] = F[i] + avgf * dt;
        }
        for(int i = 0; i < N; i++){
            F[i] /= F[N-1];
        }
        return F;
    }
    
    Ctn<T> quantile(Ctn<T> F){
         //precondition error checking
        // assert( F.size() == x.size() )
         //looping vars
         int i_t = 0;
    
         //initialize Quantile
         Ctn<T> Q(np, 0);
    
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
    
    T eval(Ctn<T> m){
       T eps_for_now = 1e-10;
       Ctn<Ctn<T>> m_ren(renormalize(m));
       assert(m_ren.size() == data_ren.size());
       assert(m_ren.at(0).size() == num_rec * nt);
       int num_splits = m_ren.size();
       T total_sum = 0.0;
       for(int i_s = 0; i_s < num_splits; i_s++){
           for(int i_r = 0; i_r < num_rec; i_r++){
               typename Ctn<T>::iterator start = m_ren.at(i_s).begin() + i_r * nt;
               typename Ctn<T>::iterator end = start + nt;
               Ctn<T> f_q(quantile(cdf(Ctn<T>(start,end))));
               Ctn<T> d_q(quantile(cdf(Ctn<T>(start,end))));
               for(int i_t = 0; i_t < nt - 1; i_t++){
                   T dt = t.at(i_t+1) - t.at(i_t);
                   T dqr = f_q.at(i_t+1) - d_q.at(i_t+1);
                   T dql = f_q.at(i_t) - d_q.at(i_t);
                   total_sum += 0.5 * dt * (dqr*dqr + dql*dql);
               }
           }
       }
       return total_sum;
    }
};

#endif
