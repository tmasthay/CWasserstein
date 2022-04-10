#ifndef WASS_H
#define WASS_H

#include "misfit.hh"
#include <iostream>
#include <cassert>
#include <vector>

using namespace std;
//template<class T> using Ctn=vector<T>;


template< class T >
class Sobolev : public Misfit<T>{
private:
    valarray<T> data;
    valarray<T> x;
    valarray<T> t;
 
    //s norm to take
    float s;
public:
    //later on, make this an abstract class by adding a renormalization routine here

    Wass(valarray<T> data, int nx, int nt, T dx, T dt, T ox){
        //set container class vars
        this->data = data;
        this->x = valarray<T>(0.0, nx);
        this->t = valarray<T>(0.0, nt);
 
        for(int ix = 0; ix < nx; ix++) this->x[ix] = ox + ix * dx;
        for(int it = 0; it < nt; it++) this->t[it] = it * dt;

        assert(this->data.size() == nx * nt);
        assert(this->x.size() == nx);
        assert(this->t.size() == nt); 
    }

    //implement slightly more general misfit for flexibility    
    T eval(const valarray<T> &m, T s_tilde){
        T total_sum = 0.0;
        if( s_tilde == 0.0 ){
            for(int ix = 0; ix < nx; ix++){
                for(int it = 0; it < nt; it++){
                    T diff = m[it + ix*nx] - data[it + ix*nx];
                    total_sum += diff*diff*dx*dt; 
                }
            }
        }
        else{
            for(int ix = 0; ix < nx; ix++){
                for(int it = 0; it < nt; it++){
                    T kernel = pow(1.0 + x[ix]*x[ix] + t[it]*t[it], s_tilde);
                    T diff = m[it + ix*nx] - data[it + ix*nx];
                    total_sum += diff*diff*kernel*dx*dt;
                }
            }
        } 
        return total_sum;
    }

    //implement virtual evaluation function -- workhorse function
    T eval(const valarray<T> &m){
        return eval(m, this->s);
    }
};

#endif
