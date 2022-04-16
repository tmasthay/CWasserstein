#ifndef SOBOLEV_H
#define SOBOLEV_H

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
    T s;

    //store X discretization
    int nx;
    T dx;
    T ox;
 
    //store T discretization
    int nt;
    T dt;
private:
    void renormalization(valarray<T> &input_data){
        T total_sum = 0.0;
        for(int i = 0; i < input_data.size(); i++){
            total_sum += input_data[i]*input_data[i];
        }
        T C = 0.0;
        const T eps = 1e-16;
        if( total_sum > eps ) C = 1.0 / total_sum;
        for(int i = 0; i < input_data.size(); i++){
            input_data[i] = C * input_data[i];
        } 
    }
public:
    //later on, make this an abstract class by adding a renormalization routine here

    Sobolev(valarray<T> &input_data, 
        T s, 
        int nx, 
        int nt, 
        T dx, 
        T dt, 
        T ox,
        bool renormalize=false){
        //set container class vars
        if( renormalize ){
            renormalization(input_data);
        }
        this->data = input_data;

        this->x = valarray<T>(0.0, nx);
        this->t = valarray<T>(0.0, nt);
        this->s = s;
 
        for(int ix = 0; ix < nx; ix++) this->x[ix] = ox + ix * dx;
        for(int it = 0; it < nt; it++) this->t[it] = it * dt;

        this->nx = nx;
        this->dx = dx;
        this->ox = ox;
 
        this->nt = nt;
        this->dt = dt;

        cerr << "data.size() == " << data.size() << " , nx == " << nx <<
            " , nt == " << nt << endl;
        assert(this->data.size() == nx * nt);
        assert(this->x.size() == nx);
        assert(this->t.size() == nt);  
    }

    //implement slightly more general misfit for flexibility    
    T eval(const valarray<T> &m, T s_tilde){
        assert( m.size() == nx * nt );
        T total_sum = 0.0;
        if( s_tilde == 0.0 ){
            for(int ix = 0; ix < nx; ix++){
                for(int it = 0; it < nt; it++){
                    T diff = m[it + ix*nt] - data[it + ix*nt];
                    //cerr << "diff^2 = " << diff*diff << endl;
                    total_sum += diff*diff*dx*dt;
                }
            }
        }
        else{
            for(int ix = 0; ix < nx; ix++){
                for(int it = 0; it < nt; it++){
                    T kernel = pow(1.0 + x[ix]*x[ix] + t[it]*t[it], s_tilde);
                    T diff = m[it + ix*nt] - data[it + ix*nt];
                    total_sum += diff*diff*kernel*dx*dt;
                }
            }
        } 
        return total_sum;
    }

    //implement virtual evaluation function -- workhorse function
    T eval(const valarray<T> &m){
        valarray<T> m2 = m;
        renormalization(m2);
        return eval(m2, this->s);
    }
};

#endif
