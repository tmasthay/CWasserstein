#ifndef WASS_ALL_H
#define WASS_ALL_H

#include "wass_slicer.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSplit2 : public WassSlicer<T>{
public:
    using WassSlicer<T>::WassSlicer;

    this->num_dists = 2;

    T renorm_op(T x, int i_dist){
        if( i_dist == 0 ) return x >= 0 ? x : 0.0;
        else return x <= 0? -x : 0.0;
    }
}
template<class T>
class WassSquare : public WassSlicer<T>{
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    T renorm_op(T x, int i_dist){
        return x * x;
    }
};

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSquare : public WassSlicer<T>{
protected:
    T c = 1.0;
    T c_inv = 1.0;
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    void set_c(T x) { this-> c = x; this->c_inv = 1.0 / this->c; }

    T renorm_op(T x, int i_dist){
        if( x >= 0 ) return x + c_inv;
        else return c_inv * exp(c * x);
    }
};

#endif
