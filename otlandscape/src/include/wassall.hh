#ifndef WASS_ALL_H
#define WASS_ALL_H

#include "wass_normalizer.hh"
#include <vector>
#include <cassert>

//template<class T> using Ctn=vector<T>;
using namespace std;

template<class T>
class WassSplit : public WassNormalizer<T>{
private:
    int num_dists = 2;
public:
    using WassNormalizer<T>::WassNormalizer;

    T renorm_op(T x, int i_dist){
        if( i_dist == 0 ) return x >= 0 ? x : 0.0;
        else return x <= 0? -x : 0.0;
    }
};

template<class T>
class WassSquare : public WassNormalizer<T>{
private:
   int num_dists = 1;
public:
    //inherit superclass constructor
    using WassNormalizer<T>::WassNormalizer;

    T renorm_op(T x, int i_dist){
        return x * x;
    }
};

//template<class T> using Ctn=vector<T>;

template<class T>
class WassLinExp : public WassNormalizer<T>{
protected:
    T c = 1.0;
    T c_inv = 1.0;
public:
    //inherit superclass constructor
    using WassNormalizer<T>::WassNormalizer;
    
    void set_sharpness(T x) { this-> c = x; this->c_inv = 1.0 / this->c; }

    T renorm_op(T x, int i_dist){
        if( x >= 0 ) return x + c_inv;
        else return c_inv * exp(c * x);
    }
};

template<class T>
class WassExp : public WassNormalizer<T>{
protected:
    T c = 1.0;
public:
    //inherit superclass constructor
    using WassNormalizer<T>::WassNormalizer;

    void set_sharpness(T x) { this-> c = x; }

    T renorm_op(T x, int i_dist){
        return exp(c*x);
    }
};
#endif
