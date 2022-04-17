#ifndef WASS_SQUARE_H
#define WASS_SQUARE_H

#include "wass_slicer.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSquare : public WassSlicer<T>{
public:
    //inherit superclass constructor
    using WassSlicer<T>::WassSlicer;

    T renorm_op(T x){
        return x * x;
    }
};

#endif
