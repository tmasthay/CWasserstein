#ifndef WASS_SPLIT_H
#define WASS_SPLIT_H

#include "wass.hh"
#include <vector>
#include <cassert>

using namespace std;
template<class T> using Ctn=vector<T>;

template<class T>
class WassSplit : public Wass<T>{
public:
    using Wass<T>::Wass;

    Ctn<Ctn<T>> renormalize(Ctn<T> f){
        Ctn<Ctn<T>> g(2,f);
        int N = f.size();
        for(int i = 0; i < N; i++){
            if( g.at(0).at(i) < 0 ) g.at(0).at(i) = 0.0;
            if( g.at(1).at(i) > 0 ) g.at(0).at(i) = 0.0;
        }
        return g;
    }
};

#endif
