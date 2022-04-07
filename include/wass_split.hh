#ifndef WASS_SPLIT_H
#define WASS_SPLIT_H

#include "wass.hh"
#include <vector>
#include <cassert>

using namespace std;
//template<class T> using Ctn=vector<T>;

template<class T>
class WassSplit : public Wass<T>{
public:
    //inherit superclass constructor
    using Wass<T>::Wass;

    //implement virtual renormalization routine from superclass
    Ctn<Ctn<T>> renormalize(Ctn<T> f){
        //splitting maps to two pdfs
        Ctn<Ctn<T>> g(2,f);

        //tmp var for efficiency
        int N = f.size();
        
        //split between positive and negative parts
        for(int i = 0; i < N; i++){
            if( g[0][i] < 0 ) g[0][i] = 0.0;
            if( g[1][i] > 0 ) g[1][i] = 0.0;
        }
        return g;
    }
};

#endif
