#ifndef WASS_NORMALIZER_H
#define WASS_NORMALIZER_H

#include "wass.hh"
#include <iostream>
#include <cassert>
#include <vector>

//using namespace std;
//template<class T> using Ctn=vector<T>;


template< class T >
class WassNormalizer : public Wass<T>{
private:
    int num_splits;
public:
    using Wass<T>::Wass;

    void set_num_splits(const int &x) { num_splits = x; };

    //Define WassNormalizer Interface
    virtual T renorm_op(T x, int split);

    //Wass interface implementation #1
    valarray< valarray<T> > renormalize(const valarray<T> &f){
         //initialize output
         valarray< valarray<T> > g = valarray< valarray<T> >(
             valarray<T>(0.0, f.size()), num_splits);

         for(int split = 0; split < num_splits; split++)
             for(int i = 0; i < f.size(); i++)
                 g[split][i] = renorm_op(f[i], split);
        return g;
    }
};

#endif
