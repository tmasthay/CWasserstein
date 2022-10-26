#include "include/cub.hh"
#include "include/wassall.hh"
#include "include/sobolev.hh"
#include <vector>
#include <valarray>
#include <iostream>
#include <cassert>
#include <rsf.h>


template<typename T>
T wasstrace(const valarray<T> &f,
    const valarray<T> &g, 
    const valarray<T> &t,
    const valarray<T> &p,
    int mode,
    T c=1.0){

    //compute distances
    float value = 0.0;
    if( mode == 0 ){
        for(int i = 0; i < t.size(); i++)
            value += pow(f[i]-g[i], 2.0);
    }
    else if( mode >= 1 ){
       if( mode == 1 ){
           WassSplit2<float> my_misfit(g, t, p, 1);
           my_misfit.set_dists(2);
           value = my_misfit.eval(f);
       }
       else if( mode == 2 ){
           WassSquare<float> my_misfit(g, t, p, 1);
           value = my_misfit.eval(f);
       }
       else if( mode == 3 ){
           WassLinExp<float> my_misfit(g, t, p, 1);
           my_misfit.set_sharpness(c);
           value = my_misfit.eval(f);
       }
       else if( mode == 4 ){
           WassLin<float> my_misfit(g, t, p, nx);
           value = my_misfit.eval(f);
       }
       else if( mode == 5 ){
           WassExp<float> my_misfit(g, t, p, nx);
           value = my_misfit.eval(f);
       }
       else{
           cerr << "Mode " << mode << " not supported in wasstrace.cc\n";
           exit(-2);
       }
    }
    else{
       cerr << "Mode " << mode << " not supported in wasstrace.cc\n";
       exit(-2);
    }
    return value;
}
